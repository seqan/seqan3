// zipstream Library License:
// --------------------------
//
// The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.
//
// This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
//
// 3. This notice may not be removed or altered from any source distribution
//
//
// Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003   (original zlib stream)
// Author: David Weese, dave.weese@gmail.com, 2014             (extension to parallel block-wise compression in bgzf format)
// Author: Rene Rahn, rene.rahn AT fu-berlin.de, 2019          (adaption to SeqAn library version 3)

#pragma once

#include <cstring>
#include <condition_variable>
#include <mutex>

#include <seqan3/contrib/parallel/buffer_queue.hpp>
#include <seqan3/contrib/stream/bgzf_stream_util.hpp>
#include <seqan3/utility/parallel/detail/reader_writer_manager.hpp>

namespace seqan3::contrib
{

// ===========================================================================
// Classes
// ===========================================================================

// --------------------------------------------------------------------------
// Class basic_bgzf_istreambuf
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istreambuf : public std::basic_streambuf<Elem, Tr>
{
public:

    typedef Tr                              traits_type;
    typedef typename traits_type::char_type char_type;
    typedef typename traits_type::int_type  int_type;
    typedef typename traits_type::off_type  off_type;
    typedef typename traits_type::pos_type  pos_type;

private:

    typedef std::basic_istream<Elem, Tr>& istream_reference;

    typedef ElemA                   char_allocator_type;
    typedef ByteT                   byte_type;
    typedef ByteAT                  byte_allocator_type;
    typedef byte_type*              byte_buffer_type;

    typedef std::vector<char_type, char_allocator_type> TBuffer;
    typedef fixed_buffer_queue<int32_t>                 TJobQueue;

    static const size_t MAX_PUTBACK = 4;

    // Allows serialized access to the underlying buffer.
    struct Serializer
    {
        istream_reference   istream;
        std::mutex          lock;
        io_error            *error;
        off_type            fileOfs;

        Serializer(istream_reference istream) :
            istream(istream),
            error(NULL),
            fileOfs(0u)
        {}

        ~Serializer()
        {
            delete error;
        }
    };

    Serializer serializer;

    struct DecompressionJob
    {
        typedef std::vector<byte_type, byte_allocator_type> TInputBuffer;

        TInputBuffer            inputBuffer;
        TBuffer                 buffer;
        off_type                fileOfs;
        int32_t                 size;
        uint32_t                compressedSize;

        std::mutex              cs;
        std::condition_variable readyEvent;
        bool                    ready;
        bool                    bgzfEofMarker;

        DecompressionJob() :
            inputBuffer(DefaultPageSize<detail::bgzf_compression>::MAX_BLOCK_SIZE, 0),
            buffer(MAX_PUTBACK + DefaultPageSize<detail::bgzf_compression>::MAX_BLOCK_SIZE / sizeof(char_type), 0),
            fileOfs(),
            size(0),
            cs(),
            readyEvent(),
            ready(true),
            bgzfEofMarker(false)
        {}

        DecompressionJob(DecompressionJob const &other) :
            inputBuffer(other.inputBuffer),
            buffer(other.buffer),
            fileOfs(other.fileOfs),
            size(other.size),
            cs(),
            readyEvent(),
            ready(other.ready),
            bgzfEofMarker(other.bgzfEofMarker)
        {}
    };

    // string of recyclable jobs
    size_t                        numThreads;
    size_t                        numJobs;
    std::vector<DecompressionJob> jobs;
    TJobQueue                     runningQueue;
    TJobQueue                     todoQueue;
    detail::reader_writer_manager runningQueueManager;  // synchronises reader, writer with running queue.
    detail::reader_writer_manager todoQueueManager;    // synchronises reader, writer with todo queue.
    int                           currentJobId;

    struct DecompressionThread
    {
        basic_bgzf_istreambuf                        *streamBuf;
        CompressionContext<detail::bgzf_compression> compressionCtx;

        void operator()()
        {
            // Active reader to consume from todo queue.
            auto reader_raii = streamBuf->todoQueueManager.register_reader();
            // Active writer to produce work for the decompression queue.
            auto writer_raii = streamBuf->runningQueueManager.register_writer();

            // wait for a new job to become available
            while (true)
            {

                int jobId = -1;
                if (streamBuf->todoQueue.wait_pop(jobId) == queue_op_status::closed)
                    return;

                DecompressionJob &job = streamBuf->jobs[jobId];
                size_t tailLen = 0;

                // typically the idle queue contains only ready jobs
                // however, if seek() fast forwards running jobs into the todoQueue
                // the caller defers the task of waiting to the decompression threads
                if (!job.ready)
                {
                    std::unique_lock<std::mutex> lock(job.cs);
                    job.readyEvent.wait(lock, [&job]{return job.ready;});
                    assert(job.ready == true);
                }

                {
                    std::lock_guard<std::mutex> scopedLock(streamBuf->serializer.lock);

                    job.bgzfEofMarker = false;
                    if (streamBuf->serializer.error != NULL)
                        return;

                    // remember start offset (for tellg later)
                    job.fileOfs = streamBuf->serializer.fileOfs;
                    job.size = -1;
                    job.compressedSize = 0;

                    // only load if not at EOF
                    if (job.fileOfs != -1)
                    {
                        // read header
                        streamBuf->serializer.istream.read(
                            (char_type*)&job.inputBuffer[0],
                            DefaultPageSize<detail::bgzf_compression>::BLOCK_HEADER_LENGTH);

                        if (!streamBuf->serializer.istream.good())
                        {
                            streamBuf->serializer.fileOfs = -1;
                            if (streamBuf->serializer.istream.eof())
                                goto eofSkip;
                            streamBuf->serializer.error = new io_error("Stream read error.");
                            return;
                        }

                        // check header
                        if (!detail::bgzf_compression::validate_header(std::span{job.inputBuffer}))
                        {
                            streamBuf->serializer.fileOfs = -1;
                            streamBuf->serializer.error = new io_error("Invalid BGZF block header.");
                            return;
                        }

                        // extract length of compressed data
                        tailLen = _bgzfUnpack16(&job.inputBuffer[0] + 16) +
                                  1u - DefaultPageSize<detail::bgzf_compression>::BLOCK_HEADER_LENGTH;

                        // read compressed data and tail
                        streamBuf->serializer.istream.read(
                            (char_type*)&job.inputBuffer[0] + DefaultPageSize<detail::bgzf_compression>::BLOCK_HEADER_LENGTH,
                            tailLen);

                        // Check if end-of-file marker is set
                        if (memcmp(reinterpret_cast<uint8_t const *>(&job.inputBuffer[0]),
                                   reinterpret_cast<uint8_t const *>(&BGZF_END_OF_FILE_MARKER[0]),
                                   28) == 0)
                        {
                            job.bgzfEofMarker = true;
                        }

                        if (!streamBuf->serializer.istream.good())
                        {
                            streamBuf->serializer.fileOfs = -1;
                            if (streamBuf->serializer.istream.eof())
                                goto eofSkip;
                            streamBuf->serializer.error = new io_error("Stream read error.");
                            return;
                        }

                        job.compressedSize = DefaultPageSize<detail::bgzf_compression>::BLOCK_HEADER_LENGTH + tailLen;
                        streamBuf->serializer.fileOfs += job.compressedSize;
                        job.ready = false;

                    eofSkip:
                        streamBuf->serializer.istream.clear(
                            streamBuf->serializer.istream.rdstate() & ~std::ios_base::failbit);
                    }

                    if (streamBuf->runningQueue.try_push(jobId) != queue_op_status::success)
                    {
                        // signal that job is ready
                        {
                            std::unique_lock<std::mutex> lock(job.cs);
                            job.ready = true;
                        }
                        job.readyEvent.notify_all();
                        return;  // Terminate this thread.
                    }
                }

                if (!job.ready)
                {
                    // decompress block
                    job.size = _decompressBlock(
                        &job.buffer[0] + MAX_PUTBACK, job.buffer.capacity(),
                        &job.inputBuffer[0], job.compressedSize, compressionCtx);

                    // signal that job is ready
                    {
                        std::unique_lock<std::mutex> lock(job.cs);
                        job.ready = true;
                    }
                    job.readyEvent.notify_all();
                }
            }
        }
    };

    std::vector<std::thread> pool;  // pool of worker threads
    TBuffer                  putbackBuffer;

public:

    basic_bgzf_istreambuf(istream_reference istream_,
                          size_t numThreads = bgzf_thread_count,
                          size_t jobsPerThread = 8) :
        serializer(istream_),
        numThreads(numThreads),
        numJobs(numThreads * jobsPerThread),
        runningQueue(numJobs),
        todoQueue(numJobs),
        runningQueueManager(detail::reader_count{1}, detail::writer_count{numThreads}, runningQueue),
        todoQueueManager(detail::reader_count{numThreads}, detail::writer_count{1}, todoQueue),
        putbackBuffer(MAX_PUTBACK)
    {
        jobs.resize(numJobs);
        currentJobId = -1;

        // Prepare todo queue.
        for (size_t i = 0; i < numJobs; ++i)
        {
            [[maybe_unused]] queue_op_status status = todoQueue.try_push(i);
            assert(status == queue_op_status::success);
        }

        // Start off the threads.
        for (size_t i = 0; i < numThreads; ++i)
            pool.emplace_back(DecompressionThread{this, CompressionContext<detail::bgzf_compression>{}});
    }

    ~basic_bgzf_istreambuf()
    {
        // Signal todoQueue that no more work is coming and close todo queue.
        todoQueueManager.writer_arrive();

        // Wait for threads to finish there active work.
        for (auto & t : pool)
        {
            if (t.joinable())
                t.join();
        }

        // Signal running queue that reader is done.
        runningQueueManager.reader_arrive();
    }

    int_type underflow()
    {
        // no need to use the next buffer?
        if (this->gptr() && this->gptr() < this->egptr())
            return traits_type::to_int_type(*this->gptr());

        size_t putback = this->gptr() - this->eback();
        if (putback > MAX_PUTBACK)
            putback = MAX_PUTBACK;

        // save at most MAX_PUTBACK characters from previous page to putback buffer
        if (putback != 0)
            std::copy(
                this->gptr() - putback,
                this->gptr(),
                &putbackBuffer[0]);

        if (currentJobId >= 0)
            todoQueue.wait_push(currentJobId);
            // appendValue(todoQueue, currentJobId);

        while (true)
        {
            if (runningQueue.wait_pop(currentJobId) == queue_op_status::closed)
            {
                currentJobId = -1;
                assert(serializer.error != NULL);
                if (serializer.error != NULL)
                    throw *serializer.error;
                return EOF;
            }

            DecompressionJob &job = jobs[currentJobId];

            // restore putback buffer
            this->setp(&job.buffer[0], &job.buffer[0] + (job.buffer.size() - 1));
            if (putback != 0)
                std::copy(
                    &putbackBuffer[0],
                    &putbackBuffer[0] + putback,
                    &job.buffer[0] + (MAX_PUTBACK - putback));

            // wait for the end of decompression
            {
                std::unique_lock<std::mutex> lock(job.cs);
                job.readyEvent.wait(lock, [&job]{return job.ready;});
            }

            size_t size = (job.size != -1)? job.size : 0;

            // reset buffer pointers
            this->setg(
                  &job.buffer[0] + (MAX_PUTBACK - putback),     // beginning of putback area
                  &job.buffer[0] + MAX_PUTBACK,                 // read position
                  &job.buffer[0] + (MAX_PUTBACK + size));       // end of buffer

            // The end of the bgzf file is reached, either if there was an error, or if the
            // end-of-file marker was reached, while the uncompressed block had zero size.
            if (job.size == -1 || (job.size == 0 && job.bgzfEofMarker))
                return EOF;
            else if (job.size > 0)
                return traits_type::to_int_type(*this->gptr());      // return next character

            throw io_error("BGZF: Invalid end condition in decompression. "
                           "Most likely due to an empty bgzf block without end-of-file marker.");
        }
    }

    pos_type seekoff(off_type ofs, std::ios_base::seekdir dir, std::ios_base::openmode openMode)
    {
        if ((openMode & (std::ios_base::in | std::ios_base::out)) == std::ios_base::in)
        {
            if (dir == std::ios_base::cur && ofs >= 0)
            {
                // forward delta seek
                while (currentJobId < 0 || this->egptr() - this->gptr() < ofs)
                {
                    ofs -= this->egptr() - this->gptr();
                    if (this->underflow() == static_cast<int_type>(EOF))
                        break;
                }

                if (currentJobId >= 0 && ofs <= this->egptr() - this->gptr())
                {
                    DecompressionJob &job = jobs[currentJobId];

                    // reset buffer pointers
                    this->setg(
                          this->eback(),            // beginning of putback area
                          this->gptr() + ofs,       // read position
                          this->egptr());           // end of buffer

                    if (this->gptr() != this->egptr())
                        return pos_type((job.fileOfs << 16) + ((this->gptr() - &job.buffer[MAX_PUTBACK])));
                    else
                        return pos_type((job.fileOfs + job.compressedSize) << 16);
                }

            }
            else if (dir == std::ios_base::beg)
            {
                // random seek
                std::streampos destFileOfs = ofs >> 16;

                // are we in the same block?
                if (currentJobId >= 0 && jobs[currentJobId].fileOfs == (off_type)destFileOfs)
                {
                    DecompressionJob &job = jobs[currentJobId];

                    // reset buffer pointers
                    this->setg(
                          this->eback(),                                        // beginning of putback area
                          &job.buffer[0] + (MAX_PUTBACK + (ofs & 0xffff)),      // read position
                          this->egptr());                                       // end of buffer
                    return ofs;
                }

                // ok, different block
                {
                    std::lock_guard<std::mutex> scopedLock(serializer.lock);

                    // remove all running jobs and put them in the idle queue unless we
                    // find our seek target

                    if (currentJobId >= 0)
                        todoQueue.wait_push(currentJobId);

                    // Note that if we are here the current job does not represent the sought block.
                    // Hence if the running queue is empty we need to explicitly unset the jobId,
                    // otherwise we would not update the serializers istream pointer to the correct position.
                    if (runningQueue.is_empty())
                        currentJobId = -1;

                    // empty is thread-safe in serializer.lock
                    while (!runningQueue.is_empty())
                    {
                        runningQueue.wait_pop(currentJobId);

                        if (jobs[currentJobId].fileOfs == (off_type)destFileOfs)
                            break;

                        // push back useless job
                        todoQueue.wait_push(currentJobId);
                        currentJobId = -1;
                    }

                    if (currentJobId == -1)
                    {
                        assert(runningQueue.is_empty());
                        serializer.istream.clear(serializer.istream.rdstate() & ~std::ios_base::eofbit);
                        if (serializer.istream.rdbuf()->pubseekpos(destFileOfs, std::ios_base::in) == destFileOfs)
                            serializer.fileOfs = destFileOfs;
                        else
                            currentJobId = -2;      // temporarily signals a seek error
                    }
                }

                // if our block wasn't in the running queue yet, it should now
                // be the first that falls out after modifying serializer.fileOfs
                if (currentJobId == -1)
                    runningQueue.wait_pop(currentJobId);
                else if (currentJobId == -2)
                    currentJobId = -1;

                if (currentJobId >= 0)
                {
                    // wait for the end of decompression
                    DecompressionJob &job = jobs[currentJobId];

                    {
                        std::unique_lock<std::mutex> lock(job.cs);
                        job.readyEvent.wait(lock, [&job]{return job.ready;});
                    }

                    assert(job.fileOfs == (off_type)destFileOfs);

                    // reset buffer pointers
                    this->setg(
                          &job.buffer[0] + MAX_PUTBACK,                     // no putback area
                          &job.buffer[0] + (MAX_PUTBACK + (ofs & 0xffff)),  // read position
                          &job.buffer[0] + (MAX_PUTBACK + job.size));       // end of buffer
                    return ofs;
                }
            }
        }
        return pos_type(off_type(-1));
    }

    pos_type seekpos(pos_type pos, std::ios_base::openmode openMode)
    {
        return seekoff(off_type(pos), std::ios_base::beg, openMode);
    }

    // returns the compressed input istream
    istream_reference get_istream()    { return serializer.istream; };
};

// --------------------------------------------------------------------------
// Class basic_bgzf_istreambase
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
    typedef std::basic_istream<Elem, Tr>&                          istream_reference;
    typedef basic_bgzf_istreambuf<Elem, Tr, ElemA, ByteT, ByteAT>  decompression_bgzf_streambuf_type;

    basic_bgzf_istreambase(istream_reference istream_)
        : m_buf(istream_)
    {
        this->init(&m_buf);
    };

    // returns the underlying decompression bgzf istream object
    decompression_bgzf_streambuf_type* rdbuf() { return &m_buf; };

    // returns the bgzf error state
    int get_zerr() const                    { return m_buf.get_zerr(); };
    // returns the uncompressed data crc
    long get_crc() const                    { return m_buf.get_crc(); };
    // returns the uncompressed data size
    long get_out_size() const               { return m_buf.get_out_size(); };
    // returns the compressed data size
    long get_in_size() const                { return m_buf.get_in_size(); };

private:
    decompression_bgzf_streambuf_type m_buf;
};

// --------------------------------------------------------------------------
// Class basic_bgzf_istream
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istream :
    public basic_bgzf_istreambase<Elem,Tr,ElemA,ByteT,ByteAT>,
    public std::basic_istream<Elem,Tr>
{
public:
    typedef basic_bgzf_istreambase<Elem,Tr,ElemA,ByteT,ByteAT> bgzf_istreambase_type;
    typedef std::basic_istream<Elem,Tr>                        istream_type;
    typedef istream_type &                                     istream_reference;
    typedef char                                               byte_type;

    basic_bgzf_istream(istream_reference istream_) :
        bgzf_istreambase_type(istream_),
        istream_type(bgzf_istreambase_type::rdbuf()),
        m_is_gzip(false),
        m_gbgzf_data_size(0)
    {};

    // returns true if it is a gzip file
    bool is_gzip() const                { return m_is_gzip; };
    // return data size check
    bool check_data_size() const        { return this->get_out_size() == m_gbgzf_data_size; };

    // return the data size in the file
    long get_gbgzf_data_size() const    { return m_gbgzf_data_size; };

protected:
    static void read_long(istream_reference in_, unsigned long& x_);

    int check_header();
    bool m_is_gzip;
    unsigned long m_gbgzf_data_size;

#ifdef _WIN32
private:
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
};

// ===========================================================================
// Typedefs
// ===========================================================================

// A typedef for basic_bgzf_istream<char>
typedef basic_bgzf_istream<char> bgzf_istream;
// A typedef for basic_bgzf_istream<wchart>
typedef basic_bgzf_istream<wchar_t> bgzf_wistream;

} // namespace seqan3::conrib
