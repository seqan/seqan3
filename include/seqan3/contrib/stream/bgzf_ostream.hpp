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
// Author: René Rahn, rene.rahn [at] fu-berlin.de, 2019        (adaptions to SeqAn library version 3)

#pragma once

#include <seqan3/contrib/parallel/serialised_resource_pool.hpp>
#include <seqan3/contrib/parallel/suspendable_queue.hpp>
#include <seqan3/contrib/stream/bgzf_stream_util.hpp>

namespace seqan3::contrib
{

// --------------------------------------------------------------------------
// Class basic_bgzf_ostreambuf
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_ostreambuf : public std::basic_streambuf<Elem, Tr>
{
private:

    typedef std::basic_ostream<Elem, Tr>&                ostream_reference;
    typedef ElemA                                        char_allocator_type;
    typedef ByteT                                        byte_type;
    typedef ByteAT                                       byte_allocator_type;
    typedef byte_type*                                   byte_buffer_type;
    typedef ConcurrentQueue<size_t, Suspendable<Limit> > job_queue_type;

public:

    typedef Tr                                traits_type;
    typedef typename traits_type::char_type   char_type;
    typedef typename traits_type::int_type    int_type;
    typedef typename traits_type::pos_type    pos_type;
    typedef typename traits_type::off_type    off_type;

    struct ScopedLock
    {
        ScopedLock(std::function<void()> complete_fn) : completion(std::move(complete_fn))
        {}

        ~ScopedLock()
        {
            completion();
        }

        std::function<void()> completion;
    };

    // One compressed block.
    struct OutputBuffer
    {
        char    buffer[DefaultPageSize<detail::bgzf_compression>::MAX_BLOCK_SIZE];
        size_t  size;
    };

    // Writes the output to the underlying stream when invoked.
    struct BufferWriter
    {
        ostream_reference ostream;

        BufferWriter(ostream_reference ostream) :
            ostream(ostream)
        {}

        bool operator() (OutputBuffer const & outputBuffer)
        {
            ostream.write(outputBuffer.buffer, outputBuffer.size);
            return ostream.good();
        }
    };

    struct CompressionJob
    {
        typedef std::vector<char_type, char_allocator_type> TBuffer;

        TBuffer         buffer;
        size_t          size;
        OutputBuffer    *outputBuffer;

        CompressionJob() :
            buffer(DefaultPageSize<detail::bgzf_compression>::VALUE / sizeof(char_type), 0),
            size(0),
            outputBuffer(NULL)
        {}
    };

    // string of recycable jobs
    size_t                                 numThreads;
    size_t                                 numJobs;
    std::vector<CompressionJob>            jobs;
    job_queue_type                         jobQueue;
    job_queue_type                         idleQueue;
    Serializer<OutputBuffer, BufferWriter> serializer;
    size_t                                 currentJobId;
    bool                                   currentJobAvail;

    struct CompressionThread
    {
        basic_bgzf_ostreambuf                        *streamBuf;
        CompressionContext<detail::bgzf_compression> compressionCtx;

        void operator()()
        {
            ScopedLock readLock{[this] () mutable { unlockReading(this->streamBuf->jobQueue); }};
            // ScopedReadLock<TJobQueue> readLock(streamBuf->jobQueue);
            ScopedLock writeLock{[this] () mutable { unlockWriting(this->streamBuf->idleQueue); }};
            // ScopedWriteLock{obQueue> writeLock{str}amBuf->idleQueue);

            // wait for a new job to become available
            bool success = true;
            while (success)
            {
                size_t jobId = -1;
                if (!popFront(jobId, streamBuf->jobQueue))
                    return;

                CompressionJob &job = streamBuf->jobs[jobId];

                // compress block with zlib
                job.outputBuffer->size = _compressBlock(
                    job.outputBuffer->buffer, sizeof(job.outputBuffer->buffer),
                    &job.buffer[0], job.size, compressionCtx);

                success = releaseValue(streamBuf->serializer, job.outputBuffer);
                appendValue(streamBuf->idleQueue, jobId);
            }
        }
    };

    // array of worker threads
    // using TFuture = decltype(std::async(CompressionThread{nullptr, CompressionContext<BgzfFile>{}, static_cast<size_t>(0)}));
    std::vector<std::thread>  pool;

    basic_bgzf_ostreambuf(ostream_reference ostream_,
                         size_t numThreads = bgzf_thread_count,
                         size_t jobsPerThread = 8) :
        numThreads(numThreads),
        numJobs(numThreads * jobsPerThread),
        jobQueue(numJobs),
        idleQueue(numJobs),
        serializer(ostream_, numThreads * jobsPerThread)
    {
        jobs.resize(numJobs);
        currentJobId = 0;

        lockWriting(jobQueue);
        lockReading(idleQueue);
        setReaderWriterCount(jobQueue, numThreads, 1);
        setReaderWriterCount(idleQueue, 1, numThreads);

        // Prepare idle queue.
        for (size_t i = 0; i < numJobs; ++i)
        {
            [[maybe_unused]] bool success = appendValue(idleQueue, i);
            assert(success);
        }

        // Start off threads.
        for (size_t i = 0; i < numThreads; ++i)
            pool.emplace_back(CompressionThread{this, CompressionContext<detail::bgzf_compression>{}});

        currentJobAvail = popFront(currentJobId, idleQueue);
        assert(currentJobAvail);

        CompressionJob &job = jobs[currentJobId];
        job.outputBuffer = aquireValue(serializer);
        this->setp(&job.buffer[0], &job.buffer[0] + (job.buffer.size() - 1));
    }

    ~basic_bgzf_ostreambuf()
    {
        // the buffer is now (after addFooter()) and flush will append the empty EOF marker
        flush(true);

        unlockWriting(jobQueue);

        // Wait for threads to finish there active work.
        for (auto & t : pool)
        {
            if (t.joinable())
                t.join();
        }

        unlockReading(idleQueue);
    }

    bool compressBuffer(size_t size)
    {
        // submit current job
        if (currentJobAvail)
        {
            jobs[currentJobId].size = size;
            appendValue(jobQueue, currentJobId);
        }

        // recycle existing idle job
        if (!(currentJobAvail = popFront(currentJobId, idleQueue)))
            return false;

        jobs[currentJobId].outputBuffer = aquireValue(serializer);

        return serializer;
    }

    int_type overflow(int_type c)
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if (c != static_cast<int_type>(EOF))
        {
            *this->pptr() = c;
            ++w;
        }
        if (compressBuffer(w))
        {
            CompressionJob &job = jobs[currentJobId];
            this->setp(&job.buffer[0], &job.buffer[0] + (job.buffer.size() - 1));
            return c;
        }
        else
        {
            return EOF;
        }
    }

    std::streamsize flush(bool flushEmptyBuffer = false)
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if ((w != 0 || flushEmptyBuffer) && compressBuffer(w))
        {
            CompressionJob &job = jobs[currentJobId];
            this->setp(&job.buffer[0], &job.buffer[0] + (job.buffer.size() - 1));
        }
        else
        {
            w = 0;
        }

        // wait for running compressor threads
        waitForMinSize(idleQueue, numJobs - 1);

        serializer.worker.ostream.flush();
        return w;
    }

    int sync()
    {
        if (this->pptr() != this->pbase())
        {
            int_type c = overflow(EOF);
            if (c == static_cast<int_type>(EOF))
                return -1;
        }
        return 0;
    }

    void addFooter()
    {
        // we flush the filled buffer here, so that an empty (EOF) buffer is flushed in the d'tor
        if (this->pptr() != this->pbase())
            overflow(EOF);
    }

    // returns a reference to the output stream
    ostream_reference get_ostream() const    { return serializer.worker.ostream; };
};

// --------------------------------------------------------------------------
// Class basic_bgzf_ostreambase
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_ostreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr>&                         ostream_reference;
    typedef basic_bgzf_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT> bgzf_streambuf_type;

    basic_bgzf_ostreambase(ostream_reference ostream_)
        : m_buf(ostream_)
    {
        this->init(&m_buf );
    };

    // returns the underlying zip ostream object
    bgzf_streambuf_type* rdbuf()            { return &m_buf; };
    // returns the bgzf error state
    int get_zerr() const                    { return m_buf.get_err(); };
    // returns the uncompressed data crc
    long get_crc() const                    { return m_buf.get_crc(); };
    // returns the compressed data size
    long get_out_size() const               { return m_buf.get_out_size(); };
    // returns the uncompressed data size
    long get_in_size() const                { return m_buf.get_in_size(); };

private:
    bgzf_streambuf_type m_buf;
};

// --------------------------------------------------------------------------
// Class basic_bgzf_ostream
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_ostream :
    public basic_bgzf_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT>,
    public std::basic_ostream<Elem,Tr>
{
public:
    typedef basic_bgzf_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT> bgzf_ostreambase_type;
    typedef std::basic_ostream<Elem,Tr>                        ostream_type;
    typedef ostream_type&                                      ostream_reference;

    basic_bgzf_ostream(ostream_reference ostream_) :
        bgzf_ostreambase_type(ostream_),
        ostream_type(bgzf_ostreambase_type::rdbuf())
    {}

    // flush inner buffer and zipper buffer
    basic_bgzf_ostream<Elem,Tr>& flush()
    {
        ostream_type::flush(); this->rdbuf()->flush(); return *this;
    };

    ~basic_bgzf_ostream()
    {
        this->rdbuf()->addFooter();
    }

private:
    static void put_long(ostream_reference out_, unsigned long x_);
#ifdef _WIN32
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
};

// ===========================================================================
// Typedefs
// ===========================================================================

// A typedef for basic_bgzf_ostream<char>
typedef basic_bgzf_ostream<char> bgzf_ostream;
// A typedef for basic_bgzf_ostream<wchar_t>
typedef basic_bgzf_ostream<wchar_t> bgzf_wostream;

} // namespace seqan3::contrib
