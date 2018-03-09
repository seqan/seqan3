#pragma once

#include <string>
#include <variant>
#include <vector>
#include <functional>
#include <experimental/filesystem>
#include <io/detail/file_base.hpp>
#include <io/container/concepts.hpp>
#include <sequence_file_format.hpp>
#include <sequence_file_format_fasta.hpp>
#include <sequence_file_format_fastq.hpp>

namespace seqan3
{

// ==================================================================
// sequence_file_out_traits
// ==================================================================

  
template <typename t>
concept bool sequence_file_out_traits_concept = requires (t v)
{
    typename t::stream_type;
    typename t::valid_format_types;
    requires detail::meets_sequence_file_format_concept<typename t::valid_format_types> (
      std::make_index_sequence<std::variant_size_v<typename t::valid_format_types>>
      {}
    );
    
    typename t::valid_compression_formats;
};


struct sequence_file_out_default_traits
{
    using stream_type = std::ofstream;
    using valid_format_types = std::variant<sequence_file_format_fasta,
											sequence_file_format_fastq,
											sequence_file_format_raw
                                     //  	sequence_file_format_embl,
                                     //  	sequence_file_format_genbank>;
	using valid_compressions = std::variant<decltype(std::ignore)>;								 
    static constexpr std::vector<std::pair<std::string, void>> valid_compression_formats{};
};



// ==================================================================
// sequence_file_out
// ==================================================================

template <typename sequence_file_out_traits = sequence_file_out_default_traits>
    requires sequence_file_traits_concept<sequence_file_out_traits>
class sequence_file_out : protected detail::file_base<sequence_file_out_traits>
{
public:
	/* types */
	using detail::file_base<sequence_file_out_traits>::stream_type;
	using detail::file_base<sequence_file_out_traits>::valid_format_types;
    
	/* constructors */
     
	// constructor with arg
	explicit sequence_file_out(std::experimental::filesystem::path & _file_name) :
		detail::file_base<sequence_file_out_traits>(_file_name) {};  

 
	// copy construction and assignment are deleted
	// implicitly because we don't want multiple access to file
	sequence_file_out() = delete;
	sequence_file_out(sequence_file_out const &) = delete;
	sequence_file_out & operator=(sequence_file_out const &) = delete;

	// move construction and assignment are defaulted
	sequence_file_out(sequence_file_out &&) = default;
	sequence_file_out & operator=(sequence_file_out &&) = default;
  
	/* options */
	struct options_type
	{
		// post-processing filters that operate on buffer before assignment to out-value
		std::function<void(std::string &)>  sequence_filter = [] (std::string & seq)  {};
		std::function<void(std::string &)>   meta_filter = [] (std::string & meta)   {};
		std::function<void(std::string &)> qual_filter = [] (std::string & qual) {};
	};
	options_type options;
  
	/* member functions */
	//TODO make the requirements stricter
	template <typename sequence_file_out_traits>
	template <typename sequence_type, typename meta_type, typename qual_type>
        requires sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_concept<std::decay_t<meta_type>> &&
                 sequence_concept<std::decay_t<qual_type>>
    void write(sequence_type && seq, 
			   meta_type && id = std::string{}, 
			   qual_type && qual = std::string{});

	template <typename sequence_file_out_traits>
    template <typename sequence_type, typename meta_type, typename qual_type>
        requires sequence_of_sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_of_sequence_concept<std::decay_t<meta_type>> &&
                 sequence_of_sequence_concept<std::decay_t<qual_type>>
    void write(sequence_type && seqs, 
			   meta_type && ids = std::vector<std::string>{},
			   quals_type && quals = std::vector<std::string>{},
			   size_t max_records = 0);


protected:
    ~sequence_file_out() = default;
	
	/* member variables */
	using detail::file_base<sequence_file_out_traits>::format;
	using detail::file_base<sequence_file_out_traits>::stream;


};


template <typename sequence_file_out_traits>
template <typename sequence_type, typename meta_type, typename qual_type>
	requires sequence_concept<std::decay_t<sequence_type>> &&
			 sequence_concept<std::decay_t<meta_type>> &&
             sequence_concept<std::decay_t<qual_type>>
inline void sequence_file_out<sequence_file_out_traits>::write(sequence_type && seq, 
															   meta_type && id, 
															   qual_type && qual)
{
	if ( format.valueless_by_exception() )
        throw ParseError("format not specified!");
	
    std::visit([&] (sequence_file_format_concept & f) 
	{ 
		f.write(seq, id, qual, stream, options); 
		
	}, format);
	
}


template <typename sequence_file_out_traits>
template <typename sequence_type, typename meta_type, typename qual_type>
        requires sequence_of_sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_of_sequence_concept<std::decay_t<meta_type>> &&
                 sequence_of_sequence_concept<std::decay_t<qual_type>>
inline void sequence_file_out<sequence_file_out_traits>::write(sequence_type && seqs, 
															   meta_type && ids,
															   qual_type && quals,
															   size_t max_records)
{
	std::visit([&] (sequence_file_format_concept & f)
    {
        f.write(seqs, ids, quals, stream, options, max_records);
    }, format);
	
}



} // namespace seqan