cmake_minimum_required(VERSION 3.6)
project(seqan3 CXX)

set(SEQAN3_LIBRARIES "pthread")

#set(CMAKE_CXX_FLAGS "-std=c++17")
set(CMAKE_CXX_FLAGS "-std=c++1z -fconcepts")

set(SOURCE_FILES
        alphabet/gaps/gapped_alphabet.hpp
        alphabet/gaps/gaps.hpp
        alphabet/nucleotide/dna4.hpp
        alphabet/nucleotide/dna4_sequence.hpp
        alphabet/quality/phred.hpp
        alphabet/alphabet.hpp
        alphabet/alphabet_sequence.hpp
        alphabet/compound_alphabet.hpp
        alphabet/quality.hpp
        #io/alignment/align_file.hpp
        #io/alignment/align_file_in.hpp
        #io/sequence/seq_file_in.hpp
        #io/sequence/seq_file_in_format.hpp
        #io/sequence/seq_file_in_format_fasta.hpp
        )


add_executable(seqan3 ${SOURCE_FILES})
add_subdirectory(alphabet)