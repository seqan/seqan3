#include <vector>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>

using namespace seqan3;

int main()
{

// TODO(comments): outputs [A, A, G] instead of AAG

{
debug_stream << "Example extend_left(seq)\n";
//! [extend_left_seq]
std::vector<dna4> genome{"GAATTAATGAAC"_dna4};
bi_fm_index index{genome}; // build the index

auto it = index.begin();         // create an iterator
it.extend_right("AAC"_dna4);     // search the sequence "AAC"
debug_stream << it.query() << '\n'; // outputs "AAC"
it.extend_left("ATG"_dna4);      // extend the query to "ATGAAT"
                                 // The rightmost character of "ATG" is extended to the left first.
debug_stream << it.query() << '\n'; // outputs "ATGAAT"
//! [extend_left_seq]
}

{
debug_stream << "Example cycle_back() and cycle_front()\n";
//! [cycle]
std::vector<dna4> genome{"GAATTAATGAAC"_dna4};
bi_fm_index index{genome};           // build the index

auto it = index.begin();             // create an iterator
// it.cycle_back();                  // cycle_back / cycle_front on begin() is undefined behaviour!
it.extend_right("AAC"_dna4);         // search the sequence "AAC"
debug_stream << it.query() << '\n';     // outputs "AAC"
debug_stream << it.last_char() << '\n'; // outputs 'C'

// it.cycle_front();                 // undefined behaviour! only cycle_back() is allowed after extend_right()
it.cycle_back();                     // search the sequence "AAT"
debug_stream << it.query() << '\n';     // outputs "AAT"
debug_stream << it.last_char() << '\n'; // outputs 'T'

it.extend_left('G'_dna4);             // search the sequence "GAAT"
debug_stream << it.query() << '\n';     // outputs "GAAC"
debug_stream << it.last_char() << '\n'; // outputs 'G'

// it.cycle_back();                  // undefined behaviour! only cycle_front() is allowed after extend_left()
it.cycle_front();                    // search the sequence "TAAT"
debug_stream << it.query() << '\n';     // outputs "TAAT"
debug_stream << it.last_char() << '\n'; // outputs 'T'

it.cycle_front();                    // search the sequence "TAAT"
debug_stream << it.query() << '\n';     // outputs "TAAT"
debug_stream << it.last_char() << '\n'; // outputs 'T'
//! [cycle]
}

{
debug_stream << "Example to_fwd_iterator()\n";
//! [to_fwd_iterator]
std::vector<dna4> genome{"GAATTAACGAAC"_dna4};
bi_fm_index index{genome};               // build the index

auto it = index.begin();                 // create an iterator
it.extend_left("AAC"_dna4);              // search the sequence "AAC"
debug_stream << it.query() << '\n';         // outputs "AAC"
auto uni_it = it.to_fwd_iterator();      // unidirectional iterator on the text "GAATTAACGAAC"
debug_stream << uni_it.query() << '\n';     // outputs "CAA"
// Undefined behaviour! Cannot be called on the forward iterator if the last extension on the bidirectional
// iterator was to the left:
// it.cycle_back();
// debug_stream << it.last_char() << '\n';

uni_it.extend_right('G'_dna4);            // search the sequence "AACG"
debug_stream << uni_it.query() << '\n';     // outputs "AACG"
debug_stream << uni_it.last_char() << '\n'; // outputs 'G'
uni_it.cycle_back();                     // returns false since there is no sequence "AACT" in the text.
//! [to_fwd_iterator]
}

{
debug_stream << "Example to_rev_iterator()\n";
//! [to_rev_iterator]
std::vector<dna4> genome{"GAATTAACGAAC"_dna4};
bi_fm_index index{genome};               // build the index

auto it = index.begin();                 // create an iterator
it.extend_right("AAC"_dna4);             // search the sequence "AAC"
debug_stream << it.query() << '\n';         // outputs "AAC"
auto uni_it = it.to_rev_iterator();      // unidirectional iterator on the text "CAAGCAATTAAG"
debug_stream << uni_it.query() << '\n';     // outputs "CAA"
// Undefined behaviour! Cannot be called on the reversed iterator if the last extension on the bidirectional
// iterator was to the right:
// it.cycle_back();
// debug_stream << it.last_char() << '\n';

uni_it.extend_right('G'_dna4);            // search the sequence "CAAG"
debug_stream << uni_it.query() << '\n';     // outputs "CAAG"
debug_stream << uni_it.last_char() << '\n'; // outputs 'G'
uni_it.cycle_back();                     // search the sequence "CAAT"
//! [to_rev_iterator]
}

    return 0;
}
