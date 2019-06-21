#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/concept.hpp>

class ab : public seqan3::alphabet_base<ab, 2>
{
private:
    // map 0 -> A and 1 -> B
    static std::array<char_type, alphabet_size> constexpr rank_to_char{'A', 'B'};

    // map every letter to rank zero, except Bs
    static std::array<rank_type, 256> constexpr char_to_rank
    {
        // initialise with an immediately evaluated lambda expression:
        []()
        {
            std::array<rank_type, 256> ret{}; // initialise all values with 0 / 'A'

            // only 'b' and 'B' result in rank 1
            ret['b'] = 1;
            ret['B'] = 1;

            return ret;
        }()
    };

    // make the base class a friend so it can access the tables:
    friend alphabet_base<ab, 2>;
};

static_assert(seqan3::Alphabet<ab>);
