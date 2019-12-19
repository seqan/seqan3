//![start]
#include <iostream>
#include <vector>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/std/ranges>

using seqan3::operator""_dna5;

template <std::ranges::forward_range urng_t> // the underlying range type
struct my_iterator : std::ranges::iterator_t<urng_t>
{
//![start]
    //![static_assert]
    static_assert(seqan3::nucleotide_alphabet<std::ranges::range_reference_t<urng_t>>,
                  "You can only iterate over ranges of nucleotides!");
    //![static_assert]

    //![solution1a]
    using base_t = std::ranges::iterator_t<urng_t>;

    // these member types are just exposed from the base type
    using value_type            = typename std::iterator_traits<base_t>::value_type;
    using pointer               = typename std::iterator_traits<base_t>::pointer;
    //![solution1a]
    //![reference]
    // If the value_type is seqan3::dna5, the reference type of the vector is
    // seqan3::dna5 & and operator* returns this so you can change the values
    // in the vector through it's iterator.
    // This won't work anymore, because now we are creating new values on access
    // so we now need to change this type to reflect that:
    using reference             = value_type;
    //![reference]

    //![solution1b]
    // this member type is explicitly set to forward_iterator_tag because we are not
    // implementing the remaining requirements
    using iterator_category     = std::forward_iterator_tag;

    // the following operators need to be explicitly defined, because the inherited
    // version has wrong return types (base_t instead of my_iterator)
    my_iterator & operator++()
    {
        base_t::operator++();       // call the implementation of the base type
        return *this;
    }

    my_iterator operator++(int)
    {
        my_iterator cpy{*this};
        ++(*this);
        return cpy;
    }

    // we do not need to define constructors, because {}-initialising this with one argument
    // calls the base-class's constructor which is what we want
    // NOTE: depending on your view/iterator, you might need/want to define your own

    // we do not need to define comparison operators, because there are comparison operators
    // for the base class and our type is implicitly convertible to the base type
    // NOTE: in our case comparing the converted-to-base iterator with end behaves as desired,
    // but it strongly depends on your view/iterator whether the same holds
    //![solution1b]

    //![dereference]
    // The operator* facilitates the access to the element.
    // This implementation calls the base class's implementation but passes the return value
    // through the seqan3::complement() function before returning it.
    reference operator*() const
    {
        return seqan3::complement(base_t::operator*());
    }
    //![dereference]

//![end]
};

// verify that your type models the concept
static_assert(std::forward_iterator<my_iterator<std::vector<seqan3::dna5>>>);

int main()
{
    std::vector<seqan3::dna5> vec{"GATTACA"_dna5};

    // instantiate the template over the underlying vector's iterator and sentinel
    // (for all standard containers the sentinel type is the same as the iterator type)
    using my_it_concrete = my_iterator<std::vector<seqan3::dna5>>;

    // create an iterator that is constructed with vec.begin() as the underlying iterator
    my_it_concrete it{vec.begin()};

    // iterate over vec, but with your custom iterator
    while (it != vec.end())
        std::cout << seqan3::to_char(*it++) << ' ';
}
//![end]
