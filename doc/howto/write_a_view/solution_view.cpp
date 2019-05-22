//![iterator]
#include <iostream>
#include <vector>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/detail/inherited_iterator_base.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

/* The iterator template */
template <std::ranges::ForwardRange urng_t>            // CRTP derivation ↓
class my_iterator : public seqan3::detail::inherited_iterator_base<my_iterator<urng_t>,
                                                                   std::ranges::iterator_t<urng_t>>
{
private:
    static_assert(NucleotideAlphabet<reference_t<urng_t>>,
                  "You can only iterate over ranges of nucleotides!");

    // the immediate base type is the CRTP-layer
    using base_t = seqan3::detail::inherited_iterator_base<my_iterator<urng_t>,
                                                           std::ranges::iterator_t<urng_t>>;

public:
    // the member types are never imported automatically, but can be explicitly inherited:
    using typename base_t::value_type;
    using typename base_t::pointer;
    using typename base_t::iterator_category;
    // this member type is overwritten as we do above:
    using reference = value_type;

    // define rule-of-six:
    my_iterator() = default;
    my_iterator(my_iterator const &) = default;
    my_iterator(my_iterator &&) = default;
    my_iterator & operator=(my_iterator const &) = default;
    my_iterator & operator=(my_iterator &&) = default;
    ~my_iterator() = default;
    // and a constructor that takes the base_type:
    my_iterator(base_t it) : base_t{std::move(it)} {}

    // we don't need to implement the ++ operators anymore!

    // only overload the operators that you actually wish to change:
    reference operator*() const noexcept
    {
        return complement(base_t::operator*());
    }
};

// The inherited_iterator_base creates the necessary code so we also model RandomAccess now!
static_assert(std::RandomAccessIterator<my_iterator<std::vector<dna5>>>);
//![iterator]

//![view_header]
/* The view class template */
template <std::ranges::View urng_t>  // CRTP derivation ↓
class my_view : public std::ranges::view_interface<my_view<urng_t>>
{
//![view_header]
//![view_private]
private:
    // this is the underlying range
    urng_t urange;
//![view_private]

//![view_member_types]
public:
    // Types of the iterators
    using iterator = my_iterator<urng_t>;
    using const_iterator = my_iterator<urng_t const>;
    //![view_member_types]

    //![view_constructors]
    // construct from a view
    my_view(urng_t urange_) : urange{std::move(urange_)}
    {}

    // construct from non-view that can be view-wrapped
    template <std::ranges::ViewableRange orng_t>
    my_view(orng_t && urange_) : urange{std::view::all(std::forward<orng_t>(urange_))}
    {}
    //![view_constructors]

    //![view_begin]
    auto begin() noexcept
    {
        return iterator{std::ranges::begin(urange)};
    }

    auto begin() const noexcept
    {
        return const_iterator{std::ranges::begin(urange)};
    }

    auto cbegin() const noexcept
    {
        return const_iterator{std::ranges::begin(urange)};
    }
    //![view_begin]

    //![view_end]
    auto end() noexcept
    {
        return std::ranges::end(urange);
    }

    auto end() const noexcept
    {
        return std::ranges::end(urange);
    }

    auto cend() const noexcept
    {
        return std::ranges::end(urange);
    }
    //![view_end]
};

//![view_deduction_guide]
// A deduction guide for the view class template
template <std::ranges::ViewableRange orng_t>
my_view(orng_t &&) -> my_view<std::ranges::all_view<orng_t>>;
//![view_deduction_guide]

//![adaptor_type_definition]
/* The adaptor object's type definition */
struct my_view_fn
{
    template <std::ranges::InputRange urng_t>
    auto operator()(urng_t && urange) const
    {
        return my_view{std::forward<urng_t>(urange)};
    }

    template <std::ranges::InputRange urng_t>
    friend auto operator|(urng_t && urange, my_view_fn const &)
    {
        return my_view{std::forward<urng_t>(urange)};
    }
};
//![adaptor_type_definition]

//![adaptor_object_definition]
/* The adaptor object's definition */
namespace view
{

inline constexpr my_view_fn my{};

}
//![adaptor_object_definition]

//![main_it]
int main()
{
    std::vector<dna5> vec{"GATTACA"_dna5};

    /* try the iterator */
    using my_it_concrete = my_iterator<std::vector<dna5>>;

     my_it_concrete it{vec.begin()};

    // now you can use operator[] on the iterator
    for (size_t i = 0; i < 7; ++i)
        std::cout << to_char(it[i]) << ' ';
//![main_it]

//![main_range]
    /* try the range */
    my_view v{vec};
    static_assert(std::ranges::RandomAccessRange<decltype(v)>);
    debug_stream << '\n' << v << '\n';
//![main_range]

//![main_adaptor]
    /* try the adaptor */
    auto v2 = vec | std::view::reverse | ::view::my;
    static_assert(std::ranges::RandomAccessRange<decltype(v2)>);
    debug_stream << v2 << '\n';
//![main_adaptor]

//![end]
}
//![end]
