// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

//![iterator]
#include <iostream>
#include <ranges>
#include <vector>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/range/detail/inherited_iterator_base.hpp>

using seqan3::operator""_dna5;

/* The iterator template */
template <std::ranges::forward_range urng_t> // CRTP derivation ↓
class my_iterator : public seqan3::detail::inherited_iterator_base<my_iterator<urng_t>, std::ranges::iterator_t<urng_t>>
{
private:
    static_assert(seqan3::nucleotide_alphabet<std::ranges::range_reference_t<urng_t>>,
                  "You can only iterate over ranges of nucleotides!");

    // the immediate base type is the CRTP-layer
    using base_t = seqan3::detail::inherited_iterator_base<my_iterator<urng_t>, std::ranges::iterator_t<urng_t>>;

public:
    // the member types are never imported automatically, but can be explicitly inherited:
    using typename base_t::difference_type;
    using typename base_t::iterator_category;
    using typename base_t::value_type;
    // this member type is overwritten as we do above:
    using reference = value_type;
    // Explicitly set the pointer to void as we return a temporary.
    using pointer = void;

    // define rule-of-six:
    my_iterator() = default;
    my_iterator(my_iterator const &) = default;
    my_iterator(my_iterator &&) = default;
    my_iterator & operator=(my_iterator const &) = default;
    my_iterator & operator=(my_iterator &&) = default;
    ~my_iterator() = default;
    // and a constructor that takes the base_type:
    my_iterator(base_t it) : base_t{std::move(it)}
    {}

    // we don't need to implement the ++ operators anymore!

    // only overload the operators that you actually wish to change:
    reference operator*() const noexcept
    {
        return seqan3::complement(base_t::operator*());
    }

    // Since the reference type changed we might as well need to override the subscript-operator.
    reference operator[](difference_type const n) const noexcept
        requires std::random_access_iterator<std::ranges::iterator_t<urng_t>>
    {
        return seqan3::complement(base_t::operator[](n));
    }

    // We delete arrow operator because of the temporary. An alternative could be to return the temporary
    // wrapped in a std::unique_ptr.
    pointer operator->() const noexcept = delete;
};

// The inherited_iterator_base creates the necessary code so we also model RandomAccess now!
static_assert(std::random_access_iterator<my_iterator<std::vector<seqan3::dna5>>>);
//![iterator]

//![view_header]
/* The view class template */
template <std::ranges::view urng_t> // CRTP derivation ↓
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
    template <std::ranges::viewable_range orng_t>
    my_view(orng_t && urange_) : urange{std::views::all(std::forward<orng_t>(urange_))}
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
template <std::ranges::viewable_range orng_t>
my_view(orng_t &&) -> my_view<std::views::all_t<orng_t>>;
//![view_deduction_guide]

//![adaptor_type_definition]
/* The adaptor object's type definition */
struct my_view_fn
{
    template <std::ranges::input_range urng_t>
    auto operator()(urng_t && urange) const
    {
        return my_view{std::forward<urng_t>(urange)};
    }

    template <std::ranges::input_range urng_t>
    friend auto operator|(urng_t && urange, my_view_fn const &)
    {
        return my_view{std::forward<urng_t>(urange)};
    }
};
//![adaptor_type_definition]

//![adaptor_object_definition]
/* The adaptor object's definition */
namespace views
{

inline constexpr my_view_fn my{};

}
//![adaptor_object_definition]

//![main_it]
int main()
{
    std::vector<seqan3::dna5> vec{"GATTACA"_dna5};

    /* try the iterator */
    using my_it_concrete = my_iterator<std::vector<seqan3::dna5>>;

    my_it_concrete it{vec.begin()};

    // now you can use operator[] on the iterator
    for (size_t i = 0; i < 7; ++i)
        std::cout << seqan3::to_char(it[i]) << ' ';
    std::cout << '\n';
    //![main_it]

    //![main_range]
    /* try the range */
    my_view v{vec};
    static_assert(std::ranges::random_access_range<decltype(v)>);
    seqan3::debug_stream << '\n' << v << '\n';
    //![main_range]

    //![main_adaptor]
    /* try the adaptor */
    auto v2 = vec | std::views::reverse | ::views::my;
    static_assert(std::ranges::random_access_range<decltype(v2)>);
    seqan3::debug_stream << v2 << '\n';
    //![main_adaptor]

    //![end]
}
//![end]
