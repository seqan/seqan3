// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_column_major_range_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

namespace seqan3::detail
{

/*!\brief Provides a range interface for alignment matrices.
 * \implements std::ranges::input_range
 * \ingroup alignment_matrix
 * \tparam derived_t The derived type implementing the data model of the matrix (see details for more information).
 *
 * \details
 *
 * This crtp-base class provides a range based interface for alignment matrices. It generates a range of ranges, where
 * both the outer range and the inner ranges model std::ranges::input_range. The interface uses a column-major-order
 * to iterate over the matrix, i.e. the outer range iterates over the columns and the inner range over the num_rows of
 * each column. In addition, the inner range type, also referred to as `alignment-column`, must model std::ranges::view.
 * It is a shallow wrapper around the actual column stored in the `derived_t`, who implements the corresponding data
 * model. The current alignment-column is created within the `derived_t` when the iterator over the outer range is
 * dereferenced. The iterator over the resulting `alignment-column` is further refined by the `derived_t` to
 * implement the actual behavior of the underlying alignment matrix.
 *
 * ### Requirements of the derived type
 *
 * This base class requires the following functions and types defined within the derived type in order to customise the
 * returned type over the underlying matrix.
 *
 *  * seqan3::detail::alignment_matrix_column_major_range_base::value_type
 *  * seqan3::detail::alignment_matrix_column_major_range_base::alignment_column_type
 *  * seqan3::detail::alignment_matrix_column_major_range_base::make_proxy
 *  * seqan3::detail::alignment_matrix_column_major_range_base::initialise_column
 *
 * The following functions are not mandatory but can be overriden by the derived type:
 *  * seqan3::detail::alignment_matrix_column_major_range_base::on_column_iterator_creation
 *  * seqan3::detail::alignment_matrix_column_major_range_base::before_column_iterator_increment
 *  * seqan3::detail::alignment_matrix_column_major_range_base::after_column_iterator_increment
 *
 * ### Example
 *
 * \include test/snippet/alignment/matrix/detail/alignment_matrix_column_major_range_base.cpp
 */
template <typename derived_t>
class alignment_matrix_column_major_range_base
{
private:
    //!\brief Befriend the derived type.
    friend derived_t;

    /*!\brief Represents a column within an alignment matrix.
     * \implements std::ranges::view
     * \implements std::ranges::input_range
     *
     * \details
     *
     * The alignment-column depends on the `column_data_view_type` of the derived type. This type is used to
     * alias the actual data storage of the underlying matrix implementation such that it does not own the data.
     */
    class alignment_column_type : public std::ranges::view_interface<alignment_column_type>
    {
    private:
        //!\brief A view aliasing the actual stored data column within the underlying matrix.
        using view_type = typename deferred_type<typename derived_t::column_data_view_type>::type;

        static_assert(std::ranges::random_access_range<view_type>, "Column view must support random access.");
        static_assert(std::ranges::sized_range<view_type>, "Column view must be a sized range.");
        static_assert(std::ranges::view<view_type>, "Column view must be a view.");

        //!\brief The sentinel type of the underlying view.
        using sentinel = std::ranges::sentinel_t<view_type>;

        /*!\brief The iterator over an alignment-column.
         * \implements std::forward_iterator
         *
         * \details
         *
         * The iterator represents the current cell within an alignment-column. When dereferenced it
         * calls the function `make_proxy` of the derived class to obtain a proxy value over the current matrix cell.
         * This proxy type depends on the derived type.
         */
        class iterator_type
        {
        public:

            /*!\name Associated types
             * \{
             */
            //!\brief A value type dependent on the derived type.
            using value_type = typename deferred_type<typename derived_t::value_type>::type;
            //!\brief The reference type dependent on the derived type.
            using reference = typename deferred_type<typename derived_t::reference>::type;
            //!\brief Pointer type.
            using pointer = void;
            //!\brief Difference type.
            using difference_type = std::ranges::range_difference_t<view_type>;
            //!\brief Iterator category.
            using iterator_category = std::forward_iterator_tag;
            //!\}

            /*!\name Constructors, destructor and assignment
             * \{
             */
            constexpr iterator_type() = default; //!< Defaulted.
            constexpr iterator_type(iterator_type const &) = default; //!< Defaulted.
            constexpr iterator_type(iterator_type &&) = default; //!< Defaulted.
            constexpr iterator_type & operator=(iterator_type const &) = default; //!< Defaulted.
            constexpr iterator_type & operator=(iterator_type &&) = default; //!< Defaulted.
            ~iterator_type() = default; //!< Defaulted.

            /*!\brief Construction from the underlying alignment-column.
             * \param host The alignment-column for this iterator.
             */
            explicit constexpr iterator_type(alignment_column_type & host) :
                host_ptr{&host},
                host_iter{host.ref.begin()}
            {
                host_ptr->me_ptr->on_column_iterator_creation(host_iter);
            }
            //!\}

            /*!\name Element access
             * \{
             */
            //!\brief Returns a proxy for the current alignment cell.
            constexpr reference operator*() const noexcept
            {
                static_assert(std::convertible_to<decltype(host_ptr->me_ptr->make_proxy(host_iter)), reference>,
                              "The returned type of make_proxy must be convertible to the reference type.");
                assert(host_ptr != nullptr);
                return host_ptr->me_ptr->make_proxy(host_iter);
            }
            //!\}

            /*!\name Arithmetic operators
             * \{
             */
            //!\brief Advances the iterator by one.
            constexpr iterator_type & operator++() noexcept
            {
                assert(host_ptr != nullptr);
                assert(host_ptr->me_ptr != nullptr);

                host_ptr->me_ptr->before_column_iterator_increment(host_iter);
                ++host_iter;
                host_ptr->me_ptr->after_column_iterator_increment(host_iter);
                return *this;
            }

            //!\brief Advances the iterator and returns previous iterator.
            constexpr iterator_type operator++(int) noexcept
            {
                iterator_type tmp{*this};
                ++(*this);
                return tmp;
            }
            //!\}

            /*!\name Comparison operators
             * \{
             */
            //!\brief Returns `true` if the host iterator reached the end, `false` otherwise.
            constexpr bool operator==(sentinel const & rhs) const noexcept
            {
                return host_iter == rhs;
            }

            //!\copydoc operator==
            friend constexpr bool operator==(sentinel const & lhs, iterator_type const & rhs) noexcept
            {
                return rhs == lhs;
            }

            //!\brief Returns `true` if both iterators are equal, `false` otherwise.
            constexpr bool operator==(iterator_type const & rhs) const noexcept
            {
                return (host_ptr == rhs.host_ptr) && (host_iter == rhs.host_iter);
            }

            //!\brief Returns `true` if the host iterator did not reach the end, `false` otherwise.
            constexpr bool operator!=(sentinel const & rhs) const noexcept
            {
                return !(*this == rhs);
            }

            //!\copydoc operator!=
            friend constexpr bool operator!=(sentinel const & lhs, iterator_type const & rhs) noexcept
            {
                return rhs != lhs;
            }

            //!\brief Returns `true` if both iterators are not equal, `false` otherwise.
            constexpr bool operator!=(iterator_type const & rhs) const noexcept
            {
                return !(*this == rhs);
            }
            //!\}

        private:
            //!\brief Pointer to the underlying alignment-column.
            alignment_column_type * host_ptr{nullptr};
            //!\brief Wrapped iterator over aliased matrix column.
            std::ranges::iterator_t<view_type> host_iter{};
        }; // class iterator_type

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr alignment_column_type() = default; //!< Defaulted.
        constexpr alignment_column_type(alignment_column_type const &) = default; //!< Defaulted.
        constexpr alignment_column_type(alignment_column_type &&) = default; //!< Defaulted.
        constexpr alignment_column_type & operator=(alignment_column_type const &) = default; //!< Defaulted.
        constexpr alignment_column_type & operator=(alignment_column_type &&) = default; //!< Defaulted.
        ~alignment_column_type() = default; //!< Defaulted.

        /*!\brief Constructs from the derived type.
         * \param[in] me  An instance of the derived type.
         * \param[in] ref The underlying view referencing the current column in the matrix.
         *
         * \details
         *
         * This constructor is called by the derived type when invoking the function `initialise_column`.
         */
        constexpr alignment_column_type(derived_t & me, view_type ref) :
            ref{std::move(ref)},
            me_ptr{&me}
        {}
        //!\}

        /*!\name Iterators
         * \brief This column type is not const-iterable.
         * \{
         */
        //!\brief Returns an iterator to the begin of the column.
        constexpr auto begin() noexcept
        {
            assert(me_ptr != nullptr);
            return iterator_type{*this};
        }

        //!\brief Deleted begin for const-qualified alignment-columns.
        constexpr auto begin() const noexcept = delete; // Not needed by the alignment algorithm

        //!\brief Deleted begin for const-qualified alignment-columns.
        constexpr auto cbegin() const noexcept = delete;  // Not needed by the alignment algorithm

        //!\brief Returns an iterator to the end of the column.
        constexpr sentinel end() noexcept
        {
            return ref.end();
        }

        //!\brief Deleted end for const-qualified alignment-columns.
        constexpr sentinel end() const noexcept = delete; // Not needed by the alignment algorithm

        //!\brief Deleted end for const-qualified alignment-columns.
        constexpr sentinel cend() const noexcept = delete; // Not needed by the alignment algorithm
        //!\}

        //!\brief Returns the size the alignment column.
        constexpr size_t size() const noexcept
        {
            return std::ranges::size(ref);
        }

    private:
        //!\brief The aliased alignment-column.
        view_type ref{};
        //!\brief Pointer to the derived type.
        derived_t * me_ptr{};
    }; // class alignment_column_type

    /*!\brief A column iterator over the alignment matrix.
     * \implements std::input_iterator
     *
     * \details
     *
     * This iterator enables the iteration over the underlying alignment matrix in column-major-order.
     * When dereferenced the iterator returns a seqan3::detail::alignment_matrix_column_major_range_base::column_type.
     * The initialisation of the alignment-column depends on the derived type.
     * This iterator models the std::input_iterator since in some cases the underlying semantics would only guarantee
     * an std::input_iterator, e.g. in the one column score matrix implementation. Both the previous and the current
     * column refer to the same data storage.
     */
    class iterator_type
    {
    public:
        /*!\name Associated types
         * \{
         */
        //!\brief The alignment-column type.
        using value_type = alignment_column_type;
        //!\brief A reference to the alignment-column type.
        using reference = value_type;
        //!\brief Pointer type.
        using pointer = void;
        //!\brief Difference type.
        using difference_type = std::ranges::range_difference_t<alignment_column_type>;
        //!\brief Iterator category.
        using iterator_category = std::input_iterator_tag;
        //!\}

        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr iterator_type() = default; //!< Defaulted.
        constexpr iterator_type(iterator_type const &) = default; //!< Defaulted.
        constexpr iterator_type(iterator_type &&) = default; //!< Defaulted.
        constexpr iterator_type & operator=(iterator_type const &) = default; //!< Defaulted.
        constexpr iterator_type & operator=(iterator_type &&) = default; //!< Defaulted.
        ~iterator_type() = default; //!< Defaulted.

        /*!\brief Construction from an instance of the derived type.
         * \param me     An instance of the derived type.
         */
        explicit constexpr iterator_type(derived_t & me) :
            me_ptr{&me},
            column_index{0}
        {}
        //!\}

        /*!\name Element access
         * \{
         */
        //!\brief Returns the current alignment-column.
        constexpr reference operator*() const noexcept
        {
            static_assert(std::convertible_to<decltype(me_ptr->initialise_column(column_index)), reference>,
                          "The returned type of initialise_column must be convertible to the reference type.");
            return me_ptr->initialise_column(column_index);
        }
        //!\}

        /*!\name Arithmetic operators
         * \{
         */
        //!\brief Increments by one.
        constexpr iterator_type & operator++() noexcept
        {
            ++column_index;
            return *this;
        }

        //!\brief Increments by one.
        constexpr void operator++(int) noexcept
        {
            ++(*this);
        }
        //!\}

        /*!\name Comparison operators
         * \{
         */
        //!\brief Returns `true` if the behind-the-end column was reached, `false` otherwise.
        constexpr bool operator==(std::ranges::default_sentinel_t const &) const noexcept
        {
            return column_index == me_ptr->num_cols;
        }

        //!\copydoc operator==
        friend constexpr bool operator==(std::ranges::default_sentinel_t const & lhs, iterator_type const & rhs)
            noexcept
        {
            return rhs == lhs;
        }

        //!\brief Returns `true` if the behind-the-end column was not reached, `false` otherwise.
        constexpr bool operator!=(std::ranges::default_sentinel_t const & rhs) const noexcept
        {
            return !(*this == rhs);
        }

        //!\copydoc operator!=
        friend constexpr bool operator!=(std::ranges::default_sentinel_t const & lhs, iterator_type const & rhs)
            noexcept
        {
            return rhs != lhs;
        }
        //!\}

    private:
        //!\brief Pointer to the derived type.
        derived_t * me_ptr;
        //!\brief The current column index.
        size_t column_index{};
    }; // class iterator_type

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    constexpr alignment_matrix_column_major_range_base() = default;
    //!\brief Defaulted.
    constexpr alignment_matrix_column_major_range_base(alignment_matrix_column_major_range_base const &) = default;
    //!\brief Defaulted.
    constexpr alignment_matrix_column_major_range_base(alignment_matrix_column_major_range_base &&) = default;
    //!\brief Defaulted.
    constexpr alignment_matrix_column_major_range_base &
    operator=(alignment_matrix_column_major_range_base const &) = default;
    //!\brief Defaulted.
    constexpr alignment_matrix_column_major_range_base &
    operator=(alignment_matrix_column_major_range_base &&) = default;
    //!\brief Defaulted.
    ~alignment_matrix_column_major_range_base() = default;
    //!\}

    /*!\name Required types
     * \brief The derived class must implement the following types.
     * \{
     */
    /*!\brief The proxy type of an alignment matrix.
     * \note This type must be defined by the derived type.
     */
    SEQAN3_DOXYGEN_ONLY(typedef /*IMPLEMENTATION_DEFINED*/ value_type;)

    /*!\brief The view over the current alignment-column; must model std::ranges::view and std::ranges::input_range.
     * \note This type must be defined by the derived type.
     */
    SEQAN3_DOXYGEN_ONLY(typedef /*IMPLEMENTATION_DEFINED*/ column_data_view_type;)
    //!\}

    /*!\name Required interfaces
     * \brief The derived class must implement the following methods.
     * \{
     */
    /*!\brief Creates the proxy value returned when dereferencing the alignment-column-iterator.
     * \param[in] host_iter The wrapped iterator to the actual memory storage.
     * \returns A proxy for the current matrix position, which must be of type
     *          seqan3::detail::alignment_matrix_column_major_range_base::value_type.
     */
    SEQAN3_DOXYGEN_ONLY(value_type make_proxy(iter_t host_iter) noexcept {})

    /*!\brief Returns the current alignment-column at the given `column_index`.
     * \param[in] column_index The current column position of the outer matrix iterator.
     * \returns A view representing the current alignment-column
     *          (seqan3::detail::alignment_matrix_column_major_range_base::alignment_column_type)
     *
     * \details
     *
     * Creates a new `seqan3::detail::alignment_matrix_column_major_range_base::alignment_column_type` initialised with
     * the current column of the alignment matrix depending on the given `column_index`. The `alignment_column_type`
     * stores a column view as defined by the derived class using the
     * seqan3::detail::alignment_matrix_column_major_range_base::column_data_view_type type definition.
     */
    SEQAN3_DOXYGEN_ONLY(alignment_column_type initialise_column(size_t column_index) {})

    /*!\brief Allows additional initialisations when calling begin on an alignment-column.
     * \tparam iter_t The iterator type of the host iterator.
     * \param[in] host_iter The wrapped iterator to the actual memory storage.
     */
    template <typename iter_t>
    constexpr void on_column_iterator_creation(iter_t SEQAN3_DOXYGEN_ONLY(host_iter)) noexcept
    {}

    /*!\brief Allows to perform additional steps before incrementing the alignment-column-iterator.
     * \tparam iter_t The iterator type of the host iterator.
     * \param[in] host_iter The wrapped iterator to the actual memory storage.
     */
    template <typename iter_t>
    constexpr void before_column_iterator_increment(iter_t SEQAN3_DOXYGEN_ONLY(host_iter)) noexcept
    {}

    /*!\brief Allows to perform additional steps after incrementing the alignment-column-iterator.
     * \tparam iter_t The iterator type of the host iterator.
     * \param[in] host_iter The wrapped iterator to the actual memory storage.
     */
    template <typename iter_t>
    constexpr void after_column_iterator_increment(iter_t SEQAN3_DOXYGEN_ONLY(host_iter)) noexcept
    {}
    //!\}

    //!\brief The type of the iterator.
    using iterator = iterator_type;
    //!\brief The type of sentinel.
    using sentinel = std::ranges::default_sentinel_t;

public:
    /*!\name Iterators
     * \brief The matrix type is not const-iterable.
     * \{
     */
    //!\brief Returns an iterator to the first column of the matrix.
    constexpr iterator begin() noexcept
    {
        return iterator{static_cast<derived_t &>(*this)};
    }

    //!\brief Deleted begin for const-qualified alignment matrix.
    constexpr iterator begin() const noexcept = delete; // not needed for the alignment algorithm

    //!\brief Deleted begin for const-qualified alignment matrix.
    constexpr iterator cbegin() const noexcept = delete; // not needed for the alignment algorithm

    //!\brief Returns a sentinel marking the end of the matrix.
    constexpr sentinel end() noexcept
    {
        return std::ranges::default_sentinel;
    }

    //!\brief Deleted end for const-qualified alignment matrix.
    constexpr sentinel end() const noexcept = delete; // not needed for the alignment algorithm

    //!\brief Deleted end for const-qualified alignment matrix.
    constexpr sentinel cend() const noexcept = delete; // not needed for the alignment algorithm
    //!\}
};
} // namespace seqan3::detail
