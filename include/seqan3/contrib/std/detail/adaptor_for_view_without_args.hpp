// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::stl::detail::adaptor_for_view_without_args
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_DETAIL_ADAPTOR_FOR_VIEW_WITHOUT_ARGS
#define SEQAN_STD_DETAIL_ADAPTOR_FOR_VIEW_WITHOUT_ARGS

#include "adaptor_base.hpp"

namespace seqan::stl::detail
{

// ============================================================================
//  adaptor_for_view_without_args
// ============================================================================

//!\brief Template for range adaptor closure objects that store no arguments and delegate to the view constructor.
template <template <typename, typename...> typename view_type>
class adaptor_for_view_without_args : public adaptor_base<adaptor_for_view_without_args<view_type>>
{
private:
    //!\brief Type of the CRTP-base.
    using base_type = adaptor_base<adaptor_for_view_without_args<view_type>>;

    //!\brief Befriend the base class so it can call impl().
    friend base_type;

    /*!\brief Call the view's constructor with the given arguments (all of the base class'es operators ultimately
     * resolve to this function call).
     * \tparam    arg_types The arguments to the view (this first one will be a range, the rest is optional).
     * \param[in] args      The arguments to the constructor.
     * \returns An instance of `view_type`.
     */
    template <typename... arg_types>
    static auto impl(arg_types &&... args)
    {
        return view_type{std::forward<arg_types>(args)...};
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args() = default;
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args(adaptor_for_view_without_args const &) noexcept = default;
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args(adaptor_for_view_without_args &&) noexcept = default;
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args & operator=(adaptor_for_view_without_args const &) noexcept = default;
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args & operator=(adaptor_for_view_without_args &&) noexcept = default;
    //!\brief Defaulted.
    ~adaptor_for_view_without_args() noexcept = default;

    //!\brief Inherit the base type's constructors.
    using base_type::base_type;
    //!\}
};

} // namespace seqan::stl::detail

#endif // SEQAN_STD_DETAIL_ADAPTOR_FOR_VIEW_WITHOUT_ARGS
