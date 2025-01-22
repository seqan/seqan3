// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::adaptor_for_view_without_args
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/range/detail/adaptor_base.hpp>

namespace seqan3::detail
{

// ============================================================================
//  adaptor_for_view_without_args
// ============================================================================

/*!\brief Template for range adaptor closure objects that store no arguments and always delegate to the
 *        view constructor.
 * \tparam view_type The view type template.
 * \ingroup core_range
 *
 * \details
 *
 * Use this adaptor template when you always want to delegate to the view's constructor and you have no arguments.
 * Since it's a one-line it's easier than specialising seqan3::detail::adaptor_base.
 *
 * # Example
 *
 * (from include/seqan3/utility/views/single_pass_input.hpp)
 *
 * This is the signature of the view type template in namespace seqan3::detail:
 *
 * \snippet include/seqan3/utility/views/single_pass_input.hpp view_def
 *
 * This is the definition of the range adaptor closure object, it will always delegate to the
 * constructor: seqan3::detail::single_pass_input_view::single_pass_input_view():
 *
 * \snippet include/seqan3/utility/views/single_pass_input.hpp adaptor_def
 */
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

} // namespace seqan3::detail
