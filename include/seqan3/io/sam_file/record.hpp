// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::sam_record.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/record.hpp>

namespace seqan3
{
/*!\brief The record type of seqan3::alignment_file_input.
 * \ingroup alignment_file
 * \implements seqan3::detail::record_like
 * \tparam field_types The types of the fields in this record as a seqan3::type_list.
 * \tparam field_ids A seqan3::fields type with seqan3::field IDs corresponding to field_types.
 */
template <typename field_types, typename field_ids>
class sam_record : public record<field_types, field_ids>
{
    //!\brief The base class.
    using base_t = record<field_types, field_ids>;

    //!\brief The underlying std::tuple class.
    using tuple_base_t = typename base_t::base_type;

    //!\copydoc seqan3::record::field_constant
    template <field f>
    using field_constant = typename base_t::template field_constant<f>;

    using base_t::get_impl;
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    sam_record() = default; //!< Defaulted.
    sam_record(sam_record const &) = default; //!< Defaulted.
    sam_record & operator=(sam_record const &) = default; //!< Defaulted.
    sam_record(sam_record &&) = default; //!< Defaulted.
    sam_record & operator=(sam_record &&) = default; //!< Defaulted.
    ~sam_record() = default; //!< Defaulted.

    //!\brief Inherit std::tuple's and seqan3::record constructors.
    using base_t::base_t;
    //!\}

    //!\brief The identifier, usually a string. (SAM Column ID: QNAME)
    decltype(auto) id() &&
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::id
    decltype(auto) id() const &&
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::id
    decltype(auto) id() &
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::id
    decltype(auto) id() const &
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief The "sequence", usually a range of nucleotides or amino acids. (SAM Column ID: SEQ)
    decltype(auto) sequence() &&
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::sequence
    decltype(auto) sequence() const &&
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::sequence
    decltype(auto) sequence() &
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::sequence
    decltype(auto) sequence() const &
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief The qualities, usually in Phred score notation. (SAM Column ID: QUAL)
    decltype(auto) base_qualities() &&
    {
        return get_impl(field_constant<seqan3::field::qual>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::base_qualities
    decltype(auto) base_qualities() const &&
    {
        return get_impl(field_constant<seqan3::field::qual>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::base_qualities
    decltype(auto) base_qualities() &
    {
        return get_impl(field_constant<seqan3::field::qual>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::base_qualities
    decltype(auto) base_qualities() const &
    {
        return get_impl(field_constant<seqan3::field::qual>{}, static_cast<tuple_base_t const &>(*this));
    }

    /*!\brief Sequence (seqan3::sam_record::sequence) relative start position (0-based), unsigned value.
     *        (SAM Column ID: POS)
     *
     * \details
     *
     * The position is the length of the soft-clipping at the start of the seqan3::sam_record::cigar_sequence if a
     * soft-clipping is present and `0` otherwise.
     */
    decltype(auto) sequence_position() &&
    {
        return get_impl(field_constant<seqan3::field::offset>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::sequence_position
    decltype(auto) sequence_position() const &&
    {
        return get_impl(field_constant<seqan3::field::offset>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::sequence_position
    decltype(auto) sequence_position() &
    {
        return get_impl(field_constant<seqan3::field::offset>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::sequence_position
    decltype(auto) sequence_position() const &
    {
        return get_impl(field_constant<seqan3::field::offset>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief The (pairwise) alignment stored in an object that models seqan3::detail::pairwise_alignment.
    decltype(auto) alignment() &&
    {
        return get_impl(field_constant<seqan3::field::alignment>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::alignment
    decltype(auto) alignment() const &&
    {
        return get_impl(field_constant<seqan3::field::alignment>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::alignment
    decltype(auto) alignment() &
    {
        return get_impl(field_constant<seqan3::field::alignment>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::alignment
    decltype(auto) alignment() const &
    {
        return get_impl(field_constant<seqan3::field::alignment>{}, static_cast<tuple_base_t const &>(*this));
    }

    /*!\brief The identifier of the (reference) sequence that seqan3::sam_record::sequence was aligned to.
     *        (SAM Column ID: RNAME)
     */
    decltype(auto) reference_id() &&
    {
        return get_impl(field_constant<seqan3::field::ref_id>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::reference_id
    decltype(auto) reference_id() const &&
    {
        return get_impl(field_constant<seqan3::field::ref_id>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::reference_id
    decltype(auto) reference_id() &
    {
        return get_impl(field_constant<seqan3::field::ref_id>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::reference_id
    decltype(auto) reference_id() const &
    {
        return get_impl(field_constant<seqan3::field::ref_id>{}, static_cast<tuple_base_t const &>(*this));
    }

    /*!\brief The (reference) "sequence" information, usually a range of nucleotides or amino acids. (Currently not
     *        implemented!)
     */
    decltype(auto) reference_sequence() = delete;

    /*!\brief (Reference) Sequence (seqan3::sam_record::reference_sequence) relative start position (0-based),
     *        unsigned value. (SAM Column ID: POS)
     */
    decltype(auto) reference_position() &&
    {
        return get_impl(field_constant<seqan3::field::ref_offset>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::reference_position
    decltype(auto) reference_position() const &&
    {
        return get_impl(field_constant<seqan3::field::ref_offset>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::reference_position
    decltype(auto) reference_position() &
    {
        return get_impl(field_constant<seqan3::field::ref_offset>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::reference_position
    decltype(auto) reference_position() const &
    {
        return get_impl(field_constant<seqan3::field::ref_offset>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief A pointer to the seqan3::sam_file_header object storing header information.
    decltype(auto) header_ptr() &&
    {
        return get_impl(field_constant<seqan3::field::header_ptr>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::header_ptr
    decltype(auto) header_ptr() const &&
    {
        return get_impl(field_constant<seqan3::field::header_ptr>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::header_ptr
    decltype(auto) header_ptr() &
    {
        return get_impl(field_constant<seqan3::field::header_ptr>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::header_ptr
    decltype(auto) header_ptr() const &
    {
        return get_impl(field_constant<seqan3::field::header_ptr>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief The alignment flag (bit information), `uint16_t` value. (SAM Column ID: FLAG)
    decltype(auto) flag() &&
    {
        return get_impl(field_constant<seqan3::field::flag>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::flag
    decltype(auto) flag() const &&
    {
        return get_impl(field_constant<seqan3::field::flag>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::flag
    decltype(auto) flag() &
    {
        return get_impl(field_constant<seqan3::field::flag>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::flag
    decltype(auto) flag() const &
    {
        return get_impl(field_constant<seqan3::field::flag>{}, static_cast<tuple_base_t const &>(*this));
    }

    /*!\brief The identifier of the (reference) sequence of the mate. (SAM Column ID: RNEXT)
     *
     * \details
     *
     * If `RNEXT` is `=`, it returns the same as seqan3::sam_record::reference_id.
     */
    decltype(auto) mate_reference_id() &&
    {
        return std::get<0>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t &&>(*this)));
    }
    //!\copydoc seqan3::sam_record::mate_reference_id
    decltype(auto) mate_reference_id() const &&
    {
        return std::get<0>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t const &&>(*this)));
    }
    //!\copydoc seqan3::sam_record::mate_reference_id
    decltype(auto) mate_reference_id() &
    {
        return std::get<0>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t &>(*this)));
    }
    //!\copydoc seqan3::sam_record::mate_reference_id
    decltype(auto) mate_reference_id() const &
    {
        return std::get<0>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t const &>(*this)));
    }

    //!\brief (Reference) Sequence relative start position (0-based) of the mate. (SAM Column ID: PNEXT)
    decltype(auto) mate_position() &&
    {
        return std::get<1>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t &&>(*this)));
    }
    //!\copydoc seqan3::sam_record::mate_position
    decltype(auto) mate_position() const &&
    {
        return std::get<1>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t const &&>(*this)));
    }
    //!\copydoc seqan3::sam_record::mate_position
    decltype(auto) mate_position() &
    {
        return std::get<1>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t &>(*this)));
    }
    //!\copydoc seqan3::sam_record::mate_position
    decltype(auto) mate_position() const &
    {
        return std::get<1>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t const &>(*this)));
    }

    //!\brief The observed template length. (SAM Column ID: TLEN)
    decltype(auto) template_length() &&
    {
        return std::get<2>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t &&>(*this)));
    }
    //!\copydoc seqan3::sam_record::template_length
    decltype(auto) template_length() const &&
    {
        return std::get<2>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t const &&>(*this)));
    }
    //!\copydoc seqan3::sam_record::template_length
    decltype(auto) template_length() &
    {
        return std::get<2>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t &>(*this)));
    }
    //!\copydoc seqan3::sam_record::template_length
    decltype(auto) template_length() const &
    {
        return std::get<2>(get_impl(field_constant<seqan3::field::mate>{}, static_cast<tuple_base_t const &>(*this)));
    }

    //!\brief The mapping quality of the alignment, usually a Phred-scaled score. (SAM Column ID: MAPQ)
    decltype(auto) mapping_quality() &&
    {
        return get_impl(field_constant<seqan3::field::mapq>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::mapping_quality
    decltype(auto) mapping_quality() const &&
    {
        return get_impl(field_constant<seqan3::field::mapq>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::mapping_quality
    decltype(auto) mapping_quality() &
    {
        return get_impl(field_constant<seqan3::field::mapq>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::mapping_quality
    decltype(auto) mapping_quality() const &
    {
        return get_impl(field_constant<seqan3::field::mapq>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief The cigar vector (std::vector<seqan3::cigar>) representing the alignment. (SAM Column ID: CIGAR)
    decltype(auto) cigar_sequence() &&
    {
        return get_impl(field_constant<seqan3::field::cigar>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::cigar_sequence
    decltype(auto) cigar_sequence() const &&
    {
        return get_impl(field_constant<seqan3::field::cigar>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::cigar_sequence
    decltype(auto) cigar_sequence() &
    {
        return get_impl(field_constant<seqan3::field::cigar>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::cigar_sequence
    decltype(auto) cigar_sequence() const &
    {
        return get_impl(field_constant<seqan3::field::cigar>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief The optional tags in the SAM format, stored in a seqan3::sam_tag_dictionary.
    decltype(auto) tags() &&
    {
        return get_impl(field_constant<seqan3::field::tags>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sam_record::tags
    decltype(auto) tags() const &&
    {
        return get_impl(field_constant<seqan3::field::tags>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sam_record::tags
    decltype(auto) tags() &
    {
        return get_impl(field_constant<seqan3::field::tags>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sam_record::tags
    decltype(auto) tags() const &
    {
        return get_impl(field_constant<seqan3::field::tags>{}, static_cast<tuple_base_t const &>(*this));
    }
};
} // namespace seqan3

namespace std
{

/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::unary_type_trait
 * \relates seqan3::sam_record
 * \see std::tuple_size_v
 */
template <typename field_types, typename field_ids>
struct tuple_size<seqan3::sam_record<field_types, field_ids>>
    : tuple_size<typename seqan3::sam_record<field_types, field_ids>::base_type>
{};

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::transformation_trait
 * \relates seqan3::sam_record
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template <size_t elem_no, typename field_types, typename field_ids>
struct tuple_element<elem_no, seqan3::sam_record<field_types, field_ids>>
    : tuple_element<elem_no, typename seqan3::sam_record<field_types, field_ids>::base_type>
{};

} // namespace std
