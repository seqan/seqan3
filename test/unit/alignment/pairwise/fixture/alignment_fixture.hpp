
#pragma once

#include <limits>

#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>

namespace seqan3::literal
{}

namespace seqan3::fixture
{
using namespace seqan3::literal;
using detail::alignment_coordinate;

static constexpr auto INF = std::numeric_limits<int>::max();

static constexpr auto NON = detail::trace_directions::none;
static constexpr auto D = detail::trace_directions::diagonal;
static constexpr auto U = detail::trace_directions::up;
static constexpr auto L = detail::trace_directions::left;
static constexpr auto DU = D | U;
static constexpr auto UL = U | L;
static constexpr auto DL = D | L;
static constexpr auto DUL = D | U | L;

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t, typename score_matrix_t, typename trace_matrix_t>
struct alignment_fixture
{
    template <typename trace_t>
    alignment_fixture(
        sequence1_t && _sequence1,
        sequence2_t && _sequence2,
        config_t _config,
        score_t _score,
        std::string _gapped_sequence1,
        std::string _gapped_sequence2,
        alignment_coordinate _begin_coordinate,
        alignment_coordinate _end_coordinate,
        std::vector<score_t> _score_vector,
        std::vector<trace_t> _trace_vector
    ) : sequence1{_sequence1},
        sequence2{_sequence2},
        config{_config},
        score{_score},
        gapped_sequence1{_gapped_sequence1},
        gapped_sequence2{_gapped_sequence2},
        begin_coordinate{_begin_coordinate},
        end_coordinate{_end_coordinate},
        score_matrix{std::move(_score_vector), sequence2.size()+1, sequence1.size()+1},
        trace_matrix{std::move(_trace_vector), sequence2.size()+1, sequence1.size()+1}
    {
    }

    sequence1_t sequence1;
    sequence2_t sequence2;

    config_t config;

    score_t score;
    std::string gapped_sequence1;
    std::string gapped_sequence2;

    alignment_coordinate begin_coordinate;
    alignment_coordinate end_coordinate;
    score_matrix_t score_matrix;
    trace_matrix_t trace_matrix;

};

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t, typename trace_t>
alignment_fixture(
    sequence1_t && _sequence1,
    sequence2_t && _sequence2,
    config_t _config,
    score_t _score,
    std::string _gapped_sequence1,
    std::string _gapped_sequence2,
    alignment_coordinate _begin_coordinate,
    alignment_coordinate _end_coordinate,
    std::vector<score_t> _score_vector,
    std::vector<trace_t> _trace_vector
)
-> alignment_fixture<
    sequence1_t,
    sequence2_t,
    config_t,
    score_t,
    detail::alignment_score_matrix<std::vector<score_t>>,
    detail::alignment_trace_matrix<std::vector<trace_t>>
>;

} // namespace seqan3::fixture
