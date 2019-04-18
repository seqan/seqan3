
#pragma once

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>

namespace seqan3::test::alignment::fixture
{

static constexpr auto INF = detail::matrix_inf<int>;

static constexpr auto NON = detail::trace_directions::none;
static constexpr auto D = detail::trace_directions::diagonal;
static constexpr auto U = detail::trace_directions::up;
static constexpr auto L = detail::trace_directions::left;
static constexpr auto DU = D | U;
static constexpr auto UL = U | L;
static constexpr auto DL = D | L;
static constexpr auto DUL = D | U | L;

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t, typename trace_t>
struct alignment_fixture
{
    sequence1_t sequence1;
    sequence2_t sequence2;

    config_t config;

    score_t score;
    std::string aligned_sequence1;
    std::string aligned_sequence2;

    alignment_coordinate front_coordinate;
    alignment_coordinate back_coordinate;

    std::vector<score_t> score_vector{};
    std::vector<trace_t> trace_vector{};

    detail::alignment_score_matrix<std::vector<score_t>> score_matrix() const
    {
        return {score_vector, sequence2.size() + 1, sequence1.size() + 1};
    };

    detail::alignment_trace_matrix<std::vector<trace_t>> trace_matrix() const
    {
        return {trace_vector, sequence2.size() + 1, sequence1.size() + 1};
    };
};

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t, typename trace_t>
alignment_fixture(
    sequence1_t sequence1,
    sequence2_t sequence2,
    config_t config,
    score_t score,
    std::string aligned_sequence1,
    std::string aligned_sequence2,
    alignment_coordinate front_coordinate,
    alignment_coordinate back_coordinate,
    std::vector<score_t> score_vector,
    std::vector<trace_t> trace_vector
)
-> alignment_fixture<sequence1_t, sequence2_t, config_t, score_t, trace_t>;

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t>
alignment_fixture(
    sequence1_t sequence1,
    sequence2_t sequence2,
    config_t config,
    score_t score,
    std::string aligned_sequence1,
    std::string aligned_sequence2,
    alignment_coordinate front_coordinate,
    alignment_coordinate back_coordinate
)
-> alignment_fixture<sequence1_t, sequence2_t, config_t, score_t, detail::trace_directions>;


} // namespace seqan3::test::alignment::fixture
