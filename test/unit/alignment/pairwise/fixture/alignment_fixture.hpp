
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

template <typename sequence1_t,
          typename sequence2_t,
          typename config_t,
          typename score_t,
          typename score_matrix_t,
          typename trace_matrix_t>
struct alignment_fixture
{
    template <typename trace_t>
    alignment_fixture(sequence1_t && _sequence1,
                      sequence2_t && _sequence2,
                      config_t _config,
                      score_t _score,
                      std::string _aligned_sequence1,
                      std::string _aligned_sequence2,
                      alignment_coordinate _front_coordinate,
                      alignment_coordinate _back_coordinate,
                      std::vector<score_t> _score_vector,
                      std::vector<trace_t> _trace_vector) :
        sequence1{_sequence1},
        sequence2{_sequence2},
        config{_config},
        score{_score},
        aligned_sequence1{_aligned_sequence1},
        aligned_sequence2{_aligned_sequence2},
        front_coordinate{_front_coordinate},
        back_coordinate{_back_coordinate},
        score_matrix{std::move(_score_vector), sequence2.size()+1, sequence1.size()+1},
        trace_matrix{std::move(_trace_vector), sequence2.size()+1, sequence1.size()+1}
    {}

    alignment_fixture(sequence1_t && _sequence1,
                      sequence2_t && _sequence2,
                      config_t _config,
                      score_t _score,
                      std::string _aligned_sequence1,
                      std::string _aligned_sequence2,
                      alignment_coordinate _front_coordinate,
                      alignment_coordinate _back_coordinate) :
        alignment_fixture(std::forward<sequence1_t>(_sequence1),
                          std::forward<sequence2_t>(_sequence2),
                          _config,
                          _score,
                          _aligned_sequence1,
                          _aligned_sequence2,
                          _front_coordinate,
                          _back_coordinate,
                          std::vector<score_t>{},
                          std::vector<detail::trace_directions>{})
    {}

    sequence1_t sequence1;
    sequence2_t sequence2;

    config_t config;

    score_t score;
    std::string aligned_sequence1;
    std::string aligned_sequence2;

    alignment_coordinate front_coordinate;
    alignment_coordinate back_coordinate;
    score_matrix_t score_matrix;
    trace_matrix_t trace_matrix;
};

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t, typename trace_t>
alignment_fixture(sequence1_t && _sequence1,
                  sequence2_t && _sequence2,
                  config_t _config,
                  score_t _score,
                  std::string _aligned_sequence1,
                  std::string _aligned_sequence2,
                  alignment_coordinate _front_coordinate,
                  alignment_coordinate _back_coordinate,
                  std::vector<score_t> _score_vector,
                  std::vector<trace_t> _trace_vector)
-> alignment_fixture<sequence1_t,
                     sequence2_t,
                     config_t,
                     score_t,
                     detail::alignment_score_matrix<std::vector<score_t>>,
                     detail::alignment_trace_matrix<std::vector<trace_t>>>;

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t>
alignment_fixture(sequence1_t && _sequence1,
                  sequence2_t && _sequence2,
                  config_t _config,
                  score_t _score,
                  std::string _aligned_sequence1,
                  std::string _aligned_sequence2,
                  alignment_coordinate _front_coordinate,
                  alignment_coordinate _back_coordinate)
-> alignment_fixture<sequence1_t,
                     sequence2_t,
                     config_t,
                     score_t,
                     detail::alignment_score_matrix<std::vector<score_t>>,
                     detail::alignment_trace_matrix<std::vector<detail::trace_directions>>>;

} // namespace seqan3::test::alignment::fixture
