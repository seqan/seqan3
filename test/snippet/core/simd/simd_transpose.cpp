#include <array>

#include <seqan3/core/simd/all.hpp>

using namespace seqan3;

using uint8x4_t = simd_type_t<uint8_t, 4>;

int main()
{
    std::array<uint8x4_t, 4> matrix{{uint8x4_t{0, 1, 2, 3},
                                     uint8x4_t{0, 1, 2, 3},
                                     uint8x4_t{0, 1, 2, 3},
                                     uint8x4_t{0, 1, 2, 3}}};
    transpose(matrix);

    debug_stream << matrix[0] << "\n"; // [0,0,0,0]
    debug_stream << matrix[1] << "\n"; // [1,1,1,1]
    debug_stream << matrix[2] << "\n"; // [2,2,2,2]
    debug_stream << matrix[3] << "\n"; // [3,3,3,3]
    return 0;
}
