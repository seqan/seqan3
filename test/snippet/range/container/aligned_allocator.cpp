#include <iostream>
#include <vector>

#include <seqan3/range/container/aligned_allocator.hpp>

size_t memory_alignment(void * value, size_t alignment)
{
   return (reinterpret_cast<size_t>(value) & (alignment - 1));
}

int main()
{
    using namespace seqan3;

    // 128-byte memory aligned and 16bit = 2byte address width for each element
    std::vector<int16_t, aligned_allocator<int16_t, 128u>> vec128{1,2,3,4,5};

    // vector has no alignment and 16bit = 2byte address width for each element
    std::vector<int16_t> vec_unaligned{1,2,3,4,5};

    // 256-byte memory aligned and 32bit = 4byte address width for each element
    std::vector<int32_t, aligned_allocator<int32_t, 256u>> vec256{1,2,3,4,5};

    for (auto && x: vec128)
        std::cout << "Item: " << x << " (" << &x << ", 128-byte aligned offset: " << memory_alignment(&x, 128u) << ")\n";

    for (auto && x: vec_unaligned)
        std::cout << "Item: " << x << " (" << &x << ", unaligned start: " << memory_alignment(&x, 128u) << ")\n";

    for (auto && x: vec256)
        std::cout << "Item: " << x << " (" << &x << ", 256-byte aligned offset: " << memory_alignment(&x, 256u) << ")\n";
}
