#pragma once

#include <iostream>
#include <string>

namespace seqan3
{

// ------------------------------------------------------------------
// concept
// ------------------------------------------------------------------

    template <typename T>
    concept bool phred_concept() {
        return requires(T t){

        };

    }
}