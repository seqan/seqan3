#include <seqan3/alphabet/structure/dssp9.hpp>

int main()
{
    using seqan3::operator""_dssp9;

    // Using the char literal to assign a single DSSP annotation:
    seqan3::dssp9 my_letter{'I'_dssp9};
    
    my_letter.assign_char('G');      // <- assigns the char explicitly
}
