#include <vector>

#include <seqan3/alignment/matrix/detail/alignment_matrix_column_major_range_base.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/span>

class my_matrix : public seqan3::detail::alignment_matrix_column_major_range_base<my_matrix>
{
public:

    // Alias the base class
    using base_t = seqan3::detail::alignment_matrix_column_major_range_base<my_matrix>;

    friend base_t;

    // Inherit the alignment column type defined in the base class. This type is returned in initialise_column.
    using typename base_t::alignment_column_type;

    // The following types are required by the base type since they cannot be inferred within the base.

    using column_data_view_type = std::span<int>;  //This type is the underlying view over the actual memory location.
    using value_type = int; // The actual value type.
    using reference = int &; // The actual reference type.

    my_matrix() = default;
    my_matrix(my_matrix const &) = default;
    my_matrix(my_matrix &&) = default;
    my_matrix & operator=(my_matrix const &) = default;
    my_matrix & operator=(my_matrix &&) = default;
    ~my_matrix() = default;

    my_matrix(size_t const num_rows, size_t const num_cols) :
        num_rows{num_rows},
        num_cols{num_cols}
    {
        data.resize(num_rows * num_cols);
    }

protected:

    std::vector<int> data{};
    size_t num_rows{};
    size_t num_cols{};

    //Required for the base class. Initialises the current column given the column index.
    alignment_column_type initialise_column(size_t const column_index) noexcept
    {
        return alignment_column_type{*this,
                                     column_data_view_type{std::addressof(data[num_rows * column_index]), num_rows}};
    }

    //Required for the base class. Initialises the proxy for the current iterator over the current column.
    template <std::random_access_iterator iter_t>
    constexpr reference make_proxy(iter_t iter) noexcept
    {
        return *iter;
    }
};

int main()
{
    my_matrix matrix{3, 5};

    // Fill the matrix with
    int val = 0;
    for (auto col : matrix) // Iterate over the columns
        for (auto & cell : col) // Iterate over the cells in one column.
            cell = val++;

    // Print the matrix column by column
    for (auto col : matrix)
        seqan3::debug_stream << col << '\n';
}
