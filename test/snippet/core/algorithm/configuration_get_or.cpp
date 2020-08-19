#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/debug_stream.hpp>

// Initial setup used in the actual example:
enum struct my_id : int
{
    bar_id,
    foo_id
};

struct bar : public seqan3::pipeable_config_element<bar, float>
{
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr bar() = default; //!< Defaulted.
    constexpr bar(bar const &) = default; //!< Defaulted.
    constexpr bar(bar &&) = default; //!< Defaulted.
    constexpr bar & operator=(bar const &) = default; //!< Defaulted.
    constexpr bar & operator=(bar &&) = default; //!< Defaulted.
    ~bar() = default; //!< Defaulted.

    //!\brief Construct from base type.
    constexpr bar(float const & val) : seqan3::pipeable_config_element<bar, float>(val) {}
    //!\}

    static constexpr my_id id{my_id::bar_id};
};

template <typename t>
struct foo : public seqan3::pipeable_config_element<foo<t>, t>
{
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr foo() = default; //!< Defaulted.
    constexpr foo(foo const &) = default; //!< Defaulted.
    constexpr foo(foo &&) = default; //!< Defaulted.
    constexpr foo & operator=(foo const &) = default; //!< Defaulted.
    constexpr foo & operator=(foo &&) = default; //!< Defaulted.
    ~foo() = default; //!< Defaulted.

    //!\brief Construct from base type.
    constexpr foo(t const & val) : seqan3::pipeable_config_element<foo<t>, t>(val) {}
    //!\}

    static constexpr my_id id{my_id::foo_id};
};

template <typename t>
foo(t) -> foo<t>;

int main()
{
    seqan3::configuration my_cfg{foo{1}}; // Only foo<int> is present.
    seqan3::debug_stream << my_cfg.get_or(foo{std::string{"hello"}}).value << '\n';  // finds foo<int> -> prints: 1
    seqan3::debug_stream << my_cfg.get_or(bar{2.4}).value << '\n';  // bar not present -> prints: 2.4
}
