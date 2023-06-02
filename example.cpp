#include <array>
#include <cstdint>
#include <utility>
#include "constexpr_huffman.h"
	
namespace input {
    inline constexpr const char cmg_lee[]{ "A_DEAD_DAD_CEDED_A_BAD_BABE_A_BEADED_ABACA_BED" };
    inline constexpr std::uint8_t numbers[]{ 2,2, 4,4,4,4, 8,8,8,8,8,8,8,8, 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16 };
    inline constexpr const char number_string[]{ "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100" };
}
	
namespace constexpr_for {
    // https://vittorioromeo.info/index/blog/cpp20_lambdas_compiletime_for.html
	
    template <auto... Xs, typename F>
    constexpr void for_values(F&& f) {
        (f.template operator()<Xs>(), ...);
    }
	
    template <auto B, auto E, typename F>
    constexpr void for_range(F&& f)
    {
        using t = std::common_type_t<decltype(B), decltype(E)>;
	
        [&f] <auto... Xs>(std::integer_sequence<t, Xs...>) {
            for_values<(B + Xs)...>(f);
        }
        (std::make_integer_sequence<t, E - B>{});
    }
}
	
int main()
{
    using output_t = huffman::metadata::output_t;
	
    static constexpr std::array flags
    {
        /* 0  */ output_t::weights,
        /* 1  */ output_t::tree,
        /* 2  */ output_t::optimized_compression,
        /* 3  */ output_t::compress_to_bytes,
        /* 4  */ output_t::compress_to_bits,
        /* 5  */ output_t::all,
	
        /* 6  */ output_t::weights | output_t::tree,
        /* 7  */ output_t::weights | output_t::tree | output_t::optimized_compression,
        /* 8  */ output_t::tree | output_t::optimized_compression,
        /* 9  */ output_t::tree | output_t::compress_to_bytes,
        /* 10 */ output_t::tree | output_t::compress_to_bits,
        /* 11 */ output_t::tree | output_t::optimized_compression | output_t::compress_to_bytes,
        /* 12 */ output_t::tree | output_t::optimized_compression | output_t::compress_to_bits
    };
	
    constexpr_for::for_range<0, std::size(flags)>([]<auto X>()
    {
        static constexpr auto metadata = []() { return huffman::metadata::create_metadata<flags[X]>(input::cmg_lee); };
        [[maybe_unused]] static constexpr auto output = huffman::start(metadata);
    });
}
