# constexpr-huffman
A header-only C++20 library that facilitates generating Huffman trees from frequency input, and compressing arrays and string literals into a stream of bits at compile-time.

# Why?
Huffman compression isn't new and there are many examples that show off how to generate Huffman trees. 
Unfortunately, those examples are often written in C or C-style C++. 
Many of those perform badly and / or use raw pointers extensively, which tends to lead to confusion, and which may be incompatible with the goal of doing everything at compile-time; e.g. reinterpret_cast. 
Some C++ examples use standard map, but this container isn't usable at compile-time.
Standard vector and sorting are usable at compile-time in C++20, and compared to using raw pointers everywhere they reduce the code complexity significantly.

If there are other constexpr Huffman libraries publicly available, they either do not meet my needs, or take too long to compile or I'm simply unaware of them.

# How to use
Check out the example file or the Compiler Explorer example below.

https://godbolt.org/z/5Kdj8fd3Y

There are two reasons the result is different from https://en.wikipedia.org/wiki/File:Huffman_coding_visualisation.svg. The first reason is the inclusion of the null terminator, and the second reason is the incorrect switch of '\_' and 'D' characters. Since 'D' comes before '\_' it should get the lower value, so [0 0] for 'D' and [0 1] for '\_'.
