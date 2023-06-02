#ifndef CONSTEXPR_HUFFMAN_H
#define CONSTEXPR_HUFFMAN_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace huffman 
{
    namespace globals 
    {
        inline constexpr std::size_t byte_max = std::numeric_limits<std::uint8_t>::max() + 1;
        inline constexpr std::size_t int_node_val = 999;
        inline constexpr std::size_t no_node_val = 999;
    }

    namespace metadata 
    {
        using namespace globals;

        enum struct input_t : std::uint8_t
        {
            data = 0,
            weights
        };

        enum struct output_t : std::uint8_t
        {
            weights = 1,
            tree = 2,
            optimized_compression = 4,
            compress_to_bytes = 8,
            compress_to_bits = 16,
            all = 31,
            skip_compressed_check = 32
        };

        constexpr auto operator|(const output_t lhs, const output_t rhs)
        {
            using output_type_t = std::underlying_type_t<output_t>;
            return output_t(static_cast<output_type_t>(lhs) | static_cast<output_type_t>(rhs));
        }

        constexpr auto operator&(const output_t lhs, const output_t rhs)
        {
            using output_type_t = std::underlying_type_t<output_t>;
            return output_t(static_cast<output_type_t>(lhs) & static_cast<output_type_t>(rhs));
        }

        constexpr auto operator~(const output_t val)
        {
            using output_type_t = std::underlying_type_t<output_t>;
            return output_t(~static_cast<output_type_t>(val));
        }

        template <typename T, std::size_t input_size, std::size_t weight_size>
        struct metadata_t
        {
            constexpr metadata_t(input_t ipt, output_t opt, std::size_t byte_count, std::size_t node_count,
                const std::span<T, input_size> input, std::array<std::size_t, weight_size>&& weights) :
                ipt(ipt),
                opt(opt),
                byte_count(byte_count),
                node_count(node_count),
                input(input),
                weights(std::forward<decltype(weights)>(weights))
            {}

            input_t ipt = input_t::data;
            output_t opt = output_t::weights;
            std::size_t byte_count = 0;
            std::size_t node_count = 0;
            std::span<T, input_size> input;
            std::array<std::size_t, weight_size> weights;
        };

        template <output_t opt, input_t ipt = input_t::data>
        constexpr auto create_metadata(const auto& input)
        {
            constexpr auto return_data = [](auto&& weights, const auto span) {
                const std::size_t byte_count = std::count_if(std::begin(weights), std::end(weights), [](const auto val) { return val != 0; });
                const std::size_t node_count = byte_count * 2 - 1;

                if (node_count < 2) {
                    throw std::runtime_error("input consists of fewer than two different input values");
                }

                return metadata_t(ipt, opt, byte_count, node_count, span, std::forward<decltype(weights)>(weights));
            };

            using data_type = std::remove_reference_t<decltype(input[0])>;
            const std::span span{ input };

            if constexpr (ipt != input_t::weights) {
                static_assert(std::is_same_v<data_type, const char> || std::is_same_v<data_type, const unsigned char>);

                constexpr auto get_weights = [](const auto span) {
                    std::array<std::size_t, byte_max> weights{};

                    for (const std::uint8_t byte : span) {
                        ++weights[byte];
                    }

                    return weights;
                };

                return return_data(get_weights(span), span);
            } else {
                static_assert(std::size(span) == byte_max && std::is_unsigned_v<data_type>);

                const auto copy_weights = [](const auto span) {
                    std::array<std::size_t, std::size(span)> weights;
                    std::copy_n(std::begin(span), std::size(span), std::begin(weights));

                    return weights;
                };

                return return_data(copy_weights(span), span);
            }
        }
    }

    namespace tree_gen 
    {
        using namespace globals;

        template <bool, typename T = std::size_t>
        struct node_t
        {
            T symbol = 0;
            T left = 0;
            T right = 0;
            T parent = 0;
        };

        template <typename T>
        struct node_t<true, T>
        {
            T symbol = 0;
            T left = 0;
            T right = 0;
        };

        struct internal_node_t
        {
            std::size_t weight = 0;
            node_t<false> node;
        };
        using i_node_t = internal_node_t;

        struct node_ptrs_t
        {
            constexpr node_ptrs_t(std::size_t size)
            {
                ptrs.reserve(size);
            }

            std::size_t index = 0;
            std::vector<i_node_t*> ptrs;
        };

        template <std::size_t byte_count, std::size_t node_count, bool get_size>
        struct tree_t
        {
            constexpr tree_t(const auto weights)
            {
                std::array<node_ptrs_t, 2> node_ptrs{ node_ptrs_t(byte_count), node_ptrs_t(byte_count) }; // leave node pointers and internal node pointers

                create_leave_nodes<byte_count>(node_buffer, node_ptrs[0].ptrs, weights);
                create_internal_nodes(node_ptrs, node_buffer, byte_count);

                if constexpr (get_size) {
                    bitcode_size = get_bitcode_size<byte_count>(node_buffer);
                }
            }

            std::size_t bitcode_size = 0;
            std::array<i_node_t, node_count> node_buffer;
        };

        template <std::size_t byte_count>
        constexpr auto get_bitcode_size(const auto& node_buffer)
        {
            std::size_t size = 0;

            for (std::size_t current_size = 0, i = 0; i < byte_count; ++i, current_size = 0) {
                std::size_t parent_index = i;
                const auto* p_node = &node_buffer[i];

                while (p_node->node.parent != no_node_val) {
                    if (parent_index == node_buffer[p_node->node.parent].node.left || parent_index == node_buffer[p_node->node.parent].node.right) {
                        ++current_size;
                    } else {
                        throw std::runtime_error("previous node is no child node");
                    }

                    parent_index = p_node->node.parent;
                    p_node = &node_buffer[p_node->node.parent];
                }

                if (current_size > size) {
                    size = current_size;
                }
            }

            if (size > 64) {
                throw std::runtime_error("the max bitcode size shouldn't exceed 64 because weights are stored in 64 uint");
            }

            return size;
        }

        constexpr auto sort_node_ptrs(auto& ptrs)
        {
            std::sort(std::begin(ptrs), std::end(ptrs), [](const auto lhs, const auto rhs) {
                if (lhs->weight == rhs->weight) {
                    return lhs->node.symbol < rhs->node.symbol; // preserve the order of equal elements
                }

                return lhs->weight < rhs->weight;
            });
        }

        constexpr auto verify_node_pointers(const auto& node_buffer, const auto& nodes)
        {
            for (const auto node : nodes) {
                if (node == nullptr) {
                    throw std::runtime_error("nullptr node detected");
                }
                if ((node->node.left != no_node_val && node->node.left >= std::size(node_buffer)) || (node->node.right != no_node_val && node->node.right >= std::size(node_buffer))) {
                    throw std::runtime_error("index is oob");
                }
            }
        }

        constexpr auto is_sorted(const auto& nodes_ptrs)
        {
            if (std::size(nodes_ptrs) < 2) {
                return true;
            }

            // the last element of a vector that was sorted in ascending order should be the largest
            const auto upper = std::upper_bound(std::begin(nodes_ptrs), std::end(nodes_ptrs) - 1, nodes_ptrs.back(), [](const auto last_node, const auto node) {
                if (node->weight == last_node->weight) {
                    if (node->node.symbol == int_node_val && last_node->node.symbol == int_node_val) {
                        return node->node.symbol > last_node->node.symbol;
                    } else {
                        return node->node.symbol >= last_node->node.symbol; // preserve the order of equal elements
                    }
                }

                return last_node->weight < node->weight;
            });

            return (upper == std::end(nodes_ptrs) - 1);
        }

        constexpr auto eval_and_collect_nodes(auto& nodes, auto& node_ptrs)
        {
            for (std::size_t i = 0; node_ptrs.index < std::size(node_ptrs.ptrs) && i < 2; ++i, ++node_ptrs.index) {
                if (std::size(nodes) >= 2 && node_ptrs.ptrs[node_ptrs.index]->weight >= nodes[1]->weight) {
                    break; // don't collect pointers if the weight of the second pointer is larger than the weight of the current pointer in the vector
                }

                nodes.emplace_back(node_ptrs.ptrs[node_ptrs.index]);
                sort_node_ptrs(nodes);
            }

            assert(is_sorted(nodes) && std::size(nodes) <= 4);
        }

        constexpr auto create_internal_nodes(auto& node_ptrs, auto& node_buffer, std::size_t buffer_index)
        {
            const auto* buffer_begin = &node_buffer[0];

            std::vector<i_node_t*> nodes;
            nodes.reserve(4);

            while (true) {
                eval_and_collect_nodes(nodes, node_ptrs[0]); // check leaves
                eval_and_collect_nodes(nodes, node_ptrs[1]); // check internal nodes

                if (std::size(nodes) < 2) {
                    break;
                }

                verify_node_pointers(node_buffer, nodes);

                nodes[0]->node.parent = buffer_index;
                nodes[1]->node.parent = buffer_index;

                node_buffer[buffer_index] = i_node_t{
                    .weight = nodes[0]->weight + nodes[1]->weight,
                    .node{
                        .symbol = int_node_val,
                        .left = static_cast<std::size_t>(nodes[0] - buffer_begin),
                        .right = static_cast<std::size_t>(nodes[1] - buffer_begin),
                        .parent = no_node_val
                    }
                };

                node_ptrs[1].ptrs.emplace_back(&node_buffer[buffer_index]);
                ++buffer_index;

                assert(is_sorted(node_ptrs[1].ptrs));
                nodes.erase(std::begin(nodes), std::begin(nodes) + 2);
            }

            if (node_ptrs[0].index + node_ptrs[1].index == std::size(node_buffer)) {
                assert(std::size(nodes) == 1 && nodes[0] == &node_buffer.back());
                assert(is_sorted(node_ptrs[1].ptrs));
            } else {
                throw std::runtime_error("incorrect number of internal nodes created");
            }
        }

        template <std::size_t byte_count>
        constexpr auto create_leave_nodes(auto& node_buffer, auto& node_ptrs, const auto weights)
        {
            std::size_t buffer_index = 0;

            for (std::size_t i = 0; i < std::size(weights); ++i) {
                if (weights[i] == 0) {
                    continue;
                }

                node_buffer[buffer_index] = i_node_t{
                    .weight = weights[i],
                    .node{
                        .symbol = i,
                        .left = no_node_val,
                        .right = no_node_val,
                        .parent = no_node_val
                    }
                };

                node_ptrs.emplace_back(&node_buffer[buffer_index]);
                ++buffer_index;
            }

            if (buffer_index != byte_count) {
                throw std::runtime_error("incorrect number of internal nodes created");
            }

            sort_node_ptrs(node_ptrs);
        }

        template <std::size_t node_count, bool skip_parent_nodes>
        constexpr auto optimize_tree_size(const auto& node_buffer)
        {
            using type = std::conditional_t<(node_count <= std::numeric_limits<std::uint8_t>::max()), std::uint8_t, std::uint16_t>;
            constexpr std::size_t max = std::numeric_limits<type>::max();

            std::array<node_t<skip_parent_nodes, type>, node_count> nodes;

            for (std::size_t i = 0; const auto & node : node_buffer) {
                nodes[i].symbol = (node.node.symbol > max) ? static_cast<type>(max) : static_cast<type>(node.node.symbol);
                nodes[i].left = (node.node.left > max) ? static_cast<type>(max) : static_cast<type>(node.node.left);
                nodes[i].right = (node.node.right > max) ? static_cast<type>(max) : static_cast<type>(node.node.right);

                if constexpr (!skip_parent_nodes) {
                    nodes[i].parent = (node.node.parent > max) ? static_cast<type>(max) : static_cast<type>(node.node.parent);
                }

                ++i;
            }

            return nodes;
        }
    }

    namespace compressor_gen 
    {
        using namespace globals;

        template <bool, std::size_t byte_count, std::size_t node_count, std::size_t bitcode_size, std::size_t weight_size>
        struct compressor_t
        {
            constexpr compressor_t(auto&&)
            {}

            std::array<std::array<std::uint8_t, bitcode_size>, weight_size> bit_data{};
        };

        template <std::size_t byte_count, std::size_t node_count, std::size_t bitcode_size, std::size_t weight_size>
        struct compressor_t<true, byte_count, node_count, bitcode_size, weight_size>
        {
            constexpr compressor_t(auto&& data) : data(std::forward<decltype(data)>(data))
            {}

            tree_gen::tree_t<byte_count, node_count, true> data;
            std::array<std::array<std::uint8_t, bitcode_size>, weight_size> bit_data{};
        };

        template <std::size_t byte_count, std::size_t node_count, bool store_tree, tree_gen::tree_t<byte_count, node_count, true> data>
        constexpr auto create_compressor(const auto weights)
        {
            const std::size_t max = 64; // a weight cannot be larger than 2^64, so the maximum of number of bits to represent a symbol should 64
            const std::size_t bitcode_size = (data.bitcode_size > max) ? max : data.bitcode_size;

            compressor_t<store_tree, byte_count, node_count, bitcode_size + 1, std::size(weights)> comp(data);
            const auto& node_buffer = data.node_buffer;

            for (std::size_t i = 0, j = 0; i < std::size(weights); ++i) {
                if (weights[i] == 0) {
                    continue;
                }

                auto& bit_count = comp.bit_data[i][0];
                auto& bit_data = comp.bit_data[i];
                std::array<std::uint8_t, bitcode_size> local_bit_data{};

                const auto* p_node = &node_buffer[j];
                std::size_t parent_index = j++;

                while (p_node->node.parent != no_node_val) {
                    if (parent_index == node_buffer[p_node->node.parent].node.left) {
                        local_bit_data[bit_count++] = 0;
                    } else if (parent_index == node_buffer[p_node->node.parent].node.right) {
                        local_bit_data[bit_count++] = 1;
                    } else {
                        throw std::runtime_error("previous node is no child node");
                    }

                    parent_index = p_node->node.parent;
                    p_node = &node_buffer[p_node->node.parent];
                }

                std::reverse(std::begin(local_bit_data), std::end(local_bit_data));
                std::copy_n(std::begin(local_bit_data) + (bitcode_size - bit_count), bit_count, std::begin(bit_data) + 1);
            }

            return comp;
        }
    }

    namespace decompress_output 
    {
        using namespace globals;

        constexpr auto check_output(const auto& decompressed, const auto input)
        {
            if (std::size(decompressed) != std::size(input)) {
                throw std::runtime_error("invalid output size after compression and decompression");
            }

            for (std::size_t i = 0; i < std::size(input); ++i) {
                if (input[i] != decompressed[i]) {
                    throw std::runtime_error("invalid output after compression and decompression");
                }
            }
        };

        template <bool return_data = true>
        constexpr auto decompress_output_bytes(const auto& compressed_bytes, const auto& node_buffer, const auto input)
        {
            using data_type = std::remove_reference_t<std::remove_cvref_t<decltype(input[0])>>;
            static_assert(std::is_same_v<data_type, char> || std::is_same_v<data_type, unsigned char>);

            std::array<data_type, std::size(input)> decompressed;
            const auto* node = &node_buffer.back().node;

            for (std::size_t i = 0; const auto byte : compressed_bytes) {
                node = (byte == 0) ? &node_buffer[node->left].node : &node_buffer[node->right].node;

                if (node->symbol != int_node_val) {
                    assert(node->left == no_node_val && node->right == no_node_val);
                    decompressed[i++] = static_cast<data_type>(node->symbol);

                    node = &node_buffer.back().node;
                }
            }

            check_output(decompressed, input);

            if constexpr (return_data) {
                return decompressed;
            }
        };

        template <bool return_data = true>
        constexpr auto decompress_output_bits(const auto& compressed_bits, const auto& node_buffer, const auto input)
        {
            using data_type = std::remove_reference_t<std::remove_cvref_t<decltype(input[0])>>;
            static_assert(std::is_same_v<data_type, char> || std::is_same_v<data_type, unsigned char>);

            std::array<data_type, std::size(input)> decompressed;
            const auto* node = &node_buffer.back().node;

            for (std::size_t i = 0, j = 0; i < std::size(compressed_bits) * 8 && j < std::size(input); ++i) {
                node = (((compressed_bits[i / 8] >> (i & 7)) & 1) == 0) ? &node_buffer[node->left].node : &node_buffer[node->right].node;

                if (node->symbol != int_node_val) {
                    assert(node->left == no_node_val && node->right == no_node_val);
                    decompressed[j++] = static_cast<data_type>(node->symbol);

                    node = &node_buffer.back().node;
                }
            }

            check_output(decompressed, input);

            if constexpr (return_data) {
                return decompressed;
            }
        };
    }

    namespace compress_input 
    {
        using namespace globals;

        template <std::size_t input, std::size_t multiple>
        constexpr std::size_t round_up()
        {
            static_assert(multiple);
            return ((input + multiple - 1) / multiple) * multiple;
        }

        constexpr auto get_conversion_size(const auto& bit_data, const auto weights)
        {
            std::size_t size = 0;

            for (std::size_t i = 0; i < std::size(weights); ++i) {
                if (weights[i] == 0) {
                    continue;
                }

                size += weights[i] * bit_data[i][0];
            }

            return size;
        };

        template <std::size_t conversion_size, bool check_output>
        constexpr auto compress_input_bytes(const auto& bit_data, const auto& node_buffer, const auto input)
        {
            std::array<std::uint8_t, conversion_size> bytes;

            for (std::size_t i = 0; const std::uint8_t byte : input) {
                const std::size_t size = bit_data[byte][0];

                for (std::size_t j = 0; j < size; ++j) {
                    bytes[i++] = bit_data[byte][j + 1];
                }
            }

            if constexpr (check_output) {
                decompress_output::decompress_output_bytes<false>(bytes, node_buffer, input);
            }

            return bytes;
        };

        template <std::size_t conversion_size, bool check_output>
        constexpr auto compress_input_bits(const auto& bit_data, const auto& node_buffer, const auto input)
        {
            std::array<std::uint8_t, round_up<conversion_size, 8>() / 8> bits{};

            for (std::size_t i = 0; const std::uint8_t byte : input) {
                const std::size_t size = bit_data[byte][0];

                for (std::size_t j = 0; j < size; ++i, ++j) {
                    if (bit_data[byte][j + 1]) {
                        bits[i / 8] |= 1 << (i % 8);
                    }
                }
            }

            if constexpr (check_output) {
                decompress_output::decompress_output_bits<false>(bits, node_buffer, input);
            }

            return bits;
        };
    }

    constexpr auto start(auto callable)
    {
        using namespace metadata;
        using namespace tree_gen;
        using namespace compressor_gen;
        using namespace compress_input;
        using namespace decompress_output;

        constexpr auto metadata = callable();
        constexpr bool skip_parent_nodes = static_cast<bool>(metadata.opt & output_t::tree) && static_cast<bool>(metadata.opt & output_t::optimized_compression);
        constexpr bool return_byte_output = static_cast<bool>(metadata.opt & output_t::compress_to_bytes);
        constexpr bool return_bit_output = static_cast<bool>(metadata.opt & output_t::compress_to_bits);
        constexpr bool return_all = metadata.opt == output_t::all;
        constexpr bool skip_check = (metadata.ipt == input_t::weights) ? true : (static_cast<bool>(metadata.opt & output_t::skip_compressed_check));
        constexpr output_t opt = metadata.opt & ~output_t::skip_compressed_check;

        static_assert(metadata.ipt != input_t::weights || (!return_byte_output && !return_bit_output && !return_all && skip_check));
        static_assert(skip_check || opt == metadata.opt);

        if constexpr (opt == output_t::weights) {
            return metadata.weights;

        } else {
            using tree = tree_t<metadata.byte_count, metadata.node_count, true>;
            constexpr auto get_tree = [metadata]() { return tree(std::span{ metadata.weights }); };

            if constexpr (opt == output_t::tree) {
                return optimize_tree_size<metadata.node_count, skip_parent_nodes>(get_tree().node_buffer);

            } else if constexpr (opt == output_t::optimized_compression) {
                constexpr auto get_compression_array = [get_tree, metadata]() { return create_compressor<metadata.byte_count, metadata.node_count, false, get_tree()>(std::span{ metadata.weights }); };
                return get_compression_array().bit_data;

            } else if constexpr (opt == output_t::compress_to_bytes) {
                constexpr auto compression_array = create_compressor<metadata.byte_count, metadata.node_count, true, get_tree()>(std::span{ metadata.weights });
                constexpr auto get_size = [compression_array, metadata]() { return get_conversion_size(compression_array.bit_data, std::span{ metadata.weights }); };
                return compress_input_bytes<get_size(), !skip_check>(
                    compression_array.bit_data, compression_array.data.node_buffer, metadata.input);

            } else if constexpr (opt == output_t::compress_to_bits) {
                constexpr auto compression_array = create_compressor<metadata.byte_count, metadata.node_count, true, get_tree()>(std::span{ metadata.weights });
                constexpr auto get_size = [compression_array, metadata]() { return get_conversion_size(compression_array.bit_data, std::span{ metadata.weights }); };
                return compress_input_bits<get_size(), !skip_check>(
                    compression_array.bit_data, compression_array.data.node_buffer, metadata.input);

            } else if constexpr (opt == output_t::all) {
                constexpr auto compression_array = create_compressor<metadata.byte_count, metadata.node_count, true, get_tree()>(std::span{ metadata.weights });
                constexpr auto get_optimized_tree = [metadata, compression_array]() { return optimize_tree_size<metadata.node_count, skip_parent_nodes>(compression_array.data.node_buffer); };
                constexpr auto size = get_conversion_size(compression_array.bit_data, std::span{ metadata.weights });

                constexpr auto byte_output = compress_input_bytes<size, !skip_check>(
                    compression_array.bit_data, compression_array.data.node_buffer, metadata.input);
                constexpr auto bit_output = compress_input_bits<size, !skip_check>(
                    compression_array.bit_data, compression_array.data.node_buffer, metadata.input);

                return std::tuple{ metadata.weights, get_optimized_tree(), compression_array.bit_data, byte_output, bit_output,
                    decompress_output_bytes(byte_output, compression_array.data.node_buffer, metadata.input),
                    decompress_output_bits(bit_output, compression_array.data.node_buffer, metadata.input)
                };

            } else if constexpr (opt == (output_t::weights | output_t::tree)) {
                constexpr auto get_optimized_tree = [metadata, get_tree]() { return optimize_tree_size<metadata.node_count, skip_parent_nodes>(get_tree().node_buffer); };
                return std::tuple{ metadata.weights, get_optimized_tree() };

            } else if constexpr (opt == (output_t::weights | output_t::tree | output_t::optimized_compression)) {
                constexpr auto compression_array = create_compressor<metadata.byte_count, metadata.node_count, true, get_tree()>(std::span{ metadata.weights });
                constexpr auto get_optimized_tree = [metadata, compression_array]() { return optimize_tree_size<metadata.node_count, skip_parent_nodes>(compression_array.data.node_buffer); };
                return std::tuple{ metadata.weights, get_optimized_tree(), compression_array.bit_data };

            } else if constexpr (opt == (output_t::tree | output_t::optimized_compression)) {
                constexpr auto compression_array = create_compressor<metadata.byte_count, metadata.node_count, true, get_tree()>(std::span{ metadata.weights });
                constexpr auto get_optimized_tree = [metadata, compression_array]() { return optimize_tree_size<metadata.node_count, skip_parent_nodes>(compression_array.data.node_buffer); };
                return std::tuple{ get_optimized_tree(), compression_array.bit_data };

            } else if constexpr (opt == (output_t::tree | output_t::compress_to_bytes)) {
                constexpr auto compression_array = create_compressor<metadata.byte_count, metadata.node_count, true, get_tree()>(std::span{ metadata.weights });
                constexpr auto get_optimized_tree = [metadata, compression_array]() { return optimize_tree_size<metadata.node_count, skip_parent_nodes>(compression_array.data.node_buffer); };
                constexpr auto get_size = [compression_array, metadata]() { return get_conversion_size(compression_array.bit_data, std::span{ metadata.weights }); };
                constexpr auto get_compressed_bytes = [get_size, compression_array, metadata]() { return compress_input_bytes<get_size(), !skip_check>(
                    compression_array.bit_data, compression_array.data.node_buffer, metadata.input); };
                return std::tuple{ get_optimized_tree(), get_compressed_bytes() };

            } else if constexpr (opt == (output_t::tree | output_t::compress_to_bits)) {
                constexpr auto compression_array = create_compressor<metadata.byte_count, metadata.node_count, true, get_tree()>(std::span{ metadata.weights });
                constexpr auto get_optimized_tree = [metadata, compression_array]() { return optimize_tree_size<metadata.node_count, skip_parent_nodes>(compression_array.data.node_buffer); };
                constexpr auto get_size = [compression_array, metadata]() { return get_conversion_size(compression_array.bit_data, std::span{ metadata.weights }); };
                constexpr auto get_compressed_bits = [get_size, compression_array, metadata]() { return compress_input_bits<get_size(), !skip_check>(
                    compression_array.bit_data, compression_array.data.node_buffer, metadata.input); };
                return std::tuple{ get_optimized_tree(), get_compressed_bits() };

            } else if constexpr (opt == (output_t::tree | output_t::optimized_compression | output_t::compress_to_bytes)) {
                constexpr auto compression_array = create_compressor<metadata.byte_count, metadata.node_count, true, get_tree()>(std::span{ metadata.weights });
                constexpr auto get_optimized_tree = [metadata, compression_array]() { return optimize_tree_size<metadata.node_count, skip_parent_nodes>(compression_array.data.node_buffer); };
                constexpr auto get_size = [compression_array, metadata]() { return get_conversion_size(compression_array.bit_data, std::span{ metadata.weights }); };
                constexpr auto get_compressed_bytes = [get_size, compression_array, metadata]() { return compress_input_bytes<get_size(), !skip_check>(
                    compression_array.bit_data, compression_array.data.node_buffer, metadata.input); };
                return std::tuple{ get_optimized_tree(), compression_array.bit_data, get_compressed_bytes() };

            } else if constexpr (opt == (output_t::tree | output_t::optimized_compression | output_t::compress_to_bits)) {
                constexpr auto compression_array = create_compressor<metadata.byte_count, metadata.node_count, true, get_tree()>(std::span{ metadata.weights });
                constexpr auto get_optimized_tree = [metadata, compression_array]() { return optimize_tree_size<metadata.node_count, skip_parent_nodes>(compression_array.data.node_buffer); };
                constexpr auto get_size = [compression_array, metadata]() { return get_conversion_size(compression_array.bit_data, std::span{ metadata.weights }); };
                constexpr auto get_compressed_bits = [get_size, compression_array, metadata]() { return compress_input_bits<get_size(), !skip_check>(
                    compression_array.bit_data, compression_array.data.node_buffer, metadata.input); };
                return std::tuple{ get_optimized_tree(), compression_array.bit_data, get_compressed_bits() };
            }
        }
    }
}

#endif