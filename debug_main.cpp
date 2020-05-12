//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

#include <thread>
#include <chrono>
#include <memory>
#include <future>
#include <functional>
#include <benchmarks/input_generator.hpp>
#include <thread_pool.hpp>
#include <kmer_index_result.hpp>

using namespace seqan3;
constexpr size_t k = 7;

int main()
{
    /*
    input_generator<seqan3::dna4> input{0001};

    std::vector<seqan3::dna4> q = input.generate_sequence(2 * k);

    auto text = input.generate_text(1000000, {q});
    auto kmer = kmer_index<seqan3::dna4, k>{text};
    auto fm = seqan3::fm_index{text};

    auto kmer_results = kmer.search(q);

    seqan3::debug_stream << "query : " << q << "\n"
                         << "kmer : " << kmer_results.to_vector() << "\n"
                         << "fm : " << seqan3::search(q, fm) << "\n";

    return 0;*/

    input_generator<seqan3::dna4> input;

    std::vector<std::vector<seqan3::dna4>> k_queries = input.generate_queries(1000, k);
    auto subk_queries = input.generate_queries(1000, k-2);
    auto nk_queries = input.generate_queries(1000, 4*k);

    auto text = input.generate_text(100000, nk_queries);
    auto kmer = kmer_index<seqan3::dna4, k>{text};
    auto fm = seqan3::fm_index{text};

    size_t i = 0;
    for (auto q : nk_queries)
    {
        seqan3::debug_stream << i << "\n";
        auto kmer_results = kmer.search(q);

        seqan3::debug_stream << "query : " << q << "\n"
                            << "kmer : " << kmer_results.to_vector() << "\n"
                            << "fm : " << seqan3::search(q, fm) << "\n";
        i++;
    }

    seqan3::debug_stream << "done.\n";
    return 0;
}

/*
 // search first n*k parts
                size_t rest_n = query.size() % k;
                std::vector<const std::vector<position_t>*> positions{};
                positions.reserve(query.size() / k + 1);

                for (auto it = query.begin(); it + k != query.end() - rest_n; it += k)
                    positions.push_back(&_data.find(hash(it))->second);


                // search last m < k part
                if (rest_n != 0)
                    for (const auto* r : search_subk(query.end() - rest_n, rest_n))
                        positions.push_back(r);


auto result = result_t(positions.at(0));

size_t i = -1;
for (auto &start_pos : *(positions.at(0)))
{
i++;
position_t previous_pos = start_pos;

for (size_t j = 1; j <= positions.size(); ++j)
{
if (j == positions.size())  // successfull hit
{
seqan3::debug_stream << "do use : " << i <<"\n";
break;
}

const auto* current = positions.at(j);
if (std::find(current->begin(), current->end(), previous_pos + k) != current->end())
{
previous_pos += k;
}
else
{
seqan3::debug_stream << "don't use : " << i <<"\n";
result.should_not_use(i);
break;
}
}
}

return result;
*/



