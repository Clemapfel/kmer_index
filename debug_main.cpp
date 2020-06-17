//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>
#include <fast_pow.hpp>

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
#include <thread_pool.hpp>

//#include <benchmarks/benchmarks.hpp>

using namespace seqan3;
using alphabet_t = dna4;
constexpr size_t k = 5;
constexpr size_t text_size = 1000000;

void force_error(std::vector<alphabet_t> query, float seed)
{
    auto input = input_generator<alphabet_t>(seed);
    auto text = input.generate_text(text_size, {});
    auto kmer = kmer::make_kmer_index<k>(text);
    auto fm = fm_index(text);

    std::vector<unsigned int> fm_result;
    for (auto _ : search(query, fm))
        fm_result.push_back(_.second);

    auto res = kmer.search<k>(query);
    auto kmer_result = res.to_vector();

    auto equal = fm_result == kmer_result;

    seqan3::debug_stream << (equal ? "EQUAL FOR " : "NOT EQUAL FOR ") << "\nQUERY " << query << " (" << query.size() << ")\n"
                         << "\nFM : " << fm_result << "\n\n KMER : " << kmer_result << "\n";
    seqan3::debug_stream << "query size = " << query.size() << "\nseed = " << seed << "\n"
                         <<  "difference (fm - kmer) = " << int(fm_result.size()) - int(kmer_result.size()) << "\n";

    std::vector<uint32_t> diff;
    std::set_difference(fm_result.begin(), fm_result.end(), kmer_result.begin(), kmer_result.end(), std::inserter(diff, diff.begin()));
    seqan3::debug_stream << diff;

    if (not equal)
        exit(1);
    else
        exit(0);
}

template<size_t... ks>
size_t pick_right_k(size_t query_size)
{
    auto _all_ks = std::vector<size_t>{ks...};
    std::sort(_all_ks.begin(), _all_ks.end(), [](size_t a, size_t b) -> bool {return a > b;});

    size_t optimal_k = _all_ks.at(0);

    if (_all_ks.size() > 1)
    {
        size_t optimal_k_i = 0;

        for (optimal_k_i; optimal_k_i < _all_ks.size(); ++optimal_k_i)
        {
            // c.f. addendum
            size_t k = _all_ks.at(optimal_k_i);
            if (query_size % k == 0)
            {
                optimal_k = k;
                break;
            }
            else if (query_size == k-1 or query_size == k-2)
            {
                optimal_k = k;
                break;
            }
            else if (query_size % k <= query_size % optimal_k) // i % n = i & n-1
            {
                optimal_k = k;
            }
        }
    }

    return optimal_k;
}

int main()
{

    for (size_t q = 3; q < 40; ++q)
    {
        seqan3::debug_stream << "query " << q << " : picked " << pick_right_k<5, 7, 11>(q) << "\n";
    }

    /*
    size_t i = 0;
    auto text = "ACGTCGT"_dna4;
    for (auto hash : text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}}))
    {
        seqan3::debug_stream << i << " | " << hash << "\n";
        i++;
    }
    /*
    for (float seed = 2100.5; seed < 2200; seed++)
    {
        auto input = input_generator<dna4>(seed);
        auto text = input.generate_text(text_size, {});
        auto kmer = kmer::make_kmer_index<k>(text);
        auto fm = fm_index(text);

        //for (size_t query_size = 2*k; query_size < 3*k+1; ++query_size)
        size_t query_size = 18;
        {
            auto query = input.generate_sequence(query_size);

            std::vector<unsigned int> fm_result;
            for (auto _ : search(query, fm))
                fm_result.push_back(_.second);

            auto kmer_result = kmer.search(query).to_vector();

            auto equal = fm_result == kmer_result;

            seqan3::debug_stream //<< "\nFM : " << fm_result << "\n\n KMER : " << kmer_result << "\n"
                    << "_____________________\n" << (equal ? "" : "NOT ") << "EQUAL FOR QUERY " << query << " (" << query.size() << ")\n";
            seqan3::debug_stream << "seed = " << seed << "\n"
                                 << "difference (fm - kmer) = " << int(fm_result.size()) - int(kmer_result.size())
                                 << "\n"
                                 << "k = " << k << "\n"
                                 << "text_size = " << text_size << "\n";

            std::vector<uint32_t> diff;
            std::set_difference(fm_result.begin(), fm_result.end(), kmer_result.begin(), kmer_result.end(),
                                std::inserter(diff, diff.begin()));
            seqan3::debug_stream << diff;
        }
    }


    return 0;

    /*
    auto query = "ACGTAACGTA"_dna4;



    auto result_test = kmer.search_test(query).to_vector();
    //auto result_true = kmer.search(query).to_vector();

    //seqan3::debug_stream << result_true << "\n\n";
    seqan3::debug_stream << result_test << "\n\n";

    //query = "ACGTAACGTAAC"_dna4;

    auto fm = fm_index(text);
    std::vector<size_t> fm_result;
    for (auto _ : search(query, fm))
        fm_result.push_back(_.second);

    seqan3::debug_stream << fm_result << "\n";
    //seqan3::debug_stream << result_test << "\n";

    //seqan3::debug_stream << "kmer_old : " << result_true.size() << "\nkmer_new : " << result_test.size() << "\nfm : " << fm_result.size();
    */
}

