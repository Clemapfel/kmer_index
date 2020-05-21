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
#include <thread_pool.hpp>

    using namespace seqan3;
    using alphabet_t = dna4;
    constexpr size_t k = 6;
    constexpr size_t text_size = 100;

    void force_error(std::vector<dna4> query, size_t seed)
    {
        auto input = input_generator<dna4>(seed);
        auto text = input.generate_sequence(text_size);
        auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);

        debug_stream << "seed : " << seed << "\n";
        debug_stream << "text : " << text << " (" << text.size() << ")\n\n";
        debug_stream << "query : " << query << "\n";

        auto fm = fm_index(text);
        auto temp = search(query, fm);
        std::vector<uint32_t> fm_res;

        for (auto pos : temp)
            fm_res.push_back(pos.second);

        debug_stream << "fm  : " << fm_res << "\n";

        auto kmer_results = kmer.search(query);
        auto kmer_res = kmer_results.to_vector();

        debug_stream << "kmer : " << kmer_res << "\n\n";

        int dif = fm_res.size() - kmer_res.size();
        debug_stream << "difference : " << dif << "\n";

        if (dif != 0)
            exit(1);

        //exit(0);
    }

    using kmer_index_t = kmer_index_element<seqan3::dna4, k, uint32_t>;

    std::vector<size_t> generate_all_hashs(std::vector<alphabet_t> suffix)
    {
        auto sigma = seqan3::alphabet_size<alphabet_t>;
        auto m = suffix.size();

        size_t suffix_hash = 0;
        for (size_t i = k - m, index = 0; i < k; ++i, ++index)
            suffix_hash += seqan3::to_rank(suffix.at(index)) * std::pow(sigma, k - i - 1);

        size_t lower_bound = 0 + suffix_hash;
        size_t upper_bound = std::pow(sigma, k) - std::pow(sigma, m) + suffix_hash;    // [1]

        size_t step_size = std::pow(sigma, k - (k - m - 1) - 1);

        std::vector<size_t> output;

        for (size_t i = lower_bound; i <= upper_bound; i += step_size)
            output.push_back(i);

        return output;
    }

    int main()
    {
        //force_error("TCC"_dna4, 15);

        debug_stream << "starting test...\n";

        auto input = input_generator<dna4>();

        for (size_t i = 0; i < 1000; ++i)
        {
            // state not reset so new text everytime
            input.reset_state(i);
            auto text = input.generate_sequence(text_size);
            auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);
            auto fm = fm_index(text);

            for (size_t size : {k-3, k-2, k, 2*k, 3*k, 3*k+k-2})
            {
                //auto query = "GGCAGCATCT"_dna4;
                auto query = input.generate_sequence(size);

                auto fm_range = search(query, fm);
                std::vector<uint32_t> fm_vec;
                for (auto p : fm_range)
                    fm_vec.push_back(p.second);

                auto kmer_range = kmer.search(query);
                auto kmer_vec = kmer_range.to_vector();

                if (kmer_vec != fm_vec)
                {
                    debug_stream << "seed : " << i << "\n";
                    debug_stream << "text : " << std::vector<dna4>(text.begin(), text.end()) << "\n\n";
                    debug_stream << "query : " << query << " (" << query.size() << ")\n";
                    debug_stream << "fm  : " << fm_vec << " (" << fm_vec.size() << ") " << "\n";
                    debug_stream << "kmer : " << kmer_vec << " (" << kmer_vec.size() << ")" << "\n";
                    return 1;
                }


            }

        }

        debug_stream << "test passed succesfully.";
    }




