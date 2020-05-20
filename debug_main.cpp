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

    void force_error(std::vector<dna4> query, size_t seed)
    {
        auto input = input_generator<dna4>(seed);
        auto text = input.generate_sequence(1e3);
        text.insert(text.begin(), query.begin(), query.end());
        auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);

        debug_stream << "text : " << text << "\n\n";
        debug_stream << "query : " << query << "\n";

        auto fm = fm_index(text);
        debug_stream << "fm  : " << search(query, fm) << "\n";

        auto kmer_results = kmer.search(query);

        debug_stream << "kmer : " << kmer_results.to_vector() << "\n";

        exit(0);
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

        size_t step_size = std::pow(sigma, k - (k - m -1) -1);

        std::vector<size_t> output;

        for (size_t i = lower_bound; i <= upper_bound; i += step_size)
            output.push_back(i);

        return output;

        /*
         * auto input = input_generator<dna4>();

        input.reset_state();
        auto text = input.generate_sequence(1e3);
        auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);

        auto suffix = "TCG"_dna4;

        auto all_kmers = kmer.get_all_kmer_with_suffix(suffix);
        std::vector<size_t> hashs_true;

        for (auto q : all_kmers)
            hashs_true.push_back(kmer.hash(q.begin()));

        std::sort(hashs_true.begin(), hashs_true.end());

        auto hashs_test = generate_all_hashs(suffix);

        seqan3::debug_stream << all_kmers << "\n";
        seqan3::debug_stream << (hashs_true == hashs_test) << "\n";

        generate_all_hashs("TCG"_dna4);
         */
    }

    int main()
    {
        //force_error("GTTAAG"_dna4, 16);

        debug_stream << "starting test...\n";

        auto input = input_generator<dna4>();

        for (size_t i = 0; i < 1000; ++i)
        {
            // state not reset so new text everytime
            input.reset_state(i);
            auto text = input.generate_sequence(1e3);
            auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);
            auto fm = fm_index(text);

            for (size_t size : {k})
            {
                //auto query = "GGCAGCATCT"_dna4;
                auto query = input.generate_sequence(size);

                auto fm_results = search(query, fm);
                size_t fm_size = 0;
                for (auto f : fm_results)
                    fm_size++;

                bool results_equal = (kmer.search(query).size() - fm_size) == 0;
                auto kmer_results =  kmer.search(query);

                debug_stream << "text (head) : " << std::vector<dna4>(text.begin(), text.begin() + 100) << "\n";
                debug_stream << "query : " << query << " (" << query.size() << ")\n";
                debug_stream << "fm  : " << search(query, fm) << " (" << fm_size << ") " << "\n";
                debug_stream << "kmer : " << kmer_results.to_vector()  << " (" << kmer_results.size() << ")" << "\n";

                if (not results_equal)
                    return 1;
            }

        }

        debug_stream << "test passed succesfully.";
    }




