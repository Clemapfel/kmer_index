// [SEQ] fm search
template<seqan3::alphabet alphabet_t, typename index_t>
static void fm_search(benchmark::State& state, size_t query_length,const index_t* index)
{
    state.counters["seed"] = seed;

    input_generator<alphabet_t> input(seed);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(search(input.generate_sequence(query_length), *index));
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    seed++;
}

// [SEQ] multi kmer search
template<seqan3::alphabet alphabet_t, size_t... ks>
static void multi_kmer_search(
        benchmark::State& state,
        size_t query_length,
        const kmer::kmer_index<alphabet_t, uint32_t, ks...>* index)
{
    state.counters["seed"] = seed;
    input_generator<alphabet_t> input(seed);

    bool error_occurred = false;
    try
    {
        size_t i = 0;
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(index->search(input.generate_sequence(query_length)));
        }
    }
    catch (...)
    {
        error_occurred = true;
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    std::vector<size_t> _all_ks = {ks...};
    size_t j = 0;
    for (auto k : _all_ks)
    {
        state.counters["k_" + std::to_string(j)] = k;
        j++;
    }
    state.counters["valid"] = not error_occurred;

    seed++;
}

