#include <iostream>
#include <string>
#include <fstream>
#include <boost/algorithm/string/join.hpp>
#include <iterator>
#include <vector>
#include <numeric>

#include "marisa.h"
#include "vp4.h"
#include "SIMLD/Range.hpp"


std::chrono::duration<double> load_dictionary(const std::string dictionary_file_name, marisa::Trie& dictionary)
{
    auto start = std::chrono::high_resolution_clock::now();

    dictionary.load(dictionary_file_name.c_str());

    return std::chrono::high_resolution_clock::now() - start;
}

std::chrono::duration<double> reconstruct_triples(const std::string index_file_name, std::vector<SIMLD::Triple>& reconstructed_triples)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream inStream(index_file_name.c_str(), std::ios::binary);

    // get total number of triples
    uint_fast64_t no_triples;
    inStream.read(reinterpret_cast<char*>(&no_triples), sizeof(no_triples));

    // get total number of triple bytes
    uint64_t no_of_bytes;
    inStream.read(reinterpret_cast<char*>(&no_of_bytes), sizeof(no_of_bytes));
  
    // read & decompress triple deltas from file
    char* encoded_deltas = new char[no_of_bytes];
    inStream.read(encoded_deltas, no_of_bytes);

    uint_fast32_t *decoded_deltas = new uint_fast32_t[(3*no_triples)];
    p4ndec256v32((unsigned char*)encoded_deltas, 3*no_triples, decoded_deltas);

    // reconstruct all triples
    reconstructed_triples.reserve(no_triples);

    for (uint_fast64_t i = 0; i<(3*(no_triples-1)); )
    { 
        SIMLD::Triple t = SIMLD::Triple(0,0,0);
        t.INTERLACED_BITS_LSB = decoded_deltas[i];
        t.INTERLACED_BITS_NSB = decoded_deltas[i + 1];
        t.INTERLACED_BITS_MSB = decoded_deltas[i + 2];

        reconstructed_triples.emplace_back(t);
        i = i + 3;
    }

    std::partial_sum(reconstructed_triples.begin(), reconstructed_triples.end(), reconstructed_triples.begin());

    return std::chrono::high_resolution_clock::now() - start;
}

std::chrono::duration<double> prepare_query(const std::string& triplePattern, const marisa::Trie& trie, SIMLD::Range& queryBox)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::istringstream iss(triplePattern); 
    std::vector<std::string> queryComponents{ std::istream_iterator<std::string>(iss), {} };
    
    std::string wildcard = "?";
    std::string subject = queryComponents[0];
    std::string predicate = queryComponents[1];
    std::string object = boost::algorithm::join(std::vector<std::string> {queryComponents.begin() + 2, queryComponents.end()}, " ");

    uint_fast32_t s_lo = 0, p_lo = 0, o_lo = 0;
    uint_fast32_t s_hi = UINT32_MAX, p_hi = UINT32_MAX, o_hi = UINT32_MAX;

    marisa::Agent agent;

    if (wildcard.compare(subject) != 0) 
    {
        agent.set_query(subject);
        trie.lookup(agent);
        s_lo = (uint_fast32_t) agent.key().id();
        s_hi = (uint_fast32_t) agent.key().id();
    }
    if (wildcard.compare(predicate) != 0) 
    {
        agent.set_query(predicate);
        trie.lookup(agent);
        p_lo = (uint_fast32_t) agent.key().id();
        p_hi = (uint_fast32_t) agent.key().id();
    }
    if (wildcard.compare(object) != 0)
    {
        agent.set_query(object);
        trie.lookup(agent);
        o_lo = (uint_fast32_t) agent.key().id();
        o_hi = (uint_fast32_t) agent.key().id();
    }

    SIMLD::Triple lowerBound = SIMLD::Triple(s_lo, p_lo, o_lo);
    SIMLD::Triple upperBound = SIMLD::Triple(s_hi, p_hi, o_hi);       
    queryBox = SIMLD::Range(lowerBound, upperBound);

    return std::chrono::high_resolution_clock::now() - start;
}

void dictionary_decoder(std::vector<SIMLD::Triple>::const_iterator start, const std::vector<SIMLD::Triple>::const_iterator end, const marisa::Trie& trie, std::vector<std::string>& result)
{
    marisa::Agent agent;
    for (auto it = start; it != end; ++it) 
    {
        std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t> dict_encoded_triple = (*it).decode();

        agent.set_query(std::get<0>(dict_encoded_triple));
        trie.reverse_lookup(agent);
        std::string subject = std::string(agent.key().str());

        agent.set_query(std::get<1>(dict_encoded_triple));
        trie.reverse_lookup(agent);
        std::string predicate = std::string(agent.key().str());

        agent.set_query(std::get<2>(dict_encoded_triple));
        trie.reverse_lookup(agent);
        std::string object = std::string(agent.key().str());

        std::stringstream mySS;
        mySS << subject << " " << predicate << " " << object << "\n";

        result.emplace_back(mySS.str());

    }
}

std::chrono::duration<double> dictionary_decode(const std::vector<SIMLD::Triple>& solutions, const marisa::Trie& trie, std::vector<std::string>& result)
{
    auto start = std::chrono::high_resolution_clock::now();

    uint64_t length = std::distance(solutions.begin(), solutions.end());
    uint64_t offset = ceil(length / 10);

    /*
     *  @todo: implement this using parFor!
     */
    std::vector<std::string> result_0, result_1, result_2, result_3, result_4, result_5, result_6, result_7, result_8, result_9;
    result_0.reserve(offset);
    result_1.reserve(offset);
    result_2.reserve(offset);
    result_3.reserve(offset);
    result_4.reserve(offset);
    result_5.reserve(offset);
    result_6.reserve(offset);
    result_7.reserve(offset);
    result_8.reserve(offset);
    result_9.reserve(offset);


    tbb::parallel_invoke(
        [&]{ dictionary_decoder(solutions.begin(),              solutions.begin() + 1 * offset, trie, result_0); },
        [&]{ dictionary_decoder(solutions.begin() + 1 * offset, solutions.begin() + 2 * offset, trie, result_1); },
        [&]{ dictionary_decoder(solutions.begin() + 2 * offset, solutions.begin() + 3 * offset, trie, result_2); },
        [&]{ dictionary_decoder(solutions.begin() + 3 * offset, solutions.begin() + 4 * offset, trie, result_3); },
        [&]{ dictionary_decoder(solutions.begin() + 4 * offset, solutions.begin() + 5 * offset, trie, result_4); },
        [&]{ dictionary_decoder(solutions.begin() + 5 * offset, solutions.begin() + 6 * offset, trie, result_5); },
        [&]{ dictionary_decoder(solutions.begin() + 6 * offset, solutions.begin() + 7 * offset, trie, result_6); },
        [&]{ dictionary_decoder(solutions.begin() + 7 * offset, solutions.begin() + 8 * offset, trie, result_7); },
        [&]{ dictionary_decoder(solutions.begin() + 8 * offset, solutions.begin() + 9 * offset, trie, result_8); },
        [&]{ dictionary_decoder(solutions.begin() + 9 * offset, solutions.end(),                trie, result_9); }
    );

    result.insert(std::end(result), std::begin(result_0), std::end(result_0));
    result.insert(std::end(result), std::begin(result_1), std::end(result_1));
    result.insert(std::end(result), std::begin(result_2), std::end(result_2));
    result.insert(std::end(result), std::begin(result_3), std::end(result_3));
    result.insert(std::end(result), std::begin(result_4), std::end(result_4));
    result.insert(std::end(result), std::begin(result_5), std::end(result_5));
    result.insert(std::end(result), std::begin(result_6), std::end(result_6));
    result.insert(std::end(result), std::begin(result_7), std::end(result_7));
    result.insert(std::end(result), std::begin(result_8), std::end(result_8));
    result.insert(std::end(result), std::begin(result_9), std::end(result_9));


   return std::chrono::high_resolution_clock::now() - start;
}

std::chrono::duration<double> print_results(const std::vector<std::string>& results)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::copy(results.begin(), results.end(), std::ostream_iterator<std::string>(std::cout));
  
    return std::chrono::high_resolution_clock::now() - start;
}



int main(int argc, char **argv) 
{
    if (argc != 3) {
        std::cerr << "Usage: simldSearch example.simld.index example.simld.dict\n";
        return 1;
    }

    std::string index_file_name = std::string(argv[1]);
    std::string dictionary_file_name = std::string(argv[2]);

    // load dictionary
    marisa::Trie dictionary; 
    std::chrono::duration<double> dict_load_time = load_dictionary(dictionary_file_name, dictionary);
    std::cout << "SIMLD: Loaded term dictionary in " << dict_load_time.count() << " s.\n";


    std::vector<SIMLD::Triple> triples;
    std::chrono::duration<double> triple_load_time = reconstruct_triples(index_file_name, triples);
    std::cout << "SIMLD: Loaded " << triples.size() << " compressed triple index in " << triple_load_time.count() << " s.\n";
    std::cout << "SIMLD: morton_min " <<triples.front() << "\n";
    std::cout << "SIMLD: morton_max " <<triples.back() << "\n";

    std::cout << "\n" << "SIMLD: Ready for query answering!" << "\n";
    std::cout << "Please type Triple Search Pattern, using '?' for wildcards. e.g ? http://www.w3.org/1999/02/22-rdf-syntax-ns#type ?" << "\n\n";


    std::string triplePattern;
    marisa::Agent agent;
    SIMLD::Range queryBox = SIMLD::Range(SIMLD::Triple(0,0,0), SIMLD::Triple(0,0,0));

    while(std::getline(std::cin, triplePattern))
    {   
        std::chrono::duration<double> query_prep_time = prepare_query(triplePattern, dictionary, queryBox);

        uint_fast64_t cardinality_estimate = 0;
        std::chrono::duration<double> card_estim_time = queryBox.estimate_cardinality_avx2(triples.begin(), triples.end(), cardinality_estimate);    
       
        std::vector<SIMLD::Triple> solutions;
        solutions.reserve(cardinality_estimate);
        std::chrono::duration<double> query_exec_time = queryBox.range_lookup_avx2(triples, solutions);
        

        std::vector<std::string> results;
        results.reserve(solutions.size());        
        std::chrono::duration<double> dict_decode_time = dictionary_decode(solutions, dictionary, results);
        std::chrono::duration<double> output_time = print_results(results);

        
        std::cout << "SIMLD: Prepared query in " << query_prep_time.count() << " s\n";
        std::cout << "SIMLD: Estimated carinality in " << card_estim_time.count() << " s\n";
        std::cout << "SIMLD: Expecting " << cardinality_estimate << " triples. \n";

        std::cout << "SIMLD: Executed AVX2 query in " << query_exec_time.count() << " s\n";
        std::cout << "SIMLD: Found " << solutions.size() << " triples." << "\n";

        std::cout << "SIMLD: Dict-decoded triples in " << dict_decode_time.count() << " s\n";
        std::cout << "SIMLD: Serialized triples in " << output_time.count() << " s\n";
        std::cout << "SIMLD: " << (query_prep_time + query_exec_time + dict_decode_time + output_time).count() << " s\n";

        solutions.clear();
    }

    return 0;
}