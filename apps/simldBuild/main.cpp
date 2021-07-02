#include <immintrin.h>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <vector>
#include <tuple>
#include <cinttypes>
#include <numeric>
#include <chrono> 

#include "vp4.h"
#include "marisa.h"
#include "SIMLD/Triple.hpp"

std::chrono::duration<double> build_term_dictionary(const std::string triple_file_name, const std::string dictionary_file_name, std::vector<std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t>>& dict_encoded_triples)
{
    auto start = std::chrono::high_resolution_clock::now();

    marisa::Trie dictionary;
    dictionary.load(dictionary_file_name.c_str());

    std::ifstream triples_file(triple_file_name);
    if(!triples_file.is_open()) throw std::runtime_error("Could not open triples file");

    std::string line;
    while(std::getline(triples_file, line))
	{
		std::stringstream ss(line);
		uint_fast32_t subject_id, predicate_id, object_id;	
        ss >> subject_id >> predicate_id >> object_id;
        dict_encoded_triples.push_back(std::make_tuple(subject_id, predicate_id, object_id));
    }

    std::cout << "SIMLD: Compressed " << dictionary.num_keys() << " unique RDF terms into "  << dictionary.io_size() << " bytes.\n";

    return std::chrono::high_resolution_clock::now() - start;
}

std::chrono::duration<double> build_triple_index(const std::vector<std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t>>& dict_encoded_triples)
{
    auto start = std::chrono::high_resolution_clock::now();

    // morton_encoding
    std::vector<SIMLD::Triple> morton_encoded_triples;

    for (auto t: dict_encoded_triples)
    {
        uint_fast32_t subject_id   = std::get<0>(t);
        uint_fast32_t predicate_id = std::get<1>(t);
        uint_fast32_t object_id    = std::get<2>(t);; 

        morton_encoded_triples.push_back(SIMLD::Triple(subject_id, predicate_id, object_id));
    }
    std::stable_sort(morton_encoded_triples.begin(), morton_encoded_triples.end());

    // differencing
    // TODO: merge the following loops!
    std::vector<SIMLD::Triple> delta_triples;
    delta_triples.push_back(morton_encoded_triples.front());

    for (size_t i = 0; i<morton_encoded_triples.size()-1; i++)
    {
        SIMLD::Triple diff = morton_encoded_triples[i+1] - morton_encoded_triples[i];
        delta_triples.push_back(diff);
    }

    std::vector<uint_fast32_t> encoder_input;
    for(std::vector<SIMLD::Triple>::iterator it = delta_triples.begin(); it != delta_triples.end(); ++it) 
    {
        encoder_input.push_back(it->INTERLACED_BITS_LSB);
        encoder_input.push_back(it->INTERLACED_BITS_NSB);
        encoder_input.push_back(it->INTERLACED_BITS_MSB);
    }


    // differential coding and compression
    unsigned char *outBuffer = new unsigned char[encoder_input.size()*sizeof(uint32_t)];
    size_t no_of_bytes = p4nenc256v32(&encoder_input[0], encoder_input.size(), outBuffer);

    // serialization
    std::ofstream p4nenc256v32_file("eswc.simld.index", std::ios::out | std::ios::binary );

    // total number of triples
    uint_fast64_t no_triples = delta_triples.size();
    p4nenc256v32_file.write(reinterpret_cast<char const*>(&no_triples), sizeof(no_triples));

    // total buffer size
    p4nenc256v32_file.write(reinterpret_cast<char const*>(&no_of_bytes), sizeof(no_of_bytes));

    // compressed triples
    p4nenc256v32_file.write(reinterpret_cast<char*>(outBuffer), no_of_bytes);

    p4nenc256v32_file.close();

    std::cout << "SIMLD: Compressed " << no_triples << " triples into "  << no_of_bytes << " bytes.\n";

    return std::chrono::high_resolution_clock::now() - start;
}

int main(int argc, char **argv) 
{
    if (argc != 3) {
        std::cerr << "Usage: dict_encoded triple infile.tsv\n";
        return 1;
    }

    std::string dictionary_file_name = std::string(argv[1]);
    std::string triple_file_name = std::string(argv[2]);

    std::vector<std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t>> dict_encoded_triples;

    std::chrono::duration<double> dict_build_time = build_term_dictionary(triple_file_name, dictionary_file_name, dict_encoded_triples);
    std::cout << "SIMLD: Constructed term dictionary in " << dict_build_time.count() << " s.\n";

    std::chrono::duration<double> index_build_time = build_triple_index(dict_encoded_triples);  
    std::cout << "SIMLD: Constructed compressed index in " << index_build_time.count() << " s.\n";

    return 0;
}

