#include <immintrin.h>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>
#include <boost/algorithm/string/join.hpp>
#include <vector>
#include <tuple>
#include <boost/unordered_set.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>
#include <cinttypes>
#include <numeric>
#include <chrono> 
#include <filesystem>
#include <clocale>
#include <locale>

#include "vp4.h"
#include "marisa.h"
#include "SIMLD/Triple.hpp"

typedef boost::tuple<std::string, std::string, std::string> triple_t;

namespace boost::tuples {
    size_t hash_value(tuple<std::string, std::string, std::string> const & t) 
    {
        std::string subject = get<0>(t);
        std::string predicate = get<1>(t);
        std::string object = get<2>(t);    

        return std::hash<std::string>{}(subject.append(predicate.append(object))) ;
    }
}

std::chrono::duration<double> parse_ntriples_document(const std::string& triple_file_path, boost::unordered_set<boost::tuple<std::string, std::string, std::string>>& triples)
{
    auto now = std::chrono::high_resolution_clock::now();

    std::ifstream triples_file(triple_file_path);
    if(!triples_file.is_open()) throw std::runtime_error("Could not open triples file");

    std::string line;
    while(std::getline(triples_file, line))
    {
        std::istringstream iss(line); 
        std::vector<std::string> terms{ std::istream_iterator<std::string>(iss), {} };
        std::string subject = terms[0];
        std::string predicate = terms[1];
        std::string object = boost::algorithm::join(std::vector<std::string> {terms.begin() + 2, terms.end() - 1}, " ");

        triples.insert(boost::make_tuple(subject, predicate, object));
    }

    triples_file.close();

    return std::chrono::high_resolution_clock::now() - now;
}

std::chrono::duration<double> compress_and_store_term_dictionary(marisa::Trie& dictionary, const boost::unordered_set<boost::tuple<std::string, std::string, std::string>>& unique_triples, const std::string& triple_file_path)
{
    auto now = std::chrono::high_resolution_clock::now();

    // extract unique RDF terms
    boost::unordered_set<std::string> unique_rdf_terms;
    BOOST_FOREACH(const triple_t& triple, unique_triples)
    {
        unique_rdf_terms.insert(boost::get<0>(triple));
        unique_rdf_terms.insert(boost::get<1>(triple));
        unique_rdf_terms.insert(boost::get<2>(triple));

    }  

    // build MARISA trie as term-id-dictionary
    marisa::Keyset keyset;
    BOOST_FOREACH(const std::string& unique_term, unique_rdf_terms)
    {
        keyset.push_back(unique_term);
    }
    
    dictionary.build(keyset, MARISA_MAX_NUM_TRIES | MARISA_TEXT_TAIL | MARISA_LABEL_ORDER | MARISA_TINY_CACHE);
   
    std::string dict_file_path = std::filesystem::current_path().string();
    dict_file_path.append("/");
    dict_file_path.append(std::filesystem::path(triple_file_path).filename());
    dict_file_path.append(".simld");
    dictionary.save(dict_file_path.c_str());

    std::cout << "SIMLD: Compressed " << unique_rdf_terms.size() << " unique RDF terms into "  << dictionary.total_size() << " bytes.\n";
    
    return std::chrono::high_resolution_clock::now() - now;
}    

std::chrono::duration<double> dict_encode_unique_triples(const marisa::Trie& dictionary, const boost::unordered_set<boost::tuple<std::string, std::string, std::string>>& unique_triples, std::vector<std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t>>& dict_encoded_triples)
{
    auto now = std::chrono::high_resolution_clock::now();

    marisa::Agent agent;
    uint_fast32_t subject_id, predicate_id, object_id;

    BOOST_FOREACH(const triple_t& triple, unique_triples)
    {
        agent.set_query(boost::get<0>(triple));
        dictionary.lookup(agent);
        subject_id = (uint_fast32_t) agent.key().id();

        agent.set_query(boost::get<1>(triple));
        dictionary.lookup(agent);
        predicate_id = (uint_fast32_t) agent.key().id();

        agent.set_query(boost::get<2>(triple));
        dictionary.lookup(agent);
        object_id = (uint_fast32_t) agent.key().id();

        dict_encoded_triples.emplace_back(subject_id, predicate_id, object_id);
    }
    
    return std::chrono::high_resolution_clock::now() - now;
}

std::chrono::duration<double> morton_encode_unique_triples(std::vector<std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t>>& dict_encoded_triples, std::vector<SIMLD::Triple>& morton_encoded_triples)
{
    auto now = std::chrono::high_resolution_clock::now();

    for (auto& t: dict_encoded_triples)
    {
        uint_fast32_t subject_id   = std::get<0>(t);
        uint_fast32_t predicate_id = std::get<1>(t);
        uint_fast32_t object_id    = std::get<2>(t);; 

//        morton_encoded_triples.emplace_back(subject_id, predicate_id, object_id);
        SIMLD::Triple morton_code = SIMLD::Triple(subject_id, predicate_id, object_id);
        morton_encoded_triples.push_back(morton_code);
    }
    std::stable_sort(morton_encoded_triples.begin(), morton_encoded_triples.end());

    return std::chrono::high_resolution_clock::now() - now;
}

std::chrono::duration<double> compress_and_store_triple_index(const std::vector<SIMLD::Triple>& morton_encoded_triples, const std::string& triple_file_path)
{
    auto now = std::chrono::high_resolution_clock::now();

    // differencing
    // TODO: merge the following loops!
    std::vector<SIMLD::Triple> delta_triples;
    delta_triples.reserve(morton_encoded_triples.size());
    delta_triples.emplace_back(morton_encoded_triples.front());

    for (size_t i = 0; i<morton_encoded_triples.size()-1; i++)
    {
        SIMLD::Triple diff = morton_encoded_triples[i+1] - morton_encoded_triples[i];
        delta_triples.emplace_back(diff);
    }

    std::vector<uint_fast32_t> encoder_input;
    encoder_input.reserve(3*delta_triples.size());

    for(std::vector<SIMLD::Triple>::iterator it = delta_triples.begin(); it != delta_triples.end(); ++it) 
    {
        encoder_input.emplace_back(it->INTERLACED_BITS_LSB);
        encoder_input.emplace_back(it->INTERLACED_BITS_NSB);
        encoder_input.emplace_back(it->INTERLACED_BITS_MSB);
    }


    // differential coding and compression
    unsigned char *outBuffer = new unsigned char[encoder_input.size()*sizeof(uint32_t)];
    size_t no_of_bytes = p4nenc256v32(&encoder_input[0], encoder_input.size(), outBuffer);

    // serialization
    std::string index_file_path = std::filesystem::current_path().string();
    index_file_path.append("/");
    index_file_path.append(std::filesystem::path(triple_file_path).filename());
    index_file_path.append(".index.simld");

    std::ofstream p4nenc256v32_file(index_file_path, std::ios::out | std::ios::binary );

    // total number of triples
    uint_fast64_t no_triples = delta_triples.size();
    p4nenc256v32_file.write(reinterpret_cast<char const*>(&no_triples), sizeof(no_triples));

    // total buffer size
    p4nenc256v32_file.write(reinterpret_cast<char const*>(&no_of_bytes), sizeof(no_of_bytes));

    // compressed triples
    p4nenc256v32_file.write(reinterpret_cast<char*>(outBuffer), no_of_bytes);

    p4nenc256v32_file.close();

    std::cout << "SIMLD: Compressed " << no_triples << " triples into "  << no_of_bytes << " bytes.\n";
    std::cout << "SIMLD:  " << index_file_path << "\n";

    return std::chrono::high_resolution_clock::now() - now;
}

int main(int argc, char **argv) 
{
    if (argc != 2) {
        std::cerr << "Usage: ./simldBuild your_triple_file.nt\n";
        return 1;
    }

    std::string triple_file_path = std::string(argv[1]);
    
    std::cout << "SIMLD: Building SIMLD index for" << triple_file_path << ".\n";

    // @TODO: use https://github.com/dice-group/rdf-parser
    boost::unordered_set<boost::tuple<std::string, std::string, std::string>> unique_triples;
    std::chrono::duration<double> triples_parse_time = parse_ntriples_document(triple_file_path, unique_triples);
    std::cout << "SIMLD: Parsed " << unique_triples.size() << " triples in "  << triples_parse_time.count() << " s.\n";  

    marisa::Trie dictionary;

    std::chrono::duration<double> dict_build_time = compress_and_store_term_dictionary(dictionary, unique_triples, triple_file_path);
    std::cout << "SIMLD: Constructed term dictionary in " << dict_build_time.count() << " s.\n";  

    std::vector<std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t>> dict_encoded_triples;
    dict_encoded_triples.reserve(unique_triples.size());
    std::chrono::duration<double> dict_encoding_time = dict_encode_unique_triples(dictionary, unique_triples, dict_encoded_triples);
    std::cout << "SIMLD: Dictionary encoding of " << dict_encoded_triples.size() << " triples took " << dict_encoding_time.count() << " s.\n";

    std::vector<SIMLD::Triple> morton_encoded_triples;
    dict_encoded_triples.reserve(dict_encoded_triples.size());
    std::chrono::duration<double> morton_encoding_time = morton_encode_unique_triples(dict_encoded_triples, morton_encoded_triples);
    std::cout << "SIMLD: Morton encoding of " << morton_encoded_triples.size() << " triples took " << morton_encoding_time.count() << " s.\n";

    std::chrono::duration<double> index_build_time = compress_and_store_triple_index(morton_encoded_triples, triple_file_path);
    std::cout << "SIMLD: Triple index compression took " << index_build_time.count() << " s.\n";

    return 0;
}