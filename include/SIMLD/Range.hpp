#include <immintrin.h>
#include <tuple>
#include <vector>
#include <immintrin.h>
#include <math.h>
#include <tbb/parallel_invoke.h>
#include "Triple.hpp"


namespace SIMLD 
{
    struct Range
    {
        Triple lowerBound = Triple(0,0,0);
        Triple upperBound = Triple(UINT32_MAX, UINT32_MAX, UINT32_MAX);

        Range (Triple low, Triple high)
        {
            lowerBound = low;
            upperBound = high;
        }

        inline bool contains(Triple element)
        {
            __m256i min = _mm256_set_epi32(lowerBound.INTERLACED_BITS_MSB, lowerBound.INTERLACED_BITS_NSB, lowerBound.INTERLACED_BITS_LSB, 0x0, 0x0, 0x0, 0x0, 0x0);
            __m256i max = _mm256_set_epi32(upperBound.INTERLACED_BITS_MSB, upperBound.INTERLACED_BITS_NSB, upperBound.INTERLACED_BITS_LSB, UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX);

            __m256i c = _mm256_xor_si256(_mm256_cmpeq_epi64(_mm256_andnot_si256(_mm256_xor_si256(min, max), min), 
                                                            _mm256_andnot_si256(_mm256_xor_si256(min, max), _mm256_set_epi32(element.INTERLACED_BITS_MSB, element.INTERLACED_BITS_NSB, element.INTERLACED_BITS_LSB, 0x0, 0x0, 0x0, 0x0, 0x0))), 
                                         _mm256_set1_epi64x(-1LL));
            
            return _mm256_testz_si256(c, c);
        }
/*
        inline bool contains_avx2(Triple e7, Triple e6, Triple e5, Triple e4, Triple e3, Triple e2, Triple e1, Triple e0)
        {


            __m256i min_2 = _mm256_set_epi32(lowerBound.SUBJECT_ID,   lowerBound.PREDICATE_ID, lowerBound.OBJECT_ID,    lowerBound.SUBJECT_ID,   lowerBound.PREDICATE_ID, lowerBound.OBJECT_ID,    lowerBound.SUBJECT_ID,   lowerBound.PREDICATE_ID);
            __m256i min_1 = _mm256_set_epi32(lowerBound.OBJECT_ID,    lowerBound.SUBJECT_ID,   lowerBound.PREDICATE_ID, lowerBound.OBJECT_ID,    lowerBound.SUBJECT_ID,   lowerBound.PREDICATE_ID, lowerBound.OBJECT_ID,    lowerBound.SUBJECT_ID);
            __m256i min_0 = _mm256_set_epi32(lowerBound.PREDICATE_ID, lowerBound.OBJECT_ID,    lowerBound.SUBJECT_ID,   lowerBound.PREDICATE_ID, lowerBound.OBJECT_ID,    lowerBound.SUBJECT_ID,   lowerBound.PREDICATE_ID, lowerBound.OBJECT_ID);

            __m256i max_2 = _mm256_set_epi32(upperBound.SUBJECT_ID,   upperBound.PREDICATE_ID, upperBound.OBJECT_ID,    upperBound.SUBJECT_ID,   upperBound.PREDICATE_ID, upperBound.OBJECT_ID,    upperBound.SUBJECT_ID,   upperBound.PREDICATE_ID);
            __m256i max_1 = _mm256_set_epi32(upperBound.OBJECT_ID,    upperBound.SUBJECT_ID,   upperBound.PREDICATE_ID, upperBound.OBJECT_ID,    upperBound.SUBJECT_ID,   upperBound.PREDICATE_ID, upperBound.OBJECT_ID,    upperBound.SUBJECT_ID);
            __m256i max_0 = _mm256_set_epi32(upperBound.PREDICATE_ID, upperBound.OBJECT_ID,    upperBound.SUBJECT_ID,   upperBound.PREDICATE_ID, upperBound.OBJECT_ID,    upperBound.SUBJECT_ID,   upperBound.PREDICATE_ID, upperBound.OBJECT_ID);




            __m256i c = _mm256_xor_si256(_mm256_cmpeq_epi64(_mm256_andnot_si256(_mm256_xor_si256(min, max), min), 
                                                            _mm256_andnot_si256(_mm256_xor_si256(min, max), _mm256_set_epi64x(element.INTERLACED_BITS_MSB, element.INTERLACED_BITS_NSB, element.INTERLACED_BITS_LSB, 0x0))), 
                                         _mm256_set1_epi64x(-1LL));
            
            return _mm256_testz_si256(c, c);

        }
*/
    };
    
    inline void filter_by_range(Range& range, std::vector<Triple>::const_iterator start, const std::vector<Triple>::const_iterator end, std::vector<Triple> &solutions)
        {
            for (auto it = start; it != end; ++it) 
            {
                if (range.contains(*it))
                {
                    solutions.push_back(*it); 
                } 
            }
        }

    std::chrono::duration<double> range_lookup_by_filter(Range& range, const std::vector<Triple>& triples, std::vector<Triple> &solutions)
        {
            auto start = std::chrono::high_resolution_clock::now();

            std::vector<Triple> sol0, sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8, sol9;

            uint64_t length = std::distance(triples.begin(), triples.end());
            uint64_t offset = ceil(length / 10);

            tbb::parallel_invoke(
                [&]{ filter_by_range(range, triples.begin(),              triples.begin() + 1 * offset, sol0); },
                [&]{ filter_by_range(range, triples.begin() + 1 * offset, triples.begin() + 2 * offset, sol1); },
                [&]{ filter_by_range(range, triples.begin() + 2 * offset, triples.begin() + 3 * offset, sol2); },
                [&]{ filter_by_range(range, triples.begin() + 3 * offset, triples.begin() + 4 * offset, sol3); },
                [&]{ filter_by_range(range, triples.begin() + 4 * offset, triples.begin() + 5 * offset, sol4); },
                [&]{ filter_by_range(range, triples.begin() + 5 * offset, triples.begin() + 6 * offset, sol5); },
                [&]{ filter_by_range(range, triples.begin() + 6 * offset, triples.begin() + 7 * offset, sol6); },
                [&]{ filter_by_range(range, triples.begin() + 7 * offset, triples.begin() + 8 * offset, sol7); },
                [&]{ filter_by_range(range, triples.begin() + 8 * offset, triples.begin() + 9 * offset, sol8); },
                [&]{ filter_by_range(range, triples.begin() + 9 * offset, triples.end()               , sol9); }
            );

            solutions.insert(std::end(solutions), std::begin(sol0), std::end(sol0));
            solutions.insert(std::end(solutions), std::begin(sol1), std::end(sol1));
            solutions.insert(std::end(solutions), std::begin(sol2), std::end(sol2));
            solutions.insert(std::end(solutions), std::begin(sol3), std::end(sol3));
            solutions.insert(std::end(solutions), std::begin(sol4), std::end(sol4));
            solutions.insert(std::end(solutions), std::begin(sol5), std::end(sol5));
            solutions.insert(std::end(solutions), std::begin(sol6), std::end(sol6));
            solutions.insert(std::end(solutions), std::begin(sol7), std::end(sol7));
            solutions.insert(std::end(solutions), std::begin(sol8), std::end(sol8));
            solutions.insert(std::end(solutions), std::begin(sol9), std::end(sol9));

            return std::chrono::high_resolution_clock::now() - start;
        }

} // namespace SIMLD