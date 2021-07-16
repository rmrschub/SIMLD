#include <cinttypes>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>
#include <tuple>
#include <bitset>
#include <vector>
#include <iterator>
#include <algorithm>
#include <immintrin.h>
#include <math.h>
#include <tbb/parallel_invoke.h>
#include "Triple.hpp"


namespace SIMLD 
{
    namespace Morton 
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

            inline bool contains(const Triple& element) const
            {
                __m256i min = _mm256_set_epi32(lowerBound.INTERLACED_BITS_MSB, lowerBound.INTERLACED_BITS_NSB, lowerBound.INTERLACED_BITS_LSB, 0x0, 0x0, 0x0, 0x0, 0x0);
                __m256i max = _mm256_set_epi32(upperBound.INTERLACED_BITS_MSB, upperBound.INTERLACED_BITS_NSB, upperBound.INTERLACED_BITS_LSB, UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX);

                __m256i c = _mm256_xor_si256(_mm256_cmpeq_epi64(_mm256_andnot_si256(_mm256_xor_si256(min, max), min), 
                                                                _mm256_andnot_si256(_mm256_xor_si256(min, max), _mm256_set_epi32(element.INTERLACED_BITS_MSB, element.INTERLACED_BITS_NSB, element.INTERLACED_BITS_LSB, 0x0, 0x0, 0x0, 0x0, 0x0))), 
                                            _mm256_set1_epi64x(-1LL));
                
                return _mm256_testz_si256(c, c);
            }

            inline uint_fast8_t avx2_contains(std::vector<Triple>::const_iterator start) const
            {
                __m256i lsb = _mm256_cmpeq_epi32(
                                    _mm256_add_epi32(
                                        _mm256_set1_epi32(0x80000000),
                                        _mm256_andnot_si256(
                                                _mm256_xor_si256(
                                                    _mm256_set1_epi32(lowerBound.INTERLACED_BITS_LSB),
                                                    _mm256_set1_epi32(upperBound.INTERLACED_BITS_LSB)
                                                ),
                                                _mm256_set1_epi32(lowerBound.INTERLACED_BITS_LSB)
                                        )
                                    ),
                                    _mm256_add_epi32(
                                        _mm256_set1_epi32(0x80000000),
                                        _mm256_andnot_si256(
                                                _mm256_xor_si256(
                                                    _mm256_set1_epi32(lowerBound.INTERLACED_BITS_LSB),
                                                    _mm256_set1_epi32(upperBound.INTERLACED_BITS_LSB)
                                                ),
                                                _mm256_set_epi32((start+7)->INTERLACED_BITS_LSB, 
                                                                (start+6)->INTERLACED_BITS_LSB, 
                                                                (start+5)->INTERLACED_BITS_LSB, 
                                                                (start+4)->INTERLACED_BITS_LSB, 
                                                                (start+3)->INTERLACED_BITS_LSB, 
                                                                (start+2)->INTERLACED_BITS_LSB, 
                                                                (start+1)->INTERLACED_BITS_LSB, 
                                                                start->INTERLACED_BITS_LSB)
                                        )
                                    )
                            );

                if (_mm256_testz_si256(lsb, lsb)) { return 0; }
            

                __m256i nsb = _mm256_cmpeq_epi32(
                                    _mm256_add_epi32(
                                        _mm256_set1_epi32(0x80000000),
                                        _mm256_andnot_si256(
                                                _mm256_xor_si256(
                                                    _mm256_set1_epi32(lowerBound.INTERLACED_BITS_NSB),
                                                    _mm256_set1_epi32(upperBound.INTERLACED_BITS_NSB)
                                                ),
                                                _mm256_set1_epi32(lowerBound.INTERLACED_BITS_NSB)
                                        )
                                    ),
                                    _mm256_add_epi32(
                                        _mm256_set1_epi32(0x80000000),
                                        _mm256_andnot_si256(
                                                _mm256_xor_si256(
                                                    _mm256_set1_epi32(lowerBound.INTERLACED_BITS_NSB),
                                                    _mm256_set1_epi32(upperBound.INTERLACED_BITS_NSB)
                                                ), 
                                                _mm256_set_epi32((start+7)->INTERLACED_BITS_NSB, 
                                                                (start+6)->INTERLACED_BITS_NSB, 
                                                                (start+5)->INTERLACED_BITS_NSB, 
                                                                (start+4)->INTERLACED_BITS_NSB, 
                                                                (start+3)->INTERLACED_BITS_NSB, 
                                                                (start+2)->INTERLACED_BITS_NSB, 
                                                                (start+1)->INTERLACED_BITS_NSB, 
                                                                start->INTERLACED_BITS_NSB)
                                        )
                                    )
                            );

                if (_mm256_testz_si256(lsb, lsb)) { return 0; }

                __m256i msb = _mm256_cmpeq_epi32(
                                    _mm256_add_epi32(
                                        _mm256_set1_epi32(0x80000000),
                                        _mm256_andnot_si256(
                                                _mm256_xor_si256(
                                                    _mm256_set1_epi32(lowerBound.INTERLACED_BITS_MSB),
                                                    _mm256_set1_epi32(upperBound.INTERLACED_BITS_MSB)
                                                ),
                                                _mm256_set1_epi32(lowerBound.INTERLACED_BITS_MSB)
                                        )
                                    ),
                                    _mm256_add_epi32(
                                        _mm256_set1_epi32(0x80000000),
                                        _mm256_andnot_si256(
                                                _mm256_xor_si256(
                                                    _mm256_set1_epi32(lowerBound.INTERLACED_BITS_MSB),
                                                    _mm256_set1_epi32(upperBound.INTERLACED_BITS_MSB)
                                                ),
                                                _mm256_set_epi32((start+7)->INTERLACED_BITS_MSB, 
                                                                (start+6)->INTERLACED_BITS_MSB, 
                                                                (start+5)->INTERLACED_BITS_MSB, 
                                                                (start+4)->INTERLACED_BITS_MSB, 
                                                                (start+3)->INTERLACED_BITS_MSB, 
                                                                (start+2)->INTERLACED_BITS_MSB, 
                                                                (start+1)->INTERLACED_BITS_MSB, 
                                                                start->INTERLACED_BITS_MSB)
                                        )
                                    )
                            );

                return _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_and_si256(lsb, _mm256_and_si256(nsb, msb))));
            }
        
            inline void filter_by_range(std::vector<Triple>::const_iterator start, const std::vector<Triple>::const_iterator end, std::vector<Triple> &solutions)
            {
                std::copy_if(start, end, std::back_inserter(solutions), [this] (auto t) { return this->contains(t); } );
            }
    
            inline std::chrono::duration<double> range_lookup(const std::vector<Triple>& triples, std::vector<Triple> &solutions)
            {
                auto start = std::chrono::high_resolution_clock::now();

                std::vector<Triple> sol0, sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8, sol9;
                sol0.reserve(triples.size());
                sol1.reserve(triples.size());
                sol2.reserve(triples.size());
                sol3.reserve(triples.size());
                sol4.reserve(triples.size());
                sol5.reserve(triples.size());
                sol6.reserve(triples.size());
                sol7.reserve(triples.size());
                sol8.reserve(triples.size());
                sol9.reserve(triples.size());


                uint64_t length = std::distance(triples.begin(), triples.end());
                uint64_t offset = ceil(length / 10);

                tbb::parallel_invoke(
                    [&]{ filter_by_range(triples.begin(),              triples.begin() + 1 * offset, sol0); },
                    [&]{ filter_by_range(triples.begin() + 1 * offset, triples.begin() + 2 * offset, sol1); },
                    [&]{ filter_by_range(triples.begin() + 2 * offset, triples.begin() + 3 * offset, sol2); },
                    [&]{ filter_by_range(triples.begin() + 3 * offset, triples.begin() + 4 * offset, sol3); },
                    [&]{ filter_by_range(triples.begin() + 4 * offset, triples.begin() + 5 * offset, sol4); },
                    [&]{ filter_by_range(triples.begin() + 5 * offset, triples.begin() + 6 * offset, sol5); },
                    [&]{ filter_by_range(triples.begin() + 6 * offset, triples.begin() + 7 * offset, sol6); },
                    [&]{ filter_by_range(triples.begin() + 7 * offset, triples.begin() + 8 * offset, sol7); },
                    [&]{ filter_by_range(triples.begin() + 8 * offset, triples.begin() + 9 * offset, sol8); },
                    [&]{ filter_by_range(triples.begin() + 9 * offset, triples.end()               , sol9); }
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

            inline void filter_by_range_avx2(std::vector<Triple>::const_iterator start, const std::vector<Triple>::const_iterator end, std::vector<Triple>& partial_solutions)
            {
                uint_fast64_t tmp = 0;
                for (auto it = start; it != end; it = it+8)
                {
                    int m = avx2_contains(it);
                    tmp += _mm_popcnt_u32(m);

                    while (m) 
                    {
                        int ffs_bit_index = __builtin_ffs(m);
                        m = (m >> ffs_bit_index) << ffs_bit_index;

                        partial_solutions.emplace_back(*(it+(ffs_bit_index-1)));
                    }
                }
            }

            inline std::chrono::duration<double> range_lookup_avx2(const std::vector<Triple>& triples, std::vector<Triple> &solutions)
            {
                auto start = std::chrono::high_resolution_clock::now();

                uint_fast64_t elems_per_thread = (std::lldiv(triples.size(), 10*8)).quot;

                std::vector<Triple> sol0, sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8, sol9;
                sol0.reserve(8* elems_per_thread);
                sol1.reserve(8* elems_per_thread);
                sol2.reserve(8* elems_per_thread);
                sol3.reserve(8* elems_per_thread);
                sol4.reserve(8* elems_per_thread);
                sol5.reserve(8* elems_per_thread);
                sol6.reserve(8* elems_per_thread);
                sol7.reserve(8* elems_per_thread);
                sol8.reserve(8* elems_per_thread);
                sol9.reserve(8* elems_per_thread);
        
                tbb::parallel_invoke(
                    [&]{ filter_by_range_avx2(triples.begin() + 0 * 8 * elems_per_thread, triples.begin() + 1 * 8 * elems_per_thread, sol0); },
                    [&]{ filter_by_range_avx2(triples.begin() + 1 * 8 * elems_per_thread, triples.begin() + 2 * 8 * elems_per_thread, sol1); },
                    [&]{ filter_by_range_avx2(triples.begin() + 2 * 8 * elems_per_thread, triples.begin() + 3 * 8 * elems_per_thread, sol2); },
                    [&]{ filter_by_range_avx2(triples.begin() + 3 * 8 * elems_per_thread, triples.begin() + 4 * 8 * elems_per_thread, sol3); },
                    [&]{ filter_by_range_avx2(triples.begin() + 4 * 8 * elems_per_thread, triples.begin() + 5 * 8 * elems_per_thread, sol4); },
                    [&]{ filter_by_range_avx2(triples.begin() + 5 * 8 * elems_per_thread, triples.begin() + 6 * 8 * elems_per_thread, sol5); },
                    [&]{ filter_by_range_avx2(triples.begin() + 6 * 8 * elems_per_thread, triples.begin() + 7 * 8 * elems_per_thread, sol6); },
                    [&]{ filter_by_range_avx2(triples.begin() + 7 * 8 * elems_per_thread, triples.begin() + 8 * 8 * elems_per_thread, sol7); },
                    [&]{ filter_by_range_avx2(triples.begin() + 8 * 8 * elems_per_thread, triples.begin() + 9 * 8 * elems_per_thread, sol8); },
                    [&]{ filter_by_range_avx2(triples.begin() + 9 * 8 * elems_per_thread, triples.begin() + 10 * 8 * elems_per_thread, sol9); }
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

                // handle remaining elements
                filter_by_range(triples.begin() + 10*8 * elems_per_thread, triples.end(), solutions);

                return std::chrono::high_resolution_clock::now() - start;
            }

            inline void worker_estimate_cardinality_avx2(std::vector<Triple>::const_iterator start, const std::vector<Triple>::const_iterator end, uint_fast64_t& cardinality_estimate)
            {
                uint_fast64_t up_to = 8*((std::lldiv(end-start, 8)).quot);
                uint_fast64_t temp = 0;

                for (auto it = start; it != (start + up_to); it = it+8) { temp += _mm_popcnt_u32(avx2_contains(it)); }
                for (auto it = (start + up_to); it < end; it++) { if (contains(*it)) { temp += 1; } }

                cardinality_estimate = temp;           
            }

            inline std::chrono::duration<double> estimate_cardinality_avx2(std::vector<Triple>::const_iterator start, const std::vector<Triple>::const_iterator end, uint_fast64_t& cardinality_estimate)
            {
                auto now = std::chrono::high_resolution_clock::now();

                uint_fast64_t elems_per_thread = (std::lldiv(end-start, 10*8)).quot;
                uint_fast64_t sol0, sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8, sol9;

                tbb::parallel_invoke(
                    [&]{ worker_estimate_cardinality_avx2(start + 0 * 8 * elems_per_thread, start + 1 * 8 * elems_per_thread, sol0); },
                    [&]{ worker_estimate_cardinality_avx2(start + 1 * 8 * elems_per_thread, start + 2 * 8 * elems_per_thread, sol1); },
                    [&]{ worker_estimate_cardinality_avx2(start + 2 * 8 * elems_per_thread, start + 3 * 8 * elems_per_thread, sol2); },
                    [&]{ worker_estimate_cardinality_avx2(start + 3 * 8 * elems_per_thread, start + 4 * 8 * elems_per_thread, sol3); },
                    [&]{ worker_estimate_cardinality_avx2(start + 4 * 8 * elems_per_thread, start + 5 * 8 * elems_per_thread, sol4); },
                    [&]{ worker_estimate_cardinality_avx2(start + 5 * 8 * elems_per_thread, start + 6 * 8 * elems_per_thread, sol5); },
                    [&]{ worker_estimate_cardinality_avx2(start + 6 * 8 * elems_per_thread, start + 7 * 8 * elems_per_thread, sol6); },
                    [&]{ worker_estimate_cardinality_avx2(start + 7 * 8 * elems_per_thread, start + 8 * 8 * elems_per_thread, sol7); },
                    [&]{ worker_estimate_cardinality_avx2(start + 8 * 8 * elems_per_thread, start + 9 * 8 * elems_per_thread, sol8); },
                    [&]{ worker_estimate_cardinality_avx2(start + 9 * 8 * elems_per_thread, start + 10 * 8 * elems_per_thread, sol9); }

                );
                cardinality_estimate = sol0 + sol1 + sol2 + sol3 + sol4 + sol5 + sol6 + sol7 + sol8 + sol9;

                for (auto it = (start + 10*8*elems_per_thread); it < end; it++) { if (contains(*it)) { cardinality_estimate += 1; } }

                return std::chrono::high_resolution_clock::now() - now;
            }

            inline std::chrono::duration<double>  point_lookup(const Triple& point, const std::vector<Triple>& triples, std::vector<Triple> &solutions)
            {
                auto now = std::chrono::high_resolution_clock::now();

                if (std::binary_search(triples.begin(), triples.end(), point))
                {
                    solutions.push_back(point);
                }

                return std::chrono::high_resolution_clock::now() - now;
            }
        };

    }   // namespace Morton 

    namespace Hilbert
    {

    }   // namespace Hilbert

    namespace Peano
    {

    }   // namespace Peano

} // namespace SIMLD