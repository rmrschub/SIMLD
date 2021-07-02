#include <iostream>
#include <sstream>
#include <tuple>
#include <vector>
#include <cinttypes>
#include <immintrin.h>
#include <emmintrin.h>


#define BMI2_X_MASK 0x49249249  // 01001001001001001001001001001001
#define BMI2_Y_MASK 0x92492492  // 10010010010010010010010010010010
#define BMI2_Z_MASK 0x24924924  // 00100100100100100100100100100100

#define INTERLACED_BITS_LSB interlaced_bits[0]
#define INTERLACED_BITS_NSB interlaced_bits[1]
#define INTERLACED_BITS_MSB interlaced_bits[2]

namespace SIMLD 
{
    struct Triple 
    {
        uint_fast32_t interlaced_bits[3]; // 0: LSB, 1: NSB, 2: MSB

        // constructor
        Triple (uint32_t subject_id, uint32_t predicate_id, uint32_t object_id)
        {
            __m256i morton = _mm256_or_si256(_mm256_set_epi32(_pdep_u32(subject_id   >> 0, BMI2_X_MASK), _pdep_u32(subject_id   >> 11, BMI2_Y_MASK), _pdep_u32(subject_id   >> 22, BMI2_Z_MASK), 0x0, 0x0, 0x0, 0x0, 0x0),
                             _mm256_or_si256(_mm256_set_epi32(_pdep_u32(predicate_id >> 0, BMI2_Y_MASK), _pdep_u32(predicate_id >> 11, BMI2_Z_MASK), _pdep_u32(predicate_id >> 21 ,BMI2_X_MASK), 0x0, 0x0, 0x0, 0x0, 0x0),
                                             _mm256_set_epi32(_pdep_u32(object_id    >> 0, BMI2_Z_MASK), _pdep_u32(object_id    >> 10, BMI2_X_MASK), _pdep_u32(object_id    >> 21, BMI2_Y_MASK), 0x0, 0x0, 0x0, 0x0, 0x0)));

            INTERLACED_BITS_LSB = _mm256_extract_epi32(morton, 7);  // y[10]x[10]z[9]y[9]x[9]...z[0]y[0]x[0]
            INTERLACED_BITS_NSB = _mm256_extract_epi32(morton, 6);  // x[21]z[20]y[20]x[20]...z[11]y[11]x[11]z[10]
            INTERLACED_BITS_MSB = _mm256_extract_epi32(morton, 5);  // z[31]y[31]x[31]...z[22]y[22]x[22]z[21]y[21]
        }

        inline std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t> decode() const
        {
            __m256i triple = _mm256_or_si256(_mm256_set_epi32(_pext_u32(INTERLACED_BITS_MSB, BMI2_Z_MASK) << 22, _pext_u32(INTERLACED_BITS_MSB, BMI2_X_MASK) << 21, _pext_u32(INTERLACED_BITS_MSB, BMI2_Y_MASK) << 21, 0x0, 0x0, 0x0, 0x0, 0x0),
                             _mm256_or_si256(_mm256_set_epi32(_pext_u32(INTERLACED_BITS_NSB, BMI2_Y_MASK) << 11, _pext_u32(INTERLACED_BITS_NSB, BMI2_Z_MASK) << 11, _pext_u32(INTERLACED_BITS_NSB, BMI2_X_MASK) << 10, 0x0, 0x0, 0x0, 0x0, 0x0),
                                             _mm256_set_epi32(_pext_u32(INTERLACED_BITS_LSB, BMI2_X_MASK) <<  0, _pext_u32(INTERLACED_BITS_LSB, BMI2_Y_MASK) <<  0, _pext_u32(INTERLACED_BITS_LSB, BMI2_Z_MASK) <<  0, 0x0, 0x0, 0x0, 0x0, 0x0)));


            return std::make_tuple(_mm256_extract_epi32(triple, 7), _mm256_extract_epi32(triple, 6), _mm256_extract_epi32(triple, 5));
        }

        inline bool operator == (const Triple& rhs) const
        {
            
            __m128i c = _mm_xor_si128(_mm_set_epi32(INTERLACED_BITS_MSB,     INTERLACED_BITS_NSB,     INTERLACED_BITS_LSB,     0x0),
                                      _mm_set_epi32(rhs.INTERLACED_BITS_MSB, rhs.INTERLACED_BITS_NSB, rhs.INTERLACED_BITS_LSB, 0x0));
            
            return _mm_testz_si128(c, c);
        }     

        inline bool operator != (const Triple& rhs) const
        {
            return !(*this == rhs);
        }

        inline bool operator < (const Triple& rhs) const
        {
            if (INTERLACED_BITS_MSB < rhs.INTERLACED_BITS_MSB) { return true; }
            else if ((INTERLACED_BITS_MSB == rhs.INTERLACED_BITS_MSB) && (INTERLACED_BITS_NSB < rhs.INTERLACED_BITS_NSB)) { return true; }
            else if ((INTERLACED_BITS_MSB == rhs.INTERLACED_BITS_MSB) && (INTERLACED_BITS_NSB == rhs.INTERLACED_BITS_NSB) && (INTERLACED_BITS_LSB < rhs.INTERLACED_BITS_LSB)) { return true; }
            
            return false; 
        };

        inline bool operator > (const Triple& rhs) const
        {
            return (rhs < *this);
        }

        inline bool operator <= (const Triple& rhs) const
        {
            if (INTERLACED_BITS_MSB <= rhs.INTERLACED_BITS_MSB) { return true; }
            else if ((INTERLACED_BITS_MSB == rhs.INTERLACED_BITS_MSB) && (INTERLACED_BITS_NSB <= rhs.INTERLACED_BITS_NSB)) { return true; }
            else if ((INTERLACED_BITS_MSB == rhs.INTERLACED_BITS_MSB) && (INTERLACED_BITS_NSB == rhs.INTERLACED_BITS_NSB) && (INTERLACED_BITS_LSB <= rhs.INTERLACED_BITS_LSB)) { return true; }
            
            return false; 
        }

        inline bool operator >= (const Triple& rhs) const
        {
            return (rhs <= *this);
        }

        inline Triple operator + (const Triple& rhs) const
        {
            Triple sum = Triple(0,0,0);
            unsigned char carry_out = 0;
            
            carry_out = _addcarry_u32(carry_out, INTERLACED_BITS_LSB, rhs.INTERLACED_BITS_LSB, &(sum.INTERLACED_BITS_LSB));
            carry_out = _addcarry_u32(carry_out, INTERLACED_BITS_NSB, rhs.INTERLACED_BITS_NSB, &(sum.INTERLACED_BITS_NSB));
            carry_out = _addcarry_u32(carry_out, INTERLACED_BITS_MSB, rhs.INTERLACED_BITS_MSB, &(sum.INTERLACED_BITS_MSB));

            return sum;
        }

        inline Triple operator - (const Triple& rhs) const
        {
            Triple diff = Triple(0,0,0);
            unsigned char borrow_from = 0;
            borrow_from = _subborrow_u32(borrow_from, INTERLACED_BITS_LSB, rhs.INTERLACED_BITS_LSB, &(diff.INTERLACED_BITS_LSB));
            borrow_from = _subborrow_u32(borrow_from, INTERLACED_BITS_NSB, rhs.INTERLACED_BITS_NSB, &(diff.INTERLACED_BITS_NSB));
            borrow_from = _subborrow_u32(borrow_from, INTERLACED_BITS_MSB, rhs.INTERLACED_BITS_MSB, &(diff.INTERLACED_BITS_MSB));
            
            return diff;
        }

    }; 

    std::ostream& operator << (std::ostream& o, Triple& triple)
    {
        std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t> coordinates = triple.decode();

        std::stringstream ss;

        ss << "(";
        ss << (uint_fast32_t)std::get<0>(coordinates) << ", ";
        ss << (uint_fast32_t)std::get<1>(coordinates) << ", ";
        ss << (uint_fast32_t)std::get<2>(coordinates) << ")";

        std::ios_base::sync_with_stdio(false);
        std::cin.tie(NULL);
        std::cout << ss.str();

        return o;
    }

} // namespace SIMLD