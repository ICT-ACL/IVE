#pragma once

#include <x86intrin.h>

#define MM_LOAD_SI_128 _mm_loadu_si128
#define MM_STORE_SI_128 _mm_storeu_si128

constexpr int cyclic_shift1 = _MM_SHUFFLE(0,3,2,1); //rotating right
constexpr int cyclic_shift2 = _MM_SHUFFLE(2,1,0,3); //rotating left
constexpr int cyclic_shift3 = _MM_SHUFFLE(1,0,3,2); //between

static const __m128i all_zero_si128 = _mm_setzero_si128();
static const __m128i all_one_si128 = _mm_set_epi32(0xffffffff, 0xffffffff,
        0xffffffff, 0xffffffff);

static const uint8_t shuffle_pi8_array[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 
    8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 
    12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 
    8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
};
static const __m128i *shuffle_mask = (__m128i*)(shuffle_pi8_array);

uint32_t* prepare_shuffling_dict_avx()
{
    uint32_t* arr = new uint32_t[2048];
    for(int i=0; i<256; ++i){
        int count=0, rest=7;
        for(int b=0; b<8; ++b){
            if(i & (1 << b)){
                // n index at pos p - move nth element to pos p
                arr[i*8 + count] = b; // move all set bits to beginning
                ++count;
            }else{
                arr[i*8 + rest] = b; // move rest at the end
                --rest;
            }
        }
    }
    return arr;
}
static const uint32_t *shuffle_mask_avx =prepare_shuffling_dict_avx();

int* prepare_byte_check_mask_dict2()
{
    int * mask = new int[65536];

    auto trans_c_s = [](const int c) -> int {
        switch (c) {
            case 0: return -1; // no match
            case 1: return 0;
            case 2: return 1;
            case 4: return 2;
            case 8: return 3;
            default: return 4; // multiple matches.
        }
    };

    for (int x = 0; x < 65536; ++x) {        
        int c0 = (x & 0xf), c1 = ((x >> 4) & 0xf);
        int c2 = ((x >> 8) & 0xf), c3 = ((x >> 12) & 0xf);
        int s0 = trans_c_s(c0), s1= trans_c_s(c1);
        int s2 = trans_c_s(c2), s3 = trans_c_s(c3);
        
        bool is_multiple_match = (s0 == 4) || (s1 == 4) ||
                (s2 == 4) || (s3 == 4);
        if (is_multiple_match) {
            mask[x] = -1;
            continue;
        }
        bool is_no_match = (s0 == -1) && (s1 == -1) &&
                (s2 == -1) && (s3 == -1);
        if (is_no_match) {
            mask[x] = -2;
            continue;
        }
        if (s0 == -1) s0 = 0; if (s1 == -1) s1 = 1;
        if (s2 == -1) s2 = 2; if (s3 == -1) s3 = 3;
        mask[x] = (s0) | (s1 << 2) | (s2 << 4) | (s3 << 6);        
    }

    return mask;
}
static const int *byte_check_mask_dict = prepare_byte_check_mask_dict2();

uint8_t * prepare_match_shuffle_dict2()
{
    uint8_t * dict = new uint8_t[4096];

    for (int x = 0; x < 256; ++x) {
        for (int i = 0; i < 4; ++i) {
            uint8_t c = (x >> (i << 1)) & 3; // c = 0, 1, 2, 3
            int pos = x * 16 + i * 4;
            for (uint8_t j = 0; j < 4; ++j)
                dict[pos + j] = c * 4 + j;
        }
    }

    return dict;
}
static const __m128i *match_shuffle_dict = (__m128i*)prepare_match_shuffle_dict2();

static const uint8_t byte_check_group_a_pi8[64] = {
    0, 0, 0, 0, 4, 4, 4, 4, 8, 8, 8, 8, 12, 12, 12, 12,
    1, 1, 1, 1, 5, 5, 5, 5, 9, 9, 9, 9, 13, 13, 13, 13,
    2, 2, 2, 2, 6, 6, 6, 6, 10, 10, 10, 10, 14, 14, 14, 14,
    3, 3, 3, 3, 7, 7, 7, 7, 11, 11, 11, 11, 15, 15, 15, 15,
};
static const uint8_t byte_check_group_b_pi8[64] = {
    0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12,
    1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13,
    2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 14,
    3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15,
};
static const __m128i *byte_check_group_a_order = (__m128i*)(byte_check_group_a_pi8);
static const __m128i *byte_check_group_b_order = (__m128i*)(byte_check_group_b_pi8);
int intersect_qfilter_bsr_b4_v2(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    __m128i base_a = _mm_lddqu_si128((__m128i*)bases_a + i);
    __m128i base_b = _mm_lddqu_si128((__m128i*)bases_b + j);
    __m128i byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
    __m128i byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
    __m128i cmp_mask, and_state;

    while (i < qs_a && j < qs_b) {
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];

        if (__builtin_expect(ms_order != -2, 0)) {
            __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
            __m128i state_b = _mm_lddqu_si128((__m128i*)(states_b + j));    
            if (ms_order > 0) {                
                __m128i sf_base_b = _mm_shuffle_epi8(base_b, match_shuffle_dict[ms_order]);
                __m128i sf_state_b = _mm_shuffle_epi8(state_b, match_shuffle_dict[ms_order]);
                cmp_mask = _mm_cmpeq_epi32(base_a, sf_base_b);                
                and_state = _mm_and_si128(state_a, sf_state_b);
                __m128i state_mask = _mm_cmpeq_epi32(and_state, all_zero_si128);
                cmp_mask = _mm_andnot_si128(state_mask, cmp_mask);
            } else {
                __m128i cmp_mask0 = _mm_cmpeq_epi32(base_a, base_b);
                __m128i state_c0 = _mm_and_si128(
                        _mm_and_si128(state_a, state_b), cmp_mask0);
                __m128i base_sf1 = _mm_shuffle_epi32(base_b, cyclic_shift1);
                __m128i state_sf1 = _mm_shuffle_epi32(state_b, cyclic_shift1);
                __m128i cmp_mask1 = _mm_cmpeq_epi32(base_a, base_sf1);
                __m128i state_c1 = _mm_and_si128(
                        _mm_and_si128(state_a, state_sf1), cmp_mask1);
                __m128i base_sf2 = _mm_shuffle_epi32(base_b, cyclic_shift2);
                __m128i state_sf2 = _mm_shuffle_epi32(state_b, cyclic_shift2);
                __m128i cmp_mask2 = _mm_cmpeq_epi32(base_a, base_sf2);
                __m128i state_c2 = _mm_and_si128(
                        _mm_and_si128(state_a, state_sf2), cmp_mask2);
                __m128i base_sf3 = _mm_shuffle_epi32(base_b, cyclic_shift3);
                __m128i state_sf3 = _mm_shuffle_epi32(state_b, cyclic_shift3);
                __m128i cmp_mask3 = _mm_cmpeq_epi32(base_a, base_sf3);
                __m128i state_c3 = _mm_and_si128(
                        _mm_and_si128(state_a, state_sf3), cmp_mask3);
                and_state = _mm_or_si128(
                        _mm_or_si128(state_c0, state_c1),
                        _mm_or_si128(state_c2, state_c3)
                        );                
                __m128i state_mask = _mm_cmpeq_epi32(and_state, all_zero_si128);
                cmp_mask = _mm_andnot_si128(state_mask, all_one_si128);
            } 

            int mask = _mm_movemask_ps((__m128)cmp_mask);
            __m128i res_b = _mm_shuffle_epi8(base_a, shuffle_mask[mask]);
            __m128i res_s = _mm_shuffle_epi8(and_state, shuffle_mask[mask]);
            _mm_storeu_si128((__m128i*)(bases_c + size_c), res_b);
            _mm_storeu_si128((__m128i*)(states_c + size_c), res_s);

            size_c += _mm_popcnt_u32(mask);        
        }

        int a_max = bases_a[i + 3];
        int b_max = bases_b[j + 3];
        if (a_max <= b_max) {
            i += 4;
            base_a = MM_LOAD_SI_128((__m128i*)(bases_a + i));
            byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
            _mm_prefetch((char*) (bases_a + i + 16), _MM_HINT_T0);
            _mm_prefetch((char*) (states_a + i + 16), _MM_HINT_T0);
        }
        if (a_max >= b_max) {
            j += 4;
            base_b = MM_LOAD_SI_128((__m128i*)(bases_b + j));
            byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
            _mm_prefetch((char*) (bases_b + j + 16), _MM_HINT_T0);
            _mm_prefetch((char*) (states_b + j + 16), _MM_HINT_T0);
        }
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

const int BSR_WIDTH = sizeof(int) * 8;
const int BSR_SHIFT = __builtin_ctzll(BSR_WIDTH);
const int BSR_MASK = BSR_WIDTH - 1;
int offline_uint_trans_bsr(int *set_a, int size_a, int *bases_a, int *states_a)
{
    int cnt = -1;
    for (int i = 0; i < size_a; ++i) {
        int u = set_a[i];
        int u_base = (u >> BSR_SHIFT);
        int u_bit = (1 << (u & BSR_MASK));
        if (cnt == -1 || bases_a[cnt] != u_base) {
            bases_a[++cnt] = u_base;
            states_a[cnt] = u_bit;
        } else {
            states_a[cnt] |= u_bit;
        }
    }
    return ++cnt;    
}

int offline_bsr_trans_uint(int *bases_a, int *states_a, int size_a, int *set_a)
{
    int cnt = 0;
    for (int i = 0; i < size_a; ++i) {
        int u_high = (bases_a[i] << BSR_SHIFT);
        int state = states_a[i];
        while (state) {
            int u = (u_high | __builtin_ctz(state));
            set_a[cnt++] = u;
            state &= (state - 1);
        }
    }
    return cnt;
}
