#ifndef VEC_SIZE_H
#define VEC_SIZE_H
#include "dtw_opt_feature_t.h"

constexpr int get_aligned_bit(int bits){
	if(bits<=8) return 8;
	if(bits<=16) return 16;
	if(bits<=32) return 32;
	if(bits<=64) return 64;
	if(bits<=128) return 128;
	if(bits<=256) return 256;
	return 512;
}

constexpr int ALIGNED_BIT_LENGTH = get_aligned_bit(FEATURE_T_BIT_LENGTH);
constexpr int vec_size = 512/ALIGNED_BIT_LENGTH;

#endif
