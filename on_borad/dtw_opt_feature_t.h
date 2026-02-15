#ifndef DTW_OPT_FEATURE_T_H
#define DTW_OPT_FEATURE_T_H

#include "ap_int.h"

const int TS_MIN = -10;
const int TS_MAX = 10;
// min and max value of time series data
// These should be edited depending on imput time series data

const int diff_max = TS_MAX - TS_MIN;

constexpr int get_bit_feature(int num){
	return (num<=1) ? 1 : 1 + get_bit_feature(num>>1);
}

constexpr int FEATURE_T_BIT_LENGTH = get_bit_feature(diff_max) + 1;

typedef ap_int<FEATURE_T_BIT_LENGTH> feature_t;

#endif
