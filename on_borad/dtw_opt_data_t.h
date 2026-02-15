#ifndef DTW_OPT_DATA_T_H
#define DTW_OPT_DATA_T_H

#include "ap_int.h"
#include "dtw_opt_feature_t.h"

/*

This header file is for defining the bit length of 'data_t' and 'feature_t'
1) diff_max * diff_max * (num_features): maximum normal value
2) INFINITY + TS_MAX * TS_MAX * (num_features/2): maximum trash value

Constraints)
1. data_t >= diff_max * diff_max * (num_features)

2. INFINITY + TS_MAX * TS_MAX * (num_features/2) is maximum trash value.
So, (INFINITY value for data_t) = 2^(bit length) - 1 - TS_MAX * TS_MAX * (num_features/2)
Therefore, data_t won't be overlflow when maximum trash value occurs.

3. INF > (diff_max * diff_max * (num_features)) - (diff_max * diff_max)
So, 2^(bit length)-1 > ((diff_max * diff_max * (num_features)) - (diff_max * diff_max)) + (TS_MAX * TS_MAX * (num_features/2))

Goal)
Given threshold, Find the minmum bit length to accommodate the threshold

*/
const int num_features = 512; // should be larger than vec_size
const int threshold = ((diff_max * diff_max * (num_features)) - (diff_max * diff_max))
+ (TS_MAX * TS_MAX * (num_features/2));

constexpr int get_bit_data(int num){ // finding min bit length which can hold threshold
	return (num<=1) ? 1 : 1 + get_bit_data(num>>1);
}

constexpr int DATA_T_BIT_LENGTH = get_bit_data(threshold + 1);
// should be greater than. So, pass threshold + 1 as an argument.

typedef ap_uint<DATA_T_BIT_LENGTH> data_t;

const data_t INF = ~((data_t)0) - (data_t)(TS_MAX * TS_MAX  * (num_features / 2));

#endif
