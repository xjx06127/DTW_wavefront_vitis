#include "ap_int.h"
#include <hls_vector.h>
#include "dtw_opt_data_t.h"
#include "dtw_opt_feature_t.h"
#include "vec_size.h"
#define min(a, b) ((a<b) ? a : b)
#define max(a, b) ((a>b) ? a : b)

extern "C" {
    void wavefront_DTW(hls::vector<feature_t, vec_size> *data_A,
    hls::vector<feature_t, vec_size> *data_B, data_t *result) {

        #pragma HLS INTERFACE m_axi port = data_A bundle = a
        #pragma HLS INTERFACE m_axi port = data_B bundle = b
        #pragma HLS INTERFACE m_axi port = result bundle = c

        const int lane = num_features + 1;

        static data_t bfr_buff[lane];
        static data_t src_buff[lane];
        static data_t trg_buff[lane];
        #pragma HLS ARRAY_PARTITION variable = bfr_buff complete	
        #pragma HLS ARRAY_PARTITION variable = src_buff complete
        #pragma HLS ARRAY_PARTITION variable = trg_buff complete

        static feature_t local_A[num_features], local_B[num_features];
        #pragma HLS ARRAY_PARTITION variable = local_A complete
        #pragma HLS ARRAY_PARTITION variable = local_B complete

		static feature_t shifted_buff_A[num_features];
        // PE(#1 ~ #lane-1) can access local_A[] via this buffer.
		#pragma HLS ARRAY_PARTITION variable = shifted_buff_A complete

        Copy_Loop:
        for (int i = 0; i < num_features / vec_size; i++) {
        #pragma HLS PIPELINE II = 1
            
            hls::vector<feature_t, vec_size> chunk_A = data_A[i];
            hls::vector<feature_t, vec_size> chunk_B = data_B[i];

            Unpack_Loop:
            for (int j = 0; j < vec_size; j++) {
            #pragma HLS UNROLL
                local_A[i * vec_size + j] = chunk_A[j];
                local_B[i * vec_size + j] = chunk_B[j];
            }
        }

        M_Init_Loop:
        for (int i = 0; i < lane; i++) {
        #pragma HLS unroll
            bfr_buff[i] = INF;
            src_buff[i] = INF;
            trg_buff[i] = INF;
        }
        bfr_buff[0] = 0;
        
        Shifted_buff_A_Init_Loop:
        for (int i=0; i < num_features; i++){
        #pragma HLS unroll
            shifted_buff_A[i] = 0;
        }

        Compute_Outer_Loop:
        // Target cycle: 2*lane-3
        for (int k = 2; k < 2 * lane - 1; k++) {
        #pragma HLS PIPELINE II = 1

		    for(int i=lane-2; i>0; i--){
			#pragma HLS unroll
				shifted_buff_A[i] = shifted_buff_A[i-1];
			}
            
			if((k-2) < num_features) shifted_buff_A[0] = local_A[k-2];

            Compute_Inner_Loop:
            for (int l = 1; l < lane; l++) {
            #pragma HLS unroll
                //j = k - l;
                //J = l;
                    
                const feature_t dist = shifted_buff_A[l-1] - local_B[l - 1];
                const data_t diff = dist*dist;
                #pragma HLS BIND_OP variable=diff op=mul impl=dsp
                const data_t min1 = min(src_buff[l], src_buff[l-1]);
                const data_t min2 = min(bfr_buff[l-1], min1);

                trg_buff[l] = diff + min2;
                #pragma HLS BIND_OP variable=trg_buff op=add impl=fabric
            }
			Transfer_Loop:
			for(int i=0; i<lane; i++){
			#pragma HLS unroll
				bfr_buff[i] = src_buff[i];
				src_buff[i] = trg_buff[i];
			}
        }
        *result = trg_buff[lane-1];
    }
}
