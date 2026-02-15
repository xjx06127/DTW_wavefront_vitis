#include <algorithm>
#include "ap_int.h"
#include <hls_vector.h>
#include "dtw_opt_data_t.h"
#include "dtw_opt_feature_t.h"
#include "vec_size.h"

#define min(a, b) ((a<b) ? a : b)
#define max(a, b) ((a>b) ? a : b)
#define DIAG 0
#define LEFT 1
#define ABOVE 2


typedef hls::vector<ap_uint<2>, num_features> route_vec_t;

extern "C" {
    void dtw_wavefront_opt6(hls::vector<feature_t, vec_size> *data_A,
    hls::vector<feature_t, vec_size> *data_B, data_t *result, route_vec_t *routes) {

        #pragma HLS INTERFACE m_axi port = data_A bundle = a depth = num_features/vec_size
        #pragma HLS INTERFACE m_axi port = data_B bundle = b depth = num_features/vec_size
        #pragma HLS INTERFACE m_axi port = result bundle = c depth = 1
        #pragma HLS INTERFACE m_axi port = routes bundle = d depth = (2*num_features-1)

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
        // PE(#1 ~ #lane-1) can access local_A[] via this buffer
		#pragma HLS ARRAY_PARTITION variable = shifted_buff_A complete

        static route_vec_t routes_tmp[2*num_features-1];
        #pragma HLS BIND_STORAGE variable=routes_tmp type=ram_2p impl=bram
        #pragma HLS AGGREGATE variable=routes_tmp compact=auto

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

			route_vec_t cur_route_vec;

            Compute_Inner_Loop:
            for (int l = 1; l < lane; l++) {
            #pragma HLS unroll
                //j = k - l;
                //J = l;
                    
                const feature_t dist = shifted_buff_A[l-1] - local_B[l - 1];
                const data_t diff = dist*dist;
                #pragma HLS BIND_OP variable=diff op=mul impl=dsp
				
				//priority: DIAG > LEFT > ABOVE
				data_t min_val = bfr_buff[l-1];
				ap_uint<2> direction = DIAG;

				if(src_buff[l-1] < min_val){
					min_val = src_buff[l-1];
					direction = LEFT;
				}
				if(src_buff[l] < min_val){
					min_val = src_buff[l];
					direction = ABOVE;
				}

                trg_buff[l] = diff + min_val;
                #pragma HLS BIND_OP variable=trg_buff op=add impl=fabric

				cur_route_vec[l-1] = direction;
            }
			
			routes_tmp[k-2] = cur_route_vec;

			Transfer_Loop:
			for(int i=0; i<lane; i++){
			#pragma HLS unroll
				bfr_buff[i] = src_buff[i];
				src_buff[i] = trg_buff[i];
			}
        }

		Write_Loop:
		for(int k=0; k<2*lane-3;k++){
		#pragma HLS PIPELINE II=1
			routes[k] = routes_tmp[k];
		}

        *result = trg_buff[lane-1];
    }
}
