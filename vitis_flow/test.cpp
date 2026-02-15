#include <math.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>
#include "ap_int.h"
#include "hls_vector.h"
#include "dtw_opt_data_t.h"
#include "dtw_opt_feature_t.h"

#define min(a, b) ((a<b) ? a : b)

#define DIAG 0
#define LEFT 1
#define ABOVE 2

typedef hls::vector<ap_uint<2>, num_features> route_vec_t;

extern "C" {
void dtw_wavefront_opt(data_t* data_A, data_t* data_B, data_t* result);
}
extern "C" {
void dtw_wavefront_opt2(data_t* data_A, data_t* data_B, data_t* result);
}
extern "C" {
void dtw_wavefront_opt3(data_t* data_A, data_t* data_B, data_t* result);
}
extern "C" {
void dtw_wavefront_opt4(feature_t* data_A, feature_t* data_B, data_t* result);
}
extern "C" {
void dtw_wavefront_opt5(feature_t* data_A, feature_t* data_B, data_t* result);
}
extern "C"{
void dtw_wavefront_opt6(feature_t* data_A, feature_t* data_B, data_t* result, route_vec_t* routes);
}

void dtw_sw(feature_t* data_A, feature_t* data_B, data_t* result, std::vector<std::pair<int,int>>& path){
    const int lane = num_features + 1;

    data_t* M = new data_t[lane*lane];
	int* Dir_sw = new int[lane*lane];
	
	for(int i=0; i<lane*lane; i++){
		Dir_sw[i] = -1;
	}

    for(int i=1; i<lane; i++){
        M[i] = M[i*lane] = INF;
    }
    M[0] = 0;

    for(int row=1; row<lane; row++){
        const feature_t value_A = data_A[row-1];
        for(int col=1; col<lane; col++){
            const data_t diag = M[(row-1)*lane+col-1];
            const data_t abve = M[(row-1)*lane+col];
            const data_t left = M[(row)*lane+col-1];

            const feature_t diff = value_A - data_B[col-1];

			data_t min_val = diag;
			int d = DIAG;
			
			if(left < min_val){
				min_val = left;
				d = LEFT;
			}
			if(abve < min_val){
				min_val = abve;
				d = ABOVE;
			}

			M[row*lane + col] = diff*diff + min_val;
			Dir_sw[row*lane + col] = d;
        }
    }
    *result = M[lane*lane-1];
	int r = num_features;
	int c = num_features;
	path.push_back({r, c});

	while(r>1 || c>1){ // when r==1 && c==1, break
		int d = Dir_sw[r*lane+c];
		if(d==DIAG){
			r--;
			c--;
		}
		else if(d==LEFT){
			c--;
		}
		else if(d==ABOVE){
			r--;
		}
		else{
			std::cout<<"ERROR: no direction value"<<std::endl;
		}
		path.push_back({r,c});
	}

    delete[] M;
	delete[] Dir_sw;
}

void reconstruct_path(route_vec_t* routes, std::vector<std::pair<int,int>>& path){
	const int lane = num_features+1;

	int* Dir = new int[lane*lane];

	for(int i=0; i<lane*lane; i++){
		Dir[i] = -1;
	}

	for(int k=2; k<2*lane-1; k++){
		route_vec_t cur_vec = routes[k-2];

		for(int l=1; l<lane; l++){
			int r = k-l;
			int c = l;

			const bool outside = (k<=l || c==0 || r>=lane);
			if(!outside){
				Dir[r*lane + c] = (int)cur_vec[l-1];
			}
		}
	}

	int r = num_features;
    int c = num_features;
    path.push_back({r, c});

    while (r > 1 || c > 1) {
        int d = Dir[r * lane + c];
        if (d == DIAG) {
            r--;
			c--;
        } else if (d == LEFT) {
            c--;
        } else if (d== ABOVE) {
            r--;
        }else{
			std::cout<<"ERROR: no direction value"<<std::endl;
		}
        path.push_back({r, c});
    }
    delete[] Dir;
}

int main(){
    feature_t data1[num_features]; 
	feature_t data2[num_features];
    feature_t data3[num_features];

	
	route_vec_t* route_out_1 = new route_vec_t[2*num_features-1];
    route_vec_t* route_out_2 = new route_vec_t[2*num_features-1];
	// 2*lane-3 elements == 2*num_features-1 elements. Each element has 512 elements(each is 2bit).

    data_t result1_hw;
    data_t result2_hw;
    data_t result1_sw;
    data_t result2_sw;

	std::vector<std::pair<int,int>> path1_hw, path1_sw;
	std::vector<std::pair<int,int>> path2_hw, path2_sw;
    // 512 <= size <=2*num_features - 1

    for(int i=0; i<num_features; i++){
        data1[i] = (feature_t)((int)(10*sin(i))); //10 * sin
        data2[i] = (feature_t)((int)(10*sin(i+1))); //shifted 10 * sin
        data3[i] = (feature_t)((int)(10*cos(i))); //10 * cos
    }

    dtw_wavefront_opt6(data1, data2, &result1_hw, route_out_1);
	reconstruct_path(route_out_1, path1_hw);
    dtw_sw(data1, data2, &result1_sw, path1_sw);

    int naive_result1=0; // for naive comparison between data1, data2
    int naive_result2=0; // for naive comparison between data1, data2

    for(int i=0; i<num_features; i++){
        naive_result1 += (data1[i]-data2[i])*(data1[i]-data2[i]);
		naive_result2 += (data1[i]-data3[i])*(data1[i]-data3[i]);
    }  

    if(result1_hw == result1_sw){
        std::cout<<"TEST1 PASSED(DTW SCORE MATCHING)"<<std::endl;
    }
    else{
        std::cout<<"TEST1(DTW SCORE MATCHING) FAILED"<<std::endl;
    }
    std::cout<<"result1_hw: "<<result1_hw<<std::endl;
    std::cout<<"result1_sw: "<<result1_sw<<std::endl;
    std::cout<<"naive_result1: "<<naive_result1<<std::endl;
    std::cout<<std::endl;

	bool match = true;
	if(path1_hw.size() != path1_sw.size()){
		match = false;
	}
	else{
		for(int i=path1_hw.size() - 1; i>=0; i--){
			if(path1_hw[i] != path1_sw[i]){
				match = false;
				std::cout << "Mismatch at " << path1_hw.size() - i << ": HW("
                          << path1_hw[i].first << "," << path1_hw[i].second << ") vs SW("
                          << path1_sw[i].first << "," << path1_sw[i].second << ")" << std::endl;
                break;
			}
			else{
				std::cout << path1_hw.size() - i << ": HW("
                          << path1_hw[i].first << "," << path1_hw[i].second << ") vs SW("
                          << path1_sw[i].first << "," << path1_sw[i].second << ")" << std::endl;
			}
		}
	}
	if(match) std::cout << "TEST1(DTW PATH MATCHING) PASSED" << std::endl;
	else std::cout << "TEST1(DTW PATH MATCHING) FAILED" << std::endl;
    std::cout << std::endl;

	dtw_wavefront_opt6(data1, data3, &result2_hw, route_out_2);
    reconstruct_path(route_out_2, path2_hw);
    dtw_sw(data1, data3, &result2_sw, path2_sw);


    if(result2_hw == result2_sw){
        std::cout<<"TEST2 PASSED(DTW SCORE MATCHING)"<<std::endl;
    }
    else{
        std::cout<<"TEST2 FAILED(DTW SCORE MATCHING)"<<std::endl;
    }
    std::cout<<"result2_hw: "<<result2_hw<<std::endl;
    std::cout<<"result2_sw: "<<result2_sw<<std::endl;
    std::cout<<"naive_result2: "<<naive_result2<<std::endl;
    std::cout<<std::endl;

    bool match2 = true;
    if (path2_hw.size() != path2_sw.size()) {match2 = false;}
    else {
        for (int i = path2_hw.size() - 1; i>=0; i--) {
            if (path2_hw[i] != path2_sw[i]) {
                match2 = false; 
                std::cout << "Mismatch at " << path2_hw.size() - i << ": HW(" 
                          << path2_hw[i].first << "," << path2_hw[i].second << ") vs SW(" 
                          << path2_sw[i].first << "," << path2_sw[i].second << ")" << std::endl;
                break;
            }
			else {
				std::cout  << path2_hw.size() - i << ": HW("
                           << path2_hw[i].first << "," << path2_hw[i].second << ") vs SW("
                           << path2_sw[i].first << "," << path2_sw[i].second << ")" << std::endl;

			}
        }
    }
	if (match2) std::cout << "TEST2(DTW PATH MATCHING) PASSED" << std::endl;
    else std::cout << "TEST2(DTW PATH MATCHING) FAILED" << std::endl;
	std::cout<<std::endl;

    std::cout<<"data_t: ap_uint<"<<DATA_T_BIT_LENGTH<<">"<<std::endl;
    std::cout<<"feature_t: ap_int<"<<FEATURE_T_BIT_LENGTH<<">"<<std::endl;
    std::cout<<std::endl;

	delete[] route_out_1;
    delete[] route_out_2;
}
