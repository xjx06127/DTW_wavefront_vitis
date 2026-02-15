#include "cmdlineparser.h"
#include <iostream>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <chrono>

#include "dtw_opt_data_t.h"
#include "dtw_opt_feature_t.h"

// XRT includes
#include "xrt/xrt_bo.h"
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"

#define min(a, b) ((a<b) ? a : b)

void dtw_sw(feature_t* data_A, feature_t* data_B, data_t* result){
    const int lane = num_features + 1;
    data_t* M = new data_t[lane*lane];
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
            M[row*lane + col] = diff*diff + min(diag, min(abve,left));
        }
    }
    *result = M[lane*lane-1];
   delete[] M;
}


int main(int argc, char** argv) {
    // Command Line Parser
    sda::utils::CmdLineParser parser;

    // Switches
    //**************//"<Full Arg>",  "<Short Arg>", "<Description>", "<Default>"
    parser.addSwitch("--xclbin_file", "-x", "input binary file string", "");
    parser.addSwitch("--device_id", "-d", "device index", "0");
    parser.parse(argc, argv);

    // Read settings
    std::string binaryFile = parser.value("xclbin_file");
    int device_index = stoi(parser.value("device_id"));

    if (argc < 3) {
        parser.printHelp();
        return EXIT_FAILURE;
    }

    std::cout << "Open the device" << device_index << std::endl;
    auto device = xrt::device(device_index);
    std::cout << "Load the xclbin " << binaryFile << std::endl;
    auto uuid = device.load_xclbin(binaryFile);

    size_t vector_size_bytes = sizeof(feature_t) * num_features;
	size_t vector_size_bytes_out = sizeof(data_t);

    auto krnl = xrt::kernel(device, uuid, "wavefront_DTW");

    std::cout << "Allocate Buffer in Global Memory\n";
    auto bo0 = xrt::bo(device, vector_size_bytes, krnl.group_id(0));
    auto bo1 = xrt::bo(device, vector_size_bytes, krnl.group_id(1));
    auto bo_out = xrt::bo(device, vector_size_bytes_out, krnl.group_id(2));

    // Map the contents of the buffer object into host memory
    auto bo0_map = bo0.map<feature_t*>();
    auto bo1_map = bo1.map<feature_t*>();
    auto bo_out_map = bo_out.map<data_t*>();
    std::fill(bo0_map, bo0_map + num_features, 0);
    std::fill(bo1_map, bo1_map + num_features, 0);
    std::fill(bo_out_map, bo_out_map + 1, 0);

	data_t naive_result = 0;
    for (int i = 0; i < num_features; ++i) {
        bo0_map[i] = (feature_t)((int)(10*sin(i)));
        bo1_map[i] = (feature_t)((int)(10*sin(i+1)));
		naive_result += (bo0_map[i]-bo1_map[i]) * (bo0_map[i]-bo1_map[i]);
    }

    // Synchronize buffer content with device side
    std::cout << "synchronize input buffer data to device global memory\n";

	auto h2d_start = std::chrono::steady_clock::now();
    bo0.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
	auto h2d_end = std::chrono::steady_clock::now();

    std::cout << "Execution of the kernel\n";
	auto krnl_start = std::chrono::steady_clock::now();
    auto run = krnl(bo0, bo1, bo_out);
    run.wait();
	auto krnl_end = std::chrono::steady_clock::now();

    // Get the output;
    std::cout << "Get the output data from the device" << std::endl;
    bo_out.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
	
	data_t test_result;
	dtw_sw(bo0_map, bo1_map, &test_result);

    // Validate our results
    if (test_result == *bo_out_map) std::cout<<"TEST PASSED"<<std::endl;
	else std::cout<<"TEST FAILED"<<std::endl;
	std::cout<<std::endl;

	std::cout<<"result_hw: "<<*bo_out_map<<std::endl;
	std::cout<<"result_sw: "<<test_result<<std::endl;
	std::cout<<"result_naive: "<<naive_result<<std::endl;
	std::cout<<std::endl;

    std::cout<<"data_t: ap_uint<"<<DATA_T_BIT_LENGTH<<">"<<std::endl;
    std::cout<<"feature_t: ap_int<"<<FEATURE_T_BIT_LENGTH<<">"<<std::endl;
	std::cout<<std::endl;

	std::chrono::duration<double> h2d_time = h2d_end - h2d_start;
	std::chrono::duration<double> krnl_time = krnl_end - krnl_start;

	std::cout<< "H2D Transfer: " << h2d_time.count() << " s"<<std::endl;
	std::cout<< "Kernel Exec: " << krnl_time.count() << " s"<<std::endl;
	std::cout<<std::endl;

	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>(krnl_end - krnl_start).count();
	// 1 subtraction, 1 multiplication, 2 min comparison, 1 addition for one cell. and there are num_features * num_features cells
	double gops = double(num_features) * num_features * 5 / exec_time;
	std::cout<< "Time: " << exec_time*1e-9<<", GOPS: "<<gops<<std::endl;

    return 0;
}
