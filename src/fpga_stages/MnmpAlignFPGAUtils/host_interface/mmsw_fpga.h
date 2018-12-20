/*
 * =====================================================================================
 *
 *       Filename:  mmsw_fpga.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/20/2018 04:03:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiayi Sheng (jys), jysheng@falcon-computing.com
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef MMSW_FPGA_H
#define MMSW_FPGA_H
// #include <CL/opencl.h>
// #include <CL/cl_ext.h>
#include "host_types.h"
// #include "mmsw.h"
#include <vector>

#define MAX_REF_SIZE (0x80000000)
#define MAX_N_SEQ 128
#define MAX_BATCH_SIZE 4096
#define INPUT_UNIT_SIZE 33464
#define OUTPUT_UNIT_SIZE 1208
#define OPT_SIZE 168
#define IDX_SIZE 64
// using namespace std;

typedef struct {
    uint64_t offset;
    uint32_t len;
    uint32_t rid;
} ref_seq_info;

typedef struct {
    uint32_t n_seq;
    uint64_t seqs[MAX_N_SEQ * 2];//even index is offset, odd index is rid | len
    uint32_t S[MAX_REF_SIZE / 4];
} mm_idx_fpga;

typedef struct {
    align_input data[MAX_BATCH_SIZE];
} fpga_inputs;

typedef struct {
    align_output data[MAX_BATCH_SIZE];
} fpga_outputs;


void mmsw_fpga_compute(char* btsm, const mm_mapopt_t* opt, const mm_idx_t* mi, std::vector<align_input>& inputs, std::vector<align_output>& outputs);

bool init_FPGA(char* bitstream);
bool init_FPGA_buffer();


void mmsw_fpga_init(const mm_mapopt_t* opt, char* btsm, const mm_idx_t* mi);
void mmsw_fpga_release();

void ref_conversion(const mm_idx_t* mi, mm_idx_fpga* idx);

#endif