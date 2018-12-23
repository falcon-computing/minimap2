/*
 * =====================================================================================
 *
 *       Filename:  mmsw_fpga.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/20/2018 05:24:34 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiayi Sheng (jys), jysheng@falcon-computing.com
 *   Organization:  
 *
 * =====================================================================================
 */

#include "mmsw_fpga.h"
#include "host_types.h"
#include "mmsw.h"
#include <stdlib.h>
#include <string.h>
#include "CL/opencl.h"

bool initialized = false;
fpga_inputs* in = NULL;
fpga_outputs* out = NULL;
mm_idx_fpga* reference = NULL; 
mm_mapopt_t_union* opt_union = NULL;

cl_platform_id platform_id;         // platform id
cl_device_id device_id;             // compute device id 
cl_context context;                 // compute context
cl_command_queue command;          // compute command queue
cl_program program;                 // compute program
cl_mem _opt_buffer;
cl_mem _ref_buffer;
cl_mem _mmi_buffer;
cl_mem _in_buffer;
cl_mem _out_buffer;
cl_kernel _mmsw_kernel;
char btsm_filename[200];


void ref_conversion(const mm_idx_t* mi, mm_idx_fpga* idx) {
    uint64_t sum_len = 0;
    idx->n_seq = mi->n_seq;
    for (int i = 0; i < mi->n_seq; ++i) {
    mm_idx_seq_t *s = &mi->seq[i];
        idx->seqs[2 * i] = s->offset;
        idx->seqs[2 * i + 1] = ((uint64_t)i << 32) | (uint64_t)s->len;
        sum_len += s->len;
    }
    memcpy(idx->S, mi->S, (sum_len + 7) / 8 * 4);
}

int load_file_to_memory(const char *filename, char **result) { 
    size_t size = 0;
    FILE *f = fopen(filename, "rb");
    if (f == NULL) { 
        printf("ERROR : Kernel binary %s not exist!\n", filename);
        *result = NULL;
        return -1; // -1 means file opening fail 
    } 
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    *result = (char *)malloc(size+1);
    if ((int)size != (int)fread(*result, sizeof(char), size, f)) { 
        free(*result);
        return -2; // -2 means file reading fail 
    } 
    fclose(f);
    (*result)[size] = 0;
    return size;
}

bool init_FPGA(char* bitstream){
     //get platform info
    char cl_platform_vendor[1001];
    char cl_platform_name[1001];
    cl_platform_vendor[0] = 0;
    cl_platform_name[0] = 0;
    int err;
    
    err = clGetPlatformIDs(1, &platform_id, NULL);
    if (err != CL_SUCCESS) {
        printf("Warning: Failed to find an OpenCL platform!Code %i\n", err);
        return false;
    }
    printf("Successfully create platform\n");
    err = clGetPlatformInfo(platform_id, CL_PLATFORM_VENDOR, 1000, (void *)cl_platform_vendor, NULL);
    if (err != CL_SUCCESS) {
        printf("Warning: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed! Code %i\n", err);
        return false;
    }
    printf("CL_PLATFORM_VENDOR %s\n", cl_platform_vendor);
    err = clGetPlatformInfo(platform_id, CL_PLATFORM_NAME, 1000, (void *)cl_platform_name, NULL);
    if (err != CL_SUCCESS) {
        printf("Warning: clGetPlatformInfo(CL_PLATFORM_NAME) failed! Code %i\n", err);
        return false;
    }
    printf("CL_PLATFORM_NAME %s\n", cl_platform_name);

    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, 1, &device_id, NULL);
    if (err != CL_SUCCESS) {
        printf("Warning: Failed to create a device group! Code %i\n", err);
        return false;
    }
    printf("Successfully create device\n");

    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    if (!context) {
        printf("Warning: Failed to create a compute context! Code %i\n", err);
        return false;
    }
    printf("Successfully create context \n");
    command = clCreateCommandQueue(context, device_id, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err);
    if (!command) {
        printf("Warning: Failed to create a command queue commands! Code %i\n", err);
        return false;
    }
    printf("Successfully create command queue \n");
    unsigned char *kernelbinary;

    int n_i = 0;
    n_i = load_file_to_memory(bitstream, (char **) &kernelbinary);
    if (n_i < 0) {
        printf("Warning : failed to load kernel from binary: %s\n", bitstream);
        return false;
    }
    printf("Successfully load kernel from binary: %s\n", bitstream);
    
    int status;
    size_t n = n_i;
    program = clCreateProgramWithBinary(context, 1, &device_id, &n, (const unsigned char **) &kernelbinary, &status, &err);
    if ((!program) || (err!=CL_SUCCESS)) {
        printf("Warning: Failed to create compute program from binary! Code %d\n", err);
        return false;
    }
    printf("Success to create compute program from binary! \n");

    // Build the program executable
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t len;
        char buffer[2048];
        printf("Warning: Failed to build program executable!\n");
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        return false;
    }
    printf("Sucess to build program executable!\n");
    _mmsw_kernel = clCreateKernel(program, "mmsw_kernel", &err);
    if((!_mmsw_kernel) || err != CL_SUCCESS){
        printf("Warning: Failed to create compute kernel for mmsw core\n");
        return false;
    }
    return true;
}

bool init_FPGA_buffer(){
    unsigned bankID[4] = {XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1, \
        XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};
    cl_mem_ext_ptr_t bank_ext[4];
    bank_ext[0].flags = bankID[0];
    bank_ext[0].param = 0;
    bank_ext[0].obj = nullptr;
    cl_int err = 0;

    bool success = true;
    
    _mmsw_kernel = clCreateKernel(program, "mmsw_kernel", &err);
    if (!_mmsw_kernel) {
        success = false;
        printf("failed to create _mmsw_kernel\n, err code is %d\n", err);
    }
    else 
        printf("create kernel\n");
    
    _opt_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX, sizeof(mm_mapopt_t_union), &bank_ext[0], &err);
    
    if (!_opt_buffer) {
        success = false;
        printf("failed to create _opt_buffer, err code is %d\n", err);
    }
    else 
        printf("create opt_buffer\n");
    
    int cl_status; 
    cl_status = clEnqueueWriteBuffer(command, _opt_buffer, CL_TRUE,\
                    0, sizeof(mm_mapopt_t_union), opt_union->serial, 0, NULL, NULL);
    if (cl_status != CL_SUCCESS) {
        success = false;
        printf("failed to write opt buffer, err code is %d\n", err);
    }
    _ref_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX, sizeof(uint64_t) * reference->n_seq * 2, &bank_ext[0], &err);

    if (!_ref_buffer) {
        success = false;
        printf("failed to create _ref_buffer, err code is %d\n", err);
    }
    else 
        printf("create ref_buffer\n");

    cl_status = clEnqueueWriteBuffer(command, _ref_buffer, CL_TRUE,\
                    0, sizeof(uint64_t) * reference->n_seq * 2, reference->seqs , 0, NULL, NULL);
    if (cl_status != CL_SUCCESS) {
        success = false;
        printf("failed to write ref buffer, err code is %d\n", err);
    }
 

    _mmi_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX, MAX_REF_SIZE, &bank_ext[0], &err);
 
    if (!_mmi_buffer) {
        success = false;
        printf("failed to create _mmi_buffer, err code is %d\n", err);
    }
    else 
        printf("create mmi_buffer\n");

    cl_status = clEnqueueWriteBuffer(command, _mmi_buffer, CL_TRUE,\
                    0, MAX_REF_SIZE, reference->S , 0, NULL, NULL);
    if (cl_status != CL_SUCCESS) {
        success = false;
        printf("failed to write mmi buffer, err code is %d\n", err);
    }
 
    _in_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX, sizeof(fpga_inputs), &bank_ext[0], &err); 

    if (!_in_buffer) {
        success = false;
        printf("failed to create _in_buffer, err code is %d\n", err);
    }
    else 
        printf("create in_buffer\n");
 
    _out_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX, sizeof(fpga_outputs), &bank_ext[0], &err); 
 
    if (!_out_buffer) {
        success = false;
        printf("failed to create _out_buffer, err code is %d\n", err);
    }
    else 
        printf("create out_buffer\n");
    
    return success;
}

void mmsw_fpga_init(const mm_mapopt_t* opt, char* btsm, const mm_idx_t* mi) {
    if(!initialized) {
        if (!(opt_union = (mm_mapopt_t_union*)malloc(sizeof(mm_mapopt_t_union)))) {
            printf("Failed to allocate for option union\n");
            exit(-1);
        }
        if (!(reference = (mm_idx_fpga*)malloc(sizeof(mm_idx_fpga)))) {
            printf("failed to allocate for ref\n");
            exit(-1);
        }
        if (!(in = (fpga_inputs*)malloc(sizeof(fpga_inputs)))) {
            printf("failed to allocate for fpga_inputs\n");
            exit(-1);
        }
        if (!(out = (fpga_outputs*)malloc(sizeof(fpga_outputs)))) {
            printf("failed to allocate for fpga_outputs\n");
            exit(-1);
        }
        opt_union->orig = *opt;
        ref_conversion(mi, reference);
#ifndef CSIM
        init_FPGA(btsm);
        init_FPGA_buffer();
#endif
        initialized = true;
    }
}


void mmsw_fpga_destroy() {
    if (initialized) {
        if (!reference) 
            free(reference);
        if (!in)
            free(in);
        if (!out)
            free(out);
        initialized = false;
    }
}

void inputs2fpga(std::vector<align_input>& inputs, fpga_inputs* in, int seg_offset, int cur_batch_size) {
    for (int i = seg_offset; i < seg_offset + cur_batch_size; i++) {
        inputs[i].data_size = sizeof(align_input) - 16 * MAX_ANKOR_NUM + 16 * inputs[i].n_a;
        in->data[i - seg_offset] = inputs[i];
    }
#ifndef CSIM
    cl_int err;
    int cl_status; 
    cl_event event;
    cl_status = clEnqueueWriteBuffer(command, _in_buffer, CL_TRUE,\
                    0, sizeof(align_input) * cur_batch_size, (uint32_t*)in , 0, NULL, &event);
    if (cl_status != CL_SUCCESS) {
        printf("failed to write in buffer\n");
        exit(-1);
    }
    clWaitForEvents(1, &event);
#endif
}

void fpga2outputs(std::vector<align_output>& outputs, fpga_outputs* out, int seg_offset, int cur_batch_size) {
#ifndef CSIM
    int cl_status; 
    cl_event event;
    cl_status = clEnqueueReadBuffer(command, _out_buffer, CL_TRUE,\
                    0, sizeof(align_output) * cur_batch_size, (uint32_t*)out , 0, NULL, &event);
    if (cl_status != CL_SUCCESS) {
        printf("failed to read out buffer\n");
        exit(-1);
    }
    clWaitForEvents(1, &event);
#endif
    for (int i = seg_offset; i < seg_offset + cur_batch_size; i++) {
        outputs.push_back(out->data[i - seg_offset]);
    }

}



void mmsw_fpga_compute(char* btsm, const mm_mapopt_t* opt, const mm_idx_t* mi, std::vector<align_input>& inputs, std::vector<align_output>& outputs) {
    mmsw_fpga_init(opt, btsm, mi);
    int seg_num = inputs.size() / MAX_BATCH_SIZE + (inputs.size() % MAX_BATCH_SIZE != 0);
    int seg_offset = 0;
    for (int i = 0; i < seg_num; i++) {
        int cur_batch_size = seg_offset + MAX_BATCH_SIZE < inputs.size() ? MAX_BATCH_SIZE : inputs.size() - seg_offset;
        inputs2fpga(inputs, in, seg_offset, cur_batch_size);
#ifdef CSIM
        // mmsw_kernel(opt_union->serial, reference->n_seq, reference->seqs, reference->S, (uint32_t*)in, (uint32_t*)out, cur_batch_size);
#else
        int err = 0;
        cl_event event;
        err |= clSetKernelArg(_mmsw_kernel, 0, sizeof(cl_mem), &_opt_buffer);
        err |= clSetKernelArg(_mmsw_kernel, 1, sizeof(uint32_t), &reference->n_seq);
        err |= clSetKernelArg(_mmsw_kernel, 2, sizeof(cl_mem), &_ref_buffer);
        err |= clSetKernelArg(_mmsw_kernel, 3, sizeof(cl_mem), &_mmi_buffer);
        err |= clSetKernelArg(_mmsw_kernel, 4, sizeof(cl_mem), &_in_buffer);
        err |= clSetKernelArg(_mmsw_kernel, 5, sizeof(cl_mem), &_out_buffer);
        err |= clSetKernelArg(_mmsw_kernel, 6, sizeof(int), &cur_batch_size);
        err |= clEnqueueTask(command, _mmsw_kernel, 0, NULL, &event);
        if(err != CL_SUCCESS) {
            fprintf(stderr, "Error: Failed to set kernel arguments %d\n", err);
            exit(-1);
        }
        clWaitForEvents(1, &event);
#endif
        fpga2outputs(outputs, out, seg_offset, cur_batch_size);
        seg_offset += cur_batch_size;
    }
}

void mmsw_fpga_release() {
#ifndef CSIM
    clReleaseMemObject(_opt_buffer);
    clReleaseMemObject(_in_buffer);
    clReleaseMemObject(_out_buffer);
    clReleaseMemObject(_ref_buffer);
    clReleaseMemObject(_mmi_buffer);
    clReleaseKernel(_mmsw_kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(command);
    clReleaseContext(context);   
#endif
    free(in);
    free(out);
    free(reference);
    free(opt_union);
}


