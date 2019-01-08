/*
 * =====================================================================================
 *
 *       Filename:  host_types.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/07/2018 04:08:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiayi Sheng (jys), jysheng@falcon-computing.com
 *   Organization:  
 *
 * =====================================================================================
 */


#ifndef HOST_TYPES_H
#define HOST_TYPES_H
#include "minimap.h"
#include "ksw2.h"
#include <stdint.h>

#define MAX_SEQ_LENGTH 256
#define MAX_ANKOR_NUM 2048

typedef union {
    mm_reg1_t orig;
    uint32_t serial[17]; 
} mm_reg1_t_union;

typedef struct {
    int32_t data_size; // the entire size of current struct (taking n_a into consideration)
    int32_t qlen;
    int32_t n_a;
    uint32_t dp_score[2];
    mm_reg1_t_union region[2];
    uint8_t qseq[2][MAX_SEQ_LENGTH];
    mm128_t a[MAX_ANKOR_NUM];
} align_input;

typedef struct {
    mm_reg1_t_union region[2];
    mm_extra_t p;
    uint32_t cigar[MAX_SEQ_LENGTH];
} align_output;

typedef union {
    mm_mapopt_t orig;
    uint32_t serial[41];
} mm_mapopt_t_union;



#endif

