/*
 * =====================================================================================
 *
 *       Filename:  mmsw_orig.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/09/2018 10:50:40 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiayi Sheng (jys), jysheng@falcon-computing.com
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef MMSW_ORIG_H
#define MMSW_ORIG_H
#include "minimap.h"
#include "host_types.h"
#include <vector>

void mmsw_orig_compute(const mm_mapopt_t* opt, const mm_idx_t* mi, std::vector<align_input>& inputs, std::vector<align_output>& outputs);
#endif

