/*
 * =====================================================================================
 *
 *       Filename:  mmsw_baseline.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/13/2018 06:02:39 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiayi Sheng (jys), jysheng@falcon-computing.com
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef MMSW_BASELINE_H
#define MMSW_BASELINE_H
#include "minimap.h"
#include "host_types.h"
#include "mmpriv.h"
#include <time.h>
#include <vector>

// timespec diff_time(timespec& start, timespec& end);
void mmsw_baseline_compute(const mm_mapopt_t* opt, const mm_idx_t* mi, std::vector<align_input>& inputs, std::vector<align_output>& outputs, double& cells, double& sw_time);
#endif
