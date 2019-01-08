#ifndef MNMP_DEBUG_H
#define MNMP_DEBUG_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <mutex>

extern std::mutex g_dbgMutex;
extern std::unordered_map<std::string, void*> g_dbgVarTable;

template<typename T>
void RECORD_VAR(std::string const i_filename, std::string const i_scope, std::string const i_name, T const i_value)
{
  std::lock_guard<std::mutex> l_lock(g_dbgMutex);
  std::string const l_queryName = i_filename+"_"+i_scope+"_"+i_name;
  if (g_dbgVarTable.count(l_queryName)) {
    T* l_min = (T*)(g_dbgVarTable[l_queryName]);
    T* l_max = l_min+1;
    T* l_total = l_min+2;
    long* l_count = (long*)(l_min+3);

    if (i_value < *l_min)
      *l_min = i_value;
    if (i_value > *l_max)
      *l_max = i_value;
    *l_total = *l_total + i_value;
    *l_count = *l_count + 1;
  }
  else {
    T *l_min = (T*)malloc(3*sizeof(T)+sizeof(long));
    T* l_max = l_min+1;
    T* l_total = l_min+2;
    long* l_count = (long*)(l_min+3);

    *l_min = i_value;
    *l_max = i_value;
    *l_total = i_value;
    *l_count = 1;
    g_dbgVarTable[l_queryName] = l_min;
  }
}

template<typename T>
void DISPLAY_VAR(std::string const i_filename, std::string const i_scope, std::string const i_name)
{
  std::string const l_queryName = i_filename+"_"+i_scope+"_"+i_name;
  if (!g_dbgVarTable.count(l_queryName))
    return;
  T* l_min = (T*)g_dbgVarTable[l_queryName];
  T* l_max = l_min+1;
  T* l_total = l_min+2;
  long* l_count = (long*)(l_min+3);
  std::cerr << "[" << i_filename << ":" << i_scope << "] "
             << "\"" << i_name << "\" - "
             << "min: " << *l_min << ", "
             << "avg: " << *l_total/(double)(*l_count) << ", "
             << "max: " << *l_max << "." << std::endl;
}

#endif
