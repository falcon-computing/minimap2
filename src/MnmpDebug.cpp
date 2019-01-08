#include <string>
#include <unordered_map>
#include <mutex>

std::mutex g_dbgMutex;
std::unordered_map<std::string, void*> g_dbgVarTable;
