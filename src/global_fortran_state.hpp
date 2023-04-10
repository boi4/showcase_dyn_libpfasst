#ifndef GLOBAL_FORTRAN_STATE
#define GLOBAL_FORTRAN_STATE
#include <map>
#include <string>
#include <cstring>

std::map<std::string, int> global_int_dict;
std::map<std::string, std::string> global_str_dict;

extern "C" {
    int get_global_int_c(char *key, int *res);
    int set_global_int_c(char *key, int *val);
    int get_global_str_c(char *key, char *res);
    int set_global_str_c(char *key, char *val);
}

#endif
