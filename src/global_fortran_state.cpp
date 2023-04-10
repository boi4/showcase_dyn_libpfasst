#include "global_fortran_state.hpp"

extern "C" {
    int get_global_int_c(char *key, int *res) {
        *res = global_int_dict[std::string(key)];
        return 0;
    }
    int set_global_int_c(char *key, int *val) {
        global_int_dict[std::string(key)] = *val;
        return 0;
    }
    int get_global_str_c(char *key, char *res) {
        strcpy(res, global_str_dict[std::string(key)].c_str());
        return 0;
    }
    int set_global_str_c(char *key, char *val) {
        global_str_dict[std::string(key)] = std::string(val);
        return 0;
    }
}
