#include "make_save_dir_for_unix.hpp"


void nmpcsim::makeSaveDirForUnix(const std::string dir_name)
{
    int error = mkdir(dir_name.c_str(), 0755);
}