#include "make_save_dir_for_win.hpp"


void nmpcsim::makeSaveDirForWin(const std::string dir_name)
{
    int error = mkdir(dir_name.c_str());
}