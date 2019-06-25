//
// Save data of the numerical simulation.
//

#ifndef MAKE_SAVE_DIR_FOR_WIN_H
#define MAKE_SAVE_DIR_FOR_WIN_H


#include <iostream>
#include <fstream>
#include <string>
#include <direct.h>


namespace nmpcsim {
void makeSaveDirForWin(const std::string dir_name);
}


#endif