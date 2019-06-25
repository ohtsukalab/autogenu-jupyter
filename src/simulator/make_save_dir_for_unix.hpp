//
// Save data of the numerical simulation.
//

#ifndef MAKE_SAVE_DIR_FOR_UNIX_H
#define MAKE_SAVE_DIR_FOR_UNIX_H


#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>


namespace nmpcsim {
void makeSaveDirForUnix(const std::string dir_name);
}


#endif