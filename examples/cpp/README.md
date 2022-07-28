## C++ examples
This shows how to use `cgmres` as a header-only library.
First, install `cgmres` library in the project root directory as
```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=YOUR_INSTALL_DESTINATION
make install 
```

Then you can run C++ examples, e.g., via
```
cd cartpole
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=YOUR_INSTALL_DESTINATION -DCMAKE_BUILD_TYPE=Release -DVECTORIZE=ON
make 
./cartpole
```

