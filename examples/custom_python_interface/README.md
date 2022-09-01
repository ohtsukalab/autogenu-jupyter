## Custon Python interface examples
Here we show how to build Python interfaces for the custom OCP definitions.   
Specifically, we consider an cartpole that containts external reference represented as a shared pointer of a struct.
Its C++ documentation is found at https://mayataka.github.io/autogenu-jupyter/classcgmres_1_1_o_c_p__cartpole_external_reference.html.

First, install `cgmres` library in the project root directory as
```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=YOUR_INSTALL_DESTINATION
make install 
```

To build the Python interfaces, run the following command
```
mkdir build 
cd build
cmake ..
make -j8
cd ..
```

Then install the python interfaces by
```
python3 install.py cartpole_external_reference
```

Then you can run the closed-loop simulation of the cartpole with the external reference via
```
python3 cartpole_external_reference.py
```