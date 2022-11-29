## Custon Python interface examples
Here we show how to build Python interfaces for the custom OCP definitions.   
Specifically, we consider an cartpole that containts external reference represented as a shared pointer of a struct.
Its C++ documentation is found at https://ohtsukalab.github.io/autogenu-jupyter/classcgmres_1_1_o_c_p__cartpole_external_reference.html.

First, install `cgmres` C++ library and `autogenu` Python module. 
In the project root directory of `autogenu-jupyter`, do 
```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=YOUR_INSTALL_DESTINATION
make install 
```
and
```
python3 -m pip install setuptools
python3 -m pip install .
```

To build the Python interfaces, run the following command in this directory
```
mkdir build 
cd build
cmake ..
make -j8
cd ..
```

Then install the python interfaces by running the following command in this directory
```
python3 install.py 
```

Then you can run the closed-loop simulation of the cartpole with the external reference via
```
python3 cartpole_external_reference.py
```