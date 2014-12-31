@if NOT EXIST build mkdir build
@cd build
@if EXIST CMakeCache.txt del CMakeCache.txt
cmake.exe .. -DCMAKE_INSTALL_PREFIX=D:\cpplibrary\descriptors -DBUILD_SHARED=OFF -DBABEL_INCLUDE_DIR=D:\openbabel-2.3.2\include\openbabel-2.0 -DBABEL_LIB_DIR=D:\openbabel-2.3.2\bin
@cd ..
