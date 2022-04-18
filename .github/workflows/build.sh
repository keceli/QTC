#!/bin/sh

arch=Linux-x86_64
cmake_root=$(pwd)/cmake-"${cmake_version}"-"${arch}"
cmake_command=cmake #"${cmake_root}/bin/cmake"
ctest_command=ctest #"${cmake_root}/bin/ctest"
toolchain_file=$(pwd)/toolchain.cmake

#test1:
echo command test
if command -v ninja 
then
  echo "Ninja found"
  echo $(command -v ninja)
else
  echo "Ninja not found"  
fi

#test2
echo which test
if which ninja >/dev/null; then
    echo exists
else
    echo does not exist
fi

sudo apt update
sudo apt-get install ninja-build

#test1:
echo command test
if command -v ninja 
then
  echo "Ninja found"
  echo $(command -v ninja)
else
  echo "Ninja not found"  
fi

#test2
echo which test
if which ninja >/dev/null; then
    echo exists
else
    echo does not exist
fi


