set -e #cd $HOME/lib
mkdir -p $HOME/local
mkdir -p $HOME/local/lib
mkdir -p $HOME/local/include

#ls -l $HOME/lib
#ls -l $HOME/include
#rm -rf $HOME/lib
#rm -rf $HOME/include

#ls -l $HOME/
#ls -l $HOME/lib/boost
#define BOOST_VERSION 105800
if grep -Fxq "#define BOOST_VERSION 105800" $HOME/include/boost/version.hpp; 
then
  echo "Found Boost!"; 
else 
  echo "NO, rebuild BOOST"; 
  cd $HOME
  mkdir -p  build/boost_install
  wget -qO- http://downloads.sourceforge.net/project/boost/boost/1.58.0/boost_1_58_0.tar.bz2 | tar xj -C build/boost_install/
  pwd && ls && ls build && ls build/boost_install && ls build/boost_install/boost_1_58_0
  cd build/boost_install/boost_1_58_0 && ./bootstrap.sh --with-libraries=program_options,thread,system 
  ./b2 install -d0 --prefix=$HOME/local 
  echo "Done, rebuild BOOST"
fi

if [ -e "${HOME}/local/lib/libgtest.a" ]
then
    echo "Found Gtest"
else
    if [ -e "/usr/src/gtest/src/gtest.cc" ]
    then
        echo "Gtest is not present, rebuilding"
        mkdir -p ./tmp
        cp -a /usr/src/gtest/ ./tmp
        mkdir ./tmp/gtest/build
        cd ./tmp/gtest/build
        cmake ..
        make
        mv libg* ${HOME}/local/lib
        cd ${HOME}/local
        echo "DONE rebuilding gtest"
     else
         echo "Gtest not cahced, or available, skipping"
  fi
fi


EIGEN3_FILE="$HOME/local/include/eigen3/Eigen/src/Core/util/Macros.h"
if grep -Fxq "#define EIGEN_WORLD_VERSION 3" ${EIGEN3_FILE} && 
grep -Fxq "#define EIGEN_MAJOR_VERSION 2" ${EIGEN3_FILE} && 
grep -Fxq "#define EIGEN_MINOR_VERSION 4" ${EIGEN3_FILE};
then
  echo "Found eigen3!"
else
  echo "NO, rebuild eigen"
  cd $HOME/local
  mkdir -p build/eigen
  wget -qO- http://bitbucket.org/eigen/eigen/get/3.2.4.tar.bz2 | tar xj -C build/eigen/ --strip-components 1 
  cd build/eigen && mkdir build_dir && cd build_dir 
  cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/local && make install && cd ../../..
  echo "DONE rebuild eigen"
fi

BAMTOOL_FILE="$HOME/local/lib/bamtools/libbamtools.a"

if [[ -e ${BAMTOOL_FILE} ]]; 
then 
  echo "Found bamtools"
else 
  echo "NO, rebuild bamtools"
  mkdir -p $HOME/build/
  cd $HOME/build/
  git clone git://github.com/pezmaster31/bamtools.git
  cd bamtools && mkdir build && cd build
  cmake -DCMAKE_INSTALL_PREFIX=$HOME/local .. &&  make && make install
  echo "DONE rebuild bamtools"
  ls -l $HOME/local/include/bamtools
  ls -l $HOME/local/include/bamtools/shared
fi

#cd $HOME/include
#ls -l $HOME/include/boost
#ls -l $HOME/include/eigen3
