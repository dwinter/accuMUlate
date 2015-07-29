set -e #cd $HOME/lib
ls -l $HOME/lib
ls -l $HOME/include
#rm -rf $HOME/lib
#rm -rf $HOME/include

ls -l $HOME/
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
  ./b2 install -d0 --prefix=$HOME/ 
  echo "Done, rebuild BOOST"
fi

EIGEN3_FILE="$HOME/include/eigen3/Eigen/src/Core/util/Macros.h"
if grep -Fxq "#define EIGEN_WORLD_VERSION 3" ${EIGEN3_FILE} && 
grep -Fxq "#define EIGEN_MAJOR_VERSION 2" ${EIGEN3_FILE} && 
grep -Fxq "#define EIGEN_MINOR_VERSION 4" ${EIGEN3_FILE};
then
  echo "Found eigen3!"
else
  echo "NO, rebuild eigen"
  cd $HOME
  mkdir -p build/eigen
  wget -qO- http://bitbucket.org/eigen/eigen/get/3.2.4.tar.bz2 | tar xj -C build/eigen/ --strip-components 1 
  cd build/eigen && mkdir build_dir && cd build_dir 
  cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/ && make install && cd ../../..
  echo "DONE rebuild eigen"
fi

BAMTOOL_FILE="$HOME/lib/libbamtools.so.2.3.0"

if [[ -e ${BAMTOOL_FILE} ]]; 
then 
  echo "Found bamtools"
else 
  echo "NO, rebuild bamtools"
  mkdir -p $HOME/build/
  cd $HOME/build/
  git clone git://github.com/pezmaster31/bamtools.git
  cd bamtools && mkdir build && cd build && cmake .. && make
  cd $HOME/build/bamtools/
  cp -r lib/ ~
  cp -r include/ ~
  echo "DONE rebuild bamtools"
fi

#cd $HOME/include
#ls -l $HOME/include/boost
#ls -l $HOME/include/eigen3
