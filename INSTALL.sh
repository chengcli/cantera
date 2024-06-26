CXX=g++
CC=gcc
cxx_flags="-std=c++17"
prefix=${HOME}/opt/
python_package=n
f90_interface=n
system_eigen=n
extra_inc_dirs="/opt/include /usr/include/yaml-cpp"  # Include both directories
system_blas_lapack=n
boost_inc_dir=$(pwd)/ext/cliboost
system_sundials=n
system_yamlcpp=n

scons build CXX=${CXX} CC=${CC} cxx_flags="${cxx_flags}" prefix=${prefix} \
 python_package=${python_package} f90_interface=${f90_interface} \
 system_eigen=${system_eigen} extra_inc_dirs="${extra_inc_dirs}" \
 system_blas_lapack=${system_blas_lapack} boost_inc_dir=${boost_inc_dir} system_sundials=${system_sundials} \
 system_yamlcpp=${system_yamlcpp} -j8

scons install
