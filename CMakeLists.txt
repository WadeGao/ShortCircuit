cmake_minimum_required(VERSION 3.4.1)

project(ShortCircuit)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++17 -fopenmp -m64 -g -m64 -DMKL_ILP64")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -std=c++17 -Wl,--no-as-needed -fopenmp -m64 -DMKL_ILP64")

# 指定根目录
set(ROOT ${PROJECT_SOURCE_DIR})

# 指定头文件路径
set(INC_DIR ${ROOT}/include/)

# 可执行文件输出目录
set(EXECUTABLE_OUTPUT_PATH ${ROOT}/bin)

#链接库，mysqlclient为库名

link_libraries(mysqlclient)
link_libraries(mkl_intel_ilp64)
link_libraries(mkl_intel_thread)
link_libraries(mkl_core)
link_libraries(iomp5)
link_libraries(pthread)
link_libraries(m)
link_libraries(dl)

#指定头文件目录
include_directories(${INC_DIR})

# 添加当前目录下所有的源文件
aux_source_directory(${ROOT}/source/ ALL_SRCS)

add_executable(${PROJECT_NAME} ${ALL_SRCS})


