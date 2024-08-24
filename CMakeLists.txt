cmake_minimum_required(VERSION 3.28)
project(skyline_matrix CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# 手动指定Eigen库路径
set(CMAKE_PREFIX_PATH "/usr/local/include/eigen3")

# 或者使用以下方式指定Eigen3_DIR
# set(Eigen3_DIR "/usr/local/include/eigen3/share/eigen3/cmake")

# 查找Eigen库
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 设置包含目录
include_directories(include)  # 添加 include 目录为头文件搜索路径

add_executable(skyline_matrix test/main.cpp
        include/Matrix_encryption.cpp
        include/Matrix_encryption.h
        include/Skyline_search.cpp
        include/Skyline_search.h)

# 链接Eigen库到可执行文件
target_link_libraries(skyline_matrix PRIVATE Eigen3::Eigen)
target_link_libraries(skyline_matrix PRIVATE -pthread)