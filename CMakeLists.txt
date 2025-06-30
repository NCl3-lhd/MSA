cmake_minimum_required(VERSION 3.10)

project(main)

# 设置头文件路径
include_directories(${PROJECT_SOURCE_DIR}/include)

# 查找 src 目录下的所有 .cpp 文件
file(GLOB SOURCES "src/*.cpp")

# 添加可执行文件
add_executable(MSA ${SOURCES})