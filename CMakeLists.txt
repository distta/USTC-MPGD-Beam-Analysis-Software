cmake_minimum_required(VERSION 3.10)
project(BeamAnalysis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 查找必要的包
find_package(PkgConfig REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(ROOT REQUIRED)

# 查找nlohmann/json
find_package(nlohmann_json REQUIRED)

# 包含目录
include_directories(include)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})

# 源文件
set(SOURCES
    src/main.cpp
    src/Analysis.cpp
    src/Detector.cpp
    src/Geometry.cpp
    src/Processor.cpp
)

# 创建可执行文件
add_executable(BeamAnalysis ${SOURCES})

# 链接库
target_link_libraries(BeamAnalysis nlohmann_json::nlohmann_json)
target_link_libraries(BeamAnalysis ${ROOT_LIBRARIES} ROOT::Minuit)

# 设置编译选项
target_compile_options(BeamAnalysis PRIVATE -Wall -Wextra)
