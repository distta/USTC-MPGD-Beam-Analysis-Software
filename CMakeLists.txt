cmake_minimum_required(VERSION 3.14)
project(BeamAnalysis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------
# 查找 ROOT
# ------------------------------
find_package(ROOT REQUIRED)
if(NOT ROOT_FOUND)
    message(FATAL_ERROR "ROOT not found!")
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
if(NOT Eigen3_FOUND)
    message(FATAL_ERROR "Eigen3 not found!")
endif()

set(SOURCES
    src/main.cpp
    src/Detector/Detector.cpp
    src/Detector/Planar.cpp
    src/Algorithm.cpp
    src/AnalysisEngine.cpp
)

add_executable(BeamAnalysis ${SOURCES})

target_include_directories(BeamAnalysis PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR} 
)

target_link_libraries(BeamAnalysis
    ${ROOT_LIBRARIES}
    stdc++fs
)

# 头文件和 LinkDef
set(DI_HDR ${CMAKE_CURRENT_SOURCE_DIR}/include/DataModel.h)
set(LINKDEF  ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)

ROOT_GENERATE_DICTIONARY(BeamAnalysis_dict ${CMAKE_CURRENT_SOURCE_DIR}/include/DataModel.h 
                        MODULE BeamAnalysis 
                        LINKDEF ${LINKDEF})
