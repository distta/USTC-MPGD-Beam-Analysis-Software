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

find_package(Threads REQUIRED)

set(SOURCES
    src/main.cpp
    src/Event/DetectorFrame.cpp
    src/Event/EventDisplay.cpp
    src/Detector/Detector.cpp
    src/Detector/Planar.cpp
    src/Detector/PlanarPad.cpp
    src/Detector/Cylinder.cpp
    src/Detector/DetectorFactory.cpp
    src/Algorithm/AlgorithmFactory.cpp
    src/Algorithm/Analyzer/WaveformProcessor.cpp
    src/Algorithm/Oscilloscope/OscilloscopeDataProcessor.cpp
    src/Algorithm/Analyzer/ClusterBuilder.cpp
    src/Algorithm/Analyzer/ClusterReconstructor.cpp
    src/Algorithm/AnalysisUtils.cpp
    src/Algorithm/Track/TrackReconstruction.cpp
    src/Algorithm/Track/TrackerAlignment.cpp
    src/Algorithm/Track/TrackPerformance.cpp
    src/Input/ChannelMap.cpp
    src/Input/ConverterFactory.cpp
    src/Input/BTAPVDatConverter.cpp
    src/Input/APV25SRSConverter.cpp
    src/Script/Base/RawDataParser.cpp
    src/Script/Base/ScriptFactory.cpp
    src/Script/ExampleScript.cpp
    src/Script/TrackAnalysisScript.cpp
    src/Script/OscilloscopeAnalysisScript.cpp
    src/Script/TimeResolutionScript.cpp
    src/Script/DUTAnalysisScript.cpp
    src/Script/DUTEfficiencyAnalysis.cpp
    src/Script/PadDUTAnalysisScript.cpp
)

add_executable(BeamAnalysis ${SOURCES})

target_include_directories(BeamAnalysis PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Detector
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Event
    ${CMAKE_CURRENT_SOURCE_DIR}/include/EventDisplay
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Input
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Script
    ${ROOT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR} 
)

target_link_libraries(BeamAnalysis
    ${ROOT_LIBRARIES}
    Threads::Threads
    stdc++fs
)


# 头文件和 LinkDef
set(DI_HDR ${CMAKE_CURRENT_SOURCE_DIR}/include/Event/DataModel.h)
set(LINKDEF  ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)

ROOT_GENERATE_DICTIONARY(BeamAnalysis_dict ${CMAKE_CURRENT_SOURCE_DIR}/include/Event/DataModel.h 
                        MODULE BeamAnalysis 
                        LINKDEF ${LINKDEF})
