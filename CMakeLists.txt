cmake_minimum_required(VERSION 2.8.4)

project(NFsim)

#find_package(XMLRPC REQUIRED abyss-server c++2)


#if (NOT XMLRPC_FOUND) 
#  message (FATAL_ERROR "XMLRPC-C not found") 
#endif (NOT XMLRPC_FOUND) 

set(CMAKE_BUILD_TYPE Release)

set(SUB_DIRS
    src/nauty24 
    src/NFutil 
    src/NFutil/MTrand 
    src/NFtest/transcription 
    src/NFtest/tlbr 
    src/NFtest/simple_system 
    src/NFtest/agentcell/cell 
    src/NFtest/agentcell 
    src/NFscheduler 
    src/NFreactions/transformations 
    src/NFreactions/reactions 
    src/NFreactions/reactantLists 
    src/NFreactions/mappings 
    src/NFreactions 
    src/NFoutput 
    src/NFinput 
    src/NFinput/TinyXML 
    src/NFinput/json 
    src/NFfunction/muParser 
    src/NFfunction 
    src/NFcore 
    src/NFcore/reactionSelector 
    src/NFcore/moleculeLists 
)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${XMLRPC_INCLUDE_DIRS} ${SUB_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR} include ${SUB_DIRS})


file(GLOB_RECURSE src_files  "${CMAKE_CURRENT_SOURCE_DIR}/src/*cpp")
file(GLOB_RECURSE c_src_files  "${CMAKE_CURRENT_SOURCE_DIR}/src/*c")


set(SRC_FILES 
    ${src_files} 
    ${c_src_files} 
)

add_executable(${PROJECT_NAME} ${SRC_FILES} )
