PROJECT(Vision3D)


CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{IMAGINEPP_ROOT}/CMake)
FIND_PACKAGE(Imagine)

SET(CMAKE_BUILD_TYPE release)  #debug or release

ImagineAddExecutable(PhotoStereo maxflow/graph.cpp parse_params.cpp polygone.cpp camera.cpp carac.cpp lsqr.cpp stereo.cpp main.cpp)
ImagineUseModules(PhotoStereo Images LinAlg)
