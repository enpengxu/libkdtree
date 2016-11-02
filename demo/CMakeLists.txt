cmake_minimum_required(VERSION 2.8)

project(kdtree_demo)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} ../include)

add_definitions(-DKD_FLOAT -DKD_DIM=3)
add_executable(kdtree_demo  src/main.cpp src/loader.cpp
  ../src/kdtree.c ../src/kdtree.cpp ../src/csv.cpp)
target_link_libraries(kdtree_demo ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
