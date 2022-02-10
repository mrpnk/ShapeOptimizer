# CMakeList.txt: CMake-Projekt für "ShapeOptimizer". Schließen Sie die Quelle ein, und definieren Sie
# projektspezifische Logik hier.
#
cmake_minimum_required (VERSION 3.9)
project (ShapeOptimizer LANGUAGES CXX)


# Create the target
add_executable (ShapeOptimizer "main.cpp" "hydro/utils.h"  "ShapeOptimizer.h" "hydro/integrate.h")

set_property(TARGET ShapeOptimizer PROPERTY CXX_STANDARD 20)


find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND AND false)
    message("OpenMP was found!")
    target_link_libraries(ShapeOptimizer PUBLIC OpenMP::OpenMP_CXX)
endif()



