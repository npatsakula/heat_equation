cmake_minimum_required(VERSION 3.16)

project(heat_equation VERSION 1.0 LANGUAGES CXX)
add_executable(heat_equation main.cpp math.hpp cli.hpp)

target_link_libraries(heat_equation fmt)
set_target_properties(heat_equation PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        INTERPROCEDURAL_OPTIMIZATION TRUE
)

target_compile_options(heat_equation PRIVATE
        $<$<CONFIG:DEBUG>: -Wall -Wextra -pedantic>
        $<$<CONFIG:Release>: -DNDEBUG -O2 -funroll-loops -ftree-vectorize>
)

target_compile_options(heat_equation PRIVATE $<$<CXX_COMPILER_ID:GNU>:
        $<$<CONFIG:DEBUG>: -ggdb -Wformat=2 -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align>>
)

target_compile_options(heat_equation PRIVATE $<$<CXX_COMPILER_ID:Clang>:
        $<$<CONFIG:DEBUG>:-Weverything -Wno-c++98-compat-pedantic -Wno-missing-prototypes>>
)
