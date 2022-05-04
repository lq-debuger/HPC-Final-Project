
include_directories(/work/mae-liq1/lib/petsc-3.16.6/include)
find_library(PETSC petsc /work/mae-liq1/lib/petsc-3.16.6/lib)
add_executable(demo demo_petsc.c)
target_link_libraries(demo ${PETSC})