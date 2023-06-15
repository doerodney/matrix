message("reading /src/matrix/CMakeLists.txt")

set(SOURCES matrix.c)

add_library(matrix STATIC ${SOURCES})
target_include_directories(matrix PRIVATE ../include) 

