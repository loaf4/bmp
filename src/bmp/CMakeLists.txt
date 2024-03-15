add_library(BMPLibrary STATIC
    ${CMAKE_CURRENT_LIST_DIR}/bmp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/analysis.cpp)

target_include_directories(BMPLibrary PUBLIC
    ${CMAKE_CURRENT_LIST_DIR})
