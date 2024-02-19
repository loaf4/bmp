add_library(BMPLibrary STATIC
    ${CMAKE_CURRENT_LIST_DIR}/bmp.cpp)

target_include_directories(BMPLibrary PUBLIC
    ${CMAKE_CURRENT_LIST_DIR})
