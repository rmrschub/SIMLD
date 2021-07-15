cmake_minimum_required(VERSION 3.0.0)
project(SIMLD VERSION 0.1.0)

INCLUDE_DIRECTORIES(
    ${PROJECT_SOURCE_DIR}/include/
    ${PROJECT_SOURCE_DIR}/extern/marisa-trie/include/
    ${PROJECT_SOURCE_DIR}/extern/TurboPFor-Integer-Compression/
    /Users/resc01-admin/VisualStudio/vcpkg/installed/x64-osx/include/
)

LINK_DIRECTORIES(
    ${PROJECT_SOURCE_DIR}/extern/marisa-trie/lib/
    ${PROJECT_SOURCE_DIR}/extern/TurboPFor-Integer-Compression/
    /Users/resc01-admin/VisualStudio/vcpkg/installed/x64-osx/include/
)

add_subdirectory(apps)