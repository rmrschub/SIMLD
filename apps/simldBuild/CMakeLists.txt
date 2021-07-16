
project(simldBuild VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra -g -mbmi2 -mavx -mavx2 -madx -O3 -mtune='intel' -m64 -mlzcnt")

add_executable(simldBuild main.cpp NTriplesParser.cpp)

target_compile_features(simldBuild PRIVATE cxx_std_17)

target_link_libraries(simldBuild PRIVATE
    libserd.a
    marisa
    bitpack.o
    bitpack_avx2.o
    bitpack_sse.o
    vp4c.o
    vp4c_avx2.o
    vp4c_sse.o
    vp4d.o
    vp4d_avx2.o
    vp4d_sse.o
    transpose.o
    transpose_sse.o
    transpose_avx2.o
    bitunpack.o
    bitunpack_sse.o
    bitunpack_avx2.o
    fp.o
    vint.o
    trlec.o
    bitutil.o
    trled.o
    eliasfano.o
    v8.o
    vsimple.o
)