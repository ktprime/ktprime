cmake_minimum_required(VERSION 3.10)
project(NumberTheoryPrograms)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(WIN32)
    add_compile_options(
        /W4
        /permissive-
        /Zc:__cplusplus
    )

    set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /RTC1")

    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /GL /Gy")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
else()
    add_compile_options(-Wall -Wextra)

    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address")

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -mtune=native")
endif()

set(SOURCE_FILES
    PrimeNumber.cpp
    Ktprime.cpp
    TwinPrime.cpp
    FastGn.cpp
    Goldbach.cpp
    SophieGermain.cpp
    FastTn.cpp
    topK.cpp
)

foreach(source_file ${SOURCE_FILES})
    get_filename_component(exe_name ${source_file} NAME_WE)
    add_executable(${exe_name} ${source_file})
    if(NOT WIN32)
        target_link_libraries(${exe_name} m)
    endif()
endforeach()
