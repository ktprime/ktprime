
cmake_minimum_required(VERSION 3.10)
project(NumberTheoryPrograms)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -O3 -march=native)


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

    target_link_libraries(${exe_name} m)

    if (${source_file} MATCHES "Fast")
        target_compile_options(${exe_name} PRIVATE -mavx2)
    endif()

    message(STATUS "Added executable: ${exe_name} from ${source_file}")
endforeach()
