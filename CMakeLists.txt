cmake_minimum_required(VERSION 2.6)
project(ulcm)

set(CMAKE_CXX_STANDARD 11)

add_executable(ulcm
    core.cpp

    ulcm.cpp
    ufgm.cpp

    ods_p.cpp
    ods_gs.cpp

    model_s.cpp
    model_ns.cpp

    main.cpp
)
