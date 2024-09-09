cmake_minimum_required(VERSION 3.16) # ITP server standard
project(QAOA-QTG C CXX)

set(CMAKE_C_STANDARD 11) # ITP server standard
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/executables)

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

set(SRC ${CMAKE_SOURCE_DIR}/src)
set(INCLUDE ${CMAKE_SOURCE_DIR}/include)

add_subdirectory(extern/nlopt)

add_executable(main main.c
        ${SRC}/knapsack.c
        ${SRC}/stategen.c
        ${SRC}/syslinks.c
        ${SRC}/combowrp.c
        ${SRC}/combo.c
        ${SRC}/qtg_count.c
        ${SRC}/general_count.c
        ${SRC}/copula_count.c
        ${SRC}/qaoa.c
)

add_executable(landscape landscape.c
        ${SRC}/knapsack.c
        ${SRC}/stategen.c
        ${SRC}/syslinks.c
        ${SRC}/combowrp.c
        ${SRC}/combo.c
        ${SRC}/qtg_count.c
        ${SRC}/general_count.c
        ${SRC}/copula_count.c
        ${SRC}/qaoa.c
)

add_executable(test unit_test.c
        ${SRC}/knapsack.c
        ${SRC}/stategen.c
        ${SRC}/syslinks.c
        ${SRC}/combowrp.c
        ${SRC}/combo.c
        ${SRC}/qtg_count.c
        ${SRC}/general_count.c
        ${SRC}/copula_count.c
        ${SRC}/qaoa.c
)

add_executable(generate ${SRC}/generator.cpp)

target_link_libraries(main PRIVATE nlopt m)
target_link_libraries(landscape PRIVATE nlopt m)
target_link_libraries(test PRIVATE nlopt m)
target_include_directories(main PRIVATE ${INCLUDE} extern/nlopt)
target_include_directories(landscape PRIVATE ${INCLUDE} extern/nlopt)
target_include_directories(test PRIVATE ${INCLUDE} extern/nlopt)

