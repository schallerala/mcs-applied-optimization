cmake_minimum_required(VERSION 3.11)

PROJECT(AOPT_EXERCISES)

include(FetchContent)

FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

if (WIN32)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build")

    # avoid linking errors, cf https://stackoverflow.com/questions/12540970/how-to-make-gtest-build-mdd-instead-of-mtd-by-default-using-cmake
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
else ()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
endif ()


FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

enable_testing()
include(GoogleTest)


set(AOPT_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")


if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# add library
add_library (AOPT INTERFACE)

add_library (AOPT::AOPT ALIAS AOPT)

#Set target properties
set(EIGEN3_INCLUDE_DIR $<BUILD_INTERFACE:${AOPT_EXTERNAL}/eigen>)

target_include_directories(AOPT
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        ${EIGEN3_INCLUDE_DIR}
)

add_subdirectory(EigenTutorial)
add_subdirectory(GridSearch)
add_subdirectory(CsvExporter)
add_subdirectory(ConvexityTests)
add_subdirectory(MassSpringSystem)
add_subdirectory(MassSpringProblemEvaluation)
add_subdirectory(OptimalityChecker)
add_subdirectory(GradientDescent)
add_subdirectory(NewtonMethods)

# FN_INDEXES: without or with spring length
set(FN_INDEXES 0 1)
# 0: dense, 1: sparse matrix
# FIXME, missing config 2 for Newton execution
set(SPARSE_MATRIX_CONFIGS 0 1)
set(GRID_SIDES 5 10 20)
set(EXECUTABLES GradientDescent NewtonMethods)
set(MAX_ITER 10000)

SET(RES_DIR results_output)

set(ALL_TARGETS_LIST)

add_custom_command(OUTPUT ${RES_DIR}
        COMMAND mkdir ${RES_DIR}
        COMMENT "Create result output folder")

foreach(GRID_SIDE IN LISTS GRID_SIDES)
    foreach(FN_INDEX IN LISTS FN_INDEXES)
        foreach(SPARSE_MATRIX IN LISTS SPARSE_MATRIX_CONFIGS)
            foreach(EXECUTABLE IN LISTS EXECUTABLES)
                SET(CUSTOM_EXECUTION_TARGET all-${EXECUTABLE}-${FN_INDEX}_${SPARSE_MATRIX}_${GRID_SIDE})
                SET(LOG_OUTPUT_FILE ${RES_DIR}/output_${CUSTOM_EXECUTION_TARGET}.log)
                add_custom_command(OUTPUT ${LOG_OUTPUT_FILE}
                        COMMAND echo "Arguments: ${FN_INDEX} ${SPARSE_MATRIX} ${GRID_SIDE}" > ${LOG_OUTPUT_FILE}
                        DEPENDS ${RES_DIR})
                add_custom_target(${CUSTOM_EXECUTION_TARGET}
                        ${EXECUTABLE} ${FN_INDEX} ${SPARSE_MATRIX} ${GRID_SIDE} ${GRID_SIDE} ${MAX_ITER} ${RES_DIR}/execution_${EXECUTABLE}_${FN_INDEX}_${SPARSE_MATRIX}_${GRID_SIDE}__ 2>&1 | tee -a ${LOG_OUTPUT_FILE}
                        DEPENDS ${EXECUTABLE} ${LOG_OUTPUT_FILE})
                list(APPEND ALL_TARGETS_LIST ${CUSTOM_EXECUTION_TARGET})
            endforeach()
        endforeach()
    endforeach()
endforeach()

# Add a few more run to use function index 2 which can be set for Newton method
foreach(GRID_SIDE IN LISTS GRID_SIDES)
    # FN_INDEX = 2: f with length with positive local hessian
    SET(FN_INDEX 2)
    foreach(SPARSE_MATRIX IN LISTS SPARSE_MATRIX_CONFIGS)
        SET(EXECUTABLE NewtonMethods)
        SET(CUSTOM_EXECUTION_TARGET all-${EXECUTABLE}-${FN_INDEX}_${SPARSE_MATRIX}_${GRID_SIDE})
        SET(LOG_OUTPUT_FILE ${RES_DIR}/output_${CUSTOM_EXECUTION_TARGET}.log)
        add_custom_command(OUTPUT ${LOG_OUTPUT_FILE}
                COMMAND echo "Arguments: ${FN_INDEX} ${SPARSE_MATRIX} ${GRID_SIDE}" > ${LOG_OUTPUT_FILE}
                DEPENDS ${RES_DIR})
        add_custom_target(${CUSTOM_EXECUTION_TARGET}
                ${EXECUTABLE} ${FN_INDEX} ${SPARSE_MATRIX} ${GRID_SIDE} ${GRID_SIDE} ${MAX_ITER} ${RES_DIR}/execution_${EXECUTABLE}_${FN_INDEX}_${SPARSE_MATRIX}_${GRID_SIDE}__ 2>&1 | tee -a ${LOG_OUTPUT_FILE}
                DEPENDS ${EXECUTABLE} ${LOG_OUTPUT_FILE})
        list(APPEND ALL_TARGETS_LIST ${CUSTOM_EXECUTION_TARGET})
    endforeach()
endforeach()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/parse_logs.awk
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_custom_command(OUTPUT result.tsv
        COMMAND cat ${RES_DIR}/*.log | awk -f parse_logs.awk > result.tsv
        DEPENDS all-ex3
        COMMENT "Parse logs to create result file for import")

add_custom_target(all-ex3
        DEPENDS ${ALL_TARGETS_LIST} result.tsv)
