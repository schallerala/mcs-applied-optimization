get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_executable(${PROJECT_NAME}
        main.cc
        )
target_link_libraries(${PROJECT_NAME}
        AOPT::AOPT
        AOPT::MassSpringSystem
        )

set_target_properties(${PROJECT_NAME} PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        )


add_executable(${PROJECT_NAME}-test
        unit_tests.cc
        )

target_link_libraries(${PROJECT_NAME}-test
        AOPT::AOPT
        AOPT::MassSpringSystem
        gtest gtest_main

        )


gtest_add_tests(TARGET ${PROJECT_NAME}-test
        EXTRA_ARGS "--gtest_color=yes"
        )

set_target_properties(${PROJECT_NAME}-test PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        )



set(FN_INDEXES 0 1)
set(CONSTRAINTS_SCENARIOS 1 2)
set(GRID_SIDES 4 8)

set(ALL_TARGETS_LIST)

foreach(FN_INDEX IN LISTS FN_INDEXES)
    foreach(CONSTRAINTS_SCENARIO IN LISTS CONSTRAINTS_SCENARIOS)
        foreach(GRID_SIDE IN LISTS GRID_SIDES)
            SET(CUSTOM_EXECUTION_TARGET ${PROJECT_NAME}-all-${FN_INDEX}-${CONSTRAINTS_SCENARIO}-${GRID_SIDE})
            add_custom_target(${CUSTOM_EXECUTION_TARGET}
                    ${PROJECT_NAME} ${FN_INDEX} ${CONSTRAINTS_SCENARIO} ${GRID_SIDE} ${GRID_SIDE} 10000 execution_${FN_INDEX}_${CONSTRAINTS_SCENARIO}_${GRID_SIDE}_ | tee output_${CUSTOM_EXECUTION_TARGET}.log
                    DEPENDS ${PROJECT_NAME})
            list(APPEND ALL_TARGETS_LIST ${CUSTOM_EXECUTION_TARGET})
        endforeach()
    endforeach()
endforeach()

add_custom_target(${PROJECT_NAME}-all
        DEPENDS ${ALL_TARGETS_LIST})
