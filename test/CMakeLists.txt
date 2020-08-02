set(TEST_FILES
        test_strassen.cpp
        test_bini.cpp
        test_main.cpp
        helpers.cpp
        test_classic.cpp test_laderman.cpp test_schonhage.cpp)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable (tests EXCLUDE_FROM_ALL ${TEST_FILES})

target_link_libraries(tests fast_matrix_multiplication)
target_link_libraries(tests ${GTEST_LIBRARIES} pthread)
