add_library(mockturtle INTERFACE)
target_include_directories(mockturtle SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(mockturtle INTERFACE kitty lorina parallel_hashmap percy json bill libabcesop abcresub)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
target_link_libraries(mockturtle INTERFACE stdc++fs)
endif()
