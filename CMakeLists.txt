cmake_minimum_required (VERSION 2.6)

# create debug symbols
#set(CMAKE_BUILD_TYPE DEBUG)

# create "compile_commands.json"
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# some c++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -std=c++1z -O2")

file(GLOB src_files
    "src/main.cpp"
)

project (gen_srt_table)
add_executable(gen_srt_table ${src_files})
