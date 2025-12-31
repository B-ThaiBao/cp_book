cmake_minimum_required(VERSION 3.20)
project(cp_book LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)
enable_testing()

include(FetchContent)
FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG v3.6.0
)
FetchContent_MakeAvailable(Catch2)

# module catch_discover_tests
list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(Catch)

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")

foreach(src ${TEST_SOURCES})
	file(RELATIVE_PATH test "${CMAKE_CURRENT_SOURCE_DIR}/tests" "${src}")

	get_filename_component(dir "${test}" DIRECTORY)
	get_filename_component(stem "${test}" NAME_WE)

	if(dir STREQUAL "")
		set(tgt "test_${stem}")
		set(out_dir "${CMAKE_BINARY_DIR}/tests")
	else()
		string(REPLACE "/" "_" dir_id "${dir}")
		set(tgt "test_${dir_id}_${stem}")
		set(out_dir "${CMAKE_BINARY_DIR}/tests/${dir}")
	endif()

	add_executable(${tgt} "${src}")

	target_include_directories(${tgt} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(${tgt} PRIVATE Catch2::Catch2WithMain)

	set_target_properties(${tgt} PROPERTIES OUTPUT_NAME "${stem}" RUNTIME_OUTPUT_DIRECTORY "${out_dir}")
	catch_discover_tests(${tgt})
endforeach()
