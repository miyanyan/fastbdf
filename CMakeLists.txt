cmake_minimum_required(VERSION 3.25)
project(fastbdf VERSION 0.1.0 LANGUAGES CXX)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_library(fastbdf INTERFACE)
add_library(fastbdf::fastbdf ALIAS fastbdf)
target_compile_features(fastbdf INTERFACE cxx_std_20)
target_include_directories(fastbdf INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

install(TARGETS fastbdf EXPORT fastbdfTargets)
install(DIRECTORY include/fastbdf DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
set(fastbdf_cmake_dir "${CMAKE_INSTALL_LIBDIR}/cmake/fastbdf")
install(EXPORT fastbdfTargets NAMESPACE fastbdf:: FILE fastbdfConfig.cmake
    DESTINATION "${fastbdf_cmake_dir}")
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/fastbdfConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/fastbdfConfigVersion.cmake"
    DESTINATION "${fastbdf_cmake_dir}")

option(FASTBDF_TESTS "Build the native tests" OFF)
option(FASTBDF_INTEROP "Build EDFlib interoperability tests and benchmark" OFF)
option(FASTBDF_MNE "Build optional MNE-Python comparison and benchmark" OFF)
if(FASTBDF_TESTS OR FASTBDF_INTEROP OR FASTBDF_MNE)
    enable_testing()
endif()
if(FASTBDF_TESTS)
    add_executable(fastbdf_tests tests/native.cpp tests/odr.cpp)
    target_link_libraries(fastbdf_tests PRIVATE fastbdf::fastbdf)
    add_test(NAME native COMMAND fastbdf_tests)
    foreach(area IN ITEMS reader writer)
        add_executable(fastbdf_${area}_tests tests/${area}.cpp)
        target_link_libraries(fastbdf_${area}_tests PRIVATE fastbdf::fastbdf)
        add_test(NAME ${area} COMMAND fastbdf_${area}_tests)
    endforeach()
endif()
if(FASTBDF_INTEROP)
    find_package(unofficial-EDFlib CONFIG REQUIRED)
    foreach(name IN ITEMS interop bench)
        add_executable(fastbdf_${name} tests/${name}.cpp)
        target_link_libraries(fastbdf_${name} PRIVATE fastbdf::fastbdf EDFlib)
        add_test(NAME ${name} COMMAND fastbdf_${name})
    endforeach()
    target_sources(fastbdf_interop PRIVATE tests/write_codes.cpp)
    set_tests_properties(bench PROPERTIES LABELS benchmark)
endif()

if(FASTBDF_MNE)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    foreach(target IN ITEMS mne_probe mne_bench)
        add_executable(fastbdf_${target} tests/${target}.cpp)
        target_link_libraries(fastbdf_${target} PRIVATE fastbdf::fastbdf)
    endforeach()
    add_test(NAME mne_compare COMMAND Python3::Interpreter
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/mne_compare.py" $<TARGET_FILE:fastbdf_mne_probe>)
    add_test(NAME python_binding COMMAND Python3::Interpreter
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python_binding.py" $<TARGET_FILE:fastbdf_mne_probe>)
    add_test(NAME mne_bench COMMAND Python3::Interpreter
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/mne_bench.py" $<TARGET_FILE:fastbdf_mne_bench>)
    set_tests_properties(mne_bench PROPERTIES LABELS benchmark)
endif()
