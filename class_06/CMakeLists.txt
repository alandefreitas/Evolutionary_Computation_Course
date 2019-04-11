# We could find boost:
#if (NOT DEFINED Boost_FOUND)
#    find_package(Boost QUIET REQUIRED)
#endif()
#if (${Boost_FOUND})
#set(SOURCE_FILES ttest.cpp)
#add_executable(class_06_ttest_example ${SOURCE_FILES})
#    target_include_directories(class_06_ttest_example PUBLIC ${Boost_INCLUDE_DIR})
#endif()

# Or we can just use the headers
add_executable(class_06_ttest_example ttest.cpp)

if (NOT DEFINED ${PythonLibs_FOUND})
    find_package(PythonLibs)
endif()

if (${PythonLibs_FOUND})
    message("${PYTHON_INCLUDE_DIRS}")
    add_executable(class_06_plotting_example plotting.cpp matplotlibcpp.h)
    target_include_directories(class_06_plotting_example PUBLIC ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(class_06_plotting_example ${PYTHON_LIBRARIES})
    if(NOT ${PYTHON_INCLUDE_DIRS}/numpy/arrayobject.h)
        target_compile_definitions(class_06_plotting_example PUBLIC -DWITHOUT_NUMPY)
    endif()
endif()