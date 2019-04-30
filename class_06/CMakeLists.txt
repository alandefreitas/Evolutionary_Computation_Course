#
# Example: t-test
#
add_executable(class_06_ttest_example ttest.cpp)

#
# Example: plotting from C++ with Python's matlibplot library
#
# Look for Python 2 Interpreter and NumPy Library
find_package(Python2 QUIET COMPONENTS Interpreter Development NumPy)

# If we found the interpreter
if (${Python2_Interpreter_FOUND})
    message("Python2_INCLUDE_DIR=${Python2_INCLUDE_DIR}")
    message("Python2_INCLUDE_DIRS=${Python2_INCLUDE_DIRS}")
    message("Python2_LIBRARIES=${Python2_LIBRARIES}")
    # Create our application
    add_executable(class_06_plotting_example plotting.cpp matplotlibcpp.h)
    # Give our application access to Python include directories
    target_include_directories(class_06_plotting_example PUBLIC ${Python2_INCLUDE_DIRS})
    # Link our application to the Python2 libraries we found
    target_link_libraries(class_06_plotting_example PUBLIC ${Python2_LIBRARIES})
    # If we did not find the NumPy library
    if(NOT Python2_NumPy_FOUND)
        # Deactivate numpy from the C++ library with the WITHOUT_NUMPY compile definition
        target_compile_definitions(class_06_plotting_example PUBLIC -DWITHOUT_NUMPY)
    endif()
else ()
    message(WARNING "We could not find the Python 2 interpreter")
    if (${UNIX})
        if (${APPLE})
            message(WARNING "Try \"brew install python2\"")
        else ()
            message(WARNING "Try \"apt install python2.7 python-pip\"")
        endif()
    else()
        message(WARNING "Download Python 2.7 from https://www.python.org/")
    endif()
endif ()