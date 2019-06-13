#
# Example: plotting from C++ with Python's matlibplot library
#
# This is important because it is not easy to understand the performance of our algorithms without being able to look at graphs
# We could create a GUI that interactively plots the results with a library such as Qt (https://www.qcustomplot.com)
# The problem with that is that we would write a lot of code for something that:
#    a) will not be useful in the future (because people most often use C++ code as libraries for other software),
#    b) is not relevant to the kernel algorithm
#    c) would create a heavy dependency with Qt or any other library we use for this GUI
# In this case, our solution is to use a python's matploglib to plot results.
# matplotlibcpp is a library that wraps matplotlib in a C++ interface.
# We usually use that for experiments and, once we are done, we turn the plotting option off.
# From that point, users of our AE library can use it and plot the results however they want./
# Note that, unless you plots are *very* simple, it's more convenient to pybind (the next example) because:
#   1) pybind finds all relevant python packages
#   2) pybind seems to work on more operating systems
#   3) pybind lets you use all plotting functions
#   4) pybind lets you use all python functions, including functions for statistics
#   5) pybind is as simple as matplotlib
#   6) pybind is more well maintained that matplotlib

# Look for Python 2 Interpreter and NumPy Library
find_package(Python2 QUIET COMPONENTS Interpreter Development NumPy)

# If we found the interpreter
if (${Python2_Interpreter_FOUND})
    # Create our application
    add_executable(class_06_matplotlibcpp_example plotting.cpp matplotlibcpp.h)
    # Give our application access to Python include directories
    target_include_directories(class_06_matplotlibcpp_example PUBLIC ${Python2_INCLUDE_DIRS} ${PYTHON_SITE_PACKAGES} ${Python2_STDLIB})
    # Link our application to the Python2 libraries we found
    target_link_libraries(class_06_matplotlibcpp_example PUBLIC ${Python2_LIBRARIES} )
    # If we did not find the NumPy library
    if(NOT Python2_NumPy_FOUND)
        message("Python2_NumPy wasn't found")
        # Deactivate numpy from the C++ header with the WITHOUT_NUMPY compile definition
        target_compile_definitions(class_06_matplotlibcpp_example PUBLIC -DWITHOUT_NUMPY)
    else()
        # link NumPy
        target_link_libraries(class_06_matplotlibcpp_example PUBLIC Python2::NumPy)
    endif()
else ()
    # If we didn't find Python
    message(WARNING "We could not find the Python 2 interpreter")
    # If it's linux or Mac
    if (${UNIX})
        # If it's mac
        if (${APPLE})
            # Show command to install python on Mac
            message(WARNING "Try \"brew install python2\"")
        else ()
            # Show command to install python on Linux
            message(WARNING "Try \"apt install python2.7 python-pip\"")
        endif()
    else()
        # Show how to install python for windows
        message(WARNING "Download Python 2.7 from https://www.python.org/")
    endif()
endif ()