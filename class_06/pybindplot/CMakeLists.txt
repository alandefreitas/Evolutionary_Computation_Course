#
# Example: plotting from C++ with Python's pybind11 library
#
# The reasons to use this example are the same as the matplotlibcpp
# The difference is that matplotlibcpp is very limited in the number of matplotlib functions it exposes to the user.
# pybind11 is a library that makes the integration between python and C++ as simple as possible
# This means:
#   1) we can access matplotlib from pybind11 with as little code as in matplotlibcpp
#   2) we can now access to all matplotlib functionalities
#   3) we can have direct access to other important python libraries that are immensely useful to us, such as scipy.stats for hypothesis testing
# However, note that calling hypothesis testing in Python is less efficient that the functions in C++, such as the one implemented in the folder "ttest"
#    1) If it is a) an important test, b) it is running inside our evolutionary cicle, and c) it is probably not going to change, then it is best to reimplement the test in C++, such as the one implemented in the folder "ttest". There are reference implementations of most tests easily available everywhere.
#    2) If it is a) a less important test, b) it is running outside our evolutionary cicle, or c) it is probably going to change, then it is best to use the Python function for a while because python has lots of functions available and they have been tested for a long time.

set(CMAKE_CXX_STANDARD 14)

# Pybind11 can also embed the Python interpreter into a C++ program
# https://pybind11.readthedocs.io/en/stable/advanced/embedding.html
# Everything needed for embedding is in embed.h
add_subdirectory(pybind11)

# Create the executable
add_executable(class_06_pybind_plot_example main.cpp)
# Linking to pybind11::embed does everything we need (pybind already looked for Python3)
target_link_libraries(class_06_pybind_plot_example PRIVATE pybind11::embed)

# We copy the file calc.py to the binary folder because it will used in our tests
file(COPY calc.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
