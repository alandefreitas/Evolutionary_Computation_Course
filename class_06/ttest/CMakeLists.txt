# The subfolder ttest contains an example of a t-test
# This test is useful to compare the results of algorithms
# It is also an example of how we can implement other statistical tests
# This is important because, very often, we need statistical
#     tests as subroutines of our evolutionary algorithms
# Calling a test from another programming language is not convenient in this case
add_executable(class_06_ttest_example ttest.cpp)

# This ttest depends on boost/math to work.
# If you have boost installed on your computer, you could only use the command
#     find_package(Boost)
# as usual in cmake.
# As we only use the headers of boost/math and many students don't have boost installed on their computers, we include the boost/math headers in this project

# This boost library includes all functions related to statistical distributions
# This makes it very easy to implement any statistical test for which you have the formula
add_executable(class_06_distributions_example distributions.cpp)