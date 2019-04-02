set(SOURCE_FILES
        # main
        main.cpp
        # knapsack problem
        knapsack.h
        knapsack.cpp
        knapsackP.h
        knapsackP.cpp
        # real function
        real_function.h
        real_function.cpp
        real_vector.h
        real_vector.cpp
        # travelling salesman problem
        route.h
        route.cpp
        tsp.h
        tsp.cpp
        # evolutionary algorithm
        EA.h
        EA.hpp
        )
add_executable(class_02_representation ${SOURCE_FILES})