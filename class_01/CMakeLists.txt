add_library(models
        # Knapsack problem
        knapsack.cpp
        knapsack.h
        knapsackP.cpp
        knapsackP.h
        # Any real function
        real_function.cpp
        real_function.h
        real_vector.cpp
        real_vector.h
        # Traveling Salesman Problem
        route.cpp
        route.h
        tsp.cpp
        tsp.h
        )

add_executable(class_01_modelling main.cpp)
target_link_libraries(class_01_modelling PUBLIC models)