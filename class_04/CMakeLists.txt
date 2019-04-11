set(SOURCE_FILES
        main.cpp
        island_EA.h
        island_EA.hpp
        )

add_executable(class_04_population_management ${SOURCE_FILES})
target_link_libraries(class_04_population_management PUBLIC models)
target_include_directories(class_04_population_management PUBLIC ../class_01)