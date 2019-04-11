set(SOURCE_FILES
        # main
        main.cpp
        # evolutionary algorithm
        naive_EA.h
        naive_EA.hpp
    )

add_executable(class_02_representation ${SOURCE_FILES})
target_link_libraries(class_02_representation PUBLIC models)
target_include_directories(class_02_representation PUBLIC ../class_01)