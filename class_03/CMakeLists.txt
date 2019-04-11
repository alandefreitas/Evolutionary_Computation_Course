set(SOURCE_FILES
        # main
        main.cpp
        # evolutionary algorithm
        selection_EA.h
        selection_EA.hpp
        )

add_executable(class_03_selection ${SOURCE_FILES})
target_link_libraries(class_03_selection PUBLIC models)
target_include_directories(class_03_selection PUBLIC ../class_01)