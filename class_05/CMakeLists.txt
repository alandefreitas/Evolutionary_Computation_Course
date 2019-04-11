set(SOURCE_FILES
        main.cpp
        parameter_control_EA.h
        parameter_control_EA.hpp
        )

add_executable(class_05_canonic_algorithms ${SOURCE_FILES})
target_link_libraries(class_05_canonic_algorithms PUBLIC models)
target_include_directories(class_05_canonic_algorithms PUBLIC ../class_01)