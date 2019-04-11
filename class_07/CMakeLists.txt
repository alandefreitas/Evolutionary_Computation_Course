set(SOURCE_FILES
        main.cpp
        EA.h
        EA.hpp
        )

add_executable(class_07_parameter_control ${SOURCE_FILES})
target_link_libraries(class_07_parameter_control PUBLIC models)
target_include_directories(class_07_parameter_control PUBLIC ../class_01)