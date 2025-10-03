add_library(common
  common_string_operations.cpp
)

target_include_directories(common PUBLIC 
  ${PROJECT_SOURCE_DIR}
)
