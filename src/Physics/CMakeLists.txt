add_library(physics
  pengrobinson.cpp
)

# BOOST
find_package(BOOST REQUIRED)

# Eigen
Find_package(Eigen3 REQUIRED)

target_include_directories(physics PUBLIC
  ${PROJECT_INCLUDE_DIRS}
  ~/local/pugixml/pugixml-1.14/src/
  ~/tools/boost_1_84_0
  ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(physics 
  ~/local/pugixml/build/libpugixml.a
  ${BOOST_LIBRARIES}
)
