add_library(physics
  PhysicsManager.cpp
  pengrobinson.cpp
)

# BOOST
find_package(BOOST REQUIRED)

# Eigen
Find_package(Eigen3 REQUIRED)

# pugixml
find_package(pugixml REQUIRED)

target_include_directories(physics PUBLIC
  ${PROJECT_INCLUDE_DIRS}
  ~/local/pugixml/pugixml-1.14/src/
  ~/tools/boost_1_84_0
  ${EIGEN3_INCLUDE_DIR}
  /usr/local/include
)

target_link_libraries(physics 
  ${BOOST_LIBRARIES}
  pugixml
  math
  common
  IO
)
