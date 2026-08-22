add_library(FILLHISTOGRAM SHARED
  FillHistogram.cxx
)

target_include_directories(FILLHISTOGRAM
  PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(FILLHISTOGRAM 
  PUBLIC 
   GBCS 
   GHISTOGRAMER 
   ${ROOT_LIBRARIES})
