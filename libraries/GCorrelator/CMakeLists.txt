add_library(GCORRELATOR SHARED
  GCorrelator.cxx
)

target_include_directories(GCORRELATOR
  PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(GCORRELATOR
  PUBLIC
    GBCS
    GHISTOGRAMER
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
)
