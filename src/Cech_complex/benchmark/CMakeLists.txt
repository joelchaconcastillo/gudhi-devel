project(Cech_complex_benchmark)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 5.0.1)
  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  if(TARGET Boost::filesystem)
    add_executable(cech_complex_benchmark cech_complex_benchmark.cpp)
    target_link_libraries(cech_complex_benchmark Boost::filesystem)
  
    if(TARGET TBB::tbb)
      target_link_libraries(cech_complex_benchmark TBB::tbb)
    endif()
  endif()
endif()
