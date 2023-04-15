project(Dowker_complex_utilities)

# CGAL and Eigen3 are required for Euclidean version of Witness
if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  if(TARGET Boost::program_options)
	  add_executable ( Dowker_complex_dowker_persistence dowker_persistence.cpp )
	  target_link_libraries(Dowker_complex_dowker_persistence Boost::program_options)

    if(TARGET TBB::tbb)
	    target_link_libraries(Dowker_complex_dowker_persistence TBB::tbb)
	    target_link_libraries(Dowker_complex_dowker_persistence TBB::tbb)
    endif()

    add_test(NAME Dowker_complex_test_torus_persistence
	    COMMAND $<TARGET_FILE:Dowker_complex_dowker_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-l" "20" "-a" "0.5")

  install(TARGETS Dowker_complex_dowker_persistence DESTINATION bin)
  endif()
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
