project(Dowker_complex_tests)

include(GUDHI_boost_test)

add_executable ( Dowker_complex_test_simple_dowker_complex test_simple_dowker_complex.cpp )
if(TARGET TBB::tbb)
	target_link_libraries(Dowker_complex_test_simple_dowker_complex TBB::tbb)
endif()

gudhi_add_boost_test(Dowker_complex_test_simple_dowker_complex)

# CGAL and Eigen3 are required for Euclidean version of Witness
if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
	add_executable ( Dowker_complex_test_euclidean_simple_dowker_complex test_euclidean_simple_dowker_complex.cpp )
  if(TARGET TBB::tbb)
	  target_link_libraries(Dowker_complex_test_euclidean_simple_dowker_complex TBB::tbb)
  endif()

  gudhi_add_boost_test(Dowker_complex_test_euclidean_simple_dowker_complex)
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
