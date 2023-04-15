/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Joel Chacon Castillo
 *
 *    Copyright (C) 2023
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_WITNESS_COMPLEX_INTERFACE_H_
#define INCLUDE_WITNESS_COMPLEX_INTERFACE_H_

#include <gudhi/Simplex_tree.h>
#include <gudhi/Dowker_complex.h>

#include "Simplex_tree_interface.h"

#include <vector>
#include <utility>  // std::pair
#include <iostream>
#include <cstddef>

namespace Gudhi {

namespace dowker_complex {

class Dowker_complex_interface {
  using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
  using Nearest_landmark_table = std::vector<Nearest_landmark_range>;

 public:
  Dowker_complex_interface(const Nearest_landmark_table& nlt) {
    dowker_complex_ = new Dowker_complex<Nearest_landmark_table>(nlt);
  }

  ~Dowker_complex_interface() {
    delete dowker_complex_;
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double  epsilon,
                           std::size_t limit_dimension) {
    dowker_complex_->create_complex(*simplex_tree, epsilon, limit_dimension);
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree,
                           double  epsilon) {
    dowker_complex_->create_complex(*simplex_tree, epsilon);
  }

 private:
  Dowker_complex<Nearest_landmark_table>* dowker_complex_;
};

}  // namespace dowker_complex

}  // namespace Gudhi

#endif  // INCLUDE_WITNESS_COMPLEX_INTERFACE_H_

