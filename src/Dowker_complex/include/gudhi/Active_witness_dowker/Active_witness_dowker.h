/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s): Joel Chacon, Jose Frias, Ignacio Segovia, Yuzhu Chen, Yulia R. Gel
 *    Contact email: joel.chaconcastillo@utdallas.edu
 */

#ifndef ACTIVE_WITNESS_ACTIVE_WITNESS_H_
#define ACTIVE_WITNESS_ACTIVE_WITNESS_H_

#include <gudhi/Active_witness_dowker/Active_witness_dowker_iterator.h>
#include <list>

namespace Gudhi {

namespace dowker_complex {

  /** \class Active_witness
   *  \brief Class representing a list of nearest neighbors to a given witness.
   *  \details Every element is a pair of a landmark identifier and the squared distance to it.
  */
template< typename Id_distance_pair,
          typename INS_range >
class Active_witness_dowker{
 public:
  typedef Active_witness_dowker<Id_distance_pair, INS_range> ActiveWitness;
  typedef typename INS_range::iterator INS_iterator;
  typedef Active_witness_dowker_iterator< ActiveWitness, Id_distance_pair, INS_iterator > iterator;
  typedef typename std::list<Id_distance_pair> Table;

  Table nearest_landmark_table_;
  INS_range    search_range_;
  INS_iterator iterator_next_;
  INS_iterator iterator_end_;

  Active_witness_dowker(const INS_range& search_range)
    : search_range_(search_range), iterator_next_(search_range_.begin()), iterator_end_(search_range_.end()) {
  }

  iterator begin() {
    return iterator(this, nearest_landmark_table_.begin());
  }

  iterator end() {
    return iterator(this);
  }
};

}  // namespace witness_complex
}  // namespace Gudhi

#endif  // ACTIVE_WITNESS_ACTIVE_WITNESS_H_
