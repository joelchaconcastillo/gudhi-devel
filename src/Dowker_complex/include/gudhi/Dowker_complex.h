/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOWKER_COMPLEX_H_
#define DOWKER_COMPLEX_H_

#include <gudhi/Active_witness_dowker/Active_witness_dowker.h>
#include <gudhi/Dowker_complex/all_faces_in.h>

#include <utility>
#include <vector>
#include <list>
#include <limits>

namespace Gudhi {

namespace dowker_complex {

/**
 * \private
 * \class Dowker_complex
 * \brief Constructs (weak) dowker complex for a given table of nearest landmarks with respect to dowkeres.
 * \ingroup dowker_complex
 *
 * \tparam Nearest_landmark_table_ needs to be a range of a range of pairs of nearest landmarks and distances.
 *         The class Nearest_landmark_table_::value_type must be a copiable range.
 *         The range of pairs must admit a member type 'iterator'. The dereference type 
 *         of the pair range iterator needs to be 'std::pair<std::size_t, double>'.
*/
template< class Nearest_landmark_table_ >
class Dowker_complex {
 private:
  typedef typename Nearest_landmark_table_::value_type               Nearest_landmark_range;
  typedef std::size_t                                                Dowker_id;
  typedef std::size_t                                                Landmark_id;
  typedef std::pair<Landmark_id, double>                             Id_distance_pair;
  typedef Active_witness_dowker<Id_distance_pair, Nearest_landmark_range>   ActiveWitness;
  typedef std::list< ActiveWitness >                                 ActiveWitnessList;
  typedef std::vector< Landmark_id >                                 typeVectorVertex;
  typedef std::vector<Nearest_landmark_range>                        Nearest_landmark_table_internal;
  typedef Landmark_id Vertex_handle;

 protected:
  Nearest_landmark_table_internal              nearest_landmark_table_;

 public:
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* @name Constructor
   */

  //@{

  Dowker_complex() {
  }

  /**
   *  \brief Initializes member variables before constructing simplicial complex.
   *  \details Records nearest landmark table.
   *  @param[in] nearest_landmark_table needs to be a range (one entry per dowker)
   *         of sorted ranges of pairs of nearest landmarks and distances.
   *         The class Nearest_landmark_table_::value_type must be a copiable range.
   *         The range of pairs must admit a member type 'iterator'. The dereference type 
   *         of the pair range iterator needs to be 'std::pair<std::size_t, double>'
   *         where the first element is the index of the landmark, and the second its
   *         (squared) distance to the dowker.
   */

  Dowker_complex(Nearest_landmark_table_ const & nearest_landmark_table)
    : nearest_landmark_table_(std::begin(nearest_landmark_table), std::end(nearest_landmark_table)) {
  }

  /** \brief Outputs the (weak) dowker complex of relaxation 'epsilon'
   *         in a simplicial complex data structure.
   *  \details The function returns true if the construction is successful and false otherwise.
   *  @param[out] complex Simplicial complex data structure compatible which is a model of
   *              SimplicialComplexForDowker concept.
   *  @param[in] epsilon Maximal squared relaxation parameter.
   *  @param[in] limit_dimension Represents the maximal dimension of the simplicial complex
   *         (default value = no limit).
   */
  template < typename SimplicialComplexForDowker >
  bool create_complex(SimplicialComplexForDowker& complex,
                      double  epsilon,
                      std::size_t limit_dimension = std::numeric_limits<std::size_t>::max()) const {
    if (complex.num_vertices() > 0) {
      std::cerr << "Dowker complex cannot create complex - complex is not empty.\n";
      return false;
    }
    if (epsilon< 0) {
      std::cerr << "Dowker complex cannot create complex - squared relaxation parameter must be non-negative.\n";
      return false;
    }
    std::cout<<"here!!!\n";
    ActiveWitnessList active_witness;
    Landmark_id k = 0; /* current dimension in iterative construction */
    for (auto&& w : nearest_landmark_table_)
      active_witness.emplace_back(w);
    while (!active_witness.empty() && k <= limit_dimension) {
      typename ActiveWitnessList::iterator aw_it = active_witness.begin();
      std::vector<Landmark_id> simplex;
      simplex.reserve(k+1);
      while (aw_it != active_witness.end()) {
        bool ok = add_all_faces_of_dimension(k,
                                             epsilon,
                                             aw_it->begin(),
                                             simplex,
                                             complex,
                                             aw_it->end());
        assert(simplex.empty());
        if (!ok)
          active_witness.erase(aw_it++);  // First increase the iterator and then erase the previous element
        else
          aw_it++;
      }
      k++;
    }
    return true;
  }

  //@}

 private:
  /** \brief Adds recursively all the faces of a certain dimension dim dowkered by the same dowker.
   * Iterator is needed to know until how far we can take landmarks to form simplexes.
   * simplex is the prefix of the simplexes to insert.
   * The output value indicates if the dowker rests active or not.
   */
  template < typename SimplicialComplexForDowker >
  bool add_all_faces_of_dimension(int dim,
                                  double epsilon,
                                  typename ActiveWitness::iterator curr_l,
                                  std::vector<Landmark_id>& simplex,
                                  SimplicialComplexForDowker& sc,
                                  typename ActiveWitness::iterator end) const {
    if (curr_l == end)
      return false;
    bool will_be_active = false;
    typename ActiveWitness::iterator l_it = curr_l;
    if (dim > 0) {
      for (; l_it != end && l_it->second <= epsilon; ++l_it) {
        simplex.push_back(l_it->first);
        if (sc.find(simplex) != sc.null_simplex()) {
          typename ActiveWitness::iterator next_it = l_it;
          will_be_active = add_all_faces_of_dimension(dim-1, epsilon, ++next_it, simplex, sc, end) || will_be_active;
        }
        assert(!simplex.empty());
        simplex.pop_back();
      }
    } else if (dim == 0) {
      for (; l_it != end && l_it->second <= epsilon; ++l_it) {
        simplex.push_back(l_it->first);
        double filtration_value = 0;
        filtration_value = l_it->second;
        if (all_faces_in(simplex, &filtration_value, sc)) {
          will_be_active = true;
          sc.insert_simplex(simplex, filtration_value);
        }
        assert(!simplex.empty());
        simplex.pop_back();
      }
    }
    return will_be_active;
  }
};

}  // namespace dowker_complex

}  // namespace Gudhi

#endif  // DOWKER_COMPLEX_H_
