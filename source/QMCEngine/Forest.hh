/******************************************************************\

Copyright 2015 Dimitrios Galanakis

This file is part of QMCSGF

QMCSGF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 3.

QMCSGF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with QMCSGF.  If not, see <http://www.gnu.org/licenses/>.

\*******************************************************************/

#ifndef __FOREST__
#define __FOREST__

#include "TSum.hh"

/*

  This class contains the trees and also functions to update them.

*/
template<typename TSumType>
class ForestType {

  typedef typename TSumType::IndexType ElementIndex;
  typedef typename TSumType::Float Float;
  typedef typename TSumType::DoubleFloat DoubleFloat;
  typedef unsigned int TreeIndex;
  struct TreeType {
    TSumType tsum[2];
  };

  const TreeIndex numtrees;
  TreeType* trees;                    // Holds the probability trees
  std::vector<TreeType*> tree_cache;  // remembers the offset of each operator

public:

  ForestType(const ElementIndex& _nterms, const TreeIndex& _ntrees) : numtrees(_ntrees), tree_cache(_nterms) {

    trees = new TreeType[numtrees];

    for (TreeIndex i = 0; i < size(); ++i) {
      trees[i].tsum[0].resize(_nterms);
      trees[i].tsum[1].resize(_nterms);
    }

  }

  /*
    size()
      returns the number of subtrees
  */
  inline TreeIndex size() const {
    return numtrees;
  }

  /*
    norm<rl>(i)
      return the normalization of offset i in direction rl
  */
  template<int rl>
  inline const DoubleFloat& norm(const TreeIndex& i) const {
    return trees[i].tsum[rl].norm();
  }

  /*
    choose<rl>(i)
      chooses one index with offset i in direction rl
  */
  template<int rl>
  inline ElementIndex choose(const TreeIndex& i) const {
    return trees[i].tsum[rl].choose();
  }

  /*
    update<rl>(o,i,v)
      changes the value correponding to index i from whatever it was
      before to the new value v. Also if the offset of the index i
      has changed, it will move it to the correct subtree.
  */
  template<int rl>
  inline void update(const TreeIndex& ioffset, const ElementIndex& i, const Float& v) {
    TreeType* i_tree = tree_cache[i];
    TreeType* f_tree = &trees[ioffset];
    f_tree->tsum[ rl].update(i, v);
    if (i_tree != f_tree) {
      Float jme = i_tree->tsum[!rl].element(i);
      i_tree->tsum[ rl].update(i, 0);
      i_tree->tsum[!rl].update(i, 0);
      f_tree->tsum[!rl].update(i, jme);
      tree_cache[i] = f_tree;
    }
  }

  /*
    set<rl>(offset,i,v)
      This is meant to be used for initialization. It sets the value
      and offset for index i assuming they were not known before.
  */
  template<int rl>
  inline void set(const TreeIndex& ioffset, const ElementIndex& i, const Float& v) {
    TreeType* tree = &trees[ioffset];
    tree_cache[i] = tree;
    tree->tsum[rl].update(i, v);
  }

  ~ForestType() {
    delete [] trees;
  }
};

#endif // __FOREST__
