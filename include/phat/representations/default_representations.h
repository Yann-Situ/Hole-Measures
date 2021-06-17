/*  Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Michael Kerber, Jan Reininghaus

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>. */

#pragma once

#include <phat/helpers/misc.h>
#include <phat/representations/Uniform_representation.h>
#include <phat/representations/Pivot_representation.h>

#include <phat/representations/sparse_pivot_column.h>
#include <phat/representations/full_pivot_column.h>
#include <phat/representations/heap_pivot_column.h>
#include <phat/representations/bit_tree_pivot_column.h>
#include <phat/representations/vector_column_rep.h>
#include <phat/representations/list_column_rep.h>
#include <phat/representations/set_column_rep.h>
#include <phat/representations/heap_column_rep.h>


namespace phat {

 // typedef Uniform_representation< std::unordered_map<index,phat::vector_column_rep>, std::unordered_map<index,index> > hash_vector;
 typedef Uniform_representation< std::vector<phat::vector_column_rep>, std::vector<index> > vector_vector;
 typedef Uniform_representation< std::vector<phat::set_column_rep>, std::vector<index> >    vector_set;
 typedef Uniform_representation< std::vector<phat::heap_column_rep>, std::vector<index> >   vector_heap;
 typedef Uniform_representation< std::vector<phat::list_column_rep>, std::vector<index> >   vector_list;
 typedef Pivot_representation< vector_vector, full_column >     full_pivot_column;
 typedef Pivot_representation< vector_vector, heap_column >     heap_pivot_column;
 typedef Pivot_representation< vector_vector, sparse_column >   sparse_pivot_column;
 typedef Pivot_representation< vector_vector, bit_tree_column > bit_tree_pivot_column;

}
