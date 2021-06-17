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

#include <phat/representations/Container_traits.h>

namespace phat {
  
template<class ColumnContainer, class DimensionContainer>
class Uniform_representation {
  
 public:

  typedef ColumnContainer Column_container;
  typedef DimensionContainer Dimension_container;

  typedef Column_container_traits<Column_container> Column_traits;
  typedef Dimension_container_traits<Dimension_container> Dimension_traits;


  typedef typename Column_traits::Column_type Column_type;

 protected:

  // We need to label them mutable, since we have no control whether
  // the representations are manipulated by the traits class
  mutable Dimension_container dims;
  mutable Column_container matrix;

  Column_traits col_traits;
  Dimension_traits dim_traits;

  thread_local_storage< column > temp_column_buffer;

 public:
  
  index _get_num_cols() const {
    return col_traits.get_size(matrix);
  }
  
  void _set_dimensions(index nr_of_rows, index nr_of_columns) {
    col_traits.resize(matrix, nr_of_columns);
    for(index idx = 0;idx < nr_of_columns;idx++) {
      col_traits.col_at(matrix,idx)._set_nr_of_rows(nr_of_rows);
      col_traits.col_at(matrix,idx).offer_thread_local_storage(&temp_column_buffer);
    }
    dim_traits.resize(dims, nr_of_columns);
  }

  dimension _get_dim( index idx ) const {
    return dim_traits.dim_at(dims, idx );
  }

  void _set_dim( index idx, dimension dim ) { 
    dim_traits.dim_at(dims,idx) = dim; 
  }

  void _get_col( index idx, column& col  ) const { 
    col_traits.col_at(matrix, idx)._get_col( col ); 
  }

  void _set_col( index idx, const column& col  ) { 
    //col_traits.col_at(matrix, idx)._set_col( col );
    matrix[idx]._set_col(col);
  }
  
  bool _is_empty( index idx ) const {
    return col_traits.col_at(matrix, idx)._is_empty();
  }

  index _get_max_index( index idx ) const { 
    return col_traits.col_at(matrix, idx)._get_max_index(); 
  }

  void _remove_max( index idx ) { 
    return col_traits.col_at(matrix, idx)._remove_max(); 
  }

  void _add_to( index source, index target ) { 
    Column_type& source_col = col_traits.col_at(matrix, source);
    Column_type& target_col = col_traits.col_at(matrix, target);
    target_col._add_to( source_col ); 
  }

  void _clear( index idx ) { 
    col_traits.col_at(matrix, idx)._clear(); 
  }
        
  void _finalize( index idx ) { 
    col_traits.col_at(matrix, idx)._finalize(); 
  }
        
  void _sync() {}

};
 
}
