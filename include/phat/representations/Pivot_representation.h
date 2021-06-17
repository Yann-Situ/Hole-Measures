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

namespace phat{

template<typename BaseRepresentation, typename PivotColumn> 
  class Pivot_representation : public BaseRepresentation {

 public:
  
  typedef BaseRepresentation Base;
  typedef PivotColumn Pivot_column;
  
 protected:
  
  // For parallization purposes, it could be more than one full column
  mutable thread_local_storage< Pivot_column > pivot_cols;
  mutable thread_local_storage< index > idx_of_pivot_cols;

  Pivot_column& get_pivot_col() const {
    return pivot_cols();
  }

  bool is_pivot_col( index idx ) const {
    return idx_of_pivot_cols() == idx;
  }
      
  void release_pivot_col() {
    index idx = idx_of_pivot_cols();
    if( idx != -1 ) {
      Base::_clear(idx);
      column col;
      pivot_cols().get_col_and_clear( col );
      Base::_set_col(idx, col);
    }
    idx_of_pivot_cols() = -1;
  }
  
  void make_pivot_col( index idx ) {
    release_pivot_col();
    idx_of_pivot_cols() = idx;
    // Remark: Converting the column to a vector first gives a huge performance
    //         penalty, due to the fact that data is copied around
    //column col;
    //Base::_get_col( idx, col );
    
    get_pivot_col().add_col( Base::col_traits.col_at(Base::matrix,idx).begin(),
			     Base::col_traits.col_at(Base::matrix,idx).end());
    
    //get_pivot_col().add_col( Base::matrix[idx]);
  }
  
 public:  

  void _set_dimensions( index nr_of_rows, index nr_of_cols ) {
    #pragma omp parallel for
    for( int tid = 0; tid < omp_get_num_threads(); tid++ ) {
      pivot_cols[ tid ].init( nr_of_rows );
      idx_of_pivot_cols[ tid ] = -1;
    }
    Base::_set_dimensions( nr_of_rows, nr_of_cols );
  }

  void _add_to( index source, index target ) {
    if( !is_pivot_col( target ) )
      make_pivot_col( target );
    //column col;
    //Base::_get_col(source,col);
    
    get_pivot_col().add_col( Base::col_traits.col_at(Base::matrix,source).begin(),
			     Base::col_traits.col_at(Base::matrix,source).end() );
    
    //get_pivot_col().add_col(Base::matrix[source]);
  }

  void _sync() { 
    #pragma omp parallel for
    for( int tid = 0; tid < omp_get_num_threads(); tid++ )
      release_pivot_col();
  } 

  void _get_col( index idx, column& col  ) const { is_pivot_col( idx ) ? get_pivot_col().get_col( col ) : Base::_get_col( idx, col ); }
        
  bool _is_empty( index idx ) const { return is_pivot_col( idx ) ? get_pivot_col().is_empty() : Base::_is_empty( idx ); }

  index _get_max_index( index idx ) const { return is_pivot_col( idx ) ? get_pivot_col().get_max_index() : Base::_get_max_index( idx ); }

  void _clear( index idx ) { is_pivot_col( idx ) ? get_pivot_col().clear() : Base::_clear( idx ); }

  void _set_col( index idx, const column& col  ) { is_pivot_col( idx ) ? get_pivot_col().set_col( col ) : Base::_set_col( idx, col ); }

  void _remove_max( index idx ) { is_pivot_col( idx ) ? get_pivot_col().remove_max() : Base::_remove_max( idx ); }
        
  void finalize( index idx ) { Base::_finalize( idx ); }

};


}
