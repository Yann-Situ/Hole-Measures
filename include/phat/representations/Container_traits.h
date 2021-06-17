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

namespace phat {
  
template<typename ColumnContainer>
  class Column_container_traits {
 
  public:  

  typedef ColumnContainer Column_container;
  typedef typename Column_container::value_type Column_type;

  Column_type& col_at(Column_container& data, index idx) const {
    return data[ idx ];
  }

  const Column_type& col_at(const Column_container& data, index idx) const {
    return data[ idx ];
  }

  index get_size(const Column_container& data) const {
    return data.size();
  }

  void resize(Column_container& data, index nr_of_columns) const {
    data.resize(nr_of_columns);
  }

};

template<typename DimensionContainer>
  class Dimension_container_traits {
  
  public:

  typedef DimensionContainer Dimension_container;

  index& dim_at(Dimension_container& data, index idx) const {
    return data[ idx ];
  }

  const index& dim_at(const Dimension_container& data, index idx) const {
    return data[ idx ];
  }

  index get_size(const Dimension_container& data) const {
    return data.size();
  }

  void resize(Dimension_container& data, index nr_of_columns) const { 
    data.resize(nr_of_columns);
  }

};


} // namespace phat
