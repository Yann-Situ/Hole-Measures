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

// This file requires c++11

#include<phat/representations/Container_traits.h>

#include<unordered_map>
#include<mutex>

namespace phat{

template<class ColumnType>
  class Column_container_traits< std::unordered_map<index,ColumnType> > {
  
  protected:

  index _size;

  mutable std::mutex _mutex;

  public:

  typedef std::unordered_map<index,ColumnType> Column_container;
  
  typedef ColumnType Column_type;

  Column_type& col_at(Column_container& data, index idx) const {
    _mutex.lock();
    typename Column_container::iterator it = data.find(idx);
    if(it==data.end()) {
      std::pair<typename Column_container::iterator, bool> result 
	= data.insert(std::make_pair(idx,Column_type()));
      it = result.first;
    }
    _mutex.unlock();
    return it->second;
  }

  /*
  const Column_type& col_at(const Column_container& data, index idx) const {
    if(data.find(idx)==data.end()) {
      data.insert(std::make_pair(idx,Column_type()));
    }
    return data.at(idx);
  }
  */

  index get_size(const Column_container& data) const {
    return _size;
  }

  void resize(Column_container& data, index nr_of_columns) {
    _size = nr_of_columns;
  }

 };

template<>
  class Dimension_container_traits< std::unordered_map<index,index> > {
  
 protected:

  index _size;

  mutable std::mutex _mutex;

  public:

  typedef std::unordered_map<index,index> Dimension_container;

  index& dim_at(Dimension_container& data, index idx) const {
    _mutex.lock();
    typename Dimension_container::iterator it = data.find(idx);
    if(it==data.end()) {
      std::pair<typename Dimension_container::iterator, bool> result 
	= data.insert(std::make_pair(idx,0));
      it = result.first;
    }
    _mutex.unlock();
    return it->second;
  }

  /*
  const index& dim_at(Dimension_container& data, index idx) const {
    if(data.find(idx)==data.end()) {
      data.insert(std::make_pair(idx,0));
    }
    return data.at(idx);
  }
  */

  index get_size(const Dimension_container& data) const {
    return _size;
  }

  void resize(Dimension_container& data, index nr_of_columns) { 
    _size = nr_of_columns;
  }

 };

} // of namespace phat
