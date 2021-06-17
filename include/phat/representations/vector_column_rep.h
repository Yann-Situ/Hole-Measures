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

namespace phat {
    class vector_column_rep {

    protected:

        std::vector<index> indices;

        thread_local_storage< column >* temp_column_buffer;

    public:

	typedef vector_column_rep Self;

	void offer_thread_local_storage(thread_local_storage< column >* tls) {
	  temp_column_buffer=tls;
	}

	std::vector<index>::const_iterator begin() const {
	  return indices.begin();
	}

	std::vector<index>::const_iterator end() const {
	  return indices.end();
	}

        // replaces(!) content of 'col' with boundary of given index
        void _get_col( column& col  ) const {
	  col.clear();
	  col = indices; 
        }
        void _set_col( const column& col  ) { 
            indices = col; 
        }

	void _set_nr_of_rows( int nr_of_rows ) {
	  // ignore
	}

        // true iff boundary of given idx is empty
        bool _is_empty() const { 
            return indices.empty(); 
        }

        // largest row index of given column idx (new name for lowestOne())
        index _get_max_index() const { 
            return indices.empty() ? -1 : indices.back(); 
        }

        // removes the maximal index of a column
        void _remove_max() {
            indices.pop_back();
        }

        // clears given column
        void _clear() { 
            indices.clear();
        }

        // syncronizes all data structures (essential for openmp stuff)
        void _sync() {}

        // adds column 'source' to column 'target'
        void _add_to( const Self& source_col ) {
	  column& temp_col = (*temp_column_buffer)();
	              
	  size_t new_size = source_col.indices.size() + this->indices.size();
	  
	  if (new_size > temp_col.size()) temp_col.resize(new_size);
          
	  std::vector<index>::iterator col_end = std::set_symmetric_difference( this->indices.begin(), this->indices.end(),
										source_col.indices.begin(), source_col.indices.end(),
										temp_col.begin() );
	  temp_col.erase(col_end, temp_col.end());
	  
          
	  this->indices.swap(temp_col);
        }
        
        // finalizes given column
        void _finalize() {
            column(indices.begin(), indices.end()).swap(indices);
        }
    };
}
