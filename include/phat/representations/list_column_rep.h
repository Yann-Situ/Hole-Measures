/*  Copyright 2013 IST Austria
    Contributed by: Jan Reininghaus

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
    class list_column_rep {

    protected:
      std::list< index > indices;

    public:

        typedef list_column_rep Self;
      
	std::list<index>::const_iterator begin() const {
	  return indices.begin();
	}

	std::list<index>::const_iterator end() const {
	  return indices.end();
	}

	void offer_thread_local_storage(thread_local_storage< column >* tls) {
	}

        // replaces(!) content of 'col' with boundary of given index
        void _get_col( column& col  ) const {
            col.clear();
            col.reserve( indices.size() );
            std::copy (indices.begin(), indices.end(), std::back_inserter(col) );
        }

        void _set_col( const column& col  ) {
            indices.clear();
            indices.resize( col.size() );
            std::copy (col.begin(), col.end(), indices.begin() );
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
            return indices.empty() ? -1 : *indices.rbegin();
        }

        // removes the maximal index of a column
        void _remove_max() {
            std::list< index >::iterator it = indices.end();
            it--;
            indices.erase( it );
        }
        
        // clears given column
        void _clear() {
            indices.clear();    
        }

        // syncronizes all data structures (essential for openmp stuff)
        void _sync() {}

        // adds column 'source' to column 'target'
        void _add_to( Self& source) {
            std::list< index >& source_col = source.indices;
            std::list< index >& target_col = this->indices;
            std::list< index > temp_col;
            target_col.swap( temp_col );
            std::set_symmetric_difference( temp_col.begin(), temp_col.end(),
                                           source_col.begin(), source_col.end(),
                                           std::back_inserter( target_col ) );
        }
        
        // finalizes given column
        void _finalize() {
        }
    };
}
