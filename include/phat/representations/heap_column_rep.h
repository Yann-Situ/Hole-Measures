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
    class heap_column_rep {

    protected:
        std::vector< index > indices;

        index inserts_since_last_prune;

        mutable thread_local_storage< column >* temp_column_buffer;

    protected:
        void _prune()
        {
            column& col = indices;
            column& temp_col = (*temp_column_buffer)();
            temp_col.clear();
            index max_index = _pop_max_index( col );
            while( max_index != -1 ) {
                temp_col.push_back( max_index );
                max_index = _pop_max_index( col );
            }
            col = temp_col;
            std::reverse( col.begin( ), col.end( ) );
            std::make_heap( col.begin( ), col.end( ) );
            inserts_since_last_prune = 0;
        }

        index _pop_max_index(column& col) const
        {
            if( col.empty( ) )
                return -1;
            else {
                index max_element = col.front( );
                std::pop_heap( col.begin( ), col.end( ) );
                col.pop_back( );
                while( !col.empty( ) && col.front( ) == max_element ) {
                    std::pop_heap( col.begin( ), col.end( ) );
                    col.pop_back( );
                    if( col.empty( ) )
                        return -1;
                    else {
                        max_element = col.front( );
                        std::pop_heap( col.begin( ), col.end( ) );
                        col.pop_back( );
                    }
                }
                return max_element;
            }
        }

    public:

	typedef heap_column_rep Self;

	std::vector<index>::const_iterator begin() const {
	  return indices.begin();
	}

	std::vector<index>::const_iterator end() const {
	  return indices.end();
	}

	void offer_thread_local_storage(thread_local_storage< column >* tls) {
	  temp_column_buffer=tls;
	}

        // replaces(!) content of 'col' with boundary of given index
        void _get_col( column& col ) const
        {
	  col.clear();
	  (*temp_column_buffer)() = indices;
	  index max_index = _pop_max_index( (*temp_column_buffer)() );
	  while( max_index != -1 ) {
	    col.push_back( max_index );
	    max_index = _pop_max_index( (*temp_column_buffer)() );
	  }
	  std::reverse( col.begin( ), col.end( ) );
        }
        void _set_col(const column& col )
        {
            indices = col;
            std::make_heap( indices.begin( ), indices.end( ) );
        }

	void _set_nr_of_rows( int nr_of_rows ) {
	  // ignore
	}

        // true iff boundary of given idx is empty
        bool _is_empty() const
        {
            return _get_max_index() == -1;
        }

        // largest row index of given column idx (new name for lowestOne())
        index _get_max_index() const
        {
            column& col = const_cast< column& >( indices );
            index max_element = _pop_max_index( col );
            col.push_back( max_element );
            std::push_heap( col.begin( ), col.end( ) );
            return max_element;
        }

        // removes the maximal index of a column
        void _remove_max()
        {
            _pop_max_index( indices );
        }

        // clears given column
        void _clear()
        {
            indices.clear( );
        }

	index _size() {
	    _prune();
	    return indices.size();
	}

        // syncronizes all data structures (essential for openmp stuff)
        void _sync( ) {}

        // adds column 'source' to column 'target'
        void _add_to( const Self& source )
        {           
            for( index idx = 0; idx < (index)source.indices.size( ); idx++ ) {
                this->indices.push_back( source.indices[ idx ] );
                std::push_heap( this->indices.begin(), this->indices.end() );
            }
            inserts_since_last_prune += source.indices.size();

            if( 2 * inserts_since_last_prune > ( index )this->indices.size() )
                _prune();
        }
        
        // finalizes given column
        void _finalize() {
            _prune();
        }

    };
}
