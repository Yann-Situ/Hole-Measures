#ifndef FILTRATION_H
#define FILTRATION_H


#include "cgal_typedef.h"
/**
 * @brief The Filtration class
 * Abstract class for persistence
 */

template <class Simplex>
class Filtration
{
public:
    Filtration(){};

    virtual int get_filter_size() = 0;
    virtual Simplex get_filter(int pos) = 0;
    virtual Delaunay::Point get_point(int pos) = 0;
    virtual double get_filtration(int pos) = 0;
    virtual const std::map<Simplex, std::vector<int>>& get_coboundary() = 0;

    Delaunay m_dela;
};


#endif // FILTRATION_H
