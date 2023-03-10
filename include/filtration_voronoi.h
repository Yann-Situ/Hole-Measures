#ifndef FILTRATION_VORONOI_H
#define FILTRATION_VORONOI_H


#include "cgal_typedef.h"
#include "filtration.h"
/**
 * @brief The FiltrationVoronoi class
 * Contains functions to create and handle the voronoi/delaunay filtration
 * of a mesh file .off.
 */

enum class CellCritType {NonCritical, Critical};

class FiltrationVoronoi : public Filtration<Delaunay::Simplex>
{
public:
    typedef Delaunay::Simplex Simplex;
    FiltrationVoronoi(Polyhedron poly);//const std::string &filename);

    void init_Dcell_filtration_values();
    void compute_infinite_voronoi_face_duals();
    void make_filter();
    void compute_delaunay_coboundary();

    bool is_finite_medial(const Simplex s);
    bool is_on_boundary(const Delaunay::Facet f);
    std::deque<Delaunay::Facet> boundary_facets();

    std::pair<Delaunay::Point, double> search_critical(int pos, bool use_max);

    Polyhedron m_poly;
    Delaunay m_dela;

    // Filtration function
    int get_filter_size() {return m_filter.size();}
    Simplex get_filter(int pos) {return m_filter.at(pos);}
    Delaunay::Point get_point(int pos) {return m_dela.dual(m_cell.at(pos));}
    double get_filtration(int pos) {return m_filtration.at(pos);}
    const std::map<Simplex, std::vector<int>>& get_coboundary() {return m_coboundary ;}

private:
    std::set<Simplex> infinite_voronoi_face_duals;   /// set of Delaunay face whose dual is infinite
    std::set<Simplex> boundary_delaunay_faces;   /// set of Delaunay faces whose dual links inside to outside

    std::map<Delaunay::Cell_handle, double> m_cell_filtration;  /// map a Dcell to its filtration value
    std::deque<double> m_filtration;                           /// map a position (in filter) to its filtration value
    std::deque<Simplex> m_filter;                    /// list of simplices in the filter
    //std::deque<Delaunay::Point> m_point;                       /// map a position of a simplex in filter -> cell -> circumcenter
    std::deque<Delaunay::Cell_handle> m_cell;                       /// map a position of a simplex in filter to its associated Dcell

    std::map<Simplex, std::vector<int>> m_coboundary; /// map a simplex to the list of its Delaunay coboudary cell filter positions
    // SoT inside;
};


#endif // FILTRATION_VORONOI_H
