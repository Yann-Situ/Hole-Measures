#ifndef DELAUNAY_HELPER_H
#define DELAUNAY_HELPER_H

#include "cgal_typedef.h"


enum class CriticalType {NonCritical, Critical};

struct CriticalInfo
{
    CriticalType c;
    double r; // distance to border
    Delaunay::Point p; // ball center associated to the simplex

    CriticalInfo(CriticalType type = CriticalType::NonCritical,
        double dist = 0.0, Delaunay::Point pt = Delaunay::Point(0.0,0.0,0.0))
        : c(type), r(dist), p(pt) {}

    void assign(CriticalType type = CriticalType::NonCritical,
        double dist = 0.0, Delaunay::Point pt = Delaunay::Point(0.0,0.0,0.0))
        { c = type; r = dist; p = pt;}
};
/**
 * @brief The DelaunayHelper class
 */
class DelaunayHelper
{
private:
    DelaunayHelper();
    ~DelaunayHelper();

public:

	static std::list<Delaunay::Simplex> D_faces(const Delaunay::Simplex s);
    static std::list<Delaunay::Simplex> D_sub_faces(const Delaunay::Cell_handle ch);
    static std::list<Delaunay::Simplex> D_sub_faces(const Delaunay::Facet f);
    static std::list<Delaunay::Simplex> D_finite_faces(const Delaunay &m_dela, const Delaunay::Cell_handle ch);
    static std::list<Delaunay::Simplex> D_finite_sub_faces(const Delaunay &m_dela, const Delaunay::Cell_handle ch);

    static double max_edge_sq(const Delaunay &m_dela, const Delaunay::Facet f);
    static bool is_infinite(
        const Delaunay &m_dela,
        const Delaunay::Simplex s);

    static CriticalInfo get_critical_info(
        const Delaunay &m_dela,
        const Delaunay::Simplex s);

    static std::list<Delaunay::Simplex> get_incident_simplices(
        const Delaunay &m_dela,
        const Delaunay::Vertex_handle v);

    static void get_finite_neighbouring_cells(
        const Delaunay &m_dela,
        std::vector<Delaunay::Cell_handle> &cells,
        const Point p);

    static void get_finite_nearest_cells(
        const Delaunay &m_dela,
        std::vector<Delaunay::Cell_handle> &cells,
        const Point p);

    static void get_finite_incident_cells(
        const Delaunay &m_dela,
        std::vector<Delaunay::Cell_handle> &cells,
        const Delaunay::Vertex_handle v);
    static void get_finite_incident_cells(
        const Delaunay &m_dela,
        std::vector<Delaunay::Cell_handle> &cells,
        const Delaunay::Facet f);
    static void get_finite_incident_cells(
        const Delaunay &m_dela,
        std::vector<Delaunay::Cell_handle> &cells,
        const Delaunay::Edge e);
};

#endif // DELAUNAY_HELPER_H
