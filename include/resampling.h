////////////////////////////////////////////////////////////////////////////////
/* This is a work in progress. We aim to find an algorithm that produce an
 * appropriate sampling on the mesh (ideally an epsilon sampling) that provides
 * guarantees on the topology of the medial axis approximation. */
////////////////////////////////////////////////////////////////////////////////
#ifndef RESAMPLING_H
#define RESAMPLING_H

#include "cgal_typedef.h"

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

// Polyhedron type
typedef CGAL::Mesh_polyhedron_3<Epick>::type M_Polyhedron;
// Domain
typedef CGAL::Polyhedral_mesh_domain_with_features_3<Epick> Mesh_domain;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

/**
 * @brief The Resampling class
 *
 */
class Resampling
{
public:
    static void add_circumcenters(Delaunay &m_dela, Polyhedron m_poly);
    static Point circumcenter(Polyhedron::Facet_iterator f);

    static void add_barycenters(Delaunay &m_dela, Polyhedron m_poly);
    static Point barycenter(Polyhedron::Facet_iterator f);
    // static Point cross(const Point a, const Point b);
    // static double len2(const Point a);
    static void write_off(Polyhedron m_poly);
    static void add_from_criteria(Delaunay &m_dela,
        SoT inside);
    static void add_from_criteria_cell(Delaunay &m_dela,
        SoT inside);

private:
    Resampling();
    ~Resampling();

};

// For sizing field
    typedef Epick::FT FT;
    typedef Point Point_3;
    typedef Mesh_domain::Index Index;

// Sizing field
class MedialFeatureField
{
public :
    MedialFeatureField(const Delaunay &m_dela, double parameter = 1.0);
    FT operator()(const Point_3& p, const int dimension, const Index& index) const;
private :
    Delaunay m_dela;
    double parameter;
    //std::map<Delaunay::Vertex_handle, FT> sqr_circumradius;
};

#endif // RESAMPLING_H
