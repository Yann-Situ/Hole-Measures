#ifndef CGAL_TYPEDEF_H
#define CGAL_TYPEDEF_H

// Polyhedron
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
// Delaunay
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
typedef Epick::Point_3 Point;
typedef Epick::Vector_3 Vector;
typedef Epick::Plane_3 Plane;
typedef Epick::Segment_3 Segment;
typedef Epick::Triangle_3 Triangle;
typedef CGAL::Polyhedron_3<Epick> Polyhedron;

typedef CGAL::Triangulation_vertex_base_with_info_3<Epick::FT, Epick>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<Epick>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                    Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<Epick, Tds, CGAL::Fast_location> Delaunay;
typedef CGAL::Side_of_triangle_mesh<Polyhedron, Epick>                  SoT;


#endif // CGAL_TYPEDEF_H
