////////////////////////////////////////////////////////////////////////////////
/* This is a work in progress. We aim to find an algorithm that produce an
 * appropriate sampling on the mesh (ideally an epsilon sampling) that provides
 * guarantees on the topology of the medial axis approximation. */
////////////////////////////////////////////////////////////////////////////////
#include "filtration_voronoi.h"
#include "resampling.h"
//#include "persistence.h"
#include <CGAL/draw_triangulation_3.h>
//#include <CGAL/draw_polyhedron.h>

// for annoying conversion between polyhedron_3 and Mesh_polyhedron_3 :
#include <CGAL/Polyhedron_copy_3.h>
typedef CGAL::Polyhedron_copy_3<Polyhedron, M_Polyhedron::HDS> Poly_copy;

#include <vector>
#include <fstream>
#include <limits>

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;// edge_size = 0.025 instead of CGAL::parameters::edge_size

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "../data/eight.off";
  std::string m_filename(filename);
  m_filename.erase(m_filename.end()-4, m_filename.end()); // remove ".off"

  // open the mesh file and import into a Polyhedron
  Polyhedron poly;
  std::ifstream input(filename);
  if (!input || !(input >> poly) || poly.empty() || !CGAL::is_triangle_mesh(poly))
  {
    std::cerr << std::string(filename) + " is not a valid input file." << std::endl;
    exit(EXIT_FAILURE);
  }

  const double facet_size_refinement = (argc > 2) ? std::stod(std::string(argv[2])) : 1.0;
  const double edge_size_refinement = (argc > 3) ? std::stod(std::string(argv[3])) : 0.025;

  FiltrationVoronoi F(poly);

  // //Resampling::add_circumcenters(S.m_dela, S.m_poly);
  // SoT inside(F.m_poly);
  // Resampling::add_from_criteria_cell(F.m_dela, inside);
  // CGAL::draw(F.m_dela);

  /*##########################################################################*/
  // https://doc.cgal.org/latest/Mesh_3/index.html


  M_Polyhedron m_poly;
  Poly_copy polyhedron_copy_modifier(poly);
  // The following line should copy 'poly' to 'm_poly'
  m_poly.delegate(polyhedron_copy_modifier);

  // Create a vector with only one element: the pointer to the polyhedron.
  std::vector<M_Polyhedron*> poly_ptrs_vector(1, &m_poly);
  // Create a polyhedral domain, with only one polyhedron,
  // and no "bounding polyhedron", so the volumetric part of the domain will be
  // empty.
  Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
  // Get sharp features
  domain.detect_features(); //includes detection of borders
  // Mesh criteria
  MedialFeatureField field(F.m_dela, facet_size_refinement);
  Mesh_criteria criteria(//edge_size = edge_size_refinement,
      //facet_angle = 25,
      facet_size = field//facet_size_refinement//,
      //facet_distance = 0.001
      );

  /* facet_size: This parameter controls the size of surface facets. Each
   * surface facet has a surface Delaunay ball which is a ball circumscribing
   * the surface facet and centered on the surface patch. The parameter
   * facet_size is either a constant or a spatially variable scalar field,
   * providing an upper bound for the radii of surface Delaunay balls.
   */
   // we can use a sclar field for the facet size to have a good

  // Mesh generation
  std::clog << "make_mesh..." << std::endl;
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
  // Output the facets of the c3t3 to an OFF file. The facets will not be
  // oriented.
  std::clog << "nb finite vertices: " << c3t3.triangulation().number_of_vertices() << std::endl;
  CGAL::draw(c3t3.triangulation());// noice
  std::ofstream off_file("out.off");
  c3t3.output_boundary_to_off(off_file);
  return off_file.fail() ? EXIT_FAILURE : EXIT_SUCCESS;
  /*##########################################################################*/

  return 0;
}
