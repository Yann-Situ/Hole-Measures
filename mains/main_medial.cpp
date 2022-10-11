#include "filtration_medial.h"
#include "persistence.h"

#include <vector>
#include <fstream>
#include <limits>



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

  FiltrationMedial F(poly);
  F.init_medial_info();
  F.init_simplex_faces();
  F.make_filter((argc <= 2) ? MedialType::Inner : MedialType::Outer);
  //F.add_critical_to_filter((argc <= 2) ? MedialType::Inner : MedialType::Outer);
  F.compute_delaunay_coboundary();

  Persistence<FiltrationMedial::Simplex> pers(F);
  pers.run_persistence();
  pers.compute_holes_from_pairs();
  pers.save_holes(m_filename, ".medial.tb");

  return 0;
}
