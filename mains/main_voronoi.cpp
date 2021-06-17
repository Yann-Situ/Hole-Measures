#include "filtration_v.h"
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

  FiltrationV F(poly);
  F.init_Dcell_filtration_values();
  F.compute_infinite_voronoi_face_duals();
  F.make_filter();
  F.compute_delaunay_coboundary();

  Persistence<FiltrationV::Simplex> pers(F);
  pers.run_persistence();
  pers.compute_holes_from_pairs();
  pers.save_tb_balls(m_filename, "_V.tb");

  return 0;
}
