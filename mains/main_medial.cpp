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

    if      (argc > 2 && std::string(argv[2]) == "-i")
    {
        std::clog << "Using Inner Medial Axis Filtration" << std::endl;
        F.make_filter(MedialType::Inner);
        F.compute_delaunay_coboundary();

        Persistence<FiltrationMedial::Simplex> pers(F);
        pers.run_persistence();
        pers.compute_holes_from_pairs();
        save_present_holes(pers.get_holes(), m_filename, ".medial.tb");
    }
    else if (argc > 2 && std::string(argv[2]) == "-o")
    {
        std::clog << "Using Outer Medial Axis Filtration" << std::endl;
        F.make_filter(MedialType::Outer);
        F.compute_delaunay_coboundary();

        Persistence<FiltrationMedial::Simplex> pers(F);
        pers.run_persistence();
        pers.compute_holes_from_pairs();
        save_present_holes(pers.get_holes(), m_filename, ".medial.tb");
    }
    else
    {
        std::clog << "Using Inner and Outer Medial Axis Filtration combined" << std::endl;

        std::clog << "- Outer Filtration:" << std::endl;
        F.make_filter(MedialType::Outer);
        F.compute_delaunay_coboundary();
        Persistence<FiltrationMedial::Simplex> pers_out(F);
        pers_out.run_persistence();
        pers_out.compute_holes_from_pairs(false);
        // don't add the first 0-hole because we don't want to pair it

        /* we perform the inner filtration after the outer one in order to add
         * the first 0-hole easily at the end of the computations:
         */
        std::clog << "- Inner Filtration:" << std::endl;
        F.make_filter(MedialType::Inner);
        F.compute_delaunay_coboundary();
        Persistence<FiltrationMedial::Simplex> pers_in(F);
        pers_in.run_persistence();
        pers_in.compute_holes_from_pairs(false);
        // don't add the first 0-hole because we don't want to pair it

        std::clog << "- Hole Deduction:" << std::endl;
        std::vector<HoleMeas> holes_in = pers_in.get_holes();
        std::vector<HoleMeas> holes_out = pers_out.get_holes();
        std::vector<HoleMeas> holes;
        alexander_deduction(holes_out);
        tb_pairing(holes_in, holes_out, holes);

        // Add the first 0-hole t-ball:
        const double inf = std::numeric_limits<double>::infinity();
        const TBball T(-F.get_filtration(0), F.get_point(0));
        const TBball B(inf, Point(inf,inf,inf));
        holes.push_back(HoleMeas(T,B,0));

        save_present_holes(holes, m_filename, ".medial.tb");
    }

    return 0;
}
