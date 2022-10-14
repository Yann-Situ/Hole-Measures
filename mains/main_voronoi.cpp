#include "input_parser.h"
#include "filtration_voronoi.h"
#include "persistence.h"

#include <fstream>
#include <limits>



int main(int argc, char* argv[])
{
    InputParser input_parser(argc, argv);

    if(input_parser.cmdOptionExists("--help") || input_parser.cmdOptionExists("-h")){
        std::clog << "Usage: main_voronoi object.off [-E] [-o output_file] [-h]" << std::endl
        << "Compute hole measures of the 3D object object.off using the voronoi "
        << " filtration.\nBy default, store only the present hole measures "
        << "(TB-balls for early-birth-late-death holes) in the object.voronoi.tb "
        << "file using the following format: " << std::endl
        << "dimension t_radius t_center.x t_center.y t_center.z b_radius "
        << "b_center.x b_center.y b_center.z" << std::endl
        << "-E, --exhaustive : store every hole measures (not only the present hole measures)." << std::endl
        << "-o output_file   : write the hole measures in output_file." << std::endl
        << "-h, --help       : display this message." << std::endl;
        exit(EXIT_SUCCESS);
    }

    const char* filename = (argc > 1) ? argv[1] : "../data/eight.off";
    // open the mesh file and import into a Polyhedron
    Polyhedron poly;
    std::ifstream input(filename);
    if (!input || !(input >> poly) || poly.empty() || !CGAL::is_triangle_mesh(poly))
    {
        std::cerr << std::string(filename) + " is not a valid input file." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string output_filename(filename);
    const std::string& output_option = input_parser.getCmdOption("-o");
    if (!output_option.empty()){
        output_filename = output_option;
    }
    else{
        output_filename.erase(output_filename.end()-4, output_filename.end()); // remove ".off" from filename
        output_filename = output_filename + ".voronoi.tb";
    }
    bool save_exhaustive_holes = false;
    if(input_parser.cmdOptionExists("--exhaustive") || input_parser.cmdOptionExists("-E")){
        save_exhaustive_holes = true;
    }

    /* Now perform the method */


    FiltrationVoronoi F(poly);
    F.init_Dcell_filtration_values();
    F.compute_infinite_voronoi_face_duals();
    F.make_filter();
    F.compute_delaunay_coboundary();

    Persistence<FiltrationVoronoi::Simplex> pers(F);
    pers.run_persistence();
    pers.compute_holes_from_pairs();
    if (save_exhaustive_holes){
        save_holes(pers.get_holes(), output_filename, "");
    }
    else {
        save_present_holes(pers.get_holes(), output_filename, "");
    }

    return 0;
}
