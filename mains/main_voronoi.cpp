#include "input_parser.h"
#include "filtration_voronoi.h"
#include "persistence.h"

#include <fstream>
#include <limits>
#include <time.h>


int main(int argc, char* argv[])
{
    InputParser input_parser(argc, argv);

    if(input_parser.cmdOptionExists("--help") || input_parser.cmdOptionExists("-h") || argc <= 1){
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

    const char* filename = argv[1];
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

    // evaluation variables
    clock_t time_last = clock();
    double time_delaunay = 0.0, time_filter = 0.0, time_persistence = 0.0;
    int filter_size = 0; //int hole_size = 0;

    std::vector<HoleMeas> holes;
    FiltrationVoronoi F(poly);
    time_delaunay += (double)(clock() - time_last)/CLOCKS_PER_SEC;

    time_last = clock();
    F.init_Dcell_filtration_values();
    F.compute_infinite_voronoi_face_duals();
    F.make_filter();
    time_filter += (double)(clock() - time_last)/CLOCKS_PER_SEC;
    filter_size += (int)F.get_filter_size();

    time_last = clock();
    F.compute_delaunay_coboundary();

    Persistence<FiltrationVoronoi::Simplex> pers(F);
    pers.run_persistence();
    pers.compute_holes_from_pairs();
    time_persistence += (double)(clock() - time_last)/CLOCKS_PER_SEC;

    holes = pers.get_holes();
    // if (save_exhaustive_holes){
    //     hole_size = save_holes(holes, output_filename, "");
    // }
    // else {
    //     hole_size = save_present_holes(holes, output_filename, "");
    // }

    //std::cout << "## Voronoi: " << std::string(filename) << std::endl;
    std::cout.precision(3);
    // std::cout   //<< "del: " << std::fixed << time_delaunay    << "    "
    //             << "fil: " << std::fixed << time_delaunay+time_filter      << "    "
    //             << "per: " << std::fixed << time_persistence << "    " << std::endl;
    // std::cout   << "tot: " << std::fixed << time_delaunay+time_filter+time_persistence << std::endl;
    //
    // std::cout   << "sampling size: " << poly.size_of_vertices() << std::endl;
    // std::cout   << "filter size  : " << filter_size << std::endl;
    // //std::cout   << "holes size   : " << hole_size << std::endl;
    // std::cout   << std::endl;

    std::cout   << filter_size << " "
                << std::fixed << time_delaunay+time_filter << " "
                << std::fixed << time_persistence << std::endl;

    return 0;
}
