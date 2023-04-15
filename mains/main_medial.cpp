#include "input_parser.h"
#include "filtration_medial.h"
#include "persistence.h"

#include <fstream>
#include <limits>
#include <time.h>

/* Input parser code found here:
 * https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c */

int main(int argc, char* argv[])
{
    InputParser input_parser(argc, argv);

    if(input_parser.cmdOptionExists("--help") || input_parser.cmdOptionExists("-h") || argc <= 1){
        std::clog << "Usage: main_medial object.off [-I] [-O] [-E] [-c] [-o output_file] [-h]" << std::endl
        << "Compute hole measures of the 3D object object.off using the medial "
        << "axis filtration.\nBy default, store only the present hole measures "
        << "(TB-balls for early-birth-late-death holes) in the object.medial.tb "
        << "file using the following format: " << std::endl
        << "dimension t_radius t_center.x t_center.y t_center.z b_radius "
        << "b_center.x b_center.y b_center.z" << std::endl
        << "-I, --in         : compute only holes measures of the inner medial axis filtration." << std::endl
        << "-O, --out        : compute only holes measures of the outer medial axis filtration." << std::endl
        << "-E, --exhaustive : store every hole measures (not only the present hole measures)." << std::endl
        << "-c, --critical   : try to add topologically critical points to the filtration (WIP)." << std::endl
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
        output_filename = output_filename + ".medial.tb";
    }
    int medial_filtration_type = 0; // 1 for outer, -1 for inner, 0 for combined
    if(input_parser.cmdOptionExists("--in") || input_parser.cmdOptionExists("-I")){
        medial_filtration_type = -1;
    }
    else if(input_parser.cmdOptionExists("--out") || input_parser.cmdOptionExists("-O")){
        medial_filtration_type = 1;
    }
    bool save_exhaustive_holes = false;
    if(input_parser.cmdOptionExists("--exhaustive") || input_parser.cmdOptionExists("-E")){
        save_exhaustive_holes = true;
    }
    bool add_critical = false;
    if(input_parser.cmdOptionExists("--critical") || input_parser.cmdOptionExists("-c")){
        add_critical = true;
    }

    /* Now perform the method */

    // evaluation variables
    clock_t time_last = clock();
    double time_delaunay = 0.0, time_filter = 0.0, time_persistence = 0.0;
    int filter_size = 0; int number_of_holes = 0;

    std::vector<HoleMeas> holes;
    FiltrationMedial F(poly);
    time_delaunay += (double)(clock() - time_last)/CLOCKS_PER_SEC;

    time_last = clock();
    F.init_finite_filtration_info();
    F.init_infinite_filtration_info();

    if (medial_filtration_type == 0)
    {
        std::clog << "Using Inner and Outer Medial Axis Filtration combined" << std::endl;

        std::clog << "- Outer Filtration:" << std::endl;
        if (add_critical){F.add_critical_to_filter(MedialType::Outer);}
        F.make_filter(MedialType::Outer);
        time_filter += (double)(clock() - time_last)/CLOCKS_PER_SEC;
        filter_size += (int)F.get_filter_size();

        time_last = clock();
        F.compute_delaunay_coboundary();
        Persistence<FiltrationMedial::Simplex> pers_out(F);
        pers_out.run_persistence();
        pers_out.compute_holes_from_pairs(false);
        // don't add the first 0-hole because we don't want to pair it
        time_persistence += (double)(clock() - time_last)/CLOCKS_PER_SEC;

        /* we perform the inner filtration after the outer one in order to add
         * the first 0-hole easily at the end of the computations:
         */

        time_last = clock();
        std::clog << "- Inner Filtration:" << std::endl;
        if (add_critical){F.add_critical_to_filter(MedialType::Inner);}
        F.make_filter(MedialType::Inner);
        time_filter += (double)(clock() - time_last)/CLOCKS_PER_SEC;
        filter_size += (int)F.get_filter_size();

        time_last = clock();
        F.compute_delaunay_coboundary();
        Persistence<FiltrationMedial::Simplex> pers_in(F);
        pers_in.run_persistence();
        pers_in.compute_holes_from_pairs(false);
        // don't add the first 0-hole because we don't want to pair it
        time_persistence += (double)(clock() - time_last)/CLOCKS_PER_SEC;

        std::clog << "- Hole Deduction:" << std::endl;
        std::vector<HoleMeas> holes_in = pers_in.get_holes();
        std::vector<HoleMeas> holes_out = pers_out.get_holes();
        alexander_deduction(holes_out);
        tb_pairing(holes_in, holes_out, holes);

        // Add the first 0-hole t-ball:
        const double inf = std::numeric_limits<double>::infinity();
        const TBball T(-F.get_filtration(0), F.get_point(0));
        const TBball B(inf, Point(inf,inf,inf));
        holes.push_back(HoleMeas(T,B,0));
    }
    else
    {
        MedialType medial_type = MedialType::Boundary;
        if (medial_filtration_type == -1)
        {
            std::clog << "Using Inner Medial Axis Filtration" << std::endl;
            medial_type = MedialType::Inner;
        }
        else //(medial_filtration_type == 1)
        {
            std::clog << "Using Outer Medial Axis Filtration" << std::endl;
            medial_type = MedialType::Outer;
        }

        if (add_critical){F.add_critical_to_filter(medial_type);std::clog << "Add critical elements to filter" << std::endl;}
        F.make_filter(medial_type);
        time_filter += (double)(clock() - time_last)/CLOCKS_PER_SEC;
        filter_size += (int)F.get_filter_size();

        time_last = clock();
        F.compute_delaunay_coboundary();
        Persistence<FiltrationMedial::Simplex> pers(F);
        pers.run_persistence();
        pers.compute_holes_from_pairs();
        time_persistence += (double)(clock() - time_last)/CLOCKS_PER_SEC;

        holes = pers.get_holes();
    }
    // if (save_exhaustive_holes){
    //     number_of_holes = save_holes(holes, output_filename, "");
    // }
    // else {
    //     number_of_holes = save_present_holes(holes, output_filename, "");
    // }

    // std::cout << "## Medial : " << std::string(filename) << std::endl;
    std::cout.precision(3);
    // std::cout   //<< "del: " << std::fixed << time_delaunay    << "    "
    //             << "fil: " << std::fixed << time_delaunay+time_filter  << "    "
    //             << "per: " << std::fixed << time_persistence << "    "
    //             << "bis: " << std::fixed << time_filter_bis  << "    " << std::endl;
    // std::cout   << "tot: " << std::fixed << time_delaunay+time_filter+time_persistence << std::endl;
    //
    // std::cout   << "sampling size: " << poly.size_of_vertices() << std::endl;
    // std::cout   << "filter size  : " << filter_size << std::endl;
    // //std::cout   << "holes size   : " << number_of_holes << std::endl;
    // std::cout   << std::endl;

    std::cout   << std::string(filename) << " "
                << poly.size_of_vertices() << " "
                << filter_size << " "
                << std::fixed << time_delaunay+time_filter << " "
                << std::fixed << time_persistence << " ";

    return 0;
}
