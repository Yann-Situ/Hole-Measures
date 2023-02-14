#include "input_parser.h"
#include "filtration_medial.h"
#include "persistence.h"

#include <fstream>
#include <limits>

/* Input parser code found here:
 * https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c */

int main(int argc, char* argv[])
{
    InputParser input_parser(argc, argv);

    if(input_parser.cmdOptionExists("--help") || input_parser.cmdOptionExists("-h") || argc <= 1){
        std::clog << "Usage: main_medial_axes object.off [-I] [-O] [-E] [-c] [-o output_file] [-h]" << std::endl
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

    std::vector<HoleMeas> holes;
    FiltrationMedial F(poly);
    F.init_medial_info();
    F.init_simplex_faces();

    if (medial_filtration_type == 0)
    {
        std::clog << "Using Inner and Outer Medial Axis Filtration combined" << std::endl;

        std::clog << "- Outer Filtration:" << std::endl;
        if (add_critical){F.add_critical_to_filter(MedialType::Outer);}
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
        if (add_critical){F.add_critical_to_filter(MedialType::Inner);}
        F.make_filter(MedialType::Inner);
        F.compute_delaunay_coboundary();
        Persistence<FiltrationMedial::Simplex> pers_in(F);
        pers_in.run_persistence();
        pers_in.compute_holes_from_pairs(false);
        // don't add the first 0-hole because we don't want to pair it

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
        F.compute_delaunay_coboundary();
        Persistence<FiltrationMedial::Simplex> pers(F);
        pers.run_persistence();
        pers.compute_holes_from_pairs();
        holes = pers.get_holes();
    }
    if (save_exhaustive_holes){
        save_holes(holes, output_filename, "");
    }
    else {
        save_present_holes(holes, output_filename, "");
    }

    return 0;
}
