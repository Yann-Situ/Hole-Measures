#include "input_parser.h"
#include "filtration_medial.h"

#include <fstream>
#include <limits>

/* Input parser code found here:
 * https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c */

int main(int argc, char* argv[])
{
    InputParser input_parser(argc, argv);

    if(input_parser.cmdOptionExists("--help") || input_parser.cmdOptionExists("-h") || argc <= 1){
        std::clog << "Usage: main_medial_count object.off [-h]" << std::endl
        << "Count the size of the Delaunay/Voronoi complexe from the 3D object object.off."
        << "Display the size of the inner medial axis, the outer medial axis and the boundary (or link) in terms of cells."
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

    FiltrationMedial F(poly);
    F.init_medial_info();
    std::array<int, 3> count = F.get_medial_count();
    // int c_boundary = F.get_medial_count(MedialType::Boundary);
    // int c_inner = F.get_medial_count(MedialType::Inner);
    // int c_outer = F.get_medial_count(MedialType::Outer);
    int c_boundary = count[0];
    int c_inner = count[1];
    int c_outer = count[2];
    int c_total = c_boundary+c_inner+c_outer;
    double boundary_prop = 1.0*c_boundary/c_total;
    double inner_prop = 1.0*c_inner/c_total;
    double outer_prop = 1.0*c_outer/c_total;
    std::cout << "cells: Bound.\t" << c_boundary <<
        "\tInner.\t" << c_inner << "\tOuter.\t" << c_outer << std::endl;
    std::cout << "Total.\t" << c_total << std::endl;
    std::cout << "props: Bound.\t" << boundary_prop <<
        "\tInner.\t" << inner_prop << "\tOuter.\t" << outer_prop << std::endl;
    std::cout << "Medial axes approach efficiency:\t" <<
        inner_prop*inner_prop + outer_prop*outer_prop << std::endl;

    return 0;
}
