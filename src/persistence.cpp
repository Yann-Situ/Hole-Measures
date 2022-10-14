/*
Persistence
computes persistence of the Voronoi filtration of a mesh.
*/

#include "persistence.h"

#include <iostream>
#include <limits>

/* @note
 * The Persistence class is defined in persistence.h because we are dealing
 * with templates.
 * see https://cpp.developpez.com/faq/cpp/?page=Les-templates#Pourquoi-mes-templates-ne-sont-ils-pas-reconnus-a-l-edition-des-liens
 * for more details about this problem
 */

std::ostream& operator<< (std::ostream &out, const HoleMeas &pair)
{
    out  << pair.dimension << " "
         <<  pair.T.r << " " << pair.T.p << " "
         <<  pair.B.r << " " << pair.B.p;
    return out;
}

/*#################### out of class functions #########################*/
/**
 * @brief alexander_deduction
 * Perform an alexander deduction on the holes:
 * (T,B,dim) -> (B,T,3-dim-1)
 */
void alexander_deduction(std::vector<HoleMeas>& holes)
{
    for (HoleMeas& hm : holes)
    {
        const TBball temp = hm.T;
        hm.T = hm.B; hm.B = temp;
        hm.dimension = 3-hm.dimension-1;
    }
    std::clog << "-- alexander deduction done" << std::endl;
}

/**
 * @brief tb_pairing
 * Arbitrary pairing for lone T-holes and lone B-holes (i.e (T,infinity,dim)
 * and (infinity,B,dim) holes).
 * output should contains all the already paired holes and the newly paired
 * holes.
 * @note t_holes should have as many lone T-holes as the lone B-holes of
 * b_holes of same dimension. If it is not the case, the output will match as
 * many possible lone balls but will leave the remaining lone balls.
 */
void tb_pairing(
    const std::vector<HoleMeas>& t_holes,
    const std::vector<HoleMeas>& b_holes,
    std::vector<HoleMeas>& output)
{
    output.clear();
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<std::vector<HoleMeas>> lone_t, lone_b;
    for (int dim = 0; dim < 3; dim ++)
    {
        std::vector<HoleMeas> t_temp;
        lone_t.push_back(t_temp);
        std::vector<HoleMeas> b_temp;
        lone_b.push_back(b_temp);
    }

    for (HoleMeas hm : t_holes)
    {
        if (hm.B.r >= inf)// lone_T (only a T value)
            lone_t[hm.dimension].push_back(hm);
        else
            output.push_back(hm);
    }
    for (HoleMeas hm : b_holes)
    {
        if (hm.T.r >= inf)// lone_B (only a B value)
            lone_b[hm.dimension].push_back(hm);
        else
            output.push_back(hm);
    }

    for (int dim = 0; dim < 3; dim ++)
    {
        if (lone_t[dim].size() != lone_b[dim].size())
        {
            std::cerr << "Error in Persistence<Simplex>::tb_pairing: "
            << lone_t[dim].size() << " lone_t but "
            << lone_b[dim].size() << " lone_b of dimension" << dim << std::endl;
        }
        int count = (lone_t[dim].size() < lone_b[dim].size()) ?
                     lone_t[dim].size() : lone_b[dim].size(); // min
        for (int i = 0; i < count; i++)
        {
            output.push_back(HoleMeas(lone_t[dim][i].T, lone_b[dim][i].B, dim));
        }
    }
    std::clog << "-- lone balls paired arbitrarily" << std::endl;
}

/**
 * @brief Persistence<Simplex>::save_holes_criteria
 * Save the hole measures in a file filename+extension that verifies
 * criteria(T.r, B.r).
 * @Precondition : holes should have been computed.
 */
void save_holes_criteria(std::vector<HoleMeas> holes, std::string filename,
    bool (*criteria)(double,double), std::string extension)
{
    // extract the TB balls
    std::ofstream file(filename + extension, std::ios::out | std::ios::trunc);
    if (!(file))
    {
        std::cerr << "Error in Persistence<Simplex>::save_holes_criteria(): "
        << "impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }
    for(HoleMeas hole : holes)
    {
        if (criteria(hole.T.r, hole.B.r))
            file << hole << std::endl;
    }
    file.close();
    std::cout << "hole measures saved in "<< filename + extension
    << ", recall that each line contains: " << std::endl
    << "dimension t_radius t_center[3] b_radius b_center[3]" << std::endl;

}

//bool valid_hole(double t, double b){return t > -b;};
//bool present_hole(double t, double b){return t > 0 && b > 0;};

/**
 * @brief Persistence<Simplex>::save_holes
 * Save the every hole measures in a file filename+extension.
 * @Precondition : holes should have been computed.
 */
void save_holes(std::vector<HoleMeas> holes, std::string filename,
    std::string extension)
{
    std::cout << "exhaustive ";
    save_holes_criteria(holes, filename,
        [](double t, double b) -> bool {return t > -b;},
        extension);
}

/**
 * @brief Persistence<Simplex>::save_tb_balls
 * Save the tb ball pairs that have negative birth date and positive death rate
 * in a file filename+extension.
 * @Precondition : holes should have been computed.
 */
void save_present_holes(std::vector<HoleMeas> holes, std::string filename,
    std::string extension)
{
    std::cout << "present ";
    save_holes_criteria(holes, filename,
        [](double t, double b) -> bool {return t > 0 && b > 0;},
        extension);
}
