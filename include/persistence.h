#ifndef PERSISTENCE_H
#define PERSISTENCE_H

#include "filtration.h"
// PHAT
#include "../include/phat/compute_persistence_pairs.h" // wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include "../include/phat/representations/default_representations.h" // main data structure (choice affects performance)
#include "../include/phat/algorithms/standard_reduction.h" // algorithm (choice affects performance)
#include "../include/phat/algorithms/chunk_reduction.h"
#include "../include/phat/algorithms/row_reduction.h"
#include "../include/phat/algorithms/twist_reduction.h"
#include <iostream>
#include <limits>

/**
 * @brief The Persistence class
 * I only made a class because I have several parts that I want to split into
 * functions and that they share heavy containers
 */
struct TBball
{
    double r; // radius
    Delaunay::Point p; // ball center

    TBball(double rad = 0.0, Delaunay::Point pt = Delaunay::Point(0.0,0.0,0.0))
        :r(rad), p(pt) {}
};
struct HoleMeas
{
    TBball T;
    TBball B;
    int dimension;
    HoleMeas(TBball _T, TBball _B, int _dimension = 0)
        :T(_T), B(_B), dimension(_dimension) {}
    friend std::ostream& operator<< (std::ostream &out, const HoleMeas &pair);
};


std::ostream& operator<< (std::ostream &out, const HoleMeas &pair)
{
    out  << pair.dimension << " "
         <<  pair.T.r << " " << pair.T.p << " "
         <<  pair.B.r << " " << pair.B.p;
    return out;
}

template <class Simplex>
class Persistence
{
public:
    Persistence(Filtration<Simplex>& f);

    void run_persistence();
    void compute_holes_from_pairs(bool add_first_0_hole = true);

    phat::persistence_pairs get_tb_pairs() {return tb_pairs;}
    std::vector<HoleMeas> get_holes() {return holes;}

private:
    phat::boundary_matrix< phat::bit_tree_pivot_column > boundary_matrix;
    phat::persistence_pairs tb_pairs; // pairs (t_index, b_index)
    std::vector<HoleMeas> holes; // the hole signed-balls
    Filtration<Simplex>& ft;
};

/**
 * @brief alexander_deduction
 * Perform an alexander deduction on the holes:
 * (T,B,dim) -> (B,T,3-dim-1)
 */
inline void alexander_deduction(std::vector<HoleMeas>& holes) {
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
inline void tb_pairing(
    const std::vector<HoleMeas>& t_holes,
    const std::vector<HoleMeas>& b_holes,
    std::vector<HoleMeas>& output) {
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
 * criteria(T.r, B.r). Return the number of such holes.
 * @Precondition : holes should have been computed.
 */
inline int save_holes_criteria(std::vector<HoleMeas> holes, std::string filename,
    bool (*criteria)(double,double),
    std::string extension = ".tb") {
    // extract the TB balls
    std::ofstream file(filename + extension, std::ios::out | std::ios::trunc);
    if (!(file))
    {
        std::cerr << "Error in Persistence<Simplex>::save_holes_criteria(): "
        << "impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }
    int n = 0;
    for(HoleMeas hole : holes)
    {
        if (criteria(hole.T.r, hole.B.r))
        {
            file << hole << std::endl;
            n += 1;
        }
    }
    file.close();
    std::clog << "hole measures saved in "<< filename + extension
              << ", recall that each line contains: " << std::endl
              << "dimension t_radius t_center[3] b_radius b_center[3]" << std::endl;
    return n;
}

/**
 * @brief Persistence<Simplex>::save_holes
 * Save the every hole measures in a file filename+extension.
 * Return the number of such holes.
 * @Precondition : holes should have been computed.
 */
inline int save_holes(std::vector<HoleMeas> holes, std::string filename,
    std::string extension = ".tb") {
    std::clog << "exhaustive ";
    return save_holes_criteria(holes, filename,
        [](double t, double b) -> bool {return t > -b;},
        extension);
}

/**
 * @brief Persistence<Simplex>::save_tb_balls
 * Save the tb ball pairs that have negative birth date and positive death rate
 * in a file filename+extension. Return the number of such holes.
 * @Precondition : holes should have been computed.
 */
inline int save_present_holes(std::vector<HoleMeas> holes, std::string filename,
    std::string extension = ".tb") {
    std::clog << "present ";
    return save_holes_criteria(holes, filename,
        [](double t, double b) -> bool {return t > 0 && b > 0;},
        extension);
}

/*############################################################################*/
/*############################# Persistence class ############################*/
/*############################################################################*/

/* @note
 * The Persistence class is defined here because we are dealing with templates.
 * see https://cpp.developpez.com/faq/cpp/?page=Les-templates#Pourquoi-mes-templates-ne-sont-ils-pas-reconnus-a-l-edition-des-liens
 * for more details about this problem
 */

/**
 * @brief Persistence constructor
 * Build the boundary matrix from a Filtration f.
 */
template <class Simplex>
Persistence<Simplex>::Persistence(Filtration<Simplex>& f)
: ft(f)
{
    const std::map<Simplex, std::vector<int>>& m_coboundary
        = f.get_coboundary();
    // make boundary matrix
    boundary_matrix.set_num_cols(f.get_filter_size());
    std::vector< phat::index > temp_col;
    for (std::size_t i = 0; i < f.get_filter_size(); i++)
    {
        const Simplex &s = f.get_filter(i);
        boundary_matrix.set_dim(i, 3-s.dimension());
        temp_col.resize(m_coboundary.at(s).size());
        int j = 0;
        for (const int coface_pos : m_coboundary.at(s))
        {
            temp_col.at(j) = coface_pos;
            j++;
        }
        std::sort(temp_col.begin(), temp_col.end()); // @note This is very important! // they should already be sorted
        boundary_matrix.set_col( i, temp_col );
    }
    std::clog << "-- the boundary matrix has " << boundary_matrix.get_num_cols()
    << " columns and " << boundary_matrix.get_num_entries() << " entries."
    << std::endl;
    //boundary_matrix.save_ascii("phat_boundary.dat");
}

/**
 * @brief Persistence<Simplex>::run_persistence
 * Run the phat algorithm to compute tb_pairs.
 */
template <class Simplex>
void Persistence<Simplex>::run_persistence()
{
    // compute persistence pairs
    tb_pairs.clear();
    phat::compute_persistence_pairs<phat::chunk_reduction>(tb_pairs, boundary_matrix);
    std::clog << "-- persistence pairs computed" << std::endl;
}


/**
 * @brief Persistence<Simplex>::compute_holes_from_pairs
 * Compute the holes from the filtration and the tb_pairs.
 * @Precondition : tb_pairs should have been computed (eg by run_persistence).
 */
template <class Simplex>
void Persistence<Simplex>::compute_holes_from_pairs(bool add_first_0_hole)
{
    std::set<int> unseen_index;
    for (int i = (add_first_0_hole) ? 0 : 1 ; i < ft.get_filter_size(); i++) {
        unseen_index.insert(i);
    }

    holes.clear();
    for(phat::index idx = 0; idx < tb_pairs.get_num_pairs(); idx++)
    {
        const int t_index = tb_pairs.get_pair(idx).first;
        const int b_index = tb_pairs.get_pair(idx).second;
        unseen_index.erase(b_index);
        unseen_index.erase(t_index);
        const int dimension = 3-ft.get_filter(t_index).dimension();
        const TBball T(-ft.get_filtration(t_index), ft.get_point(t_index));
        const TBball B(ft.get_filtration(b_index), ft.get_point(b_index));
        holes.push_back(HoleMeas(T,B,dimension));
    }

    const double inf = std::numeric_limits<double>::infinity();
    /* Due to the persistence algorithm (hidden in compute_persistence_pairs),
     * the unseen_index represent the holes that never died and are still
     * present at the end of the filtration.
     * The first 0-hole is also in this case. */
    for (int n : unseen_index)
    {
        const int dimension = 3-ft.get_filter(n).dimension();
        const TBball T(-ft.get_filtration(n), ft.get_point(n));

        /* actually there is no point associated with the B-measure: putting
         * p_infinity is theoretically wrong, but used for implementation
         * simplicity.*/
        const TBball B(inf, Point(inf,inf,inf));
        holes.push_back(HoleMeas(T,B,dimension));
    }
    std::clog << "-- holes computed" << std::endl;
}


//#include "../src/persistence.tpp"

#endif // PERSISTENCE_H
