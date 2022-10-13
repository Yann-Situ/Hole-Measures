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

void alexander_deduction(std::vector<HoleMeas>& holes);
void tb_pairing(
    const std::vector<HoleMeas>& t_holes,
    const std::vector<HoleMeas>& b_holes,
    std::vector<HoleMeas>& output);
void save_holes_criteria(std::vector<HoleMeas> holes, std::string filename,
    bool (*criteria)(double,double),
    std::string extension = ".tb");
void save_holes(std::vector<HoleMeas> holes, std::string filename,
    std::string extension = ".tb");
void save_present_holes(std::vector<HoleMeas> holes, std::string filename,
    std::string extension = ".tb");
//void save_tb_pairs_bis(std::string filename ,std::string extension = "_V_bis.tb");

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
