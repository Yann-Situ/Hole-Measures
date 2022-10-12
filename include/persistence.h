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
    void compute_holes_from_pairs();
    void save_holes(std::string filename ,std::string extension = ".tb",
        bool add_first_0_hole = false);
    void save_present_holes(std::string filename ,std::string extension = ".tb",
        bool add_first_0_hole = false);
    //void save_tb_pairs_bis(std::string filename ,std::string extension = "_V_bis.tb");
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
void Persistence<Simplex>::compute_holes_from_pairs()
{
    std::set<int> c;
    for (int i = 0; i < ft.get_filter_size(); i++) {
        c.insert(i);
    }

    holes.clear();
    for(phat::index idx = 0; idx < tb_pairs.get_num_pairs(); idx++)
    {
        const int t_index = tb_pairs.get_pair(idx).first;
        const int b_index = tb_pairs.get_pair(idx).second;
        c.erase(b_index);
        c.erase(t_index);
        const int dimension = 3-ft.get_filter(t_index).dimension();
        const TBball T(-ft.get_filtration(t_index), ft.get_point(t_index));
        const TBball B(ft.get_filtration(b_index), ft.get_point(b_index));
        holes.push_back(HoleMeas(T,B,dimension));
    }
    std::clog << "set of non closing cells : ";
    std::cout << "c = { ";
    for (int n : c)
        std::cout << n << ' ';
    std::cout << "}\n";
    std::cout << "v = { ";
    for (int n : c)
        std::cout << ft.get_filtration(n) << ' ';
    std::cout << "}\n";
    std::clog << "-- holes computed" << std::endl;
}

/**
 * @brief Persistence<Simplex>::save_holes
 * Save the every hole measures in a file filename+extension.
 * @Precondition : holes should have been computed.
 */
template <class Simplex>
void Persistence<Simplex>::save_holes(std::string filename,
    std::string extension, bool add_first_0_hole)
{
    // extract the TB balls
    std::cout << "Writting the TB file, recall that each line contains: " << std::endl
              << "dimension t_radius t_center[3] b_radius b_center[3]" << std::endl;
    std::ofstream file(filename + extension, std::ios::out | std::ios::trunc);
    if (!(file))
    {
        std::cerr << "Error in Persistence<Simplex>::save_holes(): impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }
    for(HoleMeas hole : holes)
    {
        if (hole.T.r > -hole.B.r)
            file << hole << std::endl;
    }
    if(add_first_0_hole) // add the unpaired cell
    {
        file << "0 "
        << - ft.get_filtration(0) << " " << ft.get_point(0) << " "
        << "inf" << std::endl;
    }

    file.close();
}

/**
 * @brief Persistence<Simplex>::save_tb_balls
 * Save the tb ball pairs that have negative birth date and positive death rate
 * in a file filename+extension.
 * @Precondition : holes should have been computed.
 */
template <class Simplex>
void Persistence<Simplex>::save_present_holes(std::string filename,
    std::string extension, bool add_first_0_hole)
{
    // extract the TB balls
    std::cout << "Writting the TB file, recall that each line contains: " << std::endl
              << "dimension t_radius t_center[3] b_radius b_center[3]" << std::endl;
    std::ofstream file(filename + extension, std::ios::out | std::ios::trunc);
    if (!(file))
    {
        std::cerr << "Error in Persistence<Simplex>::save_present_holes(): impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }
    for(HoleMeas hole : holes)
    {
        if (hole.T.r > 0 && hole.B.r > 0)
        {
            file << hole << std::endl;
        }
    }
    if(add_first_0_hole) // add the unpaired cell
    {
        file << "0 "
        << - ft.get_filtration(0) << " " << ft.get_point(0) << " "
        << "inf" << std::endl;
    }

    file.close();
}

/*
template <class Simplex>
void Persistence<Simplex>::save_tb_pairs_bis(std::string filename ,std::string extension)
{
    auto m_filtration = ft.get_filtration();
    auto m_filter = ft.get_filter();
    auto m_cell = ft.get_m_cell();

    // extract the TB balls
    std::cout << "Writting the TB file, recall that each line contains: " << std::endl
              << "dimension t_radius t_center[3] b_radius b_center[3]" << std::endl;
    std::ofstream file(filename + extension, std::ios::out | std::ios::trunc);
    if (!(file))
    {
        std::cerr << "Error in Persistence<Simplex>::compute_pairs(): impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }
    for( phat::index idx = 0; idx < tb_pairs.get_num_pairs(); idx++ )
    {
        //        std::cout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;
        int t_pos = tb_pairs.get_pair(idx).first;
        int b_pos = tb_pairs.get_pair(idx).second;

        const double t_value = m_filtration.at(t_pos);
        const double b_value = m_filtration.at(b_pos);
        if (t_value < 0 && b_value > 0)
        {
            const int dim = 3-m_filter.at(t_pos).dimension();
            std::tuple<Delaunay::Point, double> t_pr;
            std::tuple<Delaunay::Point, double> b_pr;
            std::clog << "______________\n" << dim << "-hole :" << std::endl;
            if (dim==0) {
                t_pr = ft.search_critical(t_pos, true);// we want a dist max
                b_pr = ft.search_critical(b_pos, false);// we want a dist min
            }
            else if (dim==1) {
                // actually we want saddles for dim1...
                t_pr = std::make_tuple(ft.m_dela.dual(m_cell.at(t_pos)), t_value);
                b_pr = std::make_tuple(ft.m_dela.dual(m_cell.at(b_pos)), b_value);
            }
            else {
                t_pr = ft.search_critical(t_pos, false);// we want a dist min
                b_pr = ft.search_critical(b_pos, true);// we want a dist max
            }

            file << dim << " "
                 << -std::get<1>(t_pr) << " " << std::get<0>(t_pr) << " "
                 << std::get<1>(b_pr) << " " << std::get<0>(b_pr) << " " << std::endl;
        }
    }
    // add the unpaired cell
    file << "0 "
         << - m_filtration.front() << " " << ft.m_dela.dual(m_cell.front()) << " "
         << "inf" << std::endl;

    file.close();
}
*/


//#include "../src/persistence.tpp"

#endif // PERSISTENCE_H
