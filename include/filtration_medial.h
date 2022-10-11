#ifndef FILTRATION_MEDIAL_H
#define FILTRATION_MEDIAL_H


#include "cgal_typedef.h"
#include "filtration.h"
/**
 * @brief The FiltrationMedial class
 * Contains functions to create and handle the voronoi/delaunay filtration
 * of a mesh file .off.
 */

//enum class CellCritType {NonCritical, Critical}; // moved to delauney_helper.h
enum class MedialType {Boundary, Inner, Outer};

struct MedialInfo
{
    MedialType m;
    double r; // distance to border
    Delaunay::Point p; // ball center associated to the simplex

    MedialInfo(MedialType medtype = MedialType::Boundary,
        double dist = 0.0, Delaunay::Point pt = Delaunay::Point(0.0,0.0,0.0))
        : m(medtype), r(dist), p(pt) {}

    void assign(MedialType medtype = MedialType::Boundary,
        double dist = 0.0, Delaunay::Point pt = Delaunay::Point(0.0,0.0,0.0))
        { m = medtype; r = dist; p = pt;}
};

class MedialSimplex
{
private:
    static int currID;
    int dim;
    int ID;

public:
    MedialSimplex() : ID(-1), dim(-1) {}
    MedialSimplex(int _dim) : ID(currID++), dim(_dim) {}
    MedialSimplex(const MedialSimplex& m) : ID(m.id()), dim(m.dimension()) {}
    //MedialSimplex(Delaunay::Simplex s) : ID(currID++), dim(s.dimension()) {}

    void assign(Delaunay::Simplex s) {ID = currID++; dim = s.dimension();};
    void assign(const MedialSimplex& m) {ID = m.id(); dim = m.dimension();};
    //void addface(int face_id) {faces_id.push_back(face_id);}
    int dimension() const {return dim;}
    int id() const {return ID;}
};
inline bool operator<(const MedialSimplex& lhs, const MedialSimplex& rhs)
{
    return lhs.id() < rhs.id();
}

class FiltrationMedial : public Filtration<MedialSimplex>
{
public:
    typedef MedialSimplex Simplex;
    typedef std::set<Simplex> Simplices;
    FiltrationMedial(Polyhedron poly);//const std::string &filename);

    void init_medial_info();
    void init_simplex_faces();
    void make_filter(MedialType med_type);
    void add_critical_to_filter(MedialType med_type);
    void compute_delaunay_coboundary();

    Polyhedron m_poly;
    Delaunay m_dela;

    // Filtration function
    int get_filter_size() {return filter.size();}
    Simplex get_filter(int pos) {return filter.at(pos).first;}
    Delaunay::Point get_point(int pos) {return filter.at(pos).second.p;}
    double get_filtration(int pos)    {return -filter.at(pos).second.r;}
    const std::map<Simplex, std::vector<int>>& get_coboundary() {return m_coboundary ;}

private:
    std::map<Simplex, std::vector<int>> m_coboundary; /// map simplex -> the list of its Delaunay coboudary cell positions
    std::vector<std::pair<Simplex, MedialInfo>> filter;
    std::map<Simplex, MedialInfo> medial_info; /// map simplex -> medial information

    std::map<Simplex, std::pair<Simplices, Simplices>> simplex_faces; // Dfaces and Dcofaces
    std::map<Delaunay::Simplex, Simplex> simplex_convert; /// map simplex -> medial information
    // std::map doesn't work with cell_handle / vertex_handle... that's a pity...
};


#endif // FILTRATION_MEDIAL_H
