#ifndef FILTRATION_MEDIAL_H
#define FILTRATION_MEDIAL_H


#include "cgal_typedef.h"
#include "filtration.h"
#include "delaunay_helper.h"

#include <fstream>
#include <limits>
/**
 * @brief The FiltrationMedial class
 * Contains functions to create and handle the voronoi/delaunay filtration
 * of a mesh file .off.
 */

//enum class CellCritType {NonCritical, Critical}; // moved to delauney_helper.h
enum MedialType {Boundary=0, Inner=1, Outer=2};

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

/**
 * @brief FiltrationMedial::init_finite_filtration_info
 * Computes the maps medial_type and distance_field and points (ball centers)
 * for finite simplices.
 * Also init the map simplex_faces, that contains the Delaunay faces and cofaces
 * (also only for finite simplices).
 * This function MUST be called before init_infinite_filtration_info.
 */
 struct before_than_face_value_order
 {
     // to compare pairs of (medial_info, Dcell)
     inline bool operator()
         (const std::pair<MedialInfo, Delaunay::Cell_handle>& sd1,
          const std::pair<MedialInfo, Delaunay::Cell_handle>& sd2)
     {
         return (sd1.first.r < sd2.first.r);
     }
 };


class FiltrationMedial : public Filtration<MedialSimplex>
{
public:
    typedef MedialSimplex Simplex;
    typedef std::set<Simplex> Simplices;
    /**
     * @brief FiltrationMedial::make_filter
     * Create the filter and the medial axis filtration by filling m_filter :
     * the simplices are sorted according to the filter induced by the
     * before_than_filtration_order, which corresponds to decreasing distance to
     * boundary.
     */
    struct before_than_filtration_order
    {
        // to compare pairs of (simplex, distance_field)
        inline bool operator()
            (const std::pair<FiltrationMedial::Simplex, MedialInfo>& sd1,
             const std::pair<FiltrationMedial::Simplex, MedialInfo>& sd2)
        {
            return (sd1.second.r > sd2.second.r) || ((sd1.second.r == sd2.second.r)
                && (sd1.first.dimension() > sd2.first.dimension())); // warning, it's the delaunay simplex dimension
        }
    };
    
    FiltrationMedial(Polyhedron poly) : m_poly(poly)
    {
        //CGAL::draw(m_poly);
        //CGAL::Side_of_triangle_mesh<Polyhedron, Epick> inside(m_poly);
        for (Polyhedron::Vertex_iterator it = m_poly.vertices_begin(); it != m_poly.vertices_end(); ++it)
        {
            m_dela.insert(it->point()); // m_dela corresponds to the 3D Delaunay of the surface points
    //        std::cout << it->point() << std::endl;
        }
        std::clog << "Computing persistence using Medial axes filtration" << std::endl;
        std::clog << "nb finite vertices: " << m_dela.number_of_vertices() << std::endl;
        std::clog << "nb edges   : " << m_dela.number_of_edges()  << " | finite: " << m_dela.number_of_finite_edges() << std::endl;
        std::clog << "nb facets  : " << m_dela.number_of_facets() << " | finite: " << m_dela.number_of_finite_facets()<< std::endl;
        std::clog << "nb cells   : " << m_dela.number_of_cells()  << " | finite: " << m_dela.number_of_finite_cells()<< std::endl;
        // CGAL::draw(m_dela);
    }

    /**
     * @brief FiltrationMedial::get_medial_count
     * Count the boundary, inner and outer medial axis size in terms of simplices.
     */
    std::array<int,3> get_medial_count() {
        std::array<int, 3> count{ {0, 0, 0} };
        for (const std::pair<Simplex, MedialInfo>& simp_med : medial_info)
        {
            count[simp_med.second.m] += 1;
        }
        return count;
    }

    void init_finite_filtration_info() {
        // initialize the data structure for telling if a point is inside a mesh
        // cf https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#InsideExample
        SoT inside(m_poly);//CGAL::Side_of_triangle_mesh<Polyhedron, Epick> inside(m_poly);
        const MedialInfo medinf_boundary{MedialType::Boundary, 0.0};

        // INFINITE_OPTIMISATION: here use number_of_finite_cells instead of number_of_cells.
        // std::vector<std::pair<MedialInfo, Delaunay::Cell_handle>> medinf_Dcell_pairs(m_dela.number_of_cells()); // (which corresponds to Voronoi vertices)
        std::vector<std::pair<MedialInfo, Delaunay::Cell_handle>> medinf_Dcell_pairs(m_dela.number_of_finite_cells()); // (which corresponds to Voronoi vertices)
        int i = 0;

        // INFINITE_OPTIMISATION: here use finite_cell_handles instead of all_cell_handles
        // for (Delaunay::Cell_handle ch : m_dela.all_cell_handles()) // also the infinite ones
        for (Delaunay::Cell_handle ch : m_dela.finite_cell_handles())
        {
            // INFINITE_OPTIMISATION: here we don't need to check if the cell is infinite anymore
            MedialInfo medinf;
            const Delaunay::Point dual_point = m_dela.dual(ch);
            const Delaunay::Point point = m_dela.point(ch, 0);
            const double radius = CGAL::sqrt(
                                             CGAL::squared_distance(dual_point, point));
            const MedialType ch_type =
            (inside(dual_point) == CGAL::ON_BOUNDED_SIDE) ?
            MedialType::Inner : MedialType::Outer;
            medinf.assign(ch_type, radius, dual_point);
            medinf_Dcell_pairs.at(i) = std::make_pair(medinf, ch);
            i++;
        }
        std::sort(medinf_Dcell_pairs.begin(), medinf_Dcell_pairs.end(),
                  before_than_face_value_order());
        // farthest from the border must be last (reverse order from the filtration)

        simplex_faces.clear();
        finite_simplex_convert.clear();
        medial_info.clear();
        medialtype_count.clear();

        for (auto medinf_Dcell_pair : medinf_Dcell_pairs) // for every Dcell (in decreasing sdf order)
        {
            const Delaunay::Cell_handle ch = medinf_Dcell_pair.second; // current Dcell
            MedialInfo medinf = medinf_Dcell_pair.first;
            // secondly: run through D-sub-faces
            const std::list<Delaunay::Simplex> D_sub_faces = DelaunayHelper::D_sub_faces(ch);
            for (Delaunay::Simplex D_sub_face : D_sub_faces)
            {
                auto simplex_find = finite_simplex_convert.find(D_sub_face);
                if (simplex_find == finite_simplex_convert.end())
                {
                    // face not found so it is a new face :
                    Simplex face;
                    face.assign(D_sub_face);
                    finite_simplex_convert[D_sub_face] = face;
                    simplex_faces[face] = std::make_pair(Simplices(), Simplices());
                    medial_info[face] = medinf;
                }
                else
                {
                    /* WARNING :
                     * pairs is in increasing df order because we want to add a
                     * Delaunay cell at the same time as the boundary cell that
                     * minimize df.
                     * As a result we don't need to change the previous medinfo */
                    auto face_info_find = medial_info.find(simplex_find->second);
                    if (face_info_find->second.m != medinf.m)
                    {
                        // Delaunay face on the two sides
                        face_info_find->second = medinf_boundary;
                    }
                }
            }
        }

        for (const std::pair<Delaunay::Simplex, Simplex>& simp_pair : finite_simplex_convert)
            // INFINITE_OPTIMISATION: here do not deal with infinite simplices anymore.
        {
            const std::list<Delaunay::Simplex> D_faces =
            DelaunayHelper::D_faces(simp_pair.first);
            for (Delaunay::Simplex D_face : D_faces)
            {
                auto it_face_find = finite_simplex_convert.find(D_face);
                if (it_face_find != finite_simplex_convert.end())
                {
                    simplex_faces[simp_pair.second].first.insert(it_face_find->second);
                    simplex_faces[it_face_find->second].second.insert(simp_pair.second);
                }
                else
                    std::cerr << "Error in init_finite_filtration_info: face not found" << std::endl;
            }
        }

        std::clog << "-- finite filtration information computed: " << simplex_faces.size() << " finite simplices" << std::endl;
    }

    /**
     * @brief FiltrationMedial::init_infinite_filtration_info
     * Create the infinite D-cell and update the information of the surrounding simplices.
     * update medial_info and simplex_faces.
     * This function MUST be called after init_finite_filtration_info.
     */
    void init_infinite_filtration_info() {
        // INFINITE_OPTIMISATION: here create the infinite simplex and update its
        // faces, which corresponds to all the finite faces of all the infinite simplices.
        // update simplex_faces, medial_info but not finite_simplex_convert
        const double inf = std::numeric_limits<double>::infinity();
        const MedialInfo medinf_boundary{MedialType::Boundary, 0.0};
        const MedialInfo infinite_medinf = MedialInfo(MedialType::Outer, inf, Point(inf,inf,inf));
        const Simplex infinite_medsimplex = Simplex(3);

        simplex_faces[infinite_medsimplex] = std::make_pair(Simplices(), Simplices());
        medial_info[infinite_medsimplex] = infinite_medinf;

        Delaunay::Vertex_handle inf_Dvertex = m_dela.infinite_vertex();
        std::list<Delaunay::Cell_handle> inf_Dcells;
        m_dela.incident_cells(inf_Dvertex,std::back_inserter(inf_Dcells));
        std::list<Delaunay::Simplex> D_finite_faces;

        for (Delaunay::Cell_handle inf_Dcell : inf_Dcells)
        {
            const std::list<Delaunay::Simplex> D_finite_sub_faces
                  = DelaunayHelper::D_finite_sub_faces(m_dela, inf_Dcell);
            for (Delaunay::Simplex D_sub_face : D_finite_sub_faces)
            {
                auto simplex_find = finite_simplex_convert.find(D_sub_face);
                if (simplex_find == finite_simplex_convert.end())
                {std::cerr << "Error in init_infinite_filtration_info: face not found" << std::endl;}
                else
                {
                    auto face_info_find = medial_info.find(simplex_find->second);
                    if (face_info_find->second.m != infinite_medinf.m)
                    {
                        // Delaunay face on the two sides
                        face_info_find->second = medinf_boundary;
                    }
                }
                if (D_sub_face.dimension() == 2) // if it is a finite facet of an infinite cell
                {
                    /* If the triangulation is a manifold, every finite facet should be
                     * encountered only once because a finite facet must have maximum a single infinite
                     * incident cell. Here we implemented the risky version, where we
                     * suppose that a facet can't have two infinite incident cell:*/
                    D_finite_faces.push_back(D_sub_face);
                }
            }
        }

        // INFINITE_OPTIMISATION: here compute the simplex_faces value of the
        // infinite simplex, and update the othe simplex_faces lists.
        for (Delaunay::Simplex D_face : D_finite_faces)
        {
            auto it_face_find = finite_simplex_convert.find(D_face);
            if (it_face_find != finite_simplex_convert.end())
            {
                simplex_faces[infinite_medsimplex].first.insert(it_face_find->second);
                simplex_faces[it_face_find->second].second.insert(infinite_medsimplex);
            }
            else
                std::cerr << "Error in init_infinite_filtration_info: face not found" << std::endl;
        }
        std::clog << "-- infinite filtration information computed: " << simplex_faces.size() << " simplices" << std::endl;

    }

    void make_filter(MedialType med_type) {
        // We add the med_type simplices
        filter.clear();
        for (const std::pair<Simplex, MedialInfo>& simp_info : medial_info)
        {
            if (simp_info.second.m == med_type)
            {
                filter.push_back(simp_info);
            }
        }
        // sort according to the filtration order
        std::sort(filter.begin(), filter.end(), before_than_filtration_order());
        std::clog << "-- filter computed: filter size is " << filter.size() << " simplices" << std::endl;
    }

    /****************************** Work in Progress ******************************/

    /**
     * @brief FiltrationMedial::add_critical_to_filter
     * Add every critical simplices (with med_type) to the simplicial complex
     * by subdividing the Delaunay critical cells.
     * @note Should be called before make_filter().
     */
    void add_critical_to_filter(MedialType med_type) {
        //std::vector<std::pair<Simplex, MedialInfo>> to_add;
        int critcount=0;
        int addcount=0;
        for (const std::pair<Delaunay::Simplex, Simplex>& simp_pair : finite_simplex_convert)
        {
            auto it_med = medial_info.find(simp_pair.second);
            if (it_med == medial_info.end())
            {
                std::cerr << "add_critical_to_filter error : medial_info simplex not find" << '\n';
                continue;
            }
            if (it_med->second.m != med_type){continue;}

            const Simplex& simp = it_med->first;
            const CriticalInfo crit =
            DelaunayHelper::get_critical_info(m_dela, simp_pair.first);
            /* struct CriticalInfo
            * {
            *    CriticalType c;  // enum CriticalType {NonCritical, Critical};
            *    double r; // distance to border
            *    Delaunay::Point p; // ball center associated to the simplex
            * }
            */
            if (crit.c == CriticalType::Critical)
            {
                critcount++;
                MedialInfo& medinfsimp = medial_info[simp];
                MedialInfo medinfcrit(med_type, crit.r, crit.p);
                if (crit.r <= medinfsimp.r)
                {
                    /* Easiest case, where we don't need to add new simplices:
                     * we just need to delay the critical cell.
                     */
                    /* In the filtration, simp arrives later than the max
                     * filtration value of its Voronoi faces. It's the case for
                     * df local min (e.g breadth point of 0-hole, thickness
                     * point of 2-hole, T and B point of 1-hole).

                     * @conjecture: in 3D, the bad cases (when the Voronoi
                     * vertices doesn't approximate topologically critical points)
                     * are exactly the cases where the topologically critical point
                     * is a local df_med saddle of index <= D-1 = 1
                     * (D = dim_medialaxes).
                     *
                     * To sum up, bad approximation cases are easy to handle for
                     * meidal axis filtration (we only need to change the
                     * medial_info of the critical simplex). In other cases,
                     * refining the mesh will converge to the rights spot.
                     */
                    medinfsimp = medinfcrit;
                    addcount++;
                }
                else
                {
                    /* the cell arrives before in the filtration than the max
                     * filtration value of its Voronoi faces:
                     * we need to create a subdivision of the current Voronoi
                     * element to add its critical point.
                     */
                }
            }
        }
        std::clog << "-- " << addcount << " adds and " << critcount << " crits" << '\n';
    }

    /**
     * @brief FiltrationMedial::compute_delaunay_coboundary
     * Fill the m_coboundary map.
     */
    void compute_delaunay_coboundary() {
        m_coboundary.clear();
        for (const auto& pair : filter) // decreasing sdf order
        {
            std::vector<int> c; // empty vector of coboundary for the moment.
            /* Coboundary corresponds to Delaunay coboundary,
             * which are duals of Voronoi boundary*/

            //std::clog << pair.first.id() << "\t";
            m_coboundary[pair.first] = c;
        }

        /* Filling the coboundary vectors, which contains the positions of
         * Delaunay coboundary cells in m_filter */
        int pos = 0;
        for (const auto& pair : filter)
        {
            for (Simplex face : simplex_faces[pair.first].first)
            {
                auto it_cobound_find = m_coboundary.find(face);
                if (it_cobound_find != m_coboundary.end()) // if the face is in the filtration
                    it_cobound_find->second.push_back(pos);

            }
            pos++;
        }
        std::clog << "-- Delaunay coboundary computed" << std::endl;
    }


    // Member data
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
    std::map<Delaunay::Simplex, Simplex> finite_simplex_convert;
    std::vector<int> medialtype_count;
    // std::map doesn't work with cell_handle / vertex_handle... that's a pity...
};

int MedialSimplex::currID = 0;




#endif // FILTRATION_MEDIAL_H
