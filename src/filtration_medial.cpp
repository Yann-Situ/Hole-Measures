/*
FiltrationMedial
computes the medial axes filtrations of a mesh.
It is based on Aldo's code in something.cpp that compute the Delaunay filtration.
*/

#include "filtration_medial.h"
#include "delaunay_helper.h"

#include <fstream>
#include <limits>

int MedialSimplex::currID = 0;

//#############################################################################

FiltrationMedial::FiltrationMedial(Polyhedron poly)
    : m_poly(poly)
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
    std::clog << "nb edges   : " << m_dela.number_of_edges() << " | finite: " << m_dela.number_of_finite_edges() << std::endl;
    std::clog << "nb facets  : " << m_dela.number_of_facets() << " | finite: " << m_dela.number_of_finite_facets()<< std::endl;
    std::clog << "nb cells   : " << m_dela.number_of_cells() << " | finite: " << m_dela.number_of_finite_cells()<< std::endl;
    // CGAL::draw(m_dela);
}


/**
 * @brief FiltrationMedial::init_medial_info
 * Computes the maps medial_type and distance_field and points (ball centers).
 */
void FiltrationMedial::init_medial_info()
{
    // initialize the data structure for telling if a point is inside a mesh
    // cf https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#InsideExample
    SoT inside(m_poly);//CGAL::Side_of_triangle_mesh<Polyhedron, Epick> inside(m_poly);

    const double inf = std::numeric_limits<double>::infinity();
    const MedialInfo medinf_boundary{MedialType::Boundary, 0.0};

    simplex_faces.clear();
    simplex_convert.clear();
    medial_info.clear();
    // (which corresponds to Voronoi vertices)
    for (Delaunay::Cell_handle ch : m_dela.all_cell_handles())
    {
        MedialInfo medinf; // firstly: assign medinf
        if (m_dela.is_infinite(ch))
        {
            medinf.assign(MedialType::Outer, inf, Point(inf,inf,inf));
        }
        else
        {
            const Delaunay::Point dual_point = m_dela.dual(ch);
            const Delaunay::Point point = m_dela.point(ch, 0);
            const double radius = CGAL::sqrt(
                CGAL::squared_distance(dual_point, point));
            const MedialType ch_type =
                (inside(dual_point) == CGAL::ON_BOUNDED_SIDE) ?
                MedialType::Inner : MedialType::Outer;
            medinf.assign(ch_type, radius, dual_point);
        }

        // secondly: run through D-sub-faces
        const std::list<Delaunay::Simplex> faces = DelaunayHelper::D_sub_faces(ch);
        for (Delaunay::Simplex face_D : faces)
        {
            Simplex face;
            auto simplex_find = simplex_convert.find(face_D);
            if (simplex_find == simplex_convert.end())
            {
                // face not found so it is a new face :
                face.assign(face_D);
                simplex_convert[face_D] = face;
                simplex_faces[face] = std::make_pair(Simplices(), Simplices());
                medial_info[face] = medinf;
            }
            else
            {
                //face.assign(simplex_find->second);
                auto face_info_find = medial_info.find(simplex_find->second);
                if (face_info_find->second.m == medinf.m)
                {
                    if (face_info_find->second.r > medinf.r)
                    {
                        // take the min radius (max filtration value)
                        face_info_find->second = medinf;
                    }
                }
                else
                {
                    // Delaunay face on the two sides
                    face_info_find->second = medinf_boundary;
                }
            }
        }
    }
    std::clog << "-- medial_type and distance_field computed: " << simplex_faces.size() << " simplices" << std::endl;
}

void FiltrationMedial::init_simplex_faces()
{
    for (const std::pair<Delaunay::Simplex, Simplex>& simp_pair : simplex_convert)
    {
        const std::list<Delaunay::Simplex> faces =
            DelaunayHelper::D_faces(simp_pair.first);
        for (Delaunay::Simplex face : faces)
        {
            auto it_face_find = simplex_convert.find(face);
            if (it_face_find != simplex_convert.end())
            {
                simplex_faces[simp_pair.second].first.insert(it_face_find->second);
                simplex_faces[it_face_find->second].second.insert(simp_pair.second);
            }
            else
                std::cerr << "Error in init_simplex_faces: face not found" << std::endl;
        }
    }
    std::clog << "-- simplices faces computed" << std::endl;
}

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

void FiltrationMedial::make_filter(MedialType med_type)
{
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


/**
 * @brief FiltrationMedial::compute_delaunay_coboundary
 * Fill the m_coboundary map.
 */
void FiltrationMedial::compute_delaunay_coboundary()
{
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
            if (it_cobound_find != m_coboundary.end())
                it_cobound_find->second.push_back(pos);

        }
        pos++;
    }
    std::clog << "-- Delaunay coboundary computed" << std::endl;
}

/****************************** Work in Progress ******************************/

/**
 * @brief FiltrationMedial::add_critical_to_filter
 * Add every critical simplices (with med_type) to the simplicial complex
 * by subdividing the Delaunay critical cells.
 * @note Should be called before make_filter().
 */

/* not working yet, because we cannot create Delaunay::simplices that are not
 * part of the delaunay triangulation... That's harder that it seems.
 */
void FiltrationMedial::add_critical_to_filter(MedialType med_type)
{
    //std::vector<std::pair<Simplex, MedialInfo>> to_add;
    int critcount=0;
    int addcount=0;
    for (const std::pair<Delaunay::Simplex, Simplex>& simp_pair : simplex_convert)
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
                /* in the filtration, simp arrives later than the max
                 * filtration value of its Voronoi faces. It's the case for
                 * df local min (e.g breadth point of 0-hole, thickness
                 * point of 2-hole, T and B point of 1-hole).

                 * @conjecture the bad cases, where the Voronoi
                 * vertices doesn't approximate topologically critical points,
                 * are exactly the cases where the topologically critical point
                 * is a local df_med saddle of index <= D-1 = 1
                 * (D = dim_medialaxes).
                 * A consequence is that we only need to add the Delaunay
                 * critical points that are df_med saddle point.
                 */
                medinfsimp = medinfcrit;
                addcount++;
            }
            else
            {
                /* the cell arrives before in the filtration than the max
                 * filtration value of its Voronoi faces */
            }
        }
    }
    std::clog << addcount << " adds and " << critcount << " crits" << '\n';
}

/* Some code I tried to implement to add simplices to the filtration/triangulation
 * in order to handle topologically critical points.
if (simp.dimension() == 1) // voronoi facet
{


    const Delaunay::Edge e(simp); // 2-hole : Voronoi facet
    const Delaunay::Point p0 = m_dela.point(e.first->vertex(e.second));
    const Delaunay::Point p1 = m_dela.point(e.first->vertex(e.third ));
}
else if (simp.dimension() == 2) // voronoi edge
{
    assert(simplex_faces[simp].second.size() == 2);
    // the two 3Dsimplex (Voronoi vertices) of the Voronoi boundary
    Simplex& p0 = simplex_faces[simp].second[0];
    Simplex& p1 = simplex_faces[simp].second[1];
    MedialInfo& medinfp0 = medial_info[p0];
    MedialInfo& medinfp1 = medial_info[p1];

    Simplex p(3); MedialInfo medinfp(med_type, crit.r, crit.p);
    // new Dsimplex with dimension 3
    medial_info[p] = medinfp;

    Simplex e0(2);
    // take the minimum radius:
    medial_info[e0] = (medinfp0.r > medinfp.r) ? medinfp : medinfp0;

    Simplex e1(2);
    // take the minimum radius:
    medial_info[e1] = (medinfp1.r > medinfp.r) ? medinfp : medinfp1;

    simplex_faces[p] = std::make_pair(Simplices{e0,e1},Simplices()); // Dfaces / Dcofaces
    simplex_faces[e0] = std::make_pair(Simplices(),Simplices()); // Dfaces / Dcofaces
    simplex_faces[e1] = std::make_pair(Simplices(),Simplices()); // Dfaces / Dcofaces
    for (Simplex Dface : simplex_faces[simp].first) // faces
    {
        std::pair<Simplices, Simplices>& Dface_pair = simplex_faces[Dface];
        Dface_find.second.erase(simp); // erase from Dface cofaces
        Dface_find.second.insert(e0);
        Dface_find.second.insert(e1); // update medial info: take the min ?
    }
    for (Simplex Dcoface : simplex_faces[simp].second) // cofaces
    {
        std::pair<Simplices, Simplices>& Dcoface_pair = simplex_faces[Dface];
        Dcoface_find.first.erase(simp); // erase from Dcoface faces
    }

}
simplex_faces.erase(simp);
medial_info.erase(it_med);
*/
