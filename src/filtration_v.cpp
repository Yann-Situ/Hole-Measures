/*
FiltrationV
computes the Voronoi filtration of a mesh.
It is based on Aldo's code in something.cpp that compute the Delaunay filtration.
*/

#include "filtration_v.h"
#include "delaunay_helper.h"

//#include <CGAL/draw_polyhedron.h>
//#include <CGAL/draw_triangulation_3.h>
#include <fstream>


inline int sign(double d){ return (d>=0.0) ? 1 : -1;}
inline double max3(double d1, double d2, double d3)
{
    if (d1 >= d2)
    {
        if (d1 >= d3) {return d1;}
        else {return d3;}
    }
    else
    {
        if (d2 >= d3) {return d2;}
        else {return d3;}
    }
}

//#############################################################################

FiltrationV::FiltrationV(Polyhedron poly)
    : m_poly(poly)
{
    //CGAL::draw(m_poly);
    //CGAL::Side_of_triangle_mesh<Polyhedron, Epick> inside(m_poly);
    for (Polyhedron::Vertex_iterator it = m_poly.vertices_begin(); it != m_poly.vertices_end(); ++it)
    {
        m_dela.insert(it->point()); // m_dela corresponds to the 3D Delaunay of the surface points
//        std::cout << it->point() << std::endl;
    }
    std::clog << "Computing persistence using Voronoi filtration" << std::endl;
    std::clog << "nb finite vertices: " << m_dela.number_of_vertices() << std::endl;
    std::clog << "nb edges   : " << m_dela.number_of_edges() << " | finite: " << m_dela.number_of_finite_edges() << std::endl;
    std::clog << "nb facets  : " << m_dela.number_of_facets() << " | finite: " << m_dela.number_of_finite_facets()<< std::endl;
    std::clog << "nb cells   : " << m_dela.number_of_cells() << " | finite: " << m_dela.number_of_finite_cells()<< std::endl;

    // CGAL::draw(m_dela);
}


/**
 * @brief FiltrationV::init_Dcell_filtration_values
 * Assign filtration values to the cells in m_cell_filtration.
 */
void FiltrationV::init_Dcell_filtration_values()
{
    m_cell_filtration.clear();

    // initialize the data structure for telling if a point is inside a mesh
    // cf https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#InsideExample
    SoT inside(m_poly);//CGAL::Side_of_triangle_mesh<Polyhedron, Epick> inside(m_poly);

    // compute filtration values at the Delaunay cells (= 3-simplices)
    // (which corresponds to Voronoi vertices)
    for (Delaunay::Cell_handle ch : m_dela.finite_cell_handles())
    {
        const Delaunay::Point dual_point = m_dela.dual(ch);
        const Delaunay::Point point = m_dela.point(ch, 0);
        double radius = CGAL::sqrt(CGAL::squared_distance(dual_point, point));
        CGAL::Bounded_side res = inside(dual_point);
        if (res == CGAL::ON_BOUNDED_SIDE)
            radius *= -1;
        m_cell_filtration[ch] = radius;
//        std::clog << "tet=" << m_dela.tetrahedron(ch) << "; center=" << m_dela.dual(ch) << "; r=" << radius << std::endl;
    }
    std::clog << "-- filtration values at Delaunay cells computed" << std::endl;
}

/**
 * @brief FiltrationV::make_filter
 * Create the filter and the filtration by filling m_filter and m_filtration
 * in the appropriate order : the simplices are sorted according to the
 * filter induced by the filtration values on the Dcells (3D-Delaunay-simplices).
 * Also initialize m_coboundary with empty lists.
 */
void FiltrationV::make_filter()
{
    // sort cells by the filtration values
    std::vector<std::pair<double, Delaunay::Cell_handle>> pairs(m_cell_filtration.size());
    int i = 0;
    for (const auto& kv : m_cell_filtration)
    {
        pairs.at(i) = std::make_pair(kv.second, kv.first);
        i++;
//        std::clog << kv.first << " has value " << kv.second << std::endl;
    }
    std::sort(pairs.rbegin(), pairs.rend()); // in reverse order

    //m_pos_in_filter.clear();
    m_filter.clear();
    m_filtration.clear();
    //m_point.clear();
    m_cell.clear();
    m_coboundary.clear();
    for (auto pair : pairs) // decreasing sdf order
    {
        const Delaunay::Cell_handle ch = pair.second; // current cell (3-simplex)
        //const Delaunay::Point point = m_dela.dual(ch); // its corresponding point of the Voronoi diagram
        const std::list<Simplex> simplices = DelaunayHelper::D_sub_faces(ch);

        for (Simplex simplex : simplices)
        {
            if (m_coboundary.count(simplex) == 0 // if that simplex is not already added
            && infinite_voronoi_face_duals.count(simplex) == 0) // if its dual is finite
            {
                //m_pos_in_filter[simplex] = 1; // useless for Voronoi filtration
                m_filter.push_front(simplex);
                /* WARNING push_front :
                 * As pairs is in decreasing sdf order, we want the last simplex
                 * of pair (with lowest sdf) to be the first of m_filter.
                 *
                 * pairs is in decreasing sdf in order to add a Delaunay cell at
                 * the same time of its boundary cell that maximize sdf. */
                m_filtration.push_front(pair.first);
                //m_point.push_front(point);
                m_cell.push_front(ch);
                std::vector<int> c; // empty vector of coboundary for the moment.
                /*coboundary corresponds to Delaunay coboundary, which are duals of Voronoi boundary*/
                m_coboundary[simplex] = c;
            }
        }
    }
    std::clog << "-- filter computed" << std::endl;
}

/**
 * @brief FiltrationV::compute_infinite_voronoi_face_duals
 * Set up the set infinite_voronoi_face_duals that should contains the duals
 * of the infinite Voronoi cells.
 * These Delaunay cells are the infinite ones and those on the convex hull.
 */
void FiltrationV::compute_infinite_voronoi_face_duals()
{
    infinite_voronoi_face_duals.clear();
    Delaunay::Vertex_handle inf_Dvertex = m_dela.infinite_vertex();
    std::list<Delaunay::Cell_handle> inf_Dcells;
    m_dela.incident_cells(inf_Dvertex,std::back_inserter(inf_Dcells));
    //std::cout << "IIV : " << inf_Dcells.size() << std::endl;
    for (Delaunay::Cell_handle inf_Dcell : inf_Dcells)
    {
        const std::list<Simplex> simplices = DelaunayHelper::D_sub_faces(inf_Dcell);
        for (Simplex simplex : simplices)
        {
            infinite_voronoi_face_duals.insert(simplex);
        }
    }
    std::clog << "-- infinite Voronoi face duals computed : " << infinite_voronoi_face_duals.size() << std::endl;
    //is_infinite (Cell_handle c) const ;
    //all_faces(Cell_handle c)
}

/**
 * @brief FiltrationV::compute_delaunay_coboundary
 * Fill the m_coboundary map.
 * @Precondition : m_coboundary should have been instantiated with empty vectors.
 */
void FiltrationV::compute_delaunay_coboundary()
{
    int pos = 0;
    /* Filling the coboundary vectors, which contains the positions of
     * Delaunay coboundary cells in m_filter */
    for (Simplex simplex : m_filter)
    {
        //m_pos_in_filter[simplex] = pos;
        const std::list<Simplex> faces = DelaunayHelper::D_faces(simplex);
        for (Simplex face : faces)
        {
            m_coboundary[face].push_back(pos);
        }
        pos++;
    }
    std::clog << "-- Delaunay coboundary computed" << std::endl;
}

/**
 * @brief FiltrationV::is_on_boundary
 * Check if the facet is on the boundary of the object.
 * @Precondition : m_cell_filtration should have been filled.
 */
bool FiltrationV::is_on_boundary(const Delaunay::Facet f)
{
    Delaunay::Cell_handle c1 = f.first;
    Delaunay::Cell_handle c2 = m_dela.mirror_facet(f).first;
    bool inf1=m_dela.is_infinite(c1); bool inf2=m_dela.is_infinite(c2);
    if (inf1 != inf2)
    {
        return true;
    }
    else if (inf1) // case where both cells are infinite
    {
        return false;
    }
    return (sign(m_cell_filtration[c1]) != sign(m_cell_filtration[c2]));
}
/**
 * @brief FiltrationV::boundary_facets
 * Return a deque of the boundary facets.
 * @Precondition : m_cell_filtration and infinite_voronoi_face_duals should
 * have been filled.
 */
std::deque<Delaunay::Facet> FiltrationV::boundary_facets()
{
    std::deque<Delaunay::Facet> b_facets;
    for (Delaunay::Facet f : m_dela.finite_facets())
    {
        if (is_on_boundary(f))
            b_facets.push_back(f);
    }
    return b_facets;
}

/**
 * @brief FiltrationV::is_finite_medial
 * Check if the simplex represents a part of the
 * finite-Voronoi-medial-axis-approx. More precisely, check if the Delaunay
 * simplex is not on the boundary and is not in infinite_voronoi_face_duals
 * (i.e on the convex hull and on infinite cells).
 * @Precondition : m_cell_filtration and infinite_voronoi_face_duals should
 * have been filled.
 */
bool FiltrationV::is_finite_medial(const Simplex s)
{
    //SoT inside(m_poly);
    if (infinite_voronoi_face_duals.count(s) > 0)
    {
        //std::clog << "is_on_boundary Warning : should not be called on infinite_voronoi_face_duals cells. return false" << std::endl;
        return false;
    }
    if (s.dimension() == 0)
    {
        return false; // it's a Delaunay vertex.
    }
    if (s.dimension() == 1)
    {
        const Delaunay::Edge e(s); // 2-hole : Voronoi facet
        Delaunay::Cell_circulator circ = m_dela.incident_cells(e), past_end(circ);
        int res_save = sign(m_cell_filtration[e.first]);
        do
        {
            Delaunay::Cell_handle c = circ;
            if (m_dela.is_infinite(c))
            {
                std::clog << "is_on_boundary Error : should not be called on infinite cells. return false" << std::endl;
                return false;
            }
            int res_new = sign(m_cell_filtration[c]);
            if (res_new != res_save)
            {
                return false;
            }
        } while (++circ != past_end);
    }
    else if (s.dimension() == 2)
    {
        const Delaunay::Facet f(s); // 1-hole : Voronoi edge
        Delaunay::Cell_handle c1 = f.first;
        Delaunay::Cell_handle c2 = m_dela.mirror_facet(f).first;
        if (m_dela.is_infinite(c1) || m_dela.is_infinite(c2))
        {
            std::clog << "is_on_boundary Error : should not be called on infinite cells. return false" << std::endl;
            return false;
        }
        int res1 = sign(m_cell_filtration[c1]);
        int res2 = sign(m_cell_filtration[c2]);
        if (res1 != res2)
        {
            return false;
        }
    }
    else if (s.dimension() == 3)
    {
        // it's a voronoi vertex, return true
    }
    else
    {
        std::clog << "is_on_boundary Error : found a cell with dim = " << s.dimension() << ". return false" << std::endl;
    }
    return true;
}


/**
 * @brief FiltrationV::search_critical
 * Given a position of a cell in the filter, compute all the subfaces of its
 * parent Delaunay cell_handle (3-cell) and return the critical elements of
 * the subface that maximize (or minimize if !use_max) the distance to border
 * value of its potential critical point.
 * @Precondition : m_cell, m_cell_filtration and infinite_voronoi_face_duals
 * should have been filled.
 */
std::pair<Delaunay::Point, double> FiltrationV::search_critical(int pos, bool use_max)
{
    int parentDcellpos = pos; Delaunay::Cell_handle parentDcell = m_cell[pos];
    const double r = m_filtration.at(pos);
    std::list<Simplex> faces = DelaunayHelper::D_sub_faces(parentDcell);

    double rmax = r; Delaunay::Point pmax = m_dela.dual(parentDcell);
    bool inside = (r<0);
    for (auto s : faces) {
        if (is_finite_medial(s))
        {
            const CriticalInfo crit = DelaunayHelper::get_critical_info(m_dela, s);
            if (crit.c == CriticalType::Critical)
            {
                if ( (use_max && rmax < crit.r) || ((!use_max) && rmax > crit.r) )
                {
                    pmax = crit.p;
                    rmax = crit.r;
                }
            }
        }
    }
    std::clog << "critical search in " << faces.size() << " cells of filtration value = " << r << std::endl;
    return std::make_pair(pmax, rmax);
}
