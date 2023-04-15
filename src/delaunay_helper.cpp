/*
DelaunayHelper
*/

#include "delaunay_helper.h"

/**
 * @brief DelaunayHelper::D_faces
 * Return the list of the (d-1)-simplices in the boundary of the d-simplexe s.
 */
std::list<Delaunay::Simplex> DelaunayHelper::D_faces(const Delaunay::Simplex s)
{
    std::list<Delaunay::Simplex> l;
    if (s.dimension() == 0)
    {

    }
    else if (s.dimension() == 1)
    {
        const Delaunay::Edge e(s); // cf https://doc.cgal.org/latest/Triangulation_3/classCGAL_1_1Triangulation__simplex__3.html#a8b2d5cb3dac0d210112cb4fd7f1ae715
        l.push_back(Delaunay::Simplex(e.first->vertex(e.second)));
        l.push_back(Delaunay::Simplex(e.first->vertex(e.third )));
    }
    else if (s.dimension() == 2)
    {
        const Delaunay::Facet f(s);
        l.push_back(Delaunay::Simplex(Delaunay::Edge(f.first, (f.second+1)%4, (f.second+2)%4)));
        l.push_back(Delaunay::Simplex(Delaunay::Edge(f.first, (f.second+2)%4, (f.second+3)%4)));
        l.push_back(Delaunay::Simplex(Delaunay::Edge(f.first, (f.second+3)%4, (f.second+1)%4)));
    }
    else if (s.dimension() == 3)
    {
        const Delaunay::Cell_handle ch(s);
        l.push_back(Delaunay::Simplex(Delaunay::Facet(ch, 0)));
        l.push_back(Delaunay::Simplex(Delaunay::Facet(ch, 1)));
        l.push_back(Delaunay::Simplex(Delaunay::Facet(ch, 2)));
        l.push_back(Delaunay::Simplex(Delaunay::Facet(ch, 3)));
    }
    return l;
}

/**
 * @brief DelaunayHelper::D_sub_faces
 * Return the list of every simplexe in the iterated boundary of a Delaunay object s.
 * Precisely return $\partial^0(s)\cup\partial^1(s)\cup ... \cup\partial^n(s)\cup$
 * where $\partial$ is the boundary operator.
 * Note that the simplexe corresponding to s is also in the list.
 */
std::list<Delaunay::Simplex> DelaunayHelper::D_sub_faces(const Delaunay::Cell_handle ch)
{
    std::list<Delaunay::Simplex> l;
    l.push_back(Delaunay::Simplex(ch->vertex(0)));
    l.push_back(Delaunay::Simplex(ch->vertex(1)));
    l.push_back(Delaunay::Simplex(ch->vertex(2)));
    l.push_back(Delaunay::Simplex(ch->vertex(3)));
    l.push_back(Delaunay::Simplex(Delaunay::Edge(ch, 2, 3)));
    l.push_back(Delaunay::Simplex(Delaunay::Edge(ch, 1, 3)));
    l.push_back(Delaunay::Simplex(Delaunay::Edge(ch, 0, 3)));
    l.push_back(Delaunay::Simplex(Delaunay::Edge(ch, 1, 2)));
    l.push_back(Delaunay::Simplex(Delaunay::Edge(ch, 0, 2)));
    l.push_back(Delaunay::Simplex(Delaunay::Edge(ch, 0, 1)));
    l.push_back(Delaunay::Simplex(Delaunay::Facet(ch, 1)));
    l.push_back(Delaunay::Simplex(Delaunay::Facet(ch, 0)));
    l.push_back(Delaunay::Simplex(Delaunay::Facet(ch, 2)));
    l.push_back(Delaunay::Simplex(Delaunay::Facet(ch, 3)));
    l.push_back(Delaunay::Simplex(ch));
    return l;
}


std::list<Delaunay::Simplex> DelaunayHelper::D_sub_faces(const Delaunay::Facet f)
{
    std::list<Delaunay::Simplex> l;
    Delaunay::Cell_handle c = f.first;
    int i = f.second;
    l.push_back(Delaunay::Simplex(c->vertex( (i+1)%4 )));
    l.push_back(Delaunay::Simplex(c->vertex( (i+2)%4 )));
    l.push_back(Delaunay::Simplex(c->vertex( (i+3)%4 )));
    l.push_back(Delaunay::Simplex(Delaunay::Edge(c, (i+1)%4, (i+2)%4)));
    l.push_back(Delaunay::Simplex(Delaunay::Edge(c, (i+2)%4, (i+3)%4)));
    l.push_back(Delaunay::Simplex(Delaunay::Edge(c, (i+3)%4, (i+1)%4)));
    l.push_back(Delaunay::Simplex(f));
    return l;
}

/**
 * @brief DelaunayHelper::D_finite_faces
 * Return the list of the finite 2-simplices in the boundary of the D-cell ch.
 */
std::list<Delaunay::Simplex> DelaunayHelper::D_finite_faces(const Delaunay &m_dela, const Delaunay::Cell_handle ch)
{
    std::list<Delaunay::Simplex> l;
    for (size_t i = 0; i < 4; i++) {
        const Delaunay::Facet f_i = Delaunay::Facet(ch, i);
        if (!m_dela.is_infinite(f_i))
            l.push_back(Delaunay::Simplex(f_i));
    }
    return l;
}

/**
 * @brief DelaunayHelper::D_finite_sub_faces
 * Return the list of the finite 2-simplices in the boundary of the D-cell ch.
 */
std::list<Delaunay::Simplex> DelaunayHelper::D_finite_sub_faces(const Delaunay &m_dela, const Delaunay::Cell_handle ch)
{
    std::list<Delaunay::Simplex> l;
    std::list<Delaunay::Edge> edges;
    for (size_t i = 0; i < 4; i++) {
        const Delaunay::Vertex_handle v_i = ch->vertex(i);
        if (!m_dela.is_infinite(v_i))
            l.push_back(Delaunay::Simplex(v_i));
    }
    edges.push_back(Delaunay::Edge(ch, 2, 3));
    edges.push_back(Delaunay::Edge(ch, 1, 3));
    edges.push_back(Delaunay::Edge(ch, 0, 3));
    edges.push_back(Delaunay::Edge(ch, 1, 2));
    edges.push_back(Delaunay::Edge(ch, 0, 2));
    edges.push_back(Delaunay::Edge(ch, 0, 1));
    for (Delaunay::Edge e_i: edges) {
        if (!m_dela.is_infinite(e_i))
            l.push_back(Delaunay::Simplex(e_i));
    }
    for (size_t i = 0; i < 4; i++) {
        const Delaunay::Facet f_i = Delaunay::Facet(ch, i);
        if (!m_dela.is_infinite(f_i))
            l.push_back(Delaunay::Simplex(f_i));
    }
    if (!m_dela.is_infinite(ch))
        l.push_back(Delaunay::Simplex(ch));
    return l;
}
//******************************************************************************

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

double DelaunayHelper::max_edge_sq(const Delaunay &m_dela, const Delaunay::Facet f)
{
    Delaunay::Cell_handle c = f.first;
    int i = f.second;
    const Segment s1 = m_dela.segment(Delaunay::Edge(c, (i+1)%4, (i+2)%4));
    const Segment s2 = m_dela.segment(Delaunay::Edge(c, (i+2)%4, (i+3)%4));
    const Segment s3 = m_dela.segment(Delaunay::Edge(c, (i+3)%4, (i+1)%4));
    return max3(s1.squared_length(), s2.squared_length(), s3.squared_length());
}

bool DelaunayHelper::is_infinite(
    const Delaunay &m_dela,
    const Delaunay::Simplex s)
{
    return
    (s.dimension() == 0 && m_dela.is_infinite(Delaunay::Vertex_handle(s))) ||
    (s.dimension() == 1 && m_dela.is_infinite(Delaunay::Edge(s))) ||
    (s.dimension() == 2 && m_dela.is_infinite(Delaunay::Facet(s))) ||
    (s.dimension() == 3 && m_dela.is_infinite(Delaunay::Cell_handle(s)));
}

/**
 * @brief DelaunayHelper::get_critical_info
 * Return a CriticalInfo of the simplex.
 *
 * A Delaunay simplex is critical if it intersects its Voronoi dual. The
 * potentially intersection is called the critical point.
 *
 * Precisely, the function returns a CriticalInfo(c,r,p) where c is the
 * CriticalType of s. if c is CriticalType::Critical, then p is the critical
 * point of the cell and r is its distance to border value.
 */
CriticalInfo DelaunayHelper::get_critical_info(
    const Delaunay &m_dela,
    const Delaunay::Simplex s)
{
    CriticalInfo crit(CriticalType::NonCritical);
    if (!DelaunayHelper::is_infinite(m_dela,s))
    {
        if (s.dimension() == 1) // 2-hole : Voronoi facet
        {
            const Delaunay::Edge e(s);
            if (m_dela.is_Gabriel(e))
            {
                //std::clog << "2-hole critical" << std::endl;
                const Delaunay::Point p0 = m_dela.point(e.first->vertex(e.second));
                const Delaunay::Point p1 = m_dela.point(e.first->vertex(e.third ));
                crit.assign(
                    CriticalType::Critical,
                    0.5*CGAL::sqrt(CGAL::squared_distance(p1, p0)),
                    p0 + 0.5*(p1-p0)
                );
            }
        }
        else if (s.dimension() == 2)  // 1-hole : Voronoi edge
        {
            // A delaunay triangle intersects its Voronoi edge iff it is
            // Gabriel and acute. Hence we need to check those two properties.
            const Delaunay::Facet f(s);
            if (m_dela.is_Gabriel(f))
            {
                Delaunay::Cell_handle ch = f.first;
                int i = f.second;
                const Delaunay::Point a = m_dela.point(ch->vertex( (i+1)%4 ));
                const Delaunay::Point b = m_dela.point(ch->vertex( (i+2)%4 ));
                const Delaunay::Point c = m_dela.point(ch->vertex( (i+3)%4 ));
                const double acsl = (c-a).squared_length();
                const double absl = (b-a).squared_length();
                const double bcsl = (c-b).squared_length();
                const double maxsl = max3(acsl,absl,bcsl);

                // this check if the triangle is acute (using pythagore) :
                if (2*maxsl <= acsl+absl+bcsl)
                {
                    const Delaunay::Point p = CGAL::circumcenter(a,b,c);
                    crit.assign(
                        CriticalType::Critical,
                        CGAL::sqrt((a-p).squared_length()),
                        p
                    );
                }
            }
        }
        else if (s.dimension() == 3) // 0-hole : Voronoi vertex
        {
            // maybe not critical, we need to check if the V-vertex is inside
            // the D-tetrahedra. To realise that we can use cgal or compute the
            // baricentric coordinate of the V-vertex.
            const Delaunay::Cell_handle ch(s);
            const Delaunay::Point dual_point = m_dela.dual(ch);
            if (m_dela.tetrahedron(ch).has_on_bounded_side(dual_point))
            {
                const Delaunay::Point p0 = m_dela.point(ch->vertex(0));
                crit.assign(
                    CriticalType::Critical,
                    CGAL::sqrt(CGAL::squared_distance(dual_point, p0)),
                    dual_point
                );
            }
        }
        else if (s.dimension() == 0) // 3-hole : Voronoi cell
        {
            const Delaunay::Vertex_handle p(s);
            crit.assign(
                CriticalType::Critical,
                0.0,
                m_dela.point(p)
            );
        }
        else
        {
            std::clog << "Error in DelaunayHelper::get_critical_info : s.dimension = " << s.dimension() << std::endl;
        }
    }
    return crit;
}

std::list<Delaunay::Simplex> DelaunayHelper::get_incident_simplices(
    const Delaunay &m_dela,
    const Delaunay::Vertex_handle v)
{
    std::list<Delaunay::Simplex> l;
    std::list<Delaunay::Cell_handle> lc;
    m_dela.incident_cells(v,std::back_inserter(lc));
    std::list<Delaunay::Facet> lf;
    m_dela.incident_facets(v,std::back_inserter(lf));
    std::list<Delaunay::Edge> le;
    m_dela.incident_edges(v,std::back_inserter(le));
    for (Delaunay::Cell_handle c : lc)
        l.push_back(Delaunay::Simplex(c));
    for (Delaunay::Facet f : lf)
        l.push_back(Delaunay::Simplex(f));
    for (Delaunay::Edge e : le)
        l.push_back(Delaunay::Simplex(e));
    return l;
}


void DelaunayHelper::get_finite_neighbouring_cells(
    const Delaunay &m_dela,
    std::vector<Delaunay::Cell_handle> &cells,
    const Point p)
{
    Delaunay::Locate_type lt;
    int li;
    int lj;
    Delaunay::Cell_handle ch = m_dela.locate(p, lt, li, lj) ;
    switch (lt) {
        case Delaunay::Locate_type::CELL:
            cells.push_back(ch);
            std::clog << "O";
            break;
        case Delaunay::Locate_type::FACET:
            get_finite_incident_cells(m_dela, cells, Delaunay::Facet(ch, li));
            std::clog << "#";
            break;
        case Delaunay::Locate_type::EDGE:
            get_finite_incident_cells(m_dela, cells, Delaunay::Edge(ch, li, lj));
            std::clog << "_";
            break;
        case Delaunay::Locate_type::VERTEX:
            get_finite_incident_cells(m_dela, cells, ch->vertex(li));
            std::clog << ".";
            break;
        case Delaunay::Locate_type::OUTSIDE_CONVEX_HULL:
            for (int i = 0 ; i < 4 ; i++)
            {
                if (!m_dela.is_infinite(ch,i))
                {
                    cells.push_back(
                        m_dela.mirror_facet(std::make_pair(ch,i)).first
                    );// should not be infinite if our surface is manifold
                }
            }
            std::clog << "X";
            break;
        default:
            std::clog << "Warning : query outside affine hull" << std::endl;
            break;
    }
}

void DelaunayHelper::get_finite_nearest_cells(
    const Delaunay &m_dela,
    std::vector<Delaunay::Cell_handle> &cells,
    const Point p)
{
    Delaunay::Vertex_handle v = m_dela.nearest_vertex(p);
    get_finite_incident_cells(m_dela, cells, v);
    //std::clog << cells.size() << " ";
}


void DelaunayHelper::get_finite_incident_cells(
    const Delaunay &m_dela,
    std::vector<Delaunay::Cell_handle> &cells,
    const Delaunay::Vertex_handle v)
{
    m_dela.finite_incident_cells(v,std::back_inserter(cells));
}


void DelaunayHelper::get_finite_incident_cells(
    const Delaunay &m_dela,
    std::vector<Delaunay::Cell_handle> &cells,
    const Delaunay::Facet f)
{
    if (!m_dela.is_infinite(f.first))
        cells.push_back(f.first);
    Delaunay::Cell_handle ch_mirror =
        m_dela.mirror_facet(f).first;
    if (!m_dela.is_infinite(ch_mirror))
        cells.push_back(ch_mirror);
}


void DelaunayHelper::get_finite_incident_cells(
    const Delaunay &m_dela,
    std::vector<Delaunay::Cell_handle> &cells,
    const Delaunay::Edge e)
{
    Delaunay::Cell_circulator circ = m_dela.incident_cells(e), past_end(circ);
    do
    {
        Delaunay::Cell_handle c = circ;
        if (!m_dela.is_infinite(c))
            cells.push_back(c);
    }while (++circ != past_end);
}
