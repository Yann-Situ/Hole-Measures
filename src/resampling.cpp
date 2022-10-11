#include "resampling.h"
#include "delaunay_helper.h"
#include <iostream>
#include <limits>

void Resampling::add_circumcenters(Delaunay &m_dela, Polyhedron m_poly)
{
    std::clog << "Add Circumcenters " << std::endl;
    for (Polyhedron::Facet_iterator  f = m_poly.facets_begin(); f != m_poly.facets_end(); ++f) {
        m_dela.insert(circumcenter(f));
    }
    std::clog << "nb vertices: " << m_dela.number_of_vertices() << std::endl;
    std::clog << "nb edges : " << m_dela.number_of_finite_edges() << std::endl;
    std::clog << "nb facets : " << m_dela.number_of_finite_facets() << std::endl;
    std::clog << "nb cells : " << m_dela.number_of_finite_cells() << std::endl;
}

Point  Resampling::circumcenter(Polyhedron::Facet_iterator f)
{
    Polyhedron::Facet::Halfedge_handle h = f->halfedge();
    return CGAL::circumcenter( h->vertex()->point(),
                 h->next()->vertex()->point(),
         h->next()->next()->vertex()->point());
}

/*############################################################################*/

void Resampling::add_barycenters(Delaunay &m_dela, Polyhedron m_poly)
{
    std::clog << "Add Barycenters " << std::endl;
    for (Polyhedron::Facet_iterator  f = m_poly.facets_begin(); f != m_poly.facets_end(); ++f) {
        //if (m_poly.is_triangle(f->halfedge())){
        m_dela.insert(Resampling::barycenter(f));
        //}
    }
    std::clog << "nb vertices: " << m_dela.number_of_vertices() << std::endl;
    std::clog << "nb edges : " << m_dela.number_of_finite_edges() << std::endl;
    std::clog << "nb facets : " << m_dela.number_of_finite_facets() << std::endl;
    std::clog << "nb cells : " << m_dela.number_of_finite_cells() << std::endl;
}

Point  Resampling::barycenter(Polyhedron::Facet_iterator f)
{
    Polyhedron::Facet::Halfedge_handle h = f->halfedge();
    return CGAL::centroid( h->vertex()->point(),
                 h->next()->vertex()->point(),
         h->next()->next()->vertex()->point());
}

/*############################################################################*/

void Resampling::write_off(Polyhedron m_poly)
{
    CGAL::set_ascii_mode( std::cout);
    std::cout << "OFF" << std::endl << m_poly.size_of_vertices() << ' '
    << m_poly.size_of_facets() << " 0" << std::endl;
    std::copy( m_poly.points_begin(), m_poly.points_end(),
    std::ostream_iterator<Point>( std::cout, "\n"));
    for (  Polyhedron::Facet_iterator i = m_poly.facets_begin(); i != m_poly.facets_end(); ++i) {
        Polyhedron::Halfedge_around_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        std::cout << CGAL::circulator_size(j) << ' ';
        do {
            std::cout << ' ' << std::distance(m_poly.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
        std::cout << std::endl;
    }
}

/*############################################################################*/

/**
 * @brief FiltrationVoronoi::criteria_1
 * Check if the height of the tetrahedra of f is above
 * sqrt(parameter)*max_{e edge of f}(||e||)
 */
bool criteria_1(const Delaunay m_dela, const Delaunay::Facet f, const double parameter)
{
    Delaunay::Cell_handle c1 = f.first;
    Delaunay::Cell_handle c2 = m_dela.mirror_facet(f).first;
    int i = f.second;
    const double thresh = parameter * DelaunayHelper::max_edge_sq(m_dela, f);

    const Plane plane = m_dela.triangle(f).supporting_plane();
    bool cond1 = m_dela.is_infinite(c1);
    bool cond2 = m_dela.is_infinite(c2);
    if (!cond1)
    {
        const Point p1 = m_dela.point(c1, i);
        const double p1_height = CGAL::squared_distance(plane, p1);
        cond1 = (p1_height > thresh);

    }
    if (!cond2)
    {
        const Point p2 = m_dela.point(c2, m_dela.mirror_facet(f).second);
        const double p2_height = CGAL::squared_distance(plane, p2);
        cond2 = (p2_height > thresh);
    }
    return cond1 && cond2;
}
double criteria_2(const Delaunay m_dela, const Delaunay::Facet f)
{
    Delaunay::Cell_handle c1 = f.first;
    Delaunay::Cell_handle c2 = m_dela.mirror_facet(f).first;
    int i = f.second;
    const double thresh = DelaunayHelper::max_edge_sq(m_dela, f);

    const Plane plane = m_dela.triangle(f).supporting_plane();
    bool cond1 = m_dela.is_infinite(c1); double v1 = 0.0;
    bool cond2 = m_dela.is_infinite(c2); double v2 = 0.0;
    if (!cond1)
    {
        const Point p1 = m_dela.point(c1, i);
        const double p1_height = CGAL::squared_distance(plane, p1);
        v1 = thresh/p1_height;

    }
    if (!cond2)
    {
        const Point p2 = m_dela.point(c2, m_dela.mirror_facet(f).second);
        const double p2_height = CGAL::squared_distance(plane, p2);
        v2 = thresh/p2_height;
    }
    return (v1 > v2) ? v1 : v2; // goal : criteria_2 < parameter
}

void Resampling::add_from_criteria(Delaunay &m_dela, SoT inside)
{
    std::clog << "Add from criteria" << std::endl;
    int max_add = 100;
    int nb_add = 0;
    int i = 0;
    std::vector<Point> to_add;
    std::deque<Triangle> b_facets;
    for (Delaunay::Facet f : m_dela.finite_facets())
        b_facets.push_back(m_dela.triangle(f));

    while (nb_add < max_add) {
        if (b_facets.empty())
        {
            if (to_add.empty())
                break;
            for (int j = 0; j<to_add.size(); j++)
            {
                const Delaunay::Vertex_handle new_vertex =
                      m_dela.insert(to_add[j]);
                std::list<Delaunay::Facet> new_facets;
                m_dela.finite_incident_facets(new_vertex,std::back_inserter(new_facets));
                for (Delaunay::Facet new_facet : new_facets)
                    b_facets.push_back(m_dela.triangle(new_facet));
            }
            nb_add+=to_add.size();
            std::clog << "Add " << to_add.size() << " points and " << b_facets.size() << " facets" << std::endl;
            to_add.clear();
        }
        Triangle t = b_facets.front();
        b_facets.pop_front();
        std::clog << i << " ";
        i++;

        Delaunay::Cell_handle c1;
        Delaunay::Vertex_handle v1;
        Delaunay::Vertex_handle v2;
        Delaunay::Vertex_handle v3;
        m_dela.is_vertex(t[0], v1);
        m_dela.is_vertex(t[1], v2);
        m_dela.is_vertex(t[2], v3);
        int i1,i2,i3;
        if (!m_dela.is_facet(v1,v2,v3, c1, i1,i2,i3))
        {
            std::clog << "facet not valid" << std::endl;
            continue;
        }
        // now c1, i1, i2, i3 should be updated (i_j rpz v_j) and c1 the cell
        int ind_f = 0+1+2+3 - (i1+i2+i3);
        Delaunay::Facet f = std::make_pair(c1, ind_f);
        Delaunay::Cell_handle c2 = m_dela.mirror_facet(f).first;

        // now check if it is on boundary
        bool is_inside_1 = !m_dela.is_infinite(c1) &&
        (inside(CGAL::centroid(m_dela.tetrahedron(c1)))== CGAL::ON_BOUNDED_SIDE);
        bool is_inside_2 = !m_dela.is_infinite(c2) &&
        (inside(CGAL::centroid(m_dela.tetrahedron(c2)))== CGAL::ON_BOUNDED_SIDE);

        if ((is_inside_1 != is_inside_2) && !criteria_1(m_dela, f, 0.01))
        {
            // it's on the boundary and doesn't verify the criteria
            to_add.push_back(CGAL::centroid(m_dela.triangle(f)));
            std::clog << "on boundary and !criteria" << std::endl;
        }
        else if (is_inside_1 != is_inside_2)
            std::clog << "on boundary" << std::endl;
        else
            std::clog << "not on boundary" << std::endl;
    }

    std::clog << "Resampling added " << nb_add << " vertices." << std::endl;
    std::clog << "nb vertices: " << m_dela.number_of_vertices() << std::endl;
    std::clog << "nb edges : " << m_dela.number_of_finite_edges() << std::endl;
    std::clog << "nb facets : " << m_dela.number_of_finite_facets() << std::endl;
    std::clog << "nb cells : " << m_dela.number_of_finite_cells() << std::endl;
}

void Resampling::add_from_criteria_cell(Delaunay &m_dela, SoT inside)
{
    std::clog << "Add from criteria" << std::endl;
    const int max_add = 100;
    int nb_add = 0;
    int i = 0;
    std::vector<Point> to_add;
    std::set<Delaunay::Cell_handle> cells;
    for (Delaunay::Cell_handle c : m_dela.finite_cell_handles())
        cells.insert(c);

    std::set<Delaunay::Cell_handle>::iterator it_cells = cells.begin();
    while (nb_add < max_add) {
        if (it_cells == cells.end())
        {
            cells.clear();
            if (to_add.empty())
                break;
            for (int j = 0; j<to_add.size(); j++)
            {
                const Delaunay::Vertex_handle new_vertex = m_dela.insert(to_add[j]);
                std::list<Delaunay::Cell_handle> new_cells;
                m_dela.finite_incident_cells(new_vertex,std::back_inserter(new_cells));
                cells.insert(new_cells.begin(), new_cells.end());
            }

            nb_add+=to_add.size();
            std::clog << "Add " << to_add.size() << " points and " << cells.size() << " facets" << std::endl;
            to_add.clear();
            it_cells = cells.begin();
        }
        std::clog << i << " ";
        i++;
        const Delaunay::Cell_handle c = *it_cells;
        it_cells++;
        if (!m_dela.is_cell(c))
        {
            std::clog << "cell not valid" << std::endl;
            continue;
        }
        const bool is_inside = (inside(CGAL::centroid(m_dela.tetrahedron(c)))
            == CGAL::ON_BOUNDED_SIDE);

        std::vector<Delaunay::Facet> b_facets;
        // now check if it has boundary facets
        for (int ind_f = 0; ind_f < 4; ind_f ++)
        {
            const Delaunay::Facet f = std::make_pair(c, ind_f);
            const Delaunay::Cell_handle c2 = m_dela.mirror_facet(f).first;
            const bool is_inside_2 = !m_dela.is_infinite(c2) &&
            (inside(CGAL::centroid(m_dela.tetrahedron(c2)))== CGAL::ON_BOUNDED_SIDE);
            if (is_inside != is_inside_2)
            {
                b_facets.push_back(f);
            }
        }

        if (b_facets.size() == 0){
            std::clog << "no boundary facets !" << std::endl;
            continue;
        }
        else if (b_facets.size() > 1){
            std::clog << "several boundary facets..." << std::endl;
            continue;
        }
        // it has only a single boundary facet
        const Delaunay::Facet f = b_facets.back();

        if (!criteria_1(m_dela, f, 0.01))
        {
            // it's on the boundary and doesn't verify the criteria
            to_add.push_back(CGAL::centroid(m_dela.triangle(f)));
            std::clog << "1 boundary facet and !criteria" << std::endl;
        }
        else
            std::clog << "1 boundary facet" << std::endl;
    }

    std::clog << "Resampling added " << nb_add << " vertices." << std::endl;
    std::clog << "nb vertices: " << m_dela.number_of_vertices() << std::endl;
    std::clog << "nb edges : " << m_dela.number_of_finite_edges() << std::endl;
    std::clog << "nb facets : " << m_dela.number_of_finite_facets() << std::endl;
    std::clog << "nb cells : " << m_dela.number_of_finite_cells() << std::endl;
}

/*############################################################################*/
MedialFeatureField::MedialFeatureField(const Delaunay &m_dela, double parameter)
: m_dela(m_dela), parameter(parameter)
{
    // for (Delaunay::Vertex_handle v : m_dela.finite_vertex_handles())
    // {
    //     v->info() = -1.0;
    // }
    //
    // for (Delaunay::Cell_handle ch : m_dela.finite_cell_handles())
    // {
    //     const Delaunay::Point dual_point = m_dela.dual(ch);
    //     const Delaunay::Point point = m_dela.point(ch, 0);
    //     FT r = CGAL::squared_distance(dual_point, point);
    //     for (int i = 0; i < 4; i++)
    //     {
    //         Delaunay::Vertex_handle v = ch->vertex(i);
    //         if (v->info() < 0.0)
    //             v->info() = r;
    //         else
    //             v->info() = (r < v->info()) ? r : v->info();
    //         // auto search = sqr_circumradius.find(v);
    //         // if (search == sqr_circumradius.end()) {
    //         //     sqr_circumradius[v]=r;
    //         // } else {
    //         //     // take minimum
    //         //     search->second = (r < search->second) ? r : search->second;
    //         // }
    //     }
    // }
    std::clog << "create medial feature field" << std::endl;
}

FT MedialFeatureField::operator()(const Point_3& p, const int dimension, const Index& index) const
{
    // Delaunay::Vertex_handle v = m_dela.nearest_vertex(p);
    // if (v->info() < 0.0)
    //     std::clog << v->info() << " ";
    // else
    //     std::clog << v->info() << " ";
    // FT r = v->info();
    // return parameter * CGAL::sqrt(v->info());
    //
    std::vector<Delaunay::Cell_handle> cells;
    DelaunayHelper::get_finite_nearest_cells(m_dela, cells, p);
    //DelaunayHelper::get_finite_neighbouring_cells(m_dela, cells, p);

    FT min_dist = std::numeric_limits< double >::max();// to modify with something as max_double
    for (Delaunay::Cell_handle cell : cells)
    {
        const Delaunay::Point dual_point = m_dela.dual(cell);
        const Delaunay::Point point = m_dela.point(cell, 0);
        const FT r = CGAL::squared_distance(dual_point, point);
        //const FT r = sqr_circumradius.at(cell);
        if (r < min_dist)
            min_dist = r;
    }
    return parameter * CGAL::sqrt(min_dist);
}
