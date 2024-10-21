#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include "triangulation.h"
#include <CGAL/draw_triangulation_2.h>

//////////////////////////////////////////////////////////

// Ορισμός του Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef CDT::Point Point;
using namespace std;

// Υλοποίηση της τριγωνοποίησης
void triangulate(const vector<int> &points_x, const vector<int> &points_y, const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints)
{
    CDT cdt;

    // Προσθήκη σημείων και ακμών
    vector<Point> points;
    for (size_t i = 0; i < points_x.size(); ++i)
    {
        points.push_back(Point(points_x[i], points_y[i]));
    }
    for (int i = 0; i < points_x.size(); i++)
        cout << points_x[i] << " ";
    for (int i = 0; i < points_y.size(); i++)
        cout << points_y[i] << " ";

    // Προσθήκη των ακμών του κυρτού περιβλήματος
    for (size_t i = 0; i < region_boundary.size(); ++i)
    {
        int next = (i + 1) % region_boundary.size();
        cdt.insert_constraint(points[region_boundary[i]], points[region_boundary[next]]);
    }

    // Προσθήκη πρόσθετων constraints
    for (const auto &constraint : additional_constraints)
    {
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }

    // Εμφάνιση της τριγωνοποίησης
    draw(cdt);
}
