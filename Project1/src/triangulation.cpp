#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include "triangulation.h"
#include <CGAL/draw_triangulation_2.h>
#include <cmath> // Χρήση της βιβλιοθήκης math για το εσωτερικό γινόμενο

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

// υπολογισμός εσωτερικού γινομένου δύο διανυσμάτων (χρησιμοποιεί την βιβλιοθήκη math)
double dotProduct(const Point &a, const Point &b, const Point &c)
{
    double ab_x = b.x() - a.x();
    double ab_y = b.y() - a.y();
    double ac_x = c.x() - a.x();
    double ac_y = c.y() - a.y();
    return (ab_x * ac_x + ab_y * ac_y);
}

// Η has_Obtuse_Angle κάνει έλεγχο αν το τρίγωνο έχει έστω και μία αμβλεία γωνία
bool has_Obtuse_Angle(const Point &a, const Point &b, const Point &c)
{
    // σε αυτό σημείο κάνει τον έλεγχο με την βοήθεια εσωτερικού γινομένου για το αν υπάρχει αμβλεία γωνία
    double dotAB_AC = dotProduct(a, b, c); // ελέχγει την γωνία που δημιουργούν τα ΑΒ-ΑΓ
    double dotBA_BC = dotProduct(b, a, c); // ελέχγει την γωνία που δημιουργούν τα ΑΒ-ΒΓ
    double dotCA_CB = dotProduct(c, a, b); // ελέχγει την γωνία που δημιουργούν τα ΑΓ-ΒΓ

    return (dotAB_AC < 0 || dotBA_BC < 0 || dotCA_CB < 0);
}

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

    int count_Obtuse_Angles = 0; // μετρητής για το πόσες αμβλείες γωνίες υπάρχουν

    // Ελέγχουμε κάθε τρίγωνο στην τριγωνοποίηση
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point a = fit->vertex(0)->point();
        Point b = fit->vertex(1)->point();
        Point c = fit->vertex(2)->point();

        // Έλεγχος για αμβλεία γωνία
        if (has_Obtuse_Angle(a, b, c))
        {
            count_Obtuse_Angles++;
        }
    }

    cout << "Number of obtuse angles in the triangulation: " << count_Obtuse_Angles << endl;

    // Εμφάνιση της τριγωνοποίησης
    draw(cdt);
}
