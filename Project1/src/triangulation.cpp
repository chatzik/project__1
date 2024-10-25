#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Delaunay_mesher_2.h>
#include "triangulation.h"
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <unordered_map>
#include <cmath>

//////////////////////////////////////////////////////////

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Face_handle Face_handle;
using namespace std;

///////////////////////////////////////////////////////////

double calculate_angle(const Point &A, const Point &B, const Point &C)
{
    double a2 = squared_distance(B, C); // απέναντι πλευρά από το A
    double b2 = squared_distance(A, C); // απέναντι πλευρά από το B
    double c2 = squared_distance(A, B); // απέναντι πλευρά από το C

    // νόμος των συνημιτόνων: cos(γωνία) = (b² + c² - a²) / (2*ρίζα(b)*ρίζα(c))
    double cos_A = (b2 + c2 - a2) / (2 * sqrt(b2) * sqrt(c2));

    // μεατροπή της γωνίας σε μοίρες
    return acos(cos_A) * 180.0 / M_PI;
}

// η collinear ελέγχει αν τα 3 σημεία που εισάγαμε είναι συγγραμικά, δηλαδή δεν σχηματίζουν τρίγωνο
// Η has_Obtuse_Angle κάνει έλεγχο αν το τρίγωνο έχει έστω και μία αμβλεία γωνία
int has_Obtuse_Angle(const Point &a, const Point &b, const Point &c)
{
    if (collinear(a, b, c))
    {
        return -1; // τα σημεία είναι συγγραμμικά, οπότε δεν σχηματίζουν τρίγωνο
    }
    // υπολογισμοσ  των γωνίων του τριγώνου
    double angle_A = calculate_angle(a, b, c);
    double angle_B = calculate_angle(b, a, c);
    double angle_C = calculate_angle(c, a, b);

    // έλεγχος για αμβλεία γωνία
    if (angle_A > 90)
    {
        cout << "Obtuse angle found at A! " << angle_A << endl;
        return 0; // Αμβλεία γωνία στο A
    }
    if (angle_B > 90)
    {
        cout << "Obtuse angle found at B! " << angle_B << endl;
        return 1; // Αμβλεία γωνία στο B
    }
    if (angle_C > 90)
    {
        cout << "Obtuse angle found at C! " << angle_C << endl;
        return 2; // Αμβλεία γωνία στο C
    }

    return -1; // Δεν υπάρχει αμβλεία γωνία
}

Point insert_Steiner(const Point &a, const Point &b) // επιστρέφει το μέσο της απέναντι πλευράς
{
    return Point((a.x() + b.x()) / 2, (a.y() + b.y()) / 2);
}

// επιστρέφει το πλήθος των αμβλείων γωνιών του γράφου
int count_Obtuse_Angles(CDT &cdt)
{
    int count = 0;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point a = fit->vertex(0)->point();
        Point b = fit->vertex(1)->point();
        Point c = fit->vertex(2)->point();
        cout << "Checking triangle with points: (" << a.x() << ", " << a.y() << "), (" << b.x() << ", " << b.y() << "), (" << c.x() << ", " << c.y() << ")" << endl;

        if (int is_Obtuse = has_Obtuse_Angle(a, b, c) != -1)
        {
            count++;
            // cout << "Obtuse angle found at vertex " << is_Obtuse << endl;
        }
    }
    return count;
}

// τριγωνοποίηση
void triangulate(const vector<int> &points_x, const vector<int> &points_y, const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints)
{
    CDT cdt;

    // προσθήκη σημείων και ακμών του γράφου
    vector<Point> points;
    for (size_t i = 0; i < points_x.size(); ++i)
    {
        points.push_back(Point(points_x[i], points_y[i]));
    }
    for (int i = 0; i < points_x.size(); i++)
        cout << points_x[i] << " ";
    for (int i = 0; i < points_y.size(); i++)
        cout << points_y[i] << " ";

    // προσθήκη των ακμών του κυρτού περιβλήματος για όριο
    for (size_t i = 0; i < region_boundary.size(); ++i)
    {
        int next = (i + 1) % region_boundary.size();
        cdt.insert_constraint(points[region_boundary[i]], points[region_boundary[next]]);
    }

    // προσθήκη πρόσθετων constraints
    for (const auto &constraint : additional_constraints)
    {
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }

    int original_graph_count = count_Obtuse_Angles(cdt); // μετρητής για το πόσες αμβλείες γωνίες υπάρχουν
    if (original_graph_count > 0)
    {
        bool found_steiner_point = true;
        // ελέγχουμε κάθε τρίγωνο στην τριγωνοποίηση
        while (found_steiner_point)
        {
            found_steiner_point = false; // Αν δεν βρεθεί, το loop θα σταματήσει
            int original_graph_count = count_Obtuse_Angles(cdt);
            // Κάνουμε iterate στα finite faces της τριγωνοποίησης
            for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
            {
                Point a = fit->vertex(0)->point();
                Point b = fit->vertex(1)->point();
                Point c = fit->vertex(2)->point();

                int is_Obtuse = has_Obtuse_Angle(a, b, c);
                if (is_Obtuse != -1) // δηλαδή η γωνία ΔΕΝ είναι οξεία ή κάθετη -> αμβλεία
                {
                    Point steiner;
                    if (is_Obtuse == 0)
                        steiner = insert_Steiner(b, c); // απέναντι πλευρά από το A
                    else if (is_Obtuse == 1)
                        steiner = insert_Steiner(a, c); // απέναντι πλευρά από το σημείο B
                    else                                // δλδ, is_Obtuse = 2
                        steiner = insert_Steiner(a, b); // απέναντι πλευρά από το C

                    // Εισαγωγή του Steiner σημείου και τερματισμός του loop
                    CDT temp_cdt = cdt;
                    temp_cdt.insert(steiner);

                    // Υπολογισμός των αμβλείων γωνιών μετά την εισαγωγή
                    int new_obtuse_count = count_Obtuse_Angles(temp_cdt);
                    if (new_obtuse_count < original_graph_count)
                    {
                        cdt.insert(steiner);
                        found_steiner_point = true;
                        cout << "Inserted Steiner point at: (" << steiner.x() << ", " << steiner.y() << ")" << endl;
                        break;
                    }
                }
            }
        }
        original_graph_count = count_Obtuse_Angles(cdt);
    }
    cout << "Number of obtuse angles in the triangulation: " << original_graph_count << endl;

    // Εμφάνιση της τριγωνοποίησης
    draw(cdt);
}
