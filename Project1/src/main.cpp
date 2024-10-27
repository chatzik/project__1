#include <iostream>
#include <string>
#include <fstream>
#include <queue>
#include <vector>
#include "json.hpp"
#include "triangulation.h"
#include <CGAL/Point_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using json = nlohmann::json;
using namespace std;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K; // Ορισμός του kernel
typedef K::Point_2 Point; 

void loadDataFromJSON(const string &filename, vector<int> &points_x, vector<int> &points_y, vector<int> &region_boundary, vector<pair<int, int>> &additional_constraints,string &instance_uid)
{
    // ανοιγμα JSON
    ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        cerr << "Error: Could not open the file " << filename <<    endl;
        return;
    }

    // φωρτωση δεδομένων
    json j;
    inputFile >> j;

    // αναγνωση των στοιχειων του αρχειου
    points_x = j["points_x"].get<vector<int>>();
    points_y = j["points_y"].get<vector<int>>();
    region_boundary = j["region_boundary"].get<vector<int>>();
    instance_uid =j["instance_uid"].get<string>();
    // cout << "Όνομα που διαβάστηκε: " << instance_uid << endl;
    // Ανάγνωση των constraints
    additional_constraints = j["additional_constraints"].get<vector<pair<int, int>>>();
}
void exportCompletionMessage(string instance_uid) {
    // δημιουργία ενος JSON 
    cout << instance_uid<<endl;
    json outputData;
    outputData["content_type"] = "CG_SHOP_2025_Solution";
    outputData["instance_uid"] = instance_uid;
    outputData["steiner_points_x"] ="sssss";
    outputData["steiner_points_y"] ="sssss";
    outputData["edges:"] ="sssss";
    

    // εγγραφη JSON αρχειου αφου το ανοιξουμε
    ofstream file("output.json");
    if (file.is_open()) {
        file << outputData.dump(6); 
        file.close();
        cout << "Το μήνυμα 'complete' αποθηκεύτηκε στο 'output.json'." <<   endl;
        }
         else {
        cerr << "Σφάλμα: Αδυναμία ανοίγματος του αρχείου για εγγραφή." <<   endl;
        }
}


int main()
{
    queue<Point> steiner_points_queue;
    // δεδομενα  από το JSON αρχεοο
    string instance_uid;
    vector<int> points_x;
    vector<int> points_y;
    vector<int> region_boundary;
    vector<pair<int, int>> additional_constraints;

    //  φωρτωση δεδομενων με την συναρτηση
    loadDataFromJSON("data.json", points_x, points_y, region_boundary, additional_constraints,instance_uid);

    triangulate(points_x, points_y, region_boundary, additional_constraints);
   // cout << "Όνομα που διαβάστηκε: " << instance_uid << endl;
    //std::cout << "Number of Steiner points added: " << steiner_points_queue.size() << std::endl;

    exportCompletionMessage(instance_uid);

    while (!steiner_points_queue.empty())
    {
        Point steiner = steiner_points_queue.front();
        steiner_points_queue.pop();

        //εκτυπωση των steiner σημειων σαν ζυγαρια χ , ψ
        cout << "Steiner point: (" << steiner.x() << ", " << steiner.y() << ")" << std::endl;
    }

    return 0;
}
