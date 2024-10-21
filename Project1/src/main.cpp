#include <iostream>
#include <fstream>
#include <vector>
#include "json.hpp"          // Χρήση της βιβλιοθήκης JSON
#include "triangulation.h"   // Επικεφαλίδα για την τριγωνοποίηση

using json = nlohmann::json;
using namespace std;

// Συνάρτηση για φόρτωση των δεδομένων από το JSON αρχείο
void loadDataFromJSON(const string& filename,vector<int>& points_x, vector<int>& points_y, vector<int>& region_boundary, vector<pair<int, int>>& additional_constraints) {
    // Άνοιγμα αρχείου JSON
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open the file " << filename << std::endl;
        return;
    }

    // Φόρτωση JSON δεδομένων
    json j;
    inputFile >> j;

    // Ανάγνωση των arrays από το JSON
    points_x = j["points_x"].get<vector<int>>();
    points_y = j["points_y"].get<vector<int>>();
    region_boundary = j["region_boundary"].get<vector<int>>();

    // Ανάγνωση των constraints
    additional_constraints = j["additional_constraints"].get<vector<pair<int, int>>>();
}

int main() {
    // Δεδομένα που θα φορτωθούν από το JSON αρχείο
    vector<int> points_x;
    vector<int> points_y;
    vector<int> region_boundary;
    vector<pair<int, int>> additional_constraints;

    // Κάλεσμα της συνάρτησης για φόρτωση δεδομένων
    loadDataFromJSON("data.json", points_x, points_y, region_boundary, additional_constraints);

    // Εκτέλεση τριγωνοποίησης
    triangulate(points_x, points_y, region_boundary, additional_constraints);

    return 0;
}
