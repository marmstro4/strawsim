#include "TH2F.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

void plot_radius_vs_vertex() {
    // Input file name
    const char* filename = "build/check.csv";

    // Vectors to store the data
    std::vector<double> x; // First column
    std::vector<double> y; // Third column

    // Open the input file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Read the file line by line
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string value;
        double col1, col2, col3, col4;

        // Parse the three columns
        if (std::getline(ss, value, ',')) col1 = std::stod(value);
        if (std::getline(ss, value, ',')) col2 = std::stod(value);
        if (std::getline(ss, value, ';')) col3 = std::stod(value);
        if (std::getline(ss, value, ';')) col4 = std::stod(value);

        // Store the first and third columns
        x.push_back(col1);
        y.push_back(col4);
    }
    file.close();

    // Check if data was read
    if (x.empty() || y.empty()) {
        std::cerr << "Error: No data to plot." << std::endl;
        return;
    }

    // Define histogram bins
    double x_min = 1, x_max = 7;  // Range for x-axis
    int x_bins = 7;
    double y_min = -2, y_max = 2;  // Range for x-axis
    int y_bins = 1000;

    // Create a 2D histogram
    //TH2F* hist = new TH2F("hist", "Radius vs Vertex",x_bins, x_min, x_max,y_bins, y_min, y_max);

    TH1F* hist = new TH1F("hist", "Radius vs Vertex",y_bins, y_min, y_max);

    // Fill the histogram with data
    for (size_t i = 0; i < x.size(); ++i) {
        hist->Fill(y[i]);
    }

    // Create a canvas to draw the histogram
    TCanvas* canvas = new TCanvas("canvas", "Radius vs Vertex Histogram", 800, 600);
    hist->Draw("colz"); // Draw the histogram with color mapping

    // Update the canvas
    canvas->Update();
}
