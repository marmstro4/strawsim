#include <iostream>
#include <fstream>
#include <string>
#include "TH1D.h"
#include "TCanvas.h"

void vertex() {
    // Open the CSV file
    std::ifstream inputFile("build/vertex.csv");

    // Check if the file was opened successfully
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file output.csv!" << std::endl;
        return;
    }

    // Create a histogram with 1000 bins, range [-10, 10]
    TH1D *hist = new TH1D("hist", "Histogram of Values from output.csv;Value;Entries", 1000, -10, 10);

    // Read numbers from the file and fill the histogram
    std::string line;
    while (std::getline(inputFile, line)) {
        try {
            // Convert the line to a double
            double value = std::stod(line);

            // Fill the histogram with the value
            hist->Fill(value);
        } catch (const std::exception &e) {
            // Handle conversion errors
            std::cerr << "Skipping invalid line: " << line << std::endl;
        }
    }

    // Close the file
    inputFile.close();

    // Draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Histogram", 800, 600);
    hist->Draw();

    // Save the histogram as an image (optional)
    canvas->SaveAs("histogram.png");

    // Note: The ROOT application needs to stay open for the canvas to remain visible
}
