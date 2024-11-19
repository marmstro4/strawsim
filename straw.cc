#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TEllipse.h>
#include <TLine.h>
#include <tuple>
#include <TGraph.h>
#include <vector>
#include <cmath>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <stdexcept>
#include "TVectorD.h"
#include <fstream>
#include <random>
#include <thread>
#include <chrono>
#include<TF1.h>

#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewCell.hh"

using namespace Garfield;
using namespace std;

std::random_device rd;  // Seed
std::mt19937 gen(rd()); // Mersenne Twister RNG
std::uniform_real_distribution<> dis(-1, 1); // Range [0.0, 1.0)

string GarPath = "/home/michaela/Downloads/garfield/garfieldpp/";
string GasPath = "/home/michaela/straw/MichaelGarfield/cell_study/gas/";


double tstep = 0.25; //choose so devisible by rebin
double SigI = 0;
double SigF = 5000;
int rebin = 1;
int SigBins = static_cast<int>(std::round((SigF-SigI)/tstep));
int avg_lim = 10;
double var_tol = 0.01;

TH1D *RecoZ = new TH1D("RecoZ", "RecoZ",10000,-100,100);

//Wire radius [cm]
double rWire = 25e-4;

//Tube radius [cm]
double rTube = 1;

//Tube speration [cm]
double VerticalOffset = 10;

//Voltage
double vWire = 2730;
double vTube = 0;

//Tracks
std::vector<double> trk_x,trk_y,trk_z,trk_t;

struct FitData {
    const std::vector<double>* zCenter;
    const std::vector<double>* yCenter;
    const std::vector<int>* hitI;
    const std::vector<double>* hitR;
};

double slantfunct(double *x, double *params) {
    double xx = x[0];          // x is the independent variable
    double A = params[0];      // Amplitude
    double mu = params[1];     // Mean
    double sigma = params[2];  // Standard deviation
    double alpha = params[3];  // Plateau control parameter
    double beta = params[4];

    // Gaussian numerator
    double gaussian = A * exp(0.5 * pow((xx - mu) / sigma, 2));

    // Denominator to introduce flatness
    double denominator = 1 + exp(-0.5 * pow((xx - mu) / sigma, 2));

    return gaussian / denominator - alpha * (xx - mu) + beta;
}

void SetGas(string GasName, MediumMagboltz* gas) {

    gas->SetPressure(3 * 760.);
    gas->SetTemperature(293.15);
    gas->LoadGasFile(GasPath + GasName);
    gas->LoadIonMobility(GarPath + "Data/IonMobility_Ar+_Ar.txt");

    //ViewMedium mediumView;
    //mediumView.SetMedium(gas);
    //mediumView.PlotElectronVelocity('e');


}

void SetWire(ComponentAnalyticField* cmp, MediumMagboltz* gas) {

  cmp->AddWire(0, 0, 2*rWire, vWire, "s");
  cmp->AddTube(rTube, vTube, 0);
  cmp->SetMedium(gas);

}

void SetSensor(Sensor* sensor, ComponentAnalyticField* cmp) {

  sensor->AddComponent(cmp);
  sensor->AddElectrode(cmp, "s");

  sensor->SetTimeWindow(SigI, tstep, SigBins);

   // Open the file
    std::ifstream infile("mdt_elx_delta.txt");

    // Declare vectors for times and values
    std::vector<double> times;
    std::vector<double> values;

    // Read the file line by line
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double time, value;
        if (!(iss >> time >> value)) continue; // Skip invalid lines

        times.push_back(1.e3 * time); // Convert time to ms
        values.push_back(value);
    }

    infile.close();

    //sensor->SetTransferFunction(times, values);

}

void SignalPlot(std::vector<double> x, std::vector<double> y) {
    TH1D *Signal = new TH1D("Signal", "Signal",SigBins,SigI,SigF);

    for (size_t i=0; i<x.size(); i++) {
        Signal->SetBinContent(x[i],y[i]);
        Signal->Draw("hist");
    }

}

std::pair<double,double> SignalProcessing(std::vector<double> x, std::vector<double> y, std::vector<double> x_basket) {

    //Get maximum charge and time when
    int amp_in = std::distance(y.begin(), max_element(std::begin(y), std::end(y)));
    double amp = y[amp_in];
    int UpTail;

    //Determine smoothness
    for (size_t i=amp_in+2*avg_lim; i<x.size(); i++) {

        double sum = 0;
        double diff_sum = 0;

        //Moving average
        for (size_t j=i-avg_lim; j<i+avg_lim; j++) {
            sum+=y[j];
        }

        double move_avg = sum/(2*avg_lim);

        //Moving variance
        for (size_t j=i-avg_lim; j<i+avg_lim; j++) {
            diff_sum+=(y[j]-move_avg)*(y[j]-move_avg);
        }

        double var = sqrt(diff_sum/(2*avg_lim-1));

        if (var<var_tol) {
            UpTail = i;
            break;
        }

        //cout<<i*tstep<<" "<<y[i]<<" "<<move_avg<<" "<<var<<endl;

    }

    //cout<<UpTail*tstep<<" "<<x_basket.front()*tstep<<endl;
    return std::make_pair(UpTail*tstep-x_basket.front()*tstep,x_basket.back()-x_basket.front()); //[ns]

}


void PlotStraws(const std::vector<double>& zCenter, const std::vector<double>& yCenter, double radius, TCanvas* canvas) {
    // Check if the inputs are valid

    // Create a canvas

    canvas->SetFixedAspectRatio(); // Ensure equal scaling on both axes
    canvas->DrawFrame(-15, 0, 15, 30); // Set the drawing frame (xMin, yMin, xMax, yMax)

    // Draw the circles (straws)
    for (size_t i = 0; i < yCenter.size(); ++i) {
        TEllipse* straw = new TEllipse(zCenter[i], yCenter[i], radius);
        straw->SetFillStyle(0); // No fill
        straw->SetLineColor(kBlack);
        straw->Draw();
    }

    // Update the canvas to render the drawing
    canvas->Update();
}

std::pair<std::vector<double>, std::vector<double>> GetStrawCenters(double strawRadius, double verticalOffset,
                     std::vector<double>& yCenter, std::vector<double>& zCenter) {
    // Calculate the number of tubes needed to cover +/-10 cm
    int nTubes = static_cast<int>(5 / strawRadius) + 1;

    // Generate the first layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
    }

    // Generate the second layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + std::sqrt(3) * strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius + strawRadius);
    }

    // Generate the third layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
    }

    return std::make_pair(yCenter, zCenter);

}

void PlotStrawsAndTrack(const std::vector<double>& zCenter, const std::vector<double>& yCenter,
                        double radius, const std::vector<int>& hitI, const std::vector<double>& hitR,
                        double m, double zOrigin) {
    TCanvas* canvas = new TCanvas("canvas", "Straws and Track", 800, 800);
    canvas->DrawFrame(-15, 0, 15, 30, "Plot and Trajectory");
    canvas->SetGrid();

    for (size_t i = 0; i < yCenter.size(); ++i) {
        TEllipse* strawCircle = new TEllipse(zCenter[i], yCenter[i], radius, radius);
        strawCircle->SetLineColor(kBlack);
        strawCircle->SetFillStyle(0);
        strawCircle->Draw();

        for (size_t j = 0; j < hitI.size(); ++j) {
            if (static_cast<int>(i) == hitI[j]) {
                TEllipse* hitCircle = new TEllipse(zCenter[i], yCenter[i], hitR[j], hitR[j]);
                hitCircle->SetLineColor(kRed);
                hitCircle->SetFillStyle(0);
                hitCircle->Draw();
            }
        }
    }

    double b = -m * zOrigin;
    double z1 = -15;
    double z2 = 15;
    double y1 = m * z1 + b;
    double y2 = m * z2 + b;

    TLine* trackLine = new TLine(z1, y1, z2, y2);
    trackLine->SetLineColor(kRed);
    trackLine->Draw();

    canvas->SaveAs("StrawsAndTrack.png");
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<double>> GetRadiiForTrack(const std::vector<double>& zCenter, const std::vector<double>& yCenter,
                      double radius, double m, double zOrigin,
                      std::vector<int>& hitN, std::vector<int>& hitI, std::vector<double>& hitR) {
    double b = -m * zOrigin;

    for (size_t i = 0; i < yCenter.size(); ++i) {
        double d = std::abs(m * zCenter[i] - yCenter[i] + b) / std::sqrt(m * m + 1);
        if (d <= radius) {
            hitI.push_back(static_cast<int>(i));
            hitR.push_back(d);
            hitN.push_back(i);
        }
    }

    return std::make_tuple(hitN, hitI, hitR);
}

double CalcFOM(const double* params, const FitData& data) {
    double m = params[0];
    double zOrigin = params[1];
    double b = -m * zOrigin;
    double fom = 0;

    for (size_t i = 0; i < data.hitI->size(); ++i) {
        int index = (*data.hitI)[i];
        double d = std::abs(m * (*data.zCenter)[index] - (*data.yCenter)[index] + b) / std::sqrt(m * m + 1);
        fom += std::abs((*data.hitR)[i] - d);
    }
    return fom;
}

std::pair<double, double> FitTrackFromRadii(const std::vector<double>& zCenter, const std::vector<double>& yCenter,
                       const std::vector<int>& hitI, const std::vector<double>& hitR) {
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    FitData data = {&zCenter, &yCenter, &hitI, &hitR};
    ROOT::Math::Functor functor([&data](const double* params) { return CalcFOM(params, data); }, 2);

    minimizer->SetFunction(functor);
    minimizer->SetVariable(0, "m", 1.0, 0.1);
    minimizer->SetVariable(1, "zOrigin", 0.0, 0.1);

    minimizer->Minimize();

    const double* results = minimizer->X();
    double m = results[0];
    double zOrigin = results[1];

    delete minimizer;

    return std::make_pair(m, zOrigin);
}

void PlotTrack(ComponentAnalyticField* cmp, DriftLineRKF* drift, TrackHeed* track, Sensor* sensor, std::vector<double> trk_y,std::vector<double> trk_z) {

  /*
  TCanvas cD("cD", "Drift View", 600, 600);
  ViewDrift driftView;
  driftView.SetCanvas(&cD);

  ViewCell cellView;
  cellView.SetCanvas(&cD);
  cellView.SetComponent(cmp);

  drift->EnablePlotting(&driftView);
  track->EnablePlotting(&driftView);
  cellView.SetComponent(cmp);
  cellView.SetCanvas(&cD);

  cellView.Plot2d();
  driftView.Plot(true, true);
  */

  TCanvas* canvas = new TCanvas("canvas", "Straw Tube Plot", 800, 600);
  std::vector<double> zCenter = {0.0},yCenter = {0.0};

  canvas->SetFixedAspectRatio(); // Ensure equal scaling on both axes
  canvas->DrawFrame(-2, -2, 2, 2); // Set the drawing frame (xMin, yMin, xMax, yMax)

    // Draw the circles (straws)
    for (size_t i = 0; i < yCenter.size(); ++i) {
        TEllipse* straw = new TEllipse(zCenter[i], yCenter[i], rTube);
        straw->SetFillStyle(0); // No fill
        straw->SetLineColor(kBlack);
        straw->Draw();
    }

    // Update the canvas to render the drawing
    canvas->Update();

  TGraph* graph = new TGraph(trk_z.size(), trk_z.data(), trk_y.data());

  // Set the marker style for the points (e.g., a circle)
  graph->SetMarkerStyle(21);  // 21: circle
  graph->SetMarkerColor(kRed);  // Red points
  graph->SetMarkerSize(3);  // Marker size

  // Draw the graph on the canvas
  graph->Draw("same");

  canvas->Update();

  //TCanvas cS = new TCanvas("cs", "", 600, 600);
  //sensor->PlotSignal("s", &cS);

  cin.get();

}


std::tuple<double, double, double, std::vector<double>, std::vector<double>, std::vector<double>> SetTrack(ComponentAnalyticField* cmp, DriftLineRKF* drift, Sensor* sensor, TrackHeed* track, double energy, double xdist) {


  string particle = "proton";

  track->SetParticle(particle);
  track->SetKineticEnergy(energy);
  track->SetSensor(sensor);
  drift->SetSensor(sensor);
  drift->EnableIonTail();

  double x0 = xdist;//rTube*dis(gen);
  double y0 = -sqrt(rTube*rTube - x0*x0);
  double DOCA= 999999, simR;
  track->NewTrack(x0, y0, 0, 0, 0, 1, 1);

  std::vector<double> trk_x,trk_y,trk_z;

  for (const auto& cluster : track->GetClusters()) {
      //cout<<cluster.x<<" "<<cluster.y<<" "<<cluster.z<<" "<<cluster.t<<endl;

      trk_x.push_back(cluster.x);
      trk_y.push_back(cluster.y);
      trk_z.push_back(cluster.z);

      simR = sqrt((cluster.y*cluster.y) + (cluster.z*cluster.z));//With reference to straw at 0 0

      if (simR<DOCA) {
          DOCA=simR;
      }

      for (const auto& electron : cluster.electrons) {
        drift->DriftElectron(electron.x, electron.y, electron.z, electron.t);
        //cout<<electron.x<<" "<<electron.y<<" "<<electron.z<<" "<<electron.t<<endl;
      }

    }

    // Signal processing
        std::vector<double> x(SigBins);
        std::vector<double> y(SigBins);
        std::vector<double> x_basket;

        // Collect signal
        for (int i = 0; i < SigBins; ++i) {
            x[i] = SigI + i;
            y[i] = -1*sensor->GetSignal("s", i);
            if (y[i] > 0.5) {
                x_basket.push_back(x[i]);
            }
            //cout<<x[i]<<" "<<y[i]<<endl;
        }

    auto [width, full_width] = SignalProcessing(x,y, x_basket);


    //SignalPlot(x,y);

    return std::make_tuple(DOCA, width, full_width, trk_x, trk_y, trk_z);

}

int main(int argc, char* argv[]) {

  TApplication app("app", &argc, argv);

  string P10 = "ar_90_co2_10_atms_3.gas";
  MediumMagboltz gas("ar", 90., "co2", 10.);
  SetGas(P10, &gas);

  ComponentAnalyticField cmp;
  SetWire(&cmp, &gas);

  Sensor sensor;
  SetSensor(&sensor, &cmp);

  double energy = 100e6;
  TrackHeed track;
  DriftLineRKF drift;

  //Event loop scanning full x dimension
  //std::vector<double> DOCA_l, width_l, full_width_l;
  std::ofstream outfile1("width.csv", std::ios::app);
  std::ofstream outfile2("full_width.csv", std::ios::app);

  for (int i = -99; i<99; i++) {
      double j = 0.7;//i/100.0;
    auto [DOCA, width, full_width, pos_x, pos_y, pos_z] = SetTrack(&cmp, &drift, &sensor ,&track, energy, j);
    //DOCA_l.push_back(DOCA);
    //width_l.push_back(width);
    //full_width_l.push_back(full_width);
    outfile1<<width<<endl;
    outfile2<<full_width<<endl;
    cout<<j<<endl;
  }

  outfile1.close();
  outfile2.close();

  //PlotTrack(&cmp, &drift, &track, &sensor, pos_y, pos_z);

  app.Run();
  return 0;
}
