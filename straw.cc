#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2D.h>
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
double remv_tol = 0.10;
double func_tol = 0.15;
double parstep = 5;

TH1D *RecoZ = new TH1D("RecoZ", "RecoZ",10000,-100,100);
TH2D *WidthRad = new TH2D("WidthRad", "WidthRad",100,-0.9,0.9,5000,0,5000);
TH2D *FullWidthRad = new TH2D("FullWidthRad", "FullWidthRad",100,-0.9,0.9,5000,0,5000);

//Wire radius [cm]
double rWire = 25e-4;

//Tube radius [cm]
double rTube = 0.5;

//Tube speration [cm]
double VerticalOffset = 10;
double vert_offset = 2.6;

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

double slopefunc(double *x, double *params) {
    double xx = x[0];          // x is the independent variable
    double A = params[0];      // Amplitude
    double B = params[1];     // Mean
    double C = params[2];  // Standard deviation
    double D = params[3];  // Plateau control parameter

    // Gaussian numerator
    return A*exp(B*xx)+C+D*xx;
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


void PlotStraws(const std::vector<double>& yCenter, const std::vector<double>& zCenter) {
    // Check if the inputs are valid

    TCanvas* canvas = new TCanvas("canvas", "Straws and Track", 800, 800);
    // Create a canvas
    canvas->SetFixedAspectRatio(); // Ensure equal scaling on both axes
    canvas->DrawFrame(-20, -20, 20, 20); // Set the drawing frame (xMin, yMin, xMax, yMax)

    // Draw the circles (straws)
    for (size_t i = 0; i < yCenter.size(); ++i) {
        TEllipse* straw = new TEllipse(zCenter[i], yCenter[i], rTube);
        straw->SetFillStyle(0); // No fill
        straw->SetLineColor(kBlack);
        straw->Draw();
    }

    // Update the canvas to render the drawing
    canvas->Update();
}

std::pair<std::vector<double>, std::vector<double>> GetStrawCenters(double strawRadius, double verticalOffset) {

    std::vector<double> yCenter, zCenter;

    // Calculate the number of tubes needed to cover +/-10 cm
    int nTubes = static_cast<int>(5 / strawRadius) + 1;

    // Generate the first layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
        //cout<<yCenter[i]<<" "<<zCenter[i]<<endl;
    }

    // Generate the second layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + std::sqrt(3) * strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius + strawRadius);
        //cout<<yCenter[i]<<" "<<zCenter[i]<<endl;
    }

    // Generate the third layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
        //cout<<yCenter[i]<<" "<<zCenter[i]<<endl;
    }

     // Generate the fourth layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius + 2* strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius +
        strawRadius);
        //cout<<yCenter[i]<<" "<<zCenter[i]<<endl;
    }

    // Generate the fourth layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius + 4* strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
        //cout<<yCenter[i]<<" "<<zCenter[i]<<endl;
    }

    // Generate the fourth layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius + 6* strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius +
        strawRadius);
        //cout<<yCenter[i]<<" "<<zCenter[i]<<endl;
    }

    // Generate the fourth layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius + 8* strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
        //cout<<yCenter[i]<<" "<<zCenter[i]<<endl;
    }

    // Generate the fourth layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius + 10* strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius +
        strawRadius);
        //cout<<yCenter[i]<<" "<<zCenter[i]<<endl;
    }

    // Generate the fourth layer
    for (int i = 0; i <= 2 * nTubes; ++i) {
        yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius + 12* strawRadius);
        zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
        //cout<<yCenter[i]<<" "<<zCenter[i]<<endl;
    }

    return std::make_pair(zCenter, yCenter);

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

void PlotTrack(std::vector<double> trk_y,std::vector<double> trk_z) {

  TGraph* graph = new TGraph(trk_z.size(), trk_z.data(), trk_y.data());

  // Set the marker style for the points (e.g., a circle)
  graph->SetMarkerStyle(21);  // 21: circle
  graph->SetMarkerColor(kRed);  // Red points
  graph->SetMarkerSize(3);  // Marker size

  // Draw the graph on the canvas
  graph->Draw("same");

}


std::tuple<double, double, double, std::vector<double>, std::vector<double>, std::vector<double>> SetTrack(ComponentAnalyticField* cmp, DriftLineRKF* drift, Sensor* sensor, TrackHeed* track, double energy, double xdist, std::vector<double> dirvect, std::vector<double> posvect) {

  string particle = "proton";

  track->SetParticle(particle);
  track->SetKineticEnergy(energy);
  track->SetSensor(sensor);
  drift->SetSensor(sensor);
  drift->EnableIonTail();

 //double x0 = xdist;//rTube*dis(gen);
  //double y0 = -sqrt(rTube*rTube - x0*x0);
  double DOCA = 999999, simR;
  track->NewTrack(posvect[0], posvect[1], posvect[2], 0, dirvect[0], dirvect[1], dirvect[2]);

  std::vector<double> trk_x,trk_y,trk_z;

  /*
  if (track->GetClusterDensity()==0) {
      trk_x.push_back(0); trk_y.push_back(0); trk_z.push_back(0);
      return std::make_tuple(0.5, 300, 500, trk_x, trk_y, trk_z);
    }
  */

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

std::vector<double> plotwidth(std::vector<double> width_l, std::vector<double> full_width_l, std::vector<double> DOCA_l) {

 std::vector<double> params;

  TGraph *GWidth = new TGraph();
  for (size_t i = 0; i < width_l.size(); ++i) {
    GWidth->SetPoint(i,width_l[i],DOCA_l[i]);
 }

 TF1 *slope = new TF1("slopefunc", slopefunc,
                      *min_element(std::begin(width_l),std::end(width_l)),
                      *max_element(std::begin(width_l),std::end(width_l)),4);

 TF1 *slope_remv = new TF1("slopefunc", slopefunc,
                      *min_element(std::begin(width_l),std::end(width_l)),
                      *max_element(std::begin(width_l),std::end(width_l)),4);

 // Fit the graph
 GWidth->Fit(slope, "RQ+");

  TCanvas *canvas = new TCanvas("canvas", "Two TGraphs", 800, 600);

  std::vector<double> width_remv;
  std::vector<double> DOCA_remv;

  for (int i = 0; i<width_l.size(); i++) {
      double func = slope->GetParameter(0)*exp(slope->GetParameter(1)*width_l[i])+slope->GetParameter(2)+slope->GetParameter(3)*width_l[i];

      //cout<<width_l[i]<<" "<<func<<endl;
      double func_diff = abs(DOCA_l[i]-func)/DOCA_l[i];

      if (func_diff<func_tol) {
          width_remv.push_back(width_l[i]);
          DOCA_remv.push_back(DOCA_l[i]);
          //cout<<"no_remv"<<endl;

    }

    }

    TGraph *GWidth_remv = new TGraph();

  for (size_t i = 0; i < width_remv.size(); ++i) {
    GWidth_remv->SetPoint(i,width_remv[i],DOCA_remv[i]);
  }

    GWidth_remv->Fit(slope_remv, "RQ+");
    //GWidth_remv->Draw("AP");
    //slope_remv->Draw("P SAME");

    params.push_back(slope_remv->GetParameter(0));
    params.push_back(slope_remv->GetParameter(1));
    params.push_back(slope_remv->GetParameter(2));
    params.push_back(slope_remv->GetParameter(3));

    return params;

}

double getDOCA(double width, std::vector<double> params) {
    return  params[0]*exp(params[1]*width)+params[2]+params[3]*width;
}

//, DriftLineRKF* drift, , TrackHeed* track


std::vector<double> GetDTCorr(Sensor* sensor, ComponentAnalyticField* cmp, double energy, std::vector<double> dirvect) {

    std::vector<double> width_l, full_width_l, j_l, DOCA_l;

    for (double i =-99; i<99; i=i+parstep) {

    if ((i==66) || (i==-66) ||(i==-52) || (i ==52)|| (i==-38)|| (i==38)|| (i==-26)|| (i==26) || (i==-16)|| (i==16)|| (i==-15)|| (i==15)|| (i==-32.5)|| (i==32.5)|| (i==-26.5)|| (i==26.5)|| (i==-19.5)|| (i==19.5)||(i==-15.5)|| (i==15.5)||(i==-2.5)|| (i==2.5)) continue;

    double j = i/100.0;

    TrackHeed track;
    DriftLineRKF drift;
    sensor->ClearSignal();

    auto [DOCA, width, full_width, pos_x, pos_y, pos_z] = SetTrack(cmp, &drift, sensor ,&track, energy, j, dirvect, {0, 0.01,0});

    width_l.push_back(width);
    full_width_l.push_back(full_width);
    j_l.push_back(j);
    DOCA_l.push_back(DOCA);
  }

    return plotwidth(width_l,full_width_l,DOCA_l);

}

double GetAvgResidual(Sensor* sensor, ComponentAnalyticField* cmp, std::vector<double> params, double energy, std::vector<double> dirvect) {

    double sum = 0;
    double step = 0.1, start = 0.1*rTube, stop = 0.9*rTube;

    for (double i=start; i<stop; i=i+step) {
        TrackHeed track;
        DriftLineRKF drift;
        sensor->ClearSignal();
        auto [DOCA, width, full_width, pos_x, pos_y, pos_z] = SetTrack(cmp, &drift, sensor ,&track, energy, i, dirvect, {0,0.01,0});

        sum+=abs(DOCA-getDOCA(width,params))/DOCA;

    }

    return sum/((stop-start)/step);

}

void SetCellArray(ComponentAnalyticField* cmp, MediumMagboltz* gas) {

    auto [ycell, zcell] = GetStrawCenters(rTube ,vert_offset);

    cmp->AddTube(rTube, vTube, 0, "tube 0");
    cmp->AddWire(0, ycell[0] , 2*rWire, vWire, "cell 0");
    cmp->SetMedium(gas);

    //PlotStraws(zcell,ycell);

}

std::vector<double> hemidir() {

    // Azimuthal angle (theta), uniformly distributed between 0 and 2pi
    double theta = 2.0 * M_PI * rand() / RAND_MAX;

    // Polar angle (phi), uniformly distributed in the upper hemisphere
    double cos_phi = (rand() / double(RAND_MAX));  // Uniformly distributed between 0 and 1
    double phi = acos(cos_phi);  // Convert to polar angle

    // Convert spherical coordinates to Cartesian coordinates
    return {sin(phi) * cos(theta),sin(phi) * sin(theta),cos(phi)  };

}

std::vector<double> sectordir() {

    // Polar angle (theta), uniformly distributed between 0 and pi
    double theta = 0.5 * M_PI * rand() / RAND_MAX + 0.25 * M_PI;

    // Convert spherical coordinates to Cartesian coordinates
    return {0, sin(theta), cos(theta)} ;

}

void PlotCellsTracks(std::vector<double> DOCA_l, std::vector<double> zCenter, std::vector<double> yCenter) {
    // Draw the circles (straws)
    for (size_t i = 0; i < yCenter.size(); ++i) {

        TEllipse* DOCA = new TEllipse(zCenter[i], yCenter[i], DOCA_l[i]);
        DOCA->SetFillStyle(0); // No fill
        DOCA->SetLineColor(kBlue);
        DOCA->Draw("same");
    }

}

double RandomBetween(double min, double max) {
    return min + (static_cast<double>(rand()) / RAND_MAX) * (max - min);
}

std::pair<double, double> GenerateRandomLineParameters() {
    // Seed the random number generator
    srand(static_cast<unsigned>(time(nullptr)));

    // Generate a random angle between pi/4 and 3pi/4
    static std::random_device rd; // Random device for seeding
    static std::mt19937 gen(rd()); // Mersenne Twister generator
    static std::uniform_real_distribution<> dis(3*M_PI / 8, 5 * M_PI / 8);

    const double angle = dis(gen);

    // Compute slope (m) of the line
    double m = std::tan(angle);

    // Compute a random y-intercept (c)
    double c = RandomBetween(-0.1, 0.1); // Example range for intercept

    return {m, c};
}

/*
std::pair<std::vector<double>, std::vector<double>> hits(double slope, double intercept, std::vector<double> ycell, std::vector<double> zcell) {

    std::vector<std::pair<double, double>> ystraws;
    std::vector<std::pair<double, double>> zstraws;
    std::vector<double> yhits, zhits;

    for (double step = 0; step<15; step=step+0.1) {

        double ypos = slope*step + intercept;

        for (int i = 0; i<zcell.size(); i++) {
            double diff = sqrt(pow(ycell[i]-ypos,2)+pow(zcell[i]-step,2));

            if (diff<rTube) {
                ystraws.push_back({ypos, ycell[i]});
                zstraws.push_back({step, zcell[i]});
            }

        }

    }

    // Sort the vector based on the second element of the pairs
    std::sort(ystraws.begin(), ystraws.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    });

    std::sort(zstraws.begin(), zstraws.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    });

    // Use std::unique to remove pairs with duplicate second elements
    auto last = std::unique(ystraws.begin(), ystraws.end(), [](const auto& a, const auto& b) {
        return a.second == b.second;
    });

    auto last2 = std::unique(zstraws.begin(), zstraws.end(), [](const auto& a, const auto& b) {
        return a.second == b.second;
    });

    // Resize the vector to remove redundant elements
    ystraws.erase(last, ystraws.end());
    zstraws.erase(last2, zstraws.end());

    // Iterate through the vector and collect the first elements
    for (const auto& pair : ystraws) {
        yhits.push_back(pair.first);
    }

    // Iterate through the vector and collect the first elements
    for (const auto& pair : zstraws) {
        zhits.push_back(pair.first);
    }

    return {yhits,zhits};

}

*/

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>,std::vector<double>> hits(double slope, double intercept, std::vector<double> ycell, std::vector<double> zcell) {

    std::vector<std::pair<double, double>> hits;
    std::vector<double> yhits,ystraw;
    std::vector<double> zhits,zstraw;

    for (int i = 0; i<zcell.size(); i++) {

        double A = 1 + slope*slope;
        double B = 2*(slope*intercept-slope*ycell[i]-zcell[i]);
        double C = ycell[i]*ycell[i]-rTube*rTube+zcell[i]*zcell[i]-2*intercept*ycell[i]+intercept*intercept;

        if (B*B-4*A*C<0) {continue;}

        if (B*B-4*A*C==0) {continue;}

        if (B*B-4*A*C>0) {

            double plus_z = (-B+sqrt(B*B-4*A*C))/(2*A);
            double minus_z = (-B-sqrt(B*B-4*A*C))/(2*A);
            double plus_y = slope*plus_z+intercept;
            double minus_y = slope*minus_z+intercept;

            //search for earliest solution
            double diffplus = sqrt(plus_z*plus_z+plus_y*plus_y);
            double diffminus = sqrt(minus_z*minus_z+minus_y*minus_y);

            //cout<<"test ("<<plus_z<<","<<plus_y<<") ("<<minus_z<<","<<minus_y<<") "<<endl;

            if (diffplus>diffminus) {
                yhits.push_back(minus_y);
                zhits.push_back(minus_z);
                ystraw.push_back(ycell[i]);
                zstraw.push_back(zcell[i]);
            }

            if (diffplus<diffminus) {
                yhits.push_back(plus_y);
                zhits.push_back(plus_z);
                ystraw.push_back(ycell[i]);
                zstraw.push_back(zcell[i]);
            }

        }


    }

    return std::make_tuple(yhits, zhits,ystraw, zstraw);

}

std::pair<double, double> FitRadii(const std::vector<double> z, const std::vector<double> y, std::vector<double> DOCA_l ) {

    // Set up the system of equations
    int n = z.size();
    double sum_z = 0, sum_y = 0, sum_zy = 0, sum_zz = 0;
    for (int i = 0; i < n; ++i) {
        sum_z += z[i];
        sum_y += y[i];
        sum_zy += z[i] * y[i];
        sum_zz += z[i] * z[i];
    }

    // Calculate slope (m) and intercept (b)
    double m = (n * sum_zy - sum_z * sum_y) / (n * sum_zz - sum_z * sum_z);
    double c = (sum_y - m * sum_z) / n;

    return std::make_pair(m,c);

}

int main(int argc, char* argv[]) {

  TApplication app("app", &argc, argv);

  //Get track dir
  std::vector<double> trk_y,trk_z;

  //auto [slope, intercept] = GenerateRandomLineParameters();
    double slope = 0.75, intercept = 0;
  for (double i = 0; i<15; i=i+0.1) {
    trk_z.push_back(i);
    trk_y.push_back(slope*i+intercept);
  }

  auto [ycell, zcell] = GetStrawCenters(rTube ,vert_offset);
  auto [yhits, zhits, ystraw, zstraw] = hits(slope, intercept, ycell, zcell);

  //Get Straw arrray
  PlotStraws(ycell, zcell);
  PlotTrack(trk_y,trk_z);

  //Setup Garfield sim for each hit
  std::vector<double> dirvect = {0, slope, 1};
  double energy = 100e6;

  string P10 = "ar_90_co2_10_atms_3.gas";
  MediumMagboltz gas("ar", 90., "co2", 10.);
  SetGas(P10, &gas);

  ComponentAnalyticField cmp;
  SetWire(&cmp, &gas);
  //SetCellArray(&cmp, &gas);

  Sensor sensor;
  SetSensor(&sensor, &cmp);

  std::vector<double> DOCA_l;
  std::vector<int> hitI;

  for (int i = 0; i<yhits.size(); i++) {

    TrackHeed track;
    DriftLineRKF drift;
    sensor.ClearSignal();

    std::vector<double> posvect = {0, yhits[i]-ystraw[i], zhits[i]-zstraw[i]};

    auto [DOCA, width, full_width, pos_x, pos_y, pos_z] = SetTrack(&cmp, &drift, &sensor ,&track, energy, rTube, dirvect, posvect);

    DOCA_l.push_back(DOCA);

  }

  PlotCellsTracks(DOCA_l, zstraw, ystraw);

  auto [fitslope, fitintercept] = FitRadii(zstraw,ystraw,DOCA_l);

  cout<<fitintercept<<" "<<0.0<<" "<<fitintercept*fitintercept<<endl;
  cout<<fitslope<<" "<<slope<<" "<<endl;
  cout<<(intercept-1*fitintercept)/fitslope<<" "<<0.0<<" "<<endl;


  app.Run();
  return 0;
}
