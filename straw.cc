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
int rowN = 6;

TH1D *RecoZ = new TH1D("RecoZ", "RecoZ",50000,-2,2);//[cm]
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

double RandomBetween(double min, double max) {
    return min + (static_cast<double>(rand()) / RAND_MAX) * (max - min);
}

void SetGas(string GasName, MediumMagboltz* gas) {

    gas->SetPressure(3 * 760.);
    gas->SetTemperature(293.15);
    gas->LoadGasFile(GasPath + GasName);
    gas->LoadIonMobility(GarPath + "Data/IonMobility_Ar+_Ar.txt");

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

    //cout<<"test 14"<<endl;

    //Get maximum charge and time when
    int amp_in = std::distance(y.begin(), max_element(std::begin(y), std::end(y)));
    double amp = y[amp_in];
    int UpTail = 1;

    //cout<<"test 15"<<endl;

    if (x_basket.size()>= 2) {

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

    }

    //cout<<"test 31"<<endl;
        //cout<<i*tstep<<" "<<y[i]<<" "<<move_avg<<" "<<var<<endl;
    return std::make_pair(UpTail*tstep-x_basket.front()*tstep,x_basket.back()-x_basket.front()); //[ns]
    }

    else {
            //cout<<"test 30"<<endl;
        return std::make_pair(5000,5000); //[ns]

    }

    //cout<<"test 16"<<endl;

    //cout<<UpTail*tstep<<" "<<x_basket.front()*tstep<<endl;


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

void PlotRadii(std::vector<double> DOCA_l, std::vector<double> radius, std::vector<double> yCenter, const std::vector<double> zCenter) {

    // Draw the circles (straws)
    for (size_t i = 0; i < yCenter.size(); ++i) {
        TEllipse* straw = new TEllipse(zCenter[i], yCenter[i], radius[i]-0.01, radius[i]);
        straw->SetLineColor(kRed);
        straw->Draw("same");
    }

    // Draw the circles (straws)
    for (size_t i = 0; i < yCenter.size(); ++i) {
        TEllipse* straw = new TEllipse(zCenter[i], yCenter[i], DOCA_l[i]-0.01, DOCA_l[i]);
        straw->SetFillColor(0);
        straw->Draw("same");
    }

}

std::pair<std::vector<double>, std::vector<double>> GetStrawCenters(double strawRadius, double verticalOffset) {

    std::vector<double> yCenter, zCenter;

    // Calculate the number of tubes needed to cover +/-10 cm
    int nTubes = static_cast<int>(5 / strawRadius) + 1;
    //int nTubesH = static_cast<int>(2 / strawRadius) + 1;
    int nTubesH = rowN;

    for (int j = 0; j<=2*nTubesH; j++) {

        for (int i = 0; i <= 2 * nTubes; ++i) {

            if (j==0) {
                yCenter.push_back(verticalOffset);
                zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
            }
        }

        if (j==1) {
            for (int i = 0; i <= 2 * nTubes; ++i) {
                yCenter.push_back(verticalOffset + std::sqrt(3) * strawRadius);
                zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius +     strawRadius);
            }
        }

        if ((j!=0) && (j!=1) && (j % 2 == 0)) {
            for (int i = 0; i <= 2 * nTubes; ++i) {
                int k = (j-2);
                yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius+ 2*k*strawRadius);
                zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
            }
        }

        if ((j!=0) && (j!=1) && (j % 2 != 0)) {
            for (int i = 0; i <= 2 * nTubes; ++i) {
                int k = (j-2);
                yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius + 2*k*strawRadius);
                zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius +
            strawRadius);
            }
        }

    }

    return std::make_pair(zCenter, yCenter);

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

  //cout<<"test 5"<<endl;

  track->SetParticle(particle);
  track->SetKineticEnergy(energy);
  track->SetSensor(sensor);
  drift->SetSensor(sensor);
  drift->EnableIonTail();

  //cout<<"test 6"<<endl;

 //double x0 = xdist;//rTube*dis(gen);
  //double y0 = -sqrt(rTube*rTube - x0*x0);
  double DOCA = 999999, simR;
  track->NewTrack(posvect[0], posvect[1], posvect[2], 0, dirvect[0], dirvect[1], dirvect[2]);

  //cout<<"test 7"<<endl;


  std::vector<double> trk_x,trk_y,trk_z;


  if (track->GetClusterDensity()==0) {
      trk_x.push_back(0); trk_y.push_back(0); trk_z.push_back(0);
      return std::make_tuple(0.5, 300, 500, trk_x, trk_y, trk_z);
    }


  //cout<<"test 8"<<endl;

  for (const auto& cluster : track->GetClusters()) {
      //cout<<cluster.x<<" "<<cluster.y<<" "<<cluster.z<<" "<<cluster.t<<endl;

      //cout<<"test 9"<<endl;

      trk_x.push_back(cluster.x);
      trk_y.push_back(cluster.y);
      trk_z.push_back(cluster.z);

      simR = sqrt((cluster.y*cluster.y) + (cluster.z*cluster.z));//With reference to straw at 0 0

      if (simR<DOCA) {
          DOCA=simR;
      }

      //cout<<"test 10"<<endl;

      for (const auto& electron : cluster.electrons) {
        drift->DriftElectron(electron.x, electron.y, electron.z, electron.t);
        //cout<<electron.x<<" "<<electron.y<<" "<<electron.z<<" "<<electron.t<<endl;
      }

    }


    // Signal processing
        std::vector<double> x(SigBins);
        std::vector<double> y(SigBins);
        std::vector<double> x_basket;

        //cout<<"test 11"<<endl;

        // Collect signal
        for (int i = 0; i < SigBins; ++i) {
            x[i] = SigI + i;
            y[i] = -1*sensor->GetSignal("s", i);
            if (y[i] > 0.5) {
                x_basket.push_back(x[i]);
            }
            //cout<<x[i]<<" "<<y[i]<<endl;
        }

    //cout<<"test 12"<<endl;

    auto [width, full_width] = SignalProcessing(x,y, x_basket);

    //cout<<"test 13"<<endl;


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

  //TCanvas *canvas = new TCanvas("canvas", "Two TGraphs", 800, 600);

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

    params.push_back(slope_remv->GetParameter(0));
    params.push_back(slope_remv->GetParameter(1));
    params.push_back(slope_remv->GetParameter(2));
    params.push_back(slope_remv->GetParameter(3));

    return params;

}

double getDOCA(double width, std::vector<double> params) {
    return  params[0]*exp(params[1]*width)+params[2]+params[3]*width;
}

std::vector<double> GetDTCorr(Sensor* sensor, ComponentAnalyticField* cmp, double energy, std::vector<double> dirvect, std::vector<double> posvect) {

    std::vector<double> width_l, full_width_l, j_l, DOCA_l;

    //cout<<"test 17"<<endl;

    for (double i = 0.0; i<5; i=i+1) {

    double j = i/100.0;

    //cout<<"test 18"<<endl;

    TrackHeed track;
    DriftLineRKF drift;
    sensor->ClearSignal();

    //cout<<"test 19"<<endl;

    auto [DOCA, width, full_width, pos_x, pos_y, pos_z] = SetTrack(cmp, &drift, sensor ,&track, energy, rTube, dirvect, posvect);

    //cout<<"test 20"<<endl;

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

std::vector<double> hemidir() {

    // Azimuthal angle (theta), uniformly distributed between 0 and 2pi
    double theta = 2.0 * M_PI * rand() / RAND_MAX;

    // Polar angle (phi), uniformly distributed in the upper hemisphere
    double cos_phi = (rand() / double(RAND_MAX));  // Uniformly distributed between 0 and 1
    double phi = acos(cos_phi);  // Convert to polar angle

    // Convert spherical coordinates to Cartesian coordinates
    return {sin(phi) * cos(theta),sin(phi) * sin(theta),cos(phi)  };

}

std::pair<double, double> generateline() {

    static std::mt19937 gen(static_cast<unsigned>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()
    ));
    //static std::random_device rd;
    //static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-2.0, 2.0);

    double m = dis(gen);

    // Compute a random y-intercept (c)
    double c = RandomBetween(-0.5, 0.5); // Example range for intercept

    //fix intercept to 0 for test
    //double c = 0;

    // Convert spherical coordinates to unit vector
    return {m, c} ;

}

std::pair<double, double> GenerateRandomLineParameters() {
    // Seed the random number generator
    srand(static_cast<unsigned>(time(nullptr)));

    // Generate a random angle between pi/4 and 3pi/4
    static std::random_device rd; // Random device for seeding
    static std::mt19937 gen(rd()); // Mersenne Twister generator
    static std::uniform_real_distribution<> dis(5*M_PI / 16, 9 * M_PI / 16);

    const double angle = dis(gen);

    // Compute slope (m) of the line
    double m = std::tan(angle);

    // Compute a random y-intercept (c)
    double c = RandomBetween(-0.1, 0.1); // Example range for intercept

    return {m, c};
}

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

            ////cout<<"test ("<<plus_z<<","<<plus_y<<") ("<<minus_z<<","<<minus_y<<") "<<endl;

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

  //auto [ycell, zcell] = GetStrawCenters(rTube ,vert_offset);
  //PlotStraws(ycell, zcell);

  // Open the file in append mode
    std::ofstream outFile("check.csv", std::ios::app);

  double energy = 100e6;

  string P10 = "ar_90_co2_10_atms_3.gas";
  MediumMagboltz gas("ar", 90., "co2", 10.);
  SetGas(P10, &gas);

  ComponentAnalyticField cmp;
  SetWire(&cmp, &gas);
  //SetCellArray(&cmp, &gas);

  Sensor sensor;
  SetSensor(&sensor, &cmp);

  //for (double k = 0.3; k<; k=k+0.1 ) {
    double k = 0.5;
      rTube = k;

    for (int j = 0; j<5000; j++) {

  //Get track dir
  //std::vector<double> trk_y,trk_z;

  auto [slope, intercept] = generateline();

  cout<<slope<<","<<intercept<<endl;




  //double slope = 1, intercept = 0;

  //for (double i = 0; i<15; i=i+0.1) {

    //trk_z.push_back(i);
    //trk_y.push_back(slope*i+intercept);
  //}

  auto [ycell, zcell] = GetStrawCenters(rTube ,vert_offset);
  auto [yhits, zhits, ystraw, zstraw] = hits(slope, intercept, ycell, zcell);

  //Get Straw arrray
  //PlotStraws(ycell, zcell);
  //PlotTrack(trk_y,trk_z);

  //Setup Garfield sim for each hit
  std::vector<double> dirvect = {0, slope, intercept};
  std::vector<double> DOCA_l, Radii;

  for (int i = 0; i<yhits.size(); i++) {

    TrackHeed track;
    DriftLineRKF drift;
    sensor.ClearSignal();

    //cout<<"test 1"<<endl;

    std::vector<double> posvect = {0, yhits[i]-ystraw[i], zhits[i]-zstraw[i]};

    //cout<<"test 2"<<endl;

    auto [DOCA, width, full_width, pos_x, pos_y, pos_z] = SetTrack(&cmp, &drift, &sensor ,&track, energy, rTube, dirvect, posvect);

    //cout<<"test 3"<<endl;

    auto params = GetDTCorr(&sensor, &cmp, energy, dirvect, posvect);

    //cout<<"test 4"<<endl;

    double radius = getDOCA(width, params);

    Radii.push_back(radius);
    DOCA_l.push_back(DOCA);

  }

  //cout<<"test 22"<<endl;

  //PlotRadii(DOCA_l, Radii, ystraw, zstraw);

  //cout<<"test 23"<<endl;

  auto [fitslope, fitintercept] = FitRadii(zstraw,ystraw,Radii);

  //cout<<"test 24"<<endl;

  //cout<<j<<" "<<-fitintercept/fitslope<<" "<<-intercept/slope<<endl;

  //if (abs(-fitintercept/fitslope)>0.0001 && (abs(-fitintercept/fitslope)<2)) {

  //outFile<<rowN<<","<<j<<","<<-fitintercept/fitslope<<endl;
  cout<<j<<","<<-fitintercept/fitslope<<","<<-intercept/slope<<";"<<abs(-fitintercept/fitslope+intercept/slope)<<endl;
  outFile<<j<<","<<-fitintercept/fitslope<<","<<-intercept/slope<<";"<<(-fitintercept/fitslope+intercept/slope)<<endl;

  //}
  //cout<<-fitintercept/fitslope<<endl;

  //RecoZ->Fill(10.0*-fitintercept/fitslope); //[cm]

  //}
  }

  //RecoZ->Draw("hist");


   //outFile.close();


  app.Run();
  return 0;
}
