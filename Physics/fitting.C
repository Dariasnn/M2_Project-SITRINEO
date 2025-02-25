#include <iostream>
#include <cmath>
#include <vector>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TRandom3.h>
#include <ctime>

// in root run .L fitting.C
// then any following function :
// trajectory() 
// trajectory_scattered()
// fitting()



void sitritrajectorybended(double *z, double *x, double *Par){
  // calculates de x position for a given z, for a charged particle in a simple uniform B field
  // execute trajectory() to get the trajectory
  // execute trajectory_error() to see the error on the impact position in the 3rd plane if the wrong B field is used (with a precision of 0.05T)
  // execute trajectory_comparison() to compare the trajectory in a uniform field and a field composed of three steps

  double thetax = Par[0];
  double p = Par[1];
  double zA = Par[2];
  double zB = Par[3];
  double q = Par[4];
  double B = Par[5];

  if(B==0.0){B=1e-10;} // Avoid division by zero

  double R = fabs(p/(q*3e-4*B)); // Radius of curvature
  double xC = R*cos(q*thetax);
  double Delta = R*R - std::pow(zB - zA - R*sin(q*thetax), 2.0);
  double l, gamma, thetaxB;
  static double xA, xB;

  // Trajectory calculation
  if (z[0] < zA) { // Before B-field
    x[0] = tan(thetax) * z[0] + 30000*tan(thetax); //15000 à mettre en variable plus tard
  } 
  else if (z[0] == zA) { // Start B field
    xA = tan(thetax) * (z[0]) + 30000*tan(thetax);
  } 

  else if (z[0] < zB) { // Inside B-field
    Delta = R * R - std::pow(z[0] - zA - R * sin(q * thetax), 2.0);
    if (Delta < 0) {
      x[0] = 0;
    } 
    else {
      x[0] = xA + q * (-xC + sqrt(Delta));
    }
  } 
  else if (z[0] == zB) { // Start B field
    Delta = R * R - std::pow(z[0] - zA - R * sin(q * thetax), 2.0);
    if (Delta < 0) {
      xB = 0;
    } 
    else {
      xB = xA + q * (-xC + sqrt(Delta));
    }
  } 
  else { // After B-field
    l = sqrt(std::pow(zA - zB, 2.0) + std::pow(xA - xB, 2.0));
    gamma = -acos(1 - (l * l) / (2 * R * R));
    thetaxB = thetax + q * gamma;
    x[0] = tan(thetaxB) * (z[0] - zB) + xB;
  }
  return;
}


void sitritrajectoryscattered(double *z, double *x, double *Par, double *planeZ, double *theta) {
  // calculates de x position for a given z, for a charged particle in a simple uniform B field that undergoes multiple scattering in each plane
  // execute trajectory_scattered() to get the trajectory
  // execute simple_vs_scattered() to see the difference when we add multiple scattering

  // Parameters from input
  double thetax = Par[0];     // Initial X angle
  double p = Par[1];          // Momentum (MeV/c)
  double zA = Par[2];         // Start of magnetic field
  double zB = Par[3];         // End of magnetic field
  double q = Par[4];          // Particle charge
  double B = Par[5];          // Magnetic field strength (T)

  double plane1 = planeZ[0];
  double plane2 = planeZ[1];
  double plane3 = planeZ[2];
  double plane4 = planeZ[3];


  // Multiple scattering angles for each plane, for now random
  double theta_1 = theta[0];
  double theta_2 = theta[1];
  double theta_3 = theta[2];
  double theta_4 = theta[3];

  if (B == 0.0){B = 1e-10;} // Avoid division by zero

  double R = fabs(p / (q * 3e-4 * B)); // Radius of curvature in the B-field
  double xC = R * cos(q * (theta_1+theta_2));     // X center of the circular trajectory
  double Delta = R*R - pow(zB - zA - R*sin(q*thetax), 2.0);
  static double l, gamma, thetaxB;

  static double xA = 0;
  static double xB = 0;
  static double x3 = 0;

  // Trajectory calculation with bending and multiple scattering
  if (z[0] < plane1) { 
      // Before first plane (no scattering)
      x[0] = tan(thetax) * (z[0]) + 30000*tan(thetax);
  } 
  else if (z[0] < plane2) {
      // After first plane 
      x[0] = tan(theta_1) * z[0]  + 30000*tan(thetax) - plane1*(tan(theta_1)-tan(thetax));
  } 
  else if (z[0] < zA) {
      double scatter = theta_2;
      x[0] = tan(scatter) * z[0] + 30000*tan(thetax) - plane2*(tan(scatter)-tan(theta_1))- plane1*(tan(theta_1)-tan(thetax));
  }
  else if (z[0] < plane3) {
      // After second plane
      double scatter = theta_2;
      if (z[0] == zA) { 
          xA = x[0];
	  //std::cout<<"xA="<<xA<<std::endl;
          //std::cout << ": theta_1=" << theta_1 << ", theta_2=" << theta_2 << std::endl;
      } 
      else if (z[0] > zA && z[0] <= zB) { // Inside B-field
      Delta = R * R - std::pow(z[0] - zA - R * sin(q * (theta_1 + scatter)), 2.0);
        if (Delta < 0) {
          x[0] = 0;
        } 
        if (z[0] == zB) { 
          xB = x[0];
	  //std::cout<<"xB="<<xB<<std::endl;
	} 
        else {
          x[0] = xA + q * (-xC + sqrt(Delta));
        }
      } 
      else {
        l = sqrt(std::pow(zA - zB, 2.0) + std::pow(xA - xB, 2.0));
        gamma = -acos(1 - (l * l) / (2 * R * R));
        thetaxB = (scatter) + q * gamma;
        x[0] = tan(thetaxB) * (z[0]-zB)+xB;
      }
    } 
  else if (z[0] == plane3){ x3 = x[0];}
  else if (z[0] < plane4) {
    // After third plane
    double scatter = theta_3;
    x[0] = tan(scatter) * (z[0] - plane3) + x3;
  } 
  else {
    // After fourth plane
    double scatter = theta_4; 
    x[0] = tan(scatter) * (z[0] - plane4) + x3 + tan(theta_3) * (plane4 - plane3);
  }
  return;
}


//******************************************************************************
// 
// Function to draw the simple trajectory
//
//******************************************************************************

void trajectory() {
  double Par[6] = {
    0.02,     // thetax (rad)
    1.32,    // p (MeV/c) - momentum
    -5000.0,   // zA (um) start of B field
    5000.0,   // zB (um) end of B field 
    1.0,      // q (charge) = +1
    0.1       // B (T) magnetic field
  };

  // Create graph for the trajectory
  TGraph *graph = new TGraph();

  double x;
  int pointIndex = 0;

  for (double z = -30000; z <= 20000; z += 10) {
    double zArray[1] = {z};
    sitritrajectorybended(zArray, &x, Par);
    graph->SetPoint(pointIndex, z, x);
    pointIndex++;
  }

  //Graph range
  graph->GetXaxis()->SetLimits(-30000, 20000); //max 4cm
  graph->GetYaxis()->SetRangeUser(-10000, 10000); //2cm hauteur pour voir les plans


  // Create canvas
  TCanvas *c1 = new TCanvas("c1", "Trajectory through Magnetic Field", 800, 600);

  c1->SetRightMargin(0.2); 

  graph->SetTitle("Particle trajectory through a magnetic field;Z (um);X (um)");
  graph->SetLineColor(kBlue);
  graph->SetLineWidth(2);
  graph->Draw("AL");

  //add lines and bfield limit
  //Planes 

  double planeHeight = 19872; // 19.872 mm in µm
  double planeThickness = 50; // 50 µm thickness

  // Plane positions
  double planeZ[4] = {
    -15000.0,        // Plane 1
    -10000.0,        // Plane 2
    10000.0,         // Plane 3
    15000.0          // Plane 4
  };


  // Draw each plane
  for (int i = 0; i < 4; i++) {
    TLine *planeLine = new TLine(
      planeZ[i] - planeThickness / 2, -planeHeight / 2,
      planeZ[i] + planeThickness / 2, planeHeight / 2
    );
    planeLine->SetLineColor(kBlack);
    planeLine->SetLineStyle(1);
    planeLine->SetLineWidth(2);
    planeLine->Draw("same");
  }

  //Bfield
  TLine *bFieldStart = new TLine(-5000, -planeHeight/2, -5000, planeHeight/2);
  bFieldStart->SetLineColor(kMagenta);
  bFieldStart->SetLineStyle(2);
  bFieldStart->Draw("same");

  TLine *bFieldEnd = new TLine(5000, -planeHeight/2, 5000, planeHeight/2);
  bFieldEnd->SetLineColor(kMagenta);
  bFieldEnd->SetLineStyle(2);
  bFieldEnd->Draw("same");

  //Legend
  TLegend *legend = new TLegend(0.81, 0.65, 0.98, 0.85);
  legend->AddEntry(graph, "Trajectory", "l");

  // Magnetic Field Region
  TLine *legendBField = new TLine(0, 0, 1, 1);
  legendBField->SetLineColor(kMagenta);
  legendBField->SetLineStyle(2);
  legendBField->SetLineWidth(2);
  legend->AddEntry(legendBField, "B Field", "l");

  // Sensing Planes
  TLine *legendPlane = new TLine(0, 0, 1, 1);
  legendPlane->SetLineColor(kBlack);
  legendPlane->SetLineStyle(1);
  legendPlane->SetLineWidth(2);
  legend->AddEntry(legendPlane, "Planes", "l");

  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(1001);
  legend->Draw();
}

//******************************************************************************
// 
// Function to draw the scattered trajectory
//
//******************************************************************************

void trajectory_scattered() {
  double Par[6] = {
    0.02,     // thetax (rad)
    1.32,    // p (MeV/c) - momentum
    -5000.0,   // zA (um) start of B field
    5000.0,   // zB (um) end of B field 
    1.0,      // q (charge) = +1
    0.1       // B (T) magnetic field
  };

  // Plane positions
  double planeZ[4] = {
    -15000.0,        // Plane 1
    -10000.0,        // Plane 2
    10000.0,         // Plane 3
    15000.0          // Plane 4
  };

  // Physical constants
  double m_e = 0.511;        // MeV (electron rest mass)
  double X0 = 21.82;         // g/cm^2 (radiation length for Silicon)
  double rho = 2.33;         // g/cm^3 (density of Silicon)
  double thickness = 50e-4;  // cm (50 µm thickness per plane)
  double E = sqrt(Par[1] * Par[1] + m_e * m_e); // Total energy
  double beta = Par[1] / E;                // Velocity fraction of c
  std::cout<<"p="<<Par[1]<<std::endl;

  TRandom3 randGen;  // Random generator for scattering
  randGen.SetSeed(time(0));  // Unique seed based on current time

  // Highland formula for Multiple Coulomb Scattering (RMS angle)
  double theta0 = (13.6 / (beta * Par[1])) * fabs(Par[4]) * sqrt(thickness / X0) * (1 + 0.038 * log(thickness / X0));

  // Multiple scattering angles for each plane
  double theta[4] = {
    randGen.Gaus(0, theta0),        // Plane 1
    randGen.Gaus(0, theta0),        // Plane 2
    randGen.Gaus(0, theta0),        // Plane 3
    randGen.Gaus(0, theta0)         // Plane 4
  };
  // Create graph for the trajectory
  TGraph *graph = new TGraph();

  double x;
  int pointIndex = 0;

  for (double z = -30000; z <= 20000; z += 100) {
    double zArray[1] = {z};
    sitritrajectoryscattered(zArray, &x, Par, planeZ, theta);
    graph->SetPoint(pointIndex, z, x);
    pointIndex++;
  }

  //Graph range
  graph->GetXaxis()->SetLimits(-30000, 20000); //max 4cm
  graph->GetYaxis()->SetRangeUser(-10000, 10000); //2cm hauteur pour voir les plans


  // Create canvas
  TCanvas *c1 = new TCanvas("c1", "Trajectory through Magnetic Field", 800, 600);

  c1->SetRightMargin(0.2); 

  graph->SetTitle("Particle trajectory undergoing multiple scattering;Z (um);X (um)");
  graph->SetLineColor(kBlue);
  graph->SetLineWidth(2);
  graph->Draw("AL");

  //add lines and bfield limit
  //Planes 

  double planeHeight = 19872; // 19.872 mm in µm
  double planeThickness = 50; // 50 µm thickness


  // Draw each plane
  for (int i = 0; i < 4; i++) {
    TLine *planeLine = new TLine(
      planeZ[i] - planeThickness / 2, -planeHeight / 2,
      planeZ[i] + planeThickness / 2, planeHeight / 2
    );
    planeLine->SetLineColor(kBlack);
    planeLine->SetLineStyle(1);
    planeLine->SetLineWidth(2);
    planeLine->Draw("same");
  }

  //Bfield
  TLine *bFieldStart = new TLine(-5000, -planeHeight/2, -5000, planeHeight/2);
  bFieldStart->SetLineColor(kMagenta);
  bFieldStart->SetLineStyle(2);
  bFieldStart->Draw("same");

  TLine *bFieldEnd = new TLine(5000, -planeHeight/2, 5000, planeHeight/2);
  bFieldEnd->SetLineColor(kMagenta);
  bFieldEnd->SetLineStyle(2);
  bFieldEnd->Draw("same");

  //Legend
  TLegend *legend = new TLegend(0.81, 0.65, 0.98, 0.85);
  legend->AddEntry(graph, "Trajectory", "l");

  // Magnetic Field Region
  TLine *legendBField = new TLine(0, 0, 1, 1);
  legendBField->SetLineColor(kMagenta);
  legendBField->SetLineStyle(2);
  legendBField->SetLineWidth(2);
  legend->AddEntry(legendBField, "B Field", "l");

  // Sensing Planes
  TLine *legendPlane = new TLine(0, 0, 1, 1);
  legendPlane->SetLineColor(kBlack);
  legendPlane->SetLineStyle(1);
  legendPlane->SetLineWidth(2);
  legend->AddEntry(legendPlane, "Planes", "l");

  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(1001);
  legend->Draw();
}


//******************************************************************************
// 
// Function to compare the uncertainties between the simple trajectory and the scattered trajectories
//
//******************************************************************************

void fitting() {

  double Par[6] = {
    0.02,     // thetax (rad)
    1.32,    // p (MeV/c) - momentum
    -5000.0,   // zA (um) start of B field
    5000.0,   // zB (um) end of B field 
    1.0,      // q (charge) = +1
    0.1,       // B (T) magnetic field
  };

  // Plane positions
  double planeZ[4] = {
    -15000.0,        // Plane 1
    -10000.0,        // Plane 2
    10000.0,        // Plane 3
    15000.0         // Plane 4
  };

  double planeHeight = 19872; // 19.872 mm in µm
  double planeThickness = 50; // 50 µm thickness

  // Physical constants
  double m_e = 0.511;        // MeV (electron rest mass)
  double X0 = 21.82;         // g/cm^2 (radiation length for Silicon)
  double rho = 2.33;         // g/cm^3 (density of Silicon)
  double thickness = 50e-4;  // cm (50 µm thickness per plane)
  double E = sqrt(Par[1] * Par[1] + m_e * m_e); // Total energy
  double beta = Par[1] / E;                // Velocity fraction of c

  std::vector<double> impact_values; //impact point
  std::vector<double> Run_number;
  std::vector<double> Run_number2;
  std::vector<double> impact_values_bended;

  // Highland formula for Multiple Coulomb Scattering (RMS angle)
  double theta0 = (13.6 / (beta * Par[1])) * fabs(Par[4]) * sqrt(thickness / X0) * (1 + 0.038 * log(thickness / X0));
  
  // Create a MultiGraph to hold both graphs
  TMultiGraph *mg = new TMultiGraph();
  TLegend *legend = new TLegend(0.81, 0.45, 0.98, 0.85);

  TRandom3 randGen;  // Random generator


  //Plotting mutliple runs to the compute the mean error
  for (double i = 0; i <= 50; i += 1) {

    bool localStopPlot = false; //to stop plotting if the particle won't go to the 3rd plane

    TGraph *graph_bended = new TGraph(); 

    double x, x2;
    int pointIndex = 0;

    randGen.SetSeed(time(0)+i);  // Unique seed based on current time
    Par[0] = randGen.Gaus(0, 0.005);
    // std::cout<<"randtheta0="<<Par[0]<<std::endl;

    // Multiple scattering angles for each plane
    double theta[4] = {
      randGen.Gaus(0, theta0),        // Plane 1
      randGen.Gaus(0, theta0),        // Plane 2
      randGen.Gaus(0, theta0),        // Plane 3
      randGen.Gaus(0, theta0)         // Plane 4
    };

    for (double z = -30000; z <= 20000; z += 10) { // 5 cm total
      if (localStopPlot) break;

      double zArray[1] = {z};
      sitritrajectoryscattered(zArray, &x, Par, planeZ, theta);

      double zArray2[1] = {z};
      sitritrajectorybended(zArray2, &x2, Par);

      if (z == 10000 && x > -planeHeight/2){ //if impact in the plane 3
        impact_values.push_back(x);
        // std::cout << "impact with plane 3 for run="<< i << " x="<< x << std::endl;
        Run_number.push_back(i);
      }

      if (z == 10000 && x2 > -planeHeight/2){ //if impact in the plane 3
        impact_values_bended.push_back(x2);
        // std::cout << "impact with plane 3 for run="<< i << " x2="<< x2 << std::endl;
        Run_number2.push_back(i);
      }

      graph_bended->SetPoint(pointIndex, z, x2);
      pointIndex++;
    }

    graph_bended->SetTitle("Simple particle trajectories;Z (um);X (um)");
    graph_bended->SetLineWidth(2);
    mg->Add(graph_bended);
    legend->AddEntry(graph_bended, Form("Run = %.1f", i), "l");

  }

  // New canvas for plotting the trajectory associated with sitritrajectorybended
  TCanvas *c4 = new TCanvas("c4", "Bended Trajectory", 800, 600);

  c4->SetRightMargin(0.2); 
  mg->GetXaxis()->SetLimits(-30000, 20000);
  mg->GetYaxis()->SetRangeUser(-10000, 10000); //2cm hauteur pour voir les plans

  mg->Draw("A pmc plc"); //to have beautiful bird colors

  c4->Update();

  // Histogram of Impact Positions
  TCanvas *c2 = new TCanvas("c2", "Impact Position Histogram", 800, 600);
  TH1F *hist = new TH1F("hist", "Histogram of Impact Positions;Impact Position (um);Number of Entries", 100, -10000, 10000);
  for (double impact : impact_values) {
      hist->Fill(impact);
  }
  hist->SetFillColor(kBlue);
  hist->SetLineColor(kBlack);
  hist->Draw();
  c2->Update();

  // Fit the histogram with a Gaussian function
  TF1 *gaussFit = new TF1("gaussFit", "gaus", -10000, 10000);
  hist->Fit(gaussFit, "R");
  gaussFit->SetLineColor(kRed);
  gaussFit->Draw("same");

  c2->Update();

  // Histogram of Impact Positions for Bended Trajectories
  TCanvas *c3 = new TCanvas("c3", "Impact Position Histogram (Bended)", 800, 600);
  TH1F *hist_bended = new TH1F("hist_bended", "Histogram of Impact Positions (Bended);Impact Position (um);Number of Entries", 100, -10000, 10000);
  for (double impact : impact_values_bended) {
      hist_bended->Fill(impact);
  }
  hist_bended->SetFillColor(kGreen);
  hist_bended->SetLineColor(kBlack);
  hist_bended->Draw();
    
  TF1 *gaussFit_bended = new TF1("gaussFit_bended", "gaus", -10000, 10000);
  hist_bended->Fit(gaussFit_bended, "R");
  gaussFit_bended->SetLineColor(kRed);
  gaussFit_bended->Draw("same");
    
  c3->Update();

}
