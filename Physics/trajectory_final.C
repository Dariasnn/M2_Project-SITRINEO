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

// in root run .L trajectory_final.C
// then any following function :
// trajectory() 
// trajectory_error() 
// trajectory3B() 
// trajectory_error3B() 
// trajectory_comparison() 
// trajectory_scattered()
// simple_vs_scattered()



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
    x[0] = tan(thetax) * z[0];
  } 
  else if (z[0] == zA) { // Start B field
    xA = tan(thetax) * (z[0]);
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


void sitritrajectorybended3B(double *z, double *x, double *Par) {
  // calculates de x position for a given z, for a charged particle in a three step B field
  // execute trajectory3B() to get the trajectory
  // execute trajectory_error3B() to see the error on the impact position in the 3rd plane if the wrong B field is used (with a precision of 0.05T)
  // execute trajectory_comparison() to compare the trajectory in a uniform field and a field composed of three steps

  // Parameters
  double thetax = Par[0];
  double p = Par[1];       // Momentum
  double zA = Par[2];      // Start of central field B
  double zB = Par[3];      // End of central field B
  double q = Par[4];       // Charge
  double B = Par[5];       // Central magnetic field
  double zA1 = Par[6];     // Start of first Bout field
  double zB2 = Par[7];    // End of second Bout field
  double Bout = Par[8];   // Magnetic field outside
  
  // Avoid division by zero
  if (Bout == 0.0) Bout = 1e-10;
  if (B == 0.0) B = 1e-10;

  double gamma1, gamma2, gamma3;
  // Store boundary positions
  static double xA1 = 0, xA = 0; 
  static double xB = 0;
  static double xB2 = 0;
  static double thetaxB1 = 0, thetaxB2 = 0, thetaxB3 = 0;
  static double l1 = 0, l2 = 0, l3 = 0;

  // Radius of curvature
  double R_in = fabs(p / (q * 3e-4 * B));     // Radius inside central field
  double R_out = fabs(p / (q * 3e-4 * Bout)); // Radius in outer fields


  // Trajectory Calculation
  if (z[0] < zA1) { 
    // Before first Bout (straight line)
    x[0] = tan(thetax) * z[0];
  } 
  else if (z[0] == zA1) { 
    // Get the value of xA1
    xA1 = tan(thetax) * z[0];
  } 
  else if (z[0] <= zB2) { 
    // ===================== //
    // === Inside Fields === //
    // ===================== //
    if (z[0] <= zA) { 
      //std::cout<< "inside the fields, first field" << std::endl;
      // ----- First Bout Section -----
      double Delta = R_out * R_out - pow(z[0] - zA1 - R_out * sin(q * thetax), 2.0);
      if (Delta < 0) {
        x[0] = 0; 
      } else {
        x[0] = xA1 + q * (-R_out * cos(q * thetax) + sqrt(Delta));
      }
      if (z[0] == zA) { 
        xA = x[0]; 
      }
    }
    else if (z[0] <= zB) {
      //std::cout<< "inside the fields, second" << std::endl;
      // ----- Central B Field Section -----
      l1 = sqrt(std::pow(zA1 - zA, 2.0) + std::pow(xA1 - xA, 2.0));
      gamma1 = -acos(1 - (l1 * l1) / (2 * R_out * R_out));
      thetaxB1 = thetax + q * gamma1;

      double Delta = R_in * R_in - pow(z[0] - zA - R_in * sin(q * thetaxB1), 2.0);
      if (Delta < 0) {
        x[0] = 0; 
      } else {
        x[0] = xA + q * (-R_in * cos(q * thetaxB1) + sqrt(Delta));
      }
      if (z[0] == zB) { 
        xB = x[0]; 
      }
    }
    else {
      //std::cout<< "inside the fields, third" << std::endl;
      // ----- Second Bout Section -----
      l2 = sqrt(std::pow(zA - zB, 2.0) + std::pow(xA - xB, 2.0));
      gamma2 = -acos(1 - (l2 * l2) / (2 * R_in * R_in));
      thetaxB2 = thetaxB1 + q * gamma2;

      double Delta = R_out * R_out - pow(z[0] - zB - R_out * sin(q * thetaxB2), 2.0);
      if (Delta < 0) {
        x[0] = 0; 
      } else {
        x[0] = xB + q * (-R_out * cos(q * thetaxB2) + sqrt(Delta));
      }
      if (z[0] == zB2) { 
        xB2 = x[0]; 
      }
    }
  } 
  else { 
    //std::cout<< "out of all fields" << z[0] << std::endl;
    // After last Bout (straight line)
    l3 = sqrt(std::pow(zB - zB2, 2.0) + std::pow(xB - xB2, 2.0));
    gamma3 = -acos(1 - (l3 * l3) / (2 * R_out * R_out));
    thetaxB3 = thetaxB2 + q * gamma3;

    x[0] = tan(thetaxB3) * (z[0] - zB2) + xB2;
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
      x[0] = tan(thetax) * (z[0]);
  } 
  else if (z[0] < plane2) {
      // After first plane 
      x[0] = tan(theta_1) * z[0] - plane1*(tan(theta_1)-tan(thetax));
  } 
  else if (z[0] < zA) {
      double scatter = theta_2;
      x[0] = tan(scatter) * z[0] - plane2*(tan(scatter)-tan(theta_1))- plane1*(tan(theta_1)-tan(thetax));
  }
  else if (z[0] < plane3) {
      // After second plane
      double scatter = theta_2;
      if (z[0] == zA) { 
          xA = x[0];
	  //std::cout<<"xA="<<xA<<std::endl;
          std::cout << ": theta_1=" << theta_1 << ", theta_2=" << theta_2 << std::endl;
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
// Function to draw the trajectory with 3B
//
//******************************************************************************

void trajectory3B() {
  double Par[9] = {
    0.02,     // thetax (rad)
    1.32,    // p (MeV/c) - momentum
    -5000.0,   // zA (um) start of B field
    5000.0,   // zB (um) end of B field 
    1.0,      // q (charge) = +1
    0.1,       // B (T) magnetic field
    -8000.0,  //start Bout field
    8000.0,   //end Bout field
    0.05,     //Bout field
  };

  // Create graph for the trajectory
  TGraph *graph = new TGraph();

  double x;
  int pointIndex = 0;

  for (double z = -30000; z <= 20000; z += 100) {
    double zArray[1] = {z};
    sitritrajectorybended3B(zArray, &x, Par);
    graph->SetPoint(pointIndex, z, x);
    pointIndex++;
  }

  //Graph range
  graph->GetXaxis()->SetLimits(-30000, 20000); //max 4cm
  graph->GetYaxis()->SetRangeUser(-10000, 10000); //2cm hauteur pour voir les plans


  // Create canvas
  TCanvas *c1 = new TCanvas("c1", "Trajectory through Magnetic Field", 800, 600);

  c1->SetRightMargin(0.2); 

  graph->SetTitle("Particle trajectory through 3 magnetic fields;Z (um);X (um)");
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
  TLine *bFieldStart1 = new TLine(-8000, -planeHeight/2, -8000, planeHeight/2);
  bFieldStart1->SetLineColor(kMagenta);
  bFieldStart1->SetLineStyle(2);
  bFieldStart1->Draw("same");

  TLine *bFieldStart2 = new TLine(-5000, -planeHeight/2, -5000, planeHeight/2);
  bFieldStart2->SetLineColor(kMagenta);
  bFieldStart2->SetLineStyle(2);
  bFieldStart2->Draw("same");

  TLine *bFieldEnd1 = new TLine(8000, -planeHeight/2, 8000, planeHeight/2);
  bFieldEnd1->SetLineColor(kMagenta);
  bFieldEnd1->SetLineStyle(2);
  bFieldEnd1->Draw("same");

  TLine *bFieldEnd2 = new TLine(5000, -planeHeight/2, 5000, planeHeight/2);
  bFieldEnd2->SetLineColor(kMagenta);
  bFieldEnd2->SetLineStyle(2);
  bFieldEnd2->Draw("same");

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
// Function to draw the trajectory comparison between B fields
//
//******************************************************************************

void trajectory_comparison() {
  double Par[9] = {
    0.02,     // thetax (rad)
    1.32,    // p (MeV/c) - momentum
    -5000.0,   // zA (um) start of B field
    5000.0,   // zB (um) end of B field 
    1.0,      // q (charge) = +1
    0.1,       // B (T) magnetic field
    -8000.0,  // start Bout field
    8000.0,   // end Bout field
    0.05,     // Bout field
  };

  double planeHeight = 19872; // 19.872 mm in µm
  double planeThickness = 50; // 50 µm thickness

  // Plane positions
  double planeZ[4] = {
    -15000.0,        // Plane 1
    -10000.0,        // Plane 2
    10000.0,        // Plane 3
    15000.0         // Plane 4
  };
  
  // Create a MultiGraph to hold both graphs
  TMultiGraph *mg = new TMultiGraph();
  TLegend *legend = new TLegend(0.81, 0.65, 0.98, 0.85);

  // Create graphs for the two trajectories
  TGraph *graph1 = new TGraph(); // Trajectory 1
  TGraph *graph2 = new TGraph(); // Trajectory 2

  double x;
  int pointIndex = 0;

  // First trajectory: sitritrajectorybended3B
  for (double z = -30000; z <= 20000; z += 100) {
    double zArray[1] = {z};
    sitritrajectorybended3B(zArray, &x, Par);
    graph1->SetPoint(pointIndex++, z, x);
  }
  graph1->SetLineColor(kBlue);
  graph1->SetLineWidth(2);
  mg->Add(graph1);
  legend->AddEntry(graph1, "Trajectory Bended 3B", "l");

  // Reset point index
  pointIndex = 0;

  // Second trajectory: sitritrajectorybended
  for (double z = -30000; z <= 20000; z += 10) {
    double zArray[1] = {z};
    sitritrajectorybended(zArray, &x, Par);
    graph2->SetPoint(pointIndex++, z, x);
  }
  graph2->SetLineColor(kRed);
  graph2->SetLineWidth(2);
  mg->Add(graph2);
  legend->AddEntry(graph2, "Trajectory Bended", "l");

  // Create canvas and draw graph
  TCanvas *c1 = new TCanvas("c1", "Trajectory Comparison", 800, 600);
  mg->GetXaxis()->SetLimits(-30000, 20000);
  mg->GetYaxis()->SetRangeUser(-10000, 10000); //2cm hauteur pour voir les plans
  mg->SetTitle("Trajectory Comparison between both B field types; Z position (um); X position (um)");
  mg->Draw("AL");

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
  TLine *bFieldStart1 = new TLine(-8000, -planeHeight/2, -8000, planeHeight/2);
  bFieldStart1->SetLineColor(kMagenta);
  bFieldStart1->SetLineStyle(2);
  bFieldStart1->SetLineWidth(2);
  bFieldStart1->Draw("same");

  int deepMagenta = TColor::GetColor("#E0115F");

  TLine *bFieldStart2 = new TLine(-5000, -planeHeight/2, -5000, planeHeight/2);
  bFieldStart2->SetLineColor(deepMagenta);
  bFieldStart2->SetLineStyle(2);
  bFieldStart2->SetLineWidth(2);
  bFieldStart2->Draw("same");

  TLine *bFieldEnd1 = new TLine(8000, -planeHeight/2, 8000, planeHeight/2);
  bFieldEnd1->SetLineColor(kMagenta);
  bFieldEnd1->SetLineStyle(2);
  bFieldEnd1->SetLineWidth(2);
  bFieldEnd1->Draw("same");

  TLine *bFieldEnd2 = new TLine(5000, -planeHeight/2, 5000, planeHeight/2);
  bFieldEnd2->SetLineColor(deepMagenta);
  bFieldEnd2->SetLineStyle(2);
  bFieldEnd2->SetLineWidth(2);
  bFieldEnd2->Draw("same");

  // Magnetic Field Region
  TLine *legendBField1 = new TLine(0, 0, 1, 1);
  TLine *legendBField2 = new TLine(0, 0, 1, 1);
  legendBField1->SetLineColor(kMagenta);
  legendBField2->SetLineColor(deepMagenta);
  legendBField1->SetLineStyle(2);
  legendBField1->SetLineWidth(2);
  legendBField2->SetLineStyle(2);
  legendBField2->SetLineWidth(2);
  legend->AddEntry(legendBField1, "external B Field", "l");
  legend->AddEntry(legendBField2, "internal B Field", "l");

  // Sensing Planes
  TLine *legendPlane = new TLine(0, 0, 1, 1);
  legendPlane->SetLineColor(kBlack);
  legendPlane->SetLineStyle(1);
  legendPlane->SetLineWidth(2);
  legend->AddEntry(legendPlane, "Planes", "l");

  
  legend->Draw();
  c1->Update();
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
// Function to draw the comparison between a scattered trajectory and a simple one
//
//******************************************************************************

void simple_vs_scattered() {
  double Par[6] = {
    0.02,     // thetax (rad)
    1.32,    // p (MeV/c) - momentum
    -5000.0,   // zA (um) start of B field
    5000.0,   // zB (um) end of B field 
    1.0,      // q (charge) = +1
    0.1,       // B (T) magnetic field
  };

  double planeHeight = 19872; // 19.872 mm in µm
  double planeThickness = 50; // 50 µm thickness

  // Plane positions
  double planeZ[4] = {
    -15000.0,        // Plane 1
    -10000.0,        // Plane 2
    10000.0,        // Plane 3
    15000.0         // Plane 4
  };

  // Physical constants
  double m_e = 0.511;        // MeV (electron rest mass)
  double X0 = 21.82;         // g/cm^2 (radiation length for Silicon)
  double rho = 2.33;         // g/cm^3 (density of Silicon)
  double thickness = 50e-4;  // cm (50 µm thickness per plane)
  double E = sqrt(Par[1] * Par[1] + m_e * m_e); // Total energy
  double beta = Par[1] / E;                // Velocity fraction of c

  TRandom3 randGen;  // Random generator for scattering

  // Highland formula for Multiple Coulomb Scattering (RMS angle)
  double theta0 = (13.6 / (beta * Par[1])) * fabs(Par[4]) * sqrt(thickness / X0) * (1 + 0.038 * log(thickness / X0));
  
  // Create a MultiGraph to hold both graphs
  TMultiGraph *mg = new TMultiGraph();
  TLegend *legend = new TLegend(0.81, 0.65, 0.98, 0.85);

  // Create graphs for the two trajectories
  TGraph *graph1 = new TGraph(); // Trajectory 1

  std::vector<double> impact_values; //impact point
  std::vector<double> Run_number;

  double x;
  int pointIndex = 0;

  //Plotting mutliple runs to the compute the mean error
  for (double i = 0; i <= 100; i += 1) {

    bool localStopPlot = false; //to stop plotting if the particle won't go to the 3rd plane

    TGraph *graph = new TGraph();
    double x, y;
    int pointIndex = 0;

    randGen.SetSeed(time(0)+i);  // Unique seed based on current time
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
     
      // Compute slope for stopping condition
      if (pointIndex > 0) {
          double slope = (x - graph->GetY()[pointIndex - 1]) / (z - graph->GetX()[pointIndex - 1]);
          //if (b == 0.4 && z > 0 && z < 5000){std::cout << "Slope: " << slope << std::endl;}
          if (fabs(slope) > 16) //tested value
          {
              localStopPlot = true; // Stop plotting for this particle, it won't make it
              break;
          }
      }
     
      if (z == 10000 && x > -planeHeight/2){ //if impact in the plane 3
        impact_values.push_back(x);
        std::cout << "impact with plane 3 for run="<< i << " x="<< x << std::endl;
        Run_number.push_back(i);
      }
     


      graph->SetPoint(pointIndex, z, x);
      pointIndex++;
    }
 
    graph->SetTitle("Particle trajectory compared to the scattered trajectories;Z (um);X (um)");
    graph->SetLineWidth(2);
    
    mg->Add(graph);

    legend->AddEntry(graph, Form("Run = %.1f T", i), "l");

  }

  // Create canvas
  TCanvas *c1 = new TCanvas("c1", "Scattered Trajectory", 800, 600);
  
  c1->SetRightMargin(0.2); 
  mg->GetXaxis()->SetLimits(-30000, 20000);
  mg->GetYaxis()->SetRangeUser(-10000, 10000); //2cm hauteur pour voir les plans

  mg->Draw("A pmc plc"); //to have beautiful bird colors

  // Reset point index
  pointIndex = 0;
  double awaited_impact = 0;

  // Simple bended trajectory
  for (double z = -30000; z <= 20000; z += 100) {
    double zArray[1] = {z};
    sitritrajectorybended(zArray, &x, Par);
    graph1->SetPoint(pointIndex++, z, x);
    if (z == 10000 && x > -planeHeight/2){ //if impact in the plane 3
      awaited_impact = x;
      std::cout << "impact with plane 3 for simple trajectory at x= "<< awaited_impact << std::endl;}
  }
  graph1->SetLineColor(kBlue);
  graph1->SetLineWidth(2);
  mg->Add(graph1);
  legend->AddEntry(graph1, "Trajectory Bended", "l");

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
  bFieldStart->SetLineWidth(2);
  bFieldStart->Draw("same");

  TLine *bFieldEnd = new TLine(5000, -planeHeight/2, 5000, planeHeight/2);
  bFieldEnd->SetLineColor(kMagenta);
  bFieldEnd->SetLineStyle(2);
  bFieldEnd->SetLineWidth(2);
  bFieldEnd->Draw("same");

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

  
  legend->Draw();
  c1->Update();

  double sum_abs_error = 0.0;
  double sum_squared_error = 0.0;
  int N = impact_values.size();

  for (double impact : impact_values) {
      double error = impact - awaited_impact;
      sum_abs_error += std::abs(error);    // Absolute error for MAE
      sum_squared_error += error * error;  // Squared error for RMSE
  }

  double MAE = sum_abs_error / N;
  double RMSE = std::sqrt(sum_squared_error / N);

  std::cout << "Mean Absolute Error (MAE): " << MAE << " µm" << std::endl;
  std::cout << "Root Mean Square Error (RMSE): " << RMSE << " µm" << std::endl;

}


//********************************************************************************
//
// Function to see different B fields and error with the simple config
//
//********************************************************************************

void trajectory_error() {
  double Par[6] = {
    0.02,     // thetax (rad)
    1.32,    // p (MeV/c) - momentum
    -5000.0,   // zA (um) start of B field
    5000.0,   // zB (um) end of B field 
    1.0,      // q (charge) = +1
    0.1       // B (T) magnetic field
  };
  //Planes 

  double planeHeight = 19872; // 19.872 mm in µm
  double planeThickness = 50; // 50 µm thickness

  // Plane positions
  double planeZ[4] = {
    -15000.0,        // Plane 1
    -10000.0,        // Plane 2
    10000.0,        // Plane 3
    15000.0         // Plane 4
  };
  

  TMultiGraph *mg = new TMultiGraph();
  TLegend *legend = new TLegend(0.81, 0.65, 0.98, 0.85);
  
  std::vector<double> Bfield_values;
  std::vector<double> impact_values; //impact point
  std::vector<double> error_dist; //to see error on impact point
  
  //changing magnetic field 
  for (double b = 0.1; b <= 0.5; b += 0.05) {
    Par[5] = b;

    bool localStopPlot = false; //to stop plotting if the particle won't go to the 3rd plane

    TGraph *graph = new TGraph();
    double x, y;
    int pointIndex = 0;

    for (double z = -30000; z <= 20000; z += 10) { // 5 cm total
      if (localStopPlot) break;

      double zArray[1] = {z};
      sitritrajectorybended(zArray, &x, Par);
     
      // Compute slope for stopping condition
      if (pointIndex > 0) {
          double slope = (x - graph->GetY()[pointIndex - 1]) / (z - graph->GetX()[pointIndex - 1]);
          //if (b == 0.4 && z > 0 && z < 5000){std::cout << "Slope: " << slope << std::endl;}
          if (fabs(slope) > 16) //tested value
          {
              localStopPlot = true; // Stop plotting for this particle, it won't make it
              break;
          }
      }
     
      if (z == 10000 && x > -planeHeight/2){ //if impact in the plane 3
        impact_values.push_back(x);
        std::cout << "impact with plane 3 for B="<< b << " T: x="<< x << std::endl;
        Bfield_values.push_back(b);
      }
     


      graph->SetPoint(pointIndex, z, x);
      //graph_y->SetPoint(pointIndex, z, y);
      pointIndex++;
    }
 
    graph->SetTitle("Particle trajectories for different B values;Z (um);X (um)");
    //graph->SetLineColor(colors[colorIndex % 7]);
    graph->SetLineWidth(2);
    
    mg->Add(graph);

    legend->AddEntry(graph, Form("B = %.1f T", b), "l");

  }

  // Create canvas
  TCanvas *c1 = new TCanvas("c1", "Trajectory through Magnetic Field", 800, 600);
  
  c1->SetRightMargin(0.2); 
  mg->GetXaxis()->SetLimits(-30000, 20000);
  mg->GetYaxis()->SetRangeUser(-10000, 10000); //2cm hauteur pour voir les plans

  mg->Draw("A pmc plc"); //to have beautiful bird colors

  //add lines and bfield limit


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
  TLine *bFieldStart = new TLine(-5000, -7000, -5000, 7000);
  bFieldStart->SetLineColor(kMagenta);
  bFieldStart->SetLineStyle(2);
  bFieldStart->Draw("same");

  TLine *bFieldEnd = new TLine(5000, -7000, 5000, 7000);
  bFieldEnd->SetLineColor(kMagenta);
  bFieldEnd->SetLineStyle(2);
  bFieldEnd->Draw("same");

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
  
  // Create a Canvas
  TCanvas *c2 = new TCanvas("c2", "Error on position for different fields", 800, 600);

  for (size_t i = 0; i < impact_values.size()-1; i++){
    double dist = fabs(impact_values[i+1] - impact_values[i]);
    error_dist.push_back(dist);
    std::cout << "error distance of "<< dist <<"µm if we use wrong Bfield"<<std::endl;
  }

  // in case there is just an upper or lower limit, set it as global error FOR NOW

  std::cout << "for B=" << Bfield_values[0] << " the error is of +" << error_dist[0] << " µm and -" << error_dist[0] << " µm"<< std::endl;

  for (size_t j = 1; j < Bfield_values.size()-1; j++){
    std::cout << "for B=" << Bfield_values[j] << " the error is of +" << error_dist[j] << " µm and -" << error_dist[j-1] << " µm"<< std::endl;}

  std::cout << "for B=" << Bfield_values[Bfield_values.size() - 1] << " the error is of +" << error_dist[error_dist.size() - 1] << " µm and -" << error_dist[error_dist.size() - 1]<< " µm"<< std::endl;

  // Create Graph with Asymmetric Errors
  int n = Bfield_values.size();
  TGraphAsymmErrors *graph2 = new TGraphAsymmErrors(n);
    
  for (size_t i = 0; i < n; i++) {
      graph2->SetPoint(i, Bfield_values[i], impact_values[i]);
      graph2->SetPointError(i, 0, 0, error_dist[i], error_dist[i-1]); // X error is 0, Y error is asymmetric
  }

  // Customize graph appearance
  graph2->SetMarkerStyle(21);
  graph2->SetMarkerSize(1.2);
  graph2->SetMarkerColor(kRed);
  graph2->SetLineColor(kBlue);
  graph2->SetLineWidth(2);
  graph2->SetTitle("Impact Position vs Magnetic Field; Magnetic Field B (T); Impact Position (um)");
    
  // Draw Graph
  graph2->Draw("AP"); // "A" for axes, "P" for points with error bars

  // Add Legend
  TLegend *legend2 = new TLegend(0.81, 0.75, 0.98, 0.88);
  legend2->AddEntry(graph2, "Impact Position with Error", "p");
  legend2->Draw();
  
  c2->SetGrid();
  c2->Update();

}

//********************************************************************************
//
// Function to see different B fields and error with the 3B config
//
//********************************************************************************

void trajectory_error3B() {
  double Par[9] = {
    0.02,     // thetax (rad)
    1.32,    // p (MeV/c) - momentum
    -5000.0,   // zA (um) start of B field
    5000.0,   // zB (um) end of B field 
    1.0,      // q (charge) = +1
    0.1,       // B (T) magnetic field
    -8000.0,  //start Bout field
    8000.0,   //end Bout field
    0.05,     //Bout field
  };
  //Planes 

  double planeHeight = 19872; // 19.872 mm in µm
  double planeThickness = 50; // 50 µm thickness

  // Plane positions
  double planeZ[4] = {
    -15000.0,        // Plane 1
    -10000.0,        // Plane 2
    10000.0,        // Plane 3
    15000.0         // Plane 4
  };
  

  TMultiGraph *mg = new TMultiGraph();
  TLegend *legend = new TLegend(0.81, 0.65, 0.98, 0.85);
  
  std::vector<double> Bfield_values;
  std::vector<double> impact_values; //impact point
  std::vector<double> error_dist; //to see error on impact point
  
  //changing magnetic field 
  for (double b = 0.1; b <= 0.5; b += 0.05) {
    Par[5] = b;

    bool localStopPlot = false; //to stop plotting if the particle won't go to the 3rd plane

    TGraph *graph = new TGraph();
    double x, y;
    int pointIndex = 0;

    for (double z = -30000; z <= 20000; z += 10) { // 5 cm total
      if (localStopPlot) break;

      double zArray[1] = {z};
      sitritrajectorybended3B(zArray, &x, Par);
     
      // Compute slope for stopping condition
      if (pointIndex > 0) {
          double slope = (x - graph->GetY()[pointIndex - 1]) / (z - graph->GetX()[pointIndex - 1]);
          //if (b == 0.4 && z > 0 && z < 8000){std::cout << "Slope: " << slope << std::endl;}
          if (fabs(slope) > 16) //tested value
          {
              localStopPlot = true; // Stop plotting for this particle, it won't make it
              break;
          }
      }
     
      if (z == 10000 && x > -planeHeight/2){ //if impact in the plane 3
        impact_values.push_back(x);
        std::cout << "impact with plane 3 for B="<< b << " T: x="<< x << std::endl;
        Bfield_values.push_back(b);
      }
     


      graph->SetPoint(pointIndex, z, x);
      //graph_y->SetPoint(pointIndex, z, y);
      pointIndex++;
    }
 
    graph->SetTitle("Particle trajectories for different B values;Z (um);X (um)");
    //graph->SetLineColor(colors[colorIndex % 7]);
    graph->SetLineWidth(2);
    
    mg->Add(graph);

    legend->AddEntry(graph, Form("B = %.1f T", b), "l");

  }

  // Create canvas
  TCanvas *c1 = new TCanvas("c1", "Trajectory through Magnetic Fields", 800, 600);
  
  c1->SetRightMargin(0.2); 
  mg->GetXaxis()->SetLimits(-30000, 20000);
  mg->GetYaxis()->SetRangeUser(-10000, 10000); //2cm hauteur pour voir les plans

  mg->Draw("A pmc plc"); //to have beautiful bird colors

  //add lines and bfield limit


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
  TLine *bFieldStart1 = new TLine(-8000, -planeHeight/2, -8000, planeHeight/2);
  bFieldStart1->SetLineColor(kMagenta);
  bFieldStart1->SetLineStyle(2);
  bFieldStart1->Draw("same");

  TLine *bFieldStart2 = new TLine(-5000, -planeHeight/2, -5000, planeHeight/2);
  bFieldStart2->SetLineColor(kMagenta);
  bFieldStart2->SetLineStyle(2);
  bFieldStart2->Draw("same");

  TLine *bFieldEnd1 = new TLine(8000, -planeHeight/2, 8000, planeHeight/2);
  bFieldEnd1->SetLineColor(kMagenta);
  bFieldEnd1->SetLineStyle(2);
  bFieldEnd1->Draw("same");

  TLine *bFieldEnd2 = new TLine(5000, -planeHeight/2, 5000, planeHeight/2);
  bFieldEnd2->SetLineColor(kMagenta);
  bFieldEnd2->SetLineStyle(2);
  bFieldEnd2->Draw("same");

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
  
  // Create a Canvas
  TCanvas *c2 = new TCanvas("c2", "Error on position for different fields", 800, 600);

  for (size_t i = 0; i < impact_values.size()-1; i++){
    double dist = fabs(impact_values[i+1] - impact_values[i]);
    error_dist.push_back(dist);
    std::cout << "error distance of "<< dist <<" µm if we use wrong Bfield"<<std::endl;
  }

  // in case there is just an upper or lower limit, set it as global error FOR NOW

  std::cout << "for B=" << Bfield_values[0] << " the error is of +" << error_dist[0] << " µm and -" << error_dist[0] << " µm"<< std::endl;

  for (size_t j = 1; j < Bfield_values.size()-1; j++){
    std::cout << "for B=" << Bfield_values[j] << " the error is of +" << error_dist[j] << " µm and -" << error_dist[j-1] << " µm"<< std::endl;}

  std::cout << "for B=" << Bfield_values[Bfield_values.size() - 1] << " the error is of +" << error_dist[error_dist.size() - 1] << " µm and -" << error_dist[error_dist.size() - 1]<< " µm"<< std::endl;

  // Create Graph with Asymmetric Errors
  int n = Bfield_values.size();
  TGraphAsymmErrors *graph2 = new TGraphAsymmErrors(n);
    
  for (size_t i = 0; i < n; i++) {
      graph2->SetPoint(i, Bfield_values[i], impact_values[i]);
      graph2->SetPointError(i, 0, 0, error_dist[i], error_dist[i-1]); // X error is 0, Y error is asymmetric
  }

  // Customize graph appearance
  graph2->SetMarkerStyle(21);
  graph2->SetMarkerSize(1.2);
  graph2->SetMarkerColor(kRed);
  graph2->SetLineColor(kBlue);
  graph2->SetLineWidth(2);
  graph2->SetTitle("Impact Position vs Magnetic Field; Magnetic Field B (T); Impact Position (um)");
    
  // Draw Graph
  graph2->Draw("AP"); // "A" for axes, "P" for points with error bars

  // Add Legend
  TLegend *legend2 = new TLegend(0.81, 0.75, 0.98, 0.88);
  legend2->AddEntry(graph2, "Impact Position with Error", "p");
  legend2->Draw();
  
  c2->SetGrid();
  c2->Update();

}