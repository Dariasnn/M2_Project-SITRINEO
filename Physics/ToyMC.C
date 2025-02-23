/*
Toy MonteCarlo of the SITRNEO tracker

Simulate the momentum distribution obtained with SITRINEO
  taking into account :
  - Fermi distribution of particle momenta emmitted by the beta source
  - a given Landau uncertainty
  - Comparaison of the Fermi distribution and the experimental distribution to have the real function distribution of the total uncertainty.

This is a ROOT macro, to be used with the ROOT package
* Compile in ROOT with with:
ROOT> .L ToyMC.C+

* Usage in ROOT:
ROOT> ToyMC(Z-of-source, beta-end-spectrum(keV), #events)

* Created by Romain Schotter & Vincent Juste, 2020
*Update by Rachel TECHI and Daria SENINA 22/02/2025
*/

#include <TRandom3.h>
#include"TF1Convolution.h"
#include"TMath.h"
#include"TH1F.h"
#include"TH3D.h"
#include"TCanvas.h"
#include"TRandom.h"
#include"TStyle.h"
#include"TGraph2D.h"
#include"TSystem.h"
#include"TFile.h"
#include <TLegend.h>
#include <TF1Convolution.h> 
#include <cmath>



using namespace std ;


//
// Display Beta spectrum
//


TH1F* extraction() {
//This function allow to extract the desired historgram from a root file

  TFile *in_file = new TFile("sitrineo_run1210.root", "READ");  // Open the root file
  if (!in_file || in_file->IsZombie()) {
      std::cerr << "Error: root file impossible to open !" << std::endl;
      return 0;
  }

  //extract the histogramme (check if it is a TH1F, TH1D, etc.)
  TH1F *histo = (TH1F*) in_file->Get("hTrackPairMomentum");
  if (!histo) {
      std::cerr << "Erreur : Histogramme not found !" << std::endl;
      in_file->Close();
      return 0;
  }

  // Creat a canva to plot the histogramme
  //TCanvas *canvas = new TCanvas("canvas", "Spectre en impulsion", 800, 600);
  //histo->Draw();

  // Hold on the plot
  //canvas->Update();
  return histo; //allow to send the histo to the function ToyMC.
}

TH1F* g_histo = nullptr;  // global variable

Double_t fitFunction(Double_t* x, Double_t* p) {
  return p[0] * g_histo->Interpolate(x[0]); //interpolation of the histogram
}

TF1* HistToFunc(TH1F* histo, const char* name) {
  g_histo = histo; 
  TF1* func = new TF1(name, fitFunction, 
                      histo->GetXaxis()->GetXmin(), 
                      histo->GetXaxis()->GetXmax(), 
                      1);
  func->SetParameter(0, 1.0); 
  return func;
}

// ---------------------------------------------------------
double EtoP(double E) {return sqrt(E*E - 511*511);} //function which allow to convert energy to momentum by using electron mass (511 keV)


// ---------------------------------------------------------
Double_t fermiFunction( Double_t *x, Double_t *par ) {
// Implementation of the Fermi function for Beta- spectrum

  Double_t Q = par[0]; // spectrum end
  Double_t Z = par[1]; // atomic number
  Double_t Norm = par[2];
  Double_t Cut_off = par[3] ; //energy cut off
  Double_t E = x[0]; // Beta energy

  if( E > Cut_off ){
    // constants
    Double_t alpha = 1./137.;
    Double_t eMass = 511; // keV/c2

    Double_t F1 = sqrt(E*E+2*eMass)*(Q-E)*(Q-E)*(E+eMass);
    Double_t g = TMath::TwoPi()*alpha*Z*(1+E/Q)*sqrt(E*E/Q/Q+2*E/Q);
    Double_t F2 = Z*Z*(E/Q+1)*(E/Q+1)*alpha*alpha+E*E/Q/Q+E/2/Q;

    Double_t p = F1*g/(exp(g)-1)*TMath::Power( F2, sqrt(1-Z*Z*alpha*alpha)-1);

    //printf("\n E=%.1f, F1=%.2f, F2=%.2f, g=%.2f\n", E, F1, F2, g);
    //printf("    g/(1-e(g))=%e, power=%.2f\n", g/(1-exp(g)), sqrt(1-Z*Z*alpha*alpha)-1);

    return Norm*p/1.e11;

  }

  else{
      return 0 ;
  }
}


// ---------------------------------------------------------

Double_t randomFermi( double Q, int Z, double cut_off) {
  // Draw an energy randomly distributed between 0 and Q
  //  according to the Fermi function with Z=atomic number.
  //  energies below cut-off are not generated (they don't go through SITRINEO)

  double par[4];
  par[0] = Q;
  par[1] = Z;
  par[2] = 1.;
  par[3] = cut_off ;
  double E[1] = { gRandom->Uniform()*Q };
  double max = 15.; // 15 for Sr90 6 for O15, 0.9 for C11, 0.2 for F18,
  while( max*gRandom->Uniform()>fermiFunction( E, par) ) {
    //printf( " essai avec E=%.2f\n", E[0]);
    E[0] = gRandom->Uniform()*Q;
  }

  return E[0];
}

// ---------------------------------------------------------

//Update by R.TECHI 18/02/2025


void ToyMC( Double_t Z = 39., Double_t Q = 2280, Int_t Nevts=1000 ){
  // Build the beta energy spectrum for a dedicated nuclei.
  // Fit the distribution with a simple function.
  //
  // Energy unit in MeV
  //
  // JB 2013/02/28
  //
  // Nuclei caracteristics
  // Beta-: H3 18, C14 156, S32 167, P33 250, Ca45 257, P32 1710, Sr39 2280
  // Beta+: C11 961, N13 1119, O15 1720, F18 650


  // init and reset of Root
   //gROOT->Reset();
   gStyle->SetOptFit(1); //show the parameters of the fit

  Char_t text[100];
  sprintf( text, "Z=%d, Q=%.0f keV, N = %d", (int)Z, Q, Nevts);


  //1D histos
  TH1F *hspectrum_F = new TH1F( "hspectrum_F", "Fermi spectrum", 100, 0., 1.2*Q);
  hspectrum_F->SetXTitle("Energy (keV)");
  hspectrum_F->SetFillStyle(1001);
  hspectrum_F->SetFillColor(3);
  hspectrum_F->SetStats(0);

  TH1F *hspectrum_G = new TH1F( "hspectrum_G", "Landau spectrum", 100, -Q, Q);
  hspectrum_G->SetXTitle("Energy (keV)");
  hspectrum_G->SetFillStyle(1001);
  hspectrum_G->SetFillColor(3);
  hspectrum_G->SetStats(0);

  TH1F *hspectrum_FG = new TH1F( "hspectrum_FG", "Fermi & Cauchy Convolution", 100, 0., 1e4); //final spectrum 
  hspectrum_FG->SetXTitle("Energy (keV)");
  hspectrum_FG->SetFillStyle(1001);
  hspectrum_FG->SetFillColor(3);
  hspectrum_FG->SetStats(0);

  TH1F *histo_exp=extraction(); // Get the experimental histo take by the acquision

  double Cutoff_initial = 100 ; // cutoff is the minimal energy to go through SITRINEO (keV)
  double Cutoff_final = 1800 ;
  int NCuts = 40 ;

  double sigma_initial = 0.1 ;
  double sigma_final = 0.4 ;
  int NSigma = 30 ;


  //contain the value of the MPV or mean or rms or FWHM or skewness or kurtosis as a function of the cutoff and sigma,
  // in order to find "good" observables to identify the cutoff and the sigma ( sigma = sigma of the Landau )
  //in the end, good observable for the cutoff is the mean because it is almost constant with sigma.

  vector<double> Sigma, cutoff, mean, mean_cutoff, MPV, MPV_cutoff, rms, fwhm, skew, kurt ;
  //search to get mean, MPV, rms, FWHM, skewness, kurtosis for different cutoff and sigma

  for( int iCut = 0 ; iCut < NCuts ; iCut++ ){//loop on the cutoff
      //if( iCut%10 == 0) cout << "Cut " << iCut << endl;

      double cut_off = Cutoff_initial + iCut*(Cutoff_final - Cutoff_initial)/NCuts ;

      for( int iSigma = 0 ; iSigma < NSigma ; iSigma++ ){//loop on sigma

            double sigma = sigma_initial + iSigma*(sigma_final - sigma_initial)/NSigma ;

            for ( int i=0; i<Nevts; i++) {//loop on the number of events
                double_t x=randomFermi(Q,(int)Z, cut_off);
                double_t y=gRandom->Landau(0.,sigma*x);
                double_t z = EtoP(x)+y ;
                hspectrum_FG->Fill(z);

            }//end loop on the number of events

            
            hspectrum_FG->Scale(1./Nevts);
            

            //find Full Width at Half Maxima
            double low_edge = 0 ; double upper_edge = 0 ; double FWHM = 0 ;
            for( int i = 0 ; i < hspectrum_FG->GetNbinsX() ; i++){
                if( i < hspectrum_FG->GetMaximumBin() ){
                    if( hspectrum_FG->GetBinContent(i) < (hspectrum_FG->GetBinContent( hspectrum_FG->GetMaximumBin() ))/2. ){
                        low_edge = hspectrum_FG->GetBinCenter(i);
                    }
                    else{ continue ;}
                }
                if( i > hspectrum_FG->GetMaximumBin() ){
                    if( hspectrum_FG->GetBinContent(i) > (hspectrum_FG->GetBinContent( hspectrum_FG->GetMaximumBin() ))/2. ){
                        upper_edge = hspectrum_FG->GetBinCenter(i);
                    }
                    else{ continue ;}
                }
            }

            FWHM = upper_edge-low_edge;

            cutoff.push_back(cut_off);
            Sigma.push_back(sigma);
            mean.push_back(hspectrum_FG->GetMean());
            mean_cutoff.push_back(hspectrum_FG->GetMean() -cut_off);
            MPV.push_back(hspectrum_FG->GetBinCenter( hspectrum_FG->GetMaximumBin()));
            MPV_cutoff.push_back(hspectrum_FG->GetBinCenter( hspectrum_FG->GetMaximumBin()) - cut_off );
            rms.push_back(hspectrum_FG->GetRMS());
            fwhm.push_back(FWHM);
            skew.push_back(hspectrum_FG->GetSkewness());
            kurt.push_back(hspectrum_FG->GetKurtosis());

            hspectrum_FG->Reset();


      }//end loop on sigma


  }//end loop on cutoff


  TGraph2D* hMean = new TGraph2D(Sigma.size(), Sigma.data(), cutoff.data(), mean.data());
  hMean->SetTitle("Mean(Sigma, Cutoff); Sigma ; Cutoff (keV); Mean (keV)");
  hMean->GetXaxis()->SetTitleOffset(0.05);
  hMean->GetYaxis()->SetTitleOffset(0.05);
  hMean->GetZaxis()->SetTitleOffset(0.05);

  TGraph2D* hMean_cutoff = new TGraph2D(Sigma.size(), Sigma.data(), cutoff.data(), mean_cutoff.data());
  hMean_cutoff->SetTitle("Mean(Sigma, Cutoff) - Cutoff; Sigma ; Cutoff (keV); Mean -cutoff (keV)");
  hMean_cutoff->GetXaxis()->SetTitleOffset(0.05);
  hMean_cutoff->GetYaxis()->SetTitleOffset(0.05);
  hMean_cutoff->GetZaxis()->SetTitleOffset(0.05);
  hMean_cutoff->SetMarkerStyle(kOpenSquare);

  TGraph2D* hMPV = new TGraph2D(Sigma.size(), Sigma.data(), cutoff.data(), MPV.data());
  hMPV->SetTitle("MPV(Sigma, Cutoff); Sigma ; Cutoff (keV); MPV (keV)");
  hMPV->GetXaxis()->SetTitleOffset(0.05);
  hMPV->GetYaxis()->SetTitleOffset(0.05);
  hMPV->GetZaxis()->SetTitleOffset(0.05);
  hMPV->SetMarkerStyle(kOpenSquare);

  TGraph2D* hMPV_cutoff = new TGraph2D(Sigma.size(), Sigma.data(), cutoff.data(), MPV_cutoff.data());
  hMPV_cutoff->SetTitle("MPV(Sigma, Cutoff) - Cutoff; Sigma ; Cutoff (keV); MPV -cutoff (keV)");
  hMPV_cutoff->GetXaxis()->SetTitleOffset(0.05);
  hMPV_cutoff->GetYaxis()->SetTitleOffset(0.05);
  hMPV_cutoff->GetZaxis()->SetTitleOffset(0.05);
  hMPV_cutoff->SetMarkerStyle(kOpenSquare);

  TGraph2D* hRMS = new TGraph2D(Sigma.size(), Sigma.data(), cutoff.data(), rms.data());
  hRMS->SetTitle("RMS(Sigma, Cutoff); Sigma (keV); Cutoff (keV); RMS (keV)");
  hRMS->GetXaxis()->SetTitleOffset(0.05);
  hRMS->GetXaxis()->SetLabelSize(0.025);
  hRMS->GetYaxis()->SetTitleOffset(0.05);
  hRMS->GetZaxis()->SetTitleOffset(0.05);

  TGraph2D* hFWHM = new TGraph2D(Sigma.size(), Sigma.data(), cutoff.data(), fwhm.data());
  hFWHM->SetTitle("FWHM(Sigma, Cutoff); Sigma (keV); Cutoff (keV); FWHM (keV)");
  hFWHM->GetXaxis()->SetTitleOffset(0.05);
  hFWHM->GetXaxis()->SetLabelSize(0.025);
  hFWHM->GetYaxis()->SetTitleOffset(0.05);
  hFWHM->GetZaxis()->SetTitleOffset(0.05);
  hFWHM->SetMarkerStyle(kOpenSquare);

  TGraph2D* hSKEW = new TGraph2D(Sigma.size(), Sigma.data(), cutoff.data(), skew.data());
  hSKEW->SetTitle("Skew(Sigma, Cutoff); Sigma (keV); Cutoff (keV); Skewness");
  hSKEW->GetXaxis()->SetTitleOffset(0.05);
  hSKEW->GetXaxis()->SetLabelSize(0.025);
  hSKEW->GetYaxis()->SetTitleOffset(0.05);
  hSKEW->GetZaxis()->SetTitleOffset(0.05);
  hSKEW->SetMarkerStyle(kOpenSquare);

  TGraph2D* hKURT = new TGraph2D(Sigma.size(), Sigma.data(), cutoff.data(), kurt.data());
  hKURT->SetTitle("Kurt(Sigma, Cutoff); Sigma (keV); Cutoff (keV); Kurtosis");
  hKURT->GetXaxis()->SetTitleOffset(0.05);
  hKURT->GetXaxis()->SetLabelSize(0.025);
  hKURT->GetYaxis()->SetTitleOffset(0.05);
  hKURT->GetZaxis()->SetTitleOffset(0.05);
  hKURT->SetMarkerStyle(kOpenSquare);

  //Generates the spectrum, resulting from the convolution of the fermi function and an uncertainty distribution function, for a given parameters and cutoff
  double cut_off = 800; //keV
  double sigma = 0.01;
  double gamma=0.01;
  double y0=0;

  TRandom3 randGen;
   for ( int i=0; i<Nevts; i++) {
    double_t x=randomFermi(Q,(int)Z, cut_off);
    double_t y=gRandom->Landau(0.,sigma*x); //if we consider a landau uncertainty

    double_t uniform = randGen.Uniform(0, 1);
    double_t y2 = y0 + gamma * tan(M_PI * (uniform - 0.5)); // Transformation inverse for a cauchy uncertainty

    double_t z = EtoP(x+y2);
    hspectrum_F->Fill(x);
    hspectrum_G->Fill(y2);
    hspectrum_FG->Fill(z);//convolute histogram

  }

  hspectrum_F->Scale(1./Nevts);
  hspectrum_G->Scale(1./Nevts);
  hspectrum_FG->Scale(1./Nevts);

  //The left part of the obtained spectrum is dominated by the Landau, so fitting the left side of the obtained spectrum with a Landau allows us to access to sigma. Thanks to sigma and the mean value, we can deduce the cutoff
  //TF1 *fit = new TF1("fit", "landau", 0, hspectrum_FG->GetMean());
  //hspectrum_FG->Fit("fit", "RN"); //allow to have information about the Khi2 test (compare hspectrum_FG and the landau distribution) and also the caracteristic of the the fit
  //double sigma_fit = fit->GetParameter(2)/fit->GetParameter(1) ;
  //std::cout << "Sigma fit relative = " << sigma_fit << std::endl ;
  //std::cout << "Sigma fit = " << fit->GetParameter(2) << std::endl ;
  //std::cout << "Mean value = " << hspectrum_FG->GetMean() << " ==> Cut off = "  << std::endl ;
  //std::cout << "Mean value Fermi = " << hspectrum_F->GetMean() << std::endl ;
  //std::cout << "Std. deviation value Landau = " << hspectrum_G->GetRMS() << std::endl ;


  double maxi_FG = hspectrum_FG->GetMaximum();
  double maxi_exp = histo_exp->GetMaximum();
  hspectrum_FG->Scale(maxi_exp/maxi_FG); // Normalize 

  TF1* func_FG = HistToFunc(hspectrum_FG, "func_FG");
  histo_exp->Fit(func_FG, "RN");//fit between MC and exp (probably doesn't work)

  //canvas with all observables
  TCanvas *c1 = new TCanvas("Cutoff", text, 0,0, 1800, 900);
  c1->Divide(4,2);
  c1->cd(1);
  hMean->Draw("surf1");
  c1->cd(2);
  hMean_cutoff->Draw("surf1");
  c1->cd(3);
  hMPV->Draw("surf1");
  c1->cd(4);
  hMPV_cutoff->Draw("surf1");
  c1->cd(5);
  hRMS->Draw("surf1");
  c1->cd(6);
  hFWHM->Draw("surf1");
  c1->cd(7);
  hSKEW->Draw("surf1");
  c1->cd(8);
  hKURT->Draw("surf1");
  c1->Update();

  //canvas with all observables
  TCanvas *c2 = new TCanvas("Theoretical_Spectrum", text, 0,0, 1600, 1000); 
  double max_FG = hspectrum_FG->GetMaximum();
  double max_exp = histo_exp->GetMaximum();
  double maxAll = std::max(max_exp, max_FG); 
  hspectrum_FG->Scale(max_exp/max_FG); // Normalize 
  hspectrum_FG->SetFillColor(0); 

  hspectrum_FG->Draw("hist");
  histo_exp->Draw("hist same");
  double xmin = hspectrum_FG->GetXaxis()->GetXmin();
  double xmax = hspectrum_FG->GetXaxis()->GetXmax();
  hspectrum_FG->GetXaxis()->SetLimits(xmin/1000,xmax/1000);// in MeV now

  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(hspectrum_FG, "Theoretical Spectrum", "l");
  legend->AddEntry(histo_exp, "Experimental Spectrum", "l");
  legend->Draw();

  hspectrum_FG->GetXaxis()->SetTitle("Momentum (MeV/c)");
  //TFile fRoot("~/Work/M2_2022/sitrineo_2mil_chi2_cuts_10_60.root", "UPDATE");
  //fRoot.Delete("hspectrum_FG;1");
  hspectrum_FG->Write();

  hspectrum_FG->SetLineColor(kRed);
  histo_exp->SetLineColor(kBlue);
  //fRoot.Close();

}

// Now the goal is to use the experimental histogram to try to find 
// a distribution function related to the uncertaincy. We can get it 
// thanks to the iterative Richardson-Lucy algorithm .

//Update by R.TECHI 18/02/2025

// Improved version of the deconvolution analysis
// The goal is to extract the uncertainty distribution by deconvolving
// the experimental spectrum with the theoretical Fermi distribution
//convertion of the histogram to function

TH1F* ApproximateDeconv(TH1F* histo_exp, TH1F* hspectrum_F) { //probably doesn't work
  // Copy experimental histogram to stock the result
  TH1F* histo_deconv = (TH1F*)hspectrum_F->Clone("histo_deconv");
  histo_deconv->Reset();  // Clean the deconvolution histogram

  int nBins = histo_exp->GetNbinsX();

  hspectrum_F->Scale(1.0 / hspectrum_F->Integral()); //normalize (area=1)

  // Loop on each bins to correct the contribution from hspectrum
    for (int i = 1; i <= nBins; ++i) {
      double exp_content = histo_exp->GetBinContent(i);

      // Substraction 
      double correction = 0.0;
      for (int j = 1; j <= nBins; ++j) {
          double kernel = hspectrum_F->GetBinContent(j);

        // shuft around the central bin
        int target_bin = i - j + nBins / 2; 
        if (target_bin >= 1 && target_bin <= nBins) {
          correction += kernel * histo_exp->GetBinContent(target_bin);
          }
      }

      // Correction (only positive value)
      double new_content = exp_content - correction;
      if (new_content < 0) new_content = 0;
      histo_deconv->SetBinContent(i, new_content);
  }

  return histo_deconv;
}

TH1F* RichardsonLucyDeconv(TH1F* histo_exp, TH1F* hspectrum_F) {
  TH1F* estimate = (TH1F*)histo_exp->Clone("estimate");
  int nBins = histo_exp->GetNbinsX();
  const int maxIterations = 10;  // Number of iteration
  double epsilon = 1e-10;        // no division by 0 !
  
  // Normalize the kernel (Fermi spectrum)
  hspectrum_F->Scale(1.0 / hspectrum_F->Integral());
  
  // Richardson-Lucy iteration
  for (int iter = 0; iter < maxIterations; iter++) {
      // Compute the current estimation convolution 
      TH1F* convolved = (TH1F*)estimate->Clone("convolved");
      convolved->Reset();
      
      for (int i = 1; i <= nBins; i++) {
          double sum = 0;
          for (int j = 1; j <= nBins; j++) {
              int shift = i - j + nBins/2;
              if (shift >= 1 && shift <= nBins) {
                  sum += hspectrum_F->GetBinContent(j) * estimate->GetBinContent(shift);
              }
          }
          convolved->SetBinContent(i, sum + epsilon);
      }
      
      // Compute the correction factor
      for (int i = 1; i <= nBins; i++) {
          double correction = histo_exp->GetBinContent(i) / convolved->GetBinContent(i);
          // Update the estimation
          double new_value = estimate->GetBinContent(i) * correction;
          estimate->SetBinContent(i, new_value);
      }
      
      delete convolved;
  }
  
  return estimate;
}


void ToyMC_deconv(Double_t Z = 39., Double_t Q = 2280, Int_t Nevts=1000){
  TCanvas *c= new TCanvas("Spectrum", "Spectrum", 0,0, 1000, 700);
  c->Divide(3,1);
  TH1F *histo_exp=extraction(); // Get the experimental histo took by the acquision
  TH1F *hspectrum_F = new TH1F( "hspectrum_F", "Fermi spectrum", 100, 0., 1.2*Q); // Creation the histogram of the Fermi distribution
  hspectrum_F->SetXTitle("Energy MeV");
  hspectrum_F->SetStats(0);
  double cut_off = 1100; //keV

  //Initialisation of the Fermi function distribution 
  for ( int i=0; i<Nevts; i++) {
    double_t x=randomFermi(Q,(int)Z, cut_off);
    x=EtoP(x);
    hspectrum_F->Fill(x);}
  hspectrum_F->Scale(1./Nevts); 

  double xmin = hspectrum_F->GetXaxis()->GetXmin();
  double xmax = hspectrum_F->GetXaxis()->GetXmax();
  hspectrum_F->GetXaxis()->SetLimits(xmin/1000,xmax/1000);

  hspectrum_F->Scale(1.0 / hspectrum_F->Integral()); //Normalization (area=1)
  histo_exp->Scale(1.0 / histo_exp->Integral()); //Normalization (area=1)

  histo_exp->Scale(hspectrum_F->Integral() / histo_exp->Integral());
  hspectrum_F->GetXaxis()->SetLimits(0, 10); //Same scale as histo_exp

  //TH1F *histo_deconv=ApproximateDeconv(histo_exp,hspectrum_F); //Approximative deconvolution 
  TH1F *histo_deconv=RichardsonLucyDeconv(histo_exp,hspectrum_F); //Deconvolution of RichardsonLucy
  histo_deconv->SetTitle("Uncertainty distribution function");

  // Plot histogram
  c->cd(1);
  hspectrum_F->Draw("HIST");
  c->cd(2);
  histo_exp->Draw("HIST");
  c->cd(3);
  histo_deconv->Draw("HIST");

  // Distribution function of the deconvolution:
  TF1* func_deconv=HistToFunc(histo_deconv,"func_histo");
  TCanvas* c1=new TCanvas("Distribution function", "Distribution function", 0,0, 1000, 700);
  func_deconv->Draw("");
}
