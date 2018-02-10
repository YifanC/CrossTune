#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#include <sstream>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TView.h>
#include <TView3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TVector3.h>
#include <TPrincipal.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>

using namespace std;

const Char_t *inputFileName = "laserResults.root";

const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;

struct calibTrackInfo
{
  Double_t x0_calib;
  Double_t y0_calib;
  Double_t z0_calib;
  Double_t theta_calib;
  Double_t phi_calib;
};

vector<Double_t> findClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB);
Double_t doCoordTransformX(const Double_t inputX);
Double_t doCoordTransformY(const Double_t inputY);
Double_t doCoordTransformZ(const Double_t inputZ);

Int_t main(Int_t argc, Char_t** argv)
{
  if(argc != 12) {
    cout << endl << endl << "Please enter ten parameters: dX_upstream dY_upstream dZ_upstream dTheta_upstream dPhi_upstream dX_downstream dY_downstream dZ_downstream dTheta_downstream dPhi_downstream" << endl << endl;
    return -1;
  }

  const Double_t dX_upstream = doCoordTransformX(atof(argv[1])) - doCoordTransformX(0.0);
  const Double_t dY_upstream = doCoordTransformY(atof(argv[2])) - doCoordTransformY(0.0);
  const Double_t dZ_upstream = doCoordTransformZ(atof(argv[3])) - doCoordTransformZ(0.0);
  const Double_t dTheta_upstream = atof(argv[4]);  // NOTE:  No coordinate change has been done!  Modify?
  const Double_t dPhi_upstream = atof(argv[5]);  // NOTE:  No coordinate change has been done!  Modify?
  
  const Double_t dX_downstream = doCoordTransformX(atof(argv[6])) - doCoordTransformX(0.0);
  const Double_t dY_downstream = doCoordTransformY(atof(argv[7])) - doCoordTransformY(0.0);
  const Double_t dZ_downstream = doCoordTransformZ(atof(argv[8])) - doCoordTransformZ(0.0);
  const Double_t dTheta_downstream = atof(argv[9]);  // NOTE:  No coordinate change has been done!  Modify?
  const Double_t dPhi_downstream = atof(argv[10]);  // NOTE:  No coordinate change has been done!  Modify?

  const int step = atoi(argv[11]);

  //Double_t theta = acos((y1-y0)/trackLength);  // Definition of Theta in Mike's Coordinates
  //Double_t phi = acos((z1-z0)/(trackLength*sin(theta)));  // Definition of Phi in Mike's Coordinates

  TH1F *crossHist = new TH1F("crossHist","",120,-15.0,15.0);
  TH2F *crossHist2D = new TH2F("crossHist2D","",40,0.0,40.0,40,0.0,40.0);

  TFile *inputFile = new TFile(inputFileName,"READ");
  
  TTreeReader readerCrossings("SpaCEtree_crossings", inputFile);
  TTreeReaderValue<Int_t> trackNum1(readerCrossings, "trackNum1.data_crossings");
  TTreeReaderValue<Int_t> trackNum2(readerCrossings, "trackNum2.data_crossings");
  TTreeReaderValue<Double_t> crossX(readerCrossings, "crossX.data_crossings");
  TTreeReaderValue<Double_t> crossY(readerCrossings, "crossY.data_crossings");
  TTreeReaderValue<Double_t> crossZ(readerCrossings, "crossZ.data_crossings");
  TTreeReaderValue<Double_t> crossDist(readerCrossings, "crossDist.data_crossings");
  TTreeReaderValue<Double_t> crossDistX(readerCrossings, "crossDistX.data_crossings");
  TTreeReaderValue<Double_t> crossDistY(readerCrossings, "crossDistY.data_crossings");
  TTreeReaderValue<Double_t> crossDistZ(readerCrossings, "crossDistZ.data_crossings");
  TTreeReaderValue<Double_t> crossX_mod(readerCrossings, "crossX_mod.data_crossings");
  TTreeReaderValue<Double_t> crossY_mod(readerCrossings, "crossY_mod.data_crossings");
  TTreeReaderValue<Double_t> crossZ_mod(readerCrossings, "crossZ_mod.data_crossings");
  TTreeReaderValue<Double_t> crossDist_mod(readerCrossings, "crossDist_mod.data_crossings");
  TTreeReaderValue<Double_t> crossDistX_mod(readerCrossings, "crossDistX_mod.data_crossings");
  TTreeReaderValue<Double_t> crossDistY_mod(readerCrossings, "crossDistY_mod.data_crossings");
  TTreeReaderValue<Double_t> crossDistZ_mod(readerCrossings, "crossDistZ_mod.data_crossings");
  TTreeReaderValue<Double_t> distWeight(readerCrossings, "distWeight.data_crossings");
  TTreeReaderValue<Int_t> crossType(readerCrossings, "crossType.data_crossings");

  TTreeReader readerTracks("SpaCEtree_tracks", inputFile);
  TTreeReaderValue<Double_t> trackX0(readerTracks, "x0_tracks");
  TTreeReaderValue<Double_t> trackY0(readerTracks, "y0_tracks");
  TTreeReaderValue<Double_t> trackZ0(readerTracks, "z0_tracks");
  TTreeReaderValue<Double_t> trackTheta(readerTracks, "theta_tracks");
  TTreeReaderValue<Double_t> trackPhi(readerTracks, "phi_tracks");

  vector<calibTrackInfo> laserTrackSet;
  
  while (readerTracks.Next())
  {
    calibTrackInfo laserTrack;

    if(*trackZ0 < Lz/2.0) {
      laserTrack.x0_calib = *trackX0 + dX_upstream;
      laserTrack.y0_calib = *trackY0 + dY_upstream;
      laserTrack.z0_calib = *trackZ0 + dZ_upstream;
      laserTrack.theta_calib = *trackTheta + dTheta_upstream;
      laserTrack.phi_calib = *trackPhi + dPhi_upstream;
    }
    else {
      laserTrack.x0_calib = *trackX0 + dX_downstream;
      laserTrack.y0_calib = *trackY0 + dY_downstream;
      laserTrack.z0_calib = *trackZ0 + dZ_downstream;
      laserTrack.theta_calib = *trackTheta + dTheta_downstream;
      laserTrack.phi_calib = *trackPhi + dPhi_downstream;
    }
    
    laserTrackSet.push_back(laserTrack);
  }

  while (readerCrossings.Next())
  {
    vector<Double_t> POAparams = findClosestPOA(laserTrackSet[*trackNum1],laserTrackSet[*trackNum2]);
    Double_t distVal = POAparams.at(0);
    Double_t distValDistorted = *crossDist_mod;

    if((distVal < 0.0) || (distVal > 0.1)) continue;
    if(distValDistorted > 0.1) continue;

    crossHist->Fill(100.0*(distValDistorted-distVal));
    crossHist2D->Fill(100.0*distVal,100.0*distValDistorted);
  }
  crossHist->Scale(1.0/crossHist->Integral());

  cout << atof(argv[1]) << " " << atof(argv[2]) << " " << atof(argv[3]) << " " << atof(argv[4]) << " " << atof(argv[5]) << " " << atof(argv[6]) << " " << atof(argv[7]) << " " << atof(argv[8]) << " " << atof(argv[9]) << " " << atof(argv[10]) << " " << crossHist->GetMean() << " " << crossHist->GetRMS() << endl;

  inputFile->Close();

  TCanvas *c_cross = new TCanvas();
  c_cross->cd();
  crossHist->GetXaxis()->SetTitle("Reco Crossing Distance - True Crossing Distance [cm]");
  crossHist->GetXaxis()->SetTitleOffset(0.95);
  crossHist->GetXaxis()->SetTitleSize(0.045);
  crossHist->GetYaxis()->SetTitle("Arb. Units");
  crossHist->GetYaxis()->SetTitleOffset(0.95);
  crossHist->GetYaxis()->SetTitleSize(0.05);
  crossHist->GetYaxis()->SetNoExponent(kTRUE);
  crossHist->SetStats(0);
  crossHist->SetLineWidth(2.0);
  crossHist->Draw("HIST");
  crossHist->Draw("AXISsame");
  crossHist->SetMinimum(0.0001);

  //cout << crossHist->GetMean() << " "  << crossHist->GetRMS() << endl;
  ostringstream stream;
  stream << "corssHist-" << step << ".png";
  std:string filename = stream.str();
  c_cross->SaveAs(filename.c_str());
  
  return 0;
}

vector<Double_t> findClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB)
{ 
  //   (xA,yA,zA)+t(xA_step,yA_step,zA_step)
  Double_t xA = calibTrackA.x0_calib, yA = calibTrackA.y0_calib, zA = calibTrackA.z0_calib;
  Double_t xB = calibTrackB.x0_calib, yB = calibTrackB.y0_calib, zB = calibTrackB.z0_calib;

  Double_t xA_step = -1.*sin(calibTrackA.theta_calib)*sin(calibTrackA.phi_calib);
  Double_t yA_step = cos(calibTrackA.theta_calib);
  Double_t zA_step = sin(calibTrackA.theta_calib)*cos(calibTrackA.phi_calib);

  Double_t xB_step = -1.*sin(calibTrackB.theta_calib)*sin(calibTrackB.phi_calib);
  Double_t yB_step = cos(calibTrackB.theta_calib);
  Double_t zB_step = sin(calibTrackB.theta_calib)*cos(calibTrackB.phi_calib);

  //perpendicular line between the two tracks
  Double_t x_prep = (yA_step*zB_step)-(yB_step*zA_step);
  Double_t y_prep = (xB_step*zA_step)-(xA_step*zB_step);
  Double_t z_prep = (xA_step*yB_step)-(xB_step*yA_step);
 
  // if cross product is zero then the lines are parallel so return distance = -2
  if (x_prep == 0 && y_prep == 0 && z_prep == 0) {
    vector<Double_t> return_vector;
    return_vector.push_back(-2.);
    return_vector.push_back((xA+xB)/2.);
    return_vector.push_back((yA+yB)/2.);
    return_vector.push_back((zA+zB)/2.);   
    return return_vector;
  }
  //normalize the perpendicular line
  Double_t mag_prep = sqrt(pow(x_prep,2)+pow(y_prep,2)+pow(z_prep,2));
  Double_t x_prep_norm = x_prep / mag_prep;
  Double_t y_prep_norm = y_prep / mag_prep;
  Double_t z_prep_norm = z_prep / mag_prep;
 
  //defined to make the math simplier
  Double_t a = y_prep_norm*(xA-xB);
  Double_t b = x_prep_norm*(yA-yB);
  Double_t c = z_prep_norm*(xA-xB);
  Double_t d = x_prep_norm*(zA-zB);
 
  Double_t g = y_prep_norm*xA_step;
  Double_t h = y_prep_norm*xB_step;
  Double_t i = x_prep_norm*yA_step;
  Double_t j = x_prep_norm*yB_step;
  Double_t k = z_prep_norm*xA_step;
  Double_t l = z_prep_norm*xB_step;
  Double_t m = x_prep_norm*zA_step;
  Double_t n = x_prep_norm*zB_step;
 
  Double_t chi = (l-n)/(h-j);
 
  //alpha: "t" for the first line //beta:  "t" for the second line
  Double_t alpha = (chi*(a-b)-c+d)/(k-m-(chi*(g-i)));
  Double_t beta = (c-d+alpha*(k-m))/(l-n);

  Double_t cpa_xA = xA+alpha*xA_step;
  Double_t cpa_yA = yA+alpha*yA_step;
  Double_t cpa_zA = zA+alpha*zA_step;

  Double_t cpa_xB = xB+beta*xB_step;
  Double_t cpa_yB = yB+beta*yB_step;
  Double_t cpa_zB = zB+beta*zB_step;
 
  //distance between the closest points on the lines\tracks
  //Double_t distance = sqrt(pow((xA+alpha*xA_step)-(xB+beta*xB_step),2)+pow((yA+alpha*yA_step)-(yB+beta*yB_step),2)+pow((zA+alpha*zA_step)-(zB+beta*zB_step),2));
  Double_t distance = sqrt(pow(cpa_xA-cpa_xB,2)+pow(cpa_yA-cpa_yB,2)+pow(cpa_zA-cpa_zB,2));

  //midpoint between points of closest approach on both lines
  Double_t x_mid = (cpa_xA+cpa_xB)/2.;
  Double_t y_mid = (cpa_yA+cpa_yB)/2.;
  Double_t z_mid = (cpa_zA+cpa_zB)/2.;
 
  //check to see if this point is outside the detector
  //  TODO: < or <=, Is 0 or Lx "outside"?
  //if (x_mid > Lx || x_mid < 0 || y_mid > Ly || y_mid < 0 || z_mid > Lz || z_mid < 0) distance = -1;
  if (cpa_xA > Lx || cpa_xA < 0 || cpa_yA > Ly || cpa_yA < 0 || cpa_zA > Lz || cpa_zA < 0 ||
      cpa_xB > Lx || cpa_xB < 0 || cpa_yB > Ly || cpa_yB < 0 || cpa_zB > Lz || cpa_zB < 0) distance = -1;
 
  vector<Double_t> return_vector;

  return_vector.push_back(distance);
  return_vector.push_back(x_mid);
  return_vector.push_back(y_mid);
  return_vector.push_back(z_mid);

  ////// NEW 12/5/2017 //////
  return_vector.push_back(cpa_xB-cpa_xA);
  return_vector.push_back(cpa_yB-cpa_yA);
  return_vector.push_back(cpa_zB-cpa_zA);
  ///////////////////////////
  
  return return_vector;
}

Double_t doCoordTransformX(const Double_t inputX)
{
  Double_t outputX;

  outputX = Lx - (Lx/2.58)*inputX/100.0;

  return outputX;
}

Double_t doCoordTransformY(const Double_t inputY)
{
  Double_t outputY;

  outputY = (Ly/(1.170+1.151))*(inputY+115.1)/100.0;

  return outputY;
}

Double_t doCoordTransformZ(const Double_t inputZ)
{
  Double_t outputZ;

  outputZ = (Lz/(10.365+0.007))*(inputZ+0.7)/100.0;
  
  return outputZ;
}
