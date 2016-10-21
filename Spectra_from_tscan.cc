//UF, calculates probability for a chain to be random 
#include <iostream>
#include <string>
#include <sstream>
#include "TH1.h"
#include "TFile.h"
#include "TF1.h"
#include <fstream>
#include <ctime>
#include <random>
#include <iomanip>
#include <algorithm>
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"


using namespace std;

void run_main();

int main() {

//Make histograms look nice
TGaxis::SetMaxDigits(3);
gStyle->SetCanvasColor(kWhite);    
gStyle->SetFrameLineWidth(2);       
gStyle->SetFrameLineColor(kWhite); //"No frame"
gStyle->SetPadLeftMargin(0.13);
gStyle->SetPadTopMargin(0.11);
gStyle->SetPadRightMargin(0.1);
gStyle->SetPadBottomMargin(0.15);
//cout << "Margins are " << gStyle->GetPadLeftMargin() << " " << gStyle->GetPadTopMargin() << " " << gStyle->GetPadRightMargin() << " " << gStyle->GetPadBottomMargin() << " " << endl;
gStyle->SetTitleSize(0.05,"xy");
//gStyle->SetTitleOffset(12,"x"); //1.2   Doesn't do anything!
//gStyle->SetTitleOffset(10,"y"); //1.0   Doesn't do anything!
gStyle->SetOptStat();
gStyle->SetLegendBorderSize(0);
gStyle->SetLineWidth(2.); //X and Y axis thickness 
gStyle->SetLabelSize(0.05,"xy"); //Size of numbers on X and Y axes  
gStyle->SetLabelOffset(0.01,"xy"); //Offset of numbers on X and Y axes
gStyle->SetLabelFont(42,"xy"); //Font on labels on X and Y axes
//gStyle->SetHistoLineWidth(5);
gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
gStyle->SetPadBorderMode(0);
//gStyle->SetPaintTextFormat("5.2f");  // What precision to put numbers if plotted with "TEXT"
gStyle->SetPaintTextFormat("g");  // What precision to put numbers if plotted with "TEXT"
gStyle->SetNumberContours(20);
//gStyle->SetTextSize(2.1);
gStyle->SetTextFont(42); 
gStyle->SetStripDecimals(0); 
gStyle->SetOptStat(1111111);




//Input files to use 
ifstream recoff("recoff.txt",ios::in);
ifstream recon("recon.txt",ios::in);
ifstream off("off.txt",ios::in);
ifstream on("on.txt",ios::in);

//Initialise histograms
TH1F* h_energy_pixel_off[1024];
TH1F* h_energy_pixel_on[1024];
TH1F* h_energy_pixel_recoff[1024];
TH1F* h_energy_pixel_recon[1024];
TH1F* h_energy_pixel_off_tot;
TH1F* h_energy_pixel_on_tot;
TH1F* h_energy_pixel_recoff_tot;
TH1F* h_energy_pixel_recon_tot;

char chis[64];
char chead[64];

//Declare histograms
for(Int_t i = 0; i < 1024; i++){ 
sprintf(chis,"h_energy_pixel_off_%d",i);
sprintf(chead,"energy_pixel_off_%d",i);
h_energy_pixel_off[i] = new TH1F(chis,chead,4096,0,40.96); 
} 
for(Int_t i = 0; i < 1024; i++){ 
sprintf(chis,"h_energy_pixel_on_%d",i);
sprintf(chead,"energy_pixel_on_%d",i);
h_energy_pixel_on[i] = new TH1F(chis,chead,4096,0,40.96); 
} 
for(Int_t i = 0; i < 1024; i++){ 
sprintf(chis,"h_energy_pixel_recoff_%d",i);
sprintf(chead,"Energy, pixel %d (Experimental data)",i);
h_energy_pixel_recoff[i] = new TH1F(chis,chead,4096,0,40.96); 
} 
for(Int_t i = 0; i < 1024; i++){ 
sprintf(chis,"h_energy_pixel_recon_%d",i);
sprintf(chead,"Energy, Beam OFF, pixel %d (Experimental data)",i);
h_energy_pixel_recon[i] = new TH1F(chis,chead,4096,0,40.96); 
} 

h_energy_pixel_recon_tot = new TH1F("h_energy_pixel_recon_tot","Energy, Total (Experimental data)",4096,0,40.960); 
h_energy_pixel_recoff_tot = new TH1F("h_energy_pixel_recoff_tot","Energy, Total, Beam OFF (Experimental data)",4096,0,40.960); 
h_energy_pixel_on_tot = new TH1F("h_energy_pixel_on_tot","on_tot",4096,0,40.960); 
h_energy_pixel_off_tot = new TH1F("h_energy_pixel_off_tot","off_tot",4096,0,40.960); 

//Create root file with subdirectories
//TFile *myoutput =  new TFile("Spectra.root", "RECREATE");
TFile *myoutput =  new TFile("Spectra_test.root", "RECREATE");
myoutput->mkdir("on");
myoutput->mkdir("off");
myoutput->mkdir("recon");
myoutput->mkdir("recoff");

//Fill spectra from the files from tscan into ordinary root histos
Int_t temp_sum = 0;
Int_t temp = 0;
Int_t n = 0;
Int_t l = 0;
Int_t m = 0; 
while (on>>temp){
if(m%4096==0 && m > 0) 
{
//cout << "Cycle " << n <<" finished" << endl;
n++; 
m=0; 
}
h_energy_pixel_on[n]->SetBinContent(m,temp);
h_energy_pixel_on_tot->AddBinContent(m,temp);
m++; 
l++;
}
//cout << "l is " << l << endl;
//cout << "temp is " << temp << endl;

l=0;
n=0;
m=0; 
while (off>>temp){
if(m%4096==0 && m > 0) 
{
//cout << "Cycle " << n <<" finished" << endl;
n++; 
m=0; 
}
h_energy_pixel_off[n]->SetBinContent(m,temp);
h_energy_pixel_off_tot->AddBinContent(m,temp);
m++; 
l++;
}
//cout << "l is " << l << endl;
//cout << "temp is " << temp << endl;  

l=0;
n=0;
m=0; 
while (recon>>temp){
if(m%4096==0 && m > 0) 
{
//cout << "Cycle " << n <<" finished" << endl;
n++; 
m=0; 
}
h_energy_pixel_recon[n]->SetBinContent(m,temp);
h_energy_pixel_recon_tot->AddBinContent(m,temp);
m++; 
l++;
}

l=0;
n=0;
m=0; 
while (recoff>>temp){
if(m%4096==0 && m > 0) 
{
//cout << "Cycle " << n <<" finished" << endl;
n++; 
m=0; 
}
h_energy_pixel_recoff[n]->SetBinContent(m,temp);
h_energy_pixel_recoff_tot->AddBinContent(m,temp);
m++; 
l++;
}
//Spectra filled!

//Write spectra to root file
for(Int_t k = 0;k<1024;k++){
myoutput->cd("off");
h_energy_pixel_off[k]->Write(); 
myoutput->cd(); 
}
for(Int_t k = 0;k<1024;k++){
myoutput->cd("on");
h_energy_pixel_on[k]->Write();
myoutput->cd(); 
}
for(Int_t k = 0;k<1024;k++){
myoutput->cd("recoff");
h_energy_pixel_recoff[k]->Write(); 
myoutput->cd(); 
}
for(Int_t k = 0;k<1024;k++){
myoutput->cd("recon");
h_energy_pixel_recon[k]->Write();   
myoutput->cd();
}

h_energy_pixel_off_tot->Write();
h_energy_pixel_on_tot->Write();
h_energy_pixel_recoff_tot->Write();
h_energy_pixel_recon_tot->Write();

myoutput->Close();


//Draw


//h_energy_pixel_off_tot->Scale(10);
TAxis *Xaxis = h_energy_pixel_off_tot->GetXaxis();
TAxis *Yaxis = h_energy_pixel_off_tot->GetYaxis();

Yaxis->SetTitle("Counts per 10 keV");
Xaxis->SetTitle("Energy (keV)");  
h_energy_pixel_off_tot->SetTitle("");
h_energy_pixel_off_tot->SetStats(0);
TGaxis::SetMaxDigits(5); 
Xaxis->SetTitleOffset(1.2);
Xaxis->SetTitleSize(0.05);
Xaxis->SetRangeUser(5010,11990);      //(5,11995);       //(9510,11490); //4505,11995
  //Xaxis->SetNdivisions(505);  

  Yaxis->SetTitleOffset(1.3);
  //Yaxis->SetNdivisions(606);
  //Yaxis->SetDecimals(1); //2 
  //Yaxis->SetRangeUser(0.001,24); 

  Xaxis->CenterTitle();
  Yaxis->CenterTitle();



  TCanvas *MyCanvas = new TCanvas("MyCanvas");
  MyCanvas->Divide(1,1);  
 

  
 MyCanvas->cd(1); 

 //gPad->SetLogy(); 



  h_energy_pixel_off_tot->SetLineWidth(1);
  h_energy_pixel_off_tot->Draw();

  h_energy_pixel_recoff_tot->SetLineColor(8);
  h_energy_pixel_recoff_tot->SetLineWidth(3);
  h_energy_pixel_recoff_tot->SetLineStyle(1);
  h_energy_pixel_recoff_tot->Draw("same");
  h_energy_pixel_off_tot->Draw("same");

 
  TCanvas *theCanvas = new TCanvas("theCanvas");
  theCanvas->cd();
  h_energy_pixel_recon[500]->Draw();
  h_energy_pixel_recon[100]->Draw("same");
  theCanvas->Update();

  /*
  MyCanvas->cd(2);
 
  gPad->SetLogy(); 


  Xaxis2 = h_energy_pixel_on_tot->GetXaxis();
  Yaxis2 = h_energy_pixel_on_tot->GetYaxis();
  Yaxis2->SetTitle("Counts per 10 keV channel");
  Xaxis2->SetTitle("Energy (keV)");  
  h_energy_pixel_on_tot->SetTitle("");
  h_energy_pixel_on_tot->SetStats(0);

  Xaxis2->SetTitleOffset(1.0);
  Xaxis2->SetTitleSize(0.05);
  //Xaxis2->SetRangeUser(5,17995); //
  Xaxis2->SetNdivisions(505);   

  //Yaxis->SetTitleOffset(1.0);

  //Yaxis->SetNdivisions(606);
  //Yaxis->SetDecimals(1); //2 
  //Yaxis->SetRangeUser(0.0001,3.399); 
  Xaxis2->SetNdivisions(505);   

  Xaxis2->CenterTitle();
  Yaxis2->CenterTitle();


  h_energy_pixel_on_tot->SetLineWidth(1);
  h_energy_pixel_on_tot->Draw();
  h_energy_pixel_recon_tot->SetLineColor(8);
  h_energy_pixel_recon_tot->SetLineWidth(3);
  h_energy_pixel_recon_tot->SetLineStyle(1);
  h_energy_pixel_recon_tot->Draw("same");
  h_energy_pixel_on_tot->Draw("same"); //Must be more visible
  */
  


 return 0; //ends "main()"
}


// ***********************************
//		Methods
// ***********************************

void run_main() {
	main();
}
