//UF, calculates probability for getting an E115-like chain of random origin
#include <iostream>
#include <string>
#include <sstream>
#include "TH1.h"
#include "TFile.h"
#include <fstream>

using namespace std;

//Functions 
double probability(int to_obs, float expected);
int fact(Int_t k);
double rate_calc(int pixel_number, int lower_limit, int upper_limit, TH1F **h_energy, double time);


//Main program starts here
int u_main() {


  //Prepare/read in root spectra
  TFile  *file = new TFile("Spectra.root");

  TH1F *h_energy_pixel_off_tot = (TH1F*)file->Get("h_energy_pixel_off_tot");
  TH1F *h_energy_pixel_recoff_tot = (TH1F*)file->Get("h_energy_pixel_recoff_tot"); 
  TH1F *h_energy_pixel_on_tot = (TH1F*)file->Get("h_energy_pixel_on_tot"); 
  TH1F *h_energy_pixel_recon_tot = (TH1F*)file->Get("h_energy_pixel_recon_tot"); 

  h_energy_pixel_recoff_tot->Draw();


  TH1F* h_energy_pixel_on[1024];
  TH1F* h_energy_pixel_off[1024];
  TH1F* h_energy_pixel_recon[1024];
  TH1F *h_energy_pixel_recoff[1024];

  //Rates per pixel
  long double rate_a1[1024]; 
  long double rate_a2[1024]; 
  long double rate_f[1024]; 

  //Relevant times
  //  double exp_time = 1550000; //Approximate number
  //double exp_time_fission = 1475000; //"Total beam off" (strålstruktur/avstängd)
  //double exp_time_alpha_on = 1433000; //"BOT" innehåller ~25% faktisk stråltid 
  double exp_time = 1433000;//1475000; //"Total beam off" 

  //  double beamoff_time = exp_time*3/4;
  //  double beamon_time = exp_time*1/4;

  Int_t temporary;
  int EVR[1024] = {0};

  int temp_value = 0;
  double fission_array[1024] = {0};  
 
  int number_of_alphas;

 int low_limit_first = 0; 
 int high_limit_first = 0; 
 int low_limit_second= 0; 
 int high_limit_second = 0; 

 int pixel = 0;
 pixel = 681;

  //Calc FISSION numbers
  //string fission_model = "one_or_zero";
  string fission_model = "one_or_average";  
  //string fission_model = "one";

  if(fission_model == "one_or_zero"){
    ifstream input("Fission-BeamOff-Pixel_modified.txt",ios::in);
    while (input){
      input>>temp_value;
      if(fission_array[temp_value] == 1) cout<<temp_value<<endl;
      fission_array[temp_value] = fission_array[temp_value] + 1;
    }
    //FULFIX!!!!!! P.G.A. att sista värdet tas med två gånger 
    fission_array[621] = 1;
  }

  if(fission_model == "one_or_average"){
    ifstream input("Fission-BeamOff-Pixel_modified.txt",ios::in);
    while (input){
      input>>temp_value;
      if(fission_array[temp_value] == 1) cout<<temp_value<<endl;
      fission_array[temp_value] = fission_array[temp_value] + 1;
    }
    //FULFIX!!!!!! P.G.A. att sista värdet tas med två gånger 
    fission_array[621] = 1;    
    for(int i = 0; i < 1024; i++){
      if(fission_array[i] == 0){
	fission_array[i] = 63.0/1024.0;    
      } 
    }
    
  }//if "one_or_average"
  
  if(fission_model == "one"){
    for(int i = 0; i < 1024; i++){
      fission_array[i] = 1;
    }
  }


  int reaction_channel = 3;
  float lifetimes = 3;
  char lookslike[7];
  float timespan[7] = {0};

  float lifetime3n1 = 0.171;  //Lifetimes in seconds, from Oganessian et al., Phys.Rev. C 87, 014302 (2013) 
  float lifetime3n2 = 0.97;
  float lifetime3n3 = 3.6;
  float lifetime3n4 = 0.54;
  float lifetime3n5 = 12.0;
  float lifetime3n6 = 27*60*60;

  float lifetime2n1 = 0.380;
  float lifetime2n2 = 5.6;
  float lifetime2n3 = 22.0;

  double random_prob_pixel[1024] = {0};
  double random_prob_sum = 0; 

  char chead[64];

  string a1;
  string a2;


  //Read in spectra, cont.
  for(Int_t i = 0; i < 1024; i++){ 
    sprintf(chead,"on/h_energy_pixel_on_%d",i);
    h_energy_pixel_on[i] = (TH1F*)file->Get(chead);   
  } 
  for(Int_t i = 0; i < 1024; i++){ 
    sprintf(chead,"off/h_energy_pixel_off_%d",i);
    h_energy_pixel_off[i] = (TH1F*)file->Get(chead);   
  } 
  for(Int_t i = 0; i < 1024; i++){ 
    sprintf(chead,"recon/h_energy_pixel_recon_%d",i);
    h_energy_pixel_recon[i] = (TH1F*)file->Get(chead);   
  } 
  for(Int_t i = 0; i < 1024; i++){ 
    sprintf(chead,"recoff/h_energy_pixel_recoff_%d",i);
    h_energy_pixel_recoff[i] = (TH1F*)file->Get(chead);   
  } 

  //h_energy_pixel_off[0]->Draw();
  //h_energy_pixel_recoff[0]->Draw();


  //Communication with user
  //cout << "What pixel?" << endl;
  //cin >> pixel_number;
  //cout << "What reaction channel?" << endl;
  //cin >> reaction_channel;
  //cout << "How many lifetimes?" << endl;
  //cin >> lifetimes;

  cout << "Number of alphas? " << endl;
  cin >> number_of_alphas;
  //  cout << "First decay is Alphalike/Escapelike/Fissionlike" << endl;
  //cin >> lookslike[0];
  cout << "How long time span for first alpha?" << endl;
  cin >> timespan[0];
  //timespan[0] = 2;
  cout << "What does it look like? <a/e> " << endl;
  cin >> lookslike[0];
  cout << "BeamOn or beamOff? <on/off>" << endl; 
  cin >> a1;
  //cout << "Second alpha/fission is Alphalike/Escapelike/Fissionlike" << endl;
  //cin >> lookslike[1];
  if(number_of_alphas == 2){
    cout << "How long time span for second alpha?" << endl;
    cin >> timespan[1];
    cout << "What does it look like? <a/e> " << endl;
    cin >> lookslike[1];
    cout << "BeamOn of beamOff? <on/off>" << endl; 
    cin >> a2;
  }
  //cout << "Thirds alpha/fission is Alphalike/Escapelike/Fissionlike" << endl;
  //cin >> lookslike[2];
  cout << "How long time span for fission?" << endl;
  cin >> timespan[2];


  //Calc number of EVR
  for (Int_t j = 0; j < 1024; j++)
    {
      temporary = 0;
      for (Int_t i = 1100; i < 1800; i++)
	{
	  temporary = temporary + h_energy_pixel_on[j]->GetBinContent(i);
	}
      EVR[j] = temporary;
    }

  //Calc rates for FIRST ALPHA based on lower and upper limits

  if(lookslike[0] == 'a'){
    low_limit_first = 900;
    high_limit_first = 1100; 
  }
  else if(lookslike[0] == 'e'){
    low_limit_first = 0;
    high_limit_first = 400;
  }
  else{
    cout << "Något är fel med inmatningen! " << lookslike[0] << endl;
  }

  if(a1 == "on"){
    for(int i = 0; i < 1024; i++){
      rate_a1[i] = rate_calc(i, low_limit_first, high_limit_first, h_energy_pixel_recon, exp_time);
    }
  }
  if(a1 == "off"){
    for(int i = 0; i < 1024; i++){
      rate_a1[i] = rate_calc(i, low_limit_first, high_limit_first, h_energy_pixel_recoff, exp_time);
    }
  }


  //Calc rates for SECOND ALPHA based on lower and upper limits

  if(number_of_alphas == 2){
    
    if(lookslike[1] == 'a'){
      low_limit_second = 900;
      high_limit_second = 1100; 
    }
    else if(lookslike[1] == 'e'){
      low_limit_second = 0;
      high_limit_second = 400;
    }
    else{
      cout << "Något är fel med inmatningen! " << lookslike[1] << endl;
    }
    
    if(a2 == "on"){
      for(int i = 0; i < 1024; i++){
	rate_a2[i] = rate_calc(i, low_limit_second, high_limit_second, h_energy_pixel_recon, exp_time);
      }
    }
    if(a2 == "off"){
      //timespan[1] = timespan[1] * 3.0/4.0;  //We only look for the alpha during the beam off time
      for(int i = 0; i < 1024; i++){
	rate_a2[i] = rate_calc(i, low_limit_second, high_limit_second, h_energy_pixel_recoff, exp_time);
      }
    }
  }//ends if number_of_alphas == 2




  //Calc rates for FISSION
  for(int i = 0; i < 1024; i++){
    rate_f[i] = fission_array[i]/exp_time;
  }
  
  
  //Calc number of expected random chains 
  random_prob_sum = 0;
  for(Int_t i = 0; i < 1024;i++){

    random_prob_pixel[i] = EVR[i]*
      (1-probability(0,rate_a1[i]*timespan[0]))*  //FIRST ALPHA
      (1-probability(0,rate_f[i]*timespan[2]));   //FISSION

    if(i == pixel){ 
      cout << "Rate a1 in pixel " << rate_a1[pixel] << endl;
      cout << "Rate f in pixel " << rate_f[pixel] << endl;
    }

    if(number_of_alphas == 2){
      random_prob_pixel[i] = random_prob_pixel[i] * (1-probability(0,rate_a2[i]*timespan[1])); 
      if(i == pixel){ 
	cout << "Rate a2 in pixel " << rate_a2[pixel] << endl;
      }
    }

    random_prob_sum = random_prob_sum + random_prob_pixel[i];    
  }

 

 cout << EVR[pixel] << "  " << probability(1,rate_a1[pixel]*timespan[0]) << " " << probability(1,rate_f[pixel]*timespan[2]) << endl;

 cout << "Random prob pixel is " << random_prob_pixel[pixel] << endl;   

  cout <<  "The number of expected chains of this type generated by stochastic fluctuations in the background is  " << random_prob_sum << endl;


  return 0; //ends "main()"

}



// ***********************************
//		Methods
// ***********************************

//Calculate the probability of observing "to_obs" events if "expected" are expected (rate * time window)
double probability(int to_obs, float expected){
  double prob;
  prob = (exp(-expected)*pow(expected,to_obs))/(fact(to_obs));
  return prob;
}


//Calculate factorial of value
int fact(Int_t k){
  int factorial;
  factorial = 1;
  for (Int_t j = 1; j <= k; j++){
    factorial = factorial*j;
  }
  //cout << "factorial is " << factorial << endl; //Function called multiple times, each one gives a printout, of course!
  return (factorial);
}


//Calculate rate of a certain type of events, between lower_limit and upper_limit
double rate_calc(int pixel_number, int lower_limit, int upper_limit, TH1F **h_energy, double time){
  int temporary = 0;
  long double rate = 0;
      temporary = 0;
      for (Int_t i = lower_limit; i < upper_limit; i++)
	{                   
	  temporary = temporary + h_energy[pixel_number]->GetBinContent(i);
	  // cout << h_energy[pixel_number]->GetBinContent(i) << endl;
	}
      rate = temporary/time;      
      if(pixel_number == 681) cout << "rate is " << rate << " and temporary is " << temporary <<endl;   
 return rate;
} //end of function "rate"










