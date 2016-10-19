#include <iostream>
#include <fstream>
#include "RandomChains.h"
#include <assert.h>
#include "TCanvas.h"
#include "TRint.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"

using namespace std;


RandomChains::RandomChains() {

	//Default value for number of pixels in the implantation detector
	nbr_pixels = 1024;
	
	//Interval for accepted superheavy nuclei alpha decays (10*keV)
	lower_limit_alphas = 900; upper_limit_alphas = 1100;

	//Interval for energy deposit of alphas that escapes the implantation detector (10*keV)
	lower_limit_escapes = 0; upper_limit_escapes = 400;

	//Interval for implanted nuclei with beam = ON (10*keV)
	lower_limit_implants = 1100; upper_limit_implants = 1800;


	
	/*
	TH1F* h_energy_pixel_reconstructed_beam_on[nbr_pixels];
	TH1F* h_energy_pixel_reconstructed_beam_off[nbr_pixels];
	TH1F* h_energy_pixel_reconstructed__beam_off = (TH1F*) calloc(sizeof(TH1F), nbr_pixels);
	assert(h_energy_pixel_reconstructed__beam_off);
	TH1F* h_energy_pixel_reconstructed__beam_on = (TH1F*) calloc(sizeof(TH1F), nbr_pixels);
	assert(h_energy_pixel_reconstructed__beam_on);
	
	*/


	Bool_t valid_input = kFALSE;

	while(!valid_input) {
		cout << "What type of run? <0/1/2> \n 0: Reproduce article numbers [default],\n 1: User customised input on Lund experimental data & \n 2: Test method with a simple method " << endl;
		cin >> run_type;

		switch(run_type) {
			case 0 : 
				cout << "You are running default" << endl;
				compute_random_chains(0);
				valid_input = kTRUE;
				break;
			case 1 : 
				cout << "User customised input: " << endl;
				compute_random_chains(1);
				valid_input = kTRUE;
				break;
			case 2 : 
				cout << "You are running the sanity check" << endl;
				compute_random_chains(2);
				valid_input = kTRUE;
				break;
			default : 
				cout << "Please insert a number, 0, 1 or 2" << endl;
		}
	}

	set_chains(run_type);

	//Calculates the background rates for alpha, escape in beam on and off. Number of implants are also calculated
	calculate_rates(run_type);

	//Calculates the number of expected random chains IN TOTAL for the specific chain/chains given as input to the program
	calculate_expected_nbr_random_chains();
	
			
}

RandomChains::~RandomChains() {
	delete this;
}

void RandomChains::set_chains(Int_t input) {
	cout << "Setting Chains ..."  << endl;

	cout << "Insert the number of decays in the chain, including a fission (if present)" << endl;
	cin >> chain_length;

	char ctemp;
	Int_t itemp;
	Int_t itemp2;

	for(Int_t i = 0; i < chain_length; i++) {
		cout << "Decay number:" << i << " What type is it? (alpha=a, escape=e or fission=f)" << endl;
		cin >> ctemp;
		decay_type.push_back(ctemp);
		cout << "Was the beam ON (type 1) or OFF (type 0) during this decay?" << endl;
		cin >> itemp2;
		beam_status.push_back(itemp2);
		cout << "What is the time span (s) that should be considered for this decay?" << endl;
		cin >> itemp;
		time_span.push_back(itemp);
	}

	cout << "Size of vector is = " << decay_type.size() << endl;
	for(vector<char>::iterator it = decay_type.begin(); it != decay_type.end(); it++) {
		cout << *it << endl;
	}
}


void RandomChains::compute_random_chains(Int_t run_type) {
	if(run_type == 2) {
		cout << "Generating test data ... " << endl;
		generate_test_data();
	}
	else {
		cout << "Reading experimental data ..." << endl;
		read_experimental_data();
	}
}

void RandomChains::generate_test_data() {
	cout << "Generating test data " << endl;
	
	for(Int_t k = 0; k < nbr_pixels; k++) {
		sprintf(ctitle,"Energy, Beam ON, pixel %d (TESTDATA)",k);
		sprintf(cname,"h_energy_pixel_reconstructed_ON_%d",k);
		h_energy_pixel_reconstructed_beam_on[k] = new TH1F(cname,ctitle, 4096, 0, 40.96);

		sprintf(ctitle,"Energy, Beam OFF, pixel %d (TESTDATA)",k);
		sprintf(cname,"h_energy_pixel_reconstructed_OFF_%d",k);
		h_energy_pixel_reconstructed_beam_off[k] = new TH1F(cname,ctitle, 4096, 0, 40.96);

		for(Int_t i = lower_limit_escapes; i < upper_limit_alphas; i++) {
			if(i < upper_limit_escapes) {
				h_energy_pixel_reconstructed_beam_on[k]->SetBinContent(i, 4);
				h_energy_pixel_reconstructed_beam_off[k]->SetBinContent(i, 3);
			}
			else if(i >= upper_limit_escapes && i < lower_limit_alphas) {
				h_energy_pixel_reconstructed_beam_on[k]->SetBinContent(i, 3);
				h_energy_pixel_reconstructed_beam_off[k]->SetBinContent(i, 2);
			}
			else {
				h_energy_pixel_reconstructed_beam_on[k]->SetBinContent(i, 2);
				h_energy_pixel_reconstructed_beam_off[k]->SetBinContent(i, 1);
			}
		}
		
		h_energy_pixel_reconstructed_beam_on[k]->SetBinContent(1500, 100);
		fissions_pixels[k] = 1;
	}


	/*
	TCanvas* theCanvas = new TCanvas("theCanvas");
	theCanvas->Divide(2,1);
	
	theCanvas->cd(1);
	h_energy_pixel_reconstructed_beam_on[100]->Draw();

	theCanvas->cd(2);
	h_energy_pixel_reconstructed_beam_off[100]->Draw();

	sprintf(ctitle,"Energy Total, Beam ON, pixel (TESTDATA)");
	sprintf(cname,"h_energy_pixel_reconstructed_ON_tot");
	*/

	h_energy_reconstructed_beam_on_tot = (TH1F*) h_energy_pixel_reconstructed_beam_on[0]->Clone(cname);

	//The total spectrum is only for plotting purposes, therefore the scaling
	h_energy_reconstructed_beam_on_tot->Scale(nbr_pixels);
	h_energy_reconstructed_beam_on_tot->SetTitle(ctitle);

	sprintf(ctitle,"Energy Total, Beam OFF, pixel (TESTDATA)");
	sprintf(cname,"h_energy_pixel_reconstructed_OFF__tot");

	h_energy_reconstructed_beam_off_tot = (TH1F*) h_energy_pixel_reconstructed_beam_off[0]->Clone(cname);

	h_energy_reconstructed_beam_off_tot->Scale(nbr_pixels);
	h_energy_reconstructed_beam_off_tot->SetTitle(ctitle);

	//theCanvas->Update();

}

void RandomChains::read_experimental_data() {
	cout << "Reading data " << endl;
}

void RandomChains::calculate_rates(Int_t run_type) {
	cout << "Calculating rates " << endl;
	if(run_type != 2) {
		cout << "Doing fission magic" << endl;
	}
	
	for(Int_t i = 0; i < decay_type.size(); i++) {
		rate_calc(decay_type.at(i), beam_status.at(i));
	}
}

void RandomChains::rate_calc(char type, Int_t beam) {
	cout << "Performing a rate calculation" << endl;

	TH1F** hist;
	if(beam) hist = h_energy_pixel_reconstructed_beam_on;
	else hist = h_energy_pixel_reconstructed_beam_off;
/*
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
*/
}

void RandomChains::calculate_expected_nbr_random_chains() {
	cout << "Calculating number of expected random chains " << endl;
}

int main() {

	RandomChains* RC = new RandomChains();
	
	RC->plot_spectra();

	return 0;
}

void run_main() {
	main();
}

void RandomChains::plot_spectra() {

#include "Style.code"

	TH1F* hist_on;
	TH1F* hist_off;
	Int_t input = -1;
	Bool_t valid_input = kFALSE;

	while(!valid_input) {
		cout << "What would you like to plot? \n" <<
			"-1: Total energy spectrum \n" <<
			"1: Energy spectrum per pixel" << endl;
		cin >> input;
		if(input < 0) {
			hist_on = (TH1F*) h_energy_reconstructed_beam_on_tot->Clone("Plot on");
			hist_off = (TH1F*) h_energy_reconstructed_beam_off_tot->Clone("Plot off");
			valid_input = kTRUE;
		}
		else if (input >= 0 && input < nbr_pixels) {
			hist_on = (TH1F*) h_energy_pixel_reconstructed_beam_on[input]->Clone("Plot on");
			hist_off = (TH1F*) h_energy_pixel_reconstructed_beam_off[input]->Clone("Plot off");
			valid_input = kTRUE;
		}

		else {
			cout << "Please insert a number" << endl;
		}
	}	
	//h_energy_pixel_off_tot->Scale(10);
	
	cout << "plotting pixel spectrum" << endl;
	TAxis *Xaxis = hist_on->GetXaxis();
	TAxis *Yaxis = hist_on->GetYaxis();

	cout << "plotting pixel spectrum" << endl;
	Yaxis->SetTitle("Counts per 10 keV");
	Xaxis->SetTitle("Energy (keV)");  
	//hist_on->SetTitle("");
	hist_on->SetStats(0);
	TGaxis::SetMaxDigits(5); 
	Xaxis->SetTitleOffset(1.2);
	Xaxis->SetTitleSize(0.05);
	Xaxis->SetRangeUser(5010,11990);      //(5,11995);       //(9510,11490); //4505,11995
	//Xaxis->SetNdivisions(505);  
	cout << "plotting pixel spectrum 1" << endl;

	Yaxis->SetTitleOffset(1.3);
	//Yaxis->SetNdivisions(606);
	//Yaxis->SetDecimals(1); //2 
	//Yaxis->SetRangeUser(0.001,24); 

	Xaxis->CenterTitle();
	Yaxis->CenterTitle();

	TCanvas *MyCanvas = new TCanvas("MyCanvas");
	MyCanvas->Divide(1,1);  

	cout << "plotting pixel spectrum 2" << endl;
	MyCanvas->cd(1); 

	//gPad->SetLogy(); 

	hist_on->SetLineWidth(2);
	hist_on->Draw();
	cout << "plotting pixel spectrum 3" << endl;

	hist_off->SetLineColor(8);
	hist_off->SetLineWidth(3);
	hist_off->SetLineStyle(1);
	hist_off->Draw("same");
	cout << "plotting pixel spectrum 4" << endl;
	hist_on->Draw("same");

}

