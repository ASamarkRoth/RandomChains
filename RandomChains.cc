#include <iostream>
#include <fstream>
#include <sstream>
#include "RandomChains.h"
#include <assert.h>
#include "TCanvas.h"
#include "TRint.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "math.h"
#include <typeinfo>

using namespace std;


RandomChains::RandomChains() {

	//Default value for number of pixels in the implantation detector
	//nbr_pixels = 1024;
	
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

	cout << "What type of run? <0/1/2> \n 0: Reproduce article numbers \n 1: User customised input on Lund experimental data \n 2: Test the program on trivial data (no input is required) " << endl;
	while(!valid_input) {
		cin >> run_type;

		switch(run_type) {
			case 0 : 
				cout << "You have chosen: \"Reproduce article numbers\"" << endl;
				prepare_data(0);
				valid_input = kTRUE;
				experiment_time = 1433000;
				break;
			case 1 : 
				cout << "You have chosen: \"user customised input\" " << endl;
				prepare_data(1);
				valid_input = kTRUE;
				experiment_time = 1433000;
				break;
			case 2 : 
				cout << "You have chosen: \"test run\"" << endl;
				prepare_data(2);
				valid_input = kTRUE;
				experiment_time = 1000000;
				break;
			default : 
				cout << "Please insert a number, 0, 1 or 2" << endl;
		}
	}

	//Sets the decay chain/chains
	set_chains(run_type);

	//Always dumps the set chains to the file "dump_input.txt" or in the case of a test "dump_test.txt"
	dump_input_to_file();
	
	//Calculate the number of implants per pixel
	calculate_implants();
	
	//Calculates the background rates for alpha, escape in beam on and off and fission.
	calculate_rates(run_type);

	//Calculates the number of expected random chains IN TOTAL for the specific chain/chains given as input to the program
	calculate_expected_nbr_random_chains();

	//Prints the result of the run in the terminal
	print_result();
			
}

RandomChains::~RandomChains() {
	delete this;
}

void RandomChains::set_chains(int input) {
	if(input == 2) {
		set_test_chains();
		return;
	}
	else if(input == 0) {
		set_article_chains();
		return;
	}
	cout << "Read decay chains from input file? <y/n>. If no, you will need to insert details of your chain step by step. "  << endl;
	char answer; 
	cin >> answer; 
	
	if(answer == 'y') {
		set_chains_from_input_file();
		return;
	}

	int l_temp;
	cout << "Insert the number of decays in the chain, including a fission (if present)" << endl;
	cin >> l_temp;
	chain_length.push_back(l_temp);

	char ctemp;
	int itemp;
	int itemp2;

	for(int i = 0; i < chain_length.front(); i++) {
		cout << "Decay number: " << i+1 << " What type is it? (alpha=a, escape=e or fission=f)" << endl;
		cin >> ctemp;
		decay_type.push_back(ctemp);
		cout << "Was the beam ON (type 1) or OFF (type 0) during this decay?" << endl;
		cin >> itemp2;
		beam_status.push_back(itemp2);
		cout << "What is the time span (s) that should be considered for this decay?" << endl;
		cin >> itemp;
		time_span.push_back(itemp);
	}


	/*
	for(vector<char>::iterator it = decay_type.begin(); it != decay_type.end(); it++) {
		cout << *it << endl;
	}
	*/
}

void RandomChains::set_test_chains() {
	//Do not modify these numbers!!!
	chain_length = {5};
	decay_type = {'a', 'e', 'a', 'e', 'f'};
	beam_status = {1, 0, 0, 1, 1};
	time_span = {1, 2, 3, 4, 5};
}

void RandomChains::set_article_chains() {
	//Do not modify these numbers!!!
	chain_length = {2, 2, 3, 3, 3, 3, 3};
	decay_type = {'a','f',  'e','f',  'a','a','f',  'a','a','f',  'a','a','f',  'a','a','f',  'e','e','f'};
	beam_status = {0,0,  0,0,  1,0,0,  1,0,0,  0,0,0,  0,0,0,  0,1,0};
	time_span = {2,10,  2,10,  2,10,50,  2,10,50,  2,10,50,  2,10,50, 2,10,50};
}


void RandomChains::prepare_data(int run_type) {
	if(run_type == 2) {
		generate_test_data();
	}
	else {
		cout << "Reading experimental data ..." << endl;
		read_experimental_data();
	}
}

void RandomChains::generate_test_data() {
	cout << "Generating test data " << endl;
	
	for(int k = 0; k < nbr_pixels; k++) {
		sprintf(ctitle,"Energy, Beam ON, pixel %d (TESTDATA)",k);
		sprintf(cname,"h_energy_pixel_reconstructed_ON_%d",k);
		h_energy_pixel_reconstructed_beam_on[k] = new TH1F(cname,ctitle, 4096, 0, 40.96);

		sprintf(ctitle,"Energy, Beam OFF, pixel %d (TESTDATA)",k);
		sprintf(cname,"h_energy_pixel_reconstructed_OFF_%d",k);
		h_energy_pixel_reconstructed_beam_off[k] = new TH1F(cname,ctitle, 4096, 0, 40.96);

		//Setting the values to insert in the test spectra:
		eon = 4;
		eoff = 3;
		non = 3;
		noff = 2;
		aon = 2;
		aoff = 1;
		imps = 100;
		fissions = 2;

		for(int i = lower_limit_escapes; i < upper_limit_alphas; i++) {
			if(i < upper_limit_escapes) {
				h_energy_pixel_reconstructed_beam_on[k]->SetBinContent(i, eon);
				h_energy_pixel_reconstructed_beam_off[k]->SetBinContent(i, eoff);
			}
			else if(i >= upper_limit_escapes && i < lower_limit_alphas) {
				h_energy_pixel_reconstructed_beam_on[k]->SetBinContent(i, non);
				h_energy_pixel_reconstructed_beam_off[k]->SetBinContent(i, noff);
			}
			else {
				h_energy_pixel_reconstructed_beam_on[k]->SetBinContent(i, aon);
				h_energy_pixel_reconstructed_beam_off[k]->SetBinContent(i, aoff);
			}
		}
		
		h_energy_pixel_reconstructed_beam_on[k]->SetBinContent(1500, imps);
		fissions_pixels[k] = fissions;
	}


	/*
	TCanvas* theCanvas = new TCanvas("theCanvas");
	theCanvas->Divide(2,1);
	
	theCanvas->cd(1);
	h_energy_pixel_reconstructed_beam_on[100]->Draw();

	theCanvas->cd(2);
	h_energy_pixel_reconstructed_beam_off[100]->Draw();

	*/
	sprintf(ctitle,"Energy Total, Beam ON, (TESTDATA)");
	sprintf(cname,"h_energy_pixel_reconstructed_ON_tot");

	h_energy_reconstructed_beam_on_tot = (TH1F*) h_energy_pixel_reconstructed_beam_on[0]->Clone(cname);

	//The total spectrum is only for plotting purposes, therefore the scaling
	h_energy_reconstructed_beam_on_tot->Scale(nbr_pixels);
	h_energy_reconstructed_beam_on_tot->SetTitle(ctitle);

	sprintf(ctitle,"Energy Total, Beam OFF, (TESTDATA)");
	sprintf(cname,"h_energy_pixel_reconstructed_OFF_tot");

	h_energy_reconstructed_beam_off_tot = (TH1F*) h_energy_pixel_reconstructed_beam_off[0]->Clone(cname);

	h_energy_reconstructed_beam_off_tot->Scale(nbr_pixels);
	h_energy_reconstructed_beam_off_tot->SetTitle(ctitle);

	//theCanvas->Update();

}

void RandomChains::read_experimental_data() {

	//Prepare/read in root spectra
	TFile  *file = new TFile("Spectra_test.root");

	h_energy_reconstructed_beam_off_tot = (TH1F*)file->Get("h_energy_pixel_recoff_tot"); 
	h_energy_reconstructed_beam_on_tot = (TH1F*)file->Get("h_energy_pixel_recon_tot"); 

	//Read in spectra, cont.
	for(int i = 0; i < 1024; i++){ 
	sprintf(ctitle,"recon/h_energy_pixel_recon_%d",i);
	h_energy_pixel_reconstructed_beam_on[i] = (TH1F*)file->Get(ctitle); 
	} 

	for(int i = 0; i < 1024; i++){ 
	sprintf(ctitle,"recoff/h_energy_pixel_recoff_%d",i);
	h_energy_pixel_reconstructed_beam_off[i] = (TH1F*)file->Get(ctitle); 
	} 

	int temp_value;
	int nbr_of_fissions = 0;
	ifstream input("Fission-BeamOff-Pixel_modified.txt",ios::in);
	input>>temp_value;
	while (input){
		//if(fissions_pixels[temp_value] == 1) cout<<temp_value<<endl;
		fissions_pixels[temp_value] += 1;
		nbr_of_fissions++;
		input>>temp_value;
	}
	//cout << "Total number of fissions are: " << nbr_of_fissions << endl;
	for(int j = 0; j < nbr_pixels; j++) {
		2;
		//cout << "Pixel " << j << " = " << fissions_pixels[j] << endl;
	}
	//If the number of fissions in a pixel is 0 then it is set to the average over the complete implantation detector. 
	for(int i = 0; i < nbr_pixels; i++){
		if(fissions_pixels[i] == 0){
			fissions_pixels[i] = (double)nbr_of_fissions/nbr_pixels;    
		} 
	}

}

void RandomChains::calculate_implants() {
	for(int i = 0; i < nbr_pixels; i++) {
		int acc_counts = 0; 
		for(int k = lower_limit_implants; k < upper_limit_implants; k++) {
			acc_counts += h_energy_pixel_reconstructed_beam_on[i]->GetBinContent(k);
		}
		nbr_implants[i] = acc_counts;
	}


}

void RandomChains::calculate_rates(int run_type) {
	cout << "Calculating rates " << endl;

	for(unsigned int i = 0; i < decay_type.size(); i++) {
		rate_calc(decay_type.at(i), beam_status.at(i));
	}

}

void RandomChains::rate_calc(char type, int beam) {

	TH1F** hist;
	if(beam) hist = h_energy_pixel_reconstructed_beam_on;
	else hist = h_energy_pixel_reconstructed_beam_off;

	//double rate_temp[nbr_pixels];
	array<long double, nbr_pixels> rate_temp;
	
	int lower_limit, upper_limit;
	if(type == 'a') {
		lower_limit = lower_limit_alphas;
		upper_limit = upper_limit_alphas;
	}
	else if(type == 'e') {
		lower_limit = lower_limit_escapes;
		upper_limit = upper_limit_escapes;
	}
	else if(type == 'f') {
		for(int i = 0; i < nbr_pixels; i++) {
			rate_temp[i] = (long double)fissions_pixels[i]/experiment_time;
		}
		rate.push_back(rate_temp);
		return;
	}
	else {
		cout << "Please input correct decay types, i.e. 'a', 'e' or 'f' " << endl;
		return;
	}

	for(int i = 0; i < nbr_pixels; i++) {
		int acc_counts = 0;
		for(int k = lower_limit; k < upper_limit; k++) {
			acc_counts += hist[i]->GetBinContent(k);
		}
		rate_temp[i] = (long double) acc_counts/experiment_time;
	}
	/*
	for(auto j = 0; j < nbr_pixels; j++) {
		cout << "Rate at j = " << j << " is " << rate_temp[j] << endl;
	}
	*/
	rate.push_back(rate_temp);
/*
//Calculate rate of a certain type of events, between lower_limit and upper_limit
double rate_calc(int pixel_number, int lower_limit, int upper_limit, TH1F **h_energy, double time){
  int temporary = 0;
  long double rate = 0;
      temporary = 0;
      for (int i = lower_limit; i < upper_limit; i++)
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

//Calculate number of expected random chains in TOTAL
void RandomChains::calculate_expected_nbr_random_chains() {
	cout << "Calculating expected number of random chains " << endl;
	int offset = 0;
	for(unsigned int j = 0; j < chain_length.size(); j++) {
		double randoms_in_pixel[nbr_pixels] = {0};
		for(int l = offset; l < offset+chain_length.at(j); l++) {
			cout << "l = " << l << endl;
			cout << "decay type = " << decay_type.at(l) << endl;
			cout << "beam status = " << beam_status.at(l) << endl;
			cout << "time span = " << time_span.at(l) << endl;
			for(int i = 0; i < nbr_pixels; i++) {
				if(l == offset) randoms_in_pixel[i] = nbr_implants[i];
				randoms_in_pixel[i] *= (1-Poisson_pmf(0,rate.at(l)[i]*time_span.at(l)));
				/*
				if(rate.at(l)[i] == 0) {
					cout << "Rate in pixel " << i << " = " << rate.at(l)[i] << endl;
					cout << "Poisson = " << Poisson_pmf(0,rate.at(l)[i]*time_span.at(l)) << endl;
				}
				*/
				cout << "rate in pixel " << i << " is  = " << (rate.at(l))[i] << endl;
				cout << "Poisson = " << Poisson_pmf(0,rate.at(l)[i]*time_span.at(l)) << endl;
			}
		}

		//Sum the number of randoms in all pixels to get the TOTAL number of random chains
		double random_chains_temp = 0; 
		for(int i = 0; i < nbr_pixels; i++) {
			random_chains_temp += randoms_in_pixel[i];
		}

		nbr_expected_random_chains.push_back(random_chains_temp);
		offset += chain_length.at(j);
	}
}


int main() {

	RandomChains* RC = new RandomChains();
	return 0;
}

void run_main() {
	main();
}



void RandomChains::plot_spectra() {

#include "Style.code"

	TH1F* hist_on;
	TH1F* hist_off;
	int input = -1;
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
	Xaxis->SetTitle("Energy (MeV)");  
	//hist_on->SetTitle("");
	hist_on->SetStats(0);
	TGaxis::SetMaxDigits(5); 
	Xaxis->SetTitleOffset(1.2);
	Xaxis->SetTitleSize(0.05);
	//Xaxis->SetRangeUser(5010,11990);      //(5,11995);       //(9510,11490); //4505,11995
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

void RandomChains::print_result() {
	cout << "**************************************************" << endl;
	if(run_type == 2) cout << "These are the result of the TEST run: " << endl;
	else cout << "These are the result of the run: " << endl;

	cout << "The total number of expected random chains are: " << endl;

	for(unsigned int j = 0; j < nbr_expected_random_chains.size(); j++) {
		cout << "For chain " << j+1 << ": " << nbr_expected_random_chains.at(j) << endl;
	}
	if(run_type == 2) print_test_result();
}

void RandomChains::print_test_result() {
	cout << "*************************************************" << endl;
	cout << "NOW TEST CALCULATION" << endl;
	cout << "By construction, the rate in every pixel will be the same" << endl;
	cout << "In the test case, all time spans are set to one second.  " << endl;

	double rate_escapes_off = (upper_limit_escapes-lower_limit_escapes)*eoff/experiment_time;
	double rate_alphas_on = (upper_limit_alphas-lower_limit_alphas)*aon/experiment_time;
	
	double rate_escapes_on = (upper_limit_escapes-lower_limit_escapes)*eon/experiment_time;
	double rate_alphas_off = (upper_limit_alphas-lower_limit_alphas)*aoff/experiment_time;

	int nbr_imps = imps;

	double fission_rate = fissions/experiment_time;

	//double test_randoms = nbr_imps*(1-Poisson_pmf(0,rate_alphas_on))*(1-Poisson_pmf(0, rate_escapes_off))*(1-Poisson_pmf(0, rate_escapes_on))*(1-Poisson_pmf(0, fission_rate)) * nbr_pixels;
	double test_randoms = nbr_imps*(1-Poisson_pmf(0,rate_alphas_on*time_span.at(0)))*(1-Poisson_pmf(0, rate_escapes_off*time_span.at(1)))*(1-Poisson_pmf(0, rate_alphas_off*time_span.at(2)))*(1-Poisson_pmf(0, rate_escapes_on*time_span.at(3)))*(1-Poisson_pmf(0, fission_rate*time_span.at(4))) * nbr_pixels;
	cout << "test_randoms = nbr_imps*(1-Poisson_pmf(0,rate_alphas_on))*(1-Poisson_pmf(0, rate_escapes_off))*(1-Poisson_pmf(0, rate_alphas_off))*(1-Poisson_pmf(0, rate_escapes_on))*(1-Poisson_pmf(0, fission_rate)) * nbr_pixels;" << endl;

	cout << "The CALCULATED total number of random chains with the test data are: " << test_randoms << endl;
	
}

void RandomChains::dump_input_to_file() {
	char output[64];
	if(run_type == 2) {
		sprintf(output, "dump_test.txt");
	}
	else if(run_type == 0) sprintf(output, "dump_article.txt");
	else sprintf(output, "dump_input.txt");
	char out[64];
	if(run_type == 0) {
		sprintf(out, "The following input was given ... ");
		cout << out << endl;
	}
	sprintf(out, "File %s was written ... ", output);
	ofstream dump;
	dump.open(output);
	dump << "Lines starting with a '#' indicates the start of a new chain. OBS, the two first lines are not read in." << endl;
	dump << "Type (alpha=a, escape=e and fission=f) 	Beam ON (=1) or OFF (=0)	Time span (s) \n";
	if(run_type == 0) {
		cout << "Lines starting with a '#' indicates the start of a new chain. OBS, the two first lines are not read in." << endl;
		cout << "Type (alpha=a, escape=e and fission=f) 	Beam ON (=1) or OFF (=0)	Time span (s) \n";
	}
	int offset = 0;
	for(unsigned int k = 0; k < chain_length.size(); k++) {
		dump << "#" << chain_length.at(k) << endl;
		if(run_type == 0) cout << "#" << chain_length.at(k) << endl;
		for(int j = offset; j < offset + chain_length.at(k); j++) {
			dump << decay_type.at(j) << " " << beam_status.at(j) << " " << time_span.at(j) << endl;
			if(run_type == 0) cout << decay_type.at(j) << " " << beam_status.at(j) << " " << time_span.at(j) << endl;
		}
		offset += chain_length.at(k);
	}
	dump.close();

	cout << out << endl;
}

void RandomChains::set_chains_from_input_file() {
	string filename;
	cout << "Enter full name of input file name: (file \"dump_input.txt\" should be available): ";
	cin >> filename;
	ifstream input_chains(filename,ios::in);
	if(!input_chains) cout << "Could not find file" << endl;
	int beam;
	double time;
	char type;
	string str;
	stringstream ss;
	int counter = 0;
	cout << "The following was read in: " << endl;
	while(getline(input_chains, str)) {
		cout << str << endl;
		if(chain_length.size() > 0) cout << "chain length = " << chain_length.at(0) << endl;
		if(counter < 2) {
			counter++;
			continue;
		}
		if(str[0] == '#') {
			string number = string(str.begin()+1, str.end());
			chain_length.push_back(stoi(number));
		}
		else {
			ss = stringstream(str);
			ss >> type >> beam >> time;
			decay_type.push_back(type);
			beam_status.push_back(beam);
			time_span.push_back(time);
		}
	}
	
	/*
	cout << "length of chain lengtj = " << chain_length.size() << endl;
	cout << "chain length = " << chain_length.at(0) << endl;
	cout << "Decay type =, beam status, time " << endl;
	for(int j = 0; j < decay_type.size(); j++) {
		cout << decay_type.at(j) << " ";
		cout << beam_status.at(j) << " ";
		cout << time_span.at(j) << " ";
	}
	*/
			
}

/* Mathematical functions */

//Poisson probability mass function
double Poisson_pmf(int nbr_to_observe, double expected_value) {
	double prob;
	prob = (exp(-expected_value)*pow(expected_value,nbr_to_observe))/(factorial(nbr_to_observe));
	return prob;
}

int factorial(int k){
	int ret;
	ret = 1;
	for (int j = 1; j <= k; j++){
	ret = ret*j;
	}
	return ret;
}

void insert_blank() {
	cout << "\n";
}
