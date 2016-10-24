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

	while(!valid_input) {
		cout << "What type of run? <0/1/2> \n 0: Reproduce article numbers \n 1: User customised input on Lund experimental data \n 2: Test the program on trivial data (no input is required) " << endl;
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

	for(int i = 0; i < nbr_pixels; i++) {
		fissions_pixels[i] = 0; 
	}

	int temp_value;
	int nbr_of_fissions = 0;
	ifstream input("Fission-BeamOff-Pixel_modified.txt",ios::in);
	input>>temp_value;
	while (input){
		//if(fissions_pixels[temp_value] == 1) cout<<temp_value<<endl;
		cout << "Temp_value = " << temp_value << endl;
		cout << "Temp_value id = " << typeid(temp_value).name() << endl;
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
	for(int i = 0; i < nbr_pixels; i++)  {
		cout << "This is where fissions are set for pixel i = " << i << " " << fissions_pixels[i] << endl;
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
				cout << "rate in pixel " << i << " is  = " << (rate.at(l))[i] << endl;
				cout << "Poisson = " << Poisson_pmf(0,rate.at(l)[i]*time_span.at(l)) << endl;
				*/
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
	RC->compare_rates();
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
 long double* u_main() {


  //Prepare/read in root spectra
  TFile  *file = new TFile("Spectra.root");

  TH1F *h_energy_pixel_off_tot = (TH1F*)file->Get("h_energy_pixel_off_tot");
  TH1F *h_energy_pixel_recoff_tot = (TH1F*)file->Get("h_energy_pixel_recoff_tot"); 
  TH1F *h_energy_pixel_on_tot = (TH1F*)file->Get("h_energy_pixel_on_tot"); 
  TH1F *h_energy_pixel_recon_tot = (TH1F*)file->Get("h_energy_pixel_recon_tot"); 

  //h_energy_pixel_recoff_tot->Draw();


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

  /*
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
  */
	lookslike[0] = 'e';
	a1 = "on";

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

  //Calc rates for FISSION
  for(int i = 0; i < 1024; i++){
    rate_f[i] = fission_array[i]/exp_time;
  }

	return rate_f;

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



void RandomChains::compare_rates() {
	long double* u_rate = u_main();
	long double compare[1024];
	for(int i = 0; i < 1024; i++) {
		compare[i] = rate.at(0)[i] - u_rate[i];
		cout << "Compare at i = " << i << " is " << compare[i] << endl;
	}
}






