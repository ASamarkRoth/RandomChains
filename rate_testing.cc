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
#include "TLegend.h"
#include "math.h"
#include <typeinfo>

using namespace std;

//Functions 
double probability(int to_obs, float expected);
int fact(Int_t k);
double rate_calc(int pixel_number, int lower_limit, int upper_limit, TH1F **h_energy, double time);
long double* u_main(int step);

void compare_rates(long double* main, long double* uf, bool print=false);

/*
Member functions and other functions of the class RandomChains are defined in this file. 

The class RandomChains handles data and computes the number of expected random chains due to random fluctuations in the background for specific decay chains. For the description of the method, see U. Forsberg et. al. / Nuclear Physics A 953 (2016) 117-138.

The user can choose between three different run types: 
	0: "Reproduce article numbers", i.e. reproduce the numbers presented in U. Forsberg et. al. / Nuclear Physics A 953 (2016) 117-138.
	1: "User customised input". This gives the user the option to:
		Step by step input the details of one chain or 
		Set the decay chain (or chains) specific from a file (with a special format, see dump_test.txt.
	2: "Test run". With this option the program is tested on trivial data. The obtained number from the program is compared to a value obtained through the calculation of a simple formula. 

Options "0" and "1" assumes the experimental data from the E115 experiment conducted 2012 lead by the Lund Nuclear Structure Group (see D. Rudolph et. al / Phys. Rev. Lett. 111, 112502 (2013). The pixel by pixel reconstructed spectra (reconstructed energies of escaping alpha particles to one of the box detectors) is read in from the .root file "Spectra.root". The pixels which detected fissions are read in from the ascii file "pixels_with_fissions.txt".

The number of pixels in the implantation detector is central in the analysis. In the header file it is set to 1024 and in the current state of the program it should not be modified. 

From a user perspective, what defines a decay chain is the following characteristics:
	chain_length = x, where x is the length of the chain.
	decay_type = 'a'=alpha, 'e'=escape and 'f'=fission
	beam_status = 1=ON and 0=OFF
	time_span = t, where t is the length of the time window during which the decay is accepted.
These characteristics can be set by the user (if run type == 1) either via an input file or step by step terminal inputs. 

The decay types are defined with the detected energies E as:
	'a'=alpha	lower_limit_alphas <= E < upper_limit_alphas
	'e'=escape	lower_limit_escapes <= E < upper_limit_escapes
Energies assigned as implants (only for beam ON) are defined as: 
	implants	lower_limit_implants <= E < upper_limit_implants

Above, lower_limit_* and upper_limit_* are the bins in the spectra which for the case of the experimental data has a width of 10 keV. These limits are set at the very beginning of the RandomChains class constructor. 

An input file after a run is automatically generated. If run type == 0 or 2 file "dump_article.txt" and "dump_test.txt" are generated, respectively. Otherwise the file "dump_input.txt" is generated. 

An example of the input file "dump_test.txt":
	Lines starting with a '#' indicates the start of a new chain. OBS, the two first lines are not read in.
	Type (alpha=a, escape=e and fission=f) 	Beam ON (=1) or OFF (=0)	Time span (s) 
	#5
	a 1 1
	e 0 2
	a 0 3
	e 1 4
	f 1 5

The first two lines are just a description. 
Lines starting with a "#" indicates the start of a new decay chain and the subsequent number indicates the chain_length of the specific chain. 
The lines with three colums (space separated) indicate the characteristics of the decay as described in the second line of the file.

The user can modify "dump_test.txt" as he/she likes (it is recommended to also modify the file name) and provide the file as an input with run type == 1.
An alternative to create a customised input is to give a step by step input through run type == 1 and after the run save the automatically generated file "dump_input.txt" with a new name.  

The result of a run is presented with an output to the terminal window. 

With ROOT the user can plot the data for either a specific pixel or the total spectrum with beam ON/OFF data overlayed. This can be achieved with the function "plot_spectra()" which requests an input for pixel number. 

*/


RandomChains::RandomChains() {
/* The constructor of class RandomChains. Here a complete run is controlled and executed. The following is done:
	1. Lower and upper limits for the decay types and implants are set.
	2. The run_type is given by the user.
	3. The experimental data is read in or the test data is generated.
	4. The chain/chains characteristics are set.
	5. The input data is dumped to a file.
	6. The number of implants per pixel is calculated.
	7. The background rates for alphas, escapes for beam ON and OFF and fission are calculated for every pixel.
	8. The TOTAL number of expected random chains due to random fluctuations in the background are calculated for the specific chain/chains given as input to the program.
	9. The result is printed in the terminal window. 
	
	*/

	// 1. Lower and upper limits for the decay types and implants are set.
	//Interval for accepted superheavy nuclei alpha decays (10*keV)
	lower_limit_alphas = 900; upper_limit_alphas = 1100;

	//Interval for energy deposit of alphas that escapes the implantation detector (10*keV)
	lower_limit_escapes = 0; upper_limit_escapes = 400;

	//Interval for implanted nuclei with beam = ON (10*keV)
	lower_limit_implants = 1100; upper_limit_implants = 1800;

	//2. The run_type is given by the user.
	Bool_t valid_input = kFALSE;

	while(!valid_input) {
		cout << "What type of run? <0/1/2> \n" << 
			"	0: Reproduce article numbers (no input required) \n "<<
			"	1: User customised input on Lund experimental data \n" << 
			"	2: Test the program on trivial data (no input required) " << endl;
		cin >> run_type;

		switch(run_type) {
		//3. The experimental data is read in or the test data is generated with the method "prepare_data()".
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

	//4. The chain/chains characteristics are set.
	set_chains(run_type);

	//5. The input data is dumped to a file.
	dump_input_to_file();

	//6. The number of implants per pixel is calculated.
	calculate_implants();

	//7. The background rates for alphas, escapes for beam ON and OFF and fission are calculated for every pixel.
	calculate_rates();

	//8. The TOTAL number of expected random chains due to random fluctuations in the background are calculated for the specific chain/chains given as input to the program.
	calculate_expected_nbr_random_chains();

	//9. The result is printed in the terminal window. 
	print_result();
		
}

RandomChains::~RandomChains() {
	/*The destructor of RandomChains. Object is deleted.*/
	delete this;
}

void RandomChains::prepare_data(int run_type) {
	/*3. The experimental data is read in or the test data is generated with the method. 
	Input argument int run_type (0, 1 or 2) is given by the user in the constructor.

	*/
		
	if(run_type == 2) {
		cout << "Generating test data " << endl;
		generate_test_data();
	}
	else {
		cout << "Reading experimental data ..." << endl;
		read_experimental_data();
	}
}

void RandomChains::generate_test_data() {
	/* The test data is generated.
	Through this method the following member data are initialised:
		TH1F* h_energy_reconstructed_beam_off_tot;
		TH1F* h_energy_reconstructed_beam_on_tot;

		TH1F* h_energy_pixel_reconstructed_beam_off[nbr_pixels];
		TH1F* h_energy_pixel_reconstructed_beam_on[nbr_pixels];
		double fissions_pixels[nbr_pixels];

		int eon;
		int eoff;
		int non;
		int noff;
		int aon;
		int aoff;
		int imps;
		int fissions;
	*/

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

}

void RandomChains::read_experimental_data() {
	/* The experimental data is read in from the files "Spectra.root" and "pixels_with_fissions.txt".
	Through this method the following member data are initialised:
		TH1F* h_energy_reconstructed_beam_off_tot;
		TH1F* h_energy_reconstructed_beam_on_tot;

		TH1F* h_energy_pixel_reconstructed_beam_off[nbr_pixels];
		TH1F* h_energy_pixel_reconstructed_beam_on[nbr_pixels];
		double fissions_pixels[nbr_pixels];

	*/

	

	//Prepare/read in root spectra
	TFile  *file = new TFile("Spectra.root");

	h_energy_reconstructed_beam_off_tot = (TH1F*)file->Get("h_energy_pixel_recoff_tot"); 
	h_energy_reconstructed_beam_on_tot = (TH1F*)file->Get("h_energy_pixel_recon_tot"); 

	for(int i = 0; i < nbr_pixels; i++){ 
	sprintf(ctitle,"recon/h_energy_pixel_recon_%d",i);
	h_energy_pixel_reconstructed_beam_on[i] = (TH1F*)file->Get(ctitle); 
	} 

	for(int i = 0; i < nbr_pixels; i++){ 
	sprintf(ctitle,"recoff/h_energy_pixel_recoff_%d",i);
	h_energy_pixel_reconstructed_beam_off[i] = (TH1F*)file->Get(ctitle); 
	} 

	for(int i = 0; i < nbr_pixels; i++){ 
	sprintf(ctitle,"on/h_energy_pixel_on_%d",i);
	h_energy_pixel_beam_on[i] = (TH1F*)file->Get(ctitle); 
	} 

	for(int i = 0; i < nbr_pixels; i++) {
		fissions_pixels[i] = 0; 
	}
	
	//The fission data are read in here and treated differently.
	int temp_value;
	int nbr_of_fissions = 0;
	ifstream input("pixels_with_fissions.txt",ios::in);
		input>>temp_value;
		while (input){
			cout << "Temp_value = " << temp_value << endl;
			//cout << "Temp_value id = " << typeid(temp_value).name() << endl;
			fissions_pixels[temp_value] += 1;
			nbr_of_fissions++;
			input>>temp_value;
		}
		cout << "Total number of fissions are: " << nbr_of_fissions << endl;

		//If the number of fissions in a pixel is 0 then it is set to the average over the complete implantation detector. 
		for(int i = 0; i < nbr_pixels; i++){
			if(fissions_pixels[i] == 0){
				fissions_pixels[i] = (double)nbr_of_fissions/nbr_pixels;    
			} 
		}
		cout << "HERE IS COMPARISON AT STEP 1, FISSION" << endl;
		compare_rates(fissions_pixels, u_main(1));

}

void RandomChains::set_chains(int input) {
	/* The chain/chains characteristics are set. 
	The input argument "int input" is the run_type and which decay chains are to be set are determined on the basis of this value.

	The following member data is initialised:
		vector<int> chain_length;
		vector<int> beam_status;
		vector<char> decay_type;
		vector<float> time_span;
	*/


	if(input == 2) {
		//In this method the test chains are set.
		set_test_chains();
		return;
	}
	else if(input == 0) {
		//In this method the article chains are set. 
		set_article_chains();
		return;
	}


	//First the user decides if he/she wants to read data from an input file
	cout << "Read decay chains from input file? <y/n>. If no, you will need to insert details of your chain step by step. "  << endl;
	char answer; 
	cin >> answer; 

	if(answer == 'y') {
		//This method handles how the chains from an input file are set.
		set_chains_from_input_file();
		return;
	}

	//The following is the step by step input that can be chosen to be provided by the user.
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

}

void RandomChains::set_test_chains() {
	/*This method defines the test chains characteristics*/
	//Do not modify these numbers!!!
	chain_length = {5};
	decay_type = {'a', 'e', 'a', 'e', 'f'};
	beam_status = {1, 0, 0, 1, 1};
	time_span = {1, 2, 3, 4, 5};
}

void RandomChains::set_article_chains() {
	/*This method defines the article chains characteristics*/
	//Do not modify these numbers!!!
	chain_length = {2, 2, 3, 3, 3, 3, 3};
	decay_type = {'a','f',  'e','f',  'a','a','f',  'a','a','f',  'a','a','f',  'a','a','f',  'e','e','f'};
	beam_status = {0,0,  0,0,  1,0,0,  1,0,0,  0,0,0,  0,0,0,  0,1,0};
	time_span = {2,10,  2,10,  2,10,50,  2,10,50,  2,10,50,  2,10,50, 2,10,50};
}

void RandomChains::dump_input_to_file() {
	/*This method dumps the set chains to a file. If the test is run the file name is "dump_test.txt" and if the reproduce article numbers is run the file name is "dump_article.txt". Otherwise the file name is "dump_input.txt"
	*/

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
	/*This method is called when the user chooses to read in input from a file.

	Through this method the following member data is initialised:
		vector<int> chain_length;
		vector<int> beam_status;
		vector<char> decay_type;
		vector<float> time_span;
	*/

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


void RandomChains::calculate_implants() {
	/* The number of implants per pixel is calculated. 
	
	Through this method the following member data is initialised:
		int nbr_implants[nbr_pixels];
	
	*/


	TH1F** hist;
	/*
	if(h_energy_pixel_beam_on!=NULL) hist = h_energy_pixel_beam_on;
	else hist = h_energy_pixel_reconstructed_beam_on;
	*/

	hist = h_energy_pixel_beam_on;
		
	for(int i = 0; i < nbr_pixels; i++) {
		int acc_counts = 0; 
		for(int k = lower_limit_implants; k < upper_limit_implants; k++) {
			acc_counts += hist[i]->GetBinContent(k);
		}
		nbr_implants[i] = acc_counts;
	}

	long double imp_temp[1024];
	for(int j = 0; j < 1024; j++) {
		imp_temp[j] = nbr_implants[j];
	}

	cout << "HERE IS COMPARISON AT STEP 2" << endl;
	compare_rates(imp_temp, u_main(2));


}

void RandomChains::calculate_rates() {
	/*This method calculates the rates in every pixel for the specific decay types, one decay at a time.

	Through this method the following member data is initialised:
		vector<array<long double, nbr_pixels>> rate;
	
	*/

	cout << "Calculating rates " << endl;

	for(unsigned int i = 0; i < decay_type.size(); i++) {
		rate_calc(decay_type.at(i), beam_status.at(i));
	}

	cout << "HERE COMPARISON STEP 3" << endl;
	long double rate_temp1[nbr_pixels];
	for(int j=0; j < nbr_pixels; j++) {
		rate_temp1[j] = rate.at(0)[j];
	}
	compare_rates(rate_temp1, u_main(3));

	
	long double* rate_temp2;
	rate_temp2 = u_main(3);
	array<long double, nbr_pixels> rate_temp;
	for(int j = 0; j<nbr_pixels; j++) {
		rate_temp[j] = rate_temp2[j];
	}
	rate.at(0) = rate_temp;

	cout << "\n Here comparison with modified shit !" << endl;
	compare_rates(rate_temp2, u_main(3));

	cout << "HERE COMPARISON STEP 4" << endl;
	long double rate_tempf[nbr_pixels];
	for(int j=0; j < nbr_pixels; j++) {
		rate_tempf[j] = rate.at(1)[j];
	}

	compare_rates(rate_tempf, u_main(4));

}

void RandomChains::rate_calc(char type, int beam) {
	/*This method calculates the rate per pixel depending on the type of decay and the beam status.
	Input arguments are:
	"char type"=decay type, i.e. 'a', 'e' or 'f'.
	"int beam"=beam status, i.e. 1/0

	Through this method the following member data is modified:
		vector<array<long double, nbr_pixels>> rate;
		
	*/


	//Based on the beam status the spectrum is determined
	TH1F** hist;
	if(beam) hist = h_energy_pixel_reconstructed_beam_on;
	else hist = h_energy_pixel_reconstructed_beam_off;

	array<long double, nbr_pixels> rate_temp;
	
	//Based on the decay type, the upper and lower bin limits are set
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

	//The rates for every pixel is calculated
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
}

void RandomChains::calculate_expected_nbr_random_chains() {
	/*This method calculates the TOTAL number of expected random chains for the input decay chain/chains.

	Through this method the following member data is modified:
		vector<double> nbr_expected_random_chains;

	*/

	cout << "Calculating expected number of random chains " << endl;
	int offset = 0;

	//looping decay chains
	for(unsigned int j = 0; j < chain_length.size(); j++) {
		long double randoms_in_pixel[nbr_pixels] = {0};

		//looping decays
		/*
		for(int l = offset; l < offset+chain_length.at(j); l++) {
			cout << "l = " << l << endl;
			cout << "decay type = " << decay_type.at(l) << endl;
			cout << "beam status = " << beam_status.at(l) << endl;
			cout << "time span = " << time_span.at(l) << endl;

			//looping pixels
			for(int i = 0; i < nbr_pixels; i++) {
				if(l == offset) randoms_in_pixel[i] = nbr_implants[i];
				randoms_in_pixel[i] *= (1-Poisson_pmf(0,rate.at(l)[i]*time_span.at(l)));
				if(rate.at(l)[i] == 0) {
					cout << "Rate in pixel " << i << " = " << rate.at(l)[i] << endl;
					cout << "Poisson = " << Poisson_pmf(0,rate.at(l)[i]*time_span.at(l)) << endl;
				}
				cout << "rate in pixel " << i << " is  = " << (rate.at(l))[i] << endl;
				cout << "Poisson = " << Poisson_pmf(0,rate.at(l)[i]*time_span.at(l)) << endl;
			}
		}
		*/

		cout << "Timespan(0) = " << time_span.at(0) << " and Timespan(1) = " << time_span.at(1) << endl;
		//Doing it as UF:
		for(int i = 0; i < nbr_pixels; i++) {
			randoms_in_pixel[i] = nbr_implants[i]*(1-Poisson_pmf(0,rate.at(0)[i]*time_span.at(0)))*(1-Poisson_pmf(0,rate.at(1)[i]*time_span.at(1)));
		}

		long double* ar_temp;
		ar_temp = u_main(3);

		for(int i = 0; i < nbr_pixels; i++) {
			2;
			//cout << "Comp. rate.at for pixel " << i << " is = " << rate.at(0)[i] - ar_temp[i] << endl;
		}

		//cout << "HERE IS COMPARISON AT STEP 5 " << j << endl;
		//compare_rates(randoms_in_pixel, u_main(5), true);


		//Sum the number of randoms in all pixels to get the TOTAL number of random chains
		double random_chains_temp = 0; 
		for(int i = 0; i < nbr_pixels; i++) {
			random_chains_temp += randoms_in_pixel[i];
		}

		nbr_expected_random_chains.push_back(random_chains_temp);
		offset += chain_length.at(j);
	}
}

void RandomChains::print_result() {
	/*This is the method that is invoked at the end of the constructor and it presents the result of the run in the terminal window. If the test was run another member function is called for further output.
	*/

	cout << "**************************************************" << endl;
	if(run_type == 2) cout << "These are the result of the TEST run: " << endl;
	else cout << "These are the result of the run: " << endl;

	cout << "The total number of expected random chains are: " << endl;

	for(unsigned int j = 0; j < nbr_expected_random_chains.size(); j++) {
		cout << "For chain " << j+1 << ": " << nbr_expected_random_chains.at(j) << endl;
	}
	if(run_type == 2) print_test_result();

/*
	for(int j = 0; j<nbr_pixels; j++) {
		cout << "For pixel " << j << " nbr_implants = " << nbr_implants[j] << endl;
	}
	*/

}

void RandomChains::print_test_result() {
	/*This method calculates the random chains in the test run with a simple formula and prints this result in the terminal window. 
	*/

	cout << "*************************************************" << endl;
	cout << "NOW TEST CALCULATION" << endl;
	cout << "By construction, the rate in every pixel will be the same" << endl;

	double rate_escapes_off = (upper_limit_escapes-lower_limit_escapes)*eoff/experiment_time;
	double rate_alphas_on = (upper_limit_alphas-lower_limit_alphas)*aon/experiment_time;
	
	double rate_escapes_on = (upper_limit_escapes-lower_limit_escapes)*eon/experiment_time;
	double rate_alphas_off = (upper_limit_alphas-lower_limit_alphas)*aoff/experiment_time;

	int nbr_imps = imps;

	double fission_rate = fissions/experiment_time;

	//Here the formula calculation is made
	double test_randoms = nbr_imps*(1-Poisson_pmf(0,rate_alphas_on*time_span.at(0)))*(1-Poisson_pmf(0, rate_escapes_off*time_span.at(1)))*(1-Poisson_pmf(0, rate_alphas_off*time_span.at(2)))*(1-Poisson_pmf(0, rate_escapes_on*time_span.at(3)))*(1-Poisson_pmf(0, fission_rate*time_span.at(4))) * nbr_pixels;
	cout << "test_randoms = nbr_imps*(1-Poisson_pmf(0,rate_alphas_on*time_span.at(0)))*(1-Poisson_pmf(0, rate_escapes_off*time_span.at(1)))*(1-Poisson_pmf(0, rate_alphas_off*time_span.at(2)))*(1-Poisson_pmf(0, rate_escapes_on*time_span.at(3)))*(1-Poisson_pmf(0, fission_rate*time_span.at(4))) * nbr_pixels;" << endl;

	cout << "The CALCULATED total number of random chains with the test data are: " << test_randoms << endl;
	
}

void RandomChains::plot_spectra() {
	/*With this method the user can plot either pixel spectra or total spectrum, assumed that an object of RandomChains has been initialised.
	*/

//This code file contains styles set for the plotting
#include "Style.code"

	TH1F* hist_on;
	TH1F* hist_off;
	int input = -1;
	Bool_t valid_input = kFALSE;

	while(!valid_input) {
		cout << "What would you like to plot? Two options: \n" <<
			"	int < 0: Total energy spectrum.\n" <<
			"	0 <= int < nbr_pixels: Energy spectrum of pixel \"int\"." << endl;
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
	
	TAxis *Xaxis = hist_on->GetXaxis();
	TAxis *Yaxis = hist_on->GetYaxis();

	Yaxis->SetTitle("Counts per 10 keV");
	Xaxis->SetTitle("Energy (MeV)");  
	char title[64];
	if(input >= 0) {
		sprintf(title, "Energy spectrum for pixel %i",  input);
	}
	else {
		sprintf(title, "Energy spectrum for all pixels");
}
hist_on->SetTitle(title);
	hist_on->SetStats(0);
	TGaxis::SetMaxDigits(5); 
	Xaxis->SetTitleOffset(1.2);
	Xaxis->SetTitleSize(0.05);
	//Xaxis->SetRangeUser(5010,11990);      //(5,11995);       //(9510,11490); //4505,11995
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

	hist_on->SetLineWidth(2);
	hist_on->Draw();

	hist_off->SetLineColor(8);
	hist_off->SetLineWidth(3);
	hist_off->SetLineStyle(1);
	hist_off->Draw("same");
	hist_on->Draw("same");

	TLegend* leg = new TLegend(0.8,0.7,0.48,0.9);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg->AddEntry(hist_on,"Beam ON","l");
	leg->AddEntry(hist_off,"Beam OFF","l");
	leg->Draw();

}

//The main function is required for c++ compilation
int main() {

	RandomChains* RC = new RandomChains();
	RC->plot_spectra();
	return 0;
}

//Some ROOT users cannot invoke a function named main from the interpreter. This is the reason for this function.
void run_main() {
	main();
}

/* Mathematical functions (non-member functions) */

//Poisson probability mass function
double Poisson_pmf(int nbr_to_observe, float expected_value) {
	double prob;
	prob = (exp(-expected_value)*pow(expected_value,nbr_to_observe))/(factorial(nbr_to_observe));
	return prob;
}

//Factorial
int factorial(int k){
	int ret;
	ret = 1;
	for (int j = 1; j <= k; j++){
	ret = ret*j;
	}
	return ret;
}


//UF, calculates probability for getting an E115-like chain of random origin
#include <iostream>
#include <string>
#include <sstream>
#include "TH1.h"
#include "TFile.h"
#include <fstream>

using namespace std;



//Main program starts here
long double* u_main(int step) {


int stepper = 0;
  //Prepare/read in root spectra
  TFile  *file = new TFile("Spectra.root");

  TH1F *h_energy_pixel_off_tot = (TH1F*)file->Get("h_energy_pixel_off_tot");
  TH1F *h_energy_pixel_recoff_tot = (TH1F*)file->Get("h_energy_pixel_recoff_tot"); 
  TH1F *h_energy_pixel_on_tot = (TH1F*)file->Get("h_energy_pixel_on_tot"); 
  TH1F *h_energy_pixel_recon_tot = (TH1F*)file->Get("h_energy_pixel_recon_tot"); 

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
  long double EVR[1024] = {0};

  int temp_value = 0;
  long double fission_array[1024] = {0};  
 
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

  int counter = 0;
  if(fission_model == "one_or_average"){
    ifstream input("pixels_with_fissions.txt",ios::in);
    while (input){
      input>>temp_value;
      //if(fission_array[temp_value] == 1) cout<<temp_value<<endl;
      //cout << "temp value UF = " << temp_value << endl;
      fission_array[temp_value] = fission_array[temp_value] + 1;
      counter++;
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

  stepper++;//1
  if(step == stepper) return fission_array;


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

  long double random_prob_pixel[1024] = {0};
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

  //Chain T1
  number_of_alphas = 1;
  lookslike[0] = 'a';
  a1 = "off";
  timespan[0] = 2;
  timespan[2] = 10;

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

    stepper++; //2
    if(step == stepper) return EVR;

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

    stepper++; //3
    if(step == stepper) return rate_a1;

  //Calc rates for FISSION
  for(int i = 0; i < 1024; i++){
    rate_f[i] = fission_array[i]/exp_time;
  }

    stepper++; //4
    if(step == stepper) return rate_f;
  cout << "Number of fissions UF = " << counter << endl;


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

    /*
    if(i == pixel){ 
      cout << "Rate a1 in pixel " << rate_a1[pixel] << endl;
      cout << "Rate f in pixel " << rate_f[pixel] << endl;
    }
    */

    if(number_of_alphas == 2){
      random_prob_pixel[i] = random_prob_pixel[i] * (1-probability(0,rate_a2[i]*timespan[1])); 
      if(i == pixel){ 
	//cout << "Rate a2 in pixel " << rate_a2[pixel] << endl;
      }
    }

    random_prob_sum = random_prob_sum + random_prob_pixel[i];    
  }

 
/*

 cout << EVR[pixel] << "  " << probability(1,rate_a1[pixel]*timespan[0]) << " " << probability(1,rate_f[pixel]*timespan[2]) << endl;

 cout << "Random prob pixel is " << random_prob_pixel[pixel] << endl;   

*/
  cout <<  "The number of expected chains of this type generated by stochastic fluctuations in the background is  " << random_prob_sum << endl;
stepper++; //5
  if(step == stepper) return random_prob_pixel;


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
/*
//Poisson probability mass function
double Poisson_pmf(int nbr_to_observe, double expected_value) {
	double prob;
	prob = (exp(-expected_value)*pow(expected_value,nbr_to_observe))/(factorial(nbr_to_observe));
	return prob;
}

//Factorial
int factorial(int k){
	int ret;
	ret = 1;
	for (int j = 1; j <= k; j++){
	ret = ret*j;
	}
	return ret;
}
*/


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
      //if(pixel_number == 681) cout << "rate is " << rate << " and temporary is " << temporary <<endl;   
 return rate;
} //end of function "rate"


void compare_rates(long double* main, long double* uf, bool print=false) {
	long double compare[1024];
	for(int i = 0; i < 1024; i++) {
		compare[i] = main[i] - uf[i];
		if(compare[i] > 0.000000000000001 || print) {
		cout << "At pixel = " << i << " main = " << main[i] << " and uf = " << uf[i] << " and diff = " << compare[i] << endl;
		}
		//cout << "At pixel = " << i << " main = " << main[i] << " and uf = " << uf[i] << " and diff = " << compare[i] << endl;
	}
}

void run_compare() {
	RandomChains* RC = new RandomChains();
}






