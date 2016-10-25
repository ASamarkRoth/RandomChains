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
			/*
			cout << "Temp_value = " << temp_value << endl;
			cout << "Temp_value id = " << typeid(temp_value).name() << endl;
			*/
			fissions_pixels[temp_value] += 1;
			nbr_of_fissions++;
			input>>temp_value;
		}
		//cout << "Total number of fissions are: " << nbr_of_fissions << endl;

		//If the number of fissions in a pixel is 0 then it is set to the average over the complete implantation detector. 
		for(int i = 0; i < nbr_pixels; i++){
			if(fissions_pixels[i] == 0){
				fissions_pixels[i] = (double)nbr_of_fissions/nbr_pixels;    
			} 
		}

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

	if(run_type == 1 || run_type == 0) hist = h_energy_pixel_beam_on;
	else hist = h_energy_pixel_reconstructed_beam_on;

	for(int i = 0; i < nbr_pixels; i++) {
		int acc_counts = 0; 
		for(int k = lower_limit_implants; k < upper_limit_implants; k++) {
			acc_counts += hist[i]->GetBinContent(k);
		}
		nbr_implants[i] = acc_counts;
	}


}

void RandomChains::calculate_rates() {
	/*This method calculates the rates in every pixel for the specific decay types, one decay at a time.

	Through this method the following member data is initialised:
		vector<array<double, nbr_pixels>> rate;
	
	*/

	cout << "Calculating rates " << endl;

	for(unsigned int i = 0; i < decay_type.size(); i++) {
		rate_calc(decay_type.at(i), beam_status.at(i));
	}

}

void RandomChains::rate_calc(char type, int beam) {
	/*This method calculates the rate per pixel depending on the type of decay and the beam status.
	Input arguments are:
	"char type"=decay type, i.e. 'a', 'e' or 'f'.
	"int beam"=beam status, i.e. 1/0

	Through this method the following member data is modified:
		vector<array<double, nbr_pixels>> rate;
		
	*/


	//Based on the beam status the spectrum is determined
	TH1F** hist;
	if(beam) hist = h_energy_pixel_reconstructed_beam_on;
	else hist = h_energy_pixel_reconstructed_beam_off;

	array<double, nbr_pixels> rate_temp;
	
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
			rate_temp[i] = (double)fissions_pixels[i]/experiment_time;
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
		rate_temp[i] = (double) acc_counts/experiment_time;
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
		double randoms_in_pixel[nbr_pixels] = {0};

		//looping decays
		for(int l = offset; l < offset+chain_length.at(j); l++) {
			/*
			cout << "l = " << l << endl;
			cout << "decay type = " << decay_type.at(l) << endl;
			cout << "beam status = " << beam_status.at(l) << endl;
			cout << "time span = " << time_span.at(l) << endl;
			*/

			//looping pixels
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

//Poisson probability mass function (For UF:s float expected value)
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

