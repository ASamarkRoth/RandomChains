/*! 
	@file RandomChains.cc
	@author Anton Roth (anton.roth@nuclear.lu.se)

	@brief Main program file
*/

/** @mainpage RandomChains
 

@section Description Program Description
	The program RandomChains handles data and computes the number of expected random chains due to random fluctuations in the background for specific decay chains for experimental data. For the description of the method, see U. Forsberg et. al. / Nuclear Physics A 953 (2016) 117-138. The main purposes of the program are:
	- Free access to the code the Lund Nuclear Structure group has used
	- Reproduce the numbers presented in the paper given above
	- Try the method on user provided data

	The complete program is located in the downloaded git repository. 

	The program consists of the following three steps:
	1. Read experimental data. Achieved with the constructor: <i> RandomChains::RandomChains(int pixels, int bins, string folder) </i>
	2. Read decay chain/chains characteristics data. Achieved with the method: <i> RandomChains::SetDecayChains(string input_chains) </i>
	3. Compute the number of expected random chains for the given decay chain/chains with the experimental data. Achieved with the method: <i> RandomChains::Run() </i>


@section Experimental data
	The Lund experimental data (see D. Rudolph et. al / Phys. Rev. Lett. 111, 112502 (2013)) is located in the folder "Lund_data" and should NOT BE MODIFIED. Within this folder the data is stored in ".csv" files. The data consists of spectra for the pixels in the implantation detector for different beam status and reconstructed or not. The structure of these files are: The spectrum histogram bin values are given after eachother separated with a ',' from the first to the last bin and then for the first to the last pixel (the pixel order does not matter). The fissions and in which pixel they have occurred are stored in another file. The data consists of four files:
	- "beam_on.csv": The spectra for all pixels for beam ON. NOT MANDATORY!
	- "recon_beam_on.csv": The reconstructed spectra (for escaped alpha particles) for all pixels for beam ON. MANDATORY!
	- "recon_beam_off.csv": The reconstructed spectra (for escaped alpha particles) for all pixels for beam OFF. MANDATORY!
	- "pixels_with_fissions.csv": The pixel number where a fission has occurred are given in this file. If two fissions have occurred within a pixel this will be presented in this file by two occurrences of this pixel number. MANDATORY!	

	To be able to read in the experimental data, the number of pixels in the implantation detector, the total number of bins in each of the spectra and the folder in which the data is stored, need to be given as input. 
       
	FOR USERS WHO WISH TO RUN THE PROGRAM ON OTHER EXPERIMENTAL DATA: Create a new folder with files with the same names as in the "Lund_data" folder and insert the new data here. In the program, the user only needs to provide the name of this created folder, the number of pixels and the total number of bins for her/his data in the constructor: <i> RandomChains::RandomChains(int pixels, int bins, string folder) </i>.

@section Decay chains
	From a user perspective, what defines a decay chain is the following characteristics:
	- chain_length = x, where x is the length of the chain
	- decay_type = 'a'=alpha, 'e'=escape and 'f'=fission
	- beam_status = 1=ON and 0=OFF
	- time_span = t, where t is the length of the time window during which the decay is accepted.

	In the program the following limits, given in bins <i> E </i>, in the spectra determines the decay type:
	'a'=alpha	<i>lower_limit_alphas <= E < upper_limit_alphas</i>
	'e'=escape	<i>lower_limit_escapes <= E < upper_limit_escapes</i>

	Implants (only for beam ON) are defined as: 
	implants	<i>lower_limit_implants <= E < upper_limit_implants</i>

	The decay chains are given to the program with an input file. In the downloaded repository examples of such input files are "dump_articles.txt" and "dump_test.txt". If no input of decay chains are provided the user can choose between two options:
	0: "Reproduce article numbers", i.e. reproduce the numbers presented in U. Forsberg et. al. / Nuclear Physics A 953 (2016) 117-138.
	2: "Test run". With this option the program is tested on trivial data. The obtained number from the program is compared to a value obtained through the calculation of a simple formula. 

	FOR USERS WHO WISH TO SPECIFY HER/HIS OWN DECAY CHAINS: Have a look at the file "dump_articles.txt". This file contains the input data if one wants to reproduce the article numbers by the Lund group. Edit this file after your own specifications and save it. In the method <i>RandomChains::SetDecayChains(string input_chains)</i> provide the name of your created file as the string <i>input_chains</i>.

	OBS: The duration of the experiment is also given in the same file.

@section Calculate the expected number of random chains
	The expected number of random chains are calculated with the method described in U. Forsberg et. al. / Nuclear Physics A 953 (2016) 117-138.

@section Files and folders
	"RandomChains.cc" - The class RandomChains is implemented here. All methods and functions can be found here.
	"RandomChains.h" - Header file for "RandomChains.cc".
	"run_file.cc" - From this file the user should control and execute the program. Examples of how this can be done already exists here.
	"Makefile" - By typing <i> make </i> within the folder this file is run and the program is built and the executable "run_file" is created.
	"run_file" - The program executable which can be run with <it> ./run_file </i>
	"dump_input.txt" - User input chains are dumped to the following file. It has the same format as the input chains files should be given.
	"dump_article.txt" - Input chains to reproduce the numbers in the Lund article are dumped to this file. 
	"dump_test.txt" - Input chains for the trivial test which can be run to verify the workings of the program are dumped to this file. 
	"Lund_data" - Folder which contains the Lund experimental data.
	"Lund_data/beam_on.csv" - Comma separated file of the spectra for all pixels in the implantation detector for beam ON. 
	"Lund_data/recon_beam_on.csv" - Comma separated file of the reconstructed spectra for all pixels in the implantation detector for beam ON. 
	"Lund_data/recon_beam_off.csv" - Comma separated file of the reconstructed spectra for all pixels in the implantation detector for beam OFF. 
	"Lund_data/pixels_with_fissions.csv" - Comma separated file of the pixels in which a fission occurred. 
*/

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

RandomChains::RandomChains(int pixels, int bins, string folder) : nbr_pixels(pixels), nbr_bins(bins) {
/* The constructor of class RandomChains. 
 	Input arguments:
		int pixels (default=1024): Number of pixels in the spectrum data
		int bins (default=4096): Total number of bins in every spectrum
		string folder (default="Lund_data"): Name of the folder which contains the experimental data
	The following method is invoked:
		void ReadExperimentalData();
	The constructor initialises the following member data:
		string folder_data;
		vector<vector<int>> data_beam_on;
		vector<vector<int>> data_reconstructed_beam_on;
		vector<vector<int>> data_reconstructed_beam_off;
	
Description: In the constructor the experimental data are read in from ".csv" files in the folder given as input argument. The number of pixels in the implantation detector and the total number of bins in each spectrum are also provided as input arguments.

*/

	folder_data = folder + "/";

	ReadExperimentalData();
		
}

void RandomChains::SetDecayChains(string input_chains) {
/* This method sets the decay chain/chains characteristics. 
	Input arguments:
		string input_chains (default=""): The name of the file from which the input should be read
	From this method the following methods are invoked:
		void set_test_chains();
		void generate_test_data();
		void set_article_chains();
		void set_chains_from_input_file(string input_file);
		void dump_input_to_file();
	In this method the following member data is initialised:
		int run_type;
	
Description: This is the main method for setting the decay chains. In this method the decay chain/chains characteristics are set with an input file provided as input argument. If no input file is given, the user gets to choose the type of run; reproduce article numbers or test the program on trivial data. Depending on the type of run different methods are invoked. To verify that this step, the input is dumped to a file named "dump_input.txt" for user provided input and "dump_article.txt" for reproduce article numbers and "dump_test.txt" for the test.
	
*/
	
	bool valid_input = false;

	if(!input_chains.empty()) {
		run_type = 1;
		valid_input = true;
	}

	// The run_type is given by the user.
	while(!valid_input) {
		cout << "What type of run? <0/1/2> \n" << 
			"	0: Reproduce article numbers (no input required) \n "<<
			"	(2: Test the program on trivial data (no input required)) " << endl;
		cin >> run_type;

		switch(run_type) {
			case 0 : 
				cout << "You have chosen: \"Reproduce article numbers\"" << endl;
				valid_input = true;
				break;
			case 2 : 
				cout << "You have chosen: \"test run\"" << endl;
				pure_beam = false;
				valid_input = true;
				break;
			default : 
				cout << "Please insert a number, 0, 1 or 2" << endl;
		}
	}
	

	// The chain/chains characteristics are set.
	if(run_type == 2) {
		//In this method the test chains are set.
		set_test_chains();
		generate_test_data();
		return;
	}
	else if(run_type == 0) {
		//In this method the article chains are set. 
		set_article_chains();
		return;
	}
	else if(run_type == 1) {
		set_chains_from_input_file(input_chains);
	}

	// The input data is dumped to a file.
	dump_input_to_file();
}

void RandomChains::Run() {
/* This method computes the expected number of random chains. 
	Input arguments:

	From this method the following methods are invoked:
		void calculate_implants();
		void rate_calc();
		void calculate_expected_nbr_random_chains();
		void print_result();

	In this method the following member data is initialised:
	
Description: With this method the expected number of random chains due to random fluctuations in the background for the decay chain/chains and experimental data provided. First the number of implants per pixel is calculated. Then the background rates in every pixel for the different decay types are calculated. This is followed by the calculation of the expected number of random chains. Finally, the results are printed in the terminal window. 
*/

	// The number of implants per pixel is calculated.
	calculate_implants();

	// The background rates for alphas, escapes for beam ON and OFF and fission are calculated for every decay and pixel for that decay.
	cout << "Calculating rates " << endl;
	for(unsigned int i = 0; i < decay_type.size(); i++) {
		rate_calc(decay_type.at(i), beam_status.at(i));
	}

	// The TOTAL number of expected random chains due to random fluctuations in the background are calculated for the specific chain/chains given as input to the program.
	calculate_expected_nbr_random_chains();

	// The result is printed in the terminal window. 
	print_result();

}

void RandomChains::ReadExperimentalData() {
/* The experimental data is read in. 
	Input arguments:

	From this method the following methods are invoked:
		void read_exp_file(string file_name);

	In this method the following member data is initialised:
		vector<vector<int>> data_beam_on;
		vector<vector<int>> data_reconstructed_beam_on;
		vector<vector<int>> data_reconstructed_beam_off;
		vector<double> fissions_pixels;
		vector<int> nbr_implants;

Description: The experimental data is read in from the folder provided in the constructor. All vectors are initialised. The spectrum data and fission data are read in with the method <i> read_exp_file(string file_name) </i>.
*/
	cout << "Reading experimental data from the relative path: " << folder_data << endl;

	data_beam_on.resize(nbr_pixels);
	data_reconstructed_beam_on.resize(nbr_pixels);
	data_reconstructed_beam_off.resize(nbr_pixels);
	fissions_pixels.resize(nbr_pixels);
	nbr_implants.resize(nbr_pixels);
	for(int j = 0; j < nbr_pixels; j++) {
		data_beam_on[j].resize(nbr_bins);
		data_reconstructed_beam_on[j].resize(nbr_bins);
		data_reconstructed_beam_off[j].resize(nbr_bins);
	}

	string read_file;

	read_file = "beam_on.csv";
	read_exp_file(read_file);

	read_file = "rec_beam_on.csv";
	read_exp_file(read_file);

	read_file = "rec_beam_off.csv";
	read_exp_file(read_file);
	
	read_file = "pixels_with_fissions.csv";
	read_exp_file(read_file);
}

void RandomChains::read_exp_file(string read_file) {
/* The experimental data files are read in. 
	Input arguments:
		string read_file: The name of the file to be read in.

	From this method the following methods are invoked:

	In this method the following member data are set:
		vector<vector<int>> data_beam_on;
		vector<vector<int>> data_reconstructed_beam_on;
		vector<vector<int>> data_reconstructed_beam_off;
		vector<double> fissions_pixels;

Description: The experimental data in the comma separated files are read in from the folder provided in the constructor. The files read in are: "beam_on.csv", "recon_beam_on.csv", "recon_beam_off.csv" and "pixels_with_fissions.csv"
*/
	
	int bin = 0; int pixel = 0;
	string val;

	ifstream ifile_stream(folder_data + read_file, ios::in);

	if(!ifile_stream) {
		cout << "File \"" << folder_data+read_file << "\" was not found " << endl;
		if(read_file != "beam_on.csv") {
			cout << "File \"" << read_file << "\" is essential for the analysis. Please add this file! " << endl;
			abort();
		}
		else if (read_file == "beam_on.csv") {
			cout << "OBS: The reconstructed data will be used instead of pure beam ON data!" << endl;
			pure_beam = false;
			return;
		}
	}

	cout << "Reading file " << folder_data+read_file << endl;

	//The fission data are read in here and treated differently.
	if(read_file == "pixels_with_fissions.csv") {
		for(int i = 0; i < nbr_pixels; i++) {
			fissions_pixels[i] = 0; 
		}
		int nbr_of_fissions = 0;
		while(getline(ifile_stream, val, ',')) {
			//cout << "number = " << nbr_of_fissions << " fission val = " << val << " Val_empty = " << val.empty() << " ";
			if(val.empty()) {
				cout << "Breaking..." << endl;
			       	break;
			}
			fissions_pixels[stoi(val)] += 1;
			nbr_of_fissions++;
		}
		cout << "Total number of fissions are: " << nbr_of_fissions << endl;

		//If the number of fissions in a pixel is 0 then it is set to the average over the complete implantation detector. 
		for(int i = 0; i < nbr_pixels; i++){
			if(fissions_pixels[i] == 0){
				fissions_pixels[i] = (double)nbr_of_fissions/nbr_pixels;    
			} 
		}

		return;
	}


	while(getline(ifile_stream, val, ',')) {
		//cout << "bin = " << bin << " val = " << val << endl;
		if(bin%nbr_bins==0 && bin > 0) {
			bin = 0;
			pixel++;
			//cout << "pixel = " << pixel << endl;
		}
		if(read_file == "beam_on.csv") data_beam_on[pixel][bin] = stoi(val);
		else if(read_file == "rec_beam_on.csv") data_reconstructed_beam_on[pixel][bin] = stoi(val);
		else if(read_file == "rec_beam_off.csv") data_reconstructed_beam_off[pixel][bin] = stoi(val);

		bin++;
	}

	if(bin%nbr_bins == 0 && (pixel+1)%nbr_pixels == 0) {
		cout << "The file was successfully read " << endl;
	}
	else {
		cout << "Something wrong with the read in ... . The following might hint on what is wrong: " << endl;
			cout << "Number of pixels read in was " << pixel+1 << endl;
			cout << " and number of bins for the last pixel was " << bin << endl;
	}

}

RandomChains::~RandomChains() {
	/*The destructor of RandomChains. Object is deleted.*/
	delete this;
}

void RandomChains::generate_test_data() {
/* The test data is generated. 
	Input arguments:

	From this method the following methods are invoked:

	In this method the following member data are set:
		int eon;
		int eoff;
		int non;
		int noff;
		int aon;
		int aoff;
		int imps;
		int fissions;
		vector<vector<int>> data_beam_on;
		vector<vector<int>> data_reconstructed_beam_on;
		vector<vector<int>> data_reconstructed_beam_off;
		vector<double> fissions_pixels;

Description: The test data is generated in this method. All the read in data is substituted.
*/
	//Clearing the data
	for(int k = 0; k < nbr_pixels; k++) {
		for(int i = 0; i < nbr_bins; i++) {
			data_reconstructed_beam_on[k][i] = 0;
			data_reconstructed_beam_off[k][i] = 0;
		}
	}

	//Setting the values to insert in the test spectra:
	eon = 4;
	eoff = 3;
	non = 3;
	noff = 2;
	aon = 2;
	aoff = 1;
	imps = 100;
	fissions = 2;

	for(int k = 0; k < nbr_pixels; k++) {
		for(int i = lower_limit_escapes; i < upper_limit_alphas; i++) {
			if(i < upper_limit_escapes) {
				data_reconstructed_beam_on[k][i] = eon;
				data_reconstructed_beam_off[k][i] = eoff;
			}
			else if(i >= upper_limit_escapes && i < lower_limit_alphas) {
				data_reconstructed_beam_on[k][i] = non;
				data_reconstructed_beam_off[k][i] = noff;
			}
			else {
				data_reconstructed_beam_on[k][i] = aon;
				data_reconstructed_beam_off[k][i] = aoff;
			}
		}

		int middle = lower_limit_implants + floor((upper_limit_implants - lower_limit_implants)/2);
		data_reconstructed_beam_on[k][middle] = imps;
		
		fissions_pixels[k] = fissions;
	}

}


void RandomChains::set_test_chains() {
/* The test chains are set. 
	Input arguments:

	From this method the following methods are invoked:

	In this method the following member data are set:
		vector<int> chain_length;
		vector<int> beam_status;
		vector<char> decay_type;
		vector<double> time_span;
		int lower_limit_alphas, upper_limit_alphas;
		int lower_limit_escapes, upper_limit_escapes;
		int lower_limit_implants, upper_limit_implants;
		double experiment_time;

Description: This method defines the test chains characteristics and limits for the different signal types. The experiment duration is also set. 
*/
	
	//Do not modify these numbers!!!
	chain_length = {5};
	decay_type = {'a', 'e', 'a', 'e', 'f'};
	beam_status = {1, 0, 0, 1, 1};
	time_span = {1, 2, 3, 4, 5};

	lower_limit_alphas = 900; upper_limit_alphas = 1100;
	lower_limit_escapes = 0; upper_limit_escapes = 400;
	lower_limit_implants = 1100; upper_limit_implants = 1800;

	experiment_time = 1000000;
}

void RandomChains::set_article_chains() {
/* The Lund article chains are set. 
	Input arguments:

	From this method the following methods are invoked:

	In this method the following member data are set:
		vector<int> chain_length;
		vector<int> beam_status;
		vector<char> decay_type;
		vector<double> time_span;
		int lower_limit_alphas, upper_limit_alphas;
		int lower_limit_escapes, upper_limit_escapes;
		int lower_limit_implants, upper_limit_implants;
		double experiment_time;

Description: This method defines the Lund article chains characteristics and limits for the different signal types. The experiment duration is also set. 
*/
	
	//Do not modify these numbers!!!
	chain_length = {2, 2, 3, 3, 3, 3, 3};
	decay_type = {'a','f',  'e','f',  'a','a','f',  'a','a','f',  'a','a','f',  'a','a','f',  'e','e','f'};
	beam_status = {0,0,  0,0,  1,0,0,  1,0,0,  0,0,0,  0,0,0,  0,1,0};
	time_span = {2,10,  2,10,  2,10,50,  2,10,50,  2,10,50,  2,10,50, 2,10,50};

	lower_limit_alphas = 900; upper_limit_alphas = 1100;
	lower_limit_escapes = 0; upper_limit_escapes = 400;
	lower_limit_implants = 1100; upper_limit_implants = 1800;

	experiment_time = 1433000;
}

void RandomChains::dump_input_to_file() {
/* This method dumps the input chains to a file. 
	Input arguments:

	From this method the following methods are invoked:

	In this method the following member data are set:

Description: This method dumps the set chains to a file. If the test is run the file name is "dump_test.txt" and if the reproduce article numbers is run the file name is "dump_article.txt". Otherwise the file name is "dump_input.txt". 
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
	dump << "Lines starting with a '#' indicates the start of a new chain. The 2nd and 4th lines are read in, here the experimental time and the bin limits for the different signal types are given. The format is very important! " << endl;
	dump << "Experiment_time(s): " << experiment_time << endl;
	dump << "alpha_low alpha_up escape_low escapes_up implants_low implants_up" << endl;
	dump << lower_limit_alphas <<" "<< upper_limit_alphas <<" "<< lower_limit_escapes <<" "<<upper_limit_escapes<<" "<<lower_limit_implants<<" "<<upper_limit_implants << endl;
	dump << "Type (alpha=a, escape=e and fission=f) 	Beam ON (=1) or OFF (=0)	Time span (s) \n";
	if(run_type == 0) {
		cout << "Lines starting with a '#' indicates the start of a new chain. The 2nd and 4th lines are read in, here the experimental time and the bin limits for the different signal types are given. The format is very important! " << endl;
		cout << "Experiment_time(s): " << experiment_time << endl;
		cout << "alpha_low alpha_up escape_low escapes_up implants_low implants_up" << endl;
		cout << lower_limit_alphas <<" "<< upper_limit_alphas <<" "<< lower_limit_escapes <<" "<<upper_limit_escapes<<" "<<lower_limit_implants<<" "<<upper_limit_implants << endl;
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

void RandomChains::set_chains_from_input_file(string input_chains) {
/* This method sets the chain/chains from a given input file. 
	Input arguments:
		string input_chains: File name of the input chains which are provided by the user in the method <i> SetDecayChains </i>

	From this method the following methods are invoked:

	In this method the following member data are set:
		vector<int> chain_length;
		vector<int> beam_status;
		vector<char> decay_type;
		vector<double> time_span;
		int lower_limit_alphas, upper_limit_alphas;
		int lower_limit_escapes, upper_limit_escapes;
		int lower_limit_implants, upper_limit_implants;
		double experiment_time;

Description: This method sets the chains from a file. It is necessary that the format of the file which is to be read in follows the correct format. The spaces are important! What is read in is printed in the terminal window and can also be found in the file dumped after a run.
*/

	string filename;

	if(input_chains.empty()) {
		cout << "Enter full name of input file name: (file \"dump_input.txt\" should be available): ";
		cin >> filename;
	}
	else filename = input_chains;
	ifstream file_stream(filename,ios::in);
	if(!file_stream) cout << "Could not find file" << endl;
	int beam;
	double time;
	char type;
	string str;
	stringstream ss;
	int counter = 0;
	cout << "The following was read in: " << endl;
	while(getline(file_stream, str)) {
		cout << str << endl;
		//if(chain_length.size() > 0) cout << "chain length = " << chain_length.at(0) << endl;
		if(counter < 5) {
			if(counter == 1) {
				if(str.find_first_of(" ") != string::npos) {
					experiment_time = stod(string(str.begin()+str.find_first_of(" "), str.end()));
				}
			}
			if(counter == 3) {
				ss.str(string());
				ss.clear();
				ss.str(str);
				ss >> lower_limit_alphas >> upper_limit_alphas >> lower_limit_escapes >> upper_limit_escapes >>lower_limit_implants >> upper_limit_implants;
			}			
			counter++;
			continue;
		}

		if(str[0] == '#') {
			string number = string(str.begin()+1, str.end());
			chain_length.push_back(stoi(number));
		}
		else {
			ss.str(string());
			ss.clear();
			ss.str(str);
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
/* Calculates the number of implants. 
	Input arguments:

	From this method the following methods are invoked:

	In this method the following member data are set:
		vector<int> nbr_implants;

Description: The number of implants in every pixel is calculated with the lower and upper limits set in <i>SetDecayChains</i>. 
*/


	vector< vector<int> > data;

	if(pure_beam) data = data_beam_on;
	else data = data_reconstructed_beam_on;

	for(int i = 0; i < nbr_pixels; i++) {
		int acc_counts = 0; 
		for(int k = lower_limit_implants; k < upper_limit_implants; k++) {
			acc_counts += data[i][k];
		}
		nbr_implants[i] = acc_counts;
	}

}

void RandomChains::calculate_rates() {
/* This method calculates the rates in every pixel for the specific decay types, one decay at a time.

 
	Input arguments:

	From this method the following methods are invoked:

	In this method the following member data are set:
		vector<int> nbr_implants;

Description: The number of implants in every pixel is calculated with the lower and upper limits set in <i>SetDecayChains</i>. 
*/
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
/* This method calculates the rates in every pixel for the specific decay types and beam status for a decay.
 
	Input arguments:
		char type=decay type, i.e. 'a', 'e' or 'f'.
		int beam=beam status, i.e. 1 or 0.

	From this method the following methods are invoked:

	In this method the following member data are set:
		vector< vector<double> > rate;

Description: Given the decay type and beam status the rate for every pixel is calculated and stored in the 2D vector <i>rate</i>. 
*/

	//Based on the beam status the spectrum is determined
	vector< vector<int> > data;
	if(beam) data = data_reconstructed_beam_on;
	else data = data_reconstructed_beam_off;


	vector<double> rate_temp;
	rate_temp.resize(nbr_pixels);
	
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

	//The rate for every pixel is calculated
	for(int i = 0; i < nbr_pixels; i++) {
		int acc_counts = 0;
		for(int k = lower_limit; k < upper_limit; k++) {
			acc_counts += data[i][k];
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
/* This method calculates the TOTAL number of expected random chains for the input decay chain/chains.
 
	Input arguments:

	From this method the following methods are invoked:

	In this method the following member data are set:
		vector<double> nbr_expected_random_chains;

Description: On the basis of the rates calculated for every decay the expected number of random chains due to random fluctuations in the background are determined per pixel and decay chain. The values of every pixel are then summed for every decay chain to a final value. 
*/

	cout << "Calculating expected number of random chains " << endl;
	int offset = 0;

	//looping decay chains
	for(unsigned int j = 0; j < chain_length.size(); j++) {
		vector<double> randoms_in_pixel(nbr_pixels,0.);

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
/* The results of the run are printed.
 
	Input arguments:

	From this method the following methods are invoked:

	In this method the following member data are set:
		vector<double> nbr_expected_random_chains;

Description: This is the method that is invoked at the end of the constructor and it presents the result of the run in the terminal window. If the test was run another member function is called for further output. 
*/

	cout << "**************************************************" << endl;
	if(run_type == 2) cout << "These are the result of the TEST run: " << endl;
	else cout << "These are the result of the run: " << endl;

	cout << "The total number of expected random chains of the same type as the given chain due to random fluctuations in the background are: " << endl;

	for(unsigned int j = 0; j < nbr_expected_random_chains.size(); j++) {
		cout << "For chain " << j+1 << ": " << nbr_expected_random_chains.at(j) << endl;
	}
	if(run_type == 2) print_test_result();

}

void RandomChains::print_test_result() {
/* The results of the test run are printed.
 
	Input arguments:

	From this method the following methods are invoked:

	In this method the following member data are set:

Description: This method calculates the random chains in the test run with a simple formula and prints this result in the terminal window.
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

/* Mathematical functions (non-member functions) */

double Poisson_pmf(int nbr_to_observe, float expected_value) {
/* Poisson probability mass function. <i>p(k) = lambda^k exp(-lambda)/k!</i>
 
	Input arguments:
		double expected_value: lambda
		int nbr_to_observe: k
	
	Returns:
		<i>p(k)</i>, i.e. the probability to observe k observations

*/
	double prob;
	prob = (exp(-expected_value)*pow(expected_value,nbr_to_observe))/(factorial(nbr_to_observe));
	return prob;
}

int factorial(int k){
/* Factorial 
	Input arguments:
		int k
	
	Returns:
		int factorial
*/
	int ret;
	ret = 1;
	for (int j = 1; j <= k; j++){
	ret = ret*j;
	}
	return ret;
}

