/*!
@file RandomChains.cc
@author Anton Roth (anton.roth@nuclear.lu.se)

@brief Main program file
*/

/** @mainpage RandomChains


@section Description_tag Program Description
        The program RandomChains handles data and computes the number
	of expected random chains due to random fluctuations in the
	background for specific decay chains for experimental
	data. For the description of the method, see <a
	href="http://www.sciencedirect.com/science/article/pii/S0375947416300768">U. Forsberg
	et. al. Nuclear Physics A 953 (2016) 117-138</a>. The main
	purposes of the program are:

        - Provide free access to the code the Lund Nuclear Structure group has
	used

        - Reproduce the numbers presented in the paper given
	above

        - Try the method on user provided chains and experimental data

	The complete program is located in the downloaded git repository.
	The program consists of the following three steps:

	-# Read experimental data. Achieved with the constructor:
           RandomChains::RandomChains(int pixels, int bins, string folder)

	-# Read decay chain/chains characteristics data. Achieved with the method:
           RandomChains::SetDecayChains(string input_chains)

	-# Compute the number of expected random chains for the given
           decay chain/chains with the experimental data. Achieved
           with the method: RandomChains::Run()

@subsection run_tag How to run RandomChains?
	The program is preferably controlled from the file <tt>run_file.cc</tt>.

	When located in the directory of RandomChains type the following in the terminal to run the program:

	-# <tt> make </tt>

	-# <tt> ./run_file </tt>

@subsection computer_tag Computer requirements
	The program has been successfully run on the following systems:
      -	Linux Ubuntu 14.04 with g++ version 4.8.4 and std=c++11



@section Experiment_tag Experimental data
        The Lund experimental data (see <a
	href="http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.111.112502">D. Rudolph
	et. al / Phys. Rev. Lett. 111, 112502 (2013)</a>) is located
	in the folder <tt>Lund_data</tt> and <b>SHOULD NOT BE
	MODIFIED</b>. Within this folder the data is stored in
	'<tt>.csv</tt>' files. The data consists of spectra for the
	pixels in the implantation detector for different beam status
	and reconstructed or not. The structure of these files are:

        The spectrum histogram bin values are given after each other
	separated with a ',' from the first to the last bin and then
	for the first to the last pixel (the pixel order does not
	matter). The fissions and in which pixel they have occurred
	are stored in another file. The data consists of four files:

	- <tt>beam_on.csv</tt>: The spectra for all pixels for beam
          ON. <b>OPTIONAL!</b>

	- <tt>rec_beam_on.csv</tt>: The reconstructed spectra (for
          escaped alpha particles) for all pixels for beam
          ON. <b>MANDATORY!</b>

	- <tt>rec_beam_off.csv</tt>: The reconstructed spectra (for
          escaped alpha particles) for all pixels for beam
          OFF. <b>MANDATORY!</b>

	- <tt>pixels_with_fissions.csv</tt>: The pixel number where a
          fission has occurred are given in this file. If two fissions
          have occurred within a pixel this will be presented in this
          file by two occurrences of this pixel
          number. <b>MANDATORY!</b>

	To be able to read in the experimental data, the number of
	pixels in the implantation detector, the total number of bins
	in each of the spectra and the folder in which the data is
	stored, need to be given as input.

@subsection users_tag FOR USERS WHO WISH TO RUN THE PROGRAM ON OTHER EXPERIMENTAL DATA:
        Create a new folder with files with the same names as in the
	<tt>Lund_data</tt> folder and insert the new data here. In the
	program, the user only needs to provide the name of this
	created folder, the number of pixels and the total number of
	bins for the data in the constructor: RandomChains::RandomChains(int pixels, int bins, string folder)

@section decay_chains_tag Decay chains
        From a user perspective, what defines a decay chain is the
	following characteristics:

        - <tt>chain_length</tt>: <em>x</em>, where <em>x</em> is the length of the chain

        - <tt>decay_type</tt>: <em>a</em>=alpha, <em>e</em>=escape and <em>f</em>=fission

        - <tt>beam_status</tt>: <em>1</em>=ON and <em>0</em>=OFF

        - <tt>time_span</tt>: <em>t</em>, where <em>t</em> is the
	length of the time window during which the decay is accepted.

	In the program the following limits, given in bins <em>E</em>, in the spectra determines the decay type:
        <table>
        <tr>
        <td><em>a</em>=alpha</td>
        <td><tt>lower_limit_alphas</tt> <= <tt>E</tt> < <tt>upper_limit_alphas</tt></td>
        </tr>
        <tr>
        <td><em>e</em>=escape</td>
	<td><tt>lower_limit_escapes</tt> <= <tt>E</tt> < <tt>upper_limit_escapes</tt></td>
        </tr>
        <tr>
        <td colspan="2">Implants (only for beam ON) are defined as:</td>
        </tr>
        <tr>
	<td>implants</td>
	<td><tt>lower_limit_implants</tt> <= <tt>E</tt> < <tt>upper_limit_implants</tt></td>
        </tr>
        </table>

	The decay chains are given to the program with an input
	file. In the downloaded repository examples of such input
	files are <tt>dump_article.txt</tt> and
	<tt>dump_test.txt</tt>. If no input of decay chains are
	provided the user can choose between two options:

        - 0: '<em>Reproduce article numbers</em>', i.e. reproduce the
	numbers presented in <a
	href="http://www.sciencedirect.com/science/article/pii/S0375947416300768">U. Forsberg
	et. al. Nuclear Physics A 953 (2016) 117-138</a>.

        - 2: '<em>Test run</em>'. With this option the program is
	tested on trivial data. The obtained number from the program
	is compared to a value obtained through the calculation of a
	simple formula.

@subsection users_chain_tag FOR USERS WHO WISH TO SPECIFY HER/HIS OWN DECAY CHAINS:
        Have a look at the file <tt>dump_articles.txt</tt>. This file
        contains the input data if one wants to reproduce the article
        numbers by the Lund group. Edit this file after your own
        specifications and save it. In the method
        RandomChains::SetDecayChains(string input_chains)
        provide the name of your created file as the string
        <tt>input_chains</tt>.

	OBS: The duration of the experiment is also given in the same file.

@section random_tag Calculate the expected number of random chains
	The expected number of random chains is calculated with the
	method described in <a
	href="http://www.sciencedirect.com/science/article/pii/S0375947416300768">U. Forsberg
	et. al. Nuclear Physics A 953 (2016) 117-138</a>.

@section files_tag Files and folders
        A list of files and folders is provided below:

	RandomChains.cc: The class RandomChains is implemented here. All methods and functions can be found here.

        RandomChains.h: Header file for RandomChains.cc.

	run_file.cc: From this file the user should control and
	execute the program. Examples of how this can be done already
	exists here.

	Makefile: By typing <tt>make</tt> within the folder this file
	is run and the program is built and the executable
	<tt>run_file</tt> is created.

	run_file: The program executable which can be run with
	<tt>./run_file</tt>

	dump_input.txt: User input chains are dumped to the following
	file. It has the same format as the input chains files should
	be given.

	dump_article.txt: Input chains to reproduce the numbers in the
	Lund article are dumped to this file.

	dump_test.txt: Input chains for the trivial test which can be
	run to verify the workings of the program are dumped to this
	file.

	Lund_data: Folder which contains the Lund experimental data.

	Lund_data/beam_on.csv: Comma separated file of the spectra for
	all pixels in the implantation detector for beam ON.

	Lund_data/rec_beam_on.csv: Comma separated file of the
	reconstructed spectra for all pixels in the implantation
	detector for beam ON.

	Lund_data/rec_beam_off.csv: Comma separated file of the
	reconstructed spectra for all pixels in the implantation
	detector for beam OFF.

	Lund_data/pixels_with_fissions.csv: Comma separated file of
	the pixels in which a fission occurred.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include "RandomChains.h"
#include <assert.h>
#include "math.h"
#include <typeinfo>

using namespace std;

// The constructor of class RandomChains.
/**
<p> In the constructor the experimental data are read in from ".csv" files in the folder given as input argument. The number of pixels in the implantation detector and the total number of bins in each spectrum are also provided as input arguments.</p>
	@param pixels number of pixels in the spectrum data
	@param bins total number of bins in every spectrum
	@param folder name of the folder which contains the experimental data
	@returns returns object of the class RandomChains

	@see ReadExperimentalData()

The following is initialised:
	- RandomChains::folder_data
	- RandomChains::data_beam_on
	- RandomChains::data_reconstructed_beam_on
	- RandomChains::data_reconstructed_beam_off


*/
RandomChains::RandomChains(int pixels, int bins, string folder) : nbr_pixels(pixels), nbr_bins(bins) {

	folder_data = folder + "/";

	ReadExperimentalData();

}

/** This method sets the decay chain/chains characteristics.
This is the main method for setting the decay chains. In this method the decay chain/chains characteristics are set with an input file provided as input argument. If no input file is given, the user gets to choose the type of run; reproduce article numbers or test the program on trivial data. Depending on the type of run different methods are invoked. To verify that this step was made successfully, the input is dumped to a file named "dump_input.txt" for user provided input and "dump_article.txt" for reproduce article numbers and "dump_test.txt" for the test.
	@param input_chains the name of the file from which the input should be read

	@see set_test_chains()
	@see generate_test_data()
	@see set_article_chains()
	@see set_chains_from_input_file(string input_file)
	@see dump_input_to_file()

	The following is initialised:
	- RandomChains::run_type

*/

void RandomChains::SetDecayChains(string input_chains) {
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

/** This method computes the expected number of random chains.
With this method the expected number of random chains due to random fluctuations in the background for the decay chain/chains and experimental data provided. First the number of implants per pixel is calculated. Then the background rates in every pixel for the different decay types are calculated. This is followed by the calculation of the expected number of random chains. Finally, the results are printed in the terminal window.
		@see RandomChains::calculate_implants()
		@see RandomChains::rate_calc()
		@see RandomChains::calculate_expected_nbr_random_chains()
		@see RandomChains::print_result()

*/
void RandomChains::Run() {

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

/** The experimental data is read in.
The experimental data is read in from the folder provided in the constructor. All vectors are initialised. The spectrum data and fission data are read in with the method <em> read_exp_file(string file_name) </em>.
		@see RandomChains::read_exp_file(string file_name)

	The following data is initialised:
		- RandomChains::data_beam_on
		- RandomChains::data_reconstructed_beam_on
		- RandomChains::data_reconstructed_beam_off
		- RandomChains::fissions_pixels
		- RandomChains::nbr_implants

*/
void RandomChains::ReadExperimentalData() {
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

/** The experimental data files are read in.
The experimental data in the comma separated files are read in from the folder provided in the constructor. The files read in are: "beam_on.csv", "recon_beam_on.csv", "recon_beam_off.csv" and "pixels_with_fissions.csv"
		@param read_file the name of the file to be read in.

	The following is initialised:
		- RandomChains::data_beam_on
		- RandomChains::data_reconstructed_beam_on
		- RandomChains::data_reconstructed_beam_off
		- RandomChains::fissions_pixels

*/
void RandomChains::read_exp_file(string read_file) {

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

/**The destructor of RandomChains. Object is deleted.*/
RandomChains::~RandomChains() {
	delete this;
}

/** The test data is generated.
The test data is generated in this method. All the read in data is substituted.

The following is initialised:
	- RandomChains::eon
	- RandomChains::eoff
	- RandomChains::non
	- RandomChains::noff
	- RandomChains::aon
	- RandomChains::aoff
	- RandomChains::imps
	- RandomChains::fissions
	- RandomChains::data_beam_on
	- RandomChains::data_reconstructed_beam_on
	- RandomChains::data_reconstructed_beam_off
	- RandomChains::fissions_pixels

*/
void RandomChains::generate_test_data() {
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


/** The test chains are set.
This method defines the test chains characteristics and limits for the different signal types. The experiment duration is also set.
	The following is initialised:

	- RandomChains::chain_length
	- RandomChains::beam_status
	- RandomChains::decay_type
	- RandomChains::time_span
	- RandomChains::lower_limit_alphas, upper_limit_alphas
	- RandomChains::lower_limit_escapes, upper_limit_escapes
	- RandomChains::lower_limit_implants, upper_limit_implants
	- RandomChains::experiment_time

*/
void RandomChains::set_test_chains() {

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

/** The Lund article chains are set.
This method defines the Lund article chains characteristics and limits for the different signal types. The experiment duration is also set.

The following is initialised:
	- RandomChains::chain_length
	- RandomChains::beam_status
	- RandomChains::decay_type
	- RandomChains::time_span
	- RandomChains::lower_limit_alphas, upper_limit_alphas
	- RandomChains::lower_limit_escapes, upper_limit_escapes
	- RandomChains::lower_limit_implants, upper_limit_implants
	- RandomChains::experiment_time

*/
void RandomChains::set_article_chains() {

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

/** This method dumps the input chains to a file.
This method dumps the set chains to a file. If the test is run the file name is "dump_test.txt" and if the reproduce article numbers is run the file name is "dump_article.txt". Otherwise the file name is "dump_input.txt".
*/
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

/** This method sets the chain/chains from a given input file.
This method sets the chains from a file. It is necessary that the format of the file which is to be read in follows the correct format. The spaces are important! What is read in is printed in the terminal window and can also be found in the file dumped after a run.
	@param input_chains file name of the input chains which are provided by the user in the method <em> SetDecayChains </em>

The following is initialised:
	- RandomChains::chain_length
	- RandomChains::beam_status
	- RandomChains::decay_type
	- RandomChains::time_span
	- RandomChains::lower_limit_alphas, upper_limit_alphas
	- RandomChains::lower_limit_escapes, upper_limit_escapes
	- RandomChains::lower_limit_implants, upper_limit_implants
	- RandomChains::experiment_time

*/
void RandomChains::set_chains_from_input_file(string input_chains) {

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


/** Calculates the number of implants.
The number of implants in every pixel is calculated with the lower and upper limits set in <em>SetDecayChains</em>.

The following is initialised:
		- RandomChains::nbr_implants
*/
void RandomChains::calculate_implants() {

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

/** This method calculates the rates in every pixel for the specific decay types, one decay at a time.
The rate in every pixel is calculated with the lower and upper limits set in <em>SetDecayChains</em>.

The following is initialised:
	- RandomChains::nbr_implants

*/
void RandomChains::calculate_rates() {
	cout << "Calculating rates " << endl;

	for(unsigned int i = 0; i < decay_type.size(); i++) {
		rate_calc(decay_type.at(i), beam_status.at(i));
	}

}

/** This method calculates the rates in every pixel for the specific decay types and beam status for a decay.
Given the decay type and beam status the rate for every pixel is calculated and stored in the 2D vector <em>rate</em>.
		@param type decay type, i.e. 'a', 'e' or 'f'.
		@param beam beam status, i.e. 1 or 0.

The following is intialised:
		- RandomChains::rate

*/
void RandomChains::rate_calc(char type, int beam) {

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

/** This method calculates the TOTAL number of expected random chains for the input decay chain/chains.
On the basis of the rates calculated for every decay the expected number of random chains due to random fluctuations in the background are determined per pixel and decay chain. The values of every pixel are then summed for every decay chain to a final value.

The following is initialised:
		- RandomChains::nbr_expected_random_chains
*/
void RandomChains::calculate_expected_nbr_random_chains() {

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

/** The results of the run are printed.
This is the method that is invoked at the end of the constructor and it presents the result of the run in the terminal window. If the test was run another member function is called for further output.

The following is initialised:
		- RandomChains::nbr_expected_random_chains

*/
void RandomChains::print_result() {

	cout << "**************************************************" << endl;
	if(run_type == 2) cout << "These are the result of the TEST run: " << endl;
	else cout << "These are the result of the run: " << endl;

	cout << "The total number of expected random chains of the same type as the given chain due to random fluctuations in the background are: " << endl;

	for(unsigned int j = 0; j < nbr_expected_random_chains.size(); j++) {
		cout << "For chain " << j+1 << ": " << nbr_expected_random_chains.at(j) << endl;
	}
	if(run_type == 2) print_test_result();

}

/** The results of the test run are printed.
This method calculates the random chains in the test run with a simple formula and prints this result in the terminal window.
*/
void RandomChains::print_test_result() {

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

/** Poisson probability mass function. <em>p(k) = lambda^k exp(-lambda)/k!</em>
	@param expected_value <i> lambda </i>
	@param nbr_to_observe <i> k </i>
	@return <em>p(k)</em>, i.e. the probability to observe k observations

*/
double Poisson_pmf(int nbr_to_observe, double expected_value) {
	double prob;
	prob = (exp(-expected_value)*pow(expected_value,nbr_to_observe))/(factorial(nbr_to_observe));
	return prob;
}

/** Factorial
	@param k
	@return factorial of <i>k</i>
*/
int factorial(int k){
	int ret;
	ret = 1;
	for (int j = 1; j <= k; j++){
	ret = ret*j;
	}
	return ret;
}
