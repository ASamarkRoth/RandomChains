/** @file RandomChains.h
@author Anton Roth (anton.roth@nuclear.lu.se)
@brief Header file for RandomChains.cc
*/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


class RandomChains {
	private:
		const int nbr_pixels; 
		const int nbr_bins;

		string folder_data;
		
		//Indicates the type of run (0,1 or 2)
		int run_type;

		//The duration of the experiment in s
		double experiment_time;

		//False if no pure beam ON spectra could be read in
		bool pure_beam = true;

		//bin limits for the different decay types
		int lower_limit_alphas, upper_limit_alphas;
		int lower_limit_escapes, upper_limit_escapes;
		int lower_limit_implants, upper_limit_implants;

		//The spectra are stored in these 2D vectors
		vector< vector<int> > data_beam_on;
		vector< vector<int> > data_reconstructed_beam_on;
		vector< vector<int> > data_reconstructed_beam_off;

		//pixels with fissions
		vector<double> fissions_pixels;

		//Number of implants for every pixel
		vector<int> nbr_implants;

		//Chain/chains characteristics
		vector<int> chain_length;
		vector<int> beam_status;
		vector<char> decay_type;
		vector<double> time_span;

		//Background rates and expected number of random chains per pixel
		vector< vector<double> > rate;
		vector<double> nbr_expected_random_chains;

		//Help variables to generate the test data and for verification
		int eon;
		int eoff;
		int non;
		int noff;
		int aon;
		int aoff;
		int imps;
		int fissions;
		
		char cname[64], ctitle[64];

		//all methods are described in "RandomChains.cc"
		void read_exp_file(string file_name);
		void generate_test_data();
		void calculate_implants();
		void calculate_rates();
		void calculate_expected_nbr_random_chains();
		void set_test_chains();
		void set_article_chains();
		void set_chains_from_input_file(string input_file);
		void rate_calc(char type, int beam);

	public:
		RandomChains(int pixels=1024, int bins=4096, string folder="Lund_data");
		void ReadExperimentalData();
		void SetDecayChains(string input_chains="");
		void Run();
		~RandomChains();
		void print_result();
		void print_test_result();
		void dump_input_to_file();
		

};

/* Mathematical functions (non-member) */
double Poisson_pmf(int nbr_to_observe, double expected_value);

int factorial(int k);

