#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include <fstream>
#include <string>
#include <vector>


using namespace std;


//template <unsigned const int nbr_pixelss, unsigned const int nbr_bins>
class RandomChains {
	private:
		/* Here the number of pixels is set!!!!! */
		const int nbr_pixels; 
		const int nbr_bins;

		string folder_data;
		
		int run_type;
		double experiment_time;

		bool pure_beam = true;

		//bin limits for the different decay types
		int lower_limit_alphas, upper_limit_alphas;
		int lower_limit_escapes, upper_limit_escapes;
		int lower_limit_implants, upper_limit_implants;

		//Spectra, summed all pixels, beam ON/OFF
		TH1F* h_energy_reconstructed_beam_off_tot;
		TH1F* h_energy_reconstructed_beam_on_tot;

		//Spectrum for every pixel, beam ON/OFF
		/*
		TH1F* h_energy_pixel_reconstructed_beam_off[nbr_pixels];
		TH1F* h_energy_pixel_reconstructed_beam_on[nbr_pixels];
		TH1F* h_energy_pixel_beam_on[nbr_pixels];
		*/

		//int* data_beam_on;
		vector<vector<int>> data_beam_on;
		vector<vector<int>> data_reconstructed_beam_on;
		vector<vector<int>> data_reconstructed_beam_off;

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

		void read_exp_file(string file_name);
		void prepare_data(int run_type);
		void generate_test_data();
		void read_experimental_data();
		void calculate_implants();
		void calculate_rates();
		void calculate_expected_nbr_random_chains();
		void set_chains(int run_type, string input_file);
		void set_test_chains();
		void set_article_chains();
		void set_chains_from_input_file(string input_file);
		void rate_calc(char type, int beam);

	public:
		RandomChains(int pixels=1024, int bins=4096, string folder="Lund_data");
		//RandomChains();
		void ReadExperimentalData();
		void Run(string input_file="");
		~RandomChains();
		void print_result();
		void print_test_result();
		//void plot_spectra();
		void run_main();
		void dump_input_to_file();
		

};

/* Mathematical functions (non-member) */
double Poisson_pmf(int nbr_to_observe, float expected_value);

int factorial(int k);

void insert_blank();

