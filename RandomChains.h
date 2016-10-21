#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include <fstream>
#include <string>
#include <vector>


using namespace std;


class RandomChains {
	private:
		/* Here the number of pixels is set!!!!! */
		static const int nbr_pixels = 1024;
		
		int run_type;
		double experiment_time;

		int lower_limit_alphas, upper_limit_alphas;
		int lower_limit_escapes, upper_limit_escapes;
		int lower_limit_implants, upper_limit_implants;

		TH1F* h_energy_reconstructed_beam_off_tot;
		TH1F* h_energy_reconstructed_beam_on_tot;

		TH1F* h_energy_pixel_reconstructed_beam_off[nbr_pixels];
		TH1F* h_energy_pixel_reconstructed_beam_on[nbr_pixels];
		double fissions_pixels[nbr_pixels];
		int nbr_implants[nbr_pixels];


		vector<int> chain_length;
		vector<int> beam_status;
		vector<char> decay_type;
		vector<double> time_span;
		vector<array<long double, nbr_pixels>> rate;
		vector<double> nbr_expected_random_chains;

		int eon;
		int eoff;
		int non;
		int noff;
		int aon;
		int aoff;
		int imps;
		int fissions;
		
		char cname[64], ctitle[64];

		void prepare_data(int run_type);
		void generate_test_data();
		void read_experimental_data();
		void calculate_implants();
		void calculate_rates(int run_type);
		void calculate_expected_nbr_random_chains();
		void set_chains(int run_type);
		void set_test_chains();
		void set_article_chains();
		void set_chains_from_input_file();
		void rate_calc(char type, int beam);

	public:
		RandomChains();
		~RandomChains();
		void print_result();
		void print_test_result();
		void plot_spectra();
		void run_main();
		void dump_input_to_file();
		

};


double Poisson_pmf(int nbr_to_observe, double expected_value);

int factorial(int k);

void insert_blank();

