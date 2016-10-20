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
		static const Int_t nbr_pixels = 1024;
		
		Int_t run_type;
		Double_t experiment_time;

		Int_t lower_limit_alphas, upper_limit_alphas;
		Int_t lower_limit_escapes, upper_limit_escapes;
		Int_t lower_limit_implants, upper_limit_implants;

		TH1F* h_energy_reconstructed_beam_off_tot;
		TH1F* h_energy_reconstructed_beam_on_tot;

		TH1F* h_energy_pixel_reconstructed_beam_off[nbr_pixels];
		TH1F* h_energy_pixel_reconstructed_beam_on[nbr_pixels];
		Double_t fissions_pixels[nbr_pixels];
		Int_t nbr_implants[nbr_pixels];


		vector<Int_t> chain_length;
		vector<Int_t> beam_status;
		vector<char> decay_type;
		vector<Double_t> time_span;
		vector<array<Double_t, nbr_pixels>> rate;
		vector<Double_t> nbr_expected_random_chains;

		Int_t eon;
		Int_t eoff;
		Int_t non;
		Int_t noff;
		Int_t aon;
		Int_t aoff;
		Int_t imps;
		Int_t fissions;
		
		char cname[64], ctitle[64];

		void prepare_data(Int_t run_type);
		void generate_test_data();
		void read_experimental_data();
		void calculate_implants();
		void calculate_rates(Int_t run_type);
		void calculate_expected_nbr_random_chains();
		void set_chains(Int_t run_type);
		void set_test_chains();
		void set_chains_from_input_file();
		void rate_calc(char type, Int_t beam);

	public:
		RandomChains();
		~RandomChains();
		void print_result();
		void print_test_result();
		void plot_spectra();
		void run_main();
		void dump_input_to_file();
		

};


Double_t Poisson_pmf(Int_t nbr_to_observe, Double_t expected_value);

Int_t factorial(Int_t k);

