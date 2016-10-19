#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include <fstream>
#include <string>
#include <vector>


using namespace std;


class RandomChains {
	private:
		Int_t run_type;
		Int_t nbr_pixels;

		Int_t lower_limit_alphas, upper_limit_alphas;
		Int_t lower_limit_escapes, upper_limit_escapes;
		Int_t lower_limit_implants, upper_limit_implants;

		TH1F* h_energy_reconstructed_beam_off_tot;
		TH1F* h_energy_reconstructed_beam_on_tot;

		//OBS the length of the arrays are declared here!
		TH1F* h_energy_pixel_reconstructed_beam_off[1024];
		TH1F* h_energy_pixel_reconstructed_beam_on[1024];
		Double_t fissions_pixels[1024];


		Int_t chain_length;
		vector<Int_t> beam_status;
		vector<char> decay_type;
		vector<Double_t> time_span;
		vector<array<Double_t, 1024>> rate;

		
		char cname[64], ctitle[64];

		void compute_random_chains(Int_t run_type);
		void generate_test_data();
		void read_experimental_data();
		void calculate_rates(Int_t run_type);
		void calculate_expected_nbr_random_chains();
		void set_chains(Int_t run_type);
		void rate_calc(char type, Int_t beam);

	public:
		RandomChains();
		~RandomChains();
		void plot_spectra();
		void run_main();
		

};

