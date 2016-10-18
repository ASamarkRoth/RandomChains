#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include <fstream>
#include <string>


using namespace std;


class RandomChains {
	private:
		int run_type;
		Int_t nbr_pixels;

		TH1F* h_energy_pixel_reconstructed_beam_off_tot;
		TH1F* h_energy_pixel_reconstructed__beam_on_tot;

		//TH1F* h_energy_pixel_reconstructed__beam_off[nbr_pixels];
		//TH1F* h_energy_pixel_reconstructed_beam_on[nbr_pixels];

		void compute_random_chains(Int_t run_type);
		void generate_test_data();
		void read_experimental_data();

	public:
		RandomChains();
		void plot_spectra();
		

};

