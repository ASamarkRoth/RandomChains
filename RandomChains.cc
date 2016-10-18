#include <iostream>
#include <fstream>
#include "RandomChains.h"

using namespace std;


RandomChains::RandomChains() {
	nbr_pixels = 1024;

	run_type = -1;

	switch(run_type) {
		case 0 : 
			cout << "You are running default" << endl;
			compute_random_chains(0);
			break;
		case 1 : 
			cout << "User customised input: " << endl;
			compute_random_chains(1);
			break;
		case 2 : 
			cout << "You are running the sanity check" << endl;
			compute_random_chains(2);
			break;
		default : 
			cout << "What type of run? <0/1/2> \n 0: Reproduce article numbers [default],\n 1: User customised input on Lund experimental data & \n 2: Test method with a simple method " << endl;
			cin >> run_type;
	}
			
}

void RandomChains::plot_spectra() {
	cout << "What would you like to plot? \n" <<
		"0: Total energy spectrum \n" <<
		"1: Energy spectrum per pixel" << endl;
}

void RandomChains::compute_random_chains(Int_t run_type) {
	if(run_type == 2) {
		cout << "Generating test data ... " << endl;
		generate_test_data();
	}
	else {
		cout << "Reading experimental data ..." << endl;
		read_experimental_data();
	}
}

void RandomChains::generate_test_data() {
	cout << "Generating test data " << endl;
}

void RandomChains::read_experimental_data() {
	cout << "Reading data " << endl;
}

int main() {

	RandomChains* RC = new RandomChains();


	return 0;
}

