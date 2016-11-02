/*!
@file RandomChains.cc
@author Anton Roth (anton.roth@nuclear.lu.se)

@brief From this file the program should be controlled.


*/
#include "RandomChains.h"

/** This function should be invoked to run the program RandomChains.
	@see RandomChains::RandomChains()
	@see RandomChains::SetDecayChains()
	@see RandomChains::Run()
*/

int main() {

	//This is how the default run is made (i.e. on Lund experimental data)
	RandomChains* RC = new RandomChains();
	RC->SetDecayChains();
	RC->Run();

	/* POSSIBLE MODIFICATIONS FOR OTHER EXPERIMENTAL DATA OR OTHER TYPES OF DECAY CHAINS CAN BE MADE ON THE BASIS OF THE FOLLOWING CODE */ 
	/*
	int nbr_of_pixels = 1024;
	int nbr_of_bins = 4096;
	string folder_with_data = "Lund_data";
	RandomChains* RC_mod = new RandomChains(nbr_of_pixels, nbr_of_bins, folder_with_data);
	
	string chains_input_file = "dump_article.txt";
	//Giving an input to the RandomChains::Run method relieves the user from giving any inputs during a run
	RC_mod->SetDecayChains(chains_input_file);
	RC_mod->Run();
	*/
	
	return 0;
}

