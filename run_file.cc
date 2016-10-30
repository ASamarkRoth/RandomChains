#include "RandomChains.h"


/*From here the RandomChains computation is controlled. */
int main() {

	//This is how the default run is made (i.e. on Lund experimental data)
	/*
	RandomChains* RC = new RandomChains();
	RC->Run();
	*/

	/* POSSIBLE MODIFICATIONS FOR OTHER EXPERIMENTAL DATA OR OTHER TYPES OF DECAY CHAINS CAN BE MADE ON THE BASIS OF THE FOLLOWING CODE */ 
	int nbr_of_pixels = 1024;
	int nbr_of_bins = 4096;
	string folder_with_data = "Lund_data";
	RandomChains* RC_mod = new RandomChains(nbr_of_pixels, nbr_of_bins, folder_with_data);
	
	string chains_input_file = "dump_article.txt";
	//Giving an input to the RandomChains::Run method relieves the user from giving any inputs during a run
	RC_mod->Run(chains_input_file);
	
	return 0;
}

