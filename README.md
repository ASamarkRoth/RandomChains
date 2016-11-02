#Program RandomChains

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

1. Read experimental data. Achieved with the constructor:
   RandomChains::RandomChains(int pixels, int bins, string folder)

2. Read decay chain/chains characteristics data. Achieved with the method:
   RandomChains::SetDecayChains(string input_chains)

3. Compute the number of expected random chains for the given
   decay chain/chains with the experimental data. Achieved
   with the method: RandomChains::Run()

<h4>How to run RandomChains? </h4>
The program is preferably controlled from the file <tt>run_file.cc</tt>.

When located in the directory of RandomChains type the following in the terminal to run the program:

1. <tt> make </tt>

2. <tt> ./run_file </tt>

<h4>System requirements</h4>
The program has been successfully run on the following systems:
- Linux Ubuntu 14.04 with g++ version 4.8.4 and std=c++11.

<h4>Documentation</h4>
A rather detailed documentation of the program can be generated with <b>Doxygen</b>, by typing <tt>doxygen</tt> in the program folder, and then be found in the created folder <i>documentation/</i>. An '.html' file is found at <i>documentation/html/index.html</i> and the <i>.pdf</i> is generated with make in the same folder as the follwoing file is created <i>documentation/latex/runman.pdf</i>.

OBS: A pregenerated documentation can be found in <i>runman.pdf</i>.
