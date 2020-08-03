#include "itensor/all.h"
#include "../headers/input.h"
#include "../headers/output.h"
#include "../headers/mcpeps.h"
#include "../headers/inner_sampling.h"
#include <ctime>
#include <cmath>
#include <complex>


int main(int argc, char *argv[]){
	int target_argc = 2;
	if(argc != target_argc){
		std::cerr << "Please provide an input file" << std::endl;
		return 1;
	}

	//Take inputs
	std::ifstream input_file_reader(argv[1]);
	if(!input_file_reader.is_open()){
		std::cerr <<"FILENAME " << argv[1] << " NOT FOUND" << endl;
		return 2;
	}

	InputClass input;
	input.Read(input_file_reader);

	int Nx = input.testInteger("Nx", 2);
	int Ny = input.testInteger("Ny", 2);
	std::string log_file = input.testString("log_file", "");
	int standard_dims = input.testInteger("D", 2);
	int max_truncation_dims = input.testInteger("Dc", 4);
	int num_trials = input.testInteger("num_trials", 10000);
	std::string out_file_name = input.testString("out_file", "");

	int num_sites = Nx*Ny*UNIT_CELL_SIZE;
	std::vector<itensor::Index> sites_vector(num_sites);
	for(int i = 0; i < num_sites; i++){
		sites_vector[i] = itensor::Index(2);
	}
	itensor::IndexSet sites(sites_vector);
	auto PEPS1 = MCKPEPS(sites, Nx, Ny, standard_dims, max_truncation_dims);
	MCKPEPS PEPS2 = PEPS1;
	PEPS2.prime();
	//auto PEPS2 = MCKPEPS(sites, Nx, Ny, 1, max_truncation_dims); //Random product state
	PEPS1.set_log_file(log_file);
	auto timestart = std::time(NULL);
	std::cerr << "Performing efficient inner product..." << std::endl;
	double inner_product = PEPS1.inner_product(PEPS2);
	double efficient_time = std::difftime(std::time(NULL), timestart);
	timestart = std::time(NULL);

	std::cerr << "Performing Monte Carlo inner product..." << std::endl;
	std::vector<double> wavefunctions;
	mc_norm(PEPS1, wavefunctions, num_trials);

	double mc_time = std::difftime(std::time(NULL), timestart);
	
	Output out;
	out.addInteger("NX", Nx);
	out.addInteger("NY", Ny);
	out.addInteger("D", standard_dims);
	out.addInteger("CHI", max_truncation_dims);
	out.addInteger("NUM_TRIALS", num_trials);
	out.addDouble("DIRECT_INNER_PRODUCT", inner_product);
	out.addVector("WAVEFUNCTIONS", wavefunctions);
	out.writeOutput(out_file_name);
	/*std::ofstream out_file(out_file_name);
	out_file << "NX: " << Nx;
	out_file << "\nNY: " << Ny;
	out_file << "\nD: " << D;
	out_file << "\nCHI: " << Dc;
	out_file << "\nNUM_TRIALS: " << num_trials;
	out_file << "\nDIRECT_INNER_PRODUCT: " << inner_product;
	out_file << "\nWAVEFUNCTIONS:"
	for(double wfn : wavefunctions){
		out_file << " " << wfn;
	}

	out_file.close();*/

	std::cerr << "Inner Product: " << inner_product << " (" << efficient_time << "s)" << std::endl;

	int num_wavefunctions = wavefunctions.size();
	double later_half_average = 0;
	for(int i = num_wavefunctions/2; i < num_wavefunctions; i++){
		later_half_average += wavefunctions[i];
	}
	later_half_average /= (num_wavefunctions - num_wavefunctions/2);
	std::cerr << "Later-Half Average: " << later_half_average << " (" << mc_time << "s)" << std::endl;
	//std::cerr << "Brute Force Inner Product (Old Method): " << brute_force_inner_product_old << " (" << brute_force_time_old << "s)" << std::endl;
	return 0;
}