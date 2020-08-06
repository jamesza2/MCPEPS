#include "itensor/all.h"
#include "../headers/input.h"
#include "../headers/output.h"
#include "../headers/mcpeps.h"
#include "../headers/inner_sampling.h"
#include <ctime>
#include <cmath>
#include <complex>


std::vector<itensor::ITensor> create_sz2_op(itensor::IndexSet &sites){
	std::vector<itensor::ITensor> ops;
	double s = 0.5*(itensor::dim(sites(1))-1);
	for(int i = 1; i <= itensor::length(sites); i++){
		itensor::Index site_index = sites(i);
		itensor::Index site_index_primed = itensor::prime(site_index);
		itensor::ITensor op(site_index, site_index_primed);
		for(int sz_index = 1; sz_index <= itensor::dim(site_index); sz_index ++){
			double sz = sz_index - s - 1.;
			op.set(site_index = sz_index, site_index_primed = sz_index, sz*sz);
		}
		ops.push_back(op);
	}
	return ops;
}

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
	int physical_dims = input.testInteger("physical_dims", 2);

	int num_sites = Nx*Ny*UNIT_CELL_SIZE;

	std::vector<itensor::Index> sites_vector(num_sites);
	for(int i = 0; i < num_sites; i++){
		sites_vector[i] = itensor::Index(physical_dims,"Site,n="+std::to_string(i+1));
	}
	itensor::IndexSet sites(sites_vector);
	auto sz2_op = create_sz2_op(sites);

	auto PEPS1 = MCKPEPS(sites, Nx, Ny, standard_dims, max_truncation_dims);
	MCKPEPS PEPS2 = PEPS1;
	PEPS2.prime();
	PEPS2.apply_spinop(sz2_op);
	//auto PEPS2 = MCKPEPS(sites, Nx, Ny, 1, max_truncation_dims); //Random product state
	PEPS1.set_log_file(log_file);
	auto timestart = std::time(NULL);
	std::cerr << "Performing efficient inner product..." << std::endl;
	double inner_product = PEPS1.inner_product(PEPS2);
	double efficient_time = std::difftime(std::time(NULL), timestart);
	timestart = std::time(NULL);

	std::cerr << "Performing Monte Carlo inner product..." << std::endl;
	std::vector<double> wavefunctions;
	std::vector<double> values;
	mc_sz2(PEPS1, wavefunctions, values, num_trials);

	double mc_time = std::difftime(std::time(NULL), timestart);
	
	Output out;
	out.addInteger("NX", Nx);
	out.addInteger("NY", Ny);
	out.addInteger("D", standard_dims);
	out.addInteger("CHI", max_truncation_dims);
	out.addInteger("NUM_TRIALS", num_trials);
	out.addDouble("DIRECT_INNER_PRODUCT", inner_product);
	out.addVector("WAVEFUNCTIONS", wavefunctions);
	out.addVector("VALUES", values);
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

	std::cerr << "Total Sz^2 (inner product): " << inner_product << " (" << efficient_time << "s)" << std::endl;

	int num_wavefunctions = values.size();
	double later_half_average = 0;
	for(int i = num_wavefunctions/2; i < num_wavefunctions; i++){
		later_half_average += values[i];
	}
	later_half_average /= (num_wavefunctions - num_wavefunctions/2);
	std::cerr << "Later-Half Average: " << later_half_average << " (" << mc_time << "s)" << std::endl;
	//std::cerr << "Brute Force Inner Product (Old Method): " << brute_force_inner_product_old << " (" << brute_force_time_old << "s)" << std::endl;
	return 0;
}