#include "itensor/all.h"
#include "../headers/input.h"
#include "../headers/output.h"
#include "../headers/mcpeps.h"
#include "../headers/auxiliary_sampling.h"
#include "../headers/inner_sampling.h"
#include "../headers/pepsop.h"
#include <ctime>
#include <cmath>
#include <complex>


/*std::vector<itensor::ITensor> create_sz2_op(itensor::IndexSet &sites){
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
}*/

itensor::ITensor create_sz2_op(int site, itensor::IndexSet &sites){
	std::vector<itensor::ITensor> ops;
	double s = 0.5*(itensor::dim(sites(1))-1);
	itensor::Index site_index = sites(site+1);
	itensor::Index site_index_primed = itensor::prime(site_index);
	itensor::ITensor op(site_index, site_index_primed);
	for(int sz_index = 1; sz_index <= itensor::dim(site_index); sz_index ++){
		double sz = sz_index - s - 1.;
		op.set(site_index = sz_index, site_index_primed = sz_index, sz*sz);
	}
	return op;
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
	std::string log_file_name = input.testString("log_file", "");
	int standard_dims = input.testInteger("D", 2);
	int max_truncation_dims = input.testInteger("Dc", 4);
	int num_trials = input.testInteger("num_trials", 10000);
	std::string out_file_name = input.testString("out_file", "");
	int physical_dims = input.testInteger("physical_dims", 4);

	std::map<std::string, double> Jvals;
	Jvals["J1"] = input.testDouble("J1", 1);
	Jvals["Jz"] = input.testDouble("Jz", 1);
	Jvals["J2"] = input.testDouble("J2", 0);
	Jvals["Jd"] = input.testDouble("Jd", 0);

	Heisenberg H(Nx, Ny, physical_dims);
	H.set_Jz(Jvals);
	PEPSop HPEPO = H.toPEPSop();

	std::string version = "_";
	version += std::to_string(Nx) + "x" + std::to_string(Ny);
	version += "_D"+std::to_string(standard_dims);
	version += "_X"+std::to_string(max_truncation_dims);
	version += "_d"+std::to_string(physical_dims);
	version += "_" + std::to_string(num_trials) + "trials";
	if(log_file_name == "AUTO"){ //../../logs/sz2_test_{Nx}x{Ny}_D{D}_X{Chi}_d{d}_{}trials
		log_file_name = "../../logs/energy_test" + version;
	}


	int num_sites = Nx*Ny*UNIT_CELL_SIZE;

	std::vector<itensor::Index> sites_vector(num_sites);
	for(int i = 0; i < num_sites; i++){
		sites_vector[i] = itensor::Index(physical_dims,"Site,n="+std::to_string(i+1));
	}
	itensor::IndexSet sites(sites_vector);

	auto PEPS1 = MCKPEPS(sites, Nx, Ny, standard_dims, max_truncation_dims);
	MCKPEPS PEPS2 = PEPS1;
	PEPS2.prime();
	PEPS1.set_log_file(log_file_name);

	
	std::cerr << "Performing efficient inner product..." << std::endl;
	auto timestart = std::time(NULL);
	double inner_product = PEPS1.inner_product(PEPS2);
	std::cerr << "Evaluating matrix element..." << std::endl;
	double total_energy = HPEPO.eval(PEPS1, PEPS2)/inner_product;
	std::cerr << "Matrix element energy: " << total_energy << std::endl;
	double efficient_time = std::difftime(std::time(NULL), timestart);
	timestart = std::time(NULL);

	std::cerr << "Performing Monte Carlo inner product..." << std::endl;
	std::vector<double> wavefunctions;
	std::vector<double> values;
	std::vector<int> spin_config(num_sites, 0);
	Randomizer r;
	randomize_in_sector(spin_config, physical_dims, r.gen, r.dist);
	for(int trial = 0; trial < num_trials; trial++){
		std::cerr << "Sampling v direction..." << std::endl;
		sample_v_direction(PEPS1, spin_config, r);
		std::cerr << "Sampling s direction..." << std::endl;
		sample_s_direction(PEPS1, spin_config, r);
		std::cerr << "Sampling l direction..." << std::endl;
		double wavefn = sample_l_direction(PEPS1, spin_config, r);
		std::cerr << "Evaluating final value..." << std::endl;
		wavefunctions.push_back(wavefn);
		values.push_back(H.eval(spin_config, spin_config));
	}

	double mc_time = std::difftime(std::time(NULL), timestart);
	
	Output out;
	out.addInteger("NX", Nx);
	out.addInteger("NY", Ny);
	out.addInteger("D", standard_dims);
	out.addInteger("CHI", max_truncation_dims);
	out.addInteger("NUM_TRIALS", num_trials);
	out.addDouble("DIRECT_INNER_PRODUCT", inner_product);
	out.addDouble("DIRECT_ENERGY", total_energy);
	out.addVector("WAVEFUNCTIONS", wavefunctions);
	out.addVector("VALUES", values);

	if(out_file_name == "AUTO"){
		out_file_name = "../../out/energy_test" + version;
		std::ifstream out_file_cand(out_file_name + "_0");
		int out_file_number = 0;
		while(out_file_cand.good()){
			out_file_cand.close();
			out_file_number ++;
			out_file_cand = std::ifstream(out_file_name + "_" + std::to_string(out_file_number));
		}
		out_file_name = out_file_name + "_" + std::to_string(out_file_number);
		out_file_cand.close();
	}
	
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

	std::cerr << "Total Energy (inner product): " << total_energy << " (" << efficient_time << "s)" << std::endl;

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