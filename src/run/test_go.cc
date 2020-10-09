#include "itensor/all.h"
#include "../headers/input.h"
#include "../headers/output.h"
#include "../headers/mcpeps.h"
#include "../headers/auxiliary_sampling.h"
#include "../headers/inner_sampling.h"
#include "../headers/pepsop.h"
#include "../headers/gradient_optimization.h"
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

int squared_distance(std::vector<int> vec1, std::vector<int> vec2){
	int sd = 0;
	if(vec1.size() != vec2.size()){std::cerr << "ERROR: VECTOR SIZES " << vec1.size() << " AND " << vec2.size() << " INCOMPATIBLE";}
	for(int i = 0; i < vec1.size(); i++){
		sd += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
	}
	return sd;
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
	double bias_dropoff = input.testDouble("bias_dropoff", 0.4);
	double update_size_init = input.testDouble("update_size_init", 0.05);
	int opt_steps = input.testInteger("optimization_steps", 100);

	std::map<std::string, double> Jvals;
	Jvals["J1"] = input.testDouble("J1", 1);
	Jvals["Jz"] = input.testDouble("Jz", 1);
	Jvals["J2"] = input.testDouble("J2", 0);
	Jvals["Jd"] = input.testDouble("Jd", 0);


	std::string version = "_";
	version += std::to_string(Nx) + "x" + std::to_string(Ny);
	version += "_J1"+std::to_string(static_cast<int>(Jvals["J1"]*1000000));
	version += "_Jz"+std::to_string(static_cast<int>(Jvals["Jz"]*1000000));
	version += "_J2"+std::to_string(static_cast<int>(Jvals["J2"]*1000000));
	if(Jvals["Jd"] != 0){
		version += "_Jd"+std::to_string(static_cast<int>(Jvals["Jd"]*1000000));
	}
	version += "_D"+std::to_string(standard_dims);
	version += "_X"+std::to_string(max_truncation_dims);
	version += "_d"+std::to_string(physical_dims);
	version += "_" + std::to_string(num_trials) + "trials";
	if(log_file_name == "AUTO"){ //../../logs/energy_test_{Nx}x{Ny}_J1_1000000_Jz_1000000_J2_0_Jd_0_D{D}_X{Chi}_d{d}_{}trials
		log_file_name = "../../logs/go_test" + version;
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
	double normsq = PEPS1.inner_product(PEPS2);
	PEPS1 /= std::sqrt(normsq);

	std::vector<int> spin_config(num_sites, 0);
	Randomizer r;
	randomize_in_sector(spin_config, physical_dims, r.gen, r.dist);
	std::vector<int> bias_config(spin_config);
	PEPS1.add_bias(bias_config, bias_dropoff);
	std::cerr << "Biased to spin config ";
	for(int sc : bias_config){std::cerr << sc << " ";}
	std::cerr << std::endl;
	
	randomize_in_sector(spin_config, physical_dims, r.gen, r.dist);

	PEPS1.set_log_file(log_file_name);

	
	std::cerr << "Performing optimization..." << std::endl;
	auto timestart = std::time(NULL);

	std::vector<double> energies;
	optimize(PEPS1, energies, Jvals, num_trials, update_size_init, opt_steps);

	double opt_time = std::difftime(std::time(NULL), timestart);
	timestart = std::time(NULL);

	std::cerr << "Performing Monte Carlo inner product..." << std::endl;
	std::vector<double> wavefunctions;
	std::vector<double> values;
	std::vector<int> squared_distances;
	
	Output out;
	out.addInteger("NX", Nx);
	out.addInteger("NY", Ny);
	out.addInteger("D", standard_dims);
	out.addInteger("CHI", max_truncation_dims);
	out.addInteger("NUM_TRIALS", num_trials);
	out.addDouble("BIAS_DROPOFF", bias_dropoff);
	out.addVector("BIAS_CONFIG", bias_config);
	out.addVector("ENERGIES", energies);

	if(out_file_name == "AUTO"){
		out_file_name = "../../out/go_test" + version;
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

	std::cerr << "Final Energy: " << energies[energies.size()-1] << " (" << opt_time << "s)" << std::endl;
	return 0;
}