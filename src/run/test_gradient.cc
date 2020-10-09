#include "itensor/all.h"
#include "../headers/input.h"
#include "../headers/output.h"
#include "../headers/mcpeps.h"
#include "../headers/gradient_optimization.h"
#include "../headers/pepsop.h"
#include "../headers/inner_sampling.h"
#include <iomanip>
#include <ctime>
#include <cmath>
#include <complex>


//Brute force inner product between PEPS1 and PEPS2, except for the tensor [i,j,k] on PEPS2
itensor::ITensor incomplete_inner(MCKPEPS &PEPS1, MCKPEPS &PEPS2, int i_omit, int j_omit, int k_omit){
	itensor::ITensor product(1);
	for(int i = 0; i < PEPS1.Nx(); i++){
		for(int j = 0; j < PEPS1.Ny(); j++){
			for(int k = 0; k < UNIT_CELL_SIZE; k++){
				if((i != i_omit) || (j != j_omit) || (k != k_omit)){
					product *= PEPS1._site_tensors[i][j][k];
					product *= PEPS2._site_tensors[i][j][k];
				}
			}
		}
	}
	return product;
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

	int Nx = input.testInteger("Nx", 3);
	int Ny = input.testInteger("Ny",3);
	std::string log_file_name = input.testString("log_file", "");
	int standard_dims = input.testInteger("D", 1);
	int max_truncation_dims = input.testInteger("Dc", 1);
	int num_trials = input.testInteger("num_trials", 10000);
	int physical_dims = input.testInteger("physical_dims", 4);
	double bias_dropoff = input.testDouble("bias_dropoff", 0.4);
	double update_size_init = input.testDouble("update_size_init", 0.05);
	int opt_steps = input.testInteger("optimization_steps", 100);

	int target_i = input.testInteger("target_i",1);
	int target_j = input.testInteger("target_j",1);
	int target_k = input.testInteger("target_k",1);

	std::map<std::string, double> Jvals;
	Jvals["J1"] = input.testDouble("J1", 1);
	Jvals["Jz"] = input.testDouble("Jz", 1);
	Jvals["J2"] = input.testDouble("J2", 0);
	Jvals["Jd"] = input.testDouble("Jd", 0);

	int print_precision = input.testInteger("precision", -1);
	if(print_precision != -1){
		std::cout << std::setprecision(print_precision);
		std::cerr << std::setprecision(print_precision);
	}


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

	Heisenberg H(Nx, Ny, physical_dims);
	H.set_J(Jvals);
	PEPSop HPEPO = H.toPEPSop();

	auto PEPS1 = MCKPEPS(sites, Nx, Ny, standard_dims, max_truncation_dims);

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

	MCKPEPS PEPS2 = PEPS1;
	PEPS2.prime();
	double normsq = PEPS1.inner_product(PEPS2);
	//Compute the exact gradient
	std::cerr << "Computing exact gradient..." << std::endl;
	itensor::ITensor exact_grad;
	for(Term t : HPEPO.terms){
		MCKPEPS PEPS_applied = PEPS2;
		t.apply(PEPS_applied);
		auto me1 = incomplete_inner(PEPS_applied, PEPS1, target_i, target_j, target_k);
		auto me2 = incomplete_inner(PEPS2, PEPS1, target_i, target_j, target_k);
		me2 *= t.eval(PEPS1, PEPS2);
		me2 /= normsq;
		exact_grad += (me1+me2);
	}
	exact_grad *= (2./normsq);


	//Get the approximated gradient
	std::cerr << "Computing approximate gradient..." << std::endl;
	std::vector<itensor::ITensor> Delta(PEPS1.size());
	std::vector<itensor::ITensor> DeltaE(PEPS1.size());
	double E = 0;
	SpinConfigPEPS scp(PEPS1, spin_config, 1);
	NoSitePEPS nsp = PEPS1.contract(scp);
	for(int sample = 0; sample < num_trials; sample++){
		get_sample(PEPS1, nsp, spin_config, H, Delta, DeltaE, E, update_size);
	}
	E /= num_trials;
	double grads_factor = 2./num_trials;
	int target_site = PEPS1.site_index_from_position(target_i, target_j, target_k);
	itensor::ITensor approx_grad = DeltaE.at(target_site)*grads_factor - Delta.at(target_site)*grads_factor*E;
	
	
	PrintData(exact_grad);
	PrintData(approx_grad);
	return 0;
}