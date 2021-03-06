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
#include <iomanip>


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
				else{
					product *= PEPS1._site_tensors[i][j][k];
				}
			}
		}
	}
	return product;
}

std::vector<double> generate_update_sizes(double update_size_init, double update_size_min, double decay, int num_trials, int num_stay=0){
	std::vector<double> update_sizes_to_return;
	for(int stay = 0; stay < num_stay; stay++){
		update_sizes_to_return.push_back(update_size_init);
	}
	double update_size = update_size_init;
	for(int trial = num_stay; trial < num_trials; trial++){
		update_sizes_to_return.push_back(update_size);
		update_size = std::max(update_size*decay, update_size_min);
	}
	return update_sizes_to_return;
}

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
	double update_size_min = input.testDouble("update_size_min", 0.001);
	int opt_steps = input.testInteger("optimization_steps", 100);
	double decay = input.testDouble("decay", 0.96);
	std::string gradient_file_name = input.testString("gradient_file", ""); //If it exists, print the gradients to that output
	int opt_steps_stay = input.testInteger("optimization_steps_no_decay", static_cast<int>(opt_steps*0.2));

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

	//std::cerr << "Evaluated input" << std::endl;

	int num_sites = Nx*Ny*UNIT_CELL_SIZE;

	std::vector<itensor::Index> sites_vector(num_sites);
	for(int i = 0; i < num_sites; i++){
		sites_vector[i] = itensor::Index(physical_dims,"Site,n="+std::to_string(i+1));
	}
	itensor::IndexSet sites(sites_vector);

	//std::cerr << "Creating " << Nx << "x" << Ny << " PEPS...";
	auto PEPS1 = MCKPEPS(sites, Nx, Ny, standard_dims, max_truncation_dims);
	//std::cerr << PEPS1.Nx() << "x" << PEPS1.Ny() << "x" << UNIT_CELL_SIZE << "...";

	std::vector<int> spin_config(num_sites, 0);
	Randomizer r;
	//std::cerr << "Creating randomized bias...";
	randomize_in_sector(spin_config, physical_dims, r.gen, r.dist);
	std::vector<int> bias_config(spin_config);
	//std::cerr << "Adding bias...";
	PEPS1.add_bias(bias_config, bias_dropoff);
	std::cerr << "Biased to spin config ";
	for(int sc : bias_config){std::cerr << sc << " ";}
	std::cerr << std::endl;
	
	randomize_in_sector(spin_config, physical_dims, r.gen, r.dist);

	PEPS1.set_log_file(log_file_name);
	MCKPEPS PEPS0 = PEPS1;

	MCKPEPS PEPS2 = PEPS1;
	PEPS2.prime();
	double normsq = PEPS1.inner_product(PEPS2);
	//std::cerr << "Dividing by norm " << normsq << std::endl;
	PEPS1 /= std::sqrt(normsq);
	
	std::cerr << "Performing optimization..." << std::endl;
	auto timestart = std::time(NULL);

	std::vector<double> energies;
	std::vector<double> fidelities;
	std::vector<double> update_sizes = generate_update_sizes(update_size_init, update_size_min, decay, opt_steps, opt_steps_stay);

	itensor::Args optimize_args = itensor::Args("GradientFile=", "");

	optimize(PEPS1, energies, update_sizes, Jvals, num_trials, opt_steps, fidelities, optimize_args);

	double opt_time = std::difftime(std::time(NULL), timestart);
	timestart = std::time(NULL);

	std::cerr << "Performing Monte Carlo inner product..." << std::endl;
	std::vector<double> wavefunctions;
	std::vector<double> values;
	std::vector<int> squared_distances;
	std::vector<double> szs;
	MCKPEPS PEPS3 = PEPS1;
	PEPS3.prime();
	for(int site = 0; site < num_sites; site++){
		PEPSop single_site = singleSiteSz(site);
		double szval = single_site.eval(PEPS1, PEPS2);
		szs.push_back(szval);
	}
	
	Output out;
	out.addInteger("NX", Nx);
	out.addInteger("NY", Ny);
	out.addInteger("D", standard_dims);
	out.addInteger("CHI", max_truncation_dims);
	out.addInteger("NUM_TRIALS", num_trials);
	out.addDouble("BIAS_DROPOFF", bias_dropoff);
	out.addVector("BIAS_CONFIG", bias_config);
	out.addVector("UPDATE_SIZES", update_sizes);
	out.addVector("ENERGIES", energies);
	out.addVector("GRADIENT_FIDELITIES", fidelities);
	out.addVector("SZS", szs);

	//Test gradient optimization by finding the gradient from a direct contraction
	bool direct_test = input.testBool("direct_test", false);
	if(direct_test){
		PEPS1 = PEPS0;
		PEPS2 = PEPS1;
		PEPS2.prime();
		std::cerr << "Computing direct gradient..." << std::endl;
		Heisenberg H(Nx, Ny, physical_dims);
		H.set_J(Jvals);
		PEPSop HPEPO = H.toPEPSop();
		std::vector<double> direct_energies;
		double update_size = update_size_init;
		for(int step = 0; step < opt_steps; step++){
			normsq = PEPS1.inner_product(PEPS2);
			//std::cerr << "Norm: " << normsq << "...";
			std::vector<itensor::ITensor> direct_grad(num_sites);
			double energy = 0;
			ArbitraryPEPS contracted = PEPS1.combine(PEPS2);
			std::vector<itensor::ITensor> me1(num_sites); //Matrix elements of <psi|H|env>
			std::vector<itensor::ITensor> me2 = contracted.environments(); //Matrix elements of <psi|env>
			for(int site = 0; site < num_sites; site++){
				me2[site] *= PEPS2.site_tensor(site);
			}
			//std::cerr << "Energy terms: ";

			for(Term t : HPEPO.terms){
				MCKPEPS PEPS_applied = PEPS2;
				t.apply(PEPS_applied);
				double energy_me = t.eval(PEPS1, PEPS2);
				//std::cerr << energy_me << " ";
				double energy_part = energy_me/normsq;
				energy += energy_part;
				ArbitraryPEPS contracted_applied = PEPS1.combine(PEPS_applied);
				std::vector<itensor::ITensor> capp_envs = contracted_applied.environments();
				for(int site = 0; site < num_sites; site++){
					me1[site] += capp_envs[site]*PEPS_applied.site_tensor(site);
				}
			}
			//std::cerr << std::endl;
			for(int site = 0; site < num_sites; site++){
				itensor::ITensor direct_gradient = (me1[site] - energy*me2[site])*2/normsq;
				direct_gradient = signelts(direct_gradient);
				PEPS1.site_tensor(site) -= update_size*r.rand()*direct_gradient;
			}
			PEPS1 /= std::sqrt(normsq);
			PEPS2 = PEPS1;
			PEPS2.prime();
			//PEPS2.print_self("PEPS2 after update");
			update_size = std::max(update_size*0.95, update_size_min);
			direct_energies.push_back(energy);
			std::cerr << "DIRECT STEP #" << step+1 << " HAS ENERGY " << energy << " AND NORM " << normsq << std::endl;
		}
		out.addVector("DIRECT_ENERGIES", direct_energies);
	}
	

	//Test gradient optimization by finding the gradient from a brute force method
	bool brute_force_test = input.testBool("brute_force_test", false);
	if(brute_force_test){
		PEPS1 = PEPS0;
		PEPS2 = PEPS1;
		PEPS2.prime();
		std::cerr << "Computing brute force gradient..." << std::endl;
		Heisenberg H(Nx, Ny, physical_dims);
		H.set_J(Jvals);
		PEPSop HPEPO = H.toPEPSop();
		std::vector<double> brute_force_energies;
		double update_size = update_size_init;
		for(int step = 0; step < opt_steps; step++){
			normsq = PEPS1.inner_product(PEPS2);
			std::vector<itensor::ITensor> exact_grad(num_sites);
			double energy = 0;
			for(Term t : HPEPO.terms){
				MCKPEPS PEPS_applied = PEPS2;
				t.apply(PEPS_applied);
				/*if(t.factor == 0){
					if(itensor::norm(PEPS_applied._site_tensors[0][0][0]) != 0){
						std::cerr << "Zero peps site norm: " << itensor::norm(PEPS_applied._site_tensors[0][0][0]);
					}
				}*/
				double energy_part = t.eval(PEPS1, PEPS2)/normsq;
				//std::cerr << "Energy part of " << t.to_string() << ": " << energy_part << std::endl;
				energy += energy_part;
				for(int site = 0; site < num_sites; site++){
					auto [i,j,k] = PEPS1.position_of_site(site);
					auto me1 = incomplete_inner(PEPS_applied, PEPS1, i, j, k);
					auto me2 = incomplete_inner(PEPS2, PEPS1, i, j, k);
					me2 *= energy_part;
					/*if(site==0){
						std::cerr << "Term " << t.to_string() << ": " <<std::endl;
						PrintData(me1);
						PrintData(me2);
						std::cerr << "Real norm: " << itensor::norm(me2) << std::endl;
						PrintData(PEPS1._site_tensors[0][0][0]);
					}*/
					exact_grad[site] += (me1-me2);
				}
				//PEPS_applied.print_self("Applied PEPS");
				//PEPS2.print_self("PEPS2");
			}
			//PrintData(exact_grad[0]);
			//PrintData(PEPS1._site_tensors[0][0][0]);
			for(int site = 0; site < num_sites; site++){
				auto [i,j,k] = PEPS1.position_of_site(site);
				if(itensor::norm(exact_grad[site]*PEPS1._site_tensors[i][j][k]) > 0.000001){
					std::cerr << "Exact grad product: " << itensor::norm(exact_grad[site]*PEPS1._site_tensors[i][j][k]) << std::endl;
				}
				exact_grad[site] = signelts(exact_grad[site]);
				PEPS1._site_tensors[i][j][k] -= update_size*r.rand()*exact_grad[site];
			}
			PEPS2 = PEPS1;
			PEPS2.prime();
			//PEPS2.print_self("PEPS2 after update");
			update_size *= 0.97;
			brute_force_energies.push_back(energy);
			std::cerr << "BFSTEP #" << step+1 << " HAS ENERGY " << energy << " AND NORM " << normsq << std::endl;
		}
		out.addVector("BRUTE_FORCE_ENERGIES", brute_force_energies);
	}

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