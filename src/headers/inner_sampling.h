#ifndef mccalculator
#define mccalculator

#include "mcpeps.h"
#include "neighbors.h"
#include <random>


double wavefunction(std::vector<int> &spin_config, MCKPEPS &state){
	SpinConfigPEPS spins(state);
	spins.set_spins(spin_config);
	return state.inner_product(spins);
}

void randomize(std::vector<int> &spin_config, int spin_max, std::mt19937 &generator, std::uniform_real_distribution<double> &distribution){
	for(int spin_index = 0; spin_index < spin_config.size(); spin_index ++){
		spin_config[spin_index] = std::floor(spin_max*distribution(generator));
	}
}
//spin_max = number of possible spin states, i.e. 2S+1
void randomize_in_sector(std::vector<int> &spin_config, int spin_max, std::mt19937 &generator, std::uniform_real_distribution<double> &distribution, int target_sz = 0){
	int total_sz = 0;
	for(int spin_index = 0; spin_index < spin_config.size(); spin_index ++){
		int spin_choice = std::floor(spin_max*distribution(generator));
		total_sz += 2*spin_choice;
		spin_config[spin_index] = spin_choice;
	}
	total_sz -= (spin_max-1)*spin_config.size();
	while((total_sz != 2*target_sz) && (total_sz != 2*target_sz+1)){// If the doubled total Sz isn't target_sz or the target_sz+1, downgrade spins until you reach the correct value
		if(total_sz - 2*target_sz > 0){
			int spin_to_flip = std::floor(spin_config.size()*distribution(generator));
			if(spin_config[spin_to_flip] > 0){
				spin_config[spin_to_flip] -= 1;
				total_sz -= 2;
			}
		}
		else{
			int spin_to_flip = std::floor(spin_config.size()*distribution(generator));
			if(spin_config[spin_to_flip] < spin_max - 1){
				spin_config[spin_to_flip] += 1;
				total_sz += 2;
			}
		}
	}
}

//Takes random NN pair of 2 spins (not necessarily pairs) and increments one spin while decremening the other
//returns ratio of total choices from old config over new config
double flip_spins(Neighbors &bonds, std::vector<int> &spin_config, int spin_max, std::mt19937 &generator, std::uniform_real_distribution<double> &distribution){
	int num_sites = spin_config.size();
	double num_choices = 0;
	for(int site_1 = 0; site_1 < num_sites; site_1 ++){
		for(int site_2 : bonds.nn_at(site_1)){
			if((spin_config[site_1] > 0) && (spin_config[site_2] < spin_max-1)){//Can flip by decrementing site_1, incrementing site_2
				num_choices += 1;
			}
			if((spin_config[site_1] < spin_max-1) && (spin_config[site_2] > 0)){//Can flip by incrementing site_1, decrementing site_2
				num_choices += 1;
			}
		}
	}
	double new_num_choices = num_choices;
	while(true){
		int site_choice_1 = std::floor(num_sites*distribution(generator));
		auto bonds_at_site_1 = bonds.nn_at(site_choice_1);
		int neighbor_choice_index = std::floor(bonds_at_site_1.size()*distribution(generator));
		int site_choice_2 = bonds_at_site_1[neighbor_choice_index];
		std::cout << "Testing flip at " << site_choice_1 << ", " << site_choice_2 << std::endl;
		/*if(site_choice_1 == site_choice_2){
			pairs_found -= 1;
			continue;
		}*/
		int old_sz_1 = spin_config[site_choice_1];
		int old_sz_2 = spin_config[site_choice_2];
		int total_site = old_sz_1 + old_sz_2;
		//If total is 0,1,2S-1 or 2S, only one possible way to split config
		if((total_site < 2) || (total_site > 2*spin_max-4)){
			continue;
		}
		int up_or_down = std::floor(2*distribution(generator))*2-1; //-1 means decrement spin_1, +1 means increment spin_1
		int new_sz_1 = spin_config[site_choice_1] + up_or_down;
		int new_sz_2 = spin_config[site_choice_2] - up_or_down;
		if((new_sz_1 < 0) || (new_sz_2 < 0) || (new_sz_1 >= spin_max) || (new_sz_2 >= spin_max)){
			continue;
		}
		spin_config[site_choice_1] = new_sz_1;
		spin_config[site_choice_2] = new_sz_2;
		if((old_sz_1 > 0) && (old_sz_2 < spin_max-1)){
			new_num_choices -= 1;
		}
		if((old_sz_1 < spin_max-1) && (old_sz_2 > 0)){
			new_num_choices -= 1;
		}
		if((new_sz_1 > 0) && (new_sz_2 < spin_max-1)){
			new_num_choices += 1;
		}
		if((new_sz_1 < spin_max-1) && (new_sz_2 > 0)){
			new_num_choices += 1;
		}

		return num_choices / new_num_choices;
	}
	
}


void mc_sz2(MCKPEPS &state, std::vector<double> &wavefunctions, std::vector<double> &values, int num_trials = 10000){
	Sz2 sz2op(state.Nx(), state.Ny(), state.physical_dims());
	mc_eval_single(state, sz2op, wavefunctions, values, num_trials);
}

void mc_eval_single(MCKPEPS &state, MCOperator &op, std::vector<double> &wavefunctions, std::vector<double> &values, int num_trials = 10000){
	std::vector<MCOperator> ops_wrapper{op};
	std::vector<std::vector<double>> values_wrapper;
	values_wrapper.push_back(std::vector<double>());
	mc_eval(state, ops_wrapper, wavefunctions, values_wrapper, num_trials);
	for(double val : values_wrapper[0]){
		values.push_back(val);
	}
}

//Computes average of a certain operator op
void mc_eval(MCKPEPS &state, std::vector<MCOperator> &ops, std::vector<double> &wavefunctions, std::vector<std::vector<double>> &values, int num_trials = 10000){
	int num_sites = state.size();
	std::vector<int> spin_config(num_sites, 0);
	std::mt19937 generator;
	std::uniform_real_distribution<double> distribution(0,1);
	randomize_in_sector(spin_config, state.physical_dims(), generator, distribution);
	double norm_estimate = 0;
	double old_wavefn = wavefunction(spin_config, state);
	wavefunctions.push_back(old_wavefn);
	for(int op_index = 0; op_index < ops.size(); op_index++){
		values[op_index].push_back(ops[op_index].eval(spin_config, spin_config));
	}
	//In each step:
	//Find new spin config by randomly flipping a few spins (say Lx pairs with 50% chance each)
	//Get wavefunction
	//Test if you will move to the new spin with metropolis probability
	for(int i = 0; i < num_trials; i++){
		std::vector<int> new_spin_config(spin_config);
		//randomize_in_sector(new_spin_config, state.physical_dims(), generator, distribution);
		double choice_ratio = flip_spins(state.bonds, new_spin_config, state.physical_dims(), generator, distribution);
		double new_wavefn = wavefunction(new_spin_config, state);
		bool switch_to_new_config = true;
		double acceptance_probability = choice_ratio*(new_wavefn*new_wavefn)/(old_wavefn*old_wavefn);
		if(acceptance_probability < distribution(generator)){
			switch_to_new_config = false;
		}
		std::cerr << "Spin config (";
		for(int spin : new_spin_config){
			std::cerr << spin << " ";
		}
		std::cerr << ") with wavefunction " << new_wavefn;
		if(switch_to_new_config){
			std::cerr << " accepted over old wavefunction ";
			old_wavefn = new_wavefn;
			spin_config = new_spin_config;
			wavefunctions.push_back(new_wavefn);
			for(int op_index = 0; op_index < ops.size(); op_index++){
				values[op_index].push_back(ops[op_index].eval(spin_config, spin_config));
			}
			//values.push_back(op.eval(spin_config, spin_config));
		}
		else{
			std::cerr << " rejected over old wavefunction ";
		}
		std::cerr << "(trial " << i << ")\n";
	}
}

#endif