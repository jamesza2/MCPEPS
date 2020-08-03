#ifndef mccalculator
#define mccalculator

#include "mcpeps.h"
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

//Takes random groups of 2 spins (not necessarily pairs) and redistributes their spins such that the sum of their sz's is the same
//Repeats num_spins_to_flip times
void flip_spins(std::vector<int> &spin_config, int spin_max, std::mt19937 &generator, std::uniform_real_distribution<double> &distribution, int num_spins_to_flip){
	for(int pairs_found = 0; pairs_found < num_spins_to_flip; pairs_found++){
		if(distribution(generator) >= 0.5){
			continue;
		}
		int site_choice_1 = std::floor(spin_max*distribution(generator));
		int site_choice_2 = std::floor(spin_max*distribution(generator));
		if(spin_config[site_choice_1] == spin_config[site_choice_2]){
			pairs_found -= 1;
			continue;
		}
		int total_site = spin_config[site_choice_1]+spin_config[site_choice_2];
		int new_sz_1 = std::floor((total_site+1)*distribution(generator));
		while((new_sz_1 >= spin_max) || (total_site - new_sz_1 >= spin_max)){
			new_sz_1 = std::floor((total_site+1)*distribution(generator));
		}
		int new_sz_2 = total_site - new_sz_1;
		spin_config[site_choice_1] = new_sz_1;
		spin_config[site_choice_2] = new_sz_2;
	}
	
}

void mc_norm(MCKPEPS &state, std::vector<double> &wavefunctions, int num_trials = 10000, int num_spins_to_flip = -1){
	int num_sites = state.size();
	std::vector<int> &spin_config(num_sites, 0);
	if(num_spins_to_flip == -1){
		num_spins_to_flip = num_sites;
	}
	std::mt19937 generator;
	std::uniform_real_distribution<double> distribution(0,1);
	randomize_in_sector(spin_config, state.physical_dims(), generator, distribution);
	double norm_estimate = 0;
	double old_wavefn = wavefunction(spin_config, state);
	wavefunctions.push_back(old_wavefn);
	//In each step:
	//Find new spin config by randomly flipping a few spins (say Lx pairs with 50% chance each)
	//Get wavefunction
	//Test if you will move to the new spin with metropolis probability
	for(int i = 0; i < num_trials; i++){
		std::vector<int> new_spin_config(spin_config);
		flip_spins(new_spin_config, state.physical_dims(), generator, distribution,num_spins_to_flip);
		double new_wavefn = wavefunction(new_spin_config, state);
		bool switch_to_new_config = true;
		if(std::abs(new_wavefn) < std::abs(old_wavefn)){
			if((new_wavefn*new_wavefn)/(old_wavefn*old_wavefn) < distribution(generator)){
				switch_to_new_config = false;
			}
		}
		std::cerr << "Spin config (";
		for(int spin : spin_config){
			std::cerr << spin << " ";
		}
		std::cerr << ") with wavefunction " << new_wavefn;
		if(switch_to_new_config){
			std::cerr << " accepted over old wavefunction ";
			old_wavefn = new_wavefn;
			spin_config = new_spin_config;
			wavefunctions.push_back(new_wavefn);
		}
		else{
			std::cerr << " rejected over old wavefunction ";
		}
		std::cerr << "(trial " << i << ")\n";
	}

}

#endif