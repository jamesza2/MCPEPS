#ifndef GRAD_OPT
#define GRAD_OPT

#include "itensor/all.h"
#include "mcpeps.h"
#include "auxiliary_sampling.h"
#include "inner_sampling.h"
#include <ctime>

itensor::ITensor signelts(itensor::ITensor site){
	auto signelt = [](itensor::Real r) {
		itensor::Real ret = 0;
		if(r < 0){
			ret -= 1;
		}
		if(r > 0){
			ret += 1;
		}
		return ret;
	};
	site.apply(signelt);
	return site;
}

//Gets samples of Delta, DeltaE, etc. for one step
void get_sample(MCKPEPS &psi, NoSitePEPS &contracted, std::vector<int> &spin_config, const Heisenberg &H, std::vector<itensor::ITensor> &Delta, std::vector<itensor::ITensor> &DeltaE, double &E, const double update_size){
	//std::cerr << "    Taking sample...";
	Randomizer r;
	//std::cerr << "v direction...";
	sample_v_direction(psi, spin_config, r);
	//std::cerr << "s direction...";
	sample_s_direction(psi, spin_config, r);
	//std::cerr << "l direction...";
	double wavefn = sample_l_direction(psi, spin_config, r);
	double real_wavefn = wavefunction(spin_config, psi);
	auto possible_mes = H.possible_matrix_elements(spin_config);
	double local_energy = 0;
	for(int me_index = 0; me_index < possible_mes.size(); me_index++){
		double new_wavefn = wavefunction(possible_mes[me_index].first, psi);
		local_energy += new_wavefn*possible_mes[me_index].second/wavefn;
	}
	//std::cerr << "Finding environments..." << std::endl;
	SpinConfigPEPS scp(psi, spin_config, 1);
	contracted = psi.contract(scp);
	auto envs = contracted.environments(psi.site_indices, spin_config);
	std::cerr << "Env0 norm: " << itensor::norm(envs[0]);
	std::cerr << " Product norm: " << itensor::elt(envs[0]*psi._site_tensors[0][0][0]);
	std::cerr << " Wavefunction: " << wavefn;
	std::cerr << " Real Wavefunction: " << real_wavefn;
	std::cerr << " Local Energy: " << local_energy << std::endl;

	for(int site_index = 0; site_index < psi.size(); site_index++){
		adapt_tensor(psi, contracted, envs[site_index], site_index);
		Delta[site_index] += envs[site_index]/wavefn;
		DeltaE[site_index] += envs[site_index]*local_energy/wavefn;
	}
	E += local_energy;
}

//Updates the PEPS for one gradient optimization step. Returns the average energy.
double update(MCKPEPS &psi, std::vector<int> &spin_config, const Heisenberg &H, const int M, const double update_size, Randomizer &r){
	//std::cerr << "Performing update..." << std::endl;
	std::vector<itensor::ITensor> Delta(psi.size());
	std::vector<itensor::ITensor> DeltaE(psi.size());
	double E = 0;
	SpinConfigPEPS scp(psi, spin_config, 1);
	NoSitePEPS nsp = psi.contract(scp);
	for(int sample = 0; sample < M; sample++){
		get_sample(psi, nsp, spin_config, H, Delta, DeltaE, E, update_size);
	}
	E /= M;
	itensor::ITensor grad;
	double grads_factor = 2./M;
	for(int site = 0; site < Delta.size(); site++){
		grad = DeltaE[site]*grads_factor - Delta[site]*grads_factor*E;
		//grad /= norm(grad);
		grad = signelts(grad);
		auto [i,j,k] = psi.position_of_site(site);
		psi._site_tensors[i][j][k] -= update_size*r.rand()*grad;
	}
	return E;
}

void optimize(MCKPEPS &psi, std::vector<double> &energies, const std::map<std::string, double> &Jvals, const int M, const double update_size_init = 0.05, const int opt_steps = 100){
	std::vector<int> spin_config(psi.size(), 0);
	Randomizer r;
	randomize_in_sector(spin_config, psi.physical_dims(), r.gen, r.dist);
	Heisenberg H(psi.Nx(), psi.Ny(), psi.physical_dims());
	H.set_J(Jvals);
	double update_size = update_size_init;
	auto timestart = std::time(NULL);
	for(int step = 0; step < opt_steps; step ++){
		double energy = update(psi, spin_config, H, M, update_size, r);
		energies.push_back(energy);
		std::cerr << "STEP#" << step+1 << " HAS ENERGY " << energy << " (" << std::difftime(std::time(NULL), timestart) << "s)" << std::endl;
		update_size *= 0.97;
		timestart = std::time(NULL);
	}
}


#endif