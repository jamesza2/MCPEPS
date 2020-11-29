#ifndef GRAD_OPT
#define GRAD_OPT

#include "itensor/all.h"
#include "mcpeps.h"
#include "auxiliary_sampling.h"
#include "inner_sampling.h"
#include <ctime>
#include <iomanip>

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

itensor::ITensor density_matrix(MCKPEPS &PEPS1, int target_site){
	MCKPEPS PEPS2 = PEPS1;
	PEPS2.prime();
	PEPS2.site_tensor.setPrime(1);
	ArbitraryPEPS contracted = PEPS1.combine(PEPS2);

}

std::vector<itensor::ITensor> direct_gradient(MCKPEPS &PEPS1, const Heisenberg &H){
	PEPSop HPEPO = H.toPEPSop();
	MCKPEPS PEPS2 = PEPS1;
	PEPS2.prime();
	//std::cerr << "Inner product with copy...";
	double normsq = PEPS1.inner_product(PEPS2);
	int num_sites = PEPS1.size();
	std::vector<itensor::ITensor> direct_grad(num_sites);
	double energy = 0;
	//std::cerr << "Combining with copy...";
	ArbitraryPEPS contracted = PEPS1.combine(PEPS2);
	std::vector<itensor::ITensor> me1(num_sites); //Matrix elements of <psi|H|env>
	//std::cerr << "Computing <psi|env>...";

	
	std::vector<itensor::ITensor> me2 = contracted.environments(); //Matrix elements of <psi|env>

	for(int site = 0; site < num_sites; site++){
		me2[site] *= PEPS2.site_tensor(site);
	}
	int target_site = PEPS1.site_index_from_position(1,1,2);
	std::cerr << std::setprecision(3) << std::fixed;
	std::cerr << "Environment of (1,1,2):\nSample#0000: ";
	auto printElt = [](itensor::Real r){
		if(r >= 0){std::cerr << "+"; r+= 0.00001;}
		std::cerr << r << " ";};
	me2[target_site].visit(printElt);

	std::cerr << "\nSite tensor at (1,1,2):\nSample#0000: ";
	auto printElt = [](itensor::Real r){
		if(r >= 0){std::cerr << "+"; r+= 0.00001;}
		std::cerr << r << " ";};
	PEPS1.site_tensor(target_site).visit(printElt);

	//std::cerr << "Computing <psi|H|env...";
	for(Term t : HPEPO.terms){
		MCKPEPS PEPS_applied = PEPS2;
		t.apply(PEPS_applied);
		double energy_me = t.eval(PEPS1, PEPS2);
		double energy_part = energy_me/normsq;
		energy += energy_part;
		ArbitraryPEPS contracted_applied = PEPS1.combine(PEPS_applied);
		std::vector<itensor::ITensor> capp_envs = contracted_applied.environments();
		for(int site = 0; site < num_sites; site++){
			me1[site] += capp_envs[site]*PEPS_applied.site_tensor(site);
		}
	}
	//std::cerr << "Computing grads...";
	std::vector<itensor::ITensor> grads;
	

	for(int site = 0; site < num_sites; site++){
		itensor::ITensor gradient = (me1[site] - energy*me2[site])*2/normsq;
		if(site == target_site){
			std::cerr << std::setprecision(3) << std::fixed;
			std::cerr << "\nDirect gradient of (1,1,2):\nSample#0000: ";
			auto printElt = [](itensor::Real r){
				if(r >= 0){std::cerr << "+"; r+= 0.00001;}
				std::cerr << r << " ";};
			gradient.visit(printElt);
		}
		//PrintData(gradient);
		gradient = signelts(gradient);
		grads.push_back(gradient);
	}
	return grads;
}

//Gets samples of Delta, DeltaE, etc. for one step
void get_sample(MCKPEPS &psi, NoSitePEPS &contracted, std::vector<int> &spin_config, const Heisenberg &H, std::vector<itensor::ITensor> &Delta, std::vector<itensor::ITensor> &DeltaE, double &E, const double update_size, Randomizer &r){
	//std::cerr << "    Taking sample...";
	//Randomizer r;
	//std::cerr << "v direction...";
	sample_v_direction(psi, spin_config, r);
	//std::cerr << "s direction...";
	sample_s_direction(psi, spin_config, r);
	//std::cerr << "l direction...";
	double wavefn = sample_l_direction(psi, spin_config, r);
	//std::cerr << "finished sampling...";
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
	/*std::cerr << "Env0 norm: " << itensor::norm(envs[0]);
	std::cerr << " Product norm: " << itensor::elt(envs[0]*adapt_tensor(contracted, psi, 0));
	std::cerr << " Wavefunction: " << wavefn;
	std::cerr << " Real Wavefunction: " << real_wavefn;
	std::cerr << " Local Energy: " << local_energy << std::endl;*/

	for(int site_index = 0; site_index < psi.size(); site_index++){
		adapt_tensor(psi, contracted, envs[site_index], site_index);
		Delta[site_index] += envs[site_index]/wavefn;
		DeltaE[site_index] += envs[site_index]*local_energy/wavefn;
	}
	/*std::cerr << "E=" << local_energy << ", W=" << wavefn << std::endl;
	PrintData(envs[0]);*/
	E += local_energy;
}

//Updates the PEPS for one gradient optimization step. Returns the average energy.
double update(MCKPEPS &psi, 
	std::vector<int> &spin_config, 
	const Heisenberg &H, 
	const int M, 
	const double update_size, 
	Randomizer &r, 
	double &gradient_fidelity,
	const itensor::Args &optimize_args)
{
	//std::cerr << "Performing update..." << std::endl;
	std::vector<itensor::ITensor> Delta(psi.size());
	std::vector<itensor::ITensor> DeltaE(psi.size());
	double E = 0;
	SpinConfigPEPS scp(psi, spin_config, 1);
	NoSitePEPS nsp = psi.contract(scp);

	int target_site = psi.site_index_from_position(1,1,2);
	std::vector<itensor::ITensor> direct_grads = direct_gradient(psi, H);
	
	std::cerr << "\nCurrent sampled gradient of (1,1,2):\n\r";
	for(int eq_step = 0; eq_step < 10; eq_step++){
		std::vector<itensor::ITensor> dud_Delta(psi.size());
		std::vector<itensor::ITensor> dud_DeltaE(psi.size());
		double dudE;
		get_sample(psi, nsp, spin_config, H, dud_Delta, dud_DeltaE, dudE, update_size, r);
	}

	for(int sample = 0; sample < M; sample++){
		//std::cerr << "Getting sample #" << sample+1 << "...";
		//std::cerr << "Update step #" << sample+1 << std::endl;
		get_sample(psi, nsp, spin_config, H, Delta, DeltaE, E, update_size, r);
		itensor::ITensor current_grad = DeltaE.at(target_site)*2./(sample+1) - Delta.at(target_site)*E*2./((sample+1)*(sample+1));
		//itensor::ITensor current_grad = Delta.at(target_site)*E*2./((sample+1)*(sample+1));
		std::cerr << "Sample#" << sample+1 << ": ";
		current_grad.visit(printElt);
		std::cerr << "\r";
		
	}

	std::cerr << std::endl;
	E /= M;
	itensor::ITensor grad;
	double grads_factor = 2./M;
	std::vector<itensor::ITensor> grads;
	for(int site = 0; site < Delta.size(); site++){
		//std::cerr << "Assembling gradient#" << site+1 << "...";
		grad = DeltaE.at(site)*grads_factor - Delta.at(site)*grads_factor*E;
		//grad /= norm(grad);
		grad = signelts(grad);
		grads.push_back(grad);
	}
	std::cerr << std::setprecision(6) << std::defaultfloat;

	for(int site = 0; site < Delta.size(); site++){
		auto [i,j,k] = psi.position_of_site(site);
		psi._site_tensors[i][j][k] -= update_size*r.rand()*grads[site];

		double fidelity = itensor::norm(direct_grads[site]*grads[site])/itensor::norm(grads[site]*grads[site]);
		gradient_fidelity += fidelity;
		psi._site_tensors[i][j][k] /= (itensor::norm(psi._site_tensors[i][j][k])*WAVEFUNCTION_NORMALIZATION_CONSTANT);
	}
	gradient_fidelity /= Delta.size();

	randomize_in_sector(spin_config, psi.physical_dims(), r.gen, r.dist, 0);
	return E;
}

//Optimizes psi, storing the vector in energies
void optimize(MCKPEPS &psi, 
	std::vector<double> &energies, 
	std::vector<double> &update_sizes, 
	const std::map<std::string, double> &Jvals, 
	const int M, 
	const int opt_steps, 
	std::vector<double> &fidelities, 
	const itensor::Args &optimize_args)
{
	std::vector<int> spin_config(psi.size(), 0);
	Randomizer r;
	randomize_in_sector(spin_config, psi.physical_dims(), r.gen, r.dist);
	Heisenberg H(psi.Nx(), psi.Ny(), psi.physical_dims());
	H.set_J(Jvals);
	double update_size = 0.05;
	//std::cerr << "Initializing updates...";
	if(update_sizes.size() > 0){update_size = update_sizes.at(0);}
	else{update_sizes.push_back(0.05);}
	//double update_size = update_size_init;
	auto timestart = std::time(NULL);
	for(int step = 0; step < opt_steps; step ++){
		//std::cerr << "starting step #" << step+1;
		double fidelity;
		double energy = update(psi, spin_config, H, M, update_size, r, fidelity, optimize_args);
		fidelities.push_back(fidelity);
		energies.push_back(energy);
		std::cerr << "STEP#" << step+1 << " HAS ENERGY " << energy << " (" << std::difftime(std::time(NULL), timestart) << "s)" << std::endl;
		if(update_sizes.size() > step+1){update_size = update_sizes.at(step+1);}
		else{
			update_size *= 0.97;
			if(step != opt_steps-1){update_sizes.push_back(update_size);}
		}
		timestart = std::time(NULL);
	}
}


#endif