#ifndef AUX_SAMPLING
#define AUX_SAMPLING

#include "itensor/all.h"
#include "mcpeps.h"
#include "auxiliary_mps.h"

class Randomizer{
	public:
		std::mt19937 gen;
		std::uniform_real_distribution<double> dist;
		Randomizer(){
			distribution = std::uniform_real_distribution<double> (0,1);
		}
		Randomizer(std::mt19937 &generator, std::uniform_real_distribution<double> &distribution){
			gen = generator;
			dist = distribution;
		}
		double rand(){
			return distribution(generator);
		}
		double rand(double min, double max){
			return min + distribution(generator)*(max-min);
		}
};

void update_site_tensor(NoSitePEPS &no_site, MCKPEPS &original, int i, int j, int k, int new_sz){
	tensor::Index site_at = original._site_indices[original.site_index_from_position(i,j,k)];
	new_tensor = adapt_tensor(no_site, original, i, j, k) * itensor::setElt(site_at = new_sz+1);
	no_site._site_tensors[i][j][k] = new_tensor;
}

//adapt a original_PEPS tensor at (i,j,k) to use the link indices of target_PEPS instead
itensor::ITensor adapt_tensor(NoSitePEPS &target_PEPS, NoSitePEPS &original_PEPS, int i, int j, int k){
	int site_index = target_PEPS.site_index_from_position(i,j,k);
	std::vector<std::pair<itensor::Index, itensor::Index>> links;
	for(int neighbor_index : target_PEPS.bonds.nn_at(site_index)){
		itensor::Index target_link = target_PEPS._link_indices[pair_to_link_index(site_index, neighbor_index)];
		itenosr::Index original_link = original_PEPS._link_indices[pair_to_link_index(site_index, neighbor_index)];
		links.push_back(std::pair<itensor::Index, itensor::Index>>(target_link, original_link));
	}
	itensor::ITensor new_tensor = original_PEPS._site_tensors[i][j][k];
	for(int l = 0; l < links.size(); l++){
		new_tensor *= itensor::delta(links[l].first, links[l].second);
	}
	//Check new tensor and no site tensor have the same indices
	if(!itensor::hasSameInds(target_PEPS._site_tensors[i][j][k].inds(), new_tensor.inds())){
		std::cout << "WARNING: Index sets of original no-site tensor and new tensor different" << std::endl;
		Print(target_PEPS._site_tensors[i][j][k]);
		Print(new_tensor);
	}
	return new_tensor;
}

double test_bond(NoSitePEPS &no_site, MCKPEPS &original, std::vector<int> &spin_config,
		int i1, int j1, int k1, int i2, int j2, int k2,
		itensor::ITensor &l_aux,
		itensor::ITensor &r_aux,
		itensor::ITensor &u_aux_1,
		itensor::ITensor &u_aux_2,
		itensor::ITensor &d_aux_1,
		itensor::ITensor &d_aux_2,
		double old_wavefunction,
		Randomizer &r){
	int old_sz_1 = spin_config[original.site_index_from_position(i1,j1,k1)];
	int old_sz_2 = spin_config[original.site_index_from_position(i2,j2,k2)];
	int up_or_down = std::floor(2*r.rand())*2-1; 
	int new_sz_1 = old_sz_1 + up_or_down;
	int new_sz_2 = old_sz_2 - up_or_down;
	double wavefn_to_return = old_wavefunction;
	if((new_sz_1 >= 0) && (new_sz_2 >= 0) && (new_sz_1 < spin_max) && (new_sz_2 < spin_max)){
		/*Contract the tensor group
				VU(iJ)---VU(iJ+1)
	VL(iJ-1)----Psi(iJ)--Psi(iJ+1)---VR(iJ+2)
				VD(iJ-1)-VD(iJ)
		with the Psi's set to new_sz_1, new_sz_2 using setElt*/
		itensor::ITensor new_product_1 = l_aux;
		new_product_1 *= u_aux_1;
		itensor::Index site_at = original._site_indices[original.site_index_from_position(i1,j1,k1)];
		new_product_1 *= (adapt_tensor(no_site, original, i1, j1, k1)*itensor::setElt(site_at = new_sz_1+1));
		new_product_1 *= d_aux_1;
		itensor::ITensor new_product_2 = r_aux;
		new_product_2 *= u_aux_2;
		site_at = psi_sites._site_indices[original.site_index_from_position(i2,j2,k2)];
		new_product_2 *= (adapt_tensor(no_site, original, i2, j2, k2)*itensor::setElt(site_at = new_sz_2+1));
		new_product_2 *= d_aux_2;
		itensor::ITensor total_product = new_product_1*new_product_2;
		double new_wavefunction = itensor::norm(total_product);
		/* Shouldn't have to worry about number of choices in the sequential case, because selection probability is always 1/2
		if((old_sz_1 > 0) && (old_sz_2 < spin_max-1)){ new_num_choices -= 1; }
		if((old_sz_1 < spin_max-1) && (old_sz_2 > 0)){ new_num_choices -= 1; }
		if((new_sz_1 > 0) && (new_sz_2 < spin_max-1)){ new_num_choices += 1; }
		if((new_sz_1 < spin_max-1) && (new_sz_2 > 0)){ new_num_choices += 1; }*/
		double transition_probability = new_wavefunction*new_wavefunction/(old_wavefunction*old_wavefunction);
		if(r.rand() < transition_probability){ //Transition is a success, change to new spin config
			//config.set_spin(i,j,1,new_sz_1, WAVEFUNCTION_NORMALIZATION_CONSTANT);
			//config.set_spin(i,j,2,new_sz_2, WAVEFUNCTION_NORMALIZATION_CONSTANT);
			update_site_tensor(no_site, original, i1,j1,k1, new_sz_1);
			update_site_tensor(no_site, original, i2,j2,k2, new_sz_1);
			spin_config[no_site.site_index_from_position(i1,j1,k1)] = new_sz_1;
			spin_config[no_site.site_index_from_position(i2,j2,k2)] = new_sz_2;
			std::cout << "Spin config ";
			for(int sp : spin_config){std::cout << sp << " ";}
			std::cout << "Accepted with wavefunction " << new_wavefunction;
			wavefn_to_return = new_wavefunction;
		}
		else{
			std::cout << "Spin config ";
			for(int sp : spin_config){std::cout << sp << " ";}
			std::cout << "Rejected with wavefunction " << new_wavefunction;
		}
	}
	//No matter if the move is accepted or rejected, we should update the auxiliary tensor accordingly
	l_aux *= u_aux_1;
	l_aux *= no_site._site_tensors[i1][j1][k1];
	l_aux *= d_aux_1;
	return wavefn_to_return;
}

double sample_v_direction(MCKPEPS &psi_sites, std::vector<int> &spin_config, Randomizer &r){
	SpinConfigPEPS config(psi_sites, spin_config, WAVEFUNCTION_NORMALIZATION_CONSTANT);
	NoSitePEPS psi = psi_sites.contract(config);
	std::list<AuxMPS> vd_list = psi.get_vd_auxiliaries();
	std::cout << "Creating up auxiliary..." << std::endl;
	double most_recent_wavefunction = -1;
	AuxMPS VUi(AuxType::VU);//Create the up auxiliary (won't have some links between certain sites though)
	for(int j = 0; j < psi._Ny; j++){
		auto left_links = itensor::commonInds(psi._site_tensors[0][j][0], psi._site_tensors[0][j][1]);
		auto [left_tensor, sing_vals, right_tensor] = itensor::svd(psi._site_tensors[0][j][0], left_links, {"MaxDim", _Dc});
		VUi.add_tensor(left_tensor);
		VUi.add_tensor(sing_vals*right_tensor);
	}
	//Begin sweeping along the rows
	auto vd_it = vd_list.begin();


	for(int i = 0; i < psi._Nx; i++){
		std::cout << "Sweeping row " << i << std::endl;
		vd_it++;
		//Create right auxiliary tensors
		std::vector<itensor::ITensor> vr_auxiliaries(psi._Ny*2+1);
		vr_auxiliaries[2*_Ny] = itensor::ITensor(1);
		for(int j = _Ny-1; j >= 0; j--){
			for(int k = 2; k >= 1; k--){
				int J = 2*j+k-1; //J is the index of the VU MPS tensors and left/right auxiliary tensors, J-1 is the index of the VD tensors
				vr_auxiliaries[J] = (vr_auxiliaries[J+1]*VUi.MPS[J])*psi._site_tensors[i][j][k];
				if(J>0){
					vr_auxiliaries[J] *= vd_it->MPS[J-1];
				}
			}
		}
		//Evaluate the old wavefunction
		double old_wavefunction = itensor::norm(vr_auxiliairies[0]);
		//double old_num_choices = config.num_choices();
		itensor::ITensor vl_auxiliary(1);
		//Sweep the horizontal (1-2) links from left to right
		for(int j = 0; j < _Ny; j++){

			int J = 2*j;
			itensor::ITensor vd_aux_1(1);
			if(j > 0){vd_aux_1 = vd_it->MPS[J-1];}
			most_recent_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,1,i,j,2, vl_auxiliary, vr_auxiliaries[J+2],VUi.MPS[J],VUi.MPS[J+1],vd_aux_1, vd_it->MPS[J], old_wavefunction, r);
			if(j < _Ny-1){
				most_recent_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,2,i,j+1,1, vl_auxiliary, vr_auxiliaries[J+3],VUi.MPS[J+1],VUi.MPS[J+2],vd_it->MPS[J], vd_it->MPS[J+1], old_wavefunction, r);
			}
			/*
			int old_sz_1 = spin_config[psi.site_index_from_position(i,j,1)];
			int old_sz_2 = spin_config[psi.site_index_from_position(i,j,2)];
			int up_or_down = std::floor(2*r.rand())*2-1; 
			int new_sz_1 = old_sz_1 + up_or_down;
			int new_sz_2 = old_sz_2 - up_or_down;
			if((new_sz_1 >= 0) && (new_sz_2 >= 0) && (new_sz_1 < spin_max) && (new_sz_2 < spin_max)){
				int J = 2*j;
				itensor::ITensor new_product_1 = vl_auxiliary;
				new_product_1 *= VUi.MPS[J];
				itensor::Index site_at = psi_sites._site_indices[psi_sites.site_index_from_position(i,j,1)];
				new_product_1 *= (psi_sites._site_tensors[i][j][1]*itensor::setElt(site_at = new_sz_1+1));
				new_product_1 *= vd_it->MPS[J-1];
				itensor::ITensor new_product_2 = vr_auxiliaries[J+2];
				new_product_2 *= VUi.MPS[J+1];
				site_at = psi_sites._site_indices[psi_sites.site_index_from_position(i,j,2)];
				new_product_2 *= (psi_sites._site_tensors[i][j][2]*itensor::setElt(site_at = new_sz_2+1));
				new_product_2 *= vd_it->MPS[J];
				itensor::ITensor total_product = new_product_1*new_product_2;
				double new_wavefunction = itensor::norm(total_product);
				double transition_probability = new_wavefunction*new_wavefunction/(old_wavefunction*old_wavefunction);
				if(r.rand() < transition_probability){ //Transition is a success, change to new spin config
					//config.set_spin(i,j,1,new_sz_1, WAVEFUNCTION_NORMALIZATION_CONSTANT);
					//config.set_spin(i,j,2,new_sz_2, WAVEFUNCTION_NORMALIZATION_CONSTANT);
					update_site_tensor(psi, psi_sites, i,j,1,new_sz_1);
					update_site_tensor(psi, psi_sites, i,j,2,new_sz_1);
				}
			}*/
		}
		if(i < psi._Nx-1){
			//Create new VU auxiliary MPS
			std::vector<itensor::ITensor> row_i_contracted;
			for(int j = 0; j < _Ny; j++){
				int J = 2*j;
				row_i_contracted.push_back(VUi[J]*psi._site_tensors[i][j][1]);
				row_i_contracted.push_back(VUi[J+1]*psi._site_tensors[i][j][2]);
			}
			//Contract the bond dimensions of the contracted row i
			for(int J = 0; J < 2*_Ny-2; J++){
				auto forward_links = itensor::commonInds(row_i_contracted[J], row_i_contracted[J+1]);
				auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_i_contracted[J], forward_links, {"MaxDim", psi._Dc});
				row_i_contracted[J+1] *= (forward_tensor*sing_vals);
				row_i_contracted[J] = back_tensor;
			}
			std::vector<itensor::ITensor>row_ip_unsplit;
			//Apply the tensors to the next row down, then split them
			for(int j = 0; j < _Ny; j++){
				int J = 2*j;
				if(j > 0){
					row_ip_unsplit[j-1] *= row_i_contracted[J];
				}
				row_ip_unsplit.push_back(row_i_contracted[J+1]*psi._site_tensors[i+1][j][0]);
			}
			for(int j = 0; j < _Ny; j++){
				int J = 2*j;
				itensor::IndexSet forward_inds;
				itensor::Index always_forward = itensor::commonIndex(row_ip_unsplit[j], psi._site_tensors[i+1][j][2]);
				if(j == _Ny-1){forward_inds = itensor::IndexSet(always_forward);}
				else{forward_inds = itensor::IndexSet(always_forward, itensor::commonIndex(row_ip_unsplit[j], row_ip_unsplit[j+1]));}
				auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_ip_unsplit[j], forward_inds, {"MaxDim", psi._Dc});
				VUi.MPS[J] = back_tensor;
				VUi.MPS[J+1] = forward_tensor*sing_vals;
			}
		}
	}	
	return most_recent_wavefunction;
}

double sample_s_direction(MCKPEPS &psi_sites, std::vector<int> &spin_config, Randomizer &r){
	SpinConfigPEPS config(psi_sites, spin_config, WAVEFUNCTION_NORMALIZATION_CONSTANT);
	NoSitePEPS psi = psi_sites.contract(config);
	std::list<AuxMPS> sd_list = psi.get_sd_auxiliaries();
	most_recent_wavefunction = -1;
	AuxMPS SUi(AuxType::SU);//Create the up auxiliary (won't have some links between certain sites though)
	for(int i = 0; i < psi._Nx; i++){
		auto left_links = itensor::commonInds(psi._site_tensors[i][0][1], psi._site_tensors[i][0][0]);
		auto [left_tensor, sing_vals, right_tensor] = itensor::svd(psi._site_tensors[i][0][1], left_links, {"MaxDim", _Dc});
		SUi.add_tensor(left_tensor);
		SUi.add_tensor(sing_vals*right_tensor);
	}
	//Begin sweeping along the rows
	auto sd_it = sd_list.begin();

	//double old_wavefunction = -1;
	//int old_num_choices = -1;
	for(int j = 0; j < psi._Ny; j++){
		sd_it++;
		//Create right auxiliary tensors
		std::vector<itensor::ITensor> sr_auxiliaries(psi._Nx*2+1);
		sr_auxiliaries[2*_Nx] = itensor::ITensor(1);
		for(int i = _Nx-1; i >= 0; i--){
			for(int k = 2; k >= 0; k-=2){
				int I = 2*i+(k/2)-1;//I is the index of the VU MPS tensors and left/right auxiliary tensors, I-1 is the index of the VD tensors
				sr_auxiliaries[I] = (sr_auxiliaries[I+1]*SUi.MPS[I])*psi._site_tensors[i][j][k];
				if(I > 0){
					sr_auxiliaries[I] *= sd_it->MPS[I-1];
				}
			}
		}
		//Evaluate the old wavefunction
		double old_wavefunction = itensor::norm(sr_auxiliairies[0]);
		//double old_num_choices = config.num_choices();
		itensor::ITensor sl_auxiliary(1);
		//Sweep the short direction (0-2) links from left to right
		for(int i = 0; i < _Nx; i++){
			int I = 2*i;
			itensor::ITensor sd_aux_1(1);
			if(i > 0){sd_aux_1 = sd_it->MPS[I-1];}
			most_recent_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,0,i,j,2, sl_auxiliary, sr_auxiliaries[I+2], SUi.MPS[I],SUi.MPS[I+1],sd_aux_1, sd_it->MPS[I], old_wavefunction, r);
			if(i < _Nx-1){
				most_recent_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,2,i+1,j,0, sl_auxiliary, sr_auxiliaries[I+3],SUi.MPS[I+1],SUi.MPS[I+2],sd_it->MPS[I], sd_it->MPS[I+1], old_wavefunction, r);
			}
		}
		if(j < psi._Ny-1){
			//Create new SU auxiliary MPS
			std::vector<itensor::ITensor> row_j_contracted;
			for(int i = 0; i < _Nx; i++){
				int I = 2*i;
				row_j_contracted.push_back(SUi[I]*psi._site_tensors[i][j][0]);
				row_j_contracted.push_back(SUi[I+1]*psi._site_tensors[i][j][2]);
			}
			//Contract the bond dimensions of the contracted row i
			for(int I = 0; I < 2*_Nx-2; I++){
				auto forward_links = itensor::commonInds(row_j_contracted[I], row_j_contracted[I+1]);
				auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_j_contracted[I], forward_links, {"MaxDim", psi._Dc});
				row_j_contracted[I+1] *= (forward_tensor*sing_vals);
				row_j_contracted[I] = back_tensor;
			}
			std::vector<itensor::ITensor>row_jp_unsplit;
			//Apply the tensors to the next row down, then split them
			for(int i = 0; i < _Nx; i++){
				int I = 2*i;
				if(i > 0){
					row_jp_unsplit[i-1] *= row_j_contracted[I];
				}
				row_jp_unsplit.push_back(row_j_contracted[I+1]*psi._site_tensors[i][j+1][1]);
			}
			for(int i = 0; i < _Ny; i++){
				int I = 2*i;
				itensor::IndexSet forward_inds;
				itensor::Index always_forward = itensor::commonIndex(row_jp_unsplit[i], psi._site_tensors[i][j+1][2]);
				if(i == _Nx-1){forward_inds = itensor::IndexSet(always_forward);}
				else{forward_inds = itensor::IndexSet(always_forward, itensor::commonIndex(row_jp_unsplit[i], row_jp_unsplit[i+1]));}
				auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_jp_unsplit[i], forward_inds, {"MaxDim", psi._Dc});
				SUi.MPS[I] = back_tensor;
				SUi.MPS[I+1] = forward_tensor*sing_vals;
			}
		}
	}	
	return most_recent_wavefunction;
}

#endif