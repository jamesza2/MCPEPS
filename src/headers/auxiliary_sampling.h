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
			dist = std::uniform_real_distribution<double> (0,1);
		}
		Randomizer(std::mt19937 &generator, std::uniform_real_distribution<double> &distribution){
			gen = generator;
			dist = distribution;
		}
		double rand(){
			return dist(gen);
		}
		double rand(double min, double max){
			return min + dist(gen)*(max-min);
		}
};



//adapt a original_PEPS tensor at (i,j,k) to use the link indices of target_PEPS instead
itensor::ITensor adapt_tensor(NoSitePEPS &target_PEPS, NoSitePEPS &original_PEPS, int i, int j, int k){
	int site_index = target_PEPS.site_index_from_position(i,j,k);
	std::vector<std::pair<itensor::Index, itensor::Index>> links;
	for(int neighbor_index : target_PEPS.bonds.nn_at(site_index)){
		itensor::Index target_link = target_PEPS.link_index(site_index, neighbor_index);
		itensor::Index original_link = original_PEPS.link_index(site_index, neighbor_index);
		links.push_back(std::pair<itensor::Index, itensor::Index>(target_link, original_link));
	}
	itensor::ITensor new_tensor = original_PEPS._site_tensors[i][j][k];
	for(int l = 0; l < links.size(); l++){
		new_tensor *= itensor::delta(links[l].first, links[l].second);
	}
	//Check new tensor contains the same indices as no site tensor except for the extra site index
	/*if(itensor::length(itensor::uniqueInds(new_tensor.inds(), target_PEPS._site_tensors[i][j][k].inds())) != 1){
		std::cout << "WARNING: Index sets of original no-site tensor and new tensor different" << std::endl;
		Print(target_PEPS._site_tensors[i][j][k]);
		Print(new_tensor);
	}*/
	return new_tensor;
}
itensor::ITensor adapt_tensor(NoSitePEPS &target_PEPS, NoSitePEPS &original_PEPS, int site){
	auto [i,j,k] = target_PEPS.position_of_site(site);
	return adapt_tensor(target_PEPS, original_PEPS, i,j,k);
}

void update_site_tensor(NoSitePEPS &no_site, MCKPEPS &original, int i, int j, int k, int new_sz, double wavefunction_normalization = 1){
	itensor::Index site_at = original.site_indices[original.site_index_from_position(i,j,k)];
	itensor::ITensor new_tensor = adapt_tensor(no_site, original, i, j, k) * itensor::setElt(site_at = new_sz+1) / wavefunction_normalization;
	no_site._site_tensors[i][j][k] = new_tensor;
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
		Randomizer &r,
		double wavefunction_normalization = 1){
	int old_sz_1 = spin_config[original.site_index_from_position(i1,j1,k1)];
	int old_sz_2 = spin_config[original.site_index_from_position(i2,j2,k2)];
	int up_or_down = std::floor(2*r.rand())*2-1; 
	int new_sz_1 = old_sz_1 + up_or_down;
	int new_sz_2 = old_sz_2 - up_or_down;
	double wavefn_to_return = old_wavefunction;
	int spin_max = original.physical_dims();

	/*
	itensor::ITensor old_product_1 = l_aux;
	old_product_1 *= u_aux_1;
	itensor::Index old_site_at = original.site_indices[original.site_index_from_position(i1,j1,k1)];
	old_product_1 *= (adapt_tensor(no_site, original, i1, j1, k1)*itensor::setElt(old_site_at = old_sz_1+1)/wavefunction_normalization);
	old_product_1 *= d_aux_1;
	itensor::ITensor old_product_2 = r_aux;
	old_product_2 *= u_aux_2;
	old_site_at = original.site_indices[original.site_index_from_position(i2,j2,k2)];
	old_product_2 *= (adapt_tensor(no_site, original, i2, j2, k2)*itensor::setElt(old_site_at = old_sz_2+1)/wavefunction_normalization);
	old_product_2 *= d_aux_2;
	itensor::ITensor old_product = old_product_1*old_product_2;
	//Print(old_product);*/

	if((new_sz_1 >= 0) && (new_sz_2 >= 0) && (new_sz_1 < spin_max) && (new_sz_2 < spin_max)){
		/*Contract the tensor group
				VU(iJ)---VU(iJ+1)
	VL(iJ-1)----Psi(iJ)--Psi(iJ+1)---VR(iJ+2)
				VD(iJ-1)-VD(iJ)
		with the Psi's set to new_sz_1, new_sz_2 using setElt*/
		itensor::ITensor new_product_1 = l_aux;
		new_product_1 *= u_aux_1;
		itensor::Index site_at = original.site_indices[original.site_index_from_position(i1,j1,k1)];
		new_product_1 *= (adapt_tensor(no_site, original, i1, j1, k1)*itensor::setElt(site_at = new_sz_1+1)/wavefunction_normalization);
		new_product_1 *= d_aux_1;
		itensor::ITensor new_product_2 = r_aux;
		new_product_2 *= u_aux_2;
		site_at = original.site_indices[original.site_index_from_position(i2,j2,k2)];
		new_product_2 *= (adapt_tensor(no_site, original, i2, j2, k2)*itensor::setElt(site_at = new_sz_2+1)/wavefunction_normalization);
		new_product_2 *= d_aux_2;
		itensor::ITensor total_product = new_product_1*new_product_2;
		//Print(total_product);
		double new_wavefunction = itensor::norm(total_product);
		/* Shouldn't have to worry about number of choices in the sequential case, because selection probability is always 1/2
		if((old_sz_1 > 0) && (old_sz_2 < spin_max-1)){ new_num_choices -= 1; }
		if((old_sz_1 < spin_max-1) && (old_sz_2 > 0)){ new_num_choices -= 1; }
		if((new_sz_1 > 0) && (new_sz_2 < spin_max-1)){ new_num_choices += 1; }
		if((new_sz_1 < spin_max-1) && (new_sz_2 > 0)){ new_num_choices += 1; }*/
		double transition_probability = new_wavefunction*new_wavefunction/(old_wavefunction*old_wavefunction);
		double choice = r.rand();
		if(choice < transition_probability){ //Transition is a success, change to new spin config
			//config.set_spin(i,j,1,new_sz_1, WAVEFUNCTION_NORMALIZATION_CONSTANT);
			//config.set_spin(i,j,2,new_sz_2, WAVEFUNCTION_NORMALIZATION_CONSTANT);
			update_site_tensor(no_site, original, i1,j1,k1, new_sz_1, wavefunction_normalization);
			update_site_tensor(no_site, original, i2,j2,k2, new_sz_2, wavefunction_normalization);
			spin_config[no_site.site_index_from_position(i1,j1,k1)] = new_sz_1;
			spin_config[no_site.site_index_from_position(i2,j2,k2)] = new_sz_2;
			std::cout << "Spin config ";
			for(int sp : spin_config){std::cout << sp << " ";}
			std::cout << "Accepted with wavefunction " << new_wavefunction;
			std::cout << " (Transition probability " << transition_probability << ",";
			std::cout << " Rolled " << choice << ")" << std::endl;
			wavefn_to_return = new_wavefunction;
		}
		else{
			spin_config[no_site.site_index_from_position(i1,j1,k1)] = new_sz_1;
			spin_config[no_site.site_index_from_position(i2,j2,k2)] = new_sz_2;
			std::cout << "Spin config ";
			for(int sp : spin_config){std::cout << sp << " ";}
			std::cout << "Rejected with wavefunction " << new_wavefunction;
			std::cout << " (Transition probability " << transition_probability << ",";
			std::cout << " Rolled " << choice << ")" << std::endl;
			spin_config[no_site.site_index_from_position(i1,j1,k1)] = old_sz_1;
			spin_config[no_site.site_index_from_position(i2,j2,k2)] = old_sz_2;
		}
	}
	else{
		//std::cout << "Spin config " << "(" << new_sz_1 << ", " << new_sz_2 << ") invalid..." << std::endl;;
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
	//std::cout << "Creating up auxiliary..." << std::endl;
	double old_wavefunction = -1;
	AuxMPS VUi(AuxType::VU);//Create the up auxiliary (won't have some links between certain sites though)
	for(int j = 0; j < psi.Ny(); j++){
		auto left_links = itensor::commonInds(psi._site_tensors[0][j][0], psi._site_tensors[0][j][1]);
		auto [left_tensor, sing_vals, right_tensor] = itensor::svd(psi._site_tensors[0][j][0], left_links, {"MaxDim", psi.Dc()});
		VUi.add_tensor(left_tensor);
		itensor::ITensor right_tensor_combined = sing_vals*right_tensor;
		VUi.add_tensor(right_tensor_combined);
	}
	//Begin sweeping along the rows
	auto vd_it = vd_list.begin();

	for(int i = 0; i < psi.Nx(); i++){
		//std::cout << "Sweeping row " << i << std::endl;
		vd_it++;
		//Create right auxiliary tensors
		std::vector<itensor::ITensor> vr_auxiliaries(psi.Ny()*2+1);
		vr_auxiliaries[2*psi.Ny()] = itensor::ITensor(1);
		for(int j = psi.Ny()-1; j >= 0; j--){
			for(int k = 2; k >= 1; k--){
				int J = 2*j+k-1; //J is the index of the VU MPS tensors and left/right auxiliary tensors, J-1 is the index of the VD tensors
				//Print(vr_auxiliaries[J+1]);
				//Print(VUi.MPS[J]);
				//Print(psi._site_tensors[i][j][k]);

				vr_auxiliaries[J] = (vr_auxiliaries[J+1]*VUi.MPS[J])*psi._site_tensors[i][j][k];
				if(J>0){
					//Print(vd_it->MPS[J-1]);
					vr_auxiliaries[J] *= vd_it->MPS[J-1];
				}
			}
		}
		//Evaluate the old wavefunction
		old_wavefunction = itensor::norm(vr_auxiliaries[0]);
		//std::cout << "Old wavefunction: " << old_wavefunction << std::endl;
		itensor::ITensor vl_auxiliary(1);
		//Sweep the horizontal (1-2) links from left to right
		for(int j = 0; j < psi.Ny(); j++){

			int J = 2*j;
			itensor::ITensor vd_aux_1(1);
			if(j > 0){vd_aux_1 = vd_it->MPS[J-1];}
			//Print(vl_auxiliary);
			old_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,1,i,j,2, vl_auxiliary, vr_auxiliaries[J+2],VUi.MPS[J],VUi.MPS[J+1],vd_aux_1, vd_it->MPS[J], old_wavefunction, r, WAVEFUNCTION_NORMALIZATION_CONSTANT);
			//std::cout << "Tensor data for " << i << ", " << j << ", 1-2..." << std::endl;
			//Print(VUi.MPS[J]);
			//Print(psi._site_tensors[i][j][1]);
			//Print(vd_aux_1);
			//Print(VUi.MPS[J+1]);
			//Print(psi._site_tensors[i][j][2]);
			//Print(vd_it->MPS[J]);
			//Print(vr_auxiliaries[J+2]);
			if(j < psi.Ny()-1){
				//Print(vl_auxiliary);
				old_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,2,i,j+1,1, vl_auxiliary, vr_auxiliaries[J+3],VUi.MPS[J+1],VUi.MPS[J+2],vd_it->MPS[J], vd_it->MPS[J+1], old_wavefunction, r, WAVEFUNCTION_NORMALIZATION_CONSTANT);
				//std::cout << "Tensor data for " << i << ", " << j << ", 2-1+..." << std::endl;
				//Print(VUi.MPS[J+1]);
				//Print(psi._site_tensors[i][j][2]);
				//Print(vd_it->MPS[J]);
				//Print(VUi.MPS[J+2]);
				//Print(psi._site_tensors[i][j+1][1]);
				//Print(vd_it->MPS[J+1]);
				//Print(vr_auxiliaries[J+3]);
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
		if(i < psi.Nx()-1){
			//std::cout << "Creating new VU auxiliary... Contracting row i..." << std::endl;
			//Create new VU auxiliary MPS
			std::vector<itensor::ITensor> row_i_contracted;
			for(int j = 0; j < psi.Ny(); j++){
				int J = 2*j;
				row_i_contracted.push_back(VUi.MPS[J]*psi._site_tensors[i][j][1]);
				row_i_contracted.push_back(VUi.MPS[J+1]*psi._site_tensors[i][j][2]);
			}
			//std::cout << "Truncating row i..." << std::endl;
			//Contract the bond dimensions of the contracted row i
			for(int J = 0; J < 2*psi.Ny()-1; J++){
				auto forward_links = itensor::commonInds(row_i_contracted[J], row_i_contracted[J+1]);
				if(itensor::length(forward_links) == itensor::length(row_i_contracted[J].inds())){
					row_i_contracted[J+1] *= row_i_contracted[J];
					row_i_contracted[J] = itensor::ITensor(1);
				}
				else{
					auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_i_contracted[J], forward_links, {"MaxDim", psi.Dc()});
					row_i_contracted[J+1] *= (forward_tensor*sing_vals);
					row_i_contracted[J] = back_tensor;
				}
				
			}

			//std::cout << "Fully truncated row i: " << std::endl;
			//for(int J = 0; J < row_i_contracted.size(); J++){Print(row_i_contracted[J]);}
			//std::cout << "Original row i+1: " << std::endl;
			//for(int j = 0; j < psi.Ny(); j++){Print(psi._site_tensors[i+1][j][0]);}
			//std::cout << "Applying row i to row i+1..." << std::endl;
			std::vector<itensor::ITensor>row_ip_unsplit;
			//Apply the tensors to the next row down, then split them
			for(int j = 0; j < psi.Ny(); j++){
				int J = 2*j;
				if(j > 0){
					row_ip_unsplit[j-1] *= row_i_contracted[J];
				}
				row_ip_unsplit.push_back(row_i_contracted[J+1]*psi._site_tensors[i+1][j][0]);
			}
			//std::cout << "Unsplit row i+1: " << std::endl;
			//for(int J = 0; J < row_ip_unsplit.size(); J++){Print(row_ip_unsplit[J]);}
			//std::cout << "Splitting row i+1..." << std::endl;
			for(int j = 0; j < psi.Ny(); j++){
				int J = 2*j;
				itensor::IndexSet forward_inds;
				itensor::Index always_forward = itensor::commonIndex(row_ip_unsplit[j], psi._site_tensors[i+1][j][2]);
				if(j == psi.Ny()-1){forward_inds = itensor::IndexSet(always_forward);}
				else{forward_inds = itensor::IndexSet(always_forward, itensor::commonIndex(row_ip_unsplit[j], row_ip_unsplit[j+1]));}
				auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_ip_unsplit[j], forward_inds, {"MaxDim", psi.Dc()});
				VUi.MPS[J] = back_tensor;
				VUi.MPS[J+1] = forward_tensor*sing_vals;
				//Print(VUi.MPS[J]);
				//Print(VUi.MPS[J+1]);
			}
		}
	}	
	return old_wavefunction;
}

double sample_s_direction(MCKPEPS &psi_sites, std::vector<int> &spin_config, Randomizer &r){
	SpinConfigPEPS config(psi_sites, spin_config, WAVEFUNCTION_NORMALIZATION_CONSTANT);
	NoSitePEPS psi = psi_sites.contract(config);
	std::list<AuxMPS> sd_list = psi.get_sd_auxiliaries();
	double old_wavefunction = -1;
	AuxMPS SUi(AuxType::SU);//Create the up auxiliary (won't have some links between certain sites though)
	for(int i = 0; i < psi.Nx(); i++){
		auto left_links = itensor::commonInds(psi._site_tensors[i][0][1], psi._site_tensors[i][0][0]);
		auto [left_tensor, sing_vals, right_tensor] = itensor::svd(psi._site_tensors[i][0][1], left_links, {"MaxDim", psi.Dc()});
		SUi.add_tensor(left_tensor);
		itensor::ITensor right_tensor_combined = sing_vals*right_tensor;
		SUi.add_tensor(right_tensor_combined);
	}
	//Begin sweeping along the rows
	auto sd_it = sd_list.begin();
	for(int j = 0; j < psi.Ny(); j++){
		//std::cout << "Sweeping column " << j << std::endl;
		sd_it++;
		//Create right auxiliary tensors
		std::vector<itensor::ITensor> sr_auxiliaries(psi.Nx()*2+1);
		sr_auxiliaries[2*psi.Nx()] = itensor::ITensor(1);
		for(int i = psi.Nx()-1; i >= 0; i--){
			for(int half_k = 1; half_k >= 0; half_k--){
				int k = 2*half_k;
				int I = 2*i+half_k;//I is the index of the VU MPS tensors and left/right auxiliary tensors, I-1 is the index of the VD tensors
				//std::cout << "Tensors at " << I << ": " << std::endl;
				//Print(sr_auxiliaries[I+1]);
				//Print(SUi.MPS[I]);
				//Print(psi._site_tensors[i][j][k]);
				sr_auxiliaries[I] = (sr_auxiliaries[I+1]*SUi.MPS[I])*psi._site_tensors[i][j][k];
				if(I > 0){
					//Print(sd_it->MPS[I-1]);
					sr_auxiliaries[I] *= sd_it->MPS[I-1];
				}
			}
		}
		//Evaluate the old wavefunction
		old_wavefunction = itensor::norm(sr_auxiliaries[0]);
		//std::cout << "original wavefunction: " << old_wavefunction << std::endl;
		//double old_num_choices = config.num_choices();
		itensor::ITensor sl_auxiliary(1);
		//Sweep the short direction (0-2) links from left to right
		for(int i = 0; i < psi.Nx(); i++){
			int I = 2*i;
			itensor::ITensor sd_aux_1(1);
			if(i > 0){sd_aux_1 = sd_it->MPS[I-1];}
			old_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,0,i,j,2, sl_auxiliary, sr_auxiliaries[I+2], SUi.MPS[I],SUi.MPS[I+1],sd_aux_1, sd_it->MPS[I], old_wavefunction, r, WAVEFUNCTION_NORMALIZATION_CONSTANT);
			if(i < psi.Nx()-1){
				old_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,2,i+1,j,0, sl_auxiliary, sr_auxiliaries[I+3],SUi.MPS[I+1],SUi.MPS[I+2],sd_it->MPS[I], sd_it->MPS[I+1], old_wavefunction, r, WAVEFUNCTION_NORMALIZATION_CONSTANT);
			}
		}
		if(j < psi.Ny()-1){
			//Create new SU auxiliary MPS
			std::vector<itensor::ITensor> row_j_contracted;
			for(int i = 0; i < psi.Nx(); i++){
				int I = 2*i;
				row_j_contracted.push_back(SUi.MPS[I]*psi._site_tensors[i][j][0]);
				row_j_contracted.push_back(SUi.MPS[I+1]*psi._site_tensors[i][j][2]);
			}
			//Truncate the bond dimensions of the contracted row i
			for(int I = 0; I < 2*psi.Nx()-1; I++){
				auto forward_links = itensor::commonInds(row_j_contracted[I], row_j_contracted[I+1]);
				if(itensor::length(forward_links) == itensor::length(row_j_contracted[I].inds())){
					row_j_contracted[I+1] *= row_j_contracted[I];
					row_j_contracted[I] = itensor::ITensor(1);
				}
				else{
					auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_j_contracted[I], forward_links, {"MaxDim", psi.Dc()});
					row_j_contracted[I+1] *= (forward_tensor*sing_vals);
					row_j_contracted[I] = back_tensor;
				}
			}
			std::vector<itensor::ITensor>row_jp_unsplit;
			//Apply the tensors to the next row down, then split them
			for(int i = 0; i < psi.Nx(); i++){
				int I = 2*i;
				if(i > 0){
					row_jp_unsplit[i-1] *= row_j_contracted[I];
				}
				row_jp_unsplit.push_back(row_j_contracted[I+1]*psi._site_tensors[i][j+1][1]);
			}
			for(int i = 0; i < psi.Nx(); i++){
				int I = 2*i;
				itensor::IndexSet forward_inds;
				itensor::Index always_forward = itensor::commonIndex(row_jp_unsplit[i], psi._site_tensors[i][j+1][2]);
				if(i == psi.Nx()-1){forward_inds = itensor::IndexSet(always_forward);}
				else{forward_inds = itensor::IndexSet(always_forward, itensor::commonIndex(row_jp_unsplit[i], row_jp_unsplit[i+1]));}
				auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_jp_unsplit[i], forward_inds, {"MaxDim", psi.Dc()});
				SUi.MPS[I] = back_tensor;
				SUi.MPS[I+1] = forward_tensor*sing_vals;
			}
		}
	}	
	return old_wavefunction;
}

double sample_l_direction(MCKPEPS &psi_sites, std::vector<int> &spin_config, Randomizer &r){
	SpinConfigPEPS config(psi_sites, spin_config, WAVEFUNCTION_NORMALIZATION_CONSTANT);
	NoSitePEPS psi = psi_sites.contract(config);
	std::list<AuxMPS> ld_list = psi.get_ld_auxiliaries();
	double old_wavefunction = -1;
	AuxMPS LUi(AuxType::LU);//Create the up auxiliary (won't have some links between certain sites though)
	itensor::ITensor blank1(1);
	itensor::ITensor blank2(1);
	LUi.add_tensor(blank1);//Very first up auxiliary is nonexistent, just two scalar 1 tensors
	LUi.add_tensor(blank2);
	auto ld_it = ld_list.begin();
	for(int h = 0; h < psi.Nx() + psi.Ny() - 1; h++){//h = i+j
		//std::cerr << "h=" << h << " sampling..." << std::endl;
		int imin = std::max(0, h - psi.Ny()+1);
		int imax = std::min(psi.Nx(), h+1);
		//Create right auxiliary tensors
		int Nd = imax-imin;
		//std::cerr << "Creating right auxiliaries..." << std::endl;
		std::vector<itensor::ITensor> lr_auxiliaries(2*Nd+1);
		lr_auxiliaries[2*Nd] = itensor::ITensor(1);
		for(int i = imax-1; i >= imin; i--){
			for(int k = 1; k >= 0; k--){
				int H = 2*(i-imin) + k;
				int j = h - i;
				/*Print(lr_auxiliaries[H+1]);
				Print(ld_it->MPS[H]);
				Print(psi._site_tensors[i][j][k]);
				Print(LUi.MPS[H]);*/
				lr_auxiliaries[H] = ((lr_auxiliaries[H+1]*ld_it->MPS[H])*psi._site_tensors[i][j][k])*LUi.MPS[H];
			}
		}
		//Print(lr_auxiliaries[0]);
		old_wavefunction = itensor::norm(lr_auxiliaries[0]);
		//std::cerr << "Testing bonds..." << std::endl;
		itensor::ITensor ll_auxiliary(1);
		for(int i = imin; i < imax; i++){
			int H = 2*(i-imin);
			int j = h-i;
			old_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,0,i,j,1,ll_auxiliary, lr_auxiliaries[H+2], LUi.MPS[H], LUi.MPS[H+1], ld_it->MPS[H], ld_it->MPS[H+1], old_wavefunction, r, WAVEFUNCTION_NORMALIZATION_CONSTANT);
			if(i < imax-1){
				old_wavefunction = test_bond(psi, psi_sites, spin_config, i,j,1,i+1,j-1,0,ll_auxiliary, lr_auxiliaries[H+3], LUi.MPS[H+1], LUi.MPS[H+2], ld_it->MPS[H+1], ld_it->MPS[H+2], old_wavefunction, r, WAVEFUNCTION_NORMALIZATION_CONSTANT);
			}
		}
		//Create new LU auxiliary MPS
		if(h < psi.Nx()+psi.Ny()-2){
			//std::cout << "Contracting old LU auxiliary..." << std::endl;
			std::vector<itensor::ITensor> row_h_contracted;
			for(int i = imin; i < imax; i++){
				//std::cout << "  i=" << i << std::endl;
				int H = 2*(i-imin);
				int j = h-i;
				row_h_contracted.push_back(LUi.MPS[H]*psi._site_tensors[i][j][0]);
				row_h_contracted.push_back(LUi.MPS[H+1]*psi._site_tensors[i][j][1]);
			}
			//std::cout << "Truncating bond dimensions..." << std::endl;
			//Truncate the bond dimensions of the contracted row h
			for(int H = 0; H < 2*Nd-1; H++){
				//std::cout << "  H=" << H << std::endl;
				auto forward_links = itensor::commonInds(row_h_contracted[H], row_h_contracted[H+1]);
				auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_h_contracted[H], forward_links, {"MaxDim", psi.Dc()});
				row_h_contracted[H] = back_tensor;
				row_h_contracted[H+1] *= (sing_vals*forward_tensor);
			}
			//std::cout << "Creating unsplit version of row h+1..." << std::endl;
			std::vector<itensor::ITensor> row_hp_unsplit;
			for(int i = imin; i < imax; i++){
				//std::cout << "  i=" << i << std::endl;
				int H = 2*(i-imin);
				int j = h-i;
				row_hp_unsplit.push_back((psi._site_tensors[i][j][2]*row_h_contracted[H])*row_h_contracted[H+1]);
			}

			LUi.clear();

			//Split the row h+1, creating scalar 1-tensors when necessary
			//std::cout << "Creating split version of row h+1..." << std::endl;
			for(int i = imin; i < imax; i++){
				//std::cerr << "  i=" << i << "...";
				int j = h-i;
				if((i == imin) && (j+1 < psi.Ny())){ //If the next diagonal has an (imin, j+1) space, add an extra 1-tensor at the front
					itensor::ITensor blank(1);
					LUi.add_tensor(blank);
				}
				if((j+1 < psi.Ny()) && (i+1 < psi.Nx())){
					//std::cerr << " gathering inds...";
					auto forward_links = itensor::commonInds(row_hp_unsplit[i-imin], psi._site_tensors[i+1][j][0]);
					if(i < imax-1){forward_links = itensor::unionInds(forward_links, itensor::commonInds(row_hp_unsplit[i-imin], row_hp_unsplit[i-imin+1]));}
					//std::cerr << "svd...";
					auto [forward_tensor, sing_vals, back_tensor] = itensor::svd(row_hp_unsplit[i-imin], forward_links, {"MaxDim", psi.Dc()});
					//std::cerr << "adding tensors...";
					LUi.add_tensor(back_tensor);
					itensor::ITensor forward_combined = forward_tensor*sing_vals;
					LUi.add_tensor(forward_combined);
				}
				else{
					LUi.add_tensor(row_hp_unsplit[i-imin]);
				}
				if((i == imax-1) && (i+1 < psi.Nx())){ //If the next diagonal has an (imax+1, jmin) space, add an extra 1-tensor at the back
					itensor::ITensor blank(1);
					LUi.add_tensor(blank);
				}
				//std::cerr << "New LU length: " << LUi.length << std::endl;
			}
		}
		ld_it++;
	}

	return old_wavefunction;
	
}

#endif