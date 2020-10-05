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
	int standard_dims = input.testInteger("D", 2);
	int max_truncation_dims = input.testInteger("Dc", 4);
	int physical_dims = input.testInteger("physical_dims", 4);

	int num_sites = Nx*Ny*UNIT_CELL_SIZE;

	std::vector<itensor::Index> sites_vector(num_sites);
	for(int i = 0; i < num_sites; i++){
		sites_vector[i] = itensor::Index(physical_dims,"Site,n="+std::to_string(i+1));
	}
	itensor::IndexSet sites(sites_vector);

	auto PEPS1 = MCKPEPS(sites, Nx, Ny, standard_dims, max_truncation_dims);

	std::vector<int> spin_config(num_sites, 0);
	Randomizer r;
	randomize_in_sector(spin_config, physical_dims, r.gen, r.dist);
	SpinConfigPEPS scp(PEPS1, spin_config, 1);
	double inner_product = PEPS1.inner_product(scp);

	NoSitePEPS nsp = PEPS1.contract(scp);
	std::cerr << "Getting environments..." << std::endl;
	std::vector<itensor::ITensor> envs = nsp.environments(PEPS1.site_indices, spin_config);

	for(int site_index = 0; site_index < envs.length(); site_index++){
		itensor::ITensor product = envs[site_index]*(*PEPS1.tensor_at(site_index));
		if(itensor::length(product.inds()) == 0){
			std::cerr << "Product #" << site_index << ": " << itensor::elt(product) << " Original: " << inner_product << std::endl;
		}
		else{
			std::cerr << "Error: Environment " << site_index << " index mismatch" << std::endl;
			Print(envs[site_index]);
			Print(*PEPS1.tensor_at(site_index));
		}
	}

	return 0;
}