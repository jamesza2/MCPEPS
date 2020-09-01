#include "itensor/all.h"
#include "../headers/input.h"
#include "../headers/output.h"
#include "../headers/mcpeps.h"
#include "../headers/inner_sampling.h"
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

std::string check_common_indices(itensor::ITensor tensor_1, std::string tensor_name_1, itensor::ITensor tensor_2, std::string tensor_name_2, int ideal_common_indices = 1){
	std::string result = "";
	itensor::IndexSet common_indices = itensor::commonInds(tensor_1,tensor_2);
	int num_common_indices = itensor::length(common_indices);
	if(num_common_indices!=ideal_common_indices){
		result += tensor_name_1 + " has incorrect common indices with " + tensor_name_2 + ", at " + std::to_string(num_common_indices) + ": ";
		for(auto ind : common_indices){
			result += itensor::id(ind) + " ";
		}
		result += "\n";
	}
	return result;
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
	int physical_dims = input.testInteger("physical_dims", 4);

	std::string version = "_";
	version += std::to_string(Nx) + "x" + std::to_string(Ny);
	version += "_D"+std::to_string(standard_dims);
	version += "_X"+std::to_string(max_truncation_dims);
	version += "_d"+std::to_string(physical_dims);
	if(log_file_name == "AUTO"){ //../../logs/aux_test_{Nx}x{Ny}_D{D}_X{Chi}_d{d}_{}trials
		log_file_name = "../../logs/aux_test" + version;
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
	PEPS1.set_log_file(log_file_name);
	PEPS1.print_self();

	NoSitePEPS PEPSC = PEPS1.contract(PEPS2);
	std::cout << "Contracted PEPS created..." << std::endl;
	PEPSC.print_self();
	std::cout << "printed..." << std::endl;
	auto vd = PEPSC.get_vd_auxiliaries();
	auto sd = PEPSC.get_sd_auxiliaries();
	auto ld = PEPSC.get_ld_auxiliaries();

	std::string final_result = "";

	std::cout << "VD auxiliary test..." << std::endl;
	//VD auxiliary test
	auto vd_it = vd.begin();
	if(vd_it->length != 0){
		final_result += "VD[0] Length test failed: " + std::to_string(vd_it->length) +"\n";
	}
	int original_i = 0;
	while(vd_it != vd.end()){
		vd_it ++;
		original_i ++;
		for(int original_j = 0; original_j < Ny; original_j++){
			auto vd_tensor_1 = vd_it->MPS[2*original_j];
			std::string site_tensor_name = "SITE[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][0]";
			std::string upper_left_site_tensor_name = "SITE[" + std::to_string(original_i-1) + "][" + std::to_string(original_j) + "][2]";
			std::string vd_tensor_name_1 = "VD[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][0]";
			final_result += check_common_indices(PEPSC.site_tensor(original_i, original_j, 0), site_tensor_name, vd_tensor_1, vd_tensor_name_1, 1);
			final_result += check_common_indices(PEPSC.site_tensor(original_i-1,original_j,2), upper_left_site_tensor_name, vd_tensor_1, vd_tensor_name_1, 1);
			
			if(original_j != Ny-1){
				auto vd_tensor_2 = vd_it->MPS[2*original_j+1];
				std::string vd_tensor_name_2 = "VD[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][1]";
				std::string upper_right_site_tensor_name = "SITE[" + std::to_string(original_i-1) + "][" + std::to_string(original_j+1) + "][1]";
				final_result += check_common_indices(vd_tensor_1, vd_tensor_name_1, vd_tensor_2, vd_tensor_name_2, 1);
				final_result += check_common_indices(vd_tensor_2, vd_tensor_name_2, PEPSC.site_tensor(original_i,original_j,0), site_tensor_name, 1);
				final_result += check_common_indices(PEPSC.site_tensor(original_i-1,original_j+1,1), upper_right_site_tensor_name, vd_tensor_2, vd_tensor_name_2, 1);
			}
		}
	}

	std::cout << "SD auxiliary test..." << std::endl;
	//SD auxiliary test
	auto sd_it = sd.begin();
	if(sd_it->length != 0){
		final_result += "SD[0] Length test failed: " + std::to_string(sd_it->length) +"\n";
	}
	int original_j = 0;
	while(sd_it != sd.end()){
		sd_it ++;
		original_j ++;
		for(int original_i = 0; original_i < Nx; original_i++){
			auto sd_tensor_1 = sd_it->MPS[2*original_i];
			std::string site_tensor_name = "SITE[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][1]";
			std::string upper_left_site_tensor_name = "SITE[" + std::to_string(original_i) + "][" + std::to_string(original_j-1) + "][2]";
			std::string sd_tensor_name_1 = "SD[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][0]";
			final_result += check_common_indices(PEPSC.site_tensor(original_i,original_j,1), site_tensor_name, sd_tensor_1, sd_tensor_name_1, 1);
			final_result += check_common_indices(PEPSC.site_tensor(original_i,original_j-1,2), upper_left_site_tensor_name, sd_tensor_1, sd_tensor_name_1, 1);
			
			if(original_i != Nx-1){
				auto sd_tensor_2 = sd_it->MPS[2*original_i+1];
				std::string sd_tensor_name_2 = "SD[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][1]";
				std::string upper_right_site_tensor_name = "SITE[" + std::to_string(original_i+1) + "][" + std::to_string(original_j-1) + "][0]";
				final_result += check_common_indices(sd_tensor_1, sd_tensor_name_1, sd_tensor_2, sd_tensor_name_2, 1);
				final_result += check_common_indices(sd_tensor_2, sd_tensor_name_2, PEPSC.site_tensor(original_i,original_j,1), site_tensor_name, 1);
				final_result += check_common_indices(PEPSC.site_tensor(original_i+1,original_j-1,0), upper_right_site_tensor_name, sd_tensor_2, sd_tensor_name_2, 1);
			}
		}
	}

	std::cout << "LD auxiliary test..." << std::endl;
	//LD auxiliary test
	auto ld_it = ld.begin();
	int original_h = 0;
	for(auto ld_it = ld.begin(); ld_it != ld.end(); ld_it++){
		int imin = std::max(0, original_h - Ny+1);
		int imax = std::min(Nx-1, original_h);
		for(int original_i = imin; original_i < imax; original_i++){
			int original_j = original_h - original_i;
			auto ld_tensor_1 = ld_it->MPS[2*(original_i-imin)];
			std::string site_tensor_name = "SITE[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][2]";
			std::string upper_left_site_tensor_name = "SITE[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][0]";
			std::string ld_tensor_name_1 = "LD[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][0]";
			final_result += check_common_indices(PEPSC.site_tensor(original_i,original_j,2), site_tensor_name, ld_tensor_1, ld_tensor_name_1, 1);
			final_result += check_common_indices(PEPSC.site_tensor(original_i,original_j,0), upper_left_site_tensor_name, ld_tensor_1, ld_tensor_name_1, 1);
			
			if(original_i != Nx-1){
				auto ld_tensor_2 = ld_it->MPS[2*original_i+1];
				std::string ld_tensor_name_2 = "LD[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][1]";
				std::string upper_right_site_tensor_name = "SITE[" + std::to_string(original_i) + "][" + std::to_string(original_j) + "][1]";
				final_result += check_common_indices(ld_tensor_1, ld_tensor_name_1, ld_tensor_2, ld_tensor_name_2, 1);
				final_result += check_common_indices(ld_tensor_2, ld_tensor_name_2, PEPSC.site_tensor(original_i,original_j,2), site_tensor_name, 1);
				final_result += check_common_indices(PEPSC.site_tensor(original_i,original_j,1), upper_right_site_tensor_name, ld_tensor_2, ld_tensor_name_2, 1);
			}
		}
		original_h ++;
	}
	std::cout << final_result;
	

	/*

	double total_Sz2;
	std::cerr << "Performing efficient inner product..." << std::endl;
	auto timestart = std::time(NULL);
	double inner_product = PEPS1.inner_product(PEPS2);
	for(int i = 0; i < num_sites; i++){
		auto Sz2_tensor = create_sz2_op(i, sites);
		MCKPEPS applied_PEPS = PEPS2;
		applied_PEPS.apply_spinop(i, Sz2_tensor);
		total_Sz2 += PEPS1.inner_product(applied_PEPS)/inner_product;
	}
	double efficient_time = std::difftime(std::time(NULL), timestart);
	timestart = std::time(NULL);

	std::cerr << "Performing Monte Carlo inner product..." << std::endl;
	std::vector<double> wavefunctions;
	std::vector<double> values;
	mc_sz2(PEPS1, wavefunctions, values, num_trials);

	double mc_time = std::difftime(std::time(NULL), timestart);

	*/
	/*
	Output out;
	out.addInteger("NX", Nx);
	out.addInteger("NY", Ny);
	out.addInteger("D", standard_dims);
	out.addInteger("CHI", max_truncation_dims);
	out.addInteger("NUM_TRIALS", num_trials);
	out.addDouble("DIRECT_INNER_PRODUCT", inner_product);
	out.addDouble("DIRECT_SZ2", total_Sz2);
	out.addVector("WAVEFUNCTIONS", wavefunctions);
	out.addVector("VALUES", values);

	if(out_file_name == "AUTO"){
		out_file_name = "../../out/sz2_test" + version;
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


	std::cerr << "Total Sz^2 (inner product): " << total_Sz2 << " (" << efficient_time << "s)" << std::endl;

	int num_wavefunctions = values.size();
	double later_half_average = 0;
	for(int i = num_wavefunctions/2; i < num_wavefunctions; i++){
		later_half_average += values[i];
	}
	later_half_average /= (num_wavefunctions - num_wavefunctions/2);
	std::cerr << "Later-Half Average: " << later_half_average << " (" << mc_time << "s)" << std::endl;
	//std::cerr << "Brute Force Inner Product (Old Method): " << brute_force_inner_product_old << " (" << brute_force_time_old << "s)" << std::endl;*/
	return 0;
}