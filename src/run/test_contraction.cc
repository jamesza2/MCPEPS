#include "itensor/all.h"
#include "../headers/input.h"
#include "../headers/mcpeps.h"
#include <ctime>
#include <cmath>
#include <complex>

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
	std::string log_file = input.testString("log_file", "");
	int standard_dims = input.testInteger("D", 2);
	int max_truncation_dims = input.testInteger("Dc", 4);

	int num_sites = Nx*Ny*UNIT_CELL_SIZE;
	std::vector<itensor::Index> sites_vector(num_sites);
	for(int i = 0; i < num_sites; i++){
		sites_vector[i] = itensor::Index(2);
	}
	itensor::IndexSet sites(sites_vector);
	auto PEPS1 = MCKPEPS(sites, Nx, Ny, standard_dims, max_truncation_dims);
	auto PEPS2 = MCKPEPS(sites, Nx, Ny, standard_dims, max_truncation_dims);
	PEPS1.set_log_file(log_file);
	double brute_force_inner_product = PEPS1.brute_force_inner_product(PEPS2);
	double inner_product = PEPS1.inner_product(PEPS2);
	std::cerr << "Inner Product: " << inner_product << std::endl;
	std::cerr << "Brute Force Inner Product: " << brute_force_inner_product << std::endl;
	return 0;
}