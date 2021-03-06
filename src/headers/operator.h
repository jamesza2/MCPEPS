#ifndef MCPEPS_OPERATOR
#define MCPEPS_OPERATOR


#include "neighbors.h"

class PEPSop;

class MCOperator{
	//A class that keeps a certain Site set in mind and has an eval method
	//Eval takes two spin configurations and calculates the matrix element
	protected:
		int _d;
		int _Nx;
		int _Ny;
		int _num_sites;

	public:

		MCOperator(int Nx, int Ny, int physical_dimension){
			_Nx = Nx;
			_Ny = Ny;
			_num_sites = Nx*Ny*3;
			_d = physical_dimension;
		}

		virtual double eval(std::vector<int> &spin_config_1, std::vector<int> &spin_config_2){
			std::cerr << "WARNING: Null evaluator called" << std::endl;
			return 0;
		}
};

class Heisenberg : public MCOperator{
	protected:
		double _J1; //NN bonds
		double _J2; //NNN bonds
		double _Jd; //d bonds
		double _Jz; //z-component anisotropy
		//For each site i, bond_matrices stores all bonds i-j
		//First int is the site j
		//Second int is the bond type (1 for NN, 2 for NNN, 3 for d)
		Neighbors bonds;
	public:
		void set_J(const std::map<std::string, double> &J_vals){
			if(J_vals.count("J1")){
				_J1 = J_vals.at("J1");
			}
			if(J_vals.count("J2")){
				_J2 = J_vals.at("J2");
			}
			if(J_vals.count("Jd")){
				_Jd = J_vals.at("Jd");
			}
			if(J_vals.count("Jz")){
				_Jz = J_vals.at("Jz");
			}
		}

		Heisenberg(int Nx, int Ny, int physical_dimension) : MCOperator(Nx, Ny, physical_dimension){
			bonds.set_dimensions(Nx, Ny);
		}

		//Given a spin config, returns all spin configs that have a possible nonzero matrix element with it
		std::vector<std::pair<std::vector<int>, double>> possible_matrix_elements(std::vector<int> &spin_config) const{
			std::vector<std::pair<std::vector<int>, double>> possible_configs;
			possible_configs.push_back(std::make_pair(spin_config,0));//The like-like config's matrix element is the sum of all Szi*Szj matrix elements, to be computed later
			std::vector<double> J{0, _J1, _J2, _Jd};
			for(int site_1 = 0; site_1 < _num_sites; site_1++){
				for(auto bond : bonds.at(site_1)){
					//Add the spin configs who differ by +/-1 at these bonds, and add the Szi*Szj matrix element to the like-like config
					if(J[bond.second] == 0){continue;}
					int site_2 = bond.first;
					double s = 0.5*(_d-1);
					double m1 = spin_config[site_1]-s;
					double m2 = spin_config[site_2]-s;
					possible_configs[0].second += m1*m2*J[bond.second]*_Jz; //SzSz matrix element on like-like config
					//S+S- part of the matrix element
					if((spin_config[site_1] < _d-1) && (spin_config[site_2] > 0)){
						std::vector<int> new_config(spin_config);
						new_config[site_1] += 1;
						new_config[site_2] -= 1;
						double matrix_element = 0.5*J[bond.second]*std::sqrt((s*(s+1) - m1*(m1+1))*(s*(s+1) - m2*(m2-1)));
						possible_configs.push_back(std::make_pair(new_config, matrix_element));
					}
					//S-S+ part of the matrix element
					if((spin_config[site_1] > 0) && (spin_config[site_2] < _d-1)){
						std::vector<int> new_config(spin_config);
						new_config[site_1] -= 1;
						new_config[site_2] += 1;
						double matrix_element = 0.5*J[bond.second]*std::sqrt((s*(s+1) - m1*(m1-1))*(s*(s+1) - m2*(m2+1)));
						possible_configs.push_back(std::make_pair(new_config, matrix_element));
					}
				}
			}
			return possible_configs;
		}


		double eval(std::vector<int> &spin_config_1, std::vector<int> &spin_config_2) override{
			//First check where spin_config_1 and spin_config_2 differ
			//Bonds can only cover all the difference points
			std::vector<int> different_sites;
			for(int site = 0; site < _num_sites; site++){
				if(spin_config_1[site] != spin_config_2[site]){
					different_sites.push_back(site);
				}
			}
			double matrix_element = 0;
			double s = 0.5*(_d-1);
			std::vector<double> J{0, _J1, _J2, _Jd};
			//For each bond:
			//Check if it covers all differences
			//Split into Sxy and Sz components
			for(int site_1 = 0; site_1 < _num_sites; site_1 ++){
				for(auto bond : bonds.at(site_1)){
					int site_2 = bond.first;
					bool auto_zero = false;
					for(int diff : different_sites){
						if((diff != site_1) && (diff != site_2)){
							auto_zero = true;
						}
					}
					if(auto_zero){continue;}
					double m1 = spin_config_1[site_1]-s;
					double m2 = spin_config_1[site_2]-s;
					double m1p = spin_config_2[site_1]-s;
					double m2p = spin_config_2[site_2]-s;
					//S1+ S2- component
					if((m1 == m1p-1) && (m2 == m2p+1)){
						matrix_element += 0.5*J[bond.second]*std::sqrt((s*(s+1) - m1*(m1+1))*(s*(s+1) - m2*(m2-1)));
					}
					//S1- S2+ component
					if((m1 == m1p+1) && (m2 == m2p-1)){
						matrix_element += 0.5*J[bond.second]*std::sqrt((s*(s+1) - m1*(m1-1))*(s*(s+1) - m2*(m2+1)));
					}
					//S1z S2z component
					if((m1 == m1p) && (m2 == m2p)){
						matrix_element += J[bond.second]*_Jz*m1*m2;
					}
				}
				return matrix_element;
			}
		}
		PEPSop toPEPSop() const;
};

class Sz2 : public MCOperator{
	public:
		Sz2(int Nx, int Ny, int physical_dimension) : MCOperator(Nx, Ny, physical_dimension) {}
		//Measures sz^2 for spin config 1, ignores spin config 2
		double eval(std::vector<int> &spin_config_1, std::vector<int> &spin_config_2) override{
			double elem = 0;
			double s = 0.5*(_d-1);
			for(int spin_1 : spin_config_1){
				double m1 = spin_1-s;
				elem += m1*m1;
			}
			return elem;
		}
};

#endif