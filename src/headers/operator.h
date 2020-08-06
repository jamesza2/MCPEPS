#ifndef MCPEPS_OPERATOR
#define MCPEPS_OPERATOR


const int UNIT_CELL_SIZE = 3;

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

		double eval(std::vector<int> &spin_config_1, std::vector<int> &spin_config_2){
			return 0;
		}
}

class Heisenberg : public MCOperator{
	protected:
		double _J1; //NN bonds
		double _J2; //NNN bonds
		double _Jd; //d bonds
		double _Jz; //z-component anisotropy
		//For each site i, bond_matrices stores all bonds i-j
		//First int is the site j
		//Second int is the bond type (1 for NN, 2 for NNN, 3 for d)
		std::vector<std::vector<std::pair<int, int>>> bond_matrices;
	public:
		void set_Jz(std::map<std::string, double> &J_vals){
			if(J_vals.count("J1")){
				_J1 = J_vals["J1"];
			}
			if(J_vals.count("J2")){
				_J2 = J_vals["J2"];
			}
			if(J_vals.count("Jd")){
				_Jd = J_vals["Jd"];
			}
			if(J_vals.count("Jz")){
				_Jz = J_vals["Jz"];
			}
		}

		Heisenberg(int Nx, int Ny, int physical_dimension) : MCOperator(Nx, Ny, physical_dimension){
			//Set up bonds
			for(int s1 = 0; s1 < _num_sites; s1++){
				bond_matrices.push_back(std::vector<std::pair<int, int>>());
			}
			set_up_bonds();
		}


		double eval(std::vector<int> &spin_config_1, std::vector<int> &spin_config_2){
			//First check where spin_config_1 and spin_config_2 differ
			//Bonds can only cover all the difference point
			std::vector<int> different_sites;
			for(int site = 0; site < _num_sites; site++){
				if(spin_config_1[site] != spin_config_2[site]){
					different_sites.push_back(site);
				}
			}
			double matrix_element = 0;
			double s = 0.5*(_d-1);
			std::vector<double> J{0, _J1, _J2, _J3};
			//For each bond:
			//Check if it covers all differences
			//Split into Sxy and Sz components
			for(int site_1 = 0; site_1 < _num_sites; site_1 ++){
				for(auto bond : bond_matrices[site_1]){
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
						matrix_element += 0.5*J[bond.second]*_Jz*m1*m2;
					}
				}
				return matrix_element;
			}
		}
	private:
		void set_up_bonds(){
			//Set up NN bonds (i,j,k)
			set_up_bond(0,1,0,0,1);
			set_up_bond(1,2,0,0,1);
			set_up_bond(2,0,0,0,1);//In-triangle bonds
			set_up_bond(2,1,0,1,1);//Rightward bond
			set_up_bond(1,0,1,-1,1);//Downward-left bond
			set_up_bond(2,0,1,0,1);//Downward bond

			//Set up NNN bonds
			set_up_bond(0,1,0,1,2);
			set_up_bond(0,2,-1,1,2);
			set_up_bond(1,0,1,0,2);
			set_up_bond(1,2,1,-1,2);
			set_up_bond(2,0,0,1,2);
			set_up_bond(2,1,1,0,2);

			//Set up d bonds
			set_up_bond(0,0,0,1,3);
			set_up_bond(1,1,1,0,3);
			set_up_bond(2,2,1,-1,3);

			//TODO: Verify bonds
		}

		int site_index_from_position(int i, int j, int k){
			return i*_Ny*UNIT_CELL_SIZE + j*UNIT_CELL_SIZE + k;
		}
		//delta_i = i2-i1
		//delta_j = j2-j1
		void set_up_bond(int k1, int k2, int delta_i, int delta_j, int bond_type){
			for(int i1 = std::max(0, -delta_i); i1 < std::min(_Nx, _Nx-delta_i); i1++){
				for(int j1 = std::max(0, -delta_j); j1 < std::min(_Ny, _Ny-delta_j); j1++){
					int site_1 = site_index_from_position(i1, j1, k1);
					int site_2 = site_index_from_position(i1+delta_i, j1+delta_j,k2);
					if(site_1 > site_2){
						int temp = site_2;
						site_2 = site_1;
						site_1 = temp;
					}
					std::pair<int, int> bond(site_2, bond_type);
					bond_matrices[site_1].push_back(bond);
				}
			}
		}
}

class Sz2 : public MCOperator{
	//Measures sz^2 for spin config 1, ignores spin config 2
	double eval(std::vector<int> &spin_config_1, std::vector<int> &spin_config_2){
		double elem = 0;
		double s = 0.5*(_d-1);
		for(int spin_1 : spin_config_1){
			double m1 = spin_1-s;
			elem += m1*m1;
		}
		return elem;
	}
}

#endif