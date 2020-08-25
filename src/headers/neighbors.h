#ifndef NEIGHBORS_LIST
#define NEIGHBORS_LIST


const int UNIT_CELL_SIZE = 3;


//A class that stores neighbor info about the kag lattice

class Neighbors{
	protected:
		int _Nx;
		int _Ny;
		int _num_sites;
	public:
		std::vector<std::vector<std::pair<int, int>>> bond_matrices;
		std::vector<std::vector<int>> nn_bonds;
		Neighbors(){}
		Neighbors(int Nx, int Ny){
			_Nx = Nx;
			_Ny = Ny;
			_num_sites = Nx*Ny*UNIT_CELL_SIZE;
			//Set up bonds
			for(int s1 = 0; s1 < _num_sites; s1++){
				nn_bonds.push_back(std::vector<int>());
				bond_matrices.push_back(std::vector<std::pair<int, int>>());
			}
			set_up_bonds();
		}
		void set_dimensions(int Nx, int Ny){
			_Nx = Nx;
			_Ny = Ny;
			_num_sites = Nx*Ny*UNIT_CELL_SIZE;
			//Set up bonds
			for(int s1 = 0; s1 < _num_sites; s1++){
				nn_bonds.push_back(std::vector<int>());
				bond_matrices.push_back(std::vector<std::pair<int, int>>());
			}
			set_up_bonds();
		}
		std::vector<std::pair<int, int>> at(int site){
			return bond_matrices[site];
		}
		std::vector<int> nn_at(int site){
			return nn_bonds[site];
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
					if(bond_type == 1){
						nn_bonds[site_1].push_back(site_2);
						nn_bonds[site_2].push_back(site_1);
					}
				}
			}
		}
};

#endif