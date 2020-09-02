#ifndef mcpeps
#define mcpeps

#include "itensor/all.h"
#include "operator.h"
#include "neighbors.h"
#include "auxiliary_mps.h"
#include <random>
#include <ctime>
#include <cmath>
#include <complex>
#include <map>
#include <list>

//Normal: Uses the previous 2 total weights to estimate energy, then attempts to make the total walker weight num_walkers
//Constant: Uses a single estimate for the energy and nothing else
//Antitrunc: Just tries to compensate for the truncation

//PEPS with no site indices
class NoSitePEPS
{
	protected:
		int _Nx;
		int _Ny;
		int _num_sites;
		int _D;
		int _Dc;
		
		std::vector<std::vector<std::vector<itensor::ITensor>>> _site_tensors;
		std::map<int, itensor::Index> _link_indices;

	public:
		std::string _log_file;
		Neighbors bonds;
		NoSitePEPS(){
			_Nx = 0;
			_Ny = 0;
			_num_sites=0;
		}
		NoSitePEPS(int input_Nx,
			int input_Ny,
			int input_max_bd,
			int input_max_truncation_bd,
			itensor::Args const& args = itensor::Args::global())
		{
			_num_sites = input_Nx*input_Ny*UNIT_CELL_SIZE;
			_Nx = input_Nx;
			_Ny = input_Ny;
			bonds.set_dimensions(input_Nx, input_Ny);
			_D = input_max_bd;
			_log_file = "";
			_Dc = input_max_truncation_bd;
			bool randomize = args.getBool("RandomizeSites", true);

			//std::cerr << "Creating link indices..." << std::endl;
			create_link_indices();
			//std::cerr << "Creating site tensors..." << std::endl;
			create_site_tensors(randomize);
		}

		NoSitePEPS(int input_Nx,
			int input_Ny,
			int input_max_bd,
			int input_max_truncation_bd,
			std::vector<std::vector<std::vector<itensor::ITensor>>> &input_site_tensors,
			std::map<int, itensor::Index> &input_link_indices)
		{
			_num_sites = input_Nx*input_Ny*UNIT_CELL_SIZE;
			_Nx = input_Nx;
			_Ny = input_Ny;
			bonds.set_dimensions(input_Nx, input_Ny);
			_D = input_max_bd;
			_log_file = "";
			_Dc = input_max_truncation_bd;
			_site_tensors = input_site_tensors;
			_link_indices = input_link_indices;
		}

		void set_log_file(std::string log_file){
			_log_file = log_file;
		}

		int size(){ return _num_sites; }
		int Dc(){ return _Dc; }
		int Chi(){ return _Dc; }
		int Nx(){ return _Nx; }
		int Ny(){ return _Ny; }

		//Prime all the link indices
		void prime(int inc = 1){
			for(auto link_iterator = _link_indices.begin(); link_iterator != _link_indices.end(); link_iterator++){
				link_iterator->second.prime(inc);
			}
			for(int i = 0; i < _site_tensors.size(); i++){
				for(int j = 0; j < _site_tensors[i].size(); j++){
					for(int k = 0; k < _site_tensors[i][j].size(); k++){
						_site_tensors[i][j][k].prime(inc, "Link");
					}
				}
			}
		}

		//Get the auxiliary MPS's for the vertical, i.e. i, direction
		//Warning: Only do this on a PEPS with no site indices, to avoid exponentially large tensors at the end
		std::list<AuxMPS> get_vd_auxiliaries(){
			bool _log = (_log_file != "");
			std::streambuf *coutbuf = std::cout.rdbuf();
			std::ofstream log_file_stream(_log_file, std::ofstream::app);
			if(_log){
				std::cout.rdbuf(log_file_stream.rdbuf());
			}
			std::list<AuxMPS> out;
			std::vector<itensor::ITensor> unsplit_MPS;
			std::vector<itensor::ITensor> previous_row; //Previous row is the row formed from (i,:,1) & (i,:,2) at the start of the contraction step
			for(int j = 0; j < _Ny; j++){
				previous_row.push_back(_site_tensors[_Nx-1][j][1]);
				previous_row.push_back(_site_tensors[_Nx-1][j][2]);
			}

			//Contract upwards and add each SVD split step into the AuxMPS list. Does not add the very last one.
			for(int i = _Nx-1; i > 0; i--){
				std::cout << "Contracting Row " << i << ":" << std::endl;
				unsplit_MPS.clear();
				for(int j = 0; j < _Ny; j++){
					unsplit_MPS.push_back((_site_tensors[i][j][0]*previous_row[2*j])*previous_row[2*j+1]);
				}
				AuxMPS aux(AuxType::VD);
				//Create auxiliary tensors
				for(int j = 0; j < _Ny-1; j++){
					
					itensor::Index left_upper_link = itensor::commonIndex(unsplit_MPS[j], _site_tensors[i-1][j][2]);
					itensor::ITensor sing_vals, right_tensor;
					itensor::ITensor left_tensor(left_upper_link);

					if(j>0){
						itensor::Index left_link = itensor::commonIndex(unsplit_MPS[j], unsplit_MPS[j-1]);
						left_tensor = itensor::ITensor(left_upper_link, left_link);
					}
					
					itensor::svd(unsplit_MPS[j], left_tensor, sing_vals, right_tensor, {"MaxDim", _Dc});
					std::cout << "SVD Error: " << itensor::sqr(itensor::norm(unsplit_MPS[j]-left_tensor*sing_vals*right_tensor)/itensor::norm(unsplit_MPS[j])) << std::endl;
					left_tensor *= sing_vals;
					std::cout << "Split tensors at i=" << i << ", j=" << j << std::endl;
					Print(left_tensor);
					Print(right_tensor);
					aux.add_tensor(left_tensor);
					aux.add_tensor(right_tensor);
				}
				aux.add_tensor(unsplit_MPS[_Ny-1]);
				out.push_front(aux);
				//Contract the auxiliary MPS into the next row
				previous_row.clear();
				for(int j = 0; j < _Ny; j++){
					previous_row.push_back(_site_tensors[i-1][j][1]);
					previous_row.push_back(_site_tensors[i-1][j][2]);
				}
				for(int aux_index = 0; aux_index < aux.length; aux_index++){
					previous_row[aux_index+1] *= aux.MPS[aux_index];
				}
				//Use SVD to reduce the bond dimension of the row
				//Combines pr_index-1 and pr_index, then splits again with reduced BD
				//Does not need to do the 0,1 link as the auxiliary MPS wasn't contracted into that
				for(int pr_index = 2; pr_index < previous_row.size(); pr_index++){
					itensor::ITensor combined_link_tensor = previous_row[pr_index-1]*previous_row[pr_index];
					auto left_links = itensor::commonInds(combined_link_tensor, previous_row[pr_index-1]);
					auto [left_tensor, singular_vals, right_tensor] = itensor::svd(combined_link_tensor, left_links, {"MaxDim",_Dc});
					previous_row[pr_index-1] = left_tensor;
					previous_row[pr_index] = singular_vals*right_tensor;
				}
			}
			out.push_front(AuxMPS(AuxType::NA)); //The first row is a dummy AuxMPS that shouldn't ever have to be called
			std::cout << "Row 0 empty auxMPS pushed " << std::endl;
			if(_log){
				std::cout.rdbuf(coutbuf);
				log_file_stream.close();
			}
			return out;
		}

		//Get the auxiliary MPS's for the short, i.e. j, direction
		std::list<AuxMPS> get_sd_auxiliaries(){
			bool _log = (_log_file != "");
			std::streambuf *coutbuf = std::cout.rdbuf();
			std::ofstream log_file_stream(_log_file, std::ofstream::app);
			if(_log){
				std::cout.rdbuf(log_file_stream.rdbuf());
			}
			std::list<AuxMPS> out;
			std::vector<itensor::ITensor> unsplit_MPS;
			std::vector<itensor::ITensor> previous_row; //Previous row is formed from (:,j,0) and (:,j,2) at the start of the truncation step
			for(int i = 0; i < _Nx; i++){
				previous_row.push_back(_site_tensors[i][_Ny-1][0]);
				previous_row.push_back(_site_tensors[i][_Ny-1][2]);
			}
			std::cout << "Creating starting row of size " << previous_row.size() << std::endl;
			//Contract to the left and add each SVD splitted column to the AuxMPS list. Does not contract the last column.
			for(int j = _Ny-1; j > 0; j--){
				/*std::cout << "Starting row:" << std::endl;
				for(auto row_element : previous_row){
					Print(row_element);
				}*/
				std::cout << "Contracting Column " << j << ":" << std::endl;
				unsplit_MPS.clear();
				for(int i = 0; i < _Nx; i++){
					unsplit_MPS.push_back((_site_tensors[i][j][1]*previous_row[2*i])*previous_row[2*i+1]);
				}
				AuxMPS aux(AuxType::SD);
				for(int i = 0; i <_Nx-1; i++){
					std::cout << "    At row " << i << ":" << std::endl;
					itensor::Index left_upper_link = itensor::commonIndex(unsplit_MPS[i], _site_tensors[i][j-1][2]);
					itensor::ITensor sing_vals, right_tensor;
					itensor::ITensor left_tensor(left_upper_link);
					if(i>0){
						itensor::Index left_link = itensor::commonIndex(unsplit_MPS[i], unsplit_MPS[i-1]);
						left_tensor = itensor::ITensor(left_upper_link, left_link);
					}
					itensor::svd(unsplit_MPS[i], left_tensor, sing_vals, right_tensor, {"MaxDim", _Dc});
					std::cout << "SVD Error: " << itensor::sqr(itensor::norm(unsplit_MPS[i]-left_tensor*sing_vals*right_tensor)/itensor::norm(unsplit_MPS[i])) << std::endl;
					left_tensor *= sing_vals;
					if(_log){
						std::cout << "Split tensors at i=" << i << ", j=" << j << std::endl;
						Print(left_tensor);
						Print(right_tensor);
					}
					aux.add_tensor(left_tensor);
					aux.add_tensor(right_tensor);
				}
				aux.add_tensor(unsplit_MPS[_Nx-1]);
				out.push_front(aux);
				//Contract into next row and truncate it
				previous_row.clear();
				for(int i = 0; i < _Nx; i++){
					previous_row.push_back(_site_tensors[i][j-1][0]);
					previous_row.push_back(_site_tensors[i][j-1][2]);
				}
				for(int aux_index = 0; aux_index < aux.length; aux_index++){
					previous_row[aux_index+1] *= aux.MPS[aux_index];
				}
				//Use SVD to reduce the bond dimension of the row
				//Combines pr_index-1 and pr_index, then splits again with reduced BD
				//Does not need to do the 0,1 link as the auxiliary MPS wasn't contracted into that
				for(int pr_index = 2; pr_index < previous_row.size(); pr_index++){
					itensor::ITensor combined_link_tensor = previous_row[pr_index-1]*previous_row[pr_index];
					auto left_links = itensor::commonInds(combined_link_tensor, previous_row[pr_index-1]);
					auto [left_tensor, singular_vals, right_tensor] = itensor::svd(combined_link_tensor, left_links, {"MaxDim",_Dc});
					previous_row[pr_index-1] = left_tensor;
					previous_row[pr_index] = singular_vals*right_tensor;
				}
			}
			out.push_front(AuxMPS(AuxType::NA)); //The first column is a dummy AuxMPS that shouldn't ever have to be called
			if(_log){
				std::cout.rdbuf(coutbuf);
				log_file_stream.close();
			}
			return out;
		}

		//Get the auxiliary MPS's for the long, i.e. i+j, direction
		std::list<AuxMPS> get_ld_auxiliaries(){
			bool _log = (_log_file != "");
			std::streambuf *coutbuf = std::cout.rdbuf();
			std::ofstream log_file_stream(_log_file, std::ofstream::app);
			if(_log){
				std::cout.rdbuf(log_file_stream.rdbuf());
			}
			std::list<AuxMPS> out;
			std::vector<itensor::ITensor> unsplit_MPS; //Unsplit MPS is formed from (:,h,2) at the start of the truncation step
			std::vector<itensor::ITensor> previous_row; //Unlike VD and SD case, previous row starts off empty
			unsplit_MPS.push_back(_site_tensors[_Nx-1][_Ny-1][2]);
			for(int h = _Nx+_Ny-2; h >= 0; h--){
				std::cout << "Contracting Diagonal " << h << ":" << std::endl;
				//Split the (:,h,2) tensors with SVD
				previous_row.clear();
				int imin = std::max(0, h-_Ny+1);
				int imax = std::min(_Nx, h+1);
				AuxMPS aux(AuxType::LD);
				for(int i = imin; i < imax; i++){
					int unsplit_index = i-imin;
					int j = h-i;
					itensor::Index left_upper_link = itensor::commonIndex(unsplit_MPS[unsplit_index], _site_tensors[i][j][0]);
					itensor::ITensor sing_vals, right_tensor;
					itensor::ITensor left_tensor(left_upper_link);
					if(unsplit_index>0){
						itensor::Index left_link = itensor::commonIndex(unsplit_MPS[unsplit_index], unsplit_MPS[unsplit_index-1]);
						left_tensor = itensor::ITensor(left_upper_link, left_link);
					}
					itensor::svd(unsplit_MPS[unsplit_index], left_tensor, sing_vals, right_tensor, {"MaxDim", _Dc});
					std::cout << "SVD Error: " << itensor::sqr(itensor::norm(unsplit_MPS[unsplit_index]-left_tensor*sing_vals*right_tensor)/itensor::norm(unsplit_MPS[unsplit_index])) << std::endl;
					left_tensor *= sing_vals;
					if(_log){
						std::cout << "Split tensors at i=" << i << ", j=" << j << std::endl;
						Print(left_tensor);
						Print(right_tensor);
					}
					aux.add_tensor(left_tensor);
					aux.add_tensor(right_tensor);
				}
				out.push_front(aux);
				if(h > 0){
					//Contract the auxiliary MPS into the (:,h,0) and (:,h,1) rows
					for(int i = imin; i < imax; i++){
						int j = h-i;
						int aux_index = 2*(i-imin);
						previous_row.push_back(_site_tensors[i][j][0]*aux.MPS[aux_index]);
						previous_row.push_back(_site_tensors[i][j][1]*aux.MPS[aux_index+1]);
					}

					//Truncate the bond dimensions in the previous row. Does all links as they all need to be truncated now.
					for(int pr_index = 1; pr_index < previous_row.size(); pr_index++){
						auto combined_pr_tensor = previous_row[pr_index-1]*previous_row[pr_index];
						auto left_links = itensor::commonInds(combined_pr_tensor, previous_row[pr_index-1]);
						auto [left_tensor, singular_vals, right_tensor] = itensor::svd(combined_pr_tensor, left_links, {"MaxDim", _Dc});
						previous_row[pr_index-1] = left_tensor;
						previous_row[pr_index] = singular_vals*right_tensor;
					}

					//Contract the (:,h,0) and (:,h,1) rows into the (:,h-1,2) row
					unsplit_MPS.clear();
					/*for(int i = std::max(0, h-_Ny); i < std::min(_Nx, h); i++){
						int j = h-i-1;
						unsplit_MPS.push_back(_site_tensors[i][j][2]);
					}*/
					for(int i = imin; i < imax; i++){
						int j = h-i;
						int pr_index = 2*(i-imin);
						if(i==0){ //A (0,j,0) index doesn't have a site in the next row to contract to
							previous_row[pr_index+1] *= previous_row[pr_index];
						}
						if(j==0){ //A (i,0,1) index doesn't have a site in the next row to contract to
							previous_row[pr_index] *= previous_row[pr_index+1];
						}
						if(i>0){
							if(i == imin){//If it's not the top site, the site it needs to contract into is already in the unsplit MPS because the (i-1,j+1,1) contraction created it
								unsplit_MPS.push_back(_site_tensors[i-1][j][2]*previous_row[pr_index]);
							}
							else{
								unsplit_MPS[unsplit_MPS.size()-1] *= previous_row[pr_index];
							}
						}
						if(j>0){
							unsplit_MPS.push_back(_site_tensors[i][j-1][2]*previous_row[pr_index+1]);
						}
					}
					for(int pr_index = 0; pr_index < previous_row.size(); pr_index++){
						Print(previous_row[pr_index]);
					}
					for(int unsplit_index = 0; unsplit_index < unsplit_MPS.size(); unsplit_index++){
						Print(unsplit_MPS[unsplit_index]);
					}
				}
			}
			//Unlike the vertical and short direction cases, the first row needs to be called
			if(_log){
				std::cout.rdbuf(coutbuf);
				log_file_stream.close();
			}
			return out;
		}

		void print_self(){
			bool _log = (_log_file != "");
			std::streambuf *coutbuf = std::cout.rdbuf();
			std::ofstream log_file_stream(_log_file, std::ofstream::app);
			if(_log){
				std::cout.rdbuf(log_file_stream.rdbuf());
			}
			std::cout << "TENSOR DATA:" << std::endl;
			for(int i = 0; i < _Nx; i++){
				for(int j = 0; j < _Ny; j++){
					for(int k = 0; k < UNIT_CELL_SIZE; k++){
						std::cout << "Tensor " << i << ", " << j << ", " << k << ":" << std::endl;
						Print(_site_tensors[i][j][k]);
					}
				}
			}
			std::cout << "LINK DATA:" << std::endl;
			for(auto const& [link_key, link_index] : _link_indices){
				auto [site_i, site_j] = link_index_to_pair(link_key);
				std::cout << site_i << "-" << site_j << " link:" << std::endl;
				Print(link_index);
			}
			if(_log){
				std::cout.rdbuf(coutbuf);
				log_file_stream.close();
			}
		}
		std::tuple<int, int, int> position_of_site(int site_index){
			int k = site_index % UNIT_CELL_SIZE;
			site_index = (site_index - k)/UNIT_CELL_SIZE;
			int j = site_index % _Ny;
			site_index = (site_index - j)/_Ny;
			int i = site_index % _Nx;
			return std::make_tuple(i, j, k);
		}

		itensor::ITensor site_tensor(int i, int j, int k){
			return _site_tensors[i][j][k];
		}

	protected:
		int pair_to_link_index(int i, int j){
			if(i < j){
				return i*_num_sites + j;
			}
			else{
				return j*_num_sites + i;
			}
		}

		std::tuple<int, int> link_index_to_pair(int link_index){
			int j = link_index % _num_sites;
			int i = (link_index - j)/_num_sites;
			return std::make_tuple(i,j);
		}

		int site_index_from_position(int i, int j, int k){
			return i*_Ny*UNIT_CELL_SIZE + j*UNIT_CELL_SIZE + k;
		}

		

		int link_index_from_position(int i1, int j1, int k1, int i2, int j2, int k2){
			return pair_to_link_index(site_index_from_position(i1, j1, k1), site_index_from_position(i2, j2, k2));
		}
		int lifp(int i1, int j1, int k1, int i2, int j2, int k2){
			return link_index_from_position(i1, j1, k1, i2, j2, k2);
		}

		void create_link_index(int i1, int j1, int k1, int i2, int j2, int k2){
			std::string link_name = "Link,n1=" + std::to_string(site_index_from_position(i1, j1, k1)+1) + ",n2=" + std::to_string(site_index_from_position(i2, j2, k2)+1);
			_link_indices[lifp(i1, j1, k1, i2, j2, k2)] = itensor::Index(_D, link_name);
		}

		void create_link_indices(){
			for(int i = 0; i < _Nx; i++){
				for(int j = 0; j < _Ny; j++){
					create_link_index(i,j,0,i,j,1);
					create_link_index(i,j,0,i,j,2);
					create_link_index(i,j,1,i,j,2);
					if(j > 0){
						create_link_index(i,j-1,2,i,j,1);
						if(i <_Nx-1){
							create_link_index(i,j,1,i+1,j-1,0);
						}
					}
					if(i > 0){
						create_link_index(i-1,j,2,i,j,0);
					}
				}
			}
		}

		virtual void create_site_tensors(bool randomize = true){
			for(int i = 0; i < _Nx; i++){
				std::vector<std::vector<itensor::ITensor>> _site_tensors_2D;
				for(int j = 0; j < _Ny; j++){
					std::vector<itensor::ITensor> _site_tensors_1D;
					for(int k = 0; k < UNIT_CELL_SIZE; k++){
						int parent_index = site_index_from_position(i,j,k); 
						//std::cerr << "Creating tensor #" << parent_index << std::endl;
						std::vector<itensor::Index> indices;
						for(int other_site = 0; other_site < _num_sites; other_site++){
							auto possible_link = _link_indices.find(pair_to_link_index(other_site, parent_index));
							if(possible_link != _link_indices.end()){
								indices.push_back(possible_link->second);
							}
						}
						//indices.push_back(sites(parent_index+1));
						//std::cerr << "  Finishing tensor creation..." << std::endl;
						itensor::ITensor new_site_tensor(indices);
						//std::cerr << "  Randomizing tensor..." << std::endl;
						if(randomize){
							new_site_tensor.randomize();
							new_site_tensor /= itensor::norm(new_site_tensor);
						}
						
						//std::cerr << "  Pushing back tensor..." << std::endl;
						_site_tensors_1D.push_back(new_site_tensor);
					}
					_site_tensors_2D.push_back(_site_tensors_1D);
				}
				_site_tensors.push_back(_site_tensors_2D);
			}
		}
};

//Kagome Lattice PEPS that uses Monte Carlo to evaluate itself
class MCKPEPS : public NoSitePEPS{
	public:
		itensor::IndexSet site_indices;
		
		MCKPEPS(itensor::IndexSet &sites,
			int input_Nx,
			int input_Ny,
			int input_max_bd,
			int input_max_truncation_bd,
			itensor::Args const& args = itensor::Args::global())
		{
			_num_sites = input_Nx*input_Ny*UNIT_CELL_SIZE;
			_Nx = input_Nx;
			_Ny = input_Ny;
			bonds.set_dimensions(input_Nx, input_Ny);
			_D = input_max_bd;
			_log_file = "";
			_Dc = input_max_truncation_bd;
			site_indices = itensor::IndexSet(sites);
			bool randomize = args.getBool("RandomizeSites", true);

			//std::cerr << "Creating link indices..." << std::endl;
			create_link_indices();
			//std::cerr << "Creating site tensors..." << std::endl;
			create_site_tensors(sites, randomize);
		}

		/*MCKPEPS(MCKPEPS &to_copy){
			_Nx = to_copy._Nx;
			_Ny = to_copy._Ny;
			_num_sites = to_copy._num_sites;
			_D = to
		}*/
		
		int physical_dims(){ return site_indices(1).dim(); }

		
		//First combines the three sites in each size-3 unit cell, then combines the unit cell tensors
		//O(D^8 + D^2Ny+2)
		double brute_force_inner_product(MCKPEPS &other){
			std::cerr << "Brute force contraction..." << std::endl;
			std::vector<std::vector<itensor::ITensor>> cell_tensors;
			for(int i = 0; i < _Nx; i++){
				std::vector<itensor::ITensor> cell_tensors_in_row;
				for(int j = 0; j < _Ny; j++){
					itensor::ITensor cell_tensor(1);
					for(int k = 0; k < UNIT_CELL_SIZE; k++){
						cell_tensor *= _site_tensors[i][j][k];
						cell_tensor *= other._site_tensors[i][j][k];
					}
					std::cerr << "Combined tensors at cell " << i << ", " << j << std::endl;
					cell_tensors_in_row.push_back(cell_tensor);
				}
				cell_tensors.push_back(cell_tensors_in_row);
			}

			itensor::ITensor brute_force_combined_tensor(1);
			for(int i = 0; i < _Nx; i++){
				for(int j = 0; j < _Ny; j++){
					brute_force_combined_tensor *= cell_tensors[i][j];
					std::cerr << "Combined cell " << i << ", " << j << std::endl;
				}
			}
			Print(brute_force_combined_tensor);
			return itensor::norm(brute_force_combined_tensor);
		}

		double brute_force_inner_product_old(MCKPEPS &other){
			std::cerr << "Brute force contraction..." << std::endl;
			itensor::ITensor brute_force_combined_tensor(1);
			for(int i = 0; i < _Nx; i++){
				for(int j = 0; j < _Ny; j++){
					for(int k = 0; k < UNIT_CELL_SIZE; k++){
						brute_force_combined_tensor *= _site_tensors[i][j][k];
						brute_force_combined_tensor *= other._site_tensors[i][j][k];
						std::cerr << "Combined tensor at site " << i << ", " << j << ", " << k << std::endl;
					}
				}
			}
			Print(brute_force_combined_tensor);
			return itensor::norm(brute_force_combined_tensor);
		}

		//Combines this PEPS with other but doesn't contract the link indices, instead returning a PEPS with no site indices
		NoSitePEPS contract(MCKPEPS &other){
			bool _log = (_log_file != "");
			std::streambuf *coutbuf = std::cout.rdbuf();
			std::ofstream log_file_stream(_log_file, std::ofstream::app);
			if(_log){
				std::cout.rdbuf(log_file_stream.rdbuf());
			}

			std::vector<std::vector<std::vector<itensor::ITensor>>> combined_tensors;
			if(_log){
				std::cout << "PEPS 1: \n";
				print_self();
				std::cout << "PEPS 2: \n";
				other.print_self();
			}
			for(int i = 0; i < _Nx; i++){
				std::vector<std::vector<itensor::ITensor>> combined_tensors_rank2;
				for(int j = 0; j < _Ny; j++){
					std::vector<itensor::ITensor> combined_tensors_rank1;
					for(int k = 0; k < UNIT_CELL_SIZE; k++){
						itensor::ITensor combined_tensor = _site_tensors[i][j][k]*other._site_tensors[i][j][k];
						combined_tensors_rank1.push_back(combined_tensor);
					}
					combined_tensors_rank2.push_back(combined_tensors_rank1);
				}
				combined_tensors.push_back(combined_tensors_rank2);
			}
			if(_log){
				std::cout << "COMBINED TENSORS:" << std::endl;
				for(int i = 0; i < _Nx; i++){
					for(int j = 0; j < _Ny; j++){
						for(int k = 0; k < UNIT_CELL_SIZE; k++){
							std::cout << "Site " << i << ", " << j << ", " << k << ":" << std::endl;
							Print(combined_tensors[i][j][k]);
						}
					}
				}
			}
			//Join the link indices of the combined tensors
			for(int site_i = 0; site_i < _num_sites; site_i++){
				for(int site_j = site_i; site_j < _num_sites; site_j++){
					int link_index = pair_to_link_index(site_i, site_j);
					if(_link_indices.count(link_index) > 0){
						//For link (site_i, site_j), create a combiner tensor and apply it to both
						auto [site_i_x, site_i_y, site_i_z] = position_of_site(site_i);
						auto [site_j_x, site_j_y, site_j_z] = position_of_site(site_j);
						itensor::Index original_link_1 = _link_indices[link_index];
						itensor::Index original_link_2 = other._link_indices[link_index];
						auto [Combiner_link, combined_index] = itensor::combiner(original_link_1, original_link_2);
						combined_tensors[site_i_x][site_i_y][site_i_z] *= Combiner_link;
						combined_tensors[site_j_x][site_j_y][site_j_z] *= Combiner_link;
					}
				}
			}

			if(_log){
				std::cout << "COMBINED TENSORS AFTER LINK JOINING:" << std::endl;
				for(int i = 0; i < _Nx; i++){
					for(int j = 0; j < _Ny; j++){
						for(int k = 0; k < UNIT_CELL_SIZE; k++){
							std::cout << "Site " << i << ", " << j << ", " << k << ":" << std::endl;
							Print(combined_tensors[i][j][k]);
						}
					}
				}
			}

			//Create a link indices map for the output PEPS
			std::map<int, itensor::Index> output_link_indices;

			for(int i = 0; i < _Nx; i++){
				for(int j = 0; j < _Ny; j++){
					add_link(i,j,0,i,j,1,output_link_indices,combined_tensors);
					add_link(i,j,0,i,j,2,output_link_indices,combined_tensors);
					add_link(i,j,1,i,j,2,output_link_indices,combined_tensors);
					if(j > 0){
						add_link(i,j-1,2,i,j,1,output_link_indices,combined_tensors);
						if(i <_Nx-1){
							add_link(i,j,1,i+1,j-1,0,output_link_indices,combined_tensors);
						}
					}
					if(i > 0){
						add_link(i-1,j,2,i,j,0,output_link_indices,combined_tensors);
					}
				}
			}

			auto nsp = NoSitePEPS(_Nx, _Ny, _D, _Dc, combined_tensors, output_link_indices);
			nsp.set_log_file(_log_file);

			if(_log){
				std::cout.rdbuf(coutbuf);
				log_file_stream.close();
			}
			return nsp;
		}

		double inner_product(MCKPEPS &other){
			bool _log = (_log_file != "");
			std::streambuf *coutbuf = std::cout.rdbuf();
			std::ofstream log_file_stream(_log_file, std::ofstream::app);
			if(_log){
				std::cout.rdbuf(log_file_stream.rdbuf());
			}
			std::cout << "Checking dimensions...\n";
			if((_Nx != other._Nx) || (_Ny != other._Ny)){
				std::cout << "ERROR: Can't contract PEPS of incompatible dimensions" << std::endl;
			}
			std::vector<std::vector<std::vector<itensor::ITensor>>> combined_tensors;
			if(_log){
				std::cout << "PEPS 1: \n";
				print_self();
				std::cout << "PEPS 2: \n";
				other.print_self();
			}
			for(int i = 0; i < _Nx; i++){
				std::vector<std::vector<itensor::ITensor>> combined_tensors_rank2;
				for(int j = 0; j < _Ny; j++){
					std::vector<itensor::ITensor> combined_tensors_rank1;
					for(int k = 0; k < UNIT_CELL_SIZE; k++){
						itensor::ITensor combined_tensor = _site_tensors[i][j][k]*other._site_tensors[i][j][k];
						combined_tensors_rank1.push_back(combined_tensor);
					}
					combined_tensors_rank2.push_back(combined_tensors_rank1);
				}
				combined_tensors.push_back(combined_tensors_rank2);
			}
			if(_log){
				std::cout << "COMBINED TENSORS:" << std::endl;
				for(int i = 0; i < _Nx; i++){
					for(int j = 0; j < _Ny; j++){
						for(int k = 0; k < UNIT_CELL_SIZE; k++){
							std::cout << "Site " << i << ", " << j << ", " << k << ":" << std::endl;
							Print(combined_tensors[i][j][k]);
						}
					}
				}
			}
			//Join the link indices of the combined tensors
			for(int site_i = 0; site_i < _num_sites; site_i++){
				for(int site_j = site_i; site_j < _num_sites; site_j++){
					int link_index = pair_to_link_index(site_i, site_j);
					if(_link_indices.count(link_index) > 0){
						//For link (site_i, site_j), create a combiner tensor and apply it to both
						auto [site_i_x, site_i_y, site_i_z] = position_of_site(site_i);
						auto [site_j_x, site_j_y, site_j_z] = position_of_site(site_j);
						itensor::Index original_link_1 = _link_indices[link_index];
						itensor::Index original_link_2 = other._link_indices[link_index];
						auto [Combiner_link, combined_index] = itensor::combiner(original_link_1, original_link_2);
						combined_tensors[site_i_x][site_i_y][site_i_z] *= Combiner_link;
						combined_tensors[site_j_x][site_j_y][site_j_z] *= Combiner_link;
					}
				}
			}

			if(_log){
				std::cout << "COMBINED TENSORS AFTER LINK JOINING:" << std::endl;
				for(int i = 0; i < _Nx; i++){
					for(int j = 0; j < _Ny; j++){
						for(int k = 0; k < UNIT_CELL_SIZE; k++){
							std::cout << "Site " << i << ", " << j << ", " << k << ":" << std::endl;
							Print(combined_tensors[i][j][k]);
						}
					}
				}
			}

			//Contract the resulting KPEPS

			for(int i = _Nx-1; i > 0; i--){
				//First layer: Contract the bottom row at k=1 and 2 
				//First contract k=1 with k=2, then contract that with k=0
				std::vector<itensor::ITensor> row_tensors;
				for(int j = 0; j < _Ny; j++){
					itensor::ITensor intermediate_tensor = combined_tensors[i][j][1]*combined_tensors[i][j][2];
					row_tensors.push_back(combined_tensors[i][j][0]*intermediate_tensor);
				}

				if(_log){
					std::cout << "ROW TENSOR: " << std::endl;
					for(int row_tensor_index = 0; row_tensor_index < row_tensors.size(); row_tensor_index++){
						std::cout << "Site " << row_tensor_index << std::endl;
						Print(row_tensors[row_tensor_index]);
					}
				}

				//Second layer: Split each row tensor into two copies using truncated SVD, then combine those copies with the tensors above
				//Need to include j=0 case later
				for(int j = 0; j < _Ny-1; j++){
					itensor::Index left_upper_link = itensor::commonIndex(row_tensors[j], combined_tensors[i-1][j][2]);
					itensor::ITensor sing_vals, right_tensor;
					itensor::ITensor left_tensor(left_upper_link);

					if(j>0){
						itensor::Index left_link = itensor::commonIndex(row_tensors[j], row_tensors[j-1]);
						left_tensor = itensor::ITensor(left_upper_link, left_link);
					}
					
					itensor::svd(row_tensors[j], left_tensor, sing_vals, right_tensor, {"MaxDim", _Dc});
					std::cout << "SVD Error: " << itensor::sqr(itensor::norm(row_tensors[j]-left_tensor*sing_vals*right_tensor)/itensor::norm(row_tensors[j])) << std::endl;
					
					left_tensor *= sing_vals;

					if(_log){
						std::cout << "Split tensors at i=" << i << ", j=" << j << std::endl;
						Print(left_tensor);
						Print(right_tensor);
					}
					combined_tensors[i-1][j][2] *= left_tensor;

					combined_tensors[i-1][j+1][1] *= right_tensor;
				}
				combined_tensors[i-1][_Ny-1][2] *= row_tensors[_Ny-1];

				if(_log){
					std::cout << "UNTRUNCATED LAYER " << i-1 << ": " << std::endl;
					for(int j = 0; j < _Ny; j++){
						std::cout << "j=" << j << std::endl;
						Print(combined_tensors[i-1][j][1]);
						Print(combined_tensors[i-1][j][2]);
					}
				}

				//Truncate indices
				for(int j = 0; j < _Ny; j++){
					itensor::ITensor left_tensor, sing_vals, right_tensor;
					//Need to truncate each link index left to right to make sure we're never truncating a 
					//tensor with two D*Dc links
					//Left_tensor keeps the link index further in the PEPS and becomes the new site tensor. 
					//First truncate the k=1 site (j=0 case special)
					auto link_up = itensor::commonIndex(combined_tensors[i-1][j][1], combined_tensors[i-1][j][0]);
					if(j == 0){
						left_tensor = itensor::ITensor(link_up);
					}
					else{
						auto links_left = itensor::commonInds(combined_tensors[i-1][j][1], combined_tensors[i-1][j-1][2]);
						left_tensor = itensor::ITensor(itensor::unionInds(links_left, link_up));
					}
					itensor::svd(combined_tensors[i-1][j][1], left_tensor, sing_vals, right_tensor, {"MaxDim", _Dc});
					std::cout << "SVD Error: " << itensor::sqr(itensor::norm(combined_tensors[i-1][j][1]-left_tensor*sing_vals*right_tensor)/itensor::norm(combined_tensors[i-1][j][1])) << std::endl;
					
					combined_tensors[i-1][j][1] = left_tensor;
					combined_tensors[i-1][j][2] *= (sing_vals*right_tensor);
					//Now truncate the k=2 site (j=0 case not special) (No truncation for j=Ny-1 case necessary)
					if(j != _Ny-1){
						link_up = itensor::commonIndex(combined_tensors[i-1][j][2], combined_tensors[i-1][j][0]);
						auto links_left_2 = itensor::commonIndex(combined_tensors[i-1][j][1], combined_tensors[i-1][j][2]);
						left_tensor = itensor::ITensor(link_up, links_left_2);
						//TODO: Check that this doesn't modify the data on combined_tensors[i-1][j][1/2]?
						itensor::svd(combined_tensors[i-1][j][2], left_tensor, sing_vals, right_tensor, {"MaxDim", _Dc});
						std::cout << "SVD Error: " << itensor::sqr(itensor::norm(combined_tensors[i-1][j][2]-left_tensor*sing_vals*right_tensor)/itensor::norm(combined_tensors[i-1][j][2])) << std::endl;
						combined_tensors[i-1][j][2] = left_tensor;
						combined_tensors[i-1][j+1][1] *= (sing_vals*right_tensor);
					}

				}
				if(_log){
					std::cout << "TRUNCATED LAYER " << i-1 << ": " << std::endl;
					for(int j = 0; j < _Ny; j++){
						std::cout << "j=" << j << std::endl;
						Print(combined_tensors[i-1][j][1]);
						Print(combined_tensors[i-1][j][2]);
					}
				}
			}
			itensor::ITensor contracted_tensor(1);
			//contracted_tensor *= combined_tensors[0][0][1];
			//contracted_tensor *= combined_tensors[0][0][2];
			//Contract the i=0 layer
			for(int j = 0; j < _Ny; j++){
				for(int k = 0; k < UNIT_CELL_SIZE; k++){
					contracted_tensor *= combined_tensors[0][j][k];
				}
			}
			Print(contracted_tensor);

			if(_log){
				std::cout.rdbuf(coutbuf);
				log_file_stream.close();
			}
			return itensor::norm(contracted_tensor);

		}

		


		//Applies a set of simple spin operators (ITensors with only site indices) to the PEPS
		//Assumes one site matches the site tensors, and the other is the primed version
		void apply_spinop(std::vector<int> sites, std::vector<itensor::ITensor> spinops){
			for(int site_index = 0; site_index < sites.size(); site_index ++){
				int site = sites[site_index];
				auto position = position_of_site(site);
				int i = std::get<0>(position);
				int j = std::get<1>(position);
				int k = std::get<2>(position);
				_site_tensors[i][j][k] *= spinops[site_index];
				_site_tensors[i][j][k].noPrime("Site");
			}
		}
		void apply_spinop(int site, itensor::ITensor spinops){
			auto position = position_of_site(site);
			int i = std::get<0>(position);
			int j = std::get<1>(position);
			int k = std::get<2>(position);
			_site_tensors[i][j][k] *= spinops;
			_site_tensors[i][j][k].noPrime("Site");
		}


		MCKPEPS get_spin_config(){//Get a random spin configuration with the same sites and dimensions
			return MCKPEPS(site_indices, _Nx, _Ny, 1, _Dc);
		}

		
	
	protected:
		void create_site_tensors(itensor::IndexSet &sites, bool randomize = true){
			for(int i = 0; i < _Nx; i++){
				std::vector<std::vector<itensor::ITensor>> _site_tensors_2D;
				for(int j = 0; j < _Ny; j++){
					std::vector<itensor::ITensor> _site_tensors_1D;
					for(int k = 0; k < UNIT_CELL_SIZE; k++){
						int parent_index = site_index_from_position(i,j,k); 
						//std::cerr << "Creating tensor #" << parent_index << std::endl;
						std::vector<itensor::Index> indices;
						for(int other_site = 0; other_site < _num_sites; other_site++){
							auto possible_link = _link_indices.find(pair_to_link_index(other_site, parent_index));
							if(possible_link != _link_indices.end()){
								//std::cerr << "Found link at " << other_site << ", " << parent_index << std::endl;
								indices.push_back(possible_link->second);
							}
						}
						indices.push_back(sites(parent_index+1));
						//std::cerr << "  Finishing tensor creation..." << std::endl;
						itensor::ITensor new_site_tensor(indices);
						//std::cerr << "  Randomizing tensor..." << std::endl;
						if(randomize){
							new_site_tensor.randomize();
							new_site_tensor /= itensor::norm(new_site_tensor);
						}
						
						//std::cerr << "  Pushing back tensor..." << std::endl;
						_site_tensors_1D.push_back(new_site_tensor);
					}
					_site_tensors_2D.push_back(_site_tensors_1D);
				}
				_site_tensors.push_back(_site_tensors_2D);
			}
		}
		void add_link(int i1, int j1, int k1, int i2, int j2, int k2, 
			std::map<int, itensor::Index> &links, 
			std::vector<std::vector<std::vector<itensor::ITensor>>> &tns){
			links[lifp(i1, j1, k1, i2, j2, k2)] = itensor::commonIndex(tns[i1][j1][k1], tns[i2][j2][k2]);
		}
};



class SpinConfigPEPS : public MCKPEPS
{
	public:
		SpinConfigPEPS(itensor::IndexSet &sites,
			int input_Nx,
			int input_Ny) : MCKPEPS(sites, input_Nx, input_Ny, 1, 1, {"RandomizeSites",false}){
		}

		SpinConfigPEPS(MCKPEPS &base_state) : MCKPEPS(base_state.site_indices, base_state.Nx(), base_state.Ny(), 1, 1, {"RandomizeSites",false}){
		}

		void set_spin(int i , int j, int k, int spin_value, double wavefunction_normalization = 1){
			int site_index = site_index_from_position(i,j,k);
			_site_tensors[i][j][k] = itensor::ITensor(_site_tensors[i][j][k].inds());
			std::vector<itensor::IndexVal> config;
			for(itensor::Index ind : _site_tensors[i][j][k].inds()){
				if(ind.dim() > 1){
					config.push_back(ind(spin_value+1));
				}
				else{
					config.push_back(ind(1));
				}
			}
			_site_tensors[i][j][k].set(config, 1.0/wavefunction_normalization);
		}

		void set_spin(int site_index, int spin_value, double wavefunction_normalization = 1){
			auto ijk = position_of_site(site_index);
			set_spin(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk),spin_value, wavefunction_normalization);
		}

		void set_spins(std::vector<int> &spin_config, double wavefunction_normalization = 1){
			for(int spin_index = 0; spin_index < spin_config.size(); spin_index ++){
				set_spin(spin_index, spin_config[spin_index], wavefunction_normalization);
			}
		}


};

#endif