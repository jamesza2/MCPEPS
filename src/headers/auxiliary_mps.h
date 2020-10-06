#ifndef AUX_MPS
#define AUX_MPS

#include "itensor/all.h"

//V = vertical
//S = short direction (increasing/decreasing j)
//L = long direction
//U = up
//D = down
enum class AuxType{NA, VU, VD, SU, SD, LU, LD, BLANK}; 

class AuxMPS{
	public:
		int position;
		AuxType type;
		std::vector<itensor::ITensor> MPS;
		int length;
		AuxMPS(){
			type=AuxType::NA;
			length = 0;
			position = -1;
		}
		AuxMPS(AuxType input_type){
			type = input_type;
			length = 0;
			position = -1;
		}
		AuxMPS(int input_length){
			type=AuxType::BLANK;
			length = input_length;
			for(int i = 0; i < length; i++){
				itensor::ITensor blank(1);
				MPS.push_back(blank);
			}
		}
		AuxMPS(const AuxMPS &to_copy){
			MPS = std::vector<itensor::ITensor>(to_copy.MPS);
			length = to_copy.length;
			type = to_copy.type;
			position = to_copy.position;
		}
		void set_type(AuxType proposed_type){
			type = proposed_type;
		}
		void set_MPS(std::vector<itensor::ITensor> &proposed_MPS){
			MPS = proposed_MPS;
			length = MPS.size();
		}
		void add_tensor(itensor::ITensor &proposed_tensor){
			MPS.push_back(proposed_tensor);
			length += 1;
		}

		void clear(){
			MPS.clear();
			length = 0;
			position = -1;
		}

		//Canonizes to the desired position
		//Tensors might have multiple sets of links with each other, that's fine
		void set_position(int desired_pos){
			std::cerr << "In position..." << std::endl;
			if(position == -1){
				for(int i = 0; i < desired_pos; i++){
					itensor::IndexSet forward_indices = itensor::commonInds(MPS[i], MPS[i+1]);
					if(itensor::length(forward_indices) == 0){continue;}
					if(itensor::length(forward_indices) == itensor::length(MPS[i].inds())){
						MPS[i+1] *= MPS[i];
						MPS[i] = itensor::ITensor(1);
						continue;
					}
					Print(MPS[i]);
					Print(MPS[i+1]);
					auto [forward, diag, back] = itensor::svd(MPS[i], forward_indices);
					MPS[i] = back;
					MPS[i+1] *= (forward*diag);
					//Perform SVD on tensor i
					//Push S*Vt onto the tensor above
				}
				for(int i = length-1; i > desired_pos; i--){
					itensor::IndexSet back_indices = itensor::commonInds(MPS[i], MPS[i-1]);
					if(itensor::length(back_indices) == 0){continue;}
					if(itensor::length(back_indices) == itensor::length(MPS[i].inds())){
						MPS[i-1] *= MPS[i];
						MPS[i] = itensor::ITensor(1);
						continue;
					}
					Print(MPS[i]);
					Print(MPS[i-1]);
					auto [back, diag, forward] = itensor::svd(MPS[i], back_indices);
					MPS[i] = forward;
					MPS[i-1] *= (back*diag);
					//Perform SVD on tensor i
					//Push U onto the tensor below
				}
			}
			else{
				if(position < desired_pos){
					for(int i = position; i < desired_pos; i++){
						itensor::IndexSet forward_indices = itensor::commonInds(MPS[i], MPS[i+1]);
						if(itensor::length(forward_indices) == 0){continue;}
						if(itensor::length(forward_indices) == itensor::length(MPS[i].inds())){
							MPS[i+1] *= MPS[i];
							MPS[i] = itensor::ITensor(1);
							continue;
						}
						auto [forward, diag, back] = itensor::svd(MPS[i], forward_indices);
						MPS[i] = back;
						MPS[i+1] *= (forward*diag);
					}
				}
				if(position > desired_pos){
					for(int i = position; i > desired_pos; i--){
						itensor::IndexSet back_indices = itensor::commonInds(MPS[i], MPS[i-1]);
						if(itensor::length(back_indices) == 0){continue;}
						if(itensor::length(back_indices) == itensor::length(MPS[i].inds())){
							MPS[i-1] *= MPS[i];
							MPS[i] = itensor::ITensor(1);
							continue;
						}
						auto [back, diag, forward] = itensor::svd(MPS[i], back_indices);
						MPS[i] = forward;
						MPS[i-1] *= (back*diag);
					}
				}
			}
			position = desired_pos;
		}

		void truncate(int truncation_index){
			std::cerr << "In truncate..." << std::endl;
			for(int i = 0; i < length-1; i++){
				set_position(i);
				std::cerr << "In truncate..." << std::endl;
				itensor::IndexSet forward_indices = itensor::commonInds(MPS[i], MPS[i+1]);
				if(itensor::length(forward_indices) == 0){continue;}
				if(itensor::length(forward_indices) == itensor::length(MPS[i].inds())){ //If MPS[i] is just an orphaned site, contract it into MPS[i+1]
					MPS[i+1] *= MPS[i];
					MPS[i] = itensor::ITensor(1);
					continue;
				}
				Print(MPS[i]);
				Print(MPS[i+1]);
				auto [forward, diag, back] = itensor::svd(MPS[i], forward_indices, {"MaxDim", truncation_index});
				MPS[i] = back;
				MPS[i+1] *= (forward*diag);
			}
		}

		void print_self(){
			for(int i = 0; i < length; i++){
				std::cerr << "Site " << i << ": ";
				Print(MPS.at(i));
			}
		}

		std::vector<itensor::ITensor>::iterator begin(){return MPS.begin();}
		std::vector<itensor::ITensor>::iterator end(){return MPS.end();}


};


#endif