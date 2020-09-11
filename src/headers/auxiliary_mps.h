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
		AuxType type;
		std::vector<itensor::ITensor> MPS;
		int length;
		AuxMPS(){
			type=AuxType::NA;
			length = 0;
		}
		AuxMPS(AuxType input_type){
			type = input_type;
			length = 0;
		}
		AuxMPS(int input_length){
			type=AuxType::BLANK;
			length = input_length;
			for(int i = 0; i < length; i++){
				MPS.push_back(itensor::ITensor(1));
			}
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
		}

		std::vector<itensor::ITensor>::iterator begin(){return MPS.begin();}
		std::vector<itensor::ITensor>::iterator end(){return MPS.end();}


};


#endif