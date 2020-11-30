#ifndef PEPS_OPERATOR
#define PEPS_OPERATOR

#include "mcpeps.h"
#include "operator.h"

enum class OpType{I,SP,SM,SZ,SX,SZ2};

//Examples include Sx, S+, Sz, etc.
//Can return the ITensor form of the operator, given input physical indices
class Spinop{
	protected:
		OpType type;
		itensor::ITensor identity(itensor::Index &ind1, itensor::Index &ind2){
			itensor::ITensor T(ind1, ind2);
			for(int i = 1; i <= itensor::dim(ind1); i++){
				T.set(ind1 = i, ind2 = i, 1.);
			}
			return T;
		}
		//Upgrades spin from ind1 to ind2
		itensor::ITensor spinplus(itensor::Index &ind1, itensor::Index &ind2){
			itensor::ITensor T(ind1, ind2);
			double s = 0.5*(itensor::dim(ind1)-1);
			for(int i = 1; i <= itensor::dim(ind1)-1; i++){
				double sz = i-s-1;
				double me = std::sqrt(s*(s+1)-sz*(sz+1));
				T.set(ind1 = i, ind2 = i+1, me);
			}
			return T;
		}
		//Downgrades spin from ind1 to ind2
		itensor::ITensor spinminus(itensor::Index &ind1, itensor::Index &ind2){
			itensor::ITensor T(ind1, ind2);
			double s = 0.5*(itensor::dim(ind1)-1);
			for(int i = 2; i <= itensor::dim(ind1); i++){
				double sz = i-s-1;
				double me = std::sqrt(s*(s+1)-sz*(sz-1));
				T.set(ind1 = i, ind2 = i-1, me);
			}
			return T;
		}
		itensor::ITensor spinz(itensor::Index &ind1, itensor::Index &ind2){
			itensor::ITensor T(ind1, ind2);
			for(int i = 1; i <= itensor::dim(ind1); i++){
				double sz = i-0.5*itensor::dim(ind1)-0.5;
				T.set(ind1 = i, ind2 = i, sz);
			}
			return T;
		}
		itensor::ITensor spinzsquared(itensor::Index &ind1, itensor::Index &ind2){
			itensor::ITensor T(ind1, ind2);
			for(int i = 1; i <= itensor::dim(ind1); i++){
				double sz = i-0.5*itensor::dim(ind1)-0.5;
				T.set(ind1 = i, ind2 = i, sz*sz);
			}
			return T;
		}

	public:
		Spinop(){type = OpType::I;}
		Spinop(OpType t){type = t;}

		itensor::ITensor tensor(itensor::Index &ind1){
			itensor::Index ind1p = itensor::prime(ind1);
			return tensor(ind1, ind1p);
		}
		itensor::ITensor tensor(itensor::Index &ind1, itensor::Index &ind2){
			if(itensor::dim(ind1) != itensor::dim(ind2)){
				std::cerr << "Error: Index lengths " << itensor::dim(ind1) << " and " << itensor::dim(ind2) << " do not match" << std::endl;
			}
			switch(type){
				case OpType::I : return identity(ind1, ind2); break;
				case OpType::SP: return spinplus(ind1, ind2); break;
				case OpType::SM: return spinminus(ind1, ind2); break;
				case OpType::SZ: return spinz(ind1, ind2); break;
				case OpType::SZ2: return spinzsquared(ind1, ind2); break;
				default: return itensor::ITensor(ind1, ind2);
			}
		}
		std::string to_string(){
			switch(type){
				case OpType::I : return "I"; break;
				case OpType::SP: return "S+"; break;
				case OpType::SM: return "S-"; break;
				case OpType::SZ: return "Sz"; break;
				case OpType::SZ2: return "Sz^2"; break;
				default: return "ERR";
			}
		}

};

//A product of spinops, e.g. 0.5*S+S- or SzSz
class Term{
	protected:
		std::list<int> sites;
		std::list<Spinop> ops;
	public:
		double factor;
		int num_sites;
		Term(){
			factor = 1;
			num_sites = 0;
		}
		Term(double f){
			factor = f;
			num_sites = 0;
		}
		void extra_factor(double ef){
			factor *= ef;
		}
		void add_site(int site_number, Spinop op){
			sites.push_back(site_number);
			ops.push_back(op);
			num_sites ++;
		}

		void add_site(int site_number, OpType type){
			Spinop s(type);
			add_site(site_number, s);
		}

		double eval(MCKPEPS &PEPS1, MCKPEPS &PEPS2){
			if(factor == 0){return 0;}
			auto sites_it = sites.begin();
			auto ops_it = ops.begin();
			std::list<itensor::ITensor> old_peps_tensors;
			for(int op_index = 0; op_index < num_sites; op_index++){
				itensor::ITensor T = ops_it->tensor(PEPS1.site_indices[*sites_it]);
				auto [i,j,k] = PEPS1.position_of_site(*sites_it);
				old_peps_tensors.push_back(itensor::ITensor(PEPS1._site_tensors[i][j][k]));
				//Print(PEPS1._site_tensors[i][j][k]);
				//Apply the op to PEPS1 
				PEPS1._site_tensors[i][j][k] *= T;
				PEPS1._site_tensors[i][j][k] *= itensor::delta(PEPS1.site_indices[*sites_it], itensor::prime(PEPS1.site_indices[*sites_it]));
				sites_it++;
				ops_it++;
				//Print(PEPS1._site_tensors[i][j][k]);
			}
			double matrix_element = PEPS1.inner_product(PEPS2)*factor;
			//Restore the original tensors
			sites_it = sites.begin();
			for(auto old_peps_it = old_peps_tensors.begin(); old_peps_it != old_peps_tensors.end(); old_peps_it++){
				auto [i,j,k] = PEPS1.position_of_site(*sites_it);
				PEPS1._site_tensors[i][j][k] = *old_peps_it;
				sites_it++;
				//Print(PEPS1._site_tensors[i][j][k]);
			}
			return matrix_element;
		}

		void apply(MCKPEPS &PEPS1){
			auto sites_it = sites.begin();
			PEPS1 *= factor;
			for(auto ops_it = ops.begin(); ops_it != ops.end(); ops_it++){
				itensor::ITensor T = ops_it->tensor(PEPS1.site_indices[*sites_it]);
				auto [i,j,k] = PEPS1.position_of_site(*sites_it);
				//std::cerr << "Original tensor: " << std::endl;
				//PrintData(PEPS1._site_tensors[i][j][k]);
				PEPS1._site_tensors[i][j][k] *= T;
				PEPS1._site_tensors[i][j][k] *= itensor::delta(PEPS1.site_indices[*sites_it], itensor::prime(PEPS1.site_indices[*sites_it]));
				//std::cerr << "Tensor after " << to_string() << std::endl;
				//PrintData(PEPS1._site_tensors[i][j][k]);
				sites_it++;
			}
			
		}
		std::string to_string(){
			auto sites_it = sites.begin(); 
			auto ops_it = ops.begin();
			std::string to_return = "[";
			while(ops_it != ops.end()){
				to_return += ops_it->to_string() + "_" + std::to_string(*sites_it) + " ";
				ops_it++;
				sites_it++;
			}
			to_return += std::to_string(factor);
			return to_return+"]";
		}
};

//A PEPS operator that takes in single-site and two-site terms
//Can evaluate using brute-force methods
class PEPSop{
	public:
		std::vector<Term> terms;
		PEPSop(){}
		void add_term(Term t){
			terms.push_back(t);
		}
		//Add the S+S- and S-S+ terms (not including the extra 0.5 factor in front)
		void add_spm(int site_1, int site_2, double factor){
			Term t1(factor);
			t1.add_site(site_1, OpType::SP);
			t1.add_site(site_2, OpType::SM);
			terms.push_back(t1);
			Term t2(factor);
			t2.add_site(site_1, OpType::SM);
			t2.add_site(site_2, OpType::SP);
			terms.push_back(t2);
		}
		//Add the SzSz term
		void add_szz(int site_1, int site_2, double factor){
			Term t1(factor);
			t1.add_site(site_1, OpType::SZ);
			t1.add_site(site_2, OpType::SZ);
			terms.push_back(t1);
		}

		void add_sz(int site, double factor){
			Term t1(factor);
			t1.add_site(site, OpType::SZ);
			terms.push_back(t1);
		}
		//Evaluates using a brute force method (i.e. taking the sum of <Psi|Hi|Psi> for all the terms Hi in the operator)
		double eval(MCKPEPS &PEPS1, MCKPEPS &PEPS2){
			double result = 0;
			for(int term_index = 0; term_index < terms.size(); term_index++){
				double term_ip = terms[term_index].eval(PEPS1, PEPS2);
				result += term_ip;
				//std::cerr << "Term value: " << term_ip << std::endl;
			}
			return result;
		}
};

PEPSop singleSiteSz(int site){
	PEPSop single_site_term;
	single_site_term.add_sz(site, 1);
	return single_site_term;
}

PEPSop Heisenberg::toPEPSop() const{
	PEPSop pop;
	std::vector<double> J{0, _J1, _J2, _Jd};
	for(int site_1 = 0; site_1 < _num_sites; site_1 ++){
		for(auto bond : bonds.at(site_1)){
			int site_2 = bond.first;
			if(J[bond.second] != 0){
				pop.add_spm(site_1, site_2, J[bond.second]*0.5);
				if(_Jz != 0){
					pop.add_szz(site_1, site_2, J[bond.second]*_Jz);
				}
			}
			
		}
	}
	return pop;
}

#endif