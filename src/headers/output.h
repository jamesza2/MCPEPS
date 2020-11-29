#ifndef OUTPUT_MAKER
#define OUTPUT_MAKER
#include <vector>
#include <string>
#include <iostream>

class Output
{
	public:
		Output(){ }
		/*Output(std::string file_name){
			_file_name = file_name;
		}*/
		/*void setFileName(std::string new_file_name){
			_file_name = new_file_name;
		}*/

		void addInteger(std::string name, int out_int){
			_output += format_name(name) + "\n" + std::to_string(out_int);
		}
		void addDouble(std::string name, double out_double){
			_output += format_name(name) + "\n" + std::to_string(out_double);
		}
		void addString(std::string name, std::string out_string){
			_output += format_name(name) + "\n" + out_string;
		}
		template<typename T>
		void addVar(std::string name, T var){
			_output += format_name(name) + "\n" + std::to_string(var);
		}
		void addVector(std::string name, std::vector<T> &out_vector){
			_output += format_name(name) + "\n";
			for(T elem : out_vector){
				_output += std::to_string(elem) + " ";
			}
		}
		template<typename T>
		void addVectorVector(std::string name, std::vector<std::vector<T>> &out_vv){
			_output += format_name(name);
			for(int i = 0; i < out_vv.size(); i++){
				_output += "\n";
				for(int j = 0; j < out_vv[i].size(); j++){
					_output += std::to_string(out_vv[i][j]) + " ";
				}
			}
		}

		void writeOutput(std::string file_name){
			std::ofstream out_file(file_name);
			out_file << _output;
			out_file.close();
		}

	private:
		//std::string _file_name;
		std::string _output;
		bool first_line = true;
		std::string format_name(std::string name){
			if(!first_line){
				return "\n" + name + ": ";
			}
			else{
				first_line = false;
				return name + ": ";
			}
		}

};

#endif