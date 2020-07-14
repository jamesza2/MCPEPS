#ifndef input_h
#define input_h
#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iostream>
using namespace std;

class InputTree
{
 public:
  string name;
  vector<InputTree*> sections;
  map<string,string> input;
  InputTree* parent;

};


class InputClass
{
 public:
  int getInteger(string myVar){
    return toInteger(GetVariable(myVar));
  }
  double getDouble(string myVar){
    return toDouble(GetVariable(myVar));
  }
  bool getBool(string myVar){
    return toBool(GetVariable(myVar));
  }

  vector<double> getVectorDouble(string myVar, char delimiter = ','){
    return toVectorDouble(GetVariable(myVar), delimiter);
  }
  vector<string> getVectorString(string myVar, char delimiter = ','){
    return toVectorString(GetVariable(myVar), delimiter);
  }
  vector<int> getVectorInteger(string myVar, char delimiter = ','){
    return toVectorInt(GetVariable(myVar), delimiter);
  }

  int toInteger(string myString)
  {
    return atoi(myString.c_str());
  }
  double toDouble(string myString)
  {
    return atof(myString.c_str());
  }

  double testDouble(string myString, double def = 0){
    if(IsVariable(myString)){
      return toDouble(GetVariable(myString));
    }
    else{
      return def;
    }
  }

  int testInteger(string myString, int def = 0){
    if(IsVariable(myString)){
      return toInteger(GetVariable(myString));
    }
    else{
      return def;
    }
  }

  bool testBool(string myString, bool def = false){
    if(IsVariable(myString)){
      return toBool(GetVariable(myString));
    }
    else{
      return def;
    }
  }

  string testString(string myString, string def = ""){
    if(IsVariable(myString)){
      return GetVariable(myString);
    }
    else{
      return def;
    }
  }

  vector<int> testVectorInt(string myString, vector<int> def){
    if(IsVariable(myString)){
      return toVectorInt(GetVariable(myString));
    }
    else{
      return def;
    }
  }

  bool toBool(string myString)
  {
    if ((myString=="True") || (myString=="true"))
      return true;
    else if ((myString=="False") || (myString=="false"))
      return false;
    cerr<<"You have a boolean which is neither true nor false"<<endl;
    exit(1);
      
  }
  vector<double> toVectorDouble(string myString, char delimiter = ',')
  {
    istringstream iss (myString);
    //vector<double> retval((istream_iterator<double>(iss)), istream_iterator<double>());
    vector<double> retval;
    string token;
    //cout << "Parsing " << myString << "..." << endl;
    while(getline(iss, token, delimiter)){
      //cout << "Found token " << token << endl;
      retval.push_back(atof(token.c_str()));
    }
    return retval;
  }
  vector<string> toVectorString(string myString, char delimiter = ','){
    istringstream iss (myString);
    vector<string> retval;
    string token;
    //cout << "Parsing " << myString << "..." << endl;
    while(getline(iss, token, delimiter)){
      //cout << "Found token " << token << endl;
      retval.push_back(token);
    }
    return retval;
  }
  vector<int> toVectorInt(string myString, char delimiter = ','){
    istringstream iss (myString);
    vector<int> retval;
    string token;
    //cout << "Parsing " << myString << "..." << endl;
    while(getline(iss, token, delimiter)){
      //cout << "Found token " << token << endl;
      retval.push_back(atoi(token.c_str()));
    }
    return retval;
  }
  bool IsVariable(string myVar)
    {
      return tree->input.count(myVar)==1;
    }
  string GetVariable(string myVar)
  {
    if (!IsVariable(myVar)){
      cerr<<"Missing variable "<<myVar<<endl;
    }
    assert(IsVariable(myVar));
    return tree->input[myVar];
  }

  bool OpenSection(string myVar,int num)
  {
    int foundSection=0;
    for (int i=0;i<tree->sections.size();i++){
      if (tree->sections[i]->name==myVar && foundSection==num){
	tree=tree->sections[i];
	return true;
      }
      else if (tree->sections[i]->name==myVar){
	foundSection++;
      }
    }
    return false;
  }

  
  bool OpenSection(string myVar)
  {
    for (int i=0;i<tree->sections.size();i++){
      if (tree->sections[i]->name==myVar){
	tree=tree->sections[i];
	return true;
      }
    }
    return false;
  }

  void CloseSection()
  {
    cerr<<"Closing the section"<<endl;
    cerr<<"The tree is currently "<<tree->name<<" "<<tree->parent->name<<endl;
    tree=tree->parent;
    cerr<<"Now the tree is  "<<tree->name<<endl;
    return;
  }
  InputTree *tree;
  void RemoveTrailingWhiteSpace(string &s)
  {
    int lastPos=s.find_last_not_of(" \t\n\v\f\r");
    if (lastPos!=string::npos){
      s=s.substr(0,lastPos+1);
    }
  }

  void RemoveLeadingWhiteSpace(string &s)
  {
    int firstPos=s.find_first_not_of(" \t\n\v\f\r");
    if (firstPos!=string::npos)
      s=s.substr(firstPos,s.size());
  }
  void RemoveWhiteSpace(string &s){
    RemoveTrailingWhiteSpace(s);
    RemoveLeadingWhiteSpace(s);
  }
  void Read(ifstream &infile)
  { 
    cerr<<"Reading input"<<endl;
    tree=new InputTree();
    tree->name="root";
    InputTree *currentTree;
    currentTree=tree;

    while (!infile.eof()){
      string line;
      infile>>line;
      RemoveWhiteSpace(line);
      if (line[line.size()-1]==':'){
      	cerr<<"Reading section "<<endl;
      	line=line.substr(0,line.size()-1);
      	InputTree *t=new InputTree();
      	t->name=line;
      	cerr<<"Pushed the tree "<<line<<endl;
      	t->parent=currentTree;
      	currentTree->sections.push_back(t);
      	currentTree=t;
      	cerr<<"Done reading section"<<endl;
      }
      else if (line=="END"){
      	cerr<<"ending section "<<endl;
      	currentTree=currentTree->parent;
      	cerr<<"Done ending section"<<endl;
      }
      else if (line.size()==0){
      }
      else {
      	int eqLoc=line.find('=');
        if(eqLoc==string::npos){
          cerr << "ERROR: Equals sign not found in line:" << endl;
          cerr << line << endl;
          throw;
        }
        assert(eqLoc!=string::npos);
      	string token=line.substr(0,eqLoc);
      	RemoveWhiteSpace(token);
      	string outVals=line.substr(eqLoc+1);
      	RemoveWhiteSpace(outVals);
      	//cerr<<"Putting in variable "<<token<< "=" << outVals<<endl;
      	currentTree->input[token]=outVals;
      }
      

    }
    cerr<<"done Reading input"<<endl;
  }
  
  


};


#endif