#include "DataModel.h"

DataModel::DataModel(): DAQDataModelBase(), postgres_helper(this) {
	// we must do this with modern ToolDAQFramework to turn off
	// typechecking with the CStore, which is necessary for use with JsonParser
	CStore = BoostStore(false,0);
}


TTree* DataModel::GetTTree(std::string name){
  
  if(m_trees.count(name)>0) return m_trees[name];
  else return 0;
  
}


void DataModel::AddTTree(std::string name,TTree *tree){
  
  m_trees[name]=tree;

}


void DataModel::DeleteTTree(std::string name){
  
  m_trees.erase(name);
}

std::ostream& operator<<(std::ostream& os, State s){
  switch(s){
  case ReplaceWater    :os << "ReplaceWater";    break;
  case Sleep    :os << "Sleep";    break;
  case PowerUP    :os << "PowerUP";    break;
  case Dark    :os << "Dark";    break;
    //case LEDW    :os << "LEDW";    break;//
    //case LED385L    :os << "LED385L";    break;//
  case DLED260J    :os << "DLED260J";    break;
  case LED260J    :os << "LED260J";    break;
  case DLED275J    :os << "DLED275J";    break;
  case LED275J    :os << "LED275J";    break;
    //case LEDR    :os << "LEDR";    break;//
    //case LEDG    :os << "LEDG";    break;//
    //case LEDB    :os << "LEDB";    break;//
    //case Dark1    :os << "Dark1";    break;  
    //case Dark2    :os << "Dark2";    break;  
    //case Dark3    :os << "Dark3";    break;  
  case DAll    :os << "DAll";    break;
  case All    :os << "All";    break;
  case Analyse     :os << "Analyse";    break;
    
  }
  return os;
}


/*
GdTree* DataModel::GetGdTree(std::string name)
{
//	std::cout << "map size " << m_gdtrees.size() << std::endl;
//	std::cout << "GETTING from map " << name << "." << std::endl;
	if (m_gdtrees.count(name))
		return m_gdtrees[name];
	else
		return 0;
}


void DataModel::AddGdTree(std::string name, GdTree *tree)
{
	std::cout << "map size " << m_gdtrees.size() << std::endl;
	std::cout << "ADDING to map " << name << ",\t" << tree << "." << std::endl;
	m_gdtrees[name]=tree;
	std::cout << m_gdtrees[name] << std::endl;
}


void DataModel::DeleteGdTree(std::string name)
{
	std::cout << "map size " << m_gdtrees.size() << std::endl;
	std::cout << "REMOVING from map " << name << "." << std::endl;
	delete m_gdtrees[name];
	m_gdtrees.erase(name);
}

int DataModel::SizeGdTree()
{
	return m_gdtrees.size();
}
*/
