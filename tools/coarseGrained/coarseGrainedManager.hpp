#ifndef COARSE_GRAINED_MANAGER_HPP
#define COARSE_GRAINED_MANAGER_HPP

#include <iostream> 
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <proteinManager.hpp>

namespace proteinManager{
namespace coarseGrainedManager{
	
	struct aminoAcidType{
		
		std::string name;
		
		std::vector<std::vector<std::string>> beadComponents;
		
		aminoAcidType(std::string name):name(name){}
		
	};
	
	struct elementType{
		
		std::string name;
		
		double mass;
		
		elementType(std::string name):name(name){}
	};
	
	struct beadType{
		std::string name;
		
		std::vector<std::string> atomComponents;
		
		beadType(std::string name):name(name){}
	};
	
	
	
	class coarseGrainedGenerator{
		
		private:
		
		//Elements types
		
		std::map<std::string,elementType*> elementsTypes;
		
		//CG model map
		
		std::map<std::string,aminoAcidType*> aminoAcidTypes;
		std::map<std::string,beadType*> beadTypes;
		
		
		public:
		
		//This function add info about the used elements e.g. mass
		void loadElementsData(std::string elementsTypesFilePath){
		
		std::stringstream ss;
		
		//Check if file exists
		std::ifstream elementTypesFile(elementsTypesFilePath);
		if(!elementTypesFile){
			ss.clear();
			ss << "File not found: " << elementsTypesFilePath;
			throw std::runtime_error(ss.str());
		}
		
		//Processing file
		std::string line;
		
		std::string stBuffer;
		double dbBuffer;
		
		while(std::getline(elementTypesFile,line)){
			
			//Empty lines or lines which starts with # are ignored
			if(line[0]=='#' or line.empty()) continue;
			
			ss.str(line);
			ss >> stBuffer >> dbBuffer;
			
			if(ss.fail()){
			ss.clear();
			ss << "Format error in file: \"" << elementsTypesFilePath
			<< "\" ,the following line couldn't be processed: \"" << line << "\"";
			throw std::runtime_error(ss.str());
			}
			
			//Check if the element has been added before
			if(elementsTypes.count(stBuffer) == 0){
			//If not, add it
			elementsTypes.insert(std::make_pair(stBuffer,new elementType(stBuffer)));
			elementsTypes[stBuffer]->mass = dbBuffer;
			} else {
			//Else an error is thrown
			ss.clear();
			ss << "The element type \"" << stBuffer << "\" has been added previously." ;
			throw std::runtime_error(ss.str());
			}
			
			ss.clear();
		}
		
		/*
		//Check
		for(auto &e : elementsTypes){
			std::cout << e.first << " " << e.second->name << " " << e.second->mass << std::endl;
		}*/
		
		}
		
		void loadCGmodel(std::string aminoAcid_Bead_Map_FilePath, std::string bead_atom_Map_FilePath){
		
		//Auxiliar variables
		std::stringstream ss;
		std::string line;
		
		std::string stBuffer1;
		std::string stBuffer2;
		
		unsigned first;
		unsigned last;
		
		//Loading of bead-atom map
		
		std::ifstream bead_atom_Map_File(bead_atom_Map_FilePath);
		//Check if file exists
		if(!bead_atom_Map_File){
			ss.clear();
			ss << "File not found: " << bead_atom_Map_FilePath;
			throw std::runtime_error(ss.str());
		}
		
		//Processing file
		while(std::getline(bead_atom_Map_File,line)){
			
			//Empty lines or lines which starts with # are ignored
			if(line[0]=='#' or line.empty()) continue;
			
			//New bead type detected
			if(line[0]=='['){
			
			first = line.find('['); first ++;
			last =  line.find(']');
			
			stBuffer1 = line.substr (first,last-first);
			stBuffer1.erase(std::remove_if(stBuffer1.begin(), stBuffer1.end(),
				[](unsigned char const c){ return std::isspace(c);}),
				stBuffer1.end());
			
			//Check if bead type has been added previously
			if(beadTypes.count(stBuffer1) > 0){
				ss.clear();
				ss << "The bead \"" << stBuffer1 << "\" has been introduced before. File: \"" << bead_atom_Map_FilePath << "\"";
				throw std::runtime_error(ss.str());
			} else {
				//Else, add it
				beadTypes.insert(std::make_pair(stBuffer1,new beadType(stBuffer1)));
			}
			
			//Add atoms for each bead
			std::getline(bead_atom_Map_File,line);
			ss.str(line);
			while(ss >> stBuffer2){
				
				//Empty lines or lines which starts with # are ignored
				if(line[0]=='#' or line.empty()) continue;
				
				if(ss.fail() or ss.str()[0] == '['){
				ss.clear();
				ss << "Format error in file: \"" << bead_atom_Map_FilePath
				<< "\" ,the following line couldn't be processed: \"" << line << "\"";
				throw std::runtime_error(ss.str());
				}
				
				beadTypes[stBuffer1]->atomComponents.push_back(stBuffer2);
			}
			
			ss.clear();
			
			} else {
			ss.clear();
			ss << "Format error in file: \"" << bead_atom_Map_FilePath
			<< "\" ,the following line couldn't be processed: \"" << line << "\"";
			throw std::runtime_error(ss.str());
			}
		}
		
		//Check if each atom has a corresponding element
		for(auto& b: beadTypes){
			for(auto& a: b.second->atomComponents){
			if(elementsTypes.count(a.substr(0,1)) == 0){
				ss.clear();
				ss << "Error in file: \"" << bead_atom_Map_FilePath << "\"" 
				<< ". The atom \"" << a << "\" in the bead \"" << b.first << "\" can not be associated with an element";
				throw std::runtime_error(ss.str());
			}
			}
		}
		
		////////////////////////////////////////////////////////////////
		
		//Loading of aminoAcid-bead map
		std::ifstream aminoAcid_Bead_Map_File(aminoAcid_Bead_Map_FilePath);
		//Check if file exists
		if(!aminoAcid_Bead_Map_File){
			ss.clear();
			ss << "File not found: " << aminoAcid_Bead_Map_FilePath;
			throw std::runtime_error(ss.str());
		}
		
		std::streampos posBuffer = -1;
		
		//Process file
		while(std::getline(aminoAcid_Bead_Map_File,line)){
			
			//Empty lines or lines which starts with # are ignored
			if(line[0]=='#' or line.empty()) continue;
			
			//New amino acid detected
			if(line[0]=='['){
			
			first = line.find('['); first ++;
			last =  line.find(']');
			
			stBuffer1 = line.substr (first,last-first);
			stBuffer1.erase(std::remove_if(stBuffer1.begin(), stBuffer1.end(),
				[](unsigned char const c){ return std::isspace(c);}),
				stBuffer1.end());
				
			//Check if an amino acid type has been added previously
			if(aminoAcidTypes.count(stBuffer1) > 0){
				ss.clear();
				ss << "The amino acid \"" << stBuffer1 << "\" has been introduced before. File: \"" << bead_atom_Map_FilePath << "\"";
				throw std::runtime_error(ss.str());
			} else {
				//Else, add it
				aminoAcidTypes.insert(std::make_pair(stBuffer1,new aminoAcidType(stBuffer1)));
			}
			
			//Each amino acid type can have several bead configurations
			//They are processed now
			while(std::getline(aminoAcid_Bead_Map_File,line)){
				
				//Empty lines or lines which starts with # are ignored
				if(line[0]=='#' or line.empty()) continue;
				
				//New amino acid found, go back
				if(line[0]=='['){
				if(posBuffer == -1){
					ss.clear();
					ss << "Format error in file: \"" << aminoAcid_Bead_Map_FilePath
					<< "\" ,the following line couldn't be processed: \"" << line << "\"";
					throw std::runtime_error(ss.str());
				}
				aminoAcid_Bead_Map_File.seekg(posBuffer);
				posBuffer = -1;
				break;
				}
				posBuffer = aminoAcid_Bead_Map_File.tellg();
				
				//Increase the number of possible configurations for the current amino acid by one
				aminoAcidTypes[stBuffer1]->beadComponents.resize(aminoAcidTypes[stBuffer1]->beadComponents.size()+1);
				
				//Processing the new configuration
				ss.str(line);
				while(ss >> stBuffer2){
				
				if(ss.fail()){
					ss.clear();
					ss << "Format error in file: \"" << aminoAcid_Bead_Map_FilePath
					<< "\" ,the following line couldn't be processed: \"" << line << "\"";
					throw std::runtime_error(ss.str());
				}
				
				//Check if each bead has been defined before
				if(beadTypes.count(stBuffer2) == 0){
					ss.clear();
					ss << "Error in file: \"" << aminoAcid_Bead_Map_FilePath << "\"" 
					<< ". The bead \"" << stBuffer2 << "\" in the amino acid \"" << stBuffer1 << "\" is not defined";
					throw std::runtime_error(ss.str());
				} else {
					aminoAcidTypes[stBuffer1]->beadComponents.back().push_back(stBuffer2);
				}
				}
				ss.clear();
			}
			
			} else {
			ss.clear();
			ss << "Format error in file: \"" << aminoAcid_Bead_Map_FilePath
			<< "\" ,the following line couldn't be processed: \"" << line << "\"";
			throw std::runtime_error(ss.str());
			}
		}
		
		/*
		//Check
		for(auto& b:beadTypes){
			std::cout << b.first << " " << b.second->name << std::endl << std::endl;
			
			for(auto& a:b.second->atomComponents){
			std::cout << a << " ";
			}
			std::cout << std::endl << std::endl;
			
		}
		
		for(auto& am:aminoAcidTypes){
			std::cout << am.first << " " << am.second->name << std::endl << std::endl;
			
			for(auto& bc:am.second->beadComponents){
			for(auto& b:bc){
				std::cout << b << " " ;
			}
			std::cout << std::endl;
			}
			std::cout << std::endl << std::endl;
		}*/
		}
		
		
		void applyCoarseGrainedMap(proteinManager::STRUCTURE& structIn,proteinManager::STRUCTURE& structOut){
		
		std::stringstream ss;
		
		bool patternMatching;
		int  atomCount;
        int  addedAtomCount = 0;
		
		proteinManager::real3 pos;
		proteinManager::real  chg;
		proteinManager::real  totalMass;
		
		for(proteinManager::MODEL& md : structIn.model()) {
			structOut.addModel(md.getModelId());
			for(proteinManager::CHAIN& ch : md.chain()) {
			structOut.model().back().addChain(ch.getChainId());
			for(proteinManager::RESIDUE& res : ch.residue()) {
				
				structOut.model().back().chain().back().addResidue(res.getResName(),res.getResSeq(),res.getResInsCode());
				proteinManager::RESIDUE& resOut = structOut.model().back().chain().back().residue().back();
				
				if(aminoAcidTypes.count(res.getResName()) == 0){
				ss.clear();
				ss << "The amino acid \"" << res.getResName() << "\"  has not been added to the current coarse grained model" << std::endl;
				throw std::runtime_error(ss.str());
				} else {
				
				for(auto const& bdCompList : aminoAcidTypes[res.getResName()]->beadComponents){
					
					patternMatching = true;
					atomCount = 0;
					
					//Bead substitution
					//Check if the pattern matchs
					for(std::string const & bd : bdCompList){
						if(patternMatching==false) {break;}
						for(std::string const & atom : beadTypes[bd]->atomComponents){
							if(!res.isAtom(atom)){patternMatching = false; break;}
							atomCount ++;
						}
					}
					if(res.atom().size() != atomCount) {patternMatching = false;}
					
					if(patternMatching){
					
					for(std::string const & bd : bdCompList){
						
						//TODO, more generic
						
                        ////////////////////////////////////////////////
                        
						pos = {0,0,0};
						chg = 0;
						totalMass = 0;
						
						for(std::string const & atom : beadTypes[bd]->atomComponents){
                            
							chg += res.atom(atom).getAtomCharge();
							totalMass += elementsTypes[atom.substr(0,1)]->mass;
						}
						
						pos = res.atom("CA").getAtomCoord();
                        
                        ////////////////////////////////////////////////
                        
						resOut.addAtom(addedAtomCount,bd);
						addedAtomCount++;
                        
						resOut.atom(bd).setAtomCoord(pos);
						resOut.atom(bd).setAtomCharge(chg);
                        
                        //Common properties
                        resOut.atom(bd).setAtomAltLoc(" ");
                        resOut.atom(bd).setAtomOccupancy(1);
                        resOut.atom(bd).setAtomTempFactor(0);
                        resOut.atom(bd).setAtomElement("");
						
					}
					
					break;
					}
				}
				
				if(patternMatching==false){
					ss.clear();
					ss << "All the substitutions are failed for the amino acid \"" << res.getResName() 
					<< "(" << res.getResSeq() << ")" << "\"" << std::endl;
					ss << res << std::endl;
					throw std::runtime_error(ss.str());
				}
				
				}
				
			}
			}
		}
		}
		
	};

}
}

#endif
