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

#include <proteinManager/proteinManager.hpp>

#include "coarseGrainedMappingSchemes.hpp"

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
		
		//CG model map
		
		std::map<std::string,aminoAcidType*> aminoAcidTypes;
		std::map<std::string,beadType*> beadTypes;
		
		public:
		
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
		
		template<class mappingSchemeType>
		void applyCoarseGrainedMap(proteinManager::STRUCTURE& structIn,proteinManager::STRUCTURE& structOut){
		
		std::stringstream ss;
		
		bool patternMatching;
		int  atomCount;
        int  addedAtomCount = 0;
        
        mappingSchemeType mSch;
		
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
                            
                            resOut.addAtom(addedAtomCount,bd);
                            addedAtomCount++;
                            
                            mSch.mappingScheme(res,resOut,bd,beadTypes[bd]->atomComponents);
                            
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
        
        template<class mappingSchemeType>
		void applyCoarseGrainedMap_IgnoreList(proteinManager::STRUCTURE& structIn,proteinManager::STRUCTURE& structOut,proteinManager::STRUCTURE& structIgnored){
		
		std::stringstream ss;
        
		bool patternMatching;
		int  atomCount;
        int  addedAtomCount = 0;
        
        bool ignoredResidue = false;
        
        mappingSchemeType mSch;
		
		for(proteinManager::MODEL& md : structIn.model()) {
            std::cout << md.getModelId() << std::endl;
			structOut.addModel(md.getModelId());
			for(proteinManager::CHAIN& ch : md.chain()) {
			structOut.model().back().addChain(ch.getChainId());
			for(proteinManager::RESIDUE& res : ch.residue()) {
                
                //std::cout << "current " << res.getModelId() << " " << res.getChainId() << " " << res.getResName() << " " << res.getResSeq() << std::endl;
                
                for(RESIDUE& resIgnored : structIgnored.model(md.getModelId()).chain(ch.getChainId()).residue() ){
                    if(ignoredResidue){ break;}
                    
                    //std::cout << std::endl << std::endl << std::endl;
                    //
                    //std::cout << res.getModelId() << " " << resIgnored.getModelId() << std::endl <<
                    //             res.getChainId() << " " << resIgnored.getChainId() << std::endl <<
                    //             res.getResName() << " " << resIgnored.getResName() << std::endl <<
                    //             res.getResSeq()  << " " << resIgnored.getResSeq()  << std::endl;
                                 
                    //std::cin.get();
                    
                    if(res.getModelId() == resIgnored.getModelId() and
                       res.getChainId() == resIgnored.getChainId() and
                       res.getResName() == resIgnored.getResName() and
                       res.getResSeq()  == resIgnored.getResSeq()){
                           //std::cout << "added " << res.getModelId() << " " << res.getChainId() << " " << res.getResName() << " " << res.getResSeq() << std::endl;
                           structOut.model(res.getModelId()).chain(res.getChainId()).addResidue(res.getResName(),res.getResSeq(),res.getResInsCode());
                           proteinManager::RESIDUE& resOut = structOut.model(res.getModelId()).chain(res.getChainId()).residue(res.getResSeq());
                           for(proteinManager::ATOM& atm : res.atom()) {
                               resOut.addAtom(atm.getAtomSerial(),atm.getAtomName());
                               
                               resOut.atom(atm.getAtomName()).setAtomCoord(atm.getAtomCoord());
                               resOut.atom(atm.getAtomName()).setAtomCharge(0);
                               resOut.atom(atm.getAtomName()).setAtomAltLoc(" ");
                               resOut.atom(atm.getAtomName()).setAtomOccupancy(1);
                               resOut.atom(atm.getAtomName()).setAtomTempFactor(0);
                               resOut.atom(atm.getAtomName()).setAtomElement("");
                               
                               
                                //std::cout << "ignored" << std::endl;
                                //std::cin.get();
                            }
                            
                            ignoredResidue = true;
                       }
                }
				
                if(!ignoredResidue){
				    structOut.model(res.getModelId()).chain(res.getChainId()).addResidue(res.getResName(),res.getResSeq(),res.getResInsCode());
				    proteinManager::RESIDUE& resOut = structOut.model(res.getModelId()).chain(res.getChainId()).residue(res.getResSeq());
				    
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
                                    
                                    resOut.addAtom(addedAtomCount,bd);
                                    addedAtomCount++;
                                    
                                    mSch.mappingScheme(res,resOut,bd,beadTypes[bd]->atomComponents);
                                    
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
                
                ignoredResidue = false;
				
			}
			}
		}
		}
		
	};

}
}

#endif
