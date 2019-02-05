#ifndef FORCE_FIELD_MANAGER_HPP
#define FORCE_FIELD_MANAGER_HPP

#include <iostream> 
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>

#include <proteinManager.hpp>

namespace proteinManager {
namespace ffManager{
    
    struct atomFFproperties{
	
	real charge;
	real c6;
	real c12;
	
	atomFFproperties(real charge,real c6,real c12):charge(charge),c6(c6),c12(c12){}
    };
    
    class forceFieldManager{
	
	private:
		
		//Data structure where the innformation about the force field is stored
		// ResName1 -----> Atom1_1 -----> atomFFProperties
		//  			   Atom1_2 -----> atomFFProperties
		//  			   Atom1_3 -----> atomFFProperties
		// ResName2 -----> Atom2_1 -----> atomFFProperties
		//  			   Atom2_2 -----> atomFFProperties
		//  			   Atom2_3 -----> atomFFProperties
	    std::map<std::string,std::shared_ptr<std::map<std::string,std::shared_ptr<atomFFproperties>>>> forceFieldData;
	    
	public:
		
		//This function is used for loading the force field data from file
	    void loadForceFieldData(std::string forceFieldDataFilePath){
			
			std::stringstream ss;
			
			//Check if file exists
			std::ifstream forceFieldDataFile(forceFieldDataFilePath);
			if(!forceFieldDataFile){
				ss.clear();
				ss << "File not found: " << forceFieldDataFilePath;
				throw std::runtime_error(ss.str());
			}
			
			//Processing file
			std::string line;
			
			std::string resNameBuffer;
			std::string atomNameBuffer;
			
			double chgBuffer;
			double c6Buffer;
			double c12Buffer;
			
			while(std::getline(forceFieldDataFile,line)){
				
				//Empty lines or lines which starts with # are ignored
				if(line[0]=='#' or line.empty()) continue;
				
				//Process line
				ss.str(line);
				ss >> resNameBuffer >> atomNameBuffer >> chgBuffer >> c6Buffer >> c12Buffer;
				
				//Check format
				if(ss.fail()){
					ss.clear();
					ss << "Format error in file: \"" << forceFieldDataFilePath
					<< "\" ,the following line couldn't be processed: \"" << line << "\"";
					throw std::runtime_error(ss.str());
				}
				
				//If the residue has not been added before, add it now
				if(forceFieldData.count(resNameBuffer) == 0){
				forceFieldData.insert(std::make_pair(resNameBuffer,std::make_shared<std::map<std::string,std::shared_ptr<atomFFproperties>>>()));
				forceFieldData[resNameBuffer]->insert(std::make_pair(atomNameBuffer,std::make_shared<atomFFproperties>(chgBuffer,c6Buffer,c12Buffer)));
				} else {
					//We care about not adding the same atom twice
					if(forceFieldData[resNameBuffer]->count(atomNameBuffer) > 0){
						ss.clear();
						ss << "Format error in file: \"" << forceFieldDataFilePath
						   << "\" ,the following atom parameters have been added previously: \"" << line << "\"";
						throw std::runtime_error(ss.str());
					} else {
						forceFieldData[resNameBuffer]->insert(std::make_pair(atomNameBuffer,std::make_shared<atomFFproperties>(chgBuffer,c6Buffer,c12Buffer)));
					}
				}
            }
		
		#ifdef DEBUG
		//Check
        real maxChg = std::numeric_limits<real>::lowest();
        real maxC6  = std::numeric_limits<real>::lowest();
        real maxC12 = std::numeric_limits<real>::lowest();
        
        real minChg = std::numeric_limits<real>::max();;
        real minC6  = std::numeric_limits<real>::max();;
        real minC12 = std::numeric_limits<real>::max();;
        
		for(auto& ff_entry : forceFieldData){
			std::cout << ff_entry.first << std::endl;
			for(auto& ff_a_entry: *ff_entry.second){
                
                if(ff_a_entry.second->charge > maxChg) maxChg = ff_a_entry.second->charge;
                if(ff_a_entry.second->c6 > maxC6) maxC6 = ff_a_entry.second->c6;
                if(ff_a_entry.second->c12 > maxC12)  maxC12 = ff_a_entry.second->c12;
                
                if(ff_a_entry.second->charge < minChg) minChg = ff_a_entry.second->charge;
                if(ff_a_entry.second->c6 < minC6 and ff_a_entry.second->c6 != 0) minC6 = ff_a_entry.second->c6;
                if(ff_a_entry.second->c12 < minC12 and ff_a_entry.second->c12 != 0)  minC12 = ff_a_entry.second->c12;
                
                
			    std::cout << "\t" << ff_a_entry.first 
			              << " " << ff_a_entry.second->charge 
				      << " c6 "  << ff_a_entry.second->c6 << " c6(sqrt) " << sqrt(ff_a_entry.second->c6)
				      << " c12 " << ff_a_entry.second->c12 << " c12(sqrt) " << sqrt(ff_a_entry.second->c12)
				      << std::endl;
			}
        }
        
        std::cout << std::endl;
        std::cout << "maxChg: " << maxChg << " maxC6: " << maxC6 << " maxC6 (sqrt) " << sqrt(maxC6) << " maxC12: " << maxC12 << " maxC12 (sqrt) " << sqrt(maxC12) << std::endl;
        std::cout << "minChg: " << minChg << " minC6: " << minC6 << " minC6 (sqrt) " << sqrt(minC6) << " minC12: " << minC12 << " minC12 (sqrt) " << sqrt(minC12) << std::endl;
        std::cout << std::endl;
		#endif
	    }
	    
	    //This function set the force field values (chg,c6,c12) for each atom in the structure structIn
	    void applyForceFieldData(proteinManager::STRUCTURE& structIn){
			
			std::stringstream ss;
			
			if(forceFieldData.empty()){
				ss.clear();
				ss << "No force field loaded" ;
				throw std::runtime_error(ss.str());
				
			}
			
			for(proteinManager::MODEL& md : structIn.model()) {
			for(proteinManager::CHAIN& ch : md.chain()) {
			for(proteinManager::RESIDUE& res : ch.residue()) {
				if(forceFieldData.count(res.getResName()) == 0){
					ss.clear();
					ss << "The residue " << res.getResName() << " ("<< res.getResSeq() << ") is not present in the current force field";
					throw std::runtime_error(ss.str());
				}
				
				for(proteinManager::ATOM& atm : res.atom()) {
					if((*forceFieldData[res.getResName()]).count(atm.getAtomName()) == 0){
						std::cerr << "WARNING: The parameters for the atom " << atm.getAtomName() 
								  << " in the residue " << res.getResName() << "("<< res.getResSeq() << ")" 
								  << " is not present in the current force field. I will assume they are equal to 0." << std::endl;
                                  
                        atm.setAtomCharge(0);
                        atm.setAtomC6(0);
                        atm.setAtomC12(0);
                                  
					} else {
						atm.setAtomCharge((*forceFieldData[res.getResName()])[atm.getAtomName()]->charge);
						atm.setAtomC6(sqrt((*forceFieldData[res.getResName()])[atm.getAtomName()]->c6));
						atm.setAtomC12(sqrt((*forceFieldData[res.getResName()])[atm.getAtomName()]->c12));
					}
					}
			}}}
		
	    }
    };
    
}
}

#endif
