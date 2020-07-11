#include <proteinManager/proteinManager.hpp>

namespace proteinManager{
    
    class radiusManager{
        
        private:
            
            std::map<std::string,real> radiusData;
            
        public:
            
            void loadRadiusData(std::string radiusDataFilePath){
			
                std::stringstream ss;
                
                //Check if file exists
                std::ifstream radiusDataFile(radiusDataFilePath);
                if(!radiusDataFile){
                    ss.clear();
                    ss << "File not found: " << radiusDataFilePath;
                    throw std::runtime_error(ss.str());
                }
                
                //Processing file
                std::string line;
                
                std::string atomNameBuffer;
                
                double radiusBuffer;
                
                while(std::getline(radiusDataFile,line)){
                    
                    //Empty lines or lines which starts with # are ignored
                    if(line[0]=='#' or line.empty()) continue;
                    
                    //Process line
                    ss.str(line);
                    ss >> atomNameBuffer >> radiusBuffer;
                    
                    //Check format
                    if(ss.fail()){
                        ss.clear();
                        ss << "Format error in file: \"" << radiusDataFilePath
                        << "\" ,the following line couldn't be processed: \"" << line << "\"";
                        throw std::runtime_error(ss.str());
                    }
                    
                    //We care about not adding the same atom twice
                    if(radiusData.count(atomNameBuffer) > 0){
                        ss.clear();
                        ss << "Format error in file: \"" << radiusDataFilePath
                        << "\" ,the following atom parameters have been added previously: \"" << line << "\"";
                        throw std::runtime_error(ss.str());
                    } else {
                        radiusData.insert(std::make_pair(atomNameBuffer,radiusBuffer));
                    }
                    
                    ss.clear();
                    
                }
                
                #ifdef DEBUG
                    
                    for(auto& rD_entry : radiusData){
                        std::cout << rD_entry.first << " " << rD_entry.second << std::endl;
                    }
                    
                #endif
            }
            
            void applyRadiusData(proteinManager::MODEL& mdlIn){
			
                std::stringstream ss;
                
                if(radiusData.empty()){
                    ss.clear();
                    ss << "No radii loaded" ;
                    throw std::runtime_error(ss.str());
                    
                }
                
                for(proteinManager::CHAIN& ch : mdlIn.chain()) {
                for(proteinManager::RESIDUE& res : ch.residue()) {
                for(proteinManager::ATOM& atm : res.atom()) {
                    
                    if(radiusData.count(atm.getAtomName()) == 0){
                        if(radiusData.count(atm.getAtomName().substr(0,1)) == 0){
                            std::cerr << "WARNING: The parameters for the atom " << atm.getAtomName() 
                                    << " in the residue " << res.getResName() << "("<< res.getResSeq() << ")" 
                                    << " is not present in the current radii data file" << std::endl;
                        } else {
                            atm.setAtomRadius(radiusData[atm.getAtomName().substr(0,1)]);
                        }
                    } else {
                        atm.setAtomRadius(radiusData[atm.getAtomName()]);
                    }
                }}}
                
            }
            
            void applyRadiusData(proteinManager::STRUCTURE& structIn){
                
                for(proteinManager::MODEL& md : structIn.model()) {
                    applyRadiusData(md);
                }
                
            }
        
    };
    
    
}
