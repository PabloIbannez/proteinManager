#include <proteinManager.hpp>

namespace proteinManager{
    
    class massesManager{
        
        private:
            
            std::map<std::string,real> massesData;
            
        public:
            
            void loadMassesData(std::string massesDataFilePath){
			
                std::stringstream ss;
                
                //Check if file exists
                std::ifstream massesDataFile(massesDataFilePath);
                if(!massesDataFile){
                    ss.clear();
                    ss << "File not found: " << massesDataFilePath;
                    throw std::runtime_error(ss.str());
                }
                
                //Processing file
                std::string line;
                
                std::string atomNameBuffer;
                
                double massBuffer;
                
                while(std::getline(massesDataFile,line)){
                    
                    //Empty lines or lines which starts with # are ignored
                    if(line[0]=='#' or line.empty()) continue;
                    
                    //Process line
                    ss.str(line);
                    ss >> atomNameBuffer >> massBuffer;
                    
                    //Check format
                    if(ss.fail()){
                        ss.clear();
                        ss << "Format error in file: \"" << massesDataFilePath
                        << "\" ,the following line couldn't be processed: \"" << line << "\"";
                        throw std::runtime_error(ss.str());
                    }
                    
                    //We care about not adding the same atom twice
                    if(massesData.count(atomNameBuffer) > 0){
                        ss.clear();
                        ss << "Format error in file: \"" << massesDataFilePath
                        << "\" ,the following atom parameters have been added previously: \"" << line << "\"";
                        throw std::runtime_error(ss.str());
                    } else {
                        massesData.insert(std::make_pair(atomNameBuffer,massBuffer));
                    }
                    
                    ss.clear();
                    
                }
                
                #ifdef DEBUG
                    
                    for(auto& mD_entry : massesData){
                        std::cout << mD_entry.first << " " << mD_entry.second << std::endl;
                    }
                    
                #endif
            }
            
            void applyMassesData(proteinManager::STRUCTURE& structIn){
			
                std::stringstream ss;
                
                if(massesData.empty()){
                    ss.clear();
                    ss << "No masses loaded" ;
                    throw std::runtime_error(ss.str());
                    
                }
                
                for(proteinManager::MODEL& md : structIn.model()) {
                for(proteinManager::CHAIN& ch : md.chain()) {
                for(proteinManager::RESIDUE& res : ch.residue()) {
                for(proteinManager::ATOM& atm : res.atom()) {
                    if(massesData.count(atm.getAtomName()) == 0){
                        if(massesData.count(atm.getAtomName().substr(0,1)) == 0){
                            std::cerr << "WARNING: The parameters for the atom " << atm.getAtomName() 
                                    << " in the residue " << res.getResName() << "("<< res.getResSeq() << ")" 
                                    << " is not present in the current masses data file" << std::endl;
                        } else {
                            atm.setAtomMass(massesData[atm.getAtomName().substr(0,1)]);
                        }
                    } else {
                        atm.setAtomMass(massesData[atm.getAtomName()]);
                    }
                }}}}
                
            }
        
    };
    
    
}
