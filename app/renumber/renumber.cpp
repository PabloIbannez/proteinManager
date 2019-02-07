#include <proteinManager.hpp>

int main(int argc, char *argv[]){
    
    std::stringstream ss;
    
    std::string inputFileName = argv[1];
    std::ifstream inputFile = std::ifstream(inputFileName);
    
    if (inputFile.fail()) {
	ss.clear();
        ss << "Error loading file \"" << inputFileName << "\".";
        throw std::runtime_error(ss.str());
    } else {
        std::cerr << "File \"" << inputFileName << "\" opened." << std::endl;
    }
    
    std::string line;
    
    int serial;
    std::string name;
    std::string altLoc;
    std::string resName;
    std::string chainID;
    int resSeq;
    std::string iCode;
    proteinManager::real3 coord;
    proteinManager::real occupancy;
    proteinManager::real tempFactor;
    std::string element;
    proteinManager::real charge;
    
    std::string resNameOld = "";
    std::vector<std::string> resComp;
    
    std::string chainID_old = "";
    
    int modelCount = 0;
    int resCount = 0;
    int atmCount = 0;
    
    while(std::getline(inputFile,line)) {
        
        if(line.compare(0,5,"MODEL") == 0)  {
            modelCount ++;
            std::cout <<  std::setw(6) << "MODEL"    <<
                          "    "                     <<
                        std::setw(4) << modelCount <<
                        std::endl;
        }
        if(line.compare(0,6,"ENDMDL") == 0) {std::cout << line << std::endl;}
        
        if(line.compare(0,3,"TER") == 0)    {
            atmCount ++;
            if(atmCount > 99999 ) { atmCount = 99999;}
            std::cout << std::left    << std::fixed <<
                         std::setw(6) << "TER"      <<
                         std::right                 <<
                         std::setw(5) << atmCount   <<
                         "      "                   <<
                         std::setw(3) << resName    <<
                         " "                        <<
                         std::setw(1) << chainID    <<
                         std::setw(4) << resCount   <<
                         std::setw(1) << iCode      <<
                         "                                                     ";
        }
        
        if(line.compare(0,3,"END") == 0)    {std::cout << line; }
        
        if(line.compare(0,4,"ATOM") == 0) {
            
            atmCount ++;
            if(atmCount > 99999 ) { atmCount = 99999;}
            
            try {
                serial     =          std::stoi(line.substr(6,5));
                name       =                    line.substr(12,4);
                name.erase(std::remove_if(name.begin(), name.end(), [](unsigned char const c) {
                    return std::isspace(c);
                }), name.end());
    
                altLoc     =                    line.substr(16,1);
                resName    =                    line.substr(17,3);
                resName.erase(std::remove_if(resName.begin(), resName.end(), [](unsigned char const c) {
                    return std::isspace(c);
                }), resName.end());
                
                if(resName != resNameOld or 
                    std::find(resComp.begin(), resComp.end(), name) != resComp.end()){ //New residue detection
                    resNameOld = resName;
                    resComp.clear();
                    resComp.push_back(name);
                    resCount++;
                } else {
                    resComp.push_back(name);
                }
    
                chainID    =                    line.substr(21,1);
                
                if(chainID != chainID_old){ //New chain detection
                    chainID_old = chainID;
                    resCount = 1;
                }
                resSeq     =          std::stoi(line.substr(22,4));
                iCode      =                    line.substr(26,1);
                coord      =    {proteinManager::real(std::stod(line.substr(30,8))),
                                 proteinManager::real(std::stod(line.substr(38,8))),
                                 proteinManager::real(std::stod(line.substr(46,8)))
                                };
                occupancy  =     proteinManager::real(std::stod(line.substr(54,6)));
                tempFactor =     proteinManager::real(std::stod(line.substr(60,6)));
                element    =                    line.substr(76,2);
    
                try {
                    charge =                   (line.substr(78,2).empty())?0:proteinManager::real(std::stod(line.substr(78,2)));
                } catch (std::invalid_argument) {
                    charge =   0;
                }
                
            } catch (const std::exception& e) { 
                ss.clear();
                ss << "ERROR ( " << e.what() << " ). Found while processing line: " << line ;
                throw std::runtime_error(ss.str());
            }
            
            ////////////////////////////////////////////////////////////
            
            std::cout  << std::left << std::fixed <<
            std::setw(6) << "ATOM"                <<
            std::right                            <<
            std::setw(5) << atmCount              <<
            " ";

            if(name.size() < 4) {
                std::cout  << std::left << std::fixed << " " <<
                std::setw(3) << name;
            } else {
                std::cout  << std::left << std::fixed <<
                std::setw(4) << name;
            }
    
            std::cout    << std::left  << std::fixed <<
            std::setw(1) << altLoc     <<
            std::setw(3) << resName    <<
            " "                        <<
            std::setw(1) << chainID    <<
            std::right                 <<
            std::setw(4) << resCount   <<
            std::setw(1) << iCode      <<
            "   "                      <<
            std::setprecision(3)       <<
            std::setw(8) << coord.x    <<
            std::setw(8) << coord.y    <<
            std::setw(8) << coord.z    <<
            std::setprecision(2)       <<
            std::setw(6) << occupancy  <<
            std::setw(6) << tempFactor <<
            "          "               <<
            std::setw(2) << element;
    
            if(int(charge)==0) {
                std::cout    << std::left << std::fixed <<
                std::setw(2) << "";
            } else {
                std::cout    << std::left << std::fixed <<
                std::setw(2) << int(charge);
            }
            }
            
            std::cout << std::endl;
    }
    
}
