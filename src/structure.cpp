/*Pablo Ibáñez Freire, pablo.ibannez@uam.es*/
#include "structure.hpp"

namespace proteinManager {

    void STRUCTURE::setOutputFormat(DATA_FORMAT output) {
        outputFormat = output;
    }
    
    void STRUCTURE::setSerialOverWrite(bool option) {
        overwrite_Serial=option;
    }
    
    void STRUCTURE::setResSeqOverWrite(bool option) {
        overwrite_ResSeq=option;
    }
    
    MODEL& STRUCTURE::model(int modelID) {
        for(MODEL& md : modelVector) {
            if(md.getModelId() == modelID) {
                return md;
            }
        }
    
        std::stringstream ss;
        ss << "No model added with the modelID \"" << modelID << "\"";
        throw std::runtime_error(ss.str());
    }
    
    boost::ptr_vector<MODEL>& STRUCTURE::model() {
        return modelVector;
    }
    
    boost::ptr_vector<CHAIN>& STRUCTURE::chain(){

        chainVector.clear();

        for(MODEL& mdl : modelVector){
            chainVector.insert(chainVector.end(),mdl.chain().begin(),mdl.chain().end());
        }

        return chainVector;
    }
    
    boost::ptr_vector<RESIDUE>& STRUCTURE::residue(){

        residueVector.clear();

        for(MODEL& mdl : modelVector){
            residueVector.insert(residueVector.end(),mdl.residue().begin(),mdl.residue().end());
        }

        return residueVector;
    }
    
    boost::ptr_vector<ATOM>& STRUCTURE::atom(){

        atomVector.clear();

        for(MODEL& mdl : modelVector){
            atomVector.insert(atomVector.end(),mdl.atom().begin(),mdl.atom().end());
        }

        return atomVector;
    }
    
    bool STRUCTURE::isModel(int modelID){
        
        for(MODEL& mdl : modelVector){
            if(mdl.getModelId() == modelID){
                return true;
            }
        }
        
        return false;
    }
    
    void STRUCTURE::loadPDB(std::string inputFileName) {
    
        std::stringstream ss;
    
        std::ifstream inputFile = std::ifstream(inputFileName);
        if (inputFile.fail()) {
            ss.clear();
            ss << "Error loading file \"" << inputFileName << "\".";
            throw std::runtime_error(ss.str());
        } else {
            std::cerr << "File \"" << inputFileName << "\" opened." << std::endl;
        }
    
        std::string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);
    
        if (extension == "pdb") {
            inputFormat = PDB;
            outputFormat = PDB;
        } else if (extension == "pdbq") {
            inputFormat = PDBQ;
            outputFormat = PDBQ;
        } else if (extension == "pqr") {
            inputFormat = PQR;
            outputFormat = PQR;
        } else if (extension == "pdrs") {
            inputFormat = PDRS;
            outputFormat = PDRS;
        } else {
            ss.clear();
            ss << "ERROR. Incorrect input file format.";
            throw std::runtime_error(ss.str());
        }
    
        bool readingInitiated = false;
        bool modelEnded = false;
    
        std::string line;
        int modelCount;
    
        int serial;
        std::string name;
        std::string altLoc;
        std::string resName;
        std::string chainID;
        int resSeq;
        std::string iCode;
        real3 coord;
        real occupancy;
        real tempFactor;
        std::string element;
        real charge;
        real radius;
        real SASA;
    
        int serialCounter;
    
        int resSeqCounter;
        int prevRes;
    
        while(std::getline(inputFile,line)) {
    
            if(line.compare(0,5,"MODEL") == 0) {
    
                if(readingInitiated == true && modelEnded == false) {
                    ss.clear();
                    ss << "ERROR (STRUCTURE). A new model starts but \"ENDMDL\" label has not been read. Error line:" << std::endl;
                    ss << line ;
                    throw std::runtime_error(ss.str());
                }
    
                if(readingInitiated == false) {
                    readingInitiated = true;
                    modelEnded = false;
                    modelCount = 0;
                    this->addModel(modelCount);
                } else {
                    readingInitiated = true;
                    modelEnded = false;
                    modelCount ++;
                    this->addModel(modelCount);
    
                }
    
                serialCounter = 1;
                resSeqCounter = 0;
                prevRes = -1;
    
            }
    
            if(line.compare(0,6,"ENDMDL") == 0) {
                modelEnded = true;
                std::cerr << "Loading of model \"" << modelVector.size()
                        << "\" complete" << std::endl;
            }
    
            if(line.compare(0,3,"TER") == 0) {
                resSeqCounter = 0;
                prevRes = -1;
                serialCounter ++ ;
            }
    
            if(line.compare(0,4,"ATOM") == 0) {
    
                if(readingInitiated == false) {
                    readingInitiated = true;
                    modelCount =   0;
                    serialCounter = 1;
                    resSeqCounter = 0;
                    prevRes = -1;
                    this->addModel(modelCount);
                }
    
                if(modelEnded == true) {
                    ss.clear();
                    ss << "ERROR (STRUCTURE). A model has ended but no \"MODEL\" label has been read. Error line:";
                    ss << line << std::endl;
                    throw std::runtime_error(ss.str());
                }
    
                switch(inputFormat) {
                    case PDB:
        
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
        
                            chainID    =                    line.substr(21,1);
                            resSeq     =          std::stoi(line.substr(22,4));
                            iCode      =                    line.substr(26,1);
                            coord      =    {real(std::stod(line.substr(30,8))),
                                            real(std::stod(line.substr(38,8))),
                                            real(std::stod(line.substr(46,8)))
                                            };
                            occupancy  =     real(std::stod(line.substr(54,6)));
                            tempFactor =     real(std::stod(line.substr(60,6)));
                            element    =                    line.substr(76,2);
        
                            try {
                                charge =                   (line.substr(78,2).empty())?0:real(std::stod(line.substr(78,2)));
                            } catch (std::invalid_argument) {
                                charge =   0;
                            }
        
                            radius     =   0;
                            SASA       =  -1;
                        } catch (const std::exception& e) {
                            ss.clear();
                            ss << "ERROR ( " << e.what() << " ). Found while processing line: " << line ;
                            throw std::runtime_error(ss.str());
                        }
        
                        break;
        
                    case PDBQ:
        
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
        
                            chainID    =                    line.substr(21,1);
                            resSeq     =          std::stoi(line.substr(22,4));
                            iCode      =                    line.substr(26,1);
                            coord      =    {real(std::stod(line.substr(30,8))),
                                            real(std::stod(line.substr(38,8))),
                                            real(std::stod(line.substr(46,8)))
                                            };
                            occupancy  =     real(std::stod(line.substr(54,6)));
                            tempFactor =     real(std::stod(line.substr(60,6)));
                            element    =                    "";
                            charge     =     real(std::stod(line.substr(70,6)));
                            radius     =   0;
                            SASA       =  -1;
                        } catch (const std::exception& e) {
                            ss.clear();
                            ss << "ERROR ( " << e.what() << " ). Found while processing line: " << line ;
                            throw std::runtime_error(ss.str());
                        }
        
                        break;
        
                    case PQR:
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
        
                            chainID    =                    line.substr(21,1);
                            resSeq     =          std::stoi(line.substr(22,4));
                            iCode      =                    line.substr(26,1);
                            coord      =    {real(std::stod(line.substr(30,8))),
                                            real(std::stod(line.substr(38,8))),
                                            real(std::stod(line.substr(46,8)))
                                            };
                            occupancy  =   1;
                            tempFactor =   0;
                            element    =  "";
                            charge     =     real(std::stod(line.substr(54,8)));
                            radius     =     real(std::stod(line.substr(62,8)));
                            SASA       =  -1;
                        } catch (const std::exception& e) {
                            ss.clear();
                            ss << "ERROR ( " << e.what() << " ). Found while processing line: " << line ;
                            throw std::runtime_error(ss.str());
                        }
        
                        break;
        
                    case PDRS:
        
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
        
                            chainID    =                    line.substr(21,1);
                            resSeq     =          std::stoi(line.substr(22,4));
                            iCode      =                    line.substr(26,1);
                            coord      =    {real(std::stod(line.substr(30,8))),
                                            real(std::stod(line.substr(38,8))),
                                            real(std::stod(line.substr(46,8)))
                                            };
                            occupancy  =   1;
                            tempFactor =   0;
                            radius     =     real(std::stod(line.substr(54,6)));
                            SASA       =     real(std::stod(line.substr(60,6)));
                        } catch (const std::exception& e) {
                            ss.clear();
                            ss << "ERROR ( " << e.what() << " ). Found while processing line: " << line ;
                            throw std::runtime_error(ss.str());
                        }
                        break;
        
                }
    
                if(overwrite_ResSeq ) {
                    if(prevRes != resSeq) {
                        resSeqCounter ++;
                        prevRes = resSeq;
                    }
    
                    resSeq = resSeqCounter;
                }
    
                if(overwrite_Serial) {
                    serial = serialCounter;
                    serialCounter ++;
                }
    
                modelVector[modelCount].addAtom(chainID,
                                                resName,resSeq,iCode,
                                                serial,name);
                
                //Common for all formats
                modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomAltLoc(altLoc);
                modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomCoord(coord);
                modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomOccupancy(occupancy);
                modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomTempFactor(tempFactor);
                
                switch(inputFormat) {
                    case PDB:
                        
                        modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomCharge(charge);
                        modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomElement(element);
        
                        break;
        
                    case PDBQ:
                        
                        modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomCharge(charge);
        
                        break;
        
                    case PQR:
                        
                        modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomCharge(charge);
                        modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomRadius(radius);
        
                        break;
        
                    case PDRS:
                        
                        modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomRadius(radius);
                        modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomSASA(SASA);
        
                       break;
        
                    }
            }
        }
    
        std::cerr << "File \"" << inputFileName << "\"" << " loaded." << std::endl;
    
    }
    
    void STRUCTURE::loadGRO(std::string inputFileName) {
    
        std::stringstream ss;
    
        std::ifstream inputFile = std::ifstream(inputFileName);
        if (inputFile.fail()) {
            ss.clear();
            ss << "Error loading file \"" << inputFileName << "\".";
            throw std::runtime_error(ss.str());
        } else {
            std::cerr << "File \"" << inputFileName << "\" opened." << std::endl;
        }
    
        std::string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);
    
        if (extension == "gro") {
            inputFormat = GRO;
            outputFormat = GRO;
        } else {
            ss.clear();
            ss << "ERROR. Incorrect input file format." ;
            throw std::runtime_error(ss.str());
        }
    
        std::regex rgBoxPattern = std::regex(R"(^\s*([-+]?[0-9]*\.?[0-9]*)\s+([-+]?[0-9]*\.?[0-9]*)\s+([-+]?[0-9]*\.?[0-9]*)\s*$)");
        std::regex rgAtomPattern = std::regex(R"(^\s+\d+[a-zA-Z]+)");
        std::smatch sm;
    
        bool readingInitiated = false;
        bool modelEnded = false;
    
        std::string line;
        int modelCount;
    
        int atomsCurrentModel;
        int atomsCount;
    
        int chainCount;
    
        int resSeqCount;
    
        int serial;
        std::string name;
        std::string altLoc;
        std::string resName;
        std::string chainID;
        int resSeq;
        std::string iCode;
        real3 coord;
        real occupancy;
        real tempFactor;
        std::string element;
        real charge;
        real radius;
    
        while(std::getline(inputFile,line)) {
    
            if((line.find("t=")  != std::string::npos) or
                    (line.find("t =") != std::string::npos) or
                    readingInitiated == false) {
    
                if(readingInitiated == true && modelEnded == false) {
                    ss.clear();
                    ss << "ERROR (STRUCTURE). A new model starts but previous one has not finished. Error line:";
                    ss << line << std::endl;
                    throw std::runtime_error(ss.str());
                }
    
                if(readingInitiated == false) {
                    readingInitiated = true;
                    modelEnded = false;
                    modelCount  = 0;
                    atomsCount  = 0;
                    chainCount  = 0;
                    resSeqCount = 0;
                    this->addModel(modelCount);
                } else {
                    readingInitiated = true;
                    modelEnded = false;
                    modelCount ++;
                    atomsCount  = 0;
                    chainCount  = 0;
                    resSeqCount = 0;
                    this->addModel(modelCount);
    
                }
    
                std::getline(inputFile,line);
                try {
                    atomsCurrentModel = std::stoi(line);
                } catch (...) {
                    ss.clear();
                    ss << "The number of atoms in the current frame/model has not been indicated correctly. Error line:" << std::endl;
                    ss << line ;
                    throw std::runtime_error(ss.str());
                }
    
                continue;
            }
    
            if(std::regex_search(line, sm, rgBoxPattern)) {
                modelEnded = true;
    
                if(atomsCurrentModel != atomsCount) {
                    ss.clear();
                    ss << "The number of atoms indecated and the number of atoms loaded in the current model don't match" ;
                    throw std::runtime_error(ss.str());
                }
    
                std::cerr << "Loading of model \"" << modelVector.size()
                        << "\" complete" << std::endl;
                continue;
            }
    
            if(std::regex_search(line, sm, rgAtomPattern)) {
    
                if(modelEnded == true) {
                    ss.clear();
                    ss << "ERROR (STRUCTURE). A model has ended but no title string has been read. Error line:" << std::endl;
                    ss << line;
                    throw std::runtime_error(ss.str());
                }
    
                serial     =          std::stoi(line.substr(15,5));
                name       =                    line.substr(10,5);
                name.erase(std::remove_if(name.begin(), name.end(), [](unsigned char const c) {
                    return std::isspace(c);
                }), name.end());
    
                altLoc     =                    "";
                resName    =                    line.substr(5,5);
                resName.erase(std::remove_if(resName.begin(), resName.end(), [](unsigned char const c) {
                    return std::isspace(c);
                }), resName.end());
    
                resSeq     =          std::stoi(line.substr(0,5));
    
                if (resSeq == 1 and resSeq != resSeqCount) {
                    resSeqCount = resSeq;
                    chainCount ++;
                } else if(resSeq == resSeqCount or resSeq == resSeqCount + 1) {
                    resSeqCount = resSeq;
                } else {
                    ss.clear();
                    ss << "The current residue number is not expected. Error line:";
                    ss << line << std::endl;
                    throw std::runtime_error(ss.str());
                }
    
                chainID    =          std::to_string(chainCount);
    
                iCode      =                    "";
                coord      =    {real(std::stod(line.substr(20,8))),
                                real(std::stod(line.substr(28,8))),
                                real(std::stod(line.substr(36,8)))
                                };
                occupancy  =     1;
                tempFactor =     0;
                element    =                    "";
                charge     =     0;
                radius     =     0;
    
                modelVector[modelCount].addAtom(chainID,
                                                resName,resSeq,iCode,
                                                serial,name);
    
                modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomAltLoc(altLoc);
                modelVector[modelCount].chain(chainID).residue(resSeq).atom().back().setAtomCoord(coord);
    
                atomsCount ++;
            }
        }
    
        std::cerr << "File \"" << inputFileName << "\"" << " loaded." << std::endl;
    
    }
    
    void STRUCTURE::addModel(int modelID) {
    
        for(MODEL& md : modelVector) {
            if(modelID == md.getModelId()) {
                std::stringstream ss;
                ss << "Error (MODEL). A model with the model identifier \"" << modelID << "\" has been added before.";
                throw std::runtime_error(ss.str());
            }
        }
    
        modelVector.push_back(new MODEL(modelID,this));
    }
    
    void STRUCTURE::addChain(int modelID,
                            std::string chainID) {
    
        for(MODEL& md : modelVector) {
            if(modelID == md.getModelId()) {
                md.addChain(chainID);
            }
        }
    
        modelVector.push_back(new MODEL(modelID,this));
        modelVector.back().addChain(chainID);
    
    }
    
    void STRUCTURE::addResidue(int modelID,
                            std::string chainID,
                            std::string resName,int resSeq,std::string iCode) {
    
        for(MODEL& md : modelVector) {
            if(modelID == md.getModelId()) {
                md.addResidue(chainID,
                            resName,resSeq,iCode);
            }
        }
    
        modelVector.push_back(new MODEL(modelID,this));
        modelVector.back().addResidue(chainID,
                                    resName,resSeq,iCode);
    }
    
    void STRUCTURE::addAtom(int modelID,
                            std::string chainID,
                            std::string resName, int resSeq, std::string iCode,
                            int serial, std::string name) {
    
        for(MODEL& md : modelVector) {
            if(modelID == md.getModelId()) {
                md.addAtom(chainID,
                        resName,resSeq,iCode,
                        serial,name);
            }
        }
    
        modelVector.push_back(new MODEL(modelID,this));
        modelVector.back().addAtom(chainID,
                                resName,resSeq,iCode,
                                serial,name);
    }
    
    void STRUCTURE::renumber(){
        
        std::vector<std::string> chainNames =  {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z",
                                                "a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z",
                                                "0","1","2","3","4","5","6","7","8","9"};
        
        int modelCount;
        int chainCount;
        int resCount  ;
        int atomCount ;
        
        modelCount = 1;
        atomCount  = 0;
        for(MODEL&   mdl : this->model()){
            mdl.setModelId(modelCount);
            modelCount ++ ;
            chainCount = 0;
        for(CHAIN&   ch  : mdl.chain()  ){
            ch.setChainId(chainNames[chainCount]);
            chainCount ++;
            if(chainCount == chainNames.size()){
                std::stringstream ss;
                ss << "ERROR (" << __FUNCTION__ << "). The number of chains is greater than the number of possible characters to name them.";
                throw std::runtime_error(ss.str());
            }
            resCount   = 1;
        for(RESIDUE& res : ch.residue() ){
            res.setResSeq(resCount);
            resCount ++;
        for(ATOM&    atm : res.atom()   ){
            atm.setAtomSerial(atomCount);
            atomCount ++;
        }}}}
    }
    
    std::ostream& operator<<(std::ostream& os, const STRUCTURE& structure) {
    
        DATA_FORMAT outputFormat = structure.outputFormat;
    
        switch(outputFormat) {
    
        case PDB:
        case PDBQ:
        case PQR:
        case PDRS:
            for(size_t i = 0; i < structure.modelVector.size(); i++) {
    
                if(structure.modelVector.size() == 1) {
                    os << structure.modelVector[i];
                } else {
                    os << std::left    << std::fixed <<
                    std::setw(6) << "MODEL"    <<
                    "    "                     <<
                    std::setw(4) << structure.modelVector[i].getModelId()  <<
                    std::endl;
                    os << structure.modelVector[i] << std::endl;
    
                    if( i < structure.modelVector.size()-1) {
                        os << std::left    << std::fixed <<
                        std::setw(6) << "ENDMDL"   <<
                        std::endl;
                    } else {
                        os << std::left    << std::fixed <<
                        std::setw(6) << "ENDMDL"   ;
                    }
                }
            }
    
            break;
    
        case SP:
        case SPQ:
        case XYZ:
            for(size_t i = 0; i < structure.modelVector.size()-1; i++) {
                os << "#" << std::endl;
                os << structure.modelVector[i] << std::endl;
            }
    
            os << "#" << std::endl;
            os << structure.modelVector[structure.modelVector.size()-1];
    
            break;
        }
    
        return os;
    }

}

