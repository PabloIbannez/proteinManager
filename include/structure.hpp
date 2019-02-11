/*Pablo Ibáñez Freire, pablo.ibannez@uam.es*/

/*
 The structure is the highest order entity. As a special feature, it is responsible for managing the input and output of data to and from files.

 USAGE:
    
    STRUCTURE pdb; //Declaration
    
    pdb.setResSeqOverWrite(true); //Overwrite options are activated
    pdb.setSerialOverWrite(true);
    
    pdb.loadPDB("asdf.pdf"); //Data is loaded from pdb file

    pdb.setOuputFormat(proteinManager::DATA_FORMAT::PQR); //We change the output format

    std::cout << pdb << std::endl; //Write

*/

#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <vector>

#include "proteinManager.hpp"

namespace proteinManager {

    class MODEL;
    class CHAIN;
    class RESIDUE;
    class ATOM;
    
    class STRUCTURE {
    
        friend MODEL;
        friend CHAIN;
        friend RESIDUE;
        friend ATOM;
    
        friend std::ostream& operator<<(std::ostream& os, const STRUCTURE& structure);
        friend std::ostream& operator<<(std::ostream& os, const MODEL& md);
        friend std::ostream& operator<<(std::ostream& os, const CHAIN& res);
        friend std::ostream& operator<<(std::ostream& os, const RESIDUE& res);
        friend std::ostream& operator<<(std::ostream& os, const ATOM& a);
    
    private:
        
        //Default input/output format is PDB
        DATA_FORMAT inputFormat = PDB;
        DATA_FORMAT outputFormat = PDB;
        
        //Atom serial and residue sequence number can been overwritten. 
        //The values that appear in the input file for these variables can be replaced by others. 
        //They are numbered as they are read, sequentially.
        //Disabled by default.
        bool overwrite_Serial = false;
        bool overwrite_ResSeq = false;
        
        //Vector where references to models of the structure are stored.
        boost::ptr_vector<MODEL> modelVector;
    
    public:
    
        /*
        ~STRUCTURE(){
            std::cerr << "Calling structure destructor" << std::endl;
        }*/
        
        //Internal options setting
        void        setOutputFormat(DATA_FORMAT output);
        void        setSerialOverWrite(bool option);
        void        setResSeqOverWrite(bool option);
        
        //Acces to models
        MODEL& model(int modelID);
        boost::ptr_vector<MODEL>& model();
        
        bool isModel(int modelID);
        
        //Input functions
        void loadPDB(std::string inputFileName);
        void loadGRO(std::string inputFileName);
        
        //Set of functions for adding different types of lower-order entities.
        void addModel(int modelID);
    
        void addChain(int modelID,
                      std::string chainID);
    
        void addResidue(int modelID,
                        std::string chainID,
                        std::string resName,int resSeq,std::string iCode);
    
        void addAtom(int modelID,
                     std::string chainID,
                     std::string resName, int resSeq, std::string iCode,
                     int serial, std::string name);
                     
        void renumber();
    
    };
}

#endif
