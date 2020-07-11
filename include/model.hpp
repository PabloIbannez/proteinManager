/*Pablo Ibáñez Freire, pablo.ibannez@uam.es*/

#ifndef MODEL_HPP
#define MODEL_HPP

#include <vector>

#include "proteinManager.hpp"


namespace proteinManager {

class STRUCTURE;
class CHAIN;
class RESIDUE;
class ATOM;

class MODEL {

    friend STRUCTURE;
    friend CHAIN;
    friend RESIDUE;
    friend ATOM;

    friend std::ostream& operator<<(std::ostream& os, const STRUCTURE& structure);
    friend std::ostream& operator<<(std::ostream& os, const MODEL& md);
    friend std::ostream& operator<<(std::ostream& os, const CHAIN& res);
    friend std::ostream& operator<<(std::ostream& os, const RESIDUE& res);
    friend std::ostream& operator<<(std::ostream& os, const ATOM& a);

  private:
  
    std::map<std::string, boost::any> modelProperties;
    STRUCTURE* parentStructure_;

    boost::ptr_vector<CHAIN> chainVector;
    
    boost::ptr_vector<RESIDUE>  residueVector;
    boost::ptr_vector<ATOM>     atomVector;

  public:

    MODEL(int modelID,STRUCTURE* parentStructure);
    MODEL(int modelID);

    /*
    ~MODEL(){
    std::cerr << "Calling model destructor" << std::endl;
    }*/

    CHAIN& chain(std::string chainID);
    boost::ptr_vector<CHAIN>& chain();

    boost::ptr_vector<ATOM>&    atom();
    boost::ptr_vector<RESIDUE>& residue();
    
    bool isChain(std::string chainID);

    STRUCTURE& getParentStructure() const;
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_MDL_PROPERTY_T(type,Name)  GET_MDL_PROPERTY_R(type,Name)
    #define GET_MDL_PROPERTY_R(type,Name)  type getModel ## Name() const;
    #define GET_MDL_PROPERTY(r, data, tuple) GET_MDL_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
    
        MDL_PROPERTY_LOOP(GET_MDL_PROPERTY)
    
    #undef GET_MDL_PROPERTY_T
    #undef GET_MDL_PROPERTY_R
    #undef GET_MDL_PROPERTY
    
    ////////////////////////////////////////////////////////////////////
    
    #define SET_MDL_PROPERTY_T(Name,type,name)  SET_MDL_PROPERTY_R(Name,type,name)
    #define SET_MDL_PROPERTY_R(Name,type,name)  void setModel ## Name( type name ) ;
    #define SET_MDL_PROPERTY(r, data, tuple) SET_MDL_PROPERTY_T(PROPNAME_CAPS(tuple),PROPTYPE(tuple),PROPNAME(tuple))
        
        MDL_PROPERTY_LOOP(SET_MDL_PROPERTY)
        
    #undef SET_MDL_PROPERTY_T
    #undef SET_MDL_PROPERTY_R
    #undef SET_MDL_PROPERTY
        
    ////////////////////////////////////////////////////////////////////

    void addChain(std::string chainID);

    void addResidue(std::string chainID,
                    std::string resName,int resSeq,std::string iCode);

    void addAtom(std::string chainID,
                 std::string resName, int resSeq, std::string iCode,
                 int serial, std::string name);


};
}

#endif
