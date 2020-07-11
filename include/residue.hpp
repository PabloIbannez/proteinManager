/*Pablo Ibáñez Freire, pablo.ibannez@uam.es*/

#ifndef RESIDUE_HPP
#define RESIDUE_HPP

#include <vector>

#include "proteinManager.hpp"

namespace proteinManager {

class STRUCTURE;
class MODEL;
class CHAIN;
class ATOM;

class RESIDUE {

    friend STRUCTURE;
    friend MODEL;
    friend CHAIN;
    friend ATOM;

    friend std::ostream& operator<<(std::ostream& os, const STRUCTURE& structure);
    friend std::ostream& operator<<(std::ostream& os, const MODEL& md);
    friend std::ostream& operator<<(std::ostream& os, const CHAIN& res);
    friend std::ostream& operator<<(std::ostream& os, const RESIDUE& res);
    friend std::ostream& operator<<(std::ostream& os, const ATOM& a);

  private:
    
    std::map<std::string, boost::any> residueProperties;
    
    CHAIN* parentChain_;

    boost::ptr_vector<ATOM> atomVector;
    std::map<std::string,int> atomMap;

  public:

    RESIDUE(std::string resName,int resSeq,std::string iCode,CHAIN* parentChain);
    RESIDUE(std::string resName, int resSeq, std::string iCode);

    /*
     ~RESIDUE(){
    std::cerr << "Calling residue destructor" << std::endl;
    }*/
    
    ATOM& atom(std::string atomName);
    boost::ptr_vector<ATOM>& atom();
    
    bool isAtom(std::string atomName);
    
    STRUCTURE& getParentStructure() const;
    MODEL& getParentModel() const ;
    CHAIN& getParentChain() const ;
    
    ////////////////////////////////////////////////////////////////////

    #define GET_MDL_PROPERTY_T(type,Name)  GET_MDL_PROPERTY_R(type,Name)
    #define GET_MDL_PROPERTY_R(type,Name)  type getModel ## Name() const;
    #define GET_MDL_PROPERTY(r, data, tuple) GET_MDL_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
    
        MDL_PROPERTY_LOOP(GET_MDL_PROPERTY)
    
    #undef GET_MDL_PROPERTY_T
    #undef GET_MDL_PROPERTY_R
    #undef GET_MDL_PROPERTY
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_CHAIN_PROPERTY_T(type,Name)  GET_CHAIN_PROPERTY_R(type,Name)
    #define GET_CHAIN_PROPERTY_R(type,Name)  type getChain ## Name() const;
    #define GET_CHAIN_PROPERTY(r, data, tuple) GET_CHAIN_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
    
        CHAIN_PROPERTY_LOOP(GET_CHAIN_PROPERTY)
    
    #undef GET_CHAIN_PROPERTY_T
    #undef GET_CHAIN_PROPERTY_R
    #undef GET_CHAIN_PROPERTY
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_RES_PROPERTY_T(type,Name)  GET_RES_PROPERTY_R(type,Name)
    #define GET_RES_PROPERTY_R(type,Name)  type getRes ## Name() const;
    #define GET_RES_PROPERTY(r, data, tuple) GET_RES_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
    
        RES_PROPERTY_LOOP(GET_RES_PROPERTY)
    
    #undef GET_RES_PROPERTY_T
    #undef GET_RES_PROPERTY_R
    #undef GET_RES_PROPERTY
    
    ////////////////////////////////////////////////////////////////////
    
    #define SET_RES_PROPERTY_T(Name,type,name)  SET_RES_PROPERTY_R(Name,type,name)
    #define SET_RES_PROPERTY_R(Name,type,name)  void setRes ## Name( type name ) ;
    #define SET_RES_PROPERTY(r, data, tuple) SET_RES_PROPERTY_T(PROPNAME_CAPS(tuple),PROPTYPE(tuple),PROPNAME(tuple))
    
        RES_PROPERTY_LOOP(SET_RES_PROPERTY)
    
    #undef SET_RES_PROPERTY_T
    #undef SET_RES_PROPERTY_R
    #undef SET_RES_PROPERTY
    
    ////////////////////////////////////////////////////////////////////
    
    void addAtom(int serial,std::string name);
};
}

#endif
