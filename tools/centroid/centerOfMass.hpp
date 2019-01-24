/* Pablo Ibáñez Freire, pablo.ibannez@uam.es */

/*
   Set of functions for calculating the center of mass of the different entities.
*/
#include <proteinManager/proteinManager.hpp>

namespace proteinManager{
    
    template <class T>
    real3 computeCenterOfMass(T& in){ return {0,0,0};}
    
    template<>
    real3 computeCenterOfMass<RESIDUE>(RESIDUE& res){
        real totalMass = 0;
        real3 centerOfMass = {0,0,0};
        for(ATOM& atm : res.atom()){
            centerOfMass += atm.getAtomCoord()*atm.getAtomMass();
            totalMass += atm.getAtomMass();
        }
        return centerOfMass/totalMass;
    }
    
    template<>
    real3 computeCenterOfMass<CHAIN>(CHAIN& ch){
        real totalMass = 0;
        real3 centerOfMass = {0,0,0};
        for(RESIDUE& res : ch.residue()){
        for(ATOM& atm : res.atom()){
            centerOfMass += atm.getAtomCoord()*atm.getAtomMass();
            totalMass += atm.getAtomMass();
        }}
        return centerOfMass/totalMass;
    }
    
    template<>
    real3 computeCenterOfMass<MODEL>(MODEL& mdl){
        real totalMass = 0;
        real3 centerOfMass = {0,0,0};
        for(CHAIN& ch : mdl.chain()){
        for(RESIDUE& res : ch.residue()){
        for(ATOM& atm : res.atom()){
            centerOfMass += atm.getAtomCoord()*atm.getAtomMass();
            totalMass += atm.getAtomMass();
        }}}
        return centerOfMass/totalMass;
    }
    
    template<>
    real3 computeCenterOfMass<STRUCTURE>(STRUCTURE& st){
        real totalMass = 0;
        real3 centerOfMass = {0,0,0};
        for(MODEL& mdl : st.model()){
        for(CHAIN& ch : mdl.chain()){
        for(RESIDUE& res : ch.residue()){
        for(ATOM& atm : res.atom()){
            centerOfMass += atm.getAtomCoord()*atm.getAtomMass();
            totalMass += atm.getAtomMass();
        }}}}
        return centerOfMass/totalMass;
    }
    
    
    
}
