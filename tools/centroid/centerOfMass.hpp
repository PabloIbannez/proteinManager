/* Pablo Ibáñez Freire, pablo.ibannez@uam.es */

/*
   Set of functions for calculating the center of mass of the different entities.
*/
#include <proteinManager/proteinManager.hpp>

namespace proteinManager{

    real3 computeCenterOfMass(RESIDUE& res){
        real totalMass = 0;
        real3 centerOfMass = {0,0,0};
        for(ATOM& atm : res.atom()){
            centerOfMass += atm.getAtomCoord()*atm.getAtomMass();
            totalMass += atm.getAtomMass();
        }
        return centerOfMass/totalMass;
    }
    
    real3 computeCenterOfMass(CHAIN& ch){
        real totalMass = 0;
        real3 centerOfMass = {0,0,0};
        for(RESIDUE& res : ch.residue()){
        for(ATOM& atm : res.atom()){
            centerOfMass += atm.getAtomCoord()*atm.getAtomMass();
            totalMass += atm.getAtomMass();
        }}
        return centerOfMass/totalMass;
    }
    
    real3 computeCenterOfMass(MODEL& mdl){
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
    
    real3 computeCenterOfMass(STRUCTURE& st){
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
