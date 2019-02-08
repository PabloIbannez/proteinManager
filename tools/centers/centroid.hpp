/* Pablo Ibáñez Freire, pablo.ibannez@uam.es */

/*
   Set of functions for calculating the centroid of the different entities.
*/

#include <proteinManager/proteinManager.hpp>

namespace proteinManager{
    

    real3 computeCentroid(RESIDUE& res){
        int atmCount = 0;
        real3 centroid = {0,0,0};
        for(ATOM& atm : res.atom()){
            centroid += atm.getAtomCoord();
            atmCount ++;
        }
        return centroid/atmCount;
    }
    
    real3 computeCentroid(CHAIN& ch){
        int atmCount = 0;
        real3 centroid = {0,0,0};
        for(RESIDUE& res : ch.residue()){
        for(ATOM& atm : res.atom()){
            centroid += atm.getAtomCoord();
            atmCount ++;
        }}
        return centroid/atmCount;
    }

    real3 computeCentroid(MODEL& mdl){
        int atmCount = 0;
        real3 centroid = {0,0,0};
        for(CHAIN& ch : mdl.chain()){
        for(RESIDUE& res : ch.residue()){
        for(ATOM& atm : res.atom()){
            centroid += atm.getAtomCoord();
            atmCount ++;
        }}}
        return centroid/atmCount;
    }
    
    real3 computeCentroid(STRUCTURE& st){
        int atmCount = 0;
        real3 centroid = {0,0,0};
        for(MODEL& mdl : st.model()){
        for(CHAIN& ch : mdl.chain()){
        for(RESIDUE& res : ch.residue()){
        for(ATOM& atm : res.atom()){
            centroid += atm.getAtomCoord();
            atmCount ++;
        }}}}
        return centroid/atmCount;
    }
}
