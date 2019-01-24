/* Pablo Ibáñez Freire, pablo.ibannez@uam.es */

/*
   Set of functions for calculating the centroid of the different entities.
*/

#include <proteinManager/proteinManager.hpp>

namespace proteinManager{
    
    template <class T>
    real3 computeCentroid(T& in){ return {0,0,0};}
    
    template<>
    real3 computeCentroid<RESIDUE>(RESIDUE& res){
        int atmCount = 0;
        real3 centroid = {0,0,0};
        for(ATOM& atm : res.atom()){
            centroid += atm.getAtomCoord();
            atmCount ++;
        }
        return centroid/atmCount;
    }
    
    template<>
    real3 computeCentroid<CHAIN>(CHAIN& ch){
        int atmCount = 0;
        real3 centroid = {0,0,0};
        for(RESIDUE& res : ch.residue()){
        for(ATOM& atm : res.atom()){
            centroid += atm.getAtomCoord();
            atmCount ++;
        }}
        return centroid/atmCount;
    }
    
    template<>
    real3 computeCentroid<MODEL>(MODEL& mdl){
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
    
    template<>
    real3 computeCentroid<STRUCTURE>(STRUCTURE& st){
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
