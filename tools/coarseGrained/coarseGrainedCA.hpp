#ifndef COARSE_GRAINED_CA_HPP
#define COARSE_GRAINED_CA_HPP

#include <proteinManager/proteinManager.hpp>

#include "../../tools/misc/misc.hpp"

namespace proteinManager{
namespace coarseGrainedCA{

    void coarseGrainedCA(STRUCTURE& structureIn,STRUCTURE& structureOut){
        
        for(MODEL&   mdl : structureIn.model()){
            structureOut.addModel(mdl.getModelId());
        for(CHAIN&   ch  : mdl.chain()       ){
            structureOut.model(mdl.getModelId()).addChain(ch.getChainId());
        for(RESIDUE& res : ch.residue()      ){
            structureOut.model(mdl.getModelId()).chain(ch.getChainId()).addResidue(res.getResName(),res.getResSeq(),res.getResInsCode());
        for(ATOM&   atm : res.atom()         ){

            if(atm.getAtomName() == "CA"){
                structureOut.model(mdl.getModelId()).chain(ch.getChainId()).residue(res.getResSeq()).addAtom(atm.getAtomSerial(),atm.getAtomName());
                structureOut.model(mdl.getModelId()).chain(ch.getChainId()).residue(res.getResSeq()).atom(atm.getAtomName()).setAtomCoord(atm.getAtomCoord());
            }

        }}}}
        
        structureOut.renumber();
        
        misc::fixHIS(structureOut);
        misc::fixCYS(structureOut);
        
        misc::setAtomsChargeEqualToResCharge(structureOut);
        misc::setAtomsNameEqualToResName(structureOut);

    }
    
    void coarseGrainedCA_SASA(STRUCTURE& structureIn,STRUCTURE& structureOut){
        
        for(MODEL&   mdl : structureIn.model()){
            structureOut.addModel(mdl.getModelId());
        for(CHAIN&   ch  : mdl.chain()       ){
            structureOut.model(mdl.getModelId()).addChain(ch.getChainId());
        for(RESIDUE& res : ch.residue()      ){
            structureOut.model(mdl.getModelId()).chain(ch.getChainId()).addResidue(res.getResName(),res.getResSeq(),res.getResInsCode());
        for(ATOM&   atm : res.atom()         ){

            if(atm.getAtomName() == "CA"){
                structureOut.model(mdl.getModelId()).chain(ch.getChainId()).residue(res.getResSeq()).addAtom(atm.getAtomSerial(),atm.getAtomName());
                structureOut.model(mdl.getModelId()).chain(ch.getChainId()).residue(res.getResSeq()).atom(atm.getAtomName()).setAtomCoord(atm.getAtomCoord());
                
                real sasa=0;
                for(ATOM&   atmSasa : res.atom()         ){
                    sasa+=atmSasa.getAtomSASA();
                }
                
                structureOut.model(mdl.getModelId()).chain(ch.getChainId()).residue(res.getResSeq()).atom(atm.getAtomName()).setAtomSASA(sasa);
            }

        }}}}
        
        structureOut.renumber();
        
        misc::fixHIS(structureOut);
        misc::fixCYS(structureOut);
        
        misc::setAtomsChargeEqualToResCharge(structureOut);
        misc::setAtomsNameEqualToResName(structureOut);

    }

}}

#endif
