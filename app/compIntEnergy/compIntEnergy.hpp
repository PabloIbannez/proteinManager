#ifndef COMP_INT_ENERGY_HPP
#define COMP_INT_ENERGY_HPP

#include <proteinManager/proteinManager.hpp>

#include "potential.hpp"

namespace proteinManager{
namespace compInt{
    
    template <class potential>
    real compIntEnergy(MODEL& mdl1,potential pot){
        
        real energy = 0;
        
        std::vector<std::shared_ptr<proteinManager::ATOM>> model1_Atom;
        
        for(proteinManager::CHAIN& ch : mdl1.chain())  {
        for(proteinManager::RESIDUE& res : ch.residue())    {
        for(proteinManager::ATOM& atm : res.atom())         {
                        model1_Atom.push_back(std::make_shared<proteinManager::ATOM>(atm));
        }}}
        
        for(int i = 0    ; i < model1_Atom.size() ; i++){
            energy += pot.energy(*model1_Atom[i]);
        }
        
        return energy;
    }
    
    template <class potential>
    real compIntEnergy(MODEL& mdl1,MODEL& mdl2,potential pot){
        
        real energy = 0;
        
        std::vector<std::shared_ptr<proteinManager::ATOM>> model1_Atom;
        std::vector<std::shared_ptr<proteinManager::ATOM>> model2_Atom;
        
        for(proteinManager::CHAIN& ch : mdl1.chain())  {
        for(proteinManager::RESIDUE& res : ch.residue())    {
        for(proteinManager::ATOM& atm : res.atom())         {
                        model1_Atom.push_back(std::make_shared<proteinManager::ATOM>(atm));
        }}}
        
        for(proteinManager::CHAIN& ch : mdl2.chain())  {
        for(proteinManager::RESIDUE& res : ch.residue())    {
        for(proteinManager::ATOM& atm : res.atom())         {
                        model2_Atom.push_back(std::make_shared<proteinManager::ATOM>(atm));
        }}}
        
        for(int i = 0    ; i < model1_Atom.size() ; i++){
            energy += pot.energy(*model1_Atom[i]);
        }
        
        for(int j = 0 ; j < model2_Atom.size() ; j++){
            energy += pot.energy(*model2_Atom[j]);
        }
        
        for(int i = 0    ; i < model1_Atom.size() ; i++){
        for(int j = i +1 ; j < model2_Atom.size() ; j++){
            energy += pot.energy(*model1_Atom[i],*model2_Atom[j]);
        }}
        
        return energy;
    }
    
}}


#endif
