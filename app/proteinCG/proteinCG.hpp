#ifndef PROTEIN_CG_HPP
#define PROTEIN_CG_HPP

#include <proteinManager/proteinManager.hpp>

#include "../../tools/massesManager/massesManager.hpp"
#include "../../tools/forceFieldManager/forceFieldManager.hpp"

#include "../../tools/coarseGrained/coarseGrainedManager.hpp"
#include "../../tools/coarseGrained/coarseGrainedMappingSchemes.hpp"

#include "../../tools/geometricTransformations/geometricTransformations.hpp"

#include "../../tools/forceFieldNumericalFitter/forceFieldNumericalFitter.cuh"
#include "../../tools/forceFieldNumericalFitter/potential.cuh"

namespace proteinManager{

class coarseGrainedTop{
    
    private:
        
        std::shared_ptr<STRUCTURE> strIn_;
        std::shared_ptr<STRUCTURE> strOut_;
        
        //Data
        
        std::string massesData = "../../tools/massesManager/massesData/atomMasses.dat";
        std::string ffData     = "../../tools/forceFieldManager/forceFieldModels/gromos.ff";
        
        std::string cgData_aminoAcid2bead = "../../tools/coarseGrained/coarseGrainedModels/RES2BEAD/aminoAcid2bead_RES2BEAD.map";
        std::string cgData_bead2atom      = "../../tools/coarseGrained/coarseGrainedModels/RES2BEAD/bead2atom_RES2BEAD.map";
        
        
        massesManager massesM;
        ffManager::forceFieldManager ffM;
        coarseGrainedManager::coarseGrainedGenerator cgM;
        
    public:
    
        coarseGrainedTop(){
            
            massesM.loadMassesData(massesData);
            ffM.loadForceFieldData(ffData);
            
            cgM.loadCGmodel(cgData_aminoAcid2bead, \
                            cgData_bead2atom);
            
        }
        
        void loadStructure(STRUCTURE& strIn){
            
            strIn_  = std::make_shared<STRUCTURE>(strIn);
            strOut_ = std::make_shared<STRUCTURE>();
            
            massesM.applyMassesData(*strIn_);
            ffM.applyForceFieldData(*strIn_);
            
            cgM.applyCoarseGrainedMap<coarseGrainedManager::coarseGrainedMappingSchemes::sasaFitting>(*strIn_,*strOut_);
            
        }
        
        STRUCTURE& getOutputStructure(){
            return *strOut_;
        }
        
        void loadChargesSurf(int modelId,std::string ChgSurfFilePath){
            
            MODEL& mdlIn = (*strOut_).model(modelId);
    
            struct chainRes{
                std::string chainId;
                int resSeq;
            };
            
            std::map<int,chainRes> oneChainMap;
            chainRes chRbuffer;
            
            int i = 1;
            for(CHAIN&   ch  : mdlIn.chain())  {
            for(RESIDUE& res : ch.residue())   {
                chRbuffer.chainId = ch.getChainId();
                chRbuffer.resSeq  = res.getResSeq();
                oneChainMap[i]=chRbuffer;
                i++;
            }}
            
            ////////////////////////////////////////////////////////////////////
            
            std::stringstream ss;
                        
            //Check if file exists
            std::ifstream inputFile(ChgSurfFilePath);
            if(!inputFile){
                ss.clear();
                ss << "File not found: " << ChgSurfFilePath;
                throw std::runtime_error(ss.str());
            }
            
            std::string line;
            
            int    resSeqBuffer;
            double chgBuffer;
            
            for(CHAIN&   ch  : mdlIn.chain())  {
            for(RESIDUE& res : ch.residue())   {
                
                if(res.atom().size() > 1){
                    ss.clear();
                    ss << "ERROR. Expected one atom per risude only";
                    throw std::runtime_error(ss.str());
                }
                res.atom()[0].setAtomCharge(0);
                res.atom()[0].setAtomSurf(false);
            }}
            
            while(std::getline(inputFile,line)){
                ss.clear();
                ss.str(line);
                
                ss >> resSeqBuffer >> chgBuffer;
                mdlIn.chain(oneChainMap[resSeqBuffer].chainId).residue(oneChainMap[resSeqBuffer].resSeq).atom()[0].setAtomCharge(chgBuffer);
                mdlIn.chain(oneChainMap[resSeqBuffer].chainId).residue(oneChainMap[resSeqBuffer].resSeq).atom()[0].setAtomSurf(true);
            }
        
        }
        
        void fit_C6(int modelId,real cutOff, real epsilon){
            
            std::cout << "C6 " << modelId << std::endl;
            
            MODEL& mdlIn  = (*strIn_).model(modelId);
            MODEL& mdlOut = (*strOut_).model(modelId);
            
            using potentialRef_C6 = ff_fitting::attractive6;
            using potentialObj_C6 = ff_fitting::attractive6;
            
            potentialRef_C6 potRef_C6(cutOff,epsilon);
            potentialObj_C6 potObj_C6(cutOff,epsilon);
            
            ff_fitting::forceFieldFitter<potentialRef_C6,potentialObj_C6>::MonteCarloIntegration mc_C6;
            mc_C6.pointsPerIntegral = 100000;
            mc_C6.samplesPerIntegral = 10;
            
            ff_fitting::forceFieldFitter<potentialRef_C6,potentialObj_C6> C6_fitter(mdlIn,mdlOut,mc_C6,potRef_C6,potObj_C6);
            C6_fitter.computeNewParametersTotalChargeConstraint();
        }
        
        void fit_C12(int modelId,real cutOff, real epsilon){
            
            std::cout << "C12 " << modelId << std::endl;
            
            MODEL& mdlIn  = (*strIn_).model(modelId);
            MODEL& mdlOut = (*strOut_).model(modelId);
        
            using potentialRef_C12 = ff_fitting::repulsive12;
            using potentialObj_C12 = ff_fitting::repulsive12;
            
            potentialRef_C12 potRef_C12(cutOff,epsilon);
            potentialObj_C12 potObj_C12(cutOff,epsilon);
            
            ff_fitting::forceFieldFitter<potentialRef_C12,potentialObj_C12>::MonteCarloIntegration mc_C12;
            mc_C12.pointsPerIntegral = 100000;
            mc_C12.samplesPerIntegral = 10;
            
            ff_fitting::forceFieldFitter<potentialRef_C12,potentialObj_C12> C12_fitter(mdlIn,mdlOut,mc_C12,potRef_C12,potObj_C12);
            C12_fitter.computeNewParametersTotalChargeConstraint();
            
        }
        
        std::shared_ptr<STRUCTURE> getStructureOut(){
            return strOut_;
        }
        
        void output(std::ostream& out){
            
            int atomCount = 0;
            for(MODEL&   mdl : (*strOut_).model()){
            for(CHAIN&   ch  : mdl.chain()  ){
            for(RESIDUE& res : ch.residue() ){
            for(ATOM&    atm : res.atom()   ){
                out << std::right
                    << std::fixed
                    << std::setprecision(4)
                    << std::setw(10)                   
                    << atomCount                       << " " 
                    << std::setw(5)                    
                    << res.getResName()                << " " 
                    << std::setw(5)                    
                    << mdl.getModelId()                << " " 
                    << std::setw(8)                    
                    << atm.getAtomSASA()               << " " 
                    << std::setw(1)                    
                    << atm.getAtomSurf()               << " " 
                    << std::setw(10)                   
                    << atm.getAtomMass()               << " " 
                    << std::setw(10)                   
                    << atm.getAtomCoord()              << " " 
                    << std::setw(10)                   
                    << std::setprecision(5)            
                    << atm.getAtomC12()                << " " 
                    << std::setw(10)                   
                    << atm.getAtomC6()                 << " "
                    << std::setw(10)                   
                    << atm.getAtomCharge()             << " "
                    << std::setw(10)                   
                    << atm.getAtomSolvE()              << std::endl;
                
                atomCount ++;
            }}}}
            
            
        }
};
}

#endif
