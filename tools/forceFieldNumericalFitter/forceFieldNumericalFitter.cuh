#include <iostream> 
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <tuple>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <cub/cub.cuh>

#include <eigen3/Eigen/Eigen>

#include <proteinManager.hpp>

#include "./integration3D/integrator3D.cuh"

#define DEBUG

namespace proteinManager{
namespace ffManager{
    
    template <class refPotType,class objPotType>
    class forceFieldFitter{
	
		private:
		
        refPotType refPot_;
        objPotType objPot_;
        
		proteinManager::STRUCTURE& refStruct_;
		proteinManager::STRUCTURE& objStruct_;
	
		std::vector<std::shared_ptr<proteinManager::ATOM>> refStructAtom;
		std::vector<std::shared_ptr<proteinManager::ATOM>> objStructAtom;
		
		//integral matrices
		
		Eigen::VectorXd vectorRef;
		Eigen::VectorXd vectorObj_fit;
        
        Eigen::VectorXd vectorObjPF;
		Eigen::MatrixXd matrixObjRef;
		Eigen::MatrixXd matrixObjObj;
        
        //Integration options
        
        enum integratorTypes {grid,grid_PF,MC};
        integratorTypes currentIntegrator;
        
        integrator::integrator3D_grid_ns::integrator3D_grid integGrid;
        
        integrator::integrator3D_MC_ns::integrator3D_MC integMC;
        int samplesPerIntegral_;
		
		public:
        
            struct GridIntegration{
                real cellSize = -1;
            };
            
            struct Grid_PF_Integration{
                std::string inputFilePath;
                real lFactor = 1;
                real fFactor = 1;
            };
        
            struct MonteCarloIntegration{
                real cutOff = -1;
                int pointsPerIntegral = -1;
                int samplesPerIntegral = -1;
            };
            
			forceFieldFitter(proteinManager::STRUCTURE& refStruct,
                             proteinManager::STRUCTURE& objStruct, 
                             GridIntegration intGrid,
                             refPotType refPot,
                             objPotType objPot
                             ):refStruct_(refStruct),objStruct_(objStruct),refPot_(refPot),objPot_(objPot){

                
				this->updateAtomVectors();
                
                #ifdef DEBUG
                std::cout << "Grid integrator selected" << std::endl;
                #endif
                std::tuple<real3,real3> box = this->computeBox();
                integGrid.init_fixedCellSize({std::get<0>(box).x,std::get<0>(box).y,std::get<0>(box).z},
                                             {std::get<1>(box).x,std::get<1>(box).y,std::get<1>(box).z},
                                             intGrid.cellSize);
                currentIntegrator = grid;
                
            }
            
            forceFieldFitter(proteinManager::STRUCTURE& refStruct,
                             proteinManager::STRUCTURE& objStruct, 
                             Grid_PF_Integration intGrid_PF):refStruct_(refStruct),objStruct_(objStruct){

                
				this->updateAtomVectors();
                
                #ifdef DEBUG
                std::cout << "Grid PF integrator selected" << std::endl;
                #endif
                integGrid.init_precompFunct(intGrid_PF.inputFilePath,intGrid_PF.lFactor,intGrid_PF.fFactor);
                currentIntegrator = grid_PF;
                
            }
            
            forceFieldFitter(proteinManager::STRUCTURE& refStruct,
                             proteinManager::STRUCTURE& objStruct, 
                             MonteCarloIntegration intMC,
                             refPotType refPot,
                             objPotType objPot
                             ):refStruct_(refStruct),objStruct_(objStruct),refPot_(refPot),objPot_(objPot){

                
				this->updateAtomVectors();
                
                #ifdef DEBUG
                std::cout << "MC integrator selected" << std::endl;
                #endif
                std::tuple<real3,real3> box = this->computeBox();
                integMC.init({std::get<0>(box).x,std::get<0>(box).y,std::get<0>(box).z},
                             {std::get<1>(box).x,std::get<1>(box).y,std::get<1>(box).z},
                             intMC.pointsPerIntegral);
                             
                samplesPerIntegral_ = intMC.samplesPerIntegral;
                currentIntegrator = MC;
                
            }
			
			void updateAtomVectors(){
				
				refStructAtom = std::vector<std::shared_ptr<proteinManager::ATOM>>();
				objStructAtom = std::vector<std::shared_ptr<proteinManager::ATOM>>();
				
				for(proteinManager::MODEL& md : refStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
								refStructAtom.push_back(std::make_shared<proteinManager::ATOM>(atm));
				}}}}
				
				for(proteinManager::MODEL& md : objStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
								objStructAtom.push_back(std::make_shared<proteinManager::ATOM>(atm));
				}}}}
				
			}
			
			std::tuple<real3,real3> computeBox(){
			
				real cutOff = std::max(refPot_.getCutOff(),objPot_.getCutOff());
				
                //Determining box
                
				real3 pos = refStruct_.model()[0].chain()[0].residue()[0].atom()[0].getAtomCoord();
				
				real3 boxMin = pos - cutOff;
				real3 boxMax = pos + cutOff;
				
				for(proteinManager::MODEL& md : refStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
					pos = atm.getAtomCoord();
					
					if(pos.x - cutOff < boxMin.x) boxMin.x = pos.x - cutOff;
					if(pos.y - cutOff < boxMin.y) boxMin.y = pos.y - cutOff;
					if(pos.z - cutOff < boxMin.z) boxMin.z = pos.z - cutOff;
					
					if(pos.x + cutOff > boxMax.x) boxMax.x = pos.x + cutOff;
					if(pos.y + cutOff > boxMax.y) boxMax.y = pos.y + cutOff;
					if(pos.z + cutOff > boxMax.z) boxMax.z = pos.z + cutOff;
						
				}}}}
				
                ////////////////////////////////////////////////////////
                
                return std::make_tuple(boxMin,boxMax);
			}
			
            
			void computeMatrices(){
				
                std::stringstream ss;
                
				real cutOff2 = std::max(refPot_.getCutOff(),objPot_.getCutOff());
                     cutOff2 = cutOff2*cutOff2;
                
                potProduct<objPotType,refPotType> potProductObjRef(objPot_,refPot_);
                potProduct<objPotType,objPotType> potProductObjObj(objPot_,objPot_);
				
				matrixObjRef.resize(objStructAtom.size(),refStructAtom.size());
				matrixObjObj.resize(objStructAtom.size(),objStructAtom.size());
                
				for(int i=0;i<objStructAtom.size();i++){
					std::cout << i << std::endl;
					for(int j=0;j<refStructAtom.size();j++){
                        real3 atmRef = refStructAtom[j]->getAtomCoord();
						real3 atmObj = objStructAtom[i]->getAtomCoord();
						
						real3 dst = atmRef - atmObj;
                        
                        potProductObjRef.pt1_.setParameters(*refStructAtom[j]);
                        
						if(dot(dst,dst) < cutOff2 and !potProductObjRef.pt1_.isNull()) {
                            
                            potProductObjRef.setParametersLinealFit(*objStructAtom[i],*refStructAtom[j]);
                            
                            if(currentIntegrator == grid ){
                                matrixObjRef(i,j) = integGrid.computeIntegral(potProductObjRef);
                            } else if (currentIntegrator == MC){
                                matrixObjRef(i,j) = integMC.computeIntegralAverage(potProductObjRef,samplesPerIntegral_).x;
                            } else {
                                ss.clear();
                                ss << "Selected integrator is not valid";
                                throw std::runtime_error(ss.str());
                            }
                                                        
						} else {
							matrixObjRef(i,j) = 0;
						}
					}
				}
				
				
				for(int i=0;i<objStructAtom.size();i++){
					std::cout << i << std::endl;
					for(int j=i;j<objStructAtom.size();j++){
						real3 atmObj1 = objStructAtom[i]->getAtomCoord();
						real3 atmObj2 = objStructAtom[j]->getAtomCoord();
						
						real3 dst = atmObj1 - atmObj2;
						
						if(dot(dst,dst) < cutOff2 ) {
                            
                            potProductObjObj.setParametersLinealFit(*objStructAtom[i],*objStructAtom[j]);
                            
                            if(currentIntegrator == grid ){
                                matrixObjObj(i,j) = integGrid.computeIntegral(potProductObjObj);
                            } else if (currentIntegrator == MC){
                                matrixObjObj(i,j) = integMC.computeIntegralAverage(potProductObjObj,samplesPerIntegral_).x;
                            } else {
                                ss.clear();
                                ss << "Selected integrator is not valid";
                                throw std::runtime_error(ss.str());
                            }

							matrixObjObj(j,i) = matrixObjObj(i,j);
                            
						} else {
							matrixObjObj(i,j) = 0;
							matrixObjObj(j,i) = 0;
						}
					}
				}
				
			}
            
            void computeNewParameters(){
                
                this->computeMatrices();
				
				vectorRef.resize(refStructAtom.size());
				vectorObj_fit.resize(objStructAtom.size());
				
				for(int i=0;i<refStructAtom.size();i++){
                    vectorRef[i] = refPot_.getInteractionParmLinealFit(*refStructAtom[i]);
				}
				
                vectorObj_fit = matrixObjObj.colPivHouseholderQr().solve(matrixObjRef*vectorRef);
                
                int i=0;
                for(proteinManager::MODEL& md : objStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
                                objPot_.setInteractionParmLinealFit(atm,vectorObj_fit[i]);
                                //Check
                                #ifdef DEBUG
                                std::cout << vectorObj_fit[i] << std::endl;
                                #endif
								i++;
				}}}}
				
				this->updateAtomVectors();
			}
            
            void computeNewParametersTotalChargeConstraint(){
                
                this->computeMatrices();
				
				matrixObjRef.conservativeResize(matrixObjRef.rows()+1,matrixObjRef.cols());
				matrixObjObj.conservativeResize(matrixObjObj.rows()+1,matrixObjObj.cols()+1);
				
				vectorRef.resize(refStructAtom.size());
				vectorObj_fit.resize(objStructAtom.size()+1);
				
				///////////////////////////////////////////////////////////////////
				
				for(int i=0;i<refStructAtom.size();i++){
                    vectorRef[i] = refPot_.getInteractionParmLinealFit(*refStructAtom[i]);
				}
				
				for(int i=0; i < matrixObjRef.cols(); i++){
					matrixObjRef(matrixObjRef.rows()-1,i) = 1;
				}
				
				//MatrixObjObj is always a symetric matrix
				for(int i=0; i < matrixObjObj.cols(); i++){
					matrixObjObj(matrixObjObj.rows()-1,i) = 1;
					matrixObjObj(i,matrixObjObj.cols()-1) = 1;
				}
				matrixObjObj(matrixObjObj.rows()-1,matrixObjObj.cols()-1) = 0;
				
				///////////////////////////////////////////////////////////////////
				
				vectorObj_fit = matrixObjObj.colPivHouseholderQr().solve(matrixObjRef*vectorRef);
                
                int i=0;
                for(proteinManager::MODEL& md : objStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
                                objPot_.setInteractionParmLinealFit(atm,vectorObj_fit[i]);
                                //Check
                                #ifdef DEBUG
                                std::cout << vectorObj_fit[i] << std::endl;
                                #endif
								i++;
				}}}}
				
				#ifdef DEBUG
				
				real totalChargeRef = 0;
				real totalChargeObj = 0;
				
				for(int i=0;i<vectorRef.size();i++){
					totalChargeRef += vectorRef[i];
				}
				
				for(int i=0;i<vectorObj_fit.size()-1;i++){
					totalChargeObj += vectorObj_fit[i];
				}
				
				std::cout << "Ref charge: " << totalChargeRef << " Obj charge: " << totalChargeObj << std::endl;
				
				real3 totalDipoleRef = {0,0,0};
				real3 totalDipoleObj = {0,0,0};
				
				for(proteinManager::MODEL& md : refStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
								totalDipoleRef += atm.getAtomCoord()*refPot_.getInteractionParmLinealFit(atm);
				}}}}
				
				for(proteinManager::MODEL& md : objStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
								totalDipoleObj += atm.getAtomCoord()*objPot_.getInteractionParmLinealFit(atm);
				}}}}
				
				std::cout << "Ref dipole: " << totalDipoleRef << " Obj dipole: " << totalDipoleObj << std::endl;
				
				#endif
				
				this->updateAtomVectors();
			}
            
            /*
            void computeVectorMatrix(){
                
                std::stringstream ss;
                
                if(currentIntegrator != grid_PF){
                    ss << "computeVectorMatrix only can be used with a grid_PF integrator";
                    throw std::runtime_error(ss.str());
                }
                
                real cutOff2 = pot_.getCutOff()*pot_.getCutOff();
                
                vector2PF.resize(objStructAtom.size());
                matrix22.resize(objStructAtom.size(),objStructAtom.size());
                
                for(int i=0;i<objStructAtom.size();i++){
					std::cout << i << std::endl;
                    
                    real3 atm1 = objStructAtom[i]->getAtomCoord();
                    
                    auto potP2PF = pot_.getPotentialProduct2PF({atm1.x,atm1.y,atm1.z});
                    
                    vector2PF[i] = integGrid.computeIntegral_PF(potP2PF);
                    
                    
                }
                
                for(int i=0;i<objStructAtom.size();i++){
					std::cout << i << std::endl;
					for(int j=i;j<objStructAtom.size();j++){
						real3 atm1 = objStructAtom[i]->getAtomCoord();
						real3 atm2 = objStructAtom[j]->getAtomCoord();
						
						real3 dst = atm1 - atm2;
						
						if(dot(dst,dst) < cutOff2) {
                            
                            auto potP22 = pot_.getPotentialProduct22({atm1.x,atm1.y,atm1.z},
                                                                     {atm2.x,atm2.y,atm2.z});

                            matrix22(i,j) = integGrid.computeIntegral(potP22);
							matrix22(j,i) = matrix22(i,j);
                            
						} else {
							matrix22(i,j) = 0;
							matrix22(j,i) = 0;
						}
					}
				}
                
                
            }
            
            void computeNewParametersUsingPF(){
                
                this->computeVectorMatrix();
                
                vector2_fit.resize(objStructAtom.size());
                
                vector2_fit = matrix22.colPivHouseholderQr().solve(vector2PF);
                
                int i=0;
                for(proteinManager::MODEL& md : objStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
                                potential::setInteractionParm(atm,vector2_fit[i]);
                                //Check
                                #ifdef DEBUG
                                std::cout << vector2_fit[i] << std::endl;
                                #endif
								i++;
				}}}}
				
				this->updateAtomVectors();
            }
			
            */
			
		
    };
}
}
