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

#include <proteinManager/proteinManager.hpp>

#include "./integration3D/integrator3D.cuh"
#include "./potential.cuh"

namespace proteinManager{
namespace ff_fitting{
    
    template <class refPotType,class objPotType>
    class forceFieldFitter{
	
		private:
		
        refPotType refPot_;
        objPotType objPot_;
        
		proteinManager::MODEL& refModel_;
		proteinManager::MODEL& objModel_;
	
		std::vector<std::shared_ptr<proteinManager::ATOM>> refModelAtom;
		std::vector<std::shared_ptr<proteinManager::ATOM>> objModelAtom;
		
		//integral matrices
		
		Eigen::VectorXd vectorRef;
		Eigen::VectorXd vectorObj_fit;
        
        Eigen::VectorXd vectorObjPF;
		Eigen::MatrixXd matrixObjRef;
		Eigen::MatrixXd matrixObjObj;
        
        //Integration options
        
        enum integratorTypes {grid,MC};
        integratorTypes currentIntegrator;
        
        integrator::integrator3D_grid_ns::integrator3D_grid integGrid;
        
        integrator::integrator3D_MC_ns::integrator3D_MC integMC;
        int samplesPerIntegral_;
		
		public:
        
            struct GridIntegration{
                real cellSize = -1;
            };
        
            struct MonteCarloIntegration{
                real cutOff = -1;
                int pointsPerIntegral = -1;
                int samplesPerIntegral = -1;
            };
            
			forceFieldFitter(proteinManager::MODEL& refModel,
                             proteinManager::MODEL& objModel, 
                             GridIntegration intGrid,
                             refPotType refPot,
                             objPotType objPot
                             ):refModel_(refModel),objModel_(objModel),refPot_(refPot),objPot_(objPot){

                
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
            
            forceFieldFitter(proteinManager::MODEL& refModel,
                             proteinManager::MODEL& objModel, 
                             MonteCarloIntegration intMC,
                             refPotType refPot,
                             objPotType objPot
                             ):refModel_(refModel),objModel_(objModel),refPot_(refPot),objPot_(objPot){

                
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
				
				refModelAtom = std::vector<std::shared_ptr<proteinManager::ATOM>>();
				objModelAtom = std::vector<std::shared_ptr<proteinManager::ATOM>>();
				
				for(proteinManager::CHAIN& ch : refModel_.chain())  {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
								refModelAtom.push_back(std::make_shared<proteinManager::ATOM>(atm));
				}}}
				
				for(proteinManager::CHAIN& ch : objModel_.chain())  {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
								objModelAtom.push_back(std::make_shared<proteinManager::ATOM>(atm));
				}}}
				
			}
			
			std::tuple<real3,real3> computeBox(){
			
				real cutOff = std::max(refPot_.getCutOff(),objPot_.getCutOff());
				
                //Determining box
                
				real3 pos = refModel_.chain()[0].residue()[0].atom()[0].getAtomCoord();
				
				real3 boxMin = pos - cutOff;
				real3 boxMax = pos + cutOff;
				
				for(proteinManager::CHAIN& ch : refModel_.chain())  {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
					pos = atm.getAtomCoord();
					
					if(pos.x - cutOff < boxMin.x) boxMin.x = pos.x - cutOff;
					if(pos.y - cutOff < boxMin.y) boxMin.y = pos.y - cutOff;
					if(pos.z - cutOff < boxMin.z) boxMin.z = pos.z - cutOff;
					
					if(pos.x + cutOff > boxMax.x) boxMax.x = pos.x + cutOff;
					if(pos.y + cutOff > boxMax.y) boxMax.y = pos.y + cutOff;
					if(pos.z + cutOff > boxMax.z) boxMax.z = pos.z + cutOff;
						
				}}}
				
                ////////////////////////////////////////////////////////
                
                return std::make_tuple(boxMin,boxMax);
			}
			
            
			void computeMatrices(){
				
                std::stringstream ss;
                
				real cutOff2 = std::max(refPot_.getCutOff(),objPot_.getCutOff());
                     cutOff2 = cutOff2*cutOff2;
                
                potProduct<objPotType,refPotType> potProductObjRef(objPot_,refPot_);
                potProduct<objPotType,objPotType> potProductObjObj(objPot_,objPot_);
				
				matrixObjRef.resize(objModelAtom.size(),refModelAtom.size());
				matrixObjObj.resize(objModelAtom.size(),objModelAtom.size());
                
				for(int i=0;i<objModelAtom.size();i++){
					std::cout << i << std::endl;
					for(int j=0;j<refModelAtom.size();j++){
                        real3 atmRef = refModelAtom[j]->getAtomCoord();
						real3 atmObj = objModelAtom[i]->getAtomCoord();
						
						real3 dst = atmRef - atmObj;
                        
                        potProductObjRef.pt1_.setParameters(*refModelAtom[j]);
                        
						if(dot(dst,dst) < cutOff2 and !potProductObjRef.pt1_.isNull()) {
                            
                            potProductObjRef.setParametersLinealFit(*objModelAtom[i],*refModelAtom[j]);
                            
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
				
				
				for(int i=0;i<objModelAtom.size();i++){
					std::cout << i << std::endl;
					for(int j=i;j<objModelAtom.size();j++){
						real3 atmObj1 = objModelAtom[i]->getAtomCoord();
						real3 atmObj2 = objModelAtom[j]->getAtomCoord();
						
						real3 dst = atmObj1 - atmObj2;
						
						if(dot(dst,dst) < cutOff2 ) {
                            
                            potProductObjObj.setParametersLinealFit(*objModelAtom[i],*objModelAtom[j]);
                            
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
				
				vectorRef.resize(refModelAtom.size());
				vectorObj_fit.resize(objModelAtom.size());
				
				for(int i=0;i<refModelAtom.size();i++){
                    vectorRef[i] = refPot_.getInteractionParmLinealFit(*refModelAtom[i]);
				}
				
                vectorObj_fit = matrixObjObj.colPivHouseholderQr().solve(matrixObjRef*vectorRef);
                
                int i=0;
				for(proteinManager::CHAIN& ch : objModel_.chain())  {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
                                objPot_.setInteractionParmLinealFit(atm,vectorObj_fit[i]);
                                //Check
                                #ifdef DEBUG
                                std::cout << vectorObj_fit[i] << std::endl;
                                #endif
								i++;
				}}}
				
				this->updateAtomVectors();
			}
            
            void computeNewParametersTotalChargeConstraint(){
                
                this->computeMatrices();
				
				matrixObjRef.conservativeResize(matrixObjRef.rows()+1,matrixObjRef.cols());
				matrixObjObj.conservativeResize(matrixObjObj.rows()+1,matrixObjObj.cols()+1);
				
				vectorRef.resize(refModelAtom.size());
				vectorObj_fit.resize(objModelAtom.size()+1);
				
				///////////////////////////////////////////////////////////////////
				
				for(int i=0;i<refModelAtom.size();i++){
                    vectorRef[i] = refPot_.getInteractionParmLinealFit(*refModelAtom[i]);
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
				for(proteinManager::CHAIN& ch : objModel_.chain())  {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
                                objPot_.setInteractionParmLinealFit(atm,vectorObj_fit[i]);
                                //Check
                                #ifdef DEBUG
                                std::cout << vectorObj_fit[i] << std::endl;
                                #endif
								i++;
				}}}
				
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
                    
                    for(proteinManager::CHAIN& ch : refModel_.chain())  {
                    for(proteinManager::RESIDUE& res : ch.residue())    {
                    for(proteinManager::ATOM& atm : res.atom())         {
                                    totalDipoleRef += atm.getAtomCoord()*refPot_.getInteractionParmLinealFit(atm);
                    }}}
                    
                    for(proteinManager::CHAIN& ch : objModel_.chain())  {
                    for(proteinManager::RESIDUE& res : ch.residue())    {
                    for(proteinManager::ATOM& atm : res.atom())         {
                                    totalDipoleObj += atm.getAtomCoord()*objPot_.getInteractionParmLinealFit(atm);
                    }}}
                    
                    std::cout << "Ref dipole: " << totalDipoleRef << " Obj dipole: " << totalDipoleObj << std::endl;
				
				#endif
				
				this->updateAtomVectors();
			}
			
		
    };
    
    template <class objPotType>
    class forceFieldFitter_PrecomputedField{
        
        private:
        
            objPotType objPot_;
            
            proteinManager::STRUCTURE& objStruct_;
            
            std::vector<std::shared_ptr<proteinManager::ATOM>> objStructAtom;
            
            //integral matrices
            
            Eigen::VectorXd vectorObj_fit;
            
            Eigen::VectorXd vectorObjPF;
            Eigen::MatrixXd matrixObjObj;
            
            integrator::integrator3D_grid_ns::integrator3D_grid integGrid;
            
            
        public:
        
            struct GridIntegrationPF{
                std::string precompFieldFilePath;
            };
            
            forceFieldFitter_PrecomputedField(proteinManager::STRUCTURE& objStruct, 
                                              GridIntegrationPF grPF,
                                              objPotType objPot
                                             ):objStruct_(objStruct),objPot_(objPot){

                
				this->updateAtomVectors();
                
                #ifdef DEBUG
                std::cout << "gridPF integrator selected" << std::endl;
                #endif
                
                integGrid.init_precompFunct(grPF.precompFieldFilePath);
            }
			
			void updateAtomVectors(){
                
				objStructAtom = std::vector<std::shared_ptr<proteinManager::ATOM>>();
				
				for(proteinManager::MODEL& md : objStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
								objStructAtom.push_back(std::make_shared<proteinManager::ATOM>(atm));
				}}}}
				
			}
            
            void computeMatrixVector(){
                
                real cutOff2 = objPot_.getCutOff();
                     cutOff2 = cutOff2*cutOff2;
                
                vectorObjPF.resize(objStructAtom.size());
                matrixObjObj.resize(objStructAtom.size(),objStructAtom.size());
                
                potPF_Product<objPotType> potProductPFObj(objPot_);
                potProduct<objPotType,objPotType> potProductObjObj(objPot_,objPot_);
                
                for(int i=0;i<objStructAtom.size();i++){
					std::cout << i << std::endl;
                    
                    potProductPFObj.setParametersLinealFit(*objStructAtom[i]);
                    
                    vectorObjPF[i] = integGrid.computeIntegral_PF(potProductPFObj);
                }
                
                for(int i=0;i<objStructAtom.size();i++){
					std::cout << i << std::endl;
					for(int j=i;j<objStructAtom.size();j++){
						real3 atmObj1 = objStructAtom[i]->getAtomCoord();
						real3 atmObj2 = objStructAtom[j]->getAtomCoord();
						
						real3 dst = atmObj1 - atmObj2;
						
						if(dot(dst,dst) < cutOff2 ) {
                            
                            potProductObjObj.setParametersLinealFit(*objStructAtom[i],*objStructAtom[j]);
                            
                            matrixObjObj(i,j) = integGrid.computeIntegral_PF(potProductPFObj);
							matrixObjObj(j,i) = matrixObjObj(i,j);
                            
						} else {
							matrixObjObj(i,j) = 0;
							matrixObjObj(j,i) = 0;
						}
					}
				}
            }
            
            void computeNewParameters(){
                
                this->computeMatrixVector();
                
                vectorObj_fit.resize(objStructAtom.size());
                
                vectorObj_fit = matrixObjObj.colPivHouseholderQr().solve(vectorObjPF);
                
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
            
            
        
    };
    
}
}
