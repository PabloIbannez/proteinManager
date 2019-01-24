#ifndef ANALYTIC_CHARGE_FITTING_HPP
#define ANALYTIC_CHARGE_FITTING_HPP

//#define DEBUG

#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <proteinManager/proteinManager.hpp>

#include <eigen3/Eigen/Eigen>

namespace proteinManager {
namespace chargeFitting{
    
    class analyticChargeFitting{
        
        private:
            
            proteinManager::STRUCTURE& refStruct_;
            proteinManager::STRUCTURE& objStruct_;
	
            std::vector<std::shared_ptr<proteinManager::ATOM>> refStructAtom;
            std::vector<std::shared_ptr<proteinManager::ATOM>> objStructAtom;
		
            //matrices
		
            Eigen::VectorXd vectorRef;
            Eigen::VectorXd vectorObj;
            
            Eigen::MatrixXd matrixRefRef;
            Eigen::MatrixXd matrixObjRef;
            Eigen::MatrixXd matrixObjObj;
            
            //position data
            
            std::vector<real3> gradientObj;
            
            //const
            
            real lbda = 0.01;
            
            //
            
            real atmContribution;
            bool atmContributionComp = false;
            
            
        public:
        
            analyticChargeFitting(proteinManager::STRUCTURE& refStruct,
                                  proteinManager::STRUCTURE& objStruct):refStruct_(refStruct),objStruct_(objStruct)
            {

				this->updateAtomVectors();
                
            }
            
            real computeResidue(){
                
                matrixRefRef.resize(refStructAtom.size(),refStructAtom.size());
                matrixObjRef.resize(objStructAtom.size(),refStructAtom.size());
				matrixObjObj.resize(objStructAtom.size(),objStructAtom.size());
                
                vectorRef.resize(refStructAtom.size());
                vectorObj.resize(objStructAtom.size());
                
                for(int i=0;i<refStructAtom.size();i++){
                    vectorRef(i) = refStructAtom[i]->getAtomCharge();
                    //std::cout << vectorRef(i) << std::endl;
                }
                
                for(int i=0;i<objStructAtom.size();i++){
                    vectorObj(i) = objStructAtom[i]->getAtomCharge();
                    //std::cout << vectorObj(i) << std::endl;
                }
                
                if(!atmContributionComp){
                    
                    for(int i=0;i<refStructAtom.size();i++){
                        for(int j=0;j<refStructAtom.size();j++){
                            real3 dst = refStructAtom[i]->getAtomCoord() - refStructAtom[j]->getAtomCoord();
                            matrixRefRef(i,j) = real(-2*M_PI)*(sqrt(dot(dst,dst)));
                        }
                    }
                    
                    atmContribution = vectorRef.transpose()*matrixRefRef*vectorRef;
                    
                    atmContributionComp = true;
                }
                
                for(int i=0;i<objStructAtom.size();i++){
					for(int j=0;j<refStructAtom.size();j++){
                        real3 dst = objStructAtom[i]->getAtomCoord() - refStructAtom[j]->getAtomCoord();
                        matrixObjRef(i,j) = real(-2*M_PI)*(sqrt(dot(dst,dst)));
                    }
                }
                
                for(int i=0;i<objStructAtom.size();i++){
					for(int j=0;j<objStructAtom.size();j++){
                        real3 dst = objStructAtom[i]->getAtomCoord() - objStructAtom[j]->getAtomCoord();
                        matrixObjObj(i,j) = real(-2*M_PI)*(sqrt(dot(dst,dst)));
                    }
                }
                
                
                
                return atmContribution+real(-2)*vectorObj.transpose()*matrixObjRef*vectorRef+vectorObj.transpose()*matrixObjObj*vectorObj;
                
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
            
            void fitCharges(){
                
                matrixObjRef.resize(objStructAtom.size(),refStructAtom.size());
				matrixObjObj.resize(objStructAtom.size(),objStructAtom.size());
                
                vectorRef.resize(refStructAtom.size());
                vectorObj.resize(objStructAtom.size());
                
                
                #ifdef DEBUG
                    
                    std::cout << "MatrixObjRef size: " << matrixObjRef.size() << std::endl;
                    std::cout << "MatrixObjObj size: " << matrixObjObj.size() << std::endl;

                    
                    std::cout << "vectorRef size: " << vectorRef.size() << std::endl;
                    std::cout << "vectorObj size: " << vectorObj.size() << std::endl;
                    
                #endif
                
                for(int i=0;i<objStructAtom.size();i++){
					//std::cout << i << std::endl;
					for(int j=0;j<refStructAtom.size();j++){
                        if(i == objStructAtom.size()-1){
                            
                            matrixObjRef(i,j) = 1;
                            
                        } else {
                            
                            real3 atmObj_1 = objStructAtom[i]->getAtomCoord();
                            real3 atmObj_2 = objStructAtom[i+1]->getAtomCoord();
                            real3 atmRef   = refStructAtom[j]->getAtomCoord();
						
                            real3 dst1 = atmObj_1 - atmRef;
                            real3 dst2 = atmObj_2 - atmRef;
                            
                            matrixObjRef(i,j) = sqrt(dot(dst1,dst1))-sqrt(dot(dst2,dst2));
                            
                        }
                    }
                }
                
                for(int i=0;i<objStructAtom.size();i++){
					//std::cout << i << std::endl;
					for(int j=0;j<objStructAtom.size();j++){
                        if(i == objStructAtom.size()-1){
                            
                            matrixObjObj(i,j) = 1;
                            
                        } else {
                            
                            real3 atmObj_1 = objStructAtom[i]->getAtomCoord();
                            real3 atmObj_2 = objStructAtom[i+1]->getAtomCoord();
                            real3 atmObj   = objStructAtom[j]->getAtomCoord();
						
                            real3 dst1 = atmObj_1 - atmObj;
                            real3 dst2 = atmObj_2 - atmObj;
                            
                            matrixObjObj(i,j) = sqrt(dot(dst1,dst1))-sqrt(dot(dst2,dst2));
                            
                        }
                    }
                }
                
                for(int i=0;i<refStructAtom.size();i++){
                    vectorRef(i) = refStructAtom[i]->getAtomCharge();
                }
                
                //Solve linear system
                
                vectorObj = matrixObjObj.colPivHouseholderQr().solve(matrixObjRef*vectorRef);
                
                //update
                
                int i=0;
                for(proteinManager::MODEL& md : objStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
                                //Check
                                #ifdef DEBUG
                                std::cout << atm.getAtomCoord() << " " << atm.getAtomCharge() << " " << vectorObj[i] << std::endl;
                                #endif 
                                atm.setAtomCharge(vectorObj[i]);
								i++;
				}}}}
				
				this->updateAtomVectors();
                
            }
            
            void fitPosition(){
                
                gradientObj.resize(objStructAtom.size());
                for(int i=0;i<objStructAtom.size();i++){
                    gradientObj[i] = {0,0,0};
                }
                
                //Computing -gradient
                
                for(int i=0;i<objStructAtom.size();i++){
                    for(int j=0;j<objStructAtom.size();j++){
                        if(i != j){
                            real3 partGrdt = objStructAtom[i]->getAtomCoord()-objStructAtom[j]->getAtomCoord();
                            partGrdt = partGrdt/sqrt(dot(partGrdt,partGrdt)); 
                            partGrdt = partGrdt*real(4.0)*M_PI*objStructAtom[i]->getAtomCharge()*objStructAtom[j]->getAtomCharge();
                            gradientObj[i] += partGrdt;
                        }
                    }
                    for(int j=0;j<refStructAtom.size();j++){
                        real3 partGrdt = objStructAtom[i]->getAtomCoord()-refStructAtom[j]->getAtomCoord();
                        partGrdt = partGrdt/sqrt(dot(partGrdt,partGrdt)); 
                        partGrdt = partGrdt*real(-4.0)*M_PI*objStructAtom[i]->getAtomCharge()*refStructAtom[j]->getAtomCharge();
                        gradientObj[i] += partGrdt;
                    }
                }
                
                int i=0;
                for(proteinManager::MODEL& md : objStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
                                real3 newPos = objStructAtom[i]->getAtomCoord()+gradientObj[i]*lbda;
                                //Check
                                #ifdef DEBUG
                                std::cout << atm.getAtomCoord() << " " << newPos << std::endl;
                                #endif 
                                atm.setAtomCoord(newPos);
								i++;
				}}}}
                
                this->updateAtomVectors();
            }
            
            void fitChargesAndPositions(int nTries, bool cmpRes = true){
                
                this->fitCharges();
                
                real res_pre;
                real res_post;
                
                if(cmpRes){
                    res_pre = this->computeResidue();
                    res_post = 0;
                }
                
                for(int i=0;i<nTries;i++){
                    this->fitCharges();
                    this->fitPosition();
                    
                    if(cmpRes){
                        res_post = this->computeResidue();
                        //std::cout << "iter nun: " << i << " Res. before: " << res_pre << " Res. after: " << res_post << " Res. var.: " << res_post-res_pre << std::endl;
                        std::cout << i << " " << res_pre << " " << res_post << " " << res_post-res_pre << std::endl;
                        res_pre = res_post;
                    }
                }
            }
        
    };
}
}

#undef _USE_MATH_DEFINES
#endif
