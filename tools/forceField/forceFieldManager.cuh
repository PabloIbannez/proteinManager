#define DEBUG

#include <iostream> 
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <tuple>

#include <proteinManager/proteinManager.hpp>
#include "forceFieldManager.hpp"
#include "fieldComputing.cuh"
#include <integration3D/src/pIntegrator.cuh>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <cub/cub.cuh>

#include <Eigen/Eigen>

namespace proteinManager {
namespace ffManager{
    
    template <class potential>
    class forceFieldFitter{
	
		private:
		
        potential pot_;
        
		proteinManager::STRUCTURE& refStruct_;
		proteinManager::STRUCTURE& objStruct_;
	
		std::vector<std::shared_ptr<proteinManager::ATOM>> refStructAtom;
		std::vector<std::shared_ptr<proteinManager::ATOM>> objStructAtom;
		
		//integral matrices
		
		Eigen::VectorXd vector1;
		Eigen::VectorXd vector2_fit;
        
        Eigen::VectorXd vector2PF;
		Eigen::MatrixXd matrix21;
		Eigen::MatrixXd matrix22;
        
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
                             potential pot):refStruct_(refStruct),objStruct_(objStruct),pot_(pot){

                
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
                             Grid_PF_Integration intGrid_PF,
                             potential pot):refStruct_(refStruct),objStruct_(objStruct),pot_(pot){

                
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
                             potential pot):refStruct_(refStruct),objStruct_(objStruct),pot_(pot){

                
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
			
				real cutOff = pot_.getCutOff();
				
                //Determining box
                
				real3 pos = refStruct_.model()[0].chain()[0].residue()[0].atom()[0].getAtomCoord()/10;
				
				real3 boxMin = pos - cutOff;
				real3 boxMax = pos + cutOff;
				
				for(proteinManager::MODEL& md : refStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
					pos = atm.getAtomCoord()/10;
					
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
                
				real cutOff2 = pot_.getCutOff()*pot_.getCutOff();
				
				matrix21.resize(objStructAtom.size(),refStructAtom.size());
				matrix22.resize(objStructAtom.size(),objStructAtom.size());
				
                
				for(int i=0;i<objStructAtom.size();i++){
					std::cout << i << std::endl;
					for(int j=0;j<refStructAtom.size();j++){
						real3 atm1 = objStructAtom[i]->getAtomCoord()/10;
						real3 atm2 = refStructAtom[j]->getAtomCoord()/10;
						
						real3 dst = atm1 - atm2;
						
						if(dot(dst,dst) < cutOff2 and potential::getInteractionParm(*refStructAtom[j]) != real(0)) {
                            
                            auto potP21 = pot_.getPotentialProduct21({atm1.x,atm1.y,atm1.z},
                                                                     {atm2.x,atm2.y,atm2.z});
                            
                            if(currentIntegrator == grid ){
                                matrix21(i,j) = integGrid.computeIntegral(potP21);
                            } else if (currentIntegrator == MC){
                                matrix21(i,j) = integMC.computeIntegralAverage(potP21,samplesPerIntegral_).x;
                            } else {
                                ss.clear();
                                ss << "Selected integrator is not valid";
                                throw std::runtime_error(ss.str());
                            }
                            
						} else {
							matrix21(i,j) = 0;
						}
					}
				}
				
				
				for(int i=0;i<objStructAtom.size();i++){
					std::cout << i << std::endl;
					for(int j=i;j<objStructAtom.size();j++){
						real3 atm1 = objStructAtom[i]->getAtomCoord()/10;
						real3 atm2 = objStructAtom[j]->getAtomCoord()/10;
						
						real3 dst = atm1 - atm2;
						
						if(dot(dst,dst) < cutOff2) {
                            
                            auto potP22 = pot_.getPotentialProduct22({atm1.x,atm1.y,atm1.z},
                                                                     {atm2.x,atm2.y,atm2.z});
                            
                            if(currentIntegrator == grid ){
                                matrix22(i,j) = integGrid.computeIntegral(potP22);
                            } else if (currentIntegrator == MC){
                                matrix22(i,j) = integMC.computeIntegralAverage(potP22,samplesPerIntegral_).x;
                            } else {
                                ss.clear();
                                ss << "Selected integrator is not valid";
                                throw std::runtime_error(ss.str());
                            }

							matrix22(j,i) = matrix22(i,j);
                            
						} else {
							matrix22(i,j) = 0;
							matrix22(j,i) = 0;
						}
					}
				}
				
			}
            
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
                    
                    real3 atm1 = objStructAtom[i]->getAtomCoord()/10;
                    
                    auto potP2PF = pot_.getPotentialProduct2PF({atm1.x,atm1.y,atm1.z});
                    
                    vector2PF[i] = integGrid.computeIntegral_PF(potP2PF);
                    
                    
                }
                
                for(int i=0;i<objStructAtom.size();i++){
					std::cout << i << std::endl;
					for(int j=i;j<objStructAtom.size();j++){
						real3 atm1 = objStructAtom[i]->getAtomCoord()/10;
						real3 atm2 = objStructAtom[j]->getAtomCoord()/10;
						
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
			
			void computeNewParameters(){
                
                this->computeMatrices();
				
				vector1.resize(refStructAtom.size());
				vector2_fit.resize(objStructAtom.size());
				
				for(int i=0;i<refStructAtom.size();i++){
                    vector1[i] = potential::getInteractionParm(*refStructAtom[i]);
				}
				
                vector2_fit = matrix22.colPivHouseholderQr().solve(matrix21*vector1);
                
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
			
			void computeNewParametersTotalChargeConstraint(){
                
                this->computeMatrices();
				
				matrix21.conservativeResize(matrix21.rows()+1,matrix21.cols());
				matrix22.conservativeResize(matrix22.rows()+1,matrix22.cols()+1);
				
				vector1.resize(refStructAtom.size());
				vector2_fit.resize(objStructAtom.size()+1);
				
				///////////////////////////////////////////////////////////////////
				
				for(int i=0;i<refStructAtom.size();i++){
                    vector1[i] = potential::getInteractionParm(*refStructAtom[i]);
				}
				
				for(int i=0; i < matrix21.cols(); i++){
					matrix21(matrix21.rows()-1,i) = 1;
				}
				
				//Matrix22 is always a symetric matrix
				for(int i=0; i < matrix22.cols(); i++){
					matrix22(matrix22.rows()-1,i) = 1;
					matrix22(i,matrix22.cols()-1) = 1;
				}
				matrix22(matrix22.rows()-1,matrix22.cols()-1) = 0;
				
				///////////////////////////////////////////////////////////////////
				
				vector2_fit = matrix22.colPivHouseholderQr().solve(matrix21*vector1);
                
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
				
				#ifdef DEBUG
				
				real totalCharge1 = 0;
				real totalCharge2 = 0;
				
				for(int i=0;i<vector1.size();i++){
					totalCharge1 += vector1[i];
				}
				
				for(int i=0;i<vector2_fit.size()-1;i++){
					totalCharge2 += vector2_fit[i];
				}
				
				std::cout << "Ref charge: " << totalCharge1 << " Obj charge: " << totalCharge2 << std::endl;
				
				real3 totalDipole1 = {0,0,0};
				real3 totalDipole2 = {0,0,0};
				
				for(proteinManager::MODEL& md : refStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
								totalDipole1 += atm.getAtomCoord()*potential::getInteractionParm(atm);
				}}}}
				
				for(proteinManager::MODEL& md : objStruct_.model()) {
				for(proteinManager::CHAIN& ch : md.chain())         {
				for(proteinManager::RESIDUE& res : ch.residue())    {
				for(proteinManager::ATOM& atm : res.atom())         {
								totalDipole2 += atm.getAtomCoord()*potential::getInteractionParm(atm);
				}}}}
				
				std::cout << "Ref dipole: " << totalDipole1 << " Obj dipole: " << totalDipole2 << std::endl;
				
				#endif
				
				this->updateAtomVectors();
			}
		
    };
}
}

template <int exponent>
struct potentialInvExponent{
    
    real cutOff_;
    real e_;
    
    potentialInvExponent(real cutOff, real maxValue){
        cutOff_ = cutOff;
        e_ = pow(1.0/maxValue,1.0/exponent);
    }
    
    real getCutOff(){return cutOff_;}
    
    static real getInteractionParm(proteinManager::ATOM& atm){
        return atm.getAtomC6();
    }
        
    static void setInteractionParm(proteinManager::ATOM& atm,real parm){
        return atm.setAtomC6(parm);
    }
    
    struct potentialAtPoint{
        
        real cutOff2_;
    
        real3 pos_;
        real  par_;
        
        potentialAtPoint(real cutOff,real3 pos,real par){
            cutOff2_ = cutOff*cutOff;
            pos_ = pos;
            par_  = par;
        }
        
        __host__ __device__ real operator()(real3 p){
        
            real r2 = dot(p-pos_,p-pos_);
            
            if(r2 > cutOff2_) return real(0);
            
            if (exponent == 1){
                return par_*(1.0/sqrt(r2));
            } else {
                return par_*powf(real(1)/r2,exponent/int(2));
            }
        }
    };
    
    //In this case are identical
    potentialAtPoint getPotentialAtPoint(real3 pos,real par){
        return potentialAtPoint(cutOff_,pos,par);
    }
    
    potentialAtPoint getPotentialAtPointRef(real3 pos,real par){
        return potentialAtPoint(cutOff_,pos,par);
    }
    
    potentialAtPoint getPotentialAtPointObj(real3 pos,real par){
        return potentialAtPoint(cutOff_,pos,par);
    }
    
    struct potentialProduct21{
        
            real cutOff2_;
            real e2_;
            real3 p1_;
            real3 p2_;
        
            potentialProduct21(real cutOff, real e,real3 p1, real3 p2){
                cutOff2_ = cutOff*cutOff;
                e2_ = e*e;
                p1_ = p1;
                p2_ = p2;
            }
        
            __host__ __device__ real operator()(real3 p){
        
                real r1_2 = dot(p-p1_,p-p1_);
                real r2_2 = dot(p-p2_,p-p2_);
                
                if(r1_2 < e2_ or r2_2 < e2_) return real(0);
                if(r1_2 > cutOff2_ or r2_2 > cutOff2_) return real(0);
                
                if (exponent == 1){
                    return (1.0/sqrt(r1_2))*(1.0/sqrt(r2_2));
                } else {
                    return powf(real(1)/r1_2,exponent/int(2))*powf(real(1)/r2_2,exponent/int(2));
                }
            }
    };
    
    potentialProduct21 getPotentialProduct21(real3 p1, real3 p2){
        return potentialProduct21(cutOff_,e_,p1,p2);
    }
    
    potentialProduct21 getPotentialProduct22(real3 p1, real3 p2){
        return potentialProduct21(cutOff_,e_,p1,p2);
    }
    
};

struct potentialCoulombDebye{
    
    enum potSel {ref,obj};
    enum potSel pS_ = ref;
    
    real cutOff_;
    real e_s_;
    real k_;
    
    real f_ = 33.238147;
    
    potentialCoulombDebye(real cutOff,real e_s,real k){
        cutOff_ = cutOff;
        e_s_ = e_s;
        k_ = k;
    }
    
    real getCutOff(){return cutOff_;}
    
    static real getInteractionParm(proteinManager::ATOM& atm){
        return atm.getAtomCharge();
    }
        
    static void setInteractionParm(proteinManager::ATOM& atm,real parm){
        return atm.setAtomCharge(parm);
    }
    
    void selectPotential(enum potSel pS){pS_ = pS;}
    
    //Coulomb
    struct potentialAtPointRef{
        
        real f_;
        
        real cutOff2_;
    
        real3 pos_;
        real  par_;
        
        potentialAtPointRef(real f,real cutOff,real3 pos,real par){
            f_ = f;
            cutOff2_ = cutOff*cutOff;
            pos_ = pos;
            par_  = par;
        }
        
        __host__ __device__ real operator()(real3 p){
        
            real r2 = dot(p-pos_,p-pos_);
            
            if(r2 > cutOff2_) return real(0);
            
            return f_*par_/sqrt(r2);
        }
    };
    
    //Debye
    struct potentialAtPointObj{
        
        real f_;
        
        real cutOff2_;
        
        real e_s_;
        real k_;
    
        real3 pos_;
        real  par_;
        
        potentialAtPointObj(real f,real cutOff,real e_s,real k,real3 pos,real par){
            f_ = f;
            e_s_ = e_s;
            k_ = k;
            cutOff2_ = cutOff*cutOff;
            pos_ = pos;
            par_  = par;
        }
        
        __host__ __device__ real operator()(real3 p){
        
            real r2 = dot(p-pos_,p-pos_);
            
            if(r2 > cutOff2_) return real(0);
            
            real r = sqrt(r2);
            
            return f_*(par_*exp(-r/k_))/(e_s_*r);
        }
    };
    
    //In this case are identical
    potentialAtPointRef getPotentialAtPointRef(real3 pos,real par){
        return potentialAtPointRef(f_,cutOff_,pos,par);
    }
    
    potentialAtPointObj getPotentialAtPointObj(real3 pos,real par){
        return potentialAtPointObj(f_,cutOff_,e_s_,k_,pos,par);
    }
    
    struct potentialProduct22{
            
            real f2_;
            
            real cutOff2_;
            real e_s_;
            real k_;
            real3 p1_;
            real3 p2_;
        
            potentialProduct22(real f,real cutOff,real e_s,real k,real3 p1, real3 p2){
                f2_ = f*f;
                cutOff2_ = cutOff*cutOff;
                e_s_ = e_s;
                k_ = k;
                p1_ = p1;
                p2_ = p2;
            }
        
            __host__ __device__ real operator()(real3 p){
        
                real r1_2 = dot(p-p1_,p-p1_);
                real r2_2 = dot(p-p2_,p-p2_);
                
                if(r1_2 > cutOff2_ or r2_2 > cutOff2_) return real(0);
                //if(r1_2 < 1 or r2_2 < 1) return real(0);
                
                real r1 = sqrt(r1_2);
                real r2 = sqrt(r2_2);
                
                return f2_*(exp(-r1/k_)/(e_s_*r1))*(exp(-r2/k_)/(e_s_*r2));
            }
    
    };
    
    potentialProduct22 getPotentialProduct22(real3 p1, real3 p2){
        return potentialProduct22(f_,cutOff_,e_s_,k_,p1,p2);
    }
    
    struct potentialProduct2PF{
            
            real f_;
            
            real cutOff2_;
            real e_s_;
            real k_;
            real3 part_;
        
            potentialProduct2PF(real f,real cutOff,real e_s,real k,real3 part){
                f_ = f;
                cutOff2_ = cutOff*cutOff;
                e_s_ = e_s;
                k_ = k;
                part_ = part;
            }
        
            __host__ __device__ real operator()(real3 p,real field){
        
                real r_2 = dot(p-part_,p-part_);
                
                if(r_2 > cutOff2_) return real(0);
                //if(r_2 < 1 ) return real(0);
                
                real r = sqrt(r_2);
                
                return f_*exp(-r/k_)/(e_s_*r)*field;
            }
    
    };
    
    potentialProduct2PF getPotentialProduct2PF(real3 part){
        return potentialProduct2PF(f_,cutOff_,e_s_,k_,part);
    }
    
};

int main(){
    
    
    proteinManager::STRUCTURE pdbRef;
    proteinManager::STRUCTURE pdbObj;
    
    pdbRef.loadPDB("./examples/1aki.pqr");
    pdbObj.loadPDB("./examples/1akiChgPosFit.pqr");
    
    //pdbRef.loadPDB("./examples/3apg.pqr");
    //pdbObj.loadPDB("./examples/3apgCG_MARRINK.pqr");
    
    /////////////////////////////////////////////////////////////////
    
    proteinManager::ffManager::forceFieldManager ffM;
    
    ffM.loadForceFieldData("./forceFieldModels/gromos.ff");
    ffM.applyForceFieldData(pdbRef);
    
    /////////////////////////////////////////////////////////////////
    
    using potential = potentialInvExponent<6>;
    potential  pot(1,200); //c12 1,20
    
    
    proteinManager::ffManager::forceFieldFitter<potential>::MonteCarloIntegration mc;
    mc.pointsPerIntegral = 100000;
    mc.samplesPerIntegral = 10;
    
    proteinManager::ffManager::forceFieldFitter<potential> ffF(pdbRef,pdbObj,mc,pot);
    
    
    //proteinManager::ffManager::forceFieldFitter<potential>::GridIntegration gInt;
    //gInt.cellSize = 0.05;
    //
    //proteinManager::ffManager::forceFieldFitter<potential> ffF(pdbRef,pdbObj,gInt,pot);
    /*
    using potential = potentialCoulombDebye;
    potential pot(3,80,0.97);
    
    proteinManager::ffManager::forceFieldFitter<potential>::Grid_PF_Integration gridPF;
    gridPF.inputFilePath = "./examples/phimap1aki.cube";
    gridPF.lFactor = 0.1;
    gridPF.fFactor = 0.593;
    
    proteinManager::ffManager::forceFieldFitter<potential> ffF(pdbRef,pdbObj,gridPF,pot);
    */
    
	ffF.computeNewParametersTotalChargeConstraint();
    //ffF.computeNewParametersUsingPF();
	
	/////////////////////////////////////////////////////////////////
    
    proteinManager::fieldComputing::fieldComputing fC;
    
    std::ofstream ref("refField.dat");
    std::ofstream obj("objField.dat");
    
    fC.init({-0.0356,-0.6309,-2.6852},{5.7644,5.6191,2.7898},0.05);
    //fC.init({-0.4225,-11.191,-6.5898},{10.6025,2.884,5.4602},0.05);
    
    fC.computeRefField<potential>(pdbRef,pot);
    //fC.outputIndex_Field(ref);
    fC.output_CUBE(ref," "," ",10,1);
    //fC.output_CUBE(ref," "," ",10,0.1);
    
    fC.setFieldValue(0);
    
    fC.computeObjField<potential>(pdbObj,pot);
    //fC.outputIndex_Field(obj);
    fC.output_CUBE(obj," "," ",10,1);
    //fC.output_CUBE(obj," "," ",10,1/0.593);
    
    return EXIT_SUCCESS;
}
