#ifndef __GEOMETRIC__TRANSFORMATIONS__
#define __GEOMETRIC__TRANSFORMATIONS__

#include <proteinManager/proteinManager.hpp>

#include <eigen3/Eigen/Eigen>
#include <random>

namespace proteinManager{
    
    namespace geometricTransformations{
        
        void translation(ATOM& atm,real3 translationVector){
            
            real3 oldCoord = atm.getAtomCoord();
            real3 newCoord;
            
            newCoord.x = oldCoord.x + translationVector.x;
            newCoord.y = oldCoord.y + translationVector.y;
            newCoord.z = oldCoord.z + translationVector.z;
            
            atm.setAtomCoord(newCoord);
        }
        
        void translation(RESIDUE& res,real3 translationVector){
            
            for(ATOM& atm : res.atom()){
                translation(atm,translationVector);
            }
        }
        
        void translation(CHAIN& ch,real3 translationVector){
            
            for(RESIDUE& res : ch.residue()){
                translation(res,translationVector);
            }
        }
        
        void translation(MODEL& mdl,real3 translationVector){
            
            for(CHAIN& ch : mdl.chain()){
                translation(ch,translationVector);
            }
        }
        
        void translation(STRUCTURE& str,real3 translationVector){
            
            for(MODEL& mdl : str.model()){
                translation(mdl,translationVector);
            }
        }
        
        ////////////////////////////////////////////////////////////////
        
        void homotheticTransformation(ATOM& atm,real3 center,real ratio){
            
            real3 oldCoord = atm.getAtomCoord();
            real3 newCoord;
            
            newCoord.x = center.x + ratio*(oldCoord.x - center.x);
            newCoord.y = center.y + ratio*(oldCoord.y - center.y);
            newCoord.z = center.z + ratio*(oldCoord.z - center.z);
            
            atm.setAtomCoord(newCoord);
        }
        
        void homotheticTransformation(RESIDUE& res,real3 center,real ratio){
            
            for(ATOM& atm : res.atom()){
                homotheticTransformation(atm,center,ratio);
            }
        }
        
        void homotheticTransformation(CHAIN& ch,real3 center,real ratio){
            
            for(RESIDUE& res : ch.residue()){
                homotheticTransformation(res,center,ratio);
            }
        }
        
        void homotheticTransformation(MODEL& mdl,real3 center,real ratio){
            
            for(CHAIN& ch : mdl.chain()){
                homotheticTransformation(ch,center,ratio);
            }
        }
        
        void homotheticTransformation(STRUCTURE& str,real3 center,real ratio){
            
            for(MODEL& mdl : str.model()){
                homotheticTransformation(mdl,center,ratio);
            }
        }
        
        ////////////////////////////////////////////////////////////////
        
        template <class T>
        void uniformScaling(T& input,real ratio){
            homotheticTransformation(input,{0,0,0},ratio);
        }
        
        ////////////////////////////////////////////////////////////////
        
        void rotation(ATOM& atm,const real3 center,const Eigen::Quaternion<real> quaternion){
            
            real3 oldCoord = atm.getAtomCoord();
            
            Eigen::Vector3f  p(oldCoord.x,oldCoord.y,oldCoord.z);
            Eigen::Vector3f  p_new;
                          
            Eigen::Vector3f  centerEigen(center.x,center.y,center.z);
            p_new = quaternion*(p-centerEigen)+centerEigen;
            
            atm.setAtomCoord({p_new.x(),p_new.y(),p_new.z()});
        }
        
        void rotation(RESIDUE& res,const real3 center,const Eigen::Quaternion<real> quaternion){
            
            for(ATOM& atm : res.atom()){
                rotation(atm,center,quaternion);
            }
        }
        
        void rotation(CHAIN& ch,const real3 center,const Eigen::Quaternion<real> quaternion){
            
            for(RESIDUE& res : ch.residue()){
                rotation(res,center,quaternion);
            }
        }
        
        void rotation(MODEL& mdl,const real3 center,const Eigen::Quaternion<real> quaternion){
            
            for(CHAIN& ch : mdl.chain()){
                rotation(ch,center,quaternion);
            }
        }
        
        void rotation(STRUCTURE& str,const real3 center,const Eigen::Quaternion<real> quaternion){
            
            for(MODEL& mdl : str.model()){
                rotation(mdl,center,quaternion);
            }
        }
        
        ////////////////////////////////////////////////////////////////
        
        void rotation(ATOM& atm,const real3 center,const real3 axis, const real angle_radians){
            
            Eigen::Vector3f  axisEigen(axis.x,axis.y,axis.z);
            
            Eigen::Quaternion<real> q;
            q = Eigen::AngleAxis<real>(angle_radians, axisEigen);
            
            rotation(atm,center,q);
        }
        
        void rotation(RESIDUE& res,const real3 center,const real3 axis, const real angle_radians){
            
            for(ATOM& atm : res.atom()){
                rotation(atm,center,axis,angle_radians);
            }
        }
        
        void rotation(CHAIN& ch,const real3 center,const real3 axis, const real angle_radians){
            
            for(RESIDUE& res : ch.residue()){
                rotation(res,center,axis,angle_radians);
            }
        }
        
        void rotation(MODEL& mdl,const real3 center,const real3 axis, const real angle_radians){
            
            for(CHAIN& ch : mdl.chain()){
                rotation(ch,center,axis,angle_radians);
            }
        }
        
        void rotation(STRUCTURE& str,const real3 center,const real3 axis, const real angle_radians){
            
            for(MODEL& mdl : str.model()){
                rotation(mdl,center,axis,angle_radians);
            }
        }
        
        ////////////////////////////////////////////////////////////////
        
        template<class entityType>
        void randomRotation(entityType& entity,const real3 center, std::mt19937& gen){
            
            std::uniform_real_distribution<real> distr(0.0,1.0);
            
            //Effective Sampling and Distance Metrics for 3D Rigid Body Path Planning
            //James J. Kuffner, 2004
            
            real s = distr(gen);
            
            real sigma1 = std::sqrt(real(1.0)-s);
            real sigma2 = std::sqrt(s);
            real theta1 = real(M_PI)*distr(gen);
            real theta2 = real(M_PI)*distr(gen);
            
            //x = sin(θ1)*σ1;
            //y = cos(θ1)*σ1;
            //z = sin(θ2)*σ2;
            //w = cos(θ2)*σ2;
            
            Eigen::Quaternion<real> randomQuaternion(std::sin(theta1)*sigma1,
                                                     std::cos(theta1)*sigma1,
                                                     std::sin(theta2)*sigma2,
                                                     std::cos(theta2)*sigma2);
            
            rotation(entity,center,randomQuaternion);
        }
        
        ////////////////////////////////////////////////////////////////
        
        template<class entityType>
        void randomRotation(entityType& entity,const real3 center,const real angle, std::mt19937& gen){
            
            std::uniform_real_distribution<real> distr(0.0,1.0);

            Eigen::Matrix<real,3,1> rndAxis;

            real theta = real(2.0)*M_PI*distr(gen);
            real phi   = acos(real(1.0)-real(2.0)*distr(gen));

            rndAxis(0) = sin(phi)*cos(theta);
            rndAxis(1) = sin(phi)*sin(theta);
            rndAxis(2) = cos(phi);
            
            Eigen::Quaternion<real> q(Eigen::AngleAxis<real>(angle,rndAxis));
            
            rotation(entity,center,q);
        }
    }
    
    
}

#endif
