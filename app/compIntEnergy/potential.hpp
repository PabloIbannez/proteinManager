#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

namespace proteinManager{
namespace compInt{
        
    struct sasaPot{
        
        real cutOff_;
        real invEpsilon_;
        real invDebyeLength_;
        //real f = 332.0636 ; //kcal·Å/(mol·e2)
        real f = 33.20636 ; //kcal·nm/(mol·e2)
        
        sasaPot(real cutOff,real epsilon,real dL){
            
            cutOff_ = cutOff;
            invEpsilon_ = 1.0/epsilon;
            invDebyeLength_ = 1.0/dL;
            
        }
        
        real energy(ATOM& atm_1){
            
            real energy;
            
            energy = (atm_1.getAtomSASA())*atm_1.getAtomSolvE();
            
            return energy;
        }
        
        real energy(ATOM& atm_i,ATOM& atm_j){
            
            real3 pi = atm_i.getAtomCoord();
            real3 pj = atm_j.getAtomCoord();
            
            real3 rji = pj-pi;

            const real r      = sqrt(dot(rji, rji));
            
            if(r > cutOff_ or r == real(0.0)) return real(0);
            
            const real invr   = real(1)/r;
            const real invr2  = invr*invr;
            const real invr6  = invr2*invr2*invr2;
            const real invr12 = invr6*invr6;
            
            real energy = 0;
            
            //vdW
            
            real C6  = atm_i.getAtomC6()*atm_j.getAtomC6();
            real C12 = atm_i.getAtomC12()*atm_j.getAtomC12();
            
            energy += C12*invr12-C6*invr6;
            
            //ele
            
            real chgProduct = atm_i.getAtomCharge()*atm_j.getAtomCharge();
            
            energy += f*chgProduct*exp(-invDebyeLength_*r)*invEpsilon_*invr;
            
            return energy;
        }
        
    };

}}

#endif
