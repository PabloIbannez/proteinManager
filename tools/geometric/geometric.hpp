#ifndef __GEOMETRIC__
#define __GEOMETRIC__

#include <proteinManager/proteinManager.hpp>

namespace proteinManager{
    
    namespace geometric{

        namespace geometric_ns{
            
            real sign(real x){
                if(x>real(0)){return real( 1);}
                if(x<real(0)){return real(-1);}
                return real(0);
            }
        
        }

        real dst(ATOM& atmi,ATOM& atmj){
            real3 dr = atmj.getAtomCoord()-atmi.getAtomCoord();
            return sqrt(dot(dr,dr));
        }
        
        real dst(ATOM& atmi,real3 rj){
            real3 dr = rj-atmi.getAtomCoord();
            return sqrt(dot(dr,dr));
        }
        
        real dst(real3 ri,real3 rj){
            real3 dr = rj-ri;
            return sqrt(dot(dr,dr));
        }
        
        real computeAngle(ATOM& atmi,ATOM& atmj,ATOM& atmk){
            
            real3 ri=atmi.getAtomCoord();
            real3 rj=atmj.getAtomCoord();
            real3 rk=atmk.getAtomCoord();
        
            real3 dr21     = ri-rj;
            real  r21      = sqrt(dot(dr21,dr21));

            dr21/=r21;
            
            real3 dr32     = rk-rj;
            real  r32    = sqrt(dot(dr32,dr32));

            dr32/=r32;

            return acos(dot(dr21,dr32));
        }
        
        real computeDihedral(ATOM& atmi,ATOM& atmj,ATOM& atmk,ATOM& atml){
            
            real3 ri=atmi.getAtomCoord();
            real3 rj=atmj.getAtomCoord();
            real3 rk=atmk.getAtomCoord();
            real3 rl=atml.getAtomCoord();
        
            real3 rij=rj-ri;
            real3 rjk=rk-rj;
            real3 rkl=rl-rk;
            
            real3 va=cross(rij,rjk);
            real3 vb=cross(rjk,rkl);
            
            real vam = sqrt(dot(va,va));
            real vbm = sqrt(dot(vb,vb));
            
            va/=vam;
            vb/=vbm;
        
            real s = dot(cross(rij,rjk),rkl);
                 s = geometric_ns::sign(s);
        
            return s*acos(dot(va,vb));
        }
        
    }
    
    
}

#endif
