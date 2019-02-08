#include <proteinManager/proteinManager.hpp>

namespace proteinManager{
    
    namespace geometricTransformations{
        
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
    }
    
    
}
