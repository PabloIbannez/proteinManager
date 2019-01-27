#ifndef ENM_MODELS_HPP
#define ENM_MODELS_HPP

namespace proteinManager{
    
    namespace enm_models{
            
            struct caOrellana{
                
                real rcut;
                
                int M = 3;
                
                real Cseq = 60;
                real Ccart = 6;
                
                int n_seq = 2;
                int n_cart = 6;
                
                const char* info = "# Approaching Elastic Network Models to Molecular Dynamics Flexibility \n \
                                    # Laura Orellana, Manuel Rueda, Carles Ferrer-Costa, José Ramón Lopez-Blanco, Pablo Chacón, and Modesto Orozco \n \
                                    # Journal of Chemical Theory and Computation 2010 6 (9), 2910-2923 \n \
                                    # DOI: 10.1021/ct100208e \n \
                                    #\n \
                                    #\n \
                                    # K units: kcal/(mol·A^2)\n";
                
                bool init(STRUCTURE& str){ 
                    
                    int resNum = 0;
                    
                    for(MODEL& mdl : str.model()){
                    for(CHAIN& ch  : mdl.chain()){
                    for(RESIDUE& res : ch.residue()){
                        resNum++;
                    }}}
                    
                    rcut = 2.9*log(resNum)-2.9;
                    
                    return true;
                }
                
                bool init(MODEL&     str){ return true;}
                bool init(CHAIN&     str){ return true;}
                bool init(RESIDUE&   str){ return true;}
                
                bond computeBond(std::shared_ptr<ATOM> atm1, std::shared_ptr<ATOM> atm2){
                    
                    real3 r12 = atm1->getAtomCoord()-atm2->getAtomCoord();
                    real r = sqrt(dot(r12,r12));
                    
                    int S12 = abs(atm1->getParentResidue()->getResSeq()-atm2->getParentResidue()->getResSeq());
                    
                    if(S12 <= M){
    
                        return {atm1,atm2,r,Cseq/(pow(S12,n_seq))};
                        
                    } else {
                        
                        if(r <= rcut){
                            return {atm1,atm2,r,pow(Ccart/r,n_cart)};
                        } else {
                            return {atm1,atm2,0,0};
                        }
                    }
                }
            };
            
        }
}

#endif
