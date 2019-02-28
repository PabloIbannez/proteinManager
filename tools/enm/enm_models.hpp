#ifndef ENM_MODELS_HPP
#define ENM_MODELS_HPP

namespace proteinManager{
    
    namespace enm_models{
        
        struct basic_nm{
                
                real rcut = 1;
                
                real k = 2000;
                
                const char* info = "# K units: kJ/(mol·nm^2)\n";
                
                bool init(STRUCTURE& str){return true;}
                bool init(MODEL&     mdl){return true;}
                bool init(CHAIN&     ch ){return true;}
                bool init(RESIDUE&   res){return true;}
                
                bond computeBond(std::shared_ptr<ATOM> atm1, std::shared_ptr<ATOM> atm2){
                    
                    real3 r12 = atm1->getAtomCoord()-atm2->getAtomCoord();
                    real r = sqrt(dot(r12,r12));
                    
                    if(r < rcut){
                        return {atm1,atm2,r,k};
                    } else {
                        return {atm1,atm2,0,0};
                    }
                }
            };
            
            struct caOrellana_A{
                
                real rcut;
                
                int M = 3;
                
                real Cseq = 60;
                real Ccart = 6;
                
                int n_seq = 2;
                int n_cart = 6;
                
                const char* info = "# Approaching Elastic Network Models to Molecular Dynamics Flexibility \n"
                                   "# Laura Orellana, Manuel Rueda, Carles Ferrer-Costa, José Ramón Lopez-Blanco, Pablo Chacón, and Modesto Orozco \n"
                                   "# Journal of Chemical Theory and Computation 2010 6 (9), 2910-2923 \n"
                                   "# DOI: 10.1021/ct100208e \n"
                                   "#\n"
                                   "#\n"
                                   "# K units: kcal/(mol·A^2)\n";
                
                bool init(STRUCTURE& str){ 
                    
                    int resNum = 0;
                    
                    for(MODEL& mdl : str.model()){
                    for(CHAIN& ch  : mdl.chain()){
                    for(RESIDUE& res : ch.residue()){
                        resNum++;
                        if(res.atom().size() > 1){
                            std::stringstream ss;
                            ss << "The current ENM model expects a structure in which "
                                  "each residue is represented by one atom only.";
                            throw std::runtime_error(ss.str());
                        }
                    }}}
                    
                    rcut = 2.9*log(resNum)-2.9;
                    
                    return true;
                }
                
                bool init(MODEL&     mdl){ 
                    
                    int resNum = 0;
                    
                    for(CHAIN& ch  : mdl.chain()){
                    for(RESIDUE& res : ch.residue()){
                        resNum++;
                        if(res.atom().size() > 1){
                            std::stringstream ss;
                            ss << "The current ENM model expects a structure in which "
                                  "each residue is represented by one atom only.";
                            throw std::runtime_error(ss.str());
                        }
                    }}
                    
                    rcut = 2.9*log(resNum)-2.9;
                    
                    return true;
                }
                
                bool init(CHAIN&     ch){ 
                    
                    int resNum = 0;
                    
                    for(RESIDUE& res : ch.residue()){
                        resNum++;
                        if(res.atom().size() > 1){
                            std::stringstream ss;
                            ss << "The current ENM model expects a structure in which "
                                  "each residue is represented by one atom only.";
                            throw std::runtime_error(ss.str());
                        }
                    }
                    
                    rcut = 2.9*log(resNum)-2.9;
                    
                    return true;
                }
                
                bool init(RESIDUE&     res){ 
                    
                    std::stringstream ss;
                    ss << "The current ENM model cannot build the EN on a residue only.";
                    throw std::runtime_error(ss.str());
                    
                    return true;
                }
                
                bond computeBond(std::shared_ptr<ATOM> atm1, std::shared_ptr<ATOM> atm2){
                    
                    real3 r12 = atm1->getAtomCoord()-atm2->getAtomCoord();
                    real r = sqrt(dot(r12,r12));
                    
                    int S12 = abs(atm1->getParentResidue()->getResSeq()-atm2->getParentResidue()->getResSeq());
                    
                    std::string chAtm1 = atm1->getParentResidue()->getParentChain()->getChainId();
                    std::string chAtm2 = atm2->getParentResidue()->getParentChain()->getChainId();
                    
                    if(S12 <= M and chAtm1 == chAtm2){
    
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
            
            struct REACH_A{
                
                real rcut = 18;
                
                real k12 = 170;
                real k13 = 1.1;
                real k14 = 6.8;
                
                real aIntra = 1460;
                real bIntra = 0.9;
                
                real aInter = 936;
                real bInter = 0.75;
                
                const char* info = "# REACH Coarse-Grained Normal Mode Analysis of Protein Dimer Interaction Dynamics,\n"
                                   "# Kei Moritsugu, Vandana Kurkal-Siebert, Jeremy C. Smith,\n"
                                   "# Biophysical Journal,Volume 97, Issue 4,2009,Pages 1158-1167\n"
                                   "# DOI: 10.1016/j.bpj.2009.05.015.\n"
                                   "#\n"
                                   "#\n"
                                   "# K units: kcal/(mol·A^2)\n";
                
                bool init(STRUCTURE& str){ 
                    
                    for(MODEL& mdl : str.model()){
                    for(CHAIN& ch  : mdl.chain()){
                    for(RESIDUE& res : ch.residue()){
                        if(res.atom().size() > 1){
                            std::stringstream ss;
                            ss << "The current ENM model expects a structure in which "
                                  "each residue is represented by one atom only.";
                            throw std::runtime_error(ss.str());
                        }
                    }}}
                    
                    return true;
                }
                
                bool init(MODEL&     mdl){ 
                    
                    for(CHAIN& ch  : mdl.chain()){
                    for(RESIDUE& res : ch.residue()){
                        if(res.atom().size() > 1){
                            std::stringstream ss;
                            ss << "The current ENM model expects a structure in which "
                                  "each residue is represented by one atom only.";
                            throw std::runtime_error(ss.str());
                        }
                    }}
                    
                    return true;
                }
                
                bool init(CHAIN&     ch){ 
                    
                    for(RESIDUE& res : ch.residue()){
                        if(res.atom().size() > 1){
                            std::stringstream ss;
                            ss << "The current ENM model expects a structure in which "
                                  "each residue is represented by one atom only.";
                            throw std::runtime_error(ss.str());
                        }
                    }
                    
                    return true;
                }
                
                bool init(RESIDUE&     res){ 
                    
                    std::stringstream ss;
                    ss << "The current ENM model cannot build the EN on a residue only.";
                    throw std::runtime_error(ss.str());
                    
                    return true;
                }
                
                bond computeBond(std::shared_ptr<ATOM> atm1, std::shared_ptr<ATOM> atm2){
                    
                    real3 r12 = atm1->getAtomCoord()-atm2->getAtomCoord();
                    real r = sqrt(dot(r12,r12));
                    
                    int S12 = abs(atm1->getParentResidue()->getResSeq()-atm2->getParentResidue()->getResSeq());
                    
                    std::string chAtm1 = atm1->getParentResidue()->getParentChain()->getChainId();
                    std::string chAtm2 = atm2->getParentResidue()->getParentChain()->getChainId();
                    
                    if(chAtm1 == chAtm2){ //intra
                        
                        if(S12 == 1){
                            return {atm1,atm2,r,k12};
                        } else if (S12 == 2){
                            return {atm1,atm2,r,k13};
                        } else if (S12 == 3){
                            return {atm1,atm2,r,k14};
                        } else {
                            
                            if(r < rcut){
                                return {atm1,atm2,r,aIntra*exp(-bIntra*r)};
                            }
                            
                            return {atm1,atm2,0,0};
                        }
                    
                    } else { //inter
                        
                        if(r < rcut){
                                return {atm1,atm2,r,aInter*exp(-bInter*r)};
                        }
                        
                        return {atm1,atm2,0,0};
                        
                    }
                }
            };
            
            struct REACH_nm{
                
                real rcut = 1.8;
                
                real k12 = 17000;
                real k13 = 100.1;
                real k14 = 600.8;
                
                real aIntra = 146000;
                real bIntra = 9;
                
                real aInter = 93600;
                real bInter = 7.5;
                
                const char* info = "# REACH Coarse-Grained Normal Mode Analysis of Protein Dimer Interaction Dynamics,\n"
                                   "# Kei Moritsugu, Vandana Kurkal-Siebert, Jeremy C. Smith,\n"
                                   "# Biophysical Journal,Volume 97, Issue 4,2009,Pages 1158-1167\n"
                                   "# DOI: 10.1016/j.bpj.2009.05.015.\n"
                                   "#\n"
                                   "#\n"
                                   "# K units: kcal/(mol·nm^2)\n";
                
                bool init(STRUCTURE& str){ 
                    
                    for(MODEL& mdl : str.model()){
                    for(CHAIN& ch  : mdl.chain()){
                    for(RESIDUE& res : ch.residue()){
                        if(res.atom().size() > 1){
                            std::stringstream ss;
                            ss << "The current ENM model expects a structure in which "
                                  "each residue is represented by one atom only.";
                            throw std::runtime_error(ss.str());
                        }
                    }}}
                    
                    return true;
                }
                
                bool init(MODEL&     mdl){ 
                    
                    for(CHAIN& ch  : mdl.chain()){
                    for(RESIDUE& res : ch.residue()){
                        if(res.atom().size() > 1){
                            std::stringstream ss;
                            ss << "The current ENM model expects a structure in which "
                                  "each residue is represented by one atom only.";
                            throw std::runtime_error(ss.str());
                        }
                    }}
                    
                    return true;
                }
                
                bool init(CHAIN&     ch){ 
                    
                    for(RESIDUE& res : ch.residue()){
                        if(res.atom().size() > 1){
                            std::stringstream ss;
                            ss << "The current ENM model expects a structure in which "
                                  "each residue is represented by one atom only.";
                            throw std::runtime_error(ss.str());
                        }
                    }
                    
                    return true;
                }
                
                bool init(RESIDUE&     res){ 
                    
                    std::stringstream ss;
                    ss << "The current ENM model cannot build the EN on a residue only.";
                    throw std::runtime_error(ss.str());
                    
                    return true;
                }
                
                bond computeBond(std::shared_ptr<ATOM> atm1, std::shared_ptr<ATOM> atm2){
                    
                    real3 r12 = atm1->getAtomCoord()-atm2->getAtomCoord();
                    real r = sqrt(dot(r12,r12));
                    
                    int S12 = abs(atm1->getParentResidue()->getResSeq()-atm2->getParentResidue()->getResSeq());
                    
                    std::string chAtm1 = atm1->getParentResidue()->getParentChain()->getChainId();
                    std::string chAtm2 = atm2->getParentResidue()->getParentChain()->getChainId();
                    
                    if(chAtm1 == chAtm2){ //intra
                        
                        if(S12 == 1){
                            return {atm1,atm2,r,k12};
                        } else if (S12 == 2){
                            return {atm1,atm2,r,k13};
                        } else if (S12 == 3){
                            return {atm1,atm2,r,k14};
                        } else {
                            
                            if(r < rcut){
                                return {atm1,atm2,r,aIntra*exp(-bIntra*r)};
                            }
                            
                            return {atm1,atm2,0,0};
                        }
                    
                    } else { //inter
                        
                        if(r < rcut){
                                return {atm1,atm2,r,aInter*exp(-bInter*r)};
                        }
                        
                        return {atm1,atm2,0,0};
                        
                    }
                }
            };
            
            struct go_dst_diffMol_nm{
                
                real rcut = 2.5;
                
                real e = 40;
                
                const char* info = "# e units: kJ, r units: nm\n";
                
                bool init(STRUCTURE& str){return true;}
                bool init(MODEL&     mdl){return true;}
                bool init(CHAIN&     ch ){return true;}
                bool init(RESIDUE&   res){return true;}
                
                bond computeBond(std::shared_ptr<ATOM> atm1, std::shared_ptr<ATOM> atm2){
                    
                    int mol1 = atm1->getModelId();
                    int mol2 = atm2->getModelId();
                    
                    if(mol1 == mol2){
                        return {atm1,atm2,0,0};
                    }
                    
                    real3 r12 = atm1->getAtomCoord()-atm2->getAtomCoord();
                    real r = sqrt(dot(r12,r12));
                    
                    if(r < rcut){
                        return {atm1,atm2,r,e};
                    } else {
                        return {atm1,atm2,0,0};
                    }
                }
            };
            
        }
}

#endif
