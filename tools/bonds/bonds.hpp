#include <proteinManager/proteinManager.hpp>

#include "../geometric/geometric.hpp"

namespace proteinManager{
    
    namespace bonds{

        struct pair{
            int i,j;
            std::string type_i,type_j;

            real r0;
        };
        
        struct angle{
            int i,j,k;
            std::string type_i,type_j,type_k;

            real theta0;
        };
        
        struct dihedral{
            int i,j,k,l;
            std::string type_i,type_j,
                        type_k,type_l;

            real phi0;
        };

        template<typename pairType>
        void residuePairBonds(STRUCTURE& structureIn,std::vector<pairType>& bondVector){
            
            for(MODEL&   mdl : structureIn.model()       ){
            for(CHAIN&   ch  : mdl.chain()               ){
            for(uint i=0 ; i<ch.residue().size()-1 ; i++ ){
                
                ATOM atm1 = ch.residue()[i  ].atom()[0];
                ATOM atm2 = ch.residue()[i+1].atom()[0];
                
                pairType pairBuffer;
                pairBuffer.i = atm1.getAtomSerial();
                pairBuffer.j = atm2.getAtomSerial();
                pairBuffer.type_i = atm1.getAtomName();
                pairBuffer.type_j = atm2.getAtomName();
                pairBuffer.r0=geometric::dst(atm1,atm2);
                pairBuffer.r0=geometric::dst(atm1,atm2);

                bondVector.push_back(pairBuffer);

            }}}

        }
        
        template<typename angleType>
        void residueAngleBonds(STRUCTURE& structureIn,std::vector<angleType>& angleVector){
            
            for(MODEL&   mdl : structureIn.model()      ){
            for(CHAIN&   ch  : mdl.chain()              ){
            for(uint i=0 ; i<ch.residue().size()-2 ; i++){

                ATOM atm1 = ch.residue()[i  ].atom()[0];
                ATOM atm2 = ch.residue()[i+1].atom()[0];
                ATOM atm3 = ch.residue()[i+2].atom()[0];

                angleType angleBuffer;

                angleBuffer.i      = atm1.getAtomSerial();
                angleBuffer.j      = atm2.getAtomSerial();
                angleBuffer.k      = atm3.getAtomSerial();
                angleBuffer.type_i = atm1.getAtomName();
                angleBuffer.type_j = atm2.getAtomName();
                angleBuffer.type_k = atm3.getAtomName();
                angleBuffer.theta0 = geometric::computeAngle(atm1,atm2,atm3);

                angleVector.push_back(angleBuffer);
            }}}
        }
        
        template<typename dihedralType>
        void residueDihedralBonds(STRUCTURE& structureIn,std::vector<dihedralType>& dihedralVector){
            
            for(MODEL&   mdl : structureIn.model()      ){
            for(CHAIN&   ch  : mdl.chain()              ){
            for(uint i=0 ; i<ch.residue().size()-3 ; i++){
                
                ATOM atm1 = ch.residue()[i  ].atom()[0];
                ATOM atm2 = ch.residue()[i+1].atom()[0];
                ATOM atm3 = ch.residue()[i+2].atom()[0];
                ATOM atm4 = ch.residue()[i+3].atom()[0];
                
                dihedralType dihedralBuffer;

                dihedralBuffer.i    = atm1.getAtomSerial();
                dihedralBuffer.j    = atm2.getAtomSerial();
                dihedralBuffer.k    = atm3.getAtomSerial();
                dihedralBuffer.l    = atm4.getAtomSerial();
                dihedralBuffer.type_i = atm1.getAtomName();
                dihedralBuffer.type_j = atm2.getAtomName();
                dihedralBuffer.type_k = atm3.getAtomName();
                dihedralBuffer.type_l = atm4.getAtomName();
                dihedralBuffer.phi0 = geometric::computeDihedral(atm1,atm2,atm3,atm4);

                dihedralVector.push_back(dihedralBuffer);
            }}}
        }

        class KaranicolasBrooksDihedral{
        
            public:
        
                struct KBdihedralInfo{
                    real K1,K2,K3,K4;
                    real a1,a2,a3,a4;
                };
        
                struct KBdihedral : public bonds::dihedral, 
                                    public KBdihedralInfo{};
        
            private:
        
                std::string KBdiheParamFilePath = "../../data/KaranicolasBrooksDiheParm.dat";
                
                std::map<std::tuple<std::string,std::string>,std::vector<real>> KBdiheParamMap; 
        
            public:
        
                KaranicolasBrooksDihedral(){
                
                    std::stringstream ss;
                    std::string line;
                    
                    std::ifstream KBdiheParamFile(KBdiheParamFilePath);
                    if(!KBdiheParamFile){
                        ss << "File not found: " << KBdiheParamFilePath;
                        throw std::runtime_error(ss.str());
                    }
                
                    std::string atm1,atm2;
                    int   n;
                    real  V;
                    real  theta0;
        
                    while(std::getline(KBdiheParamFile,line)){
                        
                        //Empty lines or lines which starts with # are ignored
                        if(line[0]=='#' or line.empty()) continue;
                        
                        //Process line
                        ss.clear();
                        ss.str(line);
                        ss >> atm1 >> atm2 >> V >> n >> theta0;
        
                        std::tuple<std::string,std::string> key=std::make_tuple(atm1,atm2);
        
                        KBdiheParamMap[key].resize(8);
                        KBdiheParamMap[key][n-1]=V;
                        KBdiheParamMap[key][n-1+4]=theta0;
                    }
        
                    #ifdef DEBUG
                    for(auto& k : KBdiheParamMap){
                    
                        std::cout << std::get<0>(k.first) << " " << std::get<1>(k.first) << " " ;
                    
                        for(auto& v : k.second){
                            if(v==0){
                                ss.clear();
                                ss << "diheMap: no value added for " 
                                   << std::get<0>(k.first)      << " " 
                                   << std::get<1>(k.first);
                                throw std::runtime_error(ss.str());
                            }
                            std::cout << v << " ";
                        }
                    
                        std::cout << std::endl;
                    }
                    #endif
                }
        
            void feedDihedralList(std::vector<KBdihedral>& diheList){
                
                for(KBdihedral& kbd : diheList){
                    std::vector<real> diheValues = KBdiheParamMap[std::make_tuple(kbd.type_j,
                                                                                  kbd.type_k)];
                    kbd.K1 = diheValues[0];
                    kbd.K2 = diheValues[1];
                    kbd.K3 = diheValues[2];
                    kbd.K4 = diheValues[3];
                    
                    kbd.a1 = diheValues[4];
                    kbd.a2 = diheValues[5];
                    kbd.a3 = diheValues[6];
                    kbd.a4 = diheValues[7];
        
                }
                
            }
        
        };

        
    }
    
    
}
