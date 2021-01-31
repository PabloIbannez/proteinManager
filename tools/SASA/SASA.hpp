#ifndef __PCASA__
#define __PCASA__

#include <proteinManager/proteinManager.hpp>

#include "../../tools/geometricTransformations/geometricTransformations.hpp"

#include "../../tools/geometric/geometric.hpp"
#include "../../tools/neighbourList/neighbourList.hpp"

namespace proteinManager{
    
    class SASA{
        
        private:

            const real rCut = 10.0;
            const int gamma = -1;

            using typePair = std::tuple<std::string,std::string>;

            std::string alphaFilePath   = "../../../SASA/data/pcasaAlpha.dat";
            std::string rndCoilFilePath = "../../../SASA/data/rndCoil.dat";
            
            std::string betaFilePath    = "../../../SASA/data/pcasaBeta.dat";
            
            std::map<std::string,real> alpha;
            std::map<std::string,real> rndCoil;
            
            std::map<typePair,real>    beta;
                
            void loadFromFile(std::string filePath,std::map<std::string,real>& m){
                
                std::string type;
                real        value;
                
                std::stringstream parser;
                std::string line;
                
                std::fstream file(filePath);
                while(std::getline(file,line)){

                    parser.clear();
                    parser.str(line);

                    parser >> type >> value;

                    if(m.count(type) == 0){
                        m[type] = value;
                    } else {
                        throw std::runtime_error("Type added previously "+type);
                    }

                }
            }
            
            void loadFromFile(std::string filePath,std::map<typePair,real>& m){
                
                std::string type1,type2;
                real        value;
                
                std::stringstream parser;
                std::string line;
                
                std::fstream file(filePath);
                while(std::getline(file,line)){
                    parser.clear();
                    parser.str(line);

                    parser >> type1 >> type2 >> value;

                    typePair pair = {type1,type2};
                    
                    if(m.count(pair) == 0){
                        m[pair] = value;
                    } else {
                        throw std::runtime_error("Type added previously " + type1 + " " + type2);
                    }

                }
                
            }
            
        public:

            SASA(){
                
                //Load alpha
                
                loadFromFile(alphaFilePath,alpha);
                loadFromFile(rndCoilFilePath,rndCoil);
                
                loadFromFile(betaFilePath,beta);
                
                /*
                for(auto& a : alpha)  {std::cout << a.first << " " << a.second << std::endl;}
                std::cout << std::endl;
                for(auto& r : rndCoil){std::cout << r.first << " " << r.second << std::endl;}
                std::cout << std::endl;
                */
                
                /*
                for(auto& b : beta){std::cout << std::get<0>(b.first) << " " 
                                              << std::get<1>(b.first) << " "                          
                                              << b.second << std::endl;}
                */

            }

            void addSASArndCoil(proteinManager::STRUCTURE& structure){
            
                for(MODEL&   mdl : structure.model()){
                for(CHAIN&   ch  : mdl.chain()      ){
                for(RESIDUE& res : ch.residue()     ){
                    
                    if(res.atom().size()>1){
                        throw std::runtime_error("SASA rndCoil, number of atoms > 1 " 
                                                  + res.getResName() + " " + std::to_string(res.getResSeq()));
                    }

                for(ATOM&    atm : res.atom()       ){
                    atm.setAtomSASA(rndCoil[atm.getAtomName()]);
                }}}}
            }
            
            void addPCASA(proteinManager::STRUCTURE& structure,bool merge_models){

                auto atomList = structure.atom();

                auto nL =  neighbourList::generateNeighbourList(atomList,rCut);

                //Complete neig
                for(int i; i<atomList.size(); i++){
                    for(int j : nL[i]){
                        if(i<j){
                            nL[j].push_back(i);
                        }
                }}

                for(int i; i<atomList.size(); i++){
                    std::string type_i = atomList[i].getAtomName();
                    real sum=0;
                    for(int j : nL[i]){
                        
                        if(!merge_models){
                            int mdl_i = atomList[i].getModelId();
                            int mdl_j = atomList[j].getModelId();

                            if(mdl_i != mdl_j){
                                continue;
                            }
                        }

                        std::string type_j = atomList[j].getAtomName();
                        real r = geometric::dst(atomList[i],atomList[j]);
                        sum += -beta[std::make_tuple(type_i,type_j)]*pow(r,gamma);
                        //std::cout << type_i << " " << type_j << " " << beta[std::make_tuple(type_i,type_j)] << " " << pow(r,gamma) << std::endl;
                    }
                    sum = alpha[type_i] + sum;
                    if(sum<0){sum=0;}
                    //std::cin.get();

                    int modelId = atomList[i].getModelId();
                    std::string ch = atomList[i].getChainId();
                    int res = atomList[i].getResSeq();
                    
                    if(structure.model(modelId).chain(ch).residue(res).atom().size()>1){
                        throw std::runtime_error("SASA PCASA, number of atoms > 1 ");
                    }

                    structure.model(modelId).chain(ch).residue(res).atom()[0].setAtomSASA(sum);
                }
            }

    };
    
    
}

#endif
