#ifndef __NATIVE_CONTACT_SOP__
#define __NATIVE_CONTACT_SOP__

#include <proteinManager/proteinManager.hpp>

#include "../../tools/geometric/geometric.hpp"
#include "../../tools/neighbourList/neighbourList.hpp"

namespace proteinManager{

class nativeContactSOP{

    public:

        struct nativeContact{
            int iMol,jMol;
            std::string chi,chj;
            int resi,resj;
            int iSerial,jSerial;
            real r0;
            real e;
        };

        struct Parameters{

            real cutOff;

            real eIntra;
            real eInter;

            bool merge_models;

        };
            
        std::vector<nativeContact> nc;

        nativeContactSOP(boost::ptr_vector<ATOM>& atomOutVector,Parameters par){

            real cutOff = par.cutOff;

            real eIntra = par.eIntra;
            real eInter = par.eInter;

            bool merge_models = par.merge_models;
            
            auto nL = neighbourList::generateNeighbourList(atomOutVector,cutOff);

            for(uint i = 0  ;i<atomOutVector.size();i++){
            std::cerr << "\r" << "Generating native contacts " << i+1 << "/" << atomOutVector.size()  << std::flush;
                for(uint j : nL[i]){

                    if(merge_models){

                    } else {
                        if(atomOutVector[i].getModelId()!=atomOutVector[j].getModelId()){continue;}
                    }

                    int iMol = atomOutVector[i].getModelId();
                    int jMol = atomOutVector[j].getModelId();

                    std::string chaini = atomOutVector[i].getChainId();
                    std::string chainj = atomOutVector[j].getChainId();
                    
                    int resi_seq = atomOutVector[i].getResSeq();
                    int resj_seq = atomOutVector[j].getResSeq();
                    
                    int iSerial = atomOutVector[i].getAtomSerial();
                    int jSerial = atomOutVector[j].getAtomSerial();

                    nativeContact nCbuffer;

                    nCbuffer.iMol = iMol;
                    nCbuffer.jMol = jMol;
                    
                    nCbuffer.chi = chaini;
                    nCbuffer.chj = chainj;
                    
                    nCbuffer.resi = resi_seq;
                    nCbuffer.resj = resj_seq;
                    
                    nCbuffer.iSerial = iSerial;
                    nCbuffer.jSerial = jSerial;

                    if(iSerial<jSerial){
                        if(iMol==jMol and chaini==chainj ){
                            if((abs(resi_seq-resj_seq)>2)){
                                nCbuffer.r0 = geometric::dst(atomOutVector[i],atomOutVector[j]);
                                nCbuffer.e = eIntra;
                                nc.push_back(nCbuffer);
                            }
                        } else {
                            nCbuffer.r0 = geometric::dst(atomOutVector[i],atomOutVector[j]);
                            nCbuffer.e = eInter;
                            nc.push_back(nCbuffer);
                        }

                    }
                }
            }
            
            std::cerr << std::endl;
        }

        std::vector<nativeContact> getNativeContactList(){
            return nc;
        }

};

}

#endif
