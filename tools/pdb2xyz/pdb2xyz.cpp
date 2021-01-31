#include <proteinManager/proteinManager.hpp>

#include "../../tools/geometricTransformations/geometricTransformations.hpp"
#include "../../tools/centers/centroid.hpp"

using namespace proteinManager;


real dst(const real3& r1,const real3& r2){
    real3 dr=r2-r1;
    return sqrt(dot(dr,dr));
}
int main(int argc, char *argv[]){
    
    STRUCTURE pdbIn;
    
    pdbIn.loadPDB(std::string(argv[1])+".pdb");
    
    real3 centroid = computeCentroid(pdbIn);
    geometricTransformations::translation(pdbIn,real(-1.0)*centroid);

    std::ofstream out(std::string(argv[1])+".xyz");

    auto& atmVector = pdbIn.atom();

    real3 pentCentroid = {0,0,0};
    uint  pentN=0;

    uint pentN_one = 0;
    std::vector<real3> pentCentroidAll;

    for(ATOM&   a : atmVector){
            
        if(a.getChainId() == "N" ){
            pentN_one++;
            pentCentroid += a.getAtomCoord();
        } else if(pentN_one!=0){
            pentCentroidAll.push_back(pentCentroid/pentN_one);
            pentCentroid = {0,0,0};
            pentN_one=0;
        }

        /*
        if(a.getModelId() == 0 or
                a.getModelId() == 1 or
                a.getModelId() == 2 or
                a.getModelId() == 3 or
                a.getModelId() == 4 ){
            if(a.getChainId() == "N" ){
                pentN++;
                pentCentroid += a.getAtomCoord();
            }
        }*/
    }

    for(int i=0;i<pentCentroidAll.size();i++){
    for(int j=i+1;j<pentCentroidAll.size();j++){
        std::cout << dst(pentCentroidAll[i],pentCentroidAll[j]) << std::endl;
    }}

    pentCentroid = pentCentroid/pentN;
    pentCentroid = pentCentroid/sqrt(dot(pentCentroid,pentCentroid));

    std::cout << pentCentroid << std::endl;

    for(ATOM&   a : atmVector){
        a.setAtomRadius(dst(a.getAtomCoord(),{0,0,0}));
    }
    
    //out << atmVector.size() << std::endl;
    out << 5486700 << std::endl;
    out << "*" << std::endl;
    for(ATOM&   a : atmVector){
        
        real3 coord = a.getAtomCoord();

        if(a.getModelId() == 0 or
           a.getModelId() == 1 or
           a.getModelId() == 2 or
           a.getModelId() == 3 or
           a.getModelId() == 4 ){
            if(a.getChainId() == "N" ){
                coord += real3({pentCentroid.x*250,pentCentroid.y*250,pentCentroid.z*250});
            }
        }
        
        if(a.getChainId() == "N" ){
            out << coord << " " << 1.671 << " " << a.getAtomRadius() << std::endl;
        }
        
        if(a.getChainId() == "A" or
           a.getChainId() == "B" or
           a.getChainId() == "C" ){
            out << coord << " " << 0.608 << " " << a.getAtomRadius() <<std::endl;
        }
        
        if(a.getChainId() == "D" or
           a.getChainId() == "E" or
           a.getChainId() == "F" ){
            out << coord << " " << 0.047 << " " << a.getAtomRadius() <<std::endl;
        }
        
        if(a.getChainId() == "G" or
           a.getChainId() == "H" or
           a.getChainId() == "I" ){
            out << coord << " " << -0.637 <<" " << a.getAtomRadius() << std::endl;
        }
        
        if(a.getChainId() == "J" or
           a.getChainId() == "K" or
           a.getChainId() == "L" ){
            out << coord << " " << -0.378 << " " << a.getAtomRadius() <<std::endl;
        }
    }
}
