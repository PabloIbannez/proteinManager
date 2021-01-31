#ifndef __NEIGHBOUR_LIST__
#define __NEIGHBOUR_LIST__

#include <proteinManager/proteinManager.hpp>

#include "../geometric/geometric.hpp"

namespace proteinManager{
    
    namespace neighbourList{
        
        template<class atomVector>
        std::vector<std::vector<int>> generateNeighbourList(atomVector& atmV,real rCut){
        
            std::vector<std::vector<int>> neigList(atmV.size());
            
            //std::cerr << "Generating neighbour list (parallel) ..." << std::endl;
            #pragma omp parallel for schedule(dynamic)
            for(uint i = 0  ;i<atmV.size();i++){
            if(omp_get_thread_num()==0){
                std::cerr << "\r" << "Generating neighbour list (parallel) " 
                          << i+1 << "/" << atmV.size()  
                          << " (number of threads: " << omp_get_num_threads() << ")"<<std::flush;
            }
            for(uint j = i+1;j<atmV.size();j++){
                if(geometric::dst(atmV[i],atmV[j])<rCut){
                    neigList[i].push_back(j);
                }
            }}
            std::cerr << "\r\e[K" << "Generating neighbour list " 
                      << atmV.size() << "/" << atmV.size() << std::endl;
        
            return neigList;
        }


    }
}

#endif
