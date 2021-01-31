#include "../../tools/geometric/geometric.hpp"
#include "../../tools/neighbourList/neighbourList.hpp"

namespace proteinManager{

    int checkClash(int modelId_i,int modelId_j,
                   STRUCTURE& linkersCargos, 
                   real cargoRadius,real atmRadius){
        
        real3 cargo_i_pos = computeCentroid(linkersCargos.model(modelId_i).chain("C"));
        real3 cargo_j_pos = computeCentroid(linkersCargos.model(modelId_j).chain("C"));
    
        auto atmLinker_i = linkersCargos.model(modelId_i).chain("L").atom();
        auto atmLinker_j = linkersCargos.model(modelId_j).chain("L").atom();
    
        if(dst(cargo_i_pos,cargo_j_pos)<(real(2.0)*cargoRadius)){return true;}
    
        //Linker i-Cargo j clash
        for(int n=0;n<int(atmLinker_i.size());n++){
            if(dst(atmLinker_i[n],cargo_j_pos)<(atmRadius+cargoRadius)){return true;}
        }
        
        //Linker j-Cargo j clash
        for(int m=0;m<int(atmLinker_j.size());m++){
            if(dst(atmLinker_j[m],cargo_i_pos)<(atmRadius+cargoRadius)){return true;}
        }
        
        //Linker i -Linker j clash
        for(int n=0;n<int(atmLinker_i.size());n++){
        for(int m=0;m<int(atmLinker_j.size());m++){
            if(dst(atmLinker_i[n],atmLinker_j[m])<(real(2.0)*atmRadius)){return true;}
        }}
    
        return false;
    }

    real square(const real& x){
        return x*x;
    }

    real computeLost(int modelId_i,int modelId_j,
                     STRUCTURE& linkersCargos,
                     real linkerCargoRadius, real cargoRadius, real atmRadius){
        
        real CCm2 = square(cargoRadius+cargoRadius);
        real CLm2 = square(atmRadius  +cargoRadius);
        real LLm2 = square(atmRadius  +atmRadius  );
    
        real lost = 0;
    
    
        /*
        for(CHAIN&   ch_i  : linkersCargos.model(modelId_i).chain()){
        for(RESIDUE& res_i : ch_i.residue()){
        for(ATOM&    atm_i : res_i.atom()){
        
            for(CHAIN&   ch_j  : linkersCargos.model(modelId_j).chain()){
            for(RESIDUE& res_j : ch_j.residue()){
            for(ATOM&    atm_j : res_j.atom()){
                
                if(atm_i.getResSeq() == 1 and atm_j.getResSeq() == 1){
                    if(dst(atm_i,atm_j)>(real(2.0)*linkerCargoRadius)){return real(0.0);}
                }
            
            
            }}}
    
    
        }}}*/
    
        //Linker i -Linker j lost
        for(RESIDUE& res_i : linkersCargos.model(modelId_i).chain("L").residue()){
        for(ATOM& atm_i : res_i.atom()){
        
            for(RESIDUE& res_j :linkersCargos.model(modelId_j).chain("L").residue()){
            for(ATOM& atm_j : res_j.atom()){
    
                if(atm_i.getResSeq() == 1 and atm_j.getResSeq() == 1){
                    if(dst(atm_i,atm_j)>(real(2.0)*linkerCargoRadius)){return real(0.0);}
                }
    
                if(atm_i.getResSeq() <=  atm_j.getResSeq()){
                    lost += square(std::max(real(0.0),LLm2-square(dst(atm_i,atm_j))));
                }
            }}
        
        }}
    
        real3 cargo_i_pos = computeCentroid(linkersCargos.model(modelId_i).chain("C"));
        real3 cargo_j_pos = computeCentroid(linkersCargos.model(modelId_j).chain("C"));
        
        //Linker i-Cargo j lost
        for(RESIDUE& res_i : linkersCargos.model(modelId_i).chain("L").residue()){
        for(ATOM& atm_i : res_i.atom()){
            lost += square(std::max(real(0.0),CLm2-square(dst(atm_i,cargo_j_pos))));
        }}
        
        //Linker j-Cargo i lost
        for(RESIDUE& res_j : linkersCargos.model(modelId_j).chain("L").residue()){
        for(ATOM& atm_j : res_j.atom()){
            lost += square(std::max(real(0.0),CLm2-square(dst(atm_j,cargo_i_pos))));
        }}
    
        //lost += square(std::max(real(0.0),CCm2-square(dst(cargo_i_pos,cargo_j_pos))));
        
        if(square(dst(cargo_i_pos,cargo_j_pos)) < CCm2){
    
            for(RESIDUE& res_i : linkersCargos.model(modelId_i).chain("C").residue()){
            for(ATOM& atm_i : res_i.atom()){
            
                for(RESIDUE& res_j :linkersCargos.model(modelId_j).chain("C").residue()){
                for(ATOM& atm_j : res_j.atom()){
    
                    if(atm_i.getResSeq() <=  atm_j.getResSeq()){
                        lost += square(std::max(real(0.0),LLm2-square(dst(atm_i,atm_j))));
                    }
                }}
            
            }}
    
        }
        
        return lost;
        
        /*
        auto atmLinker_i = linkersCargos.model(modelId_i).chain("L").atom();
        auto atmLinker_j = linkersCargos.model(modelId_j).chain("L").atom();
    
        if(dst(atmLinker_i[0],atmLinker_j[0])>(real(2.0)*linkerCargoRadius)){return lost;}
        
        real3 cargo_i_pos = computeCentroid(linkersCargos.model(modelId_i).chain("C"));
        real3 cargo_j_pos = computeCentroid(linkersCargos.model(modelId_j).chain("C"));
    
        //Cargo i - Cargo j lost
        
        lost += square(std::max(real(0.0),CCm2-square(dst(cargo_i_pos,cargo_j_pos))));
    
        //Linker i-Cargo j lost
        for(int n=0;n<int(atmLinker_i.size());n++){
            lost += square(std::max(real(0.0),CLm2-square(dst(atmLinker_i[n],cargo_j_pos))));
        }
        
        //Linker j-Cargo j clash
        for(int m=0;m<int(atmLinker_j.size());m++){
            lost += square(std::max(real(0.0),CLm2-square(dst(atmLinker_j[m],cargo_i_pos))));
        }
        
        //Linker i -Linker j clash
        for(int n=0;n<int(atmLinker_i.size());n++){
        for(int m=0;m<int(atmLinker_j.size());m++){
            lost += square(std::max(real(0.0),LLm2-square(dst(atmLinker_i[n],atmLinker_j[m]))));
        }}
    
        return lost;
        */
    }

    void updateLostMatrix(Eigen::MatrixXf& lostMatrix,STRUCTURE& linkersCargos,real linkerCargoRadius, real cargoRadius,real atmRadius){
    
        auto modelVector = linkersCargos.model();
        int nCargos = modelVector.size();
        
        for(int i=0;i<nCargos;i++){
        for(int j=0;j<nCargos;j++){
            lostMatrix(i,j)=0;
        }}
    
        for(int i=0  ;i<nCargos;i++){
        for(int j=i+1;j<nCargos;j++){
            real lost = computeLost(modelVector[i].getModelId(),
                                    modelVector[j].getModelId(),
                                    linkersCargos,linkerCargoRadius,
                                    cargoRadius,atmRadius);
            lostMatrix(i,j)=lost;
            lostMatrix(j,i)=lost;
        }}
    
    }

    void updateLostMatrix(std::map<int,int>& modelId2index,int modelId,Eigen::MatrixXf& lostMatrix,STRUCTURE& linkersCargos,real linkerCargoRadius, real cargoRadius,real atmRadius){
    
    
        for(MODEL& mdl : linkersCargos.model()){
            
            if(mdl.getModelId()==modelId){continue;}
    
            real lost = computeLost(modelId,mdl.getModelId(),linkersCargos,linkerCargoRadius,cargoRadius,atmRadius);
            lostMatrix(modelId2index[modelId],modelId2index[mdl.getModelId()])=lost;
            lostMatrix(modelId2index[mdl.getModelId()],modelId2index[modelId])=lost;
        }
    }

    int computeNumClashes(int modelId,STRUCTURE& linkersCargos, real cargoRadius,real atmRadius){
    
        auto modelVector = linkersCargos.model();
        int nCargos = modelVector.size();
    
        int nClashes=0;
        for(int i=0  ;i<nCargos;i++){
    
            if(modelVector[i].getModelId() == modelId){continue;}
            else {
                if(checkClash(modelId,modelVector[i].getModelId(),linkersCargos,cargoRadius,atmRadius)){nClashes++;}
            }
        }
    
        return nClashes;
    }
    
    int computeNumClashes(STRUCTURE& linkersCargos, real cargoRadius,real atmRadius){
    
        auto modelVector = linkersCargos.model();
        int nCargos = modelVector.size();
    
        int nClashes=0;
        for(int i=0  ;i<nCargos;i++){
        for(int j=i+1;j<nCargos;j++){
            if(checkClash(modelVector[i].getModelId(),modelVector[j].getModelId(),linkersCargos,cargoRadius,atmRadius)){nClashes++;}
        }}
    
        return nClashes;
    }

    void updateClashesMatrix(Eigen::MatrixXi& clashesMatrix,STRUCTURE& linkersCargos, real cargoRadius,real atmRadius){
    
        auto modelVector = linkersCargos.model();
        int nCargos = modelVector.size();
        
        for(int i=0;i<nCargos;i++){
        for(int j=0;j<nCargos;j++){
            clashesMatrix(i,j)=0;
        }}
    
        for(int i=0  ;i<nCargos;i++){
        for(int j=i+1;j<nCargos;j++){
            if(checkClash(modelVector[i].getModelId(),modelVector[j].getModelId(),linkersCargos,cargoRadius,atmRadius)){
                clashesMatrix(i,j)=1;
                clashesMatrix(j,i)=1;
            }
        }}
    
    }

    void updateClashesMatrix(int modelId,Eigen::MatrixXi& clashesMatrix,STRUCTURE& linkersCargos, real cargoRadius,real atmRadius){
    
        auto modelVector = linkersCargos.model();
        int nCargos = modelVector.size();
        
        int index = -1;
    
        for(int i=0;i<nCargos;i++){
            if(modelVector[i].getModelId()==modelId){index=i;};
        }
    
    
        for(int i=0  ;i<nCargos;i++){
            if(i==index){continue;}
            if(checkClash(modelId,modelVector[i].getModelId(),linkersCargos,cargoRadius,atmRadius)){
                clashesMatrix(index,i)=1;
                clashesMatrix(i,index)=1;
            } else {
                clashesMatrix(index,i)=0;
                clashesMatrix(i,index)=0;
    
            }
        }
    }

    bool checkIfCargoInsideContainer(STRUCTURE& linkersCargos,real3 containerCenter,real containerRadius){
        
        for(MODEL&   mdl : linkersCargos.model()){
        for(RESIDUE& res : mdl.chain("C").residue()        ){
        for(ATOM&    atm : res.atom()          ){
            if(dst(atm,containerCenter)>containerRadius){return false;}
        }}}
    
        return true;
    }

    bool checkIfCargoInsideContainer(int modelId,STRUCTURE& linkersCargos,real3 containerCenter,real containerRadius){
        
        for(RESIDUE& res : linkersCargos.model(modelId).chain("C").residue()        ){
        for(ATOM&    atm : res.atom()          ){
            if(dst(atm,containerCenter)>containerRadius){return false;}
        }}
    
        return true;
    }

    bool checkIfCargoInsideCone(int modelId,STRUCTURE& linkersCargos,real3 containerCenter,real containerRadius){
        
        real3 first = linkersCargos.model(modelId).chain("L").atom().front().getAtomCoord();
        real3 last  = linkersCargos.model(modelId).chain("L").atom().back().getAtomCoord();
    
        real3 v1 = containerCenter-first;
        real3 v2 = last-first;
    
        v1 = v1/sqrt(dot(v1,v1));
        v2 = v2/sqrt(dot(v2,v2));
    
        real angle = acos(dot(v1,v2));
        
        if(angle > M_PI/real(4)){return false;}
        
        return true;
    }


    void copyModelCoord(MODEL& dst,MODEL& src){
        
        for(CHAIN&   ch  : src.chain()       ){
        for(RESIDUE& res : ch.residue()      ){
        for(ATOM&   atm :  res.atom()        ){
            dst.chain(ch.getChainId()).residue(res.getResSeq()).atom(atm.getAtomName()).setAtomCoord(atm.getAtomCoord());
        }}}
    
    }

    void updateWeightVector(std::vector<real>& weightVector,Eigen::MatrixXf& lostMatrix){
    
        weightVector.resize(lostMatrix.rows());
    
        std::fill(weightVector.begin(),weightVector.end(),real(0.0));
    
        for(int i=0;i<lostMatrix.rows();i++){
        
            real w = 0;
            for(int j=0;j<lostMatrix.cols();j++){
                w+=lostMatrix(i,j);
            }
            
            weightVector[i]=w;
        }
    
        real min = *std::min_element(weightVector.begin(),weightVector.end());
        for(int i=0;i<weightVector.size();i++){
            if(weightVector[i]==0){weightVector[i]=min/real(10.0);}
        }
    }

    void minimizeCargoLost(std::mt19937& rndg,STRUCTURE& linkersCargos,real3 containerCenter, real containerRadius,real atmRadius){
    
        constexpr real maxAngle = M_PI/16;
    
        int nCargos = linkersCargos.model().size();
        Eigen::MatrixXf lostMatrix(nCargos,nCargos);
        Eigen::MatrixXf lostMatrixRef(nCargos,nCargos);
    
        std::uniform_int_distribution<> lCdistrib(0,linkersCargos.model().size()-1);
        std::uniform_real_distribution<real> angleDist(-maxAngle,maxAngle);
    
        std::vector<real> weightVector;
    
        MODEL refModel(0);
        
        for(CHAIN&   ch  : linkersCargos.model()[0].chain()       ){
            refModel.addChain(ch.getChainId());
        for(RESIDUE& res : ch.residue()      ){
            refModel.chain(ch.getChainId()).addResidue(res.getResName(),res.getResSeq(),res.getResInsCode());
        for(ATOM&   atm :  res.atom()        ){
            refModel.chain(ch.getChainId()).residue(res.getResSeq()).addAtom(atm.getAtomSerial(),atm.getAtomName());
        }}}
    
    
        real  cargoRadius;
        real linkerCargoRadius;
        {
            auto cargoVector = linkersCargos.model()[0].chain("C").atom();
            real3 cargoCentroid = computeCentroid(linkersCargos.model()[0].chain("C"));
    
            cargoRadius=0;
            for(auto& cargoAtm : cargoVector){
                cargoRadius=std::max(cargoRadius,dst(cargoAtm,cargoCentroid));
            }
            cargoRadius+=atmRadius;
            
            auto linkerVector = linkersCargos.model()[0].chain("L").atom();
            linkerCargoRadius=dst(linkerVector[0],cargoCentroid)+cargoRadius;
        }
    
        std::cerr << "Cargo radius: " << cargoRadius << std::endl;
        std::cerr << "Linker cargo radius: " << linkerCargoRadius << std::endl;
        
        std::vector<int>  index2modelId;
        std::map<int,int> modelId2index;
        {
            auto linkersCargosModelVector = linkersCargos.model();
    
            for(auto& m : linkersCargosModelVector){
                index2modelId.push_back(m.getModelId());
                modelId2index[m.getModelId()]=index2modelId.size()-1;
            }
    
        }
    
        int steps=0;
        updateLostMatrix(lostMatrix,linkersCargos,linkerCargoRadius,cargoRadius,atmRadius);
        real lost = lostMatrix.sum()/2;
    
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        while(lost>real(0.0)){
            if(steps%100==0){
                std::chrono::steady_clock::time_point check = std::chrono::steady_clock::now();
    
                real time = std::chrono::duration_cast<std::chrono::seconds> (check - start).count();
    
                std::cerr << "\r\e[K" <<  "Minimize cargo lost, step: " << steps << ", current lost: " << lost << " FPS: "<< steps/time << std::flush;
    
                if(steps%10000==0){
                    std::ofstream minLostBack("minLostBack.pdb");
                    minLostBack << linkersCargos << std::endl;
                    minLostBack.close();
                }
    
                //Print state
                /*
                auto linkersCargosVector = linkersCargos.atom();
                std::cout << "#" << std::endl;
                for(auto& lc : linkersCargosVector){
                    std::cout << lc.getAtomCoord() << " 3.8 "  << ((lc.getChainId()=="C")?255:512) << std::endl;
                }*/
            }
    
            steps++;
    
            //Random move
            
            
            //Select random linker-cargo
            
    
            //updateWeightVector(weightVector,lostMatrix);
    
            //std::discrete_distribution<> d(weightVector.begin(),weightVector.end());
            
            int cargoLinkerModelId = index2modelId[lCdistrib(rndg)];
            //int cargoLinkerModelId = index2modelId[d(rndg)];
            
            //Copy to ref
            
            copyModelCoord(refModel,linkersCargos.model(cargoLinkerModelId));
            lostMatrixRef = lostMatrix;
    
            //Rotation
                
            real3 lOrig = linkersCargos.model(cargoLinkerModelId).chain("L").atom().front().getAtomCoord();
    
            geometricTransformations::randomRotation(linkersCargos.model(cargoLinkerModelId),lOrig,angleDist(rndg),rndg);
    
            if(checkIfCargoInsideContainer(cargoLinkerModelId,linkersCargos,containerCenter,containerRadius)){
            //if(checkIfCargoInsideCone(cargoLinkerModelId,linkersCargos,containerCenter,containerRadius)){
                
                updateLostMatrix(modelId2index,cargoLinkerModelId,lostMatrix,linkersCargos,linkerCargoRadius,cargoRadius,atmRadius);
                real lostNew = lostMatrix.sum()/2;
                
                if(lostNew <= lost){
                    //Accept
                    lost = lostNew;
    
                } else {
                    //Reject
                    copyModelCoord(linkersCargos.model(cargoLinkerModelId),refModel);
                    lostMatrix = lostMatrixRef;
                }
                
            } else {
                //Reject
                copyModelCoord(linkersCargos.model(cargoLinkerModelId),refModel);
                lostMatrix = lostMatrixRef;
            }
    
        }
    
        std::cerr << std::endl;
    }

    void removeCargoClashes(std::mt19937& rndg,STRUCTURE& linkersCargos,real3 containerCenter, real containerRadius,real atmRadius){
    
        //constexpr real maxAngle = M_PI/16;
        constexpr real maxAngle = M_PI;
    
        int nCargos = linkersCargos.model().size();
        Eigen::MatrixXi clashesMatrix(nCargos,nCargos);
        Eigen::MatrixXi clashesMatrixRef(nCargos,nCargos);
    
        std::uniform_int_distribution<> lCdistrib(0,linkersCargos.model().size()-1);
        std::uniform_real_distribution<real> angleDist(-maxAngle,maxAngle);
    
        MODEL refModel(0);
        
        for(CHAIN&   ch  : linkersCargos.model()[0].chain()       ){
            refModel.addChain(ch.getChainId());
        for(RESIDUE& res : ch.residue()      ){
            refModel.chain(ch.getChainId()).addResidue(res.getResName(),res.getResSeq(),res.getResInsCode());
        for(ATOM&   atm :  res.atom()        ){
            refModel.chain(ch.getChainId()).residue(res.getResSeq()).addAtom(atm.getAtomSerial(),atm.getAtomName());
        }}}
    
    
        real  cargoRadius;
        {
            auto cargoVector = linkersCargos.model()[0].chain("C").atom();
            real3 cargoCentroid = computeCentroid(linkersCargos.model()[0].chain("C"));
    
            cargoRadius=0;
            for(auto& cargoAtm : cargoVector){
                cargoRadius=std::max(cargoRadius,dst(cargoAtm,cargoCentroid));
            }
        }
    
        std::vector<int>  index2modelId;
        std::map<int,int> modelId2index;
        {
            auto linkersCargosModelVector = linkersCargos.model();
    
            for(auto& m : linkersCargosModelVector){
                index2modelId.push_back(m.getModelId());
                modelId2index[m.getModelId()]=index2modelId.size()-1;
            }
    
        }
    
        int steps=0;
        updateClashesMatrix(clashesMatrix,linkersCargos,cargoRadius,atmRadius);
        int nClashes = clashesMatrix.sum()/2;
        while(nClashes>0){
            if(steps%100==0){
                std::cerr << "\r" <<  "Remove cargo clashes, step: " << steps << ", current number of clashes: " << nClashes << std::flush;
    
                //Print state
                /*
                auto linkersCargosVector = linkersCargos.atom();
                std::cout << "#" << std::endl;
                for(auto& lc : linkersCargosVector){
                    std::cout << lc.getAtomCoord() << " 3.8 "  << ((lc.getChainId()=="C")?255:512) << std::endl;
                }*/
            }
    
            steps++;
    
            //Random move
            
            
            //Select random linker-cargo
            
            int cargoLinkerModelId = index2modelId[lCdistrib(rndg)];
            
            //Copy to ref
            
            copyModelCoord(refModel,linkersCargos.model(cargoLinkerModelId));
            clashesMatrixRef = clashesMatrix;
    
            //Rotation
                
            real3 lOrig = linkersCargos.model(cargoLinkerModelId).chain("L").atom().front().getAtomCoord();
    
            geometricTransformations::randomRotation(linkersCargos.model(cargoLinkerModelId),lOrig,angleDist(rndg),rndg);
    
            if(checkIfCargoInsideContainer(linkersCargos,containerCenter,containerRadius)){
                updateClashesMatrix(cargoLinkerModelId,clashesMatrix,linkersCargos,cargoRadius,atmRadius);
                int nClashesNew = clashesMatrix.sum()/2;
                if(nClashesNew <= nClashes){
                    //Accept
                    nClashes = nClashesNew;
    
                } else {
                    //Reject
                    copyModelCoord(linkersCargos.model(cargoLinkerModelId),refModel);
                    clashesMatrix = clashesMatrixRef;
                }
            } else {
                //Reject
                    copyModelCoord(linkersCargos.model(cargoLinkerModelId),refModel);
                    clashesMatrix = clashesMatrixRef;
            }
    
        }
    
        std::cerr << std::endl;
    }

}
