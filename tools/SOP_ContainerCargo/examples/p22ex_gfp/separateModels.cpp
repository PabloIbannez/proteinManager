#include <proteinManager/proteinManager.hpp>

int main(int argc, char *argv[]){
    proteinManager::STRUCTURE pdbIn;
                                          
    pdbIn.loadPDB(std::string(argv[1])+".pdb");
    pdbIn.renumber();
    //pdbIn.setOutputFormat(proteinManager::DATA_FORMAT::PDB);
    
    int modelCount = 0;
    for(proteinManager::MODEL& md : pdbIn.model()){
        std::ofstream outModel(std::string(argv[1])+"_model_"+std::to_string(modelCount)+".pdb");
        outModel << md << std::endl;
        outModel.close();
        modelCount++;
    }
}
