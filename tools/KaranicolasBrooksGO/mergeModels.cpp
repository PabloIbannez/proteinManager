#include <proteinManager/proteinManager.hpp>

int main(int argc, char *argv[]){
    proteinManager::STRUCTURE pdbIn;
                                          
    std::vector<std::string> charList =  {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z",
                                          "a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z",
                                          "0","1","2","3","4","5","6","7","8","9"};
    
    pdbIn.loadPDB(argv[1]);
    pdbIn.renumber();
    //pdbIn.setOutputFormat(proteinManager::DATA_FORMAT::PDB);
    
    int chainCount = 0;
    for(proteinManager::MODEL& md : pdbIn.model())      {
    for(proteinManager::CHAIN& ch : md.chain())         {
        ch.setChainId(charList[chainCount]);
        chainCount++;
    }}

    for(auto& ch : pdbIn.chain()){
        std::cout << ch << std::endl;
    }
}
