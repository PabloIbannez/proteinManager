#include <proteinManager/proteinManager.hpp>

using namespace proteinManager;

int main(int argc, char *argv[]){
    
    STRUCTURE proteinIn;
    
    proteinIn.loadPDB(argv[1]);
    
    for(MODEL& md : proteinIn.model()){
    for(CHAIN& ch : md.chain()){
    for(RESIDUE& res : ch.residue()){
    for(ATOM& atm : res.atom()){
        std::cout << atm.getAtomName() << " " << atm.getAtomSASA() << std::endl;
    }}}}
    
    return EXIT_SUCCESS;
}
