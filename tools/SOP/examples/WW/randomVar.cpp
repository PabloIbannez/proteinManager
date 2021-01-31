#include <proteinManager/proteinManager.hpp>
#include <random>

using namespace proteinManager;

int main(int argc,char* argv[]){

    STRUCTURE pdbIn;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(-std::atof(argv[2]), std::atof(argv[2])); 

    pdbIn.loadPDB(argv[1]);

    for(MODEL& md : pdbIn.model()){
    for(CHAIN& ch : md.chain()){
    for(RESIDUE& res : ch.residue()){
    for(ATOM& atm : res.atom()){
        
        real3 atmCoord = atm.getAtomCoord();
        
        atmCoord.x += distr(gen);
        atmCoord.y += distr(gen);
        atmCoord.z += distr(gen);
        
        atm.setAtomCoord(atmCoord);
        
        //std::cout << atm << std::endl;

    }}}}

    std::cout << pdbIn << std::endl;

    return EXIT_SUCCESS;
}
