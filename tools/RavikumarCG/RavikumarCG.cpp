#include <cmath>
#include <fstream>
#include <vector>

#include <proteinManager/proteinManager.hpp>

struct particle{
    
    int id;
    std::string name;
    int mol;
    proteinManager::real3 pos;
};

struct bond{
  int i,j; 
  double k, r0;
};

std::ostream & operator << (std::ostream &out, bond b){
  out<<b.i<<" "<<b.j<<" "<<b.k<<" "<<b.r0;
  return out;
}

int main(int argc, char* argv[]){
    
    proteinManager::STRUCTURE pdbIn;
    proteinManager::STRUCTURE pdbExt;
    
    std::vector<particle> particleVector;
    
    pdbIn.loadPDB(argv[1]);
    
    if(std::string(argv[2]) != "none"){
	pdbExt.loadPDB(argv[2]);
    }
        
    std::ofstream topFile;
    std::ofstream bondsFile;
    
    topFile.open(argv[3]);
    bondsFile.open(argv[4]);
    
    int               id = 0;
    std::string       type;
    int               mol;
    proteinManager::real3 pos;
    proteinManager::real  chg;
    proteinManager::real  SASA;

    for(proteinManager::MODEL& md : pdbIn.model()) {
	for(proteinManager::CHAIN& ch : md.chain()) {
	    for(proteinManager::RESIDUE& res : ch.residue()) {
		
		type = res.getResName();
		mol = md.getModelId();
		
		try{
		    pos = res.atom("CA").getAtomCoord()/10.0;
		} catch(const std::runtime_error& e){
		    std::cerr << e.what() << std::endl;
		    std::cerr << "WARNING: I will consider the center of mass of this residue as the CG bead position." << std::endl << std::endl << std::endl;
		    
		    pos = {0,0,0};
		    for(auto a : res.atom()){
			pos += a.getAtomCoord()/10.0;
		    }
		    
		    pos /= res.atom().size();
		}
		particleVector.push_back({id,res.getResName(),mol,pos});
		
		SASA = 0;
		for(proteinManager::ATOM& atm : res.atom()){
		    SASA += atm.getAtomSASA();
		}
		
		
		if      (res.getResName() == "LYS" or res.getResName() == "ARG"){ chg =  1; }
		else if (res.getResName() == "ASP" or res.getResName() == "GLU"){ chg = -1; }
		else if (res.getResName() == "HIS"){ chg = 0.5; }
		else    { chg = 0; }
		
		topFile   << std::left << std::fixed << std::setw(6) 
			  << id         << " " 
			  << type       << " " 
			  << mol        << " " 
			  << pos.x      << " " 
			  << pos.y      << " " 
			  << pos.z      << " " 
			  << chg        << " " 
			  << SASA/100.0 << std::endl;
		id++;
	    }
	}
    }
    
    //ENM
    
    std::vector<bond> bondsVector;
    
    double r_max = 1;
    double k = 400;
    
    int i,j;
    for(int i = 0; i < particleVector.size() ; i++){
	for(int j = i + 1; j < particleVector.size() ; j++){
	    
	    double r;
	    
	    r = (particleVector[i].pos.x - particleVector[j].pos.x)*(particleVector[i].pos.x - particleVector[j].pos.x) +
	        (particleVector[i].pos.y - particleVector[j].pos.y)*(particleVector[i].pos.y - particleVector[j].pos.y) +
		(particleVector[i].pos.z - particleVector[j].pos.z)*(particleVector[i].pos.z - particleVector[j].pos.z);
	    
	    r = std::sqrt(r);
	    
	    if(r < r_max and particleVector[i].mol == particleVector[j].mol) {
		bondsVector.push_back({i,j,k,r});
	    }
			   
	}
    }
    
    bondsFile << bondsVector.size() << std::endl;
    for(bond b: bondsVector){
	bondsFile << b << std::endl;
    }
    
    if(std::string(argv[2]) == "none"){
	return EXIT_SUCCESS;
    }
    
    ////////////////////////////////////////////////////////////////////
    
    std::vector<particle> fixedParticleVector;
    
    for(proteinManager::MODEL& md : pdbExt.model()) {
	for(proteinManager::CHAIN& ch : md.chain()) {
	    for(proteinManager::RESIDUE& res : ch.residue()) {
		
		type = res.getResName();
		mol = md.getModelId();
		
		pos = res.atom("CA").getAtomCoord()/10.0;
		fixedParticleVector.push_back({id,res.getResName(),mol,pos});
		
		SASA = 0;
		for(proteinManager::ATOM& atm : res.atom()){
		    SASA += atm.getAtomSASA();
		}
		
		
		if      (res.getResName() == "LYS" or res.getResName() == "ARG"){ chg =  1; }
		else if (res.getResName() == "ASP" or res.getResName() == "GLU"){ chg = -1; }
		else if (res.getResName() == "HIS"){ chg = 0.5; }
		else    { chg = 0; }
		
		topFile   << std::left << std::fixed << std::setw(6) 
			  << id         << " " 
			  << type       << " " 
			  << -1         << " " 
			  << pos.x      << " " 
			  << pos.y      << " " 
			  << pos.z      << " " 
			  << chg        << " " 
			  << SASA/100.0 << std::endl;
		id++;
	    }
	}
    }
    
    bondsFile << fixedParticleVector.size() << std::endl;
    for(particle p: fixedParticleVector){
	bondsFile << p.id << " " << p.pos.x << " " << p.pos.y << " " << p.pos.z << " " << k << " " << 0 << std::endl;
	
    }
    
    return EXIT_SUCCESS;
}
