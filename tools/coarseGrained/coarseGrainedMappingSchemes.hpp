#ifndef COARSE_GRAINED_MAPPING_SCHEMES_HPP
#define COARSE_GRAINED_MAPPING_SCHEMES_HPP

namespace proteinManager{
namespace coarseGrainedManager{
namespace coarseGrainedMappingSchemes{
            
    struct ca{
                                    
        static void mappingScheme(RESIDUE& resIn, RESIDUE& resOut, std::string const & beadName,std::vector<std::string>& beadComponents){
            
            ////////////////////////////////////////////////

            real3 pos = {0,0,0};
            real chg = 0;
            real totalMass = 0;
            
            for(std::string const & atom : beadComponents){
                
                chg += resIn.atom(atom).getAtomCharge();
                totalMass += resIn.atom(atom).getAtomMass();
            }
            
            pos = resIn.atom("CA").getAtomCoord();
            
            ////////////////////////////////////////////////
            
            resOut.atom(beadName).setAtomCoord(pos);
            resOut.atom(beadName).setAtomCharge(chg);
            resOut.atom(beadName).setAtomMass(totalMass);
            
            //Common properties
            resOut.atom(beadName).setAtomAltLoc(" ");
            resOut.atom(beadName).setAtomOccupancy(1);
            resOut.atom(beadName).setAtomTempFactor(0);
            resOut.atom(beadName).setAtomElement("");
        }
                                    
    };
}}}

#endif
