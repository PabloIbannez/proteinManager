#ifndef COARSE_GRAINED_MAPPING_SCHEMES_HPP
#define COARSE_GRAINED_MAPPING_SCHEMES_HPP

#include <random>

namespace proteinManager{
namespace coarseGrainedManager{
namespace coarseGrainedMappingSchemes{
            
    struct ca{
                                    
        void mappingScheme(RESIDUE& resIn, RESIDUE& resOut, std::string const & beadName,std::vector<std::string>& beadComponents){
            
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
    
    struct sasa{
                                    
        void mappingScheme(RESIDUE& resIn, RESIDUE& resOut, std::string const & beadName,std::vector<std::string>& beadComponents){
            
            ////////////////////////////////////////////////

            real3 pos = resIn.atom("CA").getAtomCoord();
            
            real totalMass  = 0;
            real totalSasa  = 0;
            real weighSolvE = 0;
            
            for(std::string const & atom : beadComponents){
                
                totalMass  += resIn.atom(atom).getAtomMass();
                totalSasa  += resIn.atom(atom).getAtomSASA();
                weighSolvE += resIn.atom(atom).getAtomSASA()*resIn.atom(atom).getAtomSolvE();
            }
            
            pos = resIn.atom("CA").getAtomCoord();
            
            ////////////////////////////////////////////////
            
            resOut.atom(beadName).setAtomCoord(pos);
            resOut.atom(beadName).setAtomMass(totalMass);
            resOut.atom(beadName).setAtomSASA(totalSasa);
            if(totalSasa == 0.0){
                resOut.atom(beadName).setAtomSolvE(0.0);
            } else {
                resOut.atom(beadName).setAtomSolvE(weighSolvE/totalSasa);
            }
            //Common properties
            resOut.atom(beadName).setAtomAltLoc(" ");
            resOut.atom(beadName).setAtomOccupancy(1);
            resOut.atom(beadName).setAtomTempFactor(0);
            resOut.atom(beadName).setAtomElement("");
        }
                                    
    };
    
    struct sasaFitting{
        
        real dr = 0.01;
        real cutOff = 5;
        
        std::vector<real2> dstFieldValue;
        
        std::random_device rndDevice;
        std::mt19937 rndGen;
        std::uniform_real_distribution<real> unifDist;
        
        real densityOfPoints = 100;
        
        sasaFitting(){
            dstFieldValue.resize(int(cutOff/dr));
            
            for(int i=0;i<dstFieldValue.size();i++){
                dstFieldValue[i].x = dr*i+dr;
                dstFieldValue[i].y = real(0);
            }
            
            rndGen.seed(rndDevice());
            unifDist = std::uniform_real_distribution<real>(0,1);
            
            #ifdef DEBUG
                std::cout << "Fitting vector size: " << dstFieldValue.size() << std::endl;
                
                for(int i=0;i<dstFieldValue.size();i++){
                    std::cout << dstFieldValue[i].x << " " << dstFieldValue[i].y << std::endl;
                }
                
            #endif
        }
        
        real3 randomPointOverSphere(real radii){
            real lambda = std::acos(real(2)*unifDist(rndGen)-real(1))-real(1.0*M_PI_2);
            real phi    = real(2.0*M_PI)*unifDist(rndGen);
            
            return {radii*std::cos(lambda)*std::cos(phi),
                    radii*std::cos(lambda)*std::sin(phi),
                    radii*std::sin(lambda)};
        }
        
        real vdw(std::vector<ATOM>& atoms,real3 point){
            
            real value = 0;
            
            for(ATOM& atm : atoms){
                real3 rap = atm.getAtomCoord()-point;
                value += atm.getAtomC6()*std::pow(1.0/dot(rap,rap),real(3.0));
                //std::cout << atm.getAtomC6() << " " << atm.getAtomCoord() << " " << point << " " << value << std::endl;
            } 
            
            return value;
        }
        
        real fittingBead(RESIDUE& resIn,std::vector<std::string>& beadComponents){
            
            #ifdef DEBUG
                std::cout << "Fitting res: " << resIn.getResSeq() << std::endl;
            #endif
            
            real3 beadCenterOfMass = {0,0,0};
            real  beadTotalMass    = 0;
            
            std::vector<ATOM> beadAtoms;
            
            for(std::string const & atom : beadComponents){
                
                beadCenterOfMass += resIn.atom(atom).getAtomCoord()*resIn.atom(atom).getAtomMass();
                beadTotalMass    += resIn.atom(atom).getAtomMass();
                
                beadAtoms.push_back(resIn.atom(atom));
            }
            
            beadCenterOfMass = beadCenterOfMass/beadTotalMass;
            
            for(ATOM& atm : beadAtoms){
                atm.setAtomCoord(atm.getAtomCoord()-beadCenterOfMass);
            }
            
            for(real2& point : dstFieldValue){
                point.y = 0;
            }
            
            for(real2& point : dstFieldValue){
                real radii = point.x;
                int numOfPoints = int(densityOfPoints*real(4.0*M_PI)*radii*radii);
                for(int i=0;i<numOfPoints;i++){
                    #ifdef DEBUG
                        std::cout << "Number of points to be generated " << numOfPoints << " " << resIn.getResName() << std::endl;
                    #endif
                    point.y += vdw(beadAtoms,randomPointOverSphere(radii));
                }
                point.y/=real(numOfPoints);
            }
            
            std::ofstream test("test.dat");
            for(real2& point : dstFieldValue){
                test << point.x << " " << point.y << std::endl;
            }
            test.close();
            
            std::cin.get();
            
        }
                                    
        void mappingScheme(RESIDUE& resIn, RESIDUE& resOut, std::string const & beadName,std::vector<std::string>& beadComponents){
            
            ////////////////////////////////////////////////

            real3 pos = resIn.atom("CA").getAtomCoord();
            
            real totalMass  = 0;
            real totalCharge = 0;
            real totalSasa  = 0;
            real weighSolvE = 0;
            
            real test = fittingBead(resIn,beadComponents);
            
            
            for(std::string const & atom : beadComponents){
                
                totalMass   += resIn.atom(atom).getAtomMass();
                totalCharge += resIn.atom(atom).getAtomCharge();
                totalSasa   += resIn.atom(atom).getAtomSASA();
                weighSolvE  += resIn.atom(atom).getAtomSASA()*resIn.atom(atom).getAtomSolvE();
            }
            
            pos = resIn.atom("CA").getAtomCoord();
            
            ////////////////////////////////////////////////
            
            resOut.atom(beadName).setAtomCoord(pos);
            resOut.atom(beadName).setAtomMass(totalMass);
            resOut.atom(beadName).setAtomCharge(totalCharge);
            resOut.atom(beadName).setAtomSASA(totalSasa);
            if(totalSasa == 0.0){
                resOut.atom(beadName).setAtomSolvE(0.0);
            } else {
                resOut.atom(beadName).setAtomSolvE(weighSolvE/totalSasa);
            }
            //Common properties
            resOut.atom(beadName).setAtomAltLoc(" ");
            resOut.atom(beadName).setAtomOccupancy(1);
            resOut.atom(beadName).setAtomTempFactor(0);
            resOut.atom(beadName).setAtomElement("");
        }
                                    
    };
}}}

#endif
