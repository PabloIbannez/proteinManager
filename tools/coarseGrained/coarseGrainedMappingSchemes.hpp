#ifndef COARSE_GRAINED_MAPPING_SCHEMES_HPP
#define COARSE_GRAINED_MAPPING_SCHEMES_HPP

#include <random>

namespace proteinManager{
namespace coarseGrainedManager{
namespace coarseGrainedMappingSchemes{
            
    struct ca{
                                    
        void mappingScheme(RESIDUE& resIn, RESIDUE& resOut, std::string const & beadName,std::vector<std::string>& beadComponents){
            
            ////////////////////////////////////////////////
            
            real3 pos = resIn.atom("CA").getAtomCoord();
            
            ////////////////////////////////////////////////
            
            resOut.atom(beadName).setAtomCoord(pos);
            
            //Common properties
            resOut.atom(beadName).setAtomAltLoc(" ");
            resOut.atom(beadName).setAtomOccupancy(1);
            resOut.atom(beadName).setAtomTempFactor(0);
            resOut.atom(beadName).setAtomElement("");
        }
                                    
    };
    
    struct centerOfMass{
                                    
        void mappingScheme(RESIDUE& resIn, RESIDUE& resOut, std::string const & beadName,std::vector<std::string>& beadComponents){
            
            ////////////////////////////////////////////////

            real3 pos = {0,0,0};
            real totalMass = 0;
            
            for(std::string const & atom : beadComponents){
                
                pos       += resIn.atom(atom).getAtomCoord()*resIn.atom(atom).getAtomMass();
                totalMass += resIn.atom(atom).getAtomMass();
            }
            
            pos = pos/totalMass;
            
            ////////////////////////////////////////////////
            
            resOut.atom(beadName).setAtomCoord(pos);
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
        real cutOff = 2;
        
        real fieldMaxValue = 10;
        real fieldMinValue = 0.025;
        
        std::random_device rndDevice;
        std::mt19937 rndGen;
        std::uniform_real_distribution<real> unifDist;
        
        real densityOfPoints = 100;
        
        sasaFitting(){
            
            rndGen.seed(rndDevice());
            unifDist = std::uniform_real_distribution<real>(0,1);

        }
        
        real3 randomPointOverSphere(real radii){
            real lambda = std::acos(real(2)*unifDist(rndGen)-real(1))-real(1.0*M_PI_2);
            real phi    = real(2.0*M_PI)*unifDist(rndGen);
            
            return {radii*std::cos(lambda)*std::cos(phi),
                    radii*std::cos(lambda)*std::sin(phi),
                    radii*std::sin(lambda)};
        }
        
        real2 fittingBead(RESIDUE& resIn,std::vector<std::string>& beadComponents,real3 pos,real beadTotalSasa){
            
            if(beadTotalSasa == 0) {return {0,0};}
            
            ////////////////////////////////////////////////////////////
            
            int vectorSize     = (int(cutOff/dr));
            std::vector<real2> vdwData(vectorSize);
            std::vector<real2> stericData(vectorSize);
            
            for(int i=0;i<vectorSize;i++){
                vdwData[i] = {dr*i+dr,.0};
                stericData[i] = {dr*i+dr,.0};
            }
            
            //Compute field
            for(int i=0;i<vectorSize;i++){
                real radii = dr*i+dr;
                int numOfPoints = int(densityOfPoints*real(4.0*M_PI)*radii*radii);
                
                /*
                #ifdef DEBUG
                    std::cout << "Number of points to be generated " << numOfPoints << " " << resIn.getResName() << std::endl;
                #endif
                */
                
                for(int j=0;j<numOfPoints;j++){
                    
                    real valueC6  = 0;
                    real valueC12 = 0;
                    
                    for(std::string const & atom : beadComponents){
                        
                        real3 rap = (resIn.atom(atom).getAtomCoord()-pos)-randomPointOverSphere(radii);
                        real  invr2   = 1.0/dot(rap,rap);
                        real  invr6   = invr2*invr2*invr2;
                        real  invr12  = invr6*invr6;
                        
                        valueC6  += resIn.atom(atom).getAtomC6()*resIn.atom(atom).getAtomSASA()*invr6;
                        valueC12 += resIn.atom(atom).getAtomC12()*resIn.atom(atom).getAtomSASA()*invr12;
                        
                        //valueC6  += resIn.atom(atom).getAtomC6()*invr6;
                        //valueC12 += resIn.atom(atom).getAtomC12()*invr12;
                        
                    }
                    
                    vdwData[i].y    += valueC6/beadTotalSasa;
                    stericData[i].y += valueC12/beadTotalSasa;
                    
                    //vdwData[i].y    += valueC6;
                    //stericData[i].y += valueC12;
                
                }
                
                vdwData[i].y/=real(numOfPoints);
                stericData[i].y/=real(numOfPoints);
            }
            
            //Check if all field values are zero or -nan
            int k;
            for(k=0;k<vectorSize;k++){
                if(vdwData[k].y > 0) {break;}
                if(stericData[k].y > 0) {break;}
            }
            if(k==vectorSize) {return {0,0};}
            
            /*
            std::cout << this->fittingPowerLaw<6>(vdwData) << std::endl;
            std::cout << this->fittingPowerLaw<12>(stericData) << std::endl;
            
            std::ofstream test("test.dat");
            for(int i=0;i<vectorSize;i++){
                test << dr*i+dr << " " << vdwData[i].y << " " << stericData[i].y << std::endl;
            }
            test.close();
            std::cin.get();
            */
            ////////////////////////////////////////////////////////////
            
            return {this->fittingPowerLaw<6>(vdwData),this->fittingPowerLaw<12>(stericData)};
        }
        
        template <int exponent>
        real fittingPowerLaw(std::vector<real2>& data){
            
            int n = 0;
            real sumLogX = 0;
            real sumLogY = 0;
            
            for(real2& point : data){
                
                if(point.y > fieldMinValue and point.y < fieldMaxValue){

                    sumLogX += std::log(point.x);
                    sumLogY += std::log(point.y);
                    n++;
                }
            }
            
            return exp((sumLogY + exponent*sumLogX)/n);
            
        }
                                    
        void mappingScheme(RESIDUE& resIn, RESIDUE& resOut, std::string const & beadName,std::vector<std::string>& beadComponents){
            
            ///////////////////////////////////////////////
            
            real totalMass  = 0;
            real totalCharge = 0;
            real totalSasa  = 0;
            
            real weighSolvE = 0;
            
            real3 pos = {0,0,0};
            
            for(std::string const & atom : beadComponents){
                
                totalMass    += resIn.atom(atom).getAtomMass();
                totalCharge  += resIn.atom(atom).getAtomCharge();
                totalSasa    += resIn.atom(atom).getAtomSASA();
                weighSolvE   += resIn.atom(atom).getAtomSASA()*resIn.atom(atom).getAtomSolvE();
                pos          += resIn.atom(atom).getAtomCoord()*resIn.atom(atom).getAtomMass();;

            }
            
            pos /= totalMass;
            //pos = resIn.atom("CA").getAtomCoord();
            
            //#ifdef DEBUG
            //    std::cout << "Fitting res: " << resIn.getResName() << " " << resIn.getResSeq() << std::endl;
            //#endif
            
            real2 c6_c12 = fittingBead(resIn,beadComponents,pos,totalSasa);
            
            //real C6 = c6_c12.x;
            //real C12 = c6_c12.y;
            //
            //real sigma;
            //real epsilon;
            //
            //if(C6 == 0 or C12 == 0){
            //    sigma = 0;
            //    epsilon = 0;
            //} else {
            //    sigma   = std::cbrt(C12/C6);
            //    epsilon = real(0.25)*(C6/(sigma*sigma*sigma))*(C6/(sigma*sigma*sigma));
            //}
            //std::cout << sigma << " " << epsilon << " "<< totalSasa << std::endl;

            resOut.atom(beadName).setAtomC6(c6_c12.x);
            resOut.atom(beadName).setAtomC12(c6_c12.y);

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
