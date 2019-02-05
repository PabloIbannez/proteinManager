/*Pablo Ibáñez Freire 2018, pablo.ibannez@uam.es*/

#include <iostream>

#include "../integrator3D.cuh"
#include <proteinManager.hpp>

using namespace proteinManager;

struct functionTest{
    
    __host__ __device__ real operator()(real3 p){
        return p.x*p.y*p.z*sin(p.x)*cos(p.y*p.z);
    }
    
};


int main(int argc, char *argv[]){
    
    std::cout << "Test start" << std::endl;
    
    integrator::integrator3D_MC_ns::integrator3D_MC integ;
    
    integ.init({-1,-1,-1},{2,2,2},std::atoi(argv[1]));
    
    functionTest f;
    
    real2 integral = integ.computeIntegralAverage(f,std::atoi(argv[2]));
    
    std::cout << "Integral: " << integral.x << " +- " << integral.y << std::endl;
    
    integ.free();
    
    std::cout << "Test end" << std::endl;
    
    return 0;
}
