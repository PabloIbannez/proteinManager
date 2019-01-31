/*Pablo Ibáñez Freire 2018, pablo.ibannez@uam.es*/

#include <iostream>

#include "../pIntegrator.cuh"

struct functionTest{
    
    real A_ = 1;
    
    void setA(real A){A_=A;}
    
    __host__ __device__ real operator()(real3 p){
        return p.x*p.y*p.z*sin(A_*p.x)*cos(A_*p.y*p.z);
    }
    
};


int main(int argc, char *argv[]){
    
    std::cout << "Test start" << std::endl;
    
    integrator::integrator3D_grid_ns::integrator3D_grid integ;
    
    integ.init({-1,-1,-1},{2,2,2},100);
    
    functionTest f;
    
    for(int i=0; i<std::atoi(argv[1]); i++){
        f.setA(i);
        std::cout << "Integral: " << integ.computeIntegral<functionTest>(f) << std::endl;
    }
    
    integ.free();
    
    std::cout << "Test end" << std::endl;
    
    
    return 0;
}
