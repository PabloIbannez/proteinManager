namespace proteinManager{
namespace ffManager{
    
    struct repulsive12{
        
        real  C12_;
        real3 partPos_;
        
        real  cutOff2_;
        real e2_;
        
        repulsive12(real cutOff,real epsilon){
            cutOff2_ = cutOff*cutOff;
            e2_ = epsilon*epsilon;
        }
        
        real getCutOff() { return sqrt(cutOff2_);}
        real getEpsilon(){ return sqrt(e2_);}
        
        bool isNull(){
            if(C12_ == 0) { return true;}
            else {return false;}
        }
        
        void setParameters(ATOM& atm){
            C12_ = atm.getAtomC12();
            partPos_ = atm.getAtomCoord();
        }
        
        void setParametersLinealFit(ATOM& atm){
            C12_ = 1;
            partPos_ = atm.getAtomCoord();
        }
        
        __host__ __device__ inline real operator()(real3 point){
        
            real r2 = dot(point-partPos_,point-partPos_);
            
            if(r2 > cutOff2_) return real(0);
            
            return C12_/(powf(r2,int(6))); // C12_/r^12
        }
    };
    
    struct attractive6{
        
        real  C6_;
        real3 partPos_;
        
        real  cutOff2_;
        real e2_;
        
        attractive6(real cutOff,real epsilon){
            cutOff2_ = cutOff*cutOff;
            e2_ = epsilon*epsilon;
        }
        
        real getCutOff() { return sqrt(cutOff2_);}
        real getEpsilon(){ return sqrt(e2_);}
        
        bool isNull(){
            if(C6_ == 0) { return true;}
            else {return false;}
        }
        
        void setParameters(ATOM& atm){
            C6_ = atm.getAtomC6();
            partPos_ = atm.getAtomCoord();
        }
        
        __host__ __device__ inline real operator()(real3 point){
        
            real r2 = dot(point-partPos_,point-partPos_);
            
            if(r2 > cutOff2_) return real(0);
            
            return C6_/(powf(r2,int(3))); // C6_/r^6
        }
    };
    
    struct coulomb{
        
        real  C_;
        real3 partPos_;
        
        real  cutOff2_;
        real e2_;
        
        coulomb(real cutOff,real epsilon){
            cutOff2_ = cutOff*cutOff;
            e2_ = epsilon*epsilon;
        }
        
        real getCutOff() { return sqrt(cutOff2_);}
        real getEpsilon(){ return sqrt(e2_);}
        
        bool isNull(){
            if(C_ == 0) { return true;}
            else {return false;}
        }
        
        void setParameters(ATOM& atm){
            C_ = atm.getAtomCharge();
            partPos_ = atm.getAtomCoord();
        }
        
        __host__ __device__ inline real operator()(real3 point){
        
            real r2 = dot(point-partPos_,point-partPos_);
            
            if(r2 > cutOff2_) return real(0);
            
            return C_/sqrtf(r2); // C_/r
        }
        
        
    };
    
    struct debye{
        
        real  C_;
        real3 partPos_;
        real  f_;
        real  dl_;
        real  e_s_;
        
        real  cutOff2_;
        real e2_;
        
        debye(real f,real dl,real e_s,real cutOff,real epsilon){
            f_ = f;
            dl_ = dl;
            e_s_ = e_s;
            cutOff2_ = cutOff*cutOff;
            e2_ = epsilon*epsilon;
        }
        
        real getCutOff() { return sqrt(cutOff2_);}
        real getEpsilon(){ return sqrt(e2_);}
        
        bool isNull(){
            if(C_ == 0) { return true;}
            else {return false;}
        }
        
        void setParameters(ATOM& atm){
            C_ = atm.getAtomCharge();
            partPos_ = atm.getAtomCoord();
        }
        
        __host__ __device__ inline real operator()(real3 point){
        
            real r2 = dot(point-partPos_,point-partPos_);
            
            if(r2 > cutOff2_) return real(0);
            
            real r = sqrtf(r2);
            return f_*(C_*exp(-r/dl_))/(e_s_*r);
        }
        
    };
    
    //potential compositions
    
    template <class potType1,class potType2>
    struct potProduct{
        
        potType1 pt1_;
        potType2 pt2_;
        
        potProduct(potType1 pt1, potType2 pt2){
            pt1_ = pt1;
            pt2_ = pt2;
        }
        
        void setParameters(ATOM& atm1,ATOM& atm2){
            pt1_.setParameters(atm1);
            pt2_.setParameters(atm2);
        }
        
        void setParametersLinealFit(ATOM& atm1,ATOM& atm2){
             pt1_.setParametersLinealFit(atm1);
             pt2_.setParametersLinealFit(atm2);
        }
        
        bool isNull(){
            return pt1_.isNull() or pt1_.isNull();
        }
        
        __host__ __device__ inline real operator()(real3 point){
            
            return pt1_(point)*pt2_(point);
            
        }
        
    };
    
    template <class potType>
    struct potPF_Product{
        
        potType pt_;
        
        potPF_Product(potType pt){
            pt_ = pt;
        }
        
        real getCutOff(){
            return pt_.getCutOff();
        }
        
        void setParameters(ATOM& atm){
            pt_.setParameters(atm);
        }
        
        __host__ __device__ inline real operator()(real3 point,real fvalue){
            
            return pt_(point)*fvalue;
            
        }
        
    };
    
}}
