#include "rGeco.h"



static void checkRInterruptFn(void* dummy) { R_CheckUserInterrupt(); }

void gMethodRInterrupt() {
    // Rcpp::Rcout << "checkInterrupt called" << std::endl;
    if (R_ToplevelExec(checkRInterruptFn, nullptr) == FALSE){
        // Rcpp::Rcout << "interrupt detected, throwing gMethodsInterrupt" << std::endl;
        throw geco::methods::gMethodsInterrupt();
    }
}    
    
void gMethodRWarnings(const std::string & message){
    Rcpp::warning(message);
}



//[[Rcpp::init]]
void initPackage(DllInfo* dll){
  geco::methods::setWarningCallback(&gMethodRWarnings);
  geco::methods::setInterruptCallback(&gMethodRInterrupt);
  geco::methods::setgMethodsSeed(std::time(NULL));
}

//[[Rcpp::export]]
void setCudaUsage(bool use_cuda){
  geco::methods::setCudaUsage(use_cuda);
}

//[[Rcpp::export]]
void setNCores(int ncores){
  geco::methods::setNCores(ncores);
}

//[[Rcpp::export]]
void setVerbosity(unsigned short verbosity){
  geco::methods::setVerbosity(verbosity);
}

//[[Rcpp::export]]
void setgMethodsSeed(float seed){
  geco::methods::setgMethodsSeed(seed);
}
