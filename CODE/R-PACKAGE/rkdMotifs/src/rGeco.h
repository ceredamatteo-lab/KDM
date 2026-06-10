#include <Rcpp.h>
// #include <geco2/include/gMethods.h>
// #include <geco2/include/base.h>

#include <gMethods.h>
#include <base.h>

namespace geco{

  template<class T> Rcpp::NumericVector convert(const gArray<T> & array){
    Rcpp::NumericVector ret(array.getSize());
    for(gSize i=0;i<array.getSize();i++){
      if(array.isNA(i)){
        ret[i]=NA_REAL;
      }else{
        ret[i]=array[i];
      }
    }
    return ret;
  }
  
  template<class T> gArray<T> convert(const Rcpp::NumericVector & array){
    gArray<T> ret(0,array.size(),true);
    for(int i=0;i<array.size();i++){
      if(!Rcpp::NumericVector::is_na(array(i))){
        ret.setValue(i,array(i),false);
      }
    }
    return ret;
  }
  
  template<class T> Rcpp::NumericMatrix convert(const gMatrix<T> & matrix){
    Rcpp::NumericMatrix ret(matrix.getRowsNum(),matrix.getColsNum());
    for(gSize r=0;r<matrix.getRowsNum();r++){
      for(gSize c=0;c<matrix.getColsNum();c++){
        if(matrix.isNA(r,c)){
          ret(r,c)=NA_REAL;
        }else{
          ret(r,c)=matrix(r,c);
        }
      }
    }
    return ret;
  }
  
  template<class T> gMatrix<T> convert(const Rcpp::NumericMatrix & matrix){
    gMatrix<T> ret(matrix.rows(),matrix.cols(),0,true);
    for(int r=0;r<matrix.rows();r++){
      for(int c=0;c<matrix.cols();c++){
        if(!Rcpp::NumericMatrix::is_na(matrix(r,c))){
          ret.setValue(r,c,matrix(r,c),false);
        }
      }
    }
    return ret;
  }

  namespace methods{
    

    void gMethodRInterrupt(){
      Rcpp::checkUserInterrupt();
    }
    
    void gMethodRWarnings(const std::string & message){
      Rcpp::warning(message);
    }
  

    template<typename T>
    geco::methods::gMDense<T> convertFromR(const Rcpp::NumericMatrix & mat){
      if(mat.nrow()>0 && mat.ncol()>0){
        std::vector<T> values((size_t)mat.nrow()*(size_t)mat.ncol());
        std::copy(mat.begin(),mat.end(),values.begin());
        return geco::methods::gMDense<T>(mat.nrow(),mat.ncol(),&values[0]);
      }else{
        return geco::methods::gMDense<T>();
      }
    }

    //ATTENZIONE, questo potrebbe non funzionare per T=float, dato che std:copy non modifica il tipo e NumericMatrix è double
    template<typename T>
    Rcpp::NumericMatrix convertFromG(const geco::methods::gMDense<T> & mat){
      Rcpp::NumericMatrix ret(mat.i_nrows,mat.i_ncols);
      const T * mem=(const T *) mat.i_values;
      std::copy(mem,mem + (mat.i_nrows * mat.i_ncols),ret.begin());
      return ret;
    }
    
    template<typename T>
    geco::methods::gMSparse<T> convertFromR(const Rcpp::S4 & mat){
      if(!mat.is("dgCMatrix")){
        Rcpp::stop("Only sparse Matrix from package \"Matrix\" and class \"dgCMatrix\" can be used" );
      }
      
      Rcpp::IntegerVector Dim = mat.slot("Dim");
      Rcpp::IntegerVector i = mat.slot("i");
      Rcpp::IntegerVector p = mat.slot("p");
      Rcpp::NumericVector x = mat.slot("x");
      
      std::vector<T> cscValues((size_t)x.size());
      std::copy(x.begin(), x.end(), cscValues.begin());
      
      std::vector<geco::methods::indexType> cscColPtr((size_t)p.size());
      std::copy(p.begin(), p.end(), cscColPtr.begin());
      
      std::vector<geco::methods::indexType> cscRowIndex((size_t)i.size());
      std::copy(i.begin(), i.end(), cscRowIndex.begin());
      
      return geco::methods::gMSparse<T>(cscValues,cscColPtr,cscRowIndex,(geco::methods::indexType)Dim[0],(geco::methods::indexType)Dim[1]);
    }
    
    template<typename T>
    Rcpp::S4 convertFromG(const geco::methods::gMSparse<T> & mat, const std::string & className=std::string("")){
      Rcpp::IntegerVector dim = Rcpp::IntegerVector::create( mat.i_nrows, mat.i_ncols );  
      Rcpp::NumericVector x((const T *) mat.i_cscValues, (const T*) mat.i_cscValues + mat.i_nnz) ;
      Rcpp::IntegerVector p((const geco::methods::indexType *) mat.i_cscColPtr, (const geco::methods::indexType *) mat.i_cscColPtr + mat.i_ncols+1) ;  
      Rcpp::IntegerVector i((const geco::methods::indexType *) mat.i_cscRowIndex, (const geco::methods::indexType *) mat.i_cscRowIndex + mat.i_nnz);
      
      std::string klass="dgCMatrix";
      if(className.length()>0){
        klass=className;
      }
      
      Rcpp::S4 s( (className.length()==0)?("dgCMatrix"):(className));
      s.slot("i")   = i;
      s.slot("p")   = p;
      s.slot("x")   = x;
      s.slot("Dim") = dim;
      return s;
    }

    
    
    template<typename T, typename D>
    geco::methods::gMDense<T> convertFromR2(const D & mat){
      if(mat.nrow()>0 && mat.ncol()>0){
        std::vector<T> values((size_t)mat.nrow()*(size_t)mat.ncol());
        std::copy(mat.begin(),mat.end(),values.begin());
        return geco::methods::gMDense<T>(mat.nrow(),mat.ncol(),&values[0]);
      }else{
        return geco::methods::gMDense<T>();
      }
    }
    
    //ATTENZIONE, questo potrebbe non funzionare per T=float, dato che std:copy non modifica il tipo e NumericMatrix è double
    template<typename T,typename D>
    D convertFromG2(const geco::methods::gMDense<T> & mat){
      D ret(mat.i_nrows,mat.i_ncols);
      const T * mem=(const T *) mat.i_values;
      std::copy(mem,mem + (mat.i_nrows * mat.i_ncols),ret.begin());
      return ret;
    }    
        

  }
}

//[[Rcpp::export]]
void initPackage(){
  geco::methods::setWarningCallback(&geco::methods::gMethodRWarnings);
  geco::methods::setInterruptCallback(&geco::methods::gMethodRInterrupt);
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
