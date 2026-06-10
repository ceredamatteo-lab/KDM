#pragma once
#include <iostream>
#include <iostream>
#include <cstring>
#include <vector>
#include <map>
#include <mutex>
#include <omp.h>

#include <bits/posix_opt.h>

#define THRUST_IGNORE_CUB_VERSION_CHECK


//#define GMETHODS_USE_32BITS_INDEX
  
#define GMETHOD_INTERRUPTIBLE  


namespace geco{
    namespace methods{
        
#if defined(GMETHODS_USE_32BITS_INDEX)
#if   UINT_MAX  >= 0xffffffff
  typedef unsigned int      indexType;
#elif defined(UINT32_MAX)
  typedef          uint32_t indexType;
#else
  #error "don't know how to typedef 'indexType' on this system"
#endif
#else
#if   ULLONG_MAX >= 0xffffffffffffffff
  typedef unsigned long long indexType;
#elif defined(UINT64_MAX)
//  typedef          uint64_t  indexType;
  typedef unsigned long long indexType;
#else
    #error "don't know how to typedef 'indexType' on this system"
#endif
#endif
  

        class gMethodException:public std::exception{
        private:
            std::string i_msg;
        public:
            
            gMethodException(const std::string & msg):i_msg(msg){
            }
            
            virtual const char* what() const noexcept{
                return i_msg.c_str();
            }
            
        };
        

        
        class gMethodsInterrupt:public std::exception{
            virtual const char* what() const noexcept{
                return "User interruption";
            }
        };
        
        typedef void (*gMethodInterruptFunction)();
        typedef void (*gMethodWarningFunction)(const std::string & message);
        void setNCores(int ncores);
        void setgMethodsSeed(unsigned int seed);
        unsigned int getgMehtodsSeed();
        void setInterruptCallback(gMethodInterruptFunction function);
        void setWarningCallback(gMethodWarningFunction function);
        void setVerbosity(unsigned short verbosity);
        void setCudaUsage(bool use_cuda);
        
        bool useCuda();
        unsigned short getVerbosity();
        
        void warning(const std::string & message);
        
        void checkInterrupt(int tn=0);
        
        
        class gMethodsOut{
        private:
            unsigned short i_level;    
            std::chrono::time_point<std::chrono::high_resolution_clock> i_start;
        public:
            gMethodsOut(unsigned short level):i_level(level){
                if(i_level<=getVerbosity()){
                    i_start=std::chrono::high_resolution_clock::now();
                }
            }
            
            template <typename Function, typename... Parameters>
            void Start(Function f, Parameters... params){
                if(i_level<=getVerbosity()){
                    i_start=(std::chrono::high_resolution_clock::now());
                    f(params...);
                    std::cout.flush();
                }
            }
            
            template <typename Function, typename... Parameters>
            void End(Function f, Parameters... params){
                if(i_level<=getVerbosity()){
                    std::chrono::duration<double> dur = std::chrono::high_resolution_clock::now() - i_start;
                    f(params...);
                    std::cout << "\t(" << dur.count() << ")" << std::endl;
                }
            }
        };
        
		class gMThreadException {
			std::exception_ptr Ptr;
			std::mutex         Lock;
		public:
			gMThreadException();
			~gMThreadException();
			void Rethrow();
			void CaptureException();
			bool exceptionOccurred();
		};
        
        template <typename T>
        class gMMemory{
        public:
            const T * i_ptr;
            bool i_owns;
            indexType i_nelem;
            
            static inline T* acquire(const indexType n_elem){
                if(n_elem == 0)  { return nullptr; }        
                T* out_memptr;
    #if ( defined(_POSIX_ADVISORY_INFO) && (_POSIX_ADVISORY_INFO >= 200112L) )
                {
                    T* memptr = nullptr;
                    const size_t n_bytes   = sizeof(T)*size_t(n_elem);
                    const size_t alignment = (n_bytes >= size_t(1024)) ? size_t(32) : size_t(16);
                    // FROM ARMADILLO
                    // TODO: investigate apparent memory leak when using alignment >= 64 (as shown on Fedora 28, glibc 2.27)
                    int status = posix_memalign((void **)&memptr, ( (alignment >= sizeof(void*)) ? alignment : sizeof(void*) ), n_bytes);
                    out_memptr = (status == 0) ? memptr : nullptr;
                }
    #elif defined(_MSC_VER)
                {
                    const size_t n_bytes   = sizeof(T)*size_t(n_elem);
                    const size_t alignment = (n_bytes >= size_t(1024)) ? size_t(32) : size_t(16);
                    
                    out_memptr = (T *) _aligned_malloc( n_bytes, alignment );
                }
    #else
                {
                    out_memptr = (T *) malloc(sizeof(T)*n_elem);
                }
    #endif
                if(out_memptr == nullptr){
                    throw gMethodException("gMethods: out of memory");
                }
            
                return out_memptr;
            }

            static inline void release(T* mem){
                if(mem == nullptr)  { return; }        
#if ( defined(_POSIX_ADVISORY_INFO) && (_POSIX_ADVISORY_INFO >= 200112L) )            
                {
                    free( (void *)(mem) );
                }
#elif defined(_MSC_VER)
                {
                    //free( (void *)(mem) );
                    _aligned_free( (void *)(mem) );
                }
#else
                {
                //delete [] mem;
                    free( (void *)(mem) );
                }
#endif
            }
            
            gMMemory(){
                i_ptr=nullptr;
                i_nelem=0;
                i_owns=false;
            }
            
            gMMemory(indexType n_elem):gMMemory(){
                if(n_elem>0){
//                     std::cout << "Memory allocated" << std::endl;
                    i_ptr=acquire(n_elem);
                    i_nelem=n_elem;
                    i_owns=true;
                }                        
            }
            
            template<typename D>
            gMMemory(const D * otherPtr,indexType n_elem):gMMemory(){
                copyFrom(otherPtr,n_elem);
            }

            gMMemory(T * otherPtr,indexType n_elem):gMMemory(){
                i_ptr=otherPtr;
                i_nelem=n_elem;
                i_owns=false;
            }
            
            gMMemory(const gMMemory<T> & other):gMMemory(){
                //std::cout << "gMMemory: copy constructed" << std::endl;
                copyFrom(other.i_ptr,other.i_nelem);
            }
            
            gMMemory<T> & operator = (gMMemory<T> other){
                //std::cout << "gMMemory: copy assigned" << std::endl;
                swap(*this,other);
                return *this;
            }
            
            gMMemory(gMMemory<T> && other):gMMemory(){
                //std::cout << "gMMemory: move constructed" << std::endl;
                swap(*this,other);
            }
            
//             gMMemory<T> & operator = (gMMemory<T> && other){
//                 std::cout << "gMMemory: move assigned" << std::endl;
//                 swap(*this,other);
//                 return *this;
//             }
            
            
            ~gMMemory(){
                if(i_ptr!=nullptr){
                    if(i_owns){
//                         std::cout << "Memory released" << std::endl;
                        release(const_cast<T *>(i_ptr));
                    }
                    i_ptr=nullptr;
                    i_nelem=0;
                    i_owns=false;
                }
            }
            
            operator T * (){
                return const_cast<T *>(i_ptr);
            }
            
            operator const T * () const {
                return i_ptr;
            }
            
            template<typename D>
            void copyFrom(const D * ptr,indexType n_elem,indexType add=0){
                if(i_owns){
                    if(i_ptr!=nullptr && i_nelem!=(n_elem+add)){
                        std::cout << "i_ptr!=nullptr && i_nelem!=n_elem (release)" << std::endl;
                        release(const_cast<T *>(i_ptr));
                        i_ptr=nullptr;
                        i_nelem=0;
                        i_owns=false;
                    }
                }else{
                    i_ptr=nullptr;
                    i_nelem=0;
                }
                
                if(n_elem>0 && ptr!=nullptr){
                    if(i_ptr==nullptr){
                        i_ptr=acquire(n_elem+add);
                        i_nelem=n_elem+add;
                        i_owns=true;
                    }
#pragma omp parallel for
                    for (indexType i=0;i<n_elem;i++){
                        const_cast<T &>(i_ptr[i])=ptr[i];
                    }
//                     std::copy(ptr,ptr+n_elem,const_cast<T*>(i_ptr));
                }
            }
            
            friend void swap(gMMemory<T> & a,gMMemory<T> & b){
                using std::swap;
                swap(a.i_ptr,b.i_ptr);
                swap(a.i_nelem,b.i_nelem);
                swap(a.i_owns,b.i_owns);
            }
            
        };
        
        template <typename T>
        class gMDense{
        public:
            gMMemory<T> i_values;
            indexType i_nrows;
            indexType i_ncols;
            
            gMDense():i_values(),i_nrows(0),i_ncols(0){
            }
            
            gMDense(indexType nrows,indexType ncols):i_values(nrows*ncols),i_nrows(nrows),i_ncols(ncols){
            }
            
            template<typename D>
            gMDense(indexType nrows,indexType ncols,const D * values):i_values(values,nrows*ncols),i_nrows(nrows),i_ncols(ncols){
            }
            
            
            //copy constructor
            gMDense(const gMDense<T> & other):i_values(other.i_values),i_nrows(other.i_nrows),i_ncols(other.i_ncols){
//                 std::cout << "gMDense: copy constructed" << std::endl;
            }
            
            //move constructor
            gMDense(gMDense<T> && other):gMDense(){
//                 std::cout << "gMDense: move constructed" << std::endl;                
                swap(*this,other);
            }
            
            //copy assignment
            gMDense<T> & operator = (gMDense<T> other){
//                 std::cout << "gMDense: copy assigned" << std::endl;
                swap(*this,other);
                return *this;
            }
            
            //move assignment
//             gMDense<T> & operator = (gMDense<T> && other){
// //                 std::cout << "gMDense: move assigned" << std::endl;
//                 swap(*this,other);
//                 return *this;
//             }
            
            friend void swap(gMDense<T> & a,gMDense<T> & b){
                using std::swap;
                swap(a.i_values,b.i_values);
                swap(a.i_nrows,b.i_nrows);
                swap(a.i_ncols,b.i_ncols);
                
            }
            
            friend std::ostream & operator << (std::ostream & out, const gMDense<T> & mat){
                for(indexType row=0;row<mat.i_nrows;row++){
                    for(indexType col=0; col<mat.i_ncols;col++){
                        out << "\t" << mat.i_values.i_ptr[row + col*mat.i_nrows];
                    }
                    out << std::endl;
                }
                return out;
            }
            
        };
        
        template<typename T>
        class gMSparse{
        public:
            gMMemory<T> i_cscValues;
            gMMemory<indexType> i_cscColPtr;
            gMMemory<indexType> i_cscRowIndex;
            indexType i_nrows;
            indexType i_ncols;
            indexType i_nnz;
            
            gMSparse():gMSparse(0,0,0){
            }
            
            gMSparse(indexType nrows,indexType ncols,indexType nnz=0):i_cscValues(nnz+1),i_cscColPtr(ncols+2),i_cscRowIndex(nnz+1),i_nrows(nrows),i_ncols(ncols),i_nnz(nnz){
                memset((indexType *) i_cscColPtr,0,sizeof(T)*(ncols+1));
                ((indexType *) i_cscColPtr)[ncols+1]=std::numeric_limits<indexType>::max();
                ((indexType *) i_cscRowIndex)[nnz]=0;                
                ((T *)i_cscValues)[nnz]=0;
            }
            
            gMSparse(const std::vector<T> & cscValues, const std::vector<indexType> & cscColPtr, const std::vector<indexType> & cscRowIndex, indexType nrows,indexType ncols):i_cscValues(cscValues.size()+1),i_cscColPtr(ncols+2),i_cscRowIndex(cscValues.size()+1),i_nrows(nrows),i_ncols(ncols),i_nnz(cscValues.size()){
                i_cscValues.copyFrom(cscValues.data(),i_nnz,1);
                i_cscValues[i_nnz]=0;
                i_cscRowIndex.copyFrom(cscRowIndex.data(),i_nnz,1);
                i_cscRowIndex[i_nnz]=0;
                i_cscColPtr.copyFrom(cscColPtr.data(),i_ncols+1,1);
                i_cscColPtr[i_ncols+1]=0;
            }
            
            gMSparse(const std::vector<std::map<indexType,T>> & valuesMap, indexType nrows):i_cscColPtr(valuesMap.size()+2),i_nrows(nrows),i_ncols(valuesMap.size()){
                i_nnz=0;
                i_cscColPtr[i_ncols+1]=std::numeric_limits<indexType>::max();
                i_cscColPtr[0]=0;                
                for(indexType c=0;c<i_ncols;c++){
                    i_cscColPtr[c+1]=i_cscColPtr[c]+valuesMap[c].size();
                    i_nnz+=valuesMap[c].size();
                }
                i_cscValues=gMMemory<T>(i_nnz+1);
                i_cscValues[i_nnz]=0;
                
                i_cscRowIndex=gMMemory<indexType>(i_nnz+1);
                i_cscRowIndex[i_nnz]=0;
                //omp_set_num_threads(70);
#pragma omp parallel for
                for(indexType c=0;c<i_ncols;c++){
                    indexType cp=i_cscColPtr[c];
                    for(auto e=valuesMap[c].begin();e!=valuesMap[c].end();e++){
                        i_cscValues[cp]=e->second;
                        i_cscRowIndex[cp++]=e->first;
                    }
                }

                
                    
                
                

//                 for(indexType c=0;c<i_ncols;c++){
//                     i_cscColPtr[c+1]=i_cscColPtr[c];
//                     indexType & cp=i_cscColPtr[c+1];
//                     for(auto e=valuesMap[c].begin();e!=valuesMap[c].end();e++){
//                         i_cscValues[cp]=e->second;
//                         i_cscRowIndex[cp]=e->first;
//                         cp++;
//                     }
//                 }
                
            }
            
            
            //copy construcotr
            gMSparse(const gMSparse<T> & other):i_cscValues(other.i_cscValues),i_cscColPtr(other.i_cscColPtr),i_cscRowIndex(other.i_cscRowIndex),i_nrows(other.i_nrows),i_ncols(other.i_ncols),i_nnz(other.i_nnz){
//                 std::cout << "gMSparse: copy constructed" << std::endl;
            }
            
            //move construcotr
            gMSparse(gMSparse<T> && other):gMSparse(){
//                 std::cout << "gMSparse: move constructed" << std::endl;
                swap(*this,other);
            }
            
            //copy assignment
            gMSparse<T> & operator = (gMSparse<T> other){
//                 std::cout << "gMSparse: copy assigned" << std::endl;
                swap(*this,other);
                return *this;
            }
            
            //move assignment is equivalent equals the copy contruct and swap assignment
//             gMSparse<T> & operator = (gMSparse<T> && other){
// //                 std::cout << "gMSparse: move assigned" << std::endl;
//                 swap(*this,other);
//                 return *this;
//             }
            
            friend void swap(gMSparse<T> & a,gMSparse<T> & b){
                using std::swap;
                swap(a.i_cscValues,b.i_cscValues);
                swap(a.i_cscColPtr,b.i_cscColPtr);
                swap(a.i_cscRowIndex,b.i_cscRowIndex);
                swap(a.i_nrows,b.i_nrows);
                swap(a.i_ncols,b.i_ncols);
                swap(a.i_nnz,b.i_nnz);
            }
            
        };
        
        template<typename T>
        gMSparse<T> toSparse(const gMDense<T> & dmat);
        
        template<typename T>
        gMDense<T> toDense(const gMSparse<T> & smat);
        
        template<typename T>
        T mse(const gMDense<T> & M1, const gMDense<T> & M2);
        
        template<typename T>
        std::vector<indexType> find(const gMDense<T> & M1, T threshold, bool greater=true);
        
        template<typename T>
        gMDense<T> select_columns(gMDense<T> M1, const std::vector<indexType> columns);
        
        template<typename T>
        std::vector<T> randperm(T n, T m);
        
        template<template<typename> class A, typename T>
        class fill_method{
        public:
            A<T> operator () (indexType nrows,indexType ncols,T value) const;
        };
        
        template<template<typename> class A,typename T>
        A<T> fill(indexType nrows,indexType ncols,T value){
            return fill_method<A,T>()(nrows,ncols,value);
        }

        template<template<typename> class A, typename T>
        class initRandom_method{
        public:
            A<T> operator () (indexType nrows,indexType ncols,T sparseness) const;
        };

        template<template<typename> class A,typename T>
        A<T> initRandom(indexType nrows,indexType ncols,T sparseness){
            return initRandom_method<A,T>()(nrows,ncols,sparseness);
        }

        
        template<template<typename> class A, template<typename> class B, typename T>
        class multiply_method   {
        public:
            B<T> operator () (const A<T> & M1,const B<T> & M2,bool transposeA ,bool transposeB) const;
        };
        
        template<template<typename> class A,template<typename> class B, typename T>
        B<T> multiply(const A<T> & M1, const B<T> & M2, bool transposeA=false ,bool transposeB = false){
            indexType Cr,Cc,k;
            if(!transposeA){
                Cr=M1.i_nrows;
                k=M1.i_ncols;
                if(!transposeB){
                    Cc=M2.i_ncols;
                    k*=(k==M2.i_nrows);
                }else{
                    Cc=M2.i_nrows;
                    k*=(k==M2.i_ncols);
                }
            }else{
                Cr=M1.i_ncols;
                k=M1.i_nrows;
                if(!transposeB){
                    Cc=M2.i_ncols;
                    k*=(k==M2.i_nrows);
                }else{
                    Cc=M2.i_nrows;
                    k*=(k==M2.i_ncols);
                }
            }
            if(k==0){
                //std::cout << "M1[" << M1.i_nrows << "," << M1.i_ncols << "]" << (transposeA?"'":"") << "x M2[" << M2.i_nrows << "," << M2.i_ncols << "]" << (transposeB?"'":"") << std::endl;
                throw gMethodException("Invalid matrix dimensions");
            }
            return multiply_method<A,B,T>()(M1,M2,transposeA,transposeB);
        }
        
        
        
    }
}
