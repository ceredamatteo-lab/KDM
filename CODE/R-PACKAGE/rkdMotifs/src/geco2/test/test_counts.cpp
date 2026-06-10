#include <geco2/gKmers.h>
#include <geco2/gKDMotifs.h>

#include <fstream>

using namespace geco;
using namespace geco::methods;
using namespace geco::methods::kdmotifs;
using namespace geco::methods::kmers;
using namespace std;

#include <omp.h>


typedef double tipo;

int main_provo(int argc, char ** argv){
    auto prime_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;;
    
    if(argc<5){
        cout << "usage: test_counts nthreads l k ncases" << endl;
        return -1;
    }
    int nthreads=atol(argv[1]);
    if(nthreads!=0){
        omp_set_num_threads(nthreads);
    }
    indexType l=atol(argv[2]);
    indexType k=atol(argv[3]);
    indexType ncases=atol(argv[4]);
    
    vector<string> sequences2(ncases,string("AAGCTACGACTACGATCACGATCAGCACTAGCATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGATAGTA"));
    //vector<string> sequences2(ncases,string("AAGCTACGACTACGATCACGATCAGCACTAGNATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGATAGTA"));
    
    cout << "nthreads=" << nthreads << endl;
    cout << "l=" << l << endl;
    cout << "k=" << k << endl;
    cout << "nseqs=" << ncases << endl;
    
    bool strict=false;

    cout << "Countin lmers SS ..." << endl;
    try{
        gLmerCounter<tipo> counter(l,false);
        
        cout << "\tdense..." << endl;
        auto S1=counter.countDense(sequences2,strict);
        cout << "\tsparse..." << endl;
        auto S3=counter.countSparse(sequences2,strict);
    }catch(gMethodException &e){
        cout << e.what() << endl;
    }catch(...){
        cout << "something bad" << endl;
    }
    
    cout << "Countin lmers DS ..." << endl;
    try{
        gLmerCounter<tipo> counter(l,true);
        
        cout << "\tdense..." << endl;
        auto S1=counter.countDense(sequences2,strict);
        cout << "\tsparse..." << endl;
        auto S3=counter.countSparse(sequences2,strict);
    }catch(gMethodException &e){
        cout << e.what() << endl;
    }catch(...){
        cout << "something bad" << endl;
    }

    
    cout << "Countin gkmers SS ..." << endl;
    try{
        gGappedKmerCounter<tipo> counter(l,k,false);
        
        cout << "\tdense..." << endl;
        auto S1=counter.countDense(sequences2,strict);
        cout << "\tsparse..." << endl;
        auto S3=counter.countSparse(sequences2,strict);
    }catch(gMethodException &e){
        cout << e.what() << endl;
    }catch(...){
        cout << "something bad" << endl;
    }
    
    cout << "Countin gkmers DS ..." << endl;
    try{
        gGappedKmerCounter<tipo> counter(l,k,true);
        
        cout << "\tdense..." << endl;
        auto S1=counter.countDense(sequences2,strict);
        cout << "\tsparse..." << endl;
        auto S3=counter.countSparse(sequences2,strict);
    }catch(gMethodException &e){
        cout << e.what() << endl;
    }catch(...){
        cout << "something bad" << endl;
    }

    cout << "Estimating kmers SS ..." << endl;
    try{
        gLmerEstimator<tipo> counter(l,k,false);
        
        cout << "\tdense..." << endl;
        auto S1=counter.countDense(sequences2,strict);
        cout << "\tsparse..." << endl;
        auto S3=counter.countSparse(sequences2,strict);
    }catch(gMethodException &e){
        cout << e.what() << endl;
    }catch(...){
        cout << "something bad" << endl;
    }
    
    cout << "Estimating kmers DS ..." << endl;
    try{
        gLmerEstimator<tipo> counter(l,k,true);
        
        cout << "\tdense..." << endl;
        auto S1=counter.countDense(sequences2,strict);
        cout << "\tsparse..." << endl;
        auto S3=counter.countSparse(sequences2,strict);
    }catch(gMethodException &e){
        cout << e.what() << endl;
    }catch(...){
        cout << "something bad" << endl;
    }
    
    return 0;
}


void outMat(const gMDense<tipo> & M){
    cout << "[" << M.i_nrows << "x" << M.i_ncols << "]" << endl;
    cout << "ka" << "\tAA" << "\tAC" << "\tAG" << "\tAT" << "\tCA" << "\tCC" << "\tCG" << "\tCT" << "\tGA" << "\tGC" << "\tGG" << "\tGT" << "\tTA" << "\tTC" << "\tTG" << "\tTT" << endl ;
    for(indexType w=0;w<M.i_ncols;w++){
        cout << "w" << w;
        tipo sum=0;
        for(indexType kmer=0;kmer<M.i_nrows;kmer++){
            tipo val=*(M.i_values.i_ptr+M.i_nrows*w + kmer);
            cout << "\t" << val;
            sum+=val;
        }
        cout << "\t-->\t" << sum << endl;
    }
}


tipo check(const gMDense<tipo> & M1,const gMDense<tipo> & M2){
    tipo err=0;

    for(indexType w=0;w<M1.i_ncols;w++){
        for(indexType kmer=0;kmer<M1.i_nrows;kmer++){
            tipo val1=*(M1.i_values.i_ptr+M1.i_nrows*w + kmer);
            tipo val2=*(M2.i_values.i_ptr+M2.i_nrows*w + kmer);
            err+=abs(val1-val2);
        }
    }
    return err;
}


int main_provo2(int argc, char ** argv){
    auto prime_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;;
    
    indexType kmerlength=8;
    indexType ngaps=2;
    bool strict=false;
    indexType win=31;
    bool dstrand=true;
    
    indexType nrep=100;
    
    bool sparse=false;

    
    string seq="AAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAA";
    //string seq="AAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAAAGCTACGACTACGATCACGATCAGCACTAGTATCAGCTACGACTAGCATCAGCTACGACTAGCATCAGCATGCATCGATCAGCTACGAATCAGCATCAGCTACGACGTGGTGATGTGTAGTGATGTAGATAAGATAGATGATGAATGTAA";
    //string seq="AAGCTAC";
    
    indexType nwin = seq.size()-win+1;
    vector<string> seqs(nwin);
    for(indexType i=0;i<nwin;i++){
        seqs[i]=seq.substr(i,win);
    }
    gLmerCounter<tipo> counter1(kmerlength,dstrand);
    
    pair<vector<string>,vector<string>> strings1=counter1.getStrings();
    
    gGappedKmerCounter<tipo> counter2(kmerlength,kmerlength-ngaps,dstrand);
    
    
/*
    indexType nmotifs=109;
    vector<string> wseqs(nmotifs);
    for(indexType i=0;i<nmotifs;i++){
        wseqs[i]=seq.substr(1,150);
    }
    gMDense<tipo> W=toDense(kdmDistr<tipo>(wseqs,kmerlength,ngaps,dstrand,false,true));
    
    

    prime_time = std::chrono::high_resolution_clock::now();
    auto S1=counter1.countWinWeighted(seq,win,W,strict,true);
    dur=std::chrono::high_resolution_clock::now()-prime_time;
    cout << "sparse\t" << S1.i_nrows << " x " << S1.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;
    
//     prime_time = std::chrono::high_resolution_clock::now();
//     auto S2=counter1.countWinWeighted(seq,win,W,strict,false);
//     dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "dense\t" << S2.i_nrows << " x " << S2.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;
    
    

    prime_time = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
    for(indexType i=0;i<nrep;i++){
        auto S=counter1.countWinWeighted(seq,win,W,strict,true);
        cout << i << endl;
    }
    dur=std::chrono::high_resolution_clock::now()-prime_time;
    cout << dur.count() << endl;
*/
        
/*
    cout <<"seqlen=" << seq.length() << endl;    
    try{
        
        cout << "lmers check---->DENSE" << endl;
        {
            prime_time = std::chrono::high_resolution_clock::now();
            auto S1=counter1.countWinDense(seq,win,strict);
            dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "win\t" << S1.i_nrows << " x " << S1.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;
            
            prime_time = std::chrono::high_resolution_clock::now();
            vector<string> seqs(nwin);
            for(indexType i=0;i<nwin;i++){
                seqs[i]=seq.substr(i,win);
            }
            auto S2=counter1.countDense(seqs,strict);
            dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "dir\t" << S2.i_nrows << " x " << S2.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;
            cout << "ERR: " << check(S1,S2) << endl;
        }

        cout << "lmers check---->SPARSE" << endl;
        {
            prime_time = std::chrono::high_resolution_clock::now();
            auto S1=counter1.countWinSparse(seq,win,strict);
            dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "win\t" << S1.i_nrows << " x " << S1.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;
            
            prime_time = std::chrono::high_resolution_clock::now();
            vector<string> seqs(nwin);
            for(indexType i=0;i<nwin;i++){
                seqs[i]=seq.substr(i,win);
            }
            auto S2=counter1.countSparse(seqs,strict);
            dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "dir\t" << S2.i_nrows << " x " << S2.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;
            cout << "ERR: " << check(toDense(S1),toDense(S2)) << endl;
        }
        
        cout << "gkmers check---->DENSE" << endl;
        {
            prime_time = std::chrono::high_resolution_clock::now();
            auto S1=counter2.countWinDense(seq,win,strict);
            dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "win\t" << S1.i_nrows << " x " << S1.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;
            
            prime_time = std::chrono::high_resolution_clock::now();
            vector<string> seqs(nwin);
            for(indexType i=0;i<nwin;i++){
                seqs[i]=seq.substr(i,win);
            }
            auto S2=counter2.countDense(seqs,strict);
            dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "dir\t" << S2.i_nrows << " x " << S2.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;
            cout << "ERR: " << check(S1,S2) << endl;
        }

        cout << "gkmers check---->SPARSE" << endl;
        {
            prime_time = std::chrono::high_resolution_clock::now();
            auto S1=counter2.countWinSparse(seq,win,strict);
            dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "win\t" << S1.i_nrows << " x " << S1.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;

            prime_time = std::chrono::high_resolution_clock::now();
            vector<string> seqs(nwin);
            for(indexType i=0;i<nwin;i++){
                seqs[i]=seq.substr(i,win);
            }
            auto S2=counter2.countSparse(seqs,strict);
            dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "dir\t" << S2.i_nrows << " x " << S2.i_ncols << "\t" << dur.count()*1000 << "ms" << endl;
            cout << "ERR: " << check(toDense(S1),toDense(S2)) << endl;
        }



        {
            cout << "lmers speed dense...";cout.flush();        
            {
                cout << "\twin...";cout.flush();
                prime_time = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
                for(int i=0;i<nrep;i++){        
                    auto S1=counter1.countWinDense(seq,win,strict);
                }
                dur=std::chrono::high_resolution_clock::now()-prime_time;
                double d1=dur.count();
                cout << "\tdir...";cout.flush();
                prime_time = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
                for(int i=0;i<nrep;i++){
                    vector<string> seqs(nwin);
                    for(indexType i=0;i<nwin;i++){
                        seqs[i]=seq.substr(i,win);
                    }
                    auto S2=counter1.countDense(seqs,strict);
                }
                dur=std::chrono::high_resolution_clock::now()-prime_time;
                double d2=dur.count();
                cout << "\t" << d1 << "\t" << d2 << endl;
            }

            cout << "lmers speed sparse...";cout.flush();        
            {
                cout << "\twin...";cout.flush();
                prime_time = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
                for(int i=0;i<nrep;i++){        
                    auto S1=counter1.countWinSparse(seq,win,strict);
                }
                dur=std::chrono::high_resolution_clock::now()-prime_time;
                double d1=dur.count();
                
                cout << "\tdir...";cout.flush();
                prime_time = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
                for(int i=0;i<nrep;i++){
                    vector<string> seqs(nwin);
                    for(indexType i=0;i<nwin;i++){
                        seqs[i]=seq.substr(i,win);
                    }
                    auto S2=counter1.countSparse(seqs,strict);
                }
                dur=std::chrono::high_resolution_clock::now()-prime_time;
                double d2=dur.count();
                cout << "\t" << d1 << "\t" << d2 << endl;
            }
            
            cout << "gkmers speed dense...";cout.flush();        
            {
                cout << "\twin...";cout.flush();
                prime_time = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
                for(int i=0;i<nrep;i++){        
                    auto S1=counter2.countWinDense(seq,win,strict);
                }
                dur=std::chrono::high_resolution_clock::now()-prime_time;
                double d1=dur.count();
                cout << "\tdir...";cout.flush();
                prime_time = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
                for(int i=0;i<nrep;i++){
                    vector<string> seqs(nwin);
                    for(indexType i=0;i<nwin;i++){
                        seqs[i]=seq.substr(i,win);
                    }
                    auto S2=counter2.countDense(seqs,strict);
                }
                dur=std::chrono::high_resolution_clock::now()-prime_time;
                double d2=dur.count();
                cout << "\t" << d1 << "\t" << d2 << endl;
            }

            cout << "gkmers speed sparse...";cout.flush();
            {
                cout << "\twin...";cout.flush();
                prime_time = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
                for(int i=0;i<nrep;i++){        
                    auto S1=counter2.countWinSparse(seq,win,strict);
                }
                dur=std::chrono::high_resolution_clock::now()-prime_time;
                double d1=dur.count();
                cout << "\tdir...";cout.flush();
                prime_time = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
                for(int i=0;i<nrep;i++){
                    vector<string> seqs(nwin);
                    for(indexType i=0;i<nwin;i++){
                        seqs[i]=seq.substr(i,win);
                    }
                    auto S2=counter2.countSparse(seqs,strict);
                }
                dur=std::chrono::high_resolution_clock::now()-prime_time;
                double d2=dur.count();
                cout << "\t" << d1 << "\t" << d2 << endl;
            }
        }

        
    }catch(gMethodException &e){
        cout << e.what() << endl;
    }catch(...){
        cout << "something bad" << endl;
    }
    
*/    
    return 0;
}


int main(int argc, char ** argv){
    //return main_vectors(argc,argv);
    return main_provo(argc,argv);
}





