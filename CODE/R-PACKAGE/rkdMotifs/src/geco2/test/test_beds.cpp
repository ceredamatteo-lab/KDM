#include <geco2/base.h>
#include <geco2/gBed.h>
#include <fstream>

using namespace std;
using namespace geco;
using namespace geco::bed;

struct keyCalc{ 
    
public:
    unsigned int i_gcbins;
    unsigned int i_rpbins;
    
    keyCalc(unsigned int gcbins,unsigned int rpbins):i_gcbins(gcbins),i_rpbins(rpbins){
    }
    
    inline unsigned int operator() (const std::string & s){
        unsigned int len=s.length();
        unsigned int gc_count=0;
        unsigned int rp_count=0;
        
        for(auto c:s){
            rp_count+=islower(c)!=0;
            char tc=tolower(c);
            if(tc=='c' || tc=='g'){
                gc_count++;
            }
        }
        unsigned int gc_bin = (gc_count * i_gcbins) / len;
        unsigned int rp_bin = (rp_count * i_rpbins) / len;
        return gc_bin * (i_rpbins+1) + rp_bin;
    }
};



// struct keyCalc{ 
// public:
//     unsigned int i_gcbins;
//     unsigned int i_rpbins;
//     
//     keyCalc(unsigned int gcbins,unsigned int rpbins):i_gcbins(gcbins),i_rpbins(rpbins){
//     }
//     
//     inline unsigned int operator() (const string & s){
//         unsigned int len=s.length();
//         unsigned int gc_count=0;
//         unsigned int rp_count=0;
//         
//         for(auto c:s){
//             rp_count+=islower(c)!=0;
//             char tc=tolower(c);
//             if(tc=='c' || tc=='g'){
//                 gc_count++;
//             }
//         }
//         //unsigned int gc_bin = (gc_count * i_gcbins) / len;
//         //unsigned int rp_bin = (rp_count * i_rpbins) / len;
//         //return ( (gc_count * i_gcbins) / len) * (i_rpbins+1) + (rp_count * i_rpbins) / len;
//         return (gc_count * i_gcbins * (i_rpbins+1) + (rp_count * i_rpbins) ) / len;
//     }
// };

unsigned int rcount(const string & s){
    unsigned int  sm=0;
    for(auto c:s){
        sm+=islower(c)!=0;
    }
    return sm;
}

unsigned int gccount(const string & s){
    unsigned int sm=0;
    for(auto c:s){
        char tc=tolower(c);
        if(tc=='c' || tc=='g'){
            sm++;
        }
    }
    return sm;
}




int main(int argc, const char *argv[]){
       
    gBedData<gNarrowPeak> myFile(argv[1]);
//gBedData<gNarrowPeak> myFile(argv[1]);
    gSize n=myFile.size(),cc;
    
//     cout << " unsorted n_peaks:\t" << myFile.size() << endl;    
//     cc=0;    
//     for(auto e:myFile){
//         cout << e << endl;
//         if(cc++>=n) break;
//     }
   
    
    myFile.sort([](const gNarrowPeak & e1, const gNarrowPeak & e2){return e1.i_start>e2.i_start;});
    cout << " sorted n_peaks:\t" << myFile.size() << endl;
    geco::sequence_retrievers::gLocal2bitSequenceRetriever retriever("/adat/database/2bit/hg38.2bit",true);
    
    n=200;
    cc=0;
    for(auto e:myFile){
        gSequence seq=gElement(e.i_chrom,gReferenceInterval(e.i_start,e.i_end),e.i_isForward,retriever,gElm).getSequence();
        cout << e << "\t" << seq << endl;
        if(cc++>=n) break;
    }
    
    
//     vector<gBedData<gNarrowPeak>> selected=removeOverlappingEntries<gNarrowPeak>(myFile);
//     cout << " not overlapping n_peaks:\t" << selected[0].size() << endl;    
//     cc=0;        
//     for(auto e:selected[0]){
//         cout << e << endl;
//         if(cc++>=n) break;
//     }
//     cout << "overlapping n_peaks:\t" << selected[1].size() << endl;
//     cc=0;        
//     for(auto e:selected[1]){
//         cout << e << endl;
//         if(cc++>=n) break;
//     }
// 
    
//     gSize nsel=2000;
//     gSize npeaks=1000;    
//     gSize gcbins=20;
//     gSize rpbins=20;
//     
//     
//     gBedData<gNarrowPeak> selected2;
//     cc=0;
// 
//     for(auto e:selected[0]){
//         selected2.push_back(e);
//         if(cc>=nsel-1) break;
//         cc++;
//     }
//     cout << "selected n_peaks:\t" << selected2.size() << endl;
//     
//     geco::sequence_retrievers::gLocal2bitSequenceRetriever retriever("/adat/database/2bit/hg38.2bit",true);
//     
// //     regionFinder<gNarrowPeak> finder(selected2,myFile,retriever,gcbins,rpbins);
// //     gBedData<gBedBaseEntry> negatives=finder(npeaks);
// //     cout << "negatives n_peaks:\t" << negatives.size() << endl;
//     
// //     regionFinder2<gNarrowPeak,gNarrowPeak> finder(selected2,myFile,retriever,gcbins,rpbins);
// //     pair<gBedData<gNarrowPeak>,gBedData<gBedBaseEntry>> results=finder(npeaks);
// //     cout << "positive n_peaks:\t" << results.first.size() << endl;
// //     cout << "negative n_peaks:\t" << results.second.size() << endl;
// //     
// //     keyCalc getKey(gcbins,rpbins);
// //     auto pos=results.first.begin();
// //     auto neg=results.second.begin();
// //     
// //     ofstream tout("/tmp/results.txt");
// //     while(pos!=results.first.end()){
// //         string s1=gElement(pos->i_chrom,gReferenceInterval(pos->i_start,pos->i_end),true,retriever,gElm).getSequence();
// //         string s2=gElement(neg->i_chrom,gReferenceInterval(neg->i_start,neg->i_end),true,retriever,gElm).getSequence();
// //         tout << *pos << "\t" << pos->i_end - pos->i_start << "\t" << getKey(s1) << "\t" << gccount(s1) << "\t" << rcount(s1) << "\t" << *neg << "\t" << neg->i_end - neg->i_start << "\t" << getKey(s2) << "\t" << gccount(s2) << "\t" << rcount(s2) << endl;
// //         pos++;
// //         neg++;
// //     }
// //     tout.close();
    
    return 0;
}

