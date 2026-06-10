#pragma once
#include "array.h"
#include "2bit_retriever.h"
#include "ioutils.h"
#include "gMethods.h"

#include <list>
#include <set>

namespace geco{
    namespace bed{
        
        class gBedBaseEntry{
        public:
            typedef typename boost::tuple<std::string,gPos,gPos,std::string,unsigned int,char> i_type;
            std::string i_chrom;
            gPos i_start;
            gPos i_end;
            std::string i_name;
            unsigned int i_score;
            bool i_isForward;
            gBedBaseEntry(){
            }
            
            gBedBaseEntry(const std::string & chrom, gPos start, gPos end,const std::string & name=std::string(),gScore score=0,char strand='.'):i_chrom(chrom),i_start(start),i_end(end),i_name(name),i_score(score){
                if(strand=='-'){
                    i_isForward=false;
                }else if(strand=='.' || strand=='+'){
                    i_isForward=true;
                }
            }
            
            gBedBaseEntry(const i_type & tup){
                i_chrom=tup.get<0>();
                i_start=tup.get<1>();
                i_end=tup.get<2>();
                i_name=tup.get<3>();
                i_score=tup.get<4>();
                i_isForward= (char)tup.get<5>()!='-';
            }
            inline std::string getSequence(const geco::sequence_retrievers::gLocal2bitSequenceRetriever & retriever) const{
                return gElement(i_chrom,i_start,i_end,i_isForward,retriever,gElm).getSequence();
            }
            
        };

        std::ostream & operator << (std::ostream & out,const gBedBaseEntry e){
            out << e.i_chrom << ":" << e.i_start << "-" << e.i_end << " " << ((e.i_isForward)?("+"):("-")) << "\t" << e.i_score;
            return out;
        }

        class gNarrowPeak:public gBedBaseEntry{
        public:
            typedef typename boost::tuple<std::string,gPos,gPos,std::string,unsigned int,char,gScore,gScore,gScore,gPos> i_type;
            
            gScore i_signal;
            gScore i_pvalue;
            gScore i_qvalue;
            gPos i_peak;
            gNarrowPeak();
            
            gNarrowPeak(const std::string & chrom, gPos start, gPos end,const std::string & name, gScore score, char strand,gScore signal, gScore pvalue, gScore qvalue, gPos peak):gBedBaseEntry(chrom,start,end,name,score,strand),i_signal(signal),i_pvalue(pvalue),i_qvalue(qvalue),i_peak(peak){
            }

            gNarrowPeak(const i_type & tup){
                i_chrom=tup. template get<0>();
                i_start=tup.template get<1>();
                i_end=tup.template get<2>();
                i_name=tup.template get<3>();
                i_score=tup.template get<4>();
                i_isForward= tup.template get<5>()!='-';
                i_signal=tup.template get<6>();
                i_pvalue=tup.template get<7>();
                i_qvalue=tup.template get<8>();
                i_peak=tup.template get<9>();
            }
        };

        std::ostream & operator << (std::ostream & out,const gNarrowPeak  e){
            out << (const gBedBaseEntry &)e << "\t" << e.i_signal;
            return out;
        }

        template<typename entryType>
        class gBedData:public geco::utils::gCsvList<entryType>{
        public:
            gBedData(){
            }
            
            gBedData(const std::string & fileName):geco::utils::gCsvList<entryType>(fileName){
            }
            
            std::vector<std::string> getSequences(const std::string & fileName) const{
                geco::sequence_retrievers::gLocal2bitSequenceRetriever retriever(fileName,true);
                std::vector<std::string> res(this->size());
                //auto reg=this->begin();
                geco::methods::gMThreadException te;
#pragma omp parallel for
                for(gSize i=0;i<this->size();i++){
                    try{
                        auto reg=this->begin();
                        std::advance(reg,i);
                        const entryType & r=*reg;
                        res[i] = (std::string) gElement(r.i_chrom,gReferenceInterval(r.i_start,r.i_end),(bool) r.i_isForward,retriever,gElm).getSequence();
                    }catch(...){
                        te.CaptureException();
                    }
                }
                te.Rethrow();
                return res;
            }
        };

        struct comparePosition{
        public:
            inline bool operator () (const gBedBaseEntry & r1,const gBedBaseEntry & r2) const {
                if(r1.i_chrom==r2.i_chrom){
                    if ( r1.i_start == r2.i_start ) {
                        return r1.i_end < r2.i_end;
                    } else {
                        return r1.i_start < r2.i_start;
                    }
                }else{
                    return r1.i_chrom<r2.i_chrom;
                }
            }
        };

        struct comparePositionA{
        public:
            gSize i_maxOverlap;
            comparePositionA(gSize maxOverlap):i_maxOverlap(maxOverlap){
            }
            inline bool operator () (const gBedBaseEntry & r1,const gBedBaseEntry & r2) const {
                if(r1.i_chrom==r2.i_chrom){
                    return r1.i_end < r2.i_start;
                }else{
                    return r1.i_chrom<r2.i_chrom;
                }
            }
        };

        struct comparePositionUnique{
        public:
            inline bool operator () (const gBedBaseEntry & r1,const gBedBaseEntry & r2) const {
                if(r1.i_chrom==r2.i_chrom){
                    return r1.i_end<r2.i_start;
                }else{
                    return r1.i_chrom<r2.i_chrom;
                }
            }
        };

        template<typename entryType>
        std::vector<gBedData<entryType>> removeOverlappingEntries(const gBedData<entryType> & data){
            gBedData<entryType> good,bad;
            std::set<entryType,comparePositionA> tres(comparePositionA(10));
            for_each(data.begin(),data.end(),[&tres,&good,&bad](const entryType & entry){
                auto res=tres.insert(entry);
                if(!res.second){
                    bad.push_back(entry);
                }else{
                    good.push_back(entry);
                }
            });
            std::vector<gBedData<entryType>> ret({good,bad});
            return ret;
        }
        
        

    }
}
