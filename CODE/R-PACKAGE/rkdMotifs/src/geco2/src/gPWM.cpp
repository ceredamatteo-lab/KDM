#include "gPWM.h"
#include "TFM-Pvalue/Matrix.h"
#include <fstream>
#include <vector>
#include <map>
#include <sstream>


using namespace std;
using namespace geco;
using namespace geco::kmers;




static const unsigned char char_map2[256]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static const gPos char_map_NA2[248]= {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,66,68,69,70,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,98,100,101,102,104,105,106,107,108,109,110,111,112,113,114,115,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255};
const gArray<unsigned char> ast2 ( char_map2,256,gArray<gPos> ( char_map_NA2,248 ) );





//BEGIN CLASS gPWM
gPWM::gPWM()
{
}

gPWM::gPWM ( const gMatrix<gScore> & PWMcounts,gScore N, const gArray<gScore> & background,gScore pcounts ) :i_background ( background )
{
    setData ( PWMcounts,N,background,pcounts );
}


void gPWM::setData ( const gMatrix<gScore> & PWMcounts, gScore N, const gArray<gScore> & background,gScore pcounts )
{
    i_PWM=gMatrix<gScore> ( 4,PWMcounts.getColsNum(),PWMcounts/N );
    for ( gSize c=0; c<i_PWM.getColsNum(); c++ ) {
        if ( i_PWM.getCol ( c ).getSum() [0]!=1 ) {
            i_PWM.setCol ( c,i_PWM.getCol ( c ) /i_PWM.getCol ( c ).getSum() );
        }
    }
    i_PSSM=i_PWM;
    for ( gSize i=0; i<i_PSSM.getColsNum(); i++ ) {
        for ( gSize b=0; b<4; b++ ) {
            gScore bf=background[b];
            gScore nf= ( i_PWM ( b,i ) *N+bf*pcounts ) / ( N+pcounts );
            i_PSSM.setValue ( b,i,log2 ( nf/bf ),false );
        }
    }
}

gPWM::gPWM ( const gMatrix<gScore> & PWMdata, const gArray<gScore> & background,gScore pcounts,gScore tolerance ) :i_background ( background )
{
    setData ( PWMdata,background,pcounts,tolerance );
}

void gPWM::setData ( const gMatrix<gScore> & PWMdata, const gArray<gScore> & background,gScore pcounts,gScore tolerance )
{
    gScore N=1000;
    for ( gSize i=0; i<PWMdata.getColsNum(); i++ ) {
        if ( abs ( PWMdata.getCol ( i ).getSum() [0]-1 ) >tolerance ) {
            throw gException ( "gPWM: Invalid data matrix" );
        }
    }
    i_PWM=gMatrix<gScore> ( 4,PWMdata.getColsNum(),PWMdata );
    for ( gSize c=0; c<i_PWM.getColsNum(); c++ ) {
        if ( PWMdata.getCol ( c ).getSum() [0]!=1 ) {
            i_PWM.setCol ( c,i_PWM.getCol ( c ) /i_PWM.getCol ( c ).getSum() );
        } else {
            i_PWM.setCol ( c,i_PWM.getCol ( c ) );
        }
    }
    i_PSSM=i_PWM;
    for ( gSize i=0; i<i_PSSM.getColsNum(); i++ ) {
        for ( gSize b=0; b<4; b++ ) {
            gScore bf=background[b];
//       i_PSSM.setValue(b,i,log2( (i_PWM(b,i)+pcounts*bf)/bf),false);
            gScore pc=bf*pcounts;
            gScore nf= ( i_PWM ( b,i ) *N+pc ) / ( N+pcounts );
            i_PSSM.setValue ( b,i,log2 ( nf/bf ),false );
        }
    }
}

gArray<gScore> gPWM::getScores ( const gSequence seq,ScoringMode mode ) const
{
    gArray<gScore> ret ( 0,seq.getLength(),true );
    gSize rlen=seq.getLength()-i_PWM.getColsNum()+1;
    gArray<gSize> ncode=ast2[seq];
    gSize mlen=i_PWM.getColsNum();
    gArray<gSize> ind=getArray<gSize> ( 0,mlen-1,1 );
    if ( mlen<=seq.getLength() ) {
        if ( mode==PWMScore ) {
            for ( gSize i=0; i<rlen; i++ ) {
                gArray<gScore> temp=i_PWM ( gArray<gSize> ( ncode,i,i+mlen ),ind );
                gScore s=1;
                for ( gSize p=0; p<temp.getSize(); p++ ) {
                    s*=temp[p];
                }
                ret.setValue ( i,s,false );
            }
        } else if ( mode==PSSMScore ) {
            for ( gSize i=0; i<rlen; i++ ) {
                ret.setValue ( i,i_PSSM ( gArray<gSize> ( ncode,i,i+mlen ),ind ).getSum() [0],false );
            }
        }
    }
    return ret;
}

gArray<gScore> gPWM::getScores ( const gElement & element,bool doublestrand,ScoringMode mode ) const
{
    gSequence seq=element.getSequence();
    gSize rlen=seq.getLength()-i_PWM.getColsNum()+1;
    gArray<gScore> temp=getScores ( seq,mode );
    if ( doublestrand ) {
        gArray<gScore> rrtemp ( 0,temp.getSize(),false );
        gArray<gScore> rtemp=gArray<gScore> ( getScores ( seq.getReverseComplement(),mode ),0,rlen ).getReverted();
        for ( gSize i=0; i<rlen; i++ ) {
            rrtemp.setValue ( i,rtemp[i],false );
        }
        for ( gSize i=0; i<temp.getSize(); i++ ) {
            temp.setValue ( i,max ( temp[i],rrtemp[i] ),false );
        }
    }
    return temp;
}

pair< gArray<gScore>,gArray<bool> > gPWM::getDsScores ( const gElement & element,ScoringMode mode ) const
{
    gSequence seq=element.getSequence();
    gSize rlen=seq.getLength()-i_PWM.getColsNum()+1;
    gArray<gScore> temp=getScores ( seq,mode );
    gArray<gScore> rrtemp ( 0,temp.getSize(),false );
    gArray<bool> isForward ( true,temp.getSize(),false );
    gArray<gScore> rtemp=gArray<gScore> ( getScores ( seq.getReverseComplement(),mode ),0,rlen ).getReverted();
    for ( gSize i=0; i<rlen; i++ ) {
        rrtemp.setValue ( i,rtemp[i],false );
    }
    for ( gSize i=0; i<temp.getSize(); i++ ) {
        if ( temp[i]<rrtemp[i] ) {
            temp.setValue ( i,rrtemp[i],false );
            isForward.setValue ( i,false,false );
        }
    }
    return pair< gArray<gScore>,gArray<bool> > ( temp,isForward );
}

gArray<gScore> gPWM::getScores ( const gElement & element,const gElementRegion & region,bool doublestrand,ScoringMode mode ) const
{
    gSequence seq=region.getSequence ( element );
    gSize rlen=seq.getLength()-i_PWM.getColsNum()+1;
    gArray<gScore> temp=getScores ( seq,mode );
    if ( doublestrand ) {
        gArray<gScore> rrtemp ( 0,temp.getSize(),false );
        gArray<gScore> rtemp=gArray<gScore> ( getScores ( seq.getReverseComplement(),mode ),0,rlen ).getReverted();
        for ( gSize i=0; i<rlen; i++ ) {
            rrtemp.setValue ( i,rtemp[i],false );
        }
//     rrtemp.setValues(gArray<gPos>(0,rlen-1,1),rtemp,0,false);
//
//     for(gSize i=0;i<rtemp.getSize();i++){
//       cout << setw(10) << i+1 << setw(20) << temp[i] << setw(20) << rrtemp[i] << endl;
//     }
        for ( gSize i=0; i<temp.getSize(); i++ ) {
            temp.setValue ( i,max ( temp[i],rrtemp[i] ),false );
        }
    }
    return temp;
}

gScore gPWM::getMaxScore ( const gSequence seq,bool doublestrand,ScoringMode mode ) const
{
    gArray<gScore> temp=getScores ( seq,mode );
    if ( doublestrand ) {
        temp.concatenate ( getScores ( seq.getReverseComplement(),mode ) );
    }
    return temp.getMax() [0];
}

gArray<gScore> gPWM::getMaxScore ( const gElement & element,const gElementRegions & regions,bool doublestrand,ScoringMode mode ) const
{
    gArray<gScore> ret ( 0,regions.size(),false );
    for ( gSize r=0; r<regions.size(); r++ ) {
        ret.setValue ( r,getMaxScore ( regions.at ( r ).getSequence ( element ),doublestrand,mode ) );
    }
    return ret;
}

gArray<gScore> gPWM::getMaxScore ( const gElementCollection & elementCollection,const gElementRegions & regions,bool doublestrand,ScoringMode mode ) const
{
    gScore * tmp=new gScore[elementCollection.getElementCount() *regions.size()];
    #pragma omp parallel for
    for ( gSize e=0; e<elementCollection.getElementCount(); e++ ) {
        gArray<gScore> cas=getMaxScore ( elementCollection.getElement ( e ),regions,doublestrand,mode );
        memcpy ( tmp+e*regions.size(),cas.getConstData(),sizeof ( gScore ) *regions.size() );
    }
    gArray<gScore> ret;
    ret.ownData ( elementCollection.getElementCount() *regions.size(),tmp );
    return ret;
}

gScore gPWM::getAvgScore ( const gSequence seq,bool doublestrand,ScoringMode mode ) const
{
    gArray<gScore> temp=getScores ( seq,mode );
    if ( doublestrand ) {
        gArray<gScore> rtemp=getScores ( seq.getReverseComplement(),mode ).getReverted();
//     temp+=rtemp;
//     temp/=2;
        for ( gSize i=0; i<temp.getSize(); i++ ) {
            temp.setValue ( i,max ( temp[i],rtemp[i] ),false );
        }
    }
    return temp.getMean() [0];

}

gArray<gScore> gPWM::getAvgScore ( const gElement & element,const gElementRegions & regions,bool doublestrand,ScoringMode mode ) const
{
    gArray<gScore> ret ( 0,regions.size(),false );
    for ( gSize r=0; r<regions.size(); r++ ) {
        ret.setValue ( r,getAvgScore ( regions.at ( r ).getSequence ( element ),doublestrand,mode ) );
    }
    return ret;
}

gArray<gScore> gPWM::getAvgScore ( const gElementCollection & elementCollection,const gElementRegions & regions,bool doublestrand,ScoringMode mode ) const
{
    gScore * tmp=new gScore[elementCollection.getElementCount() *regions.size()];
    #pragma omp parallel for
    for ( gSize e=0; e<elementCollection.getElementCount(); e++ ) {
        gArray<gScore> cas=getAvgScore ( elementCollection.getElement ( e ),regions,doublestrand,mode );
        memcpy ( tmp+e*regions.size(),cas.getConstData(),sizeof ( gScore ) *regions.size() );
    }
    gArray<gScore> ret;
    ret.ownData ( elementCollection.getElementCount() *regions.size(),tmp );
    return ret;
}

pair<gScore,gScore> gPWM::getScoreRange() const
{
    pair<gScore,gScore> ret;
    ret.first=0;
    ret.second=0;
    for ( gPos k=0; k<=getLength(); k++ ) {
        ret.first+=i_PSSM.getCol ( k ).getMin() [0];
        ret.second+=i_PSSM.getCol ( k ).getMax() [0];
    }
    return ret;
}

gArray<gScore> gPWM::getPosIC() const
{

    gArray<gScore> ics ( 0,i_PWM.getColsNum(),true );
    for ( gPos i=0; i<i_PWM.getColsNum(); i++ ) {
        gScore ic=0;
        gArray<gScore> col=i_PWM.getCol ( i );
        for ( gSize b=0; b<3; b++ ) {
            if ( col[b]>0 ) {
                ic+=col[b]*log2 ( col[b] );
            }
        }
        ics.setValue ( i,ic+2,false );
    }
    return ics;
}

gScore gPWM::getIC() const
{
    return getPosIC().getMean() [0];
}

map<gScore,gScore> gPWM::computeDistribution ( const gArray<gScore> & range,gScore granularity ) const
{
    map<gScore,gScore> ret;

    //setup TFMpvalue matrix object
    Matrix M ( i_background[0],i_background[1],i_background[2],i_background[3] );

    M.length = i_PWM.getColsNum();
    M.mat = new double*[4];
    for ( int r = 0; r < 4; r++ ) {
        M.mat[r]=new double[M.length];
        for ( int c = 0; c < i_PWM.getColsNum(); c++ ) {
            M.mat[r][c]=i_PSSM ( r,c );
        }
    }

    M.computesIntegerMatrix ( granularity );

    gScore sMin= ( range.isNA ( 0 ) ) ? ( getScoreRange().first ) :range[0];
    gScore sMax= ( range.isNA ( 1 ) ) ? ( getScoreRange().second ) :range[1];
    long long int _min=sMin*M.granularity+M.offset;
    long long int _max=sMax*M.granularity+M.offset;

    map<long long, double> *nbocc = M.calcDistribWithMapMinMax ( _min,_max );
    map<long long, double>::iterator iter;
    double sum = 0;
    map<long long, double>::reverse_iterator riter = nbocc[M.length-1].rbegin();

    while ( riter != nbocc[M.length-1].rend() ) {
        sum += riter->second;
        nbocc[M.length][riter->first] = sum;
        riter++;
    }
    iter = nbocc[M.length].begin();
    while ( iter != nbocc[M.length].end() && iter->first <= _max ) {
        ret.insert ( pair<gScore,gScore> ( ( ( iter->first )-M.offset ) /M.granularity,iter->second ) );
        //cout << (((iter->first)-M.offset)/M.granularity) << " " << (iter->second) << " " << nbocc[M.length-1][iter->first] << endl;
        iter ++;
    }
    delete [] nbocc;
    return ret;
}

gArray<gScore> gPWM::computePvalues ( const gArray<gScore> & scores,gScore initialGranularity,bool forcedGranularity,gScore maxGranularity ) const
{
    gArray<gScore> pValues ( 0,scores.getSize(),false );
    //setup TFMpvalue matrix object
    Matrix M ( i_background[0],i_background[1],i_background[2],i_background[3] );
    M.length = i_PWM.getColsNum();
    M.mat = new double*[4];
    for ( int r = 0; r < 4; r++ ) {
        M.mat[r]=new double[M.length];
        for ( int c = 0; c < i_PWM.getColsNum(); c++ ) {
            M.mat[r][c]=i_PSSM ( r,c );
        }
    }

    for ( gSize s=0; s<scores.getSize(); s++ ) {
        double requestedScore=scores[s];
        long long max;
        long long min;
        double ppv;
        double pv;
        long long score;

        for ( double granularity = initialGranularity; granularity >= maxGranularity; granularity /= 10 ) {
            M.computesIntegerMatrix ( granularity );
            max = requestedScore*M.granularity + M.offset + M.errorMax+1;
            min = requestedScore*M.granularity + M.offset - M.errorMax-1;
            score = requestedScore*M.granularity + M.offset;
            M.lookForPvalue ( score,min,max,&ppv,&pv );
            if ( ppv == pv ) {
                if ( !forcedGranularity ) {
                    break;
                }
            }
        }
        pValues.setValue ( s,pv,false );
    }

    return pValues;
}

gScore gPWM::computeThreshold ( gScore alpha,gScore initialGranularity,bool forcedGranularity,gScore maxGranularity ) const
{

    Matrix M ( i_background[0],i_background[1],i_background[2],i_background[3] );
    M.length = i_PWM.getColsNum();
    M.mat = new double*[4];
    for ( int r = 0; r < 4; r++ ) {
        M.mat[r]=new double[M.length];
        for ( int c = 0; c < i_PWM.getColsNum(); c++ ) {
            M.mat[r][c]=i_PSSM ( r,c );
        }
    }
    M.computesIntegerMatrix ( initialGranularity );

    long long decrgr=10;
    long long max = M.maxScore+ceil ( M.errorMax+0.5 );
    long long min = M.minScore;
    double pv;
    long long score;

    for ( double granularity = initialGranularity; granularity >= maxGranularity; granularity /= decrgr ) {
        M.computesIntegerMatrix ( granularity );
        double ppv;
        score = M.lookForScore ( min,max,alpha,&pv,&ppv );

        min = ( score - ceil ( M.errorMax+0.5 ) ) * decrgr;
        max = ( score + ceil ( M.errorMax+0.5 ) ) * decrgr;

        if ( pv == ppv ) {
            if ( !forcedGranularity ) {
                break;
            }
        }

    }
    gScore ret= ( ( score-M.offset ) /M.granularity );

    return ret;
}

gSize gPWM::getLength() const
{
    return i_PWM.getColsNum();
}

gMatrix<gScore> gPWM::getPWM() const
{
    return i_PWM;
}

gMatrix<gScore> gPWM::getPSSM() const
{
    return i_PSSM;
}

gArray<gScore> gPWM::getBackground() const
{
    return i_background;
}

//END CLASS gPWM

//BEGIN CLASS gPWMScoreDistribution
gPWMScoreDistribution::gPWMScoreDistribution()
{
}

gPWMScoreDistribution::gPWMScoreDistribution ( const gPWM & pwm,gScore minAlpha,gScore granularity )
{
    gArray<gScore> range ( 0,2,true );
    range.setValue ( 0,pwm.computeThreshold ( minAlpha,granularity,false,granularity ),false );
    cerr << "\t" << range[0] << "\t";
    cerr.flush();

    map<gScore,gScore> ret;

    //setup TFMpvalue matrix object
    Matrix M ( pwm.getBackground() [0],pwm.getBackground() [1],pwm.getBackground() [2],pwm.getBackground() [3] );

    M.length = pwm.getLength();
    M.mat = new double*[4];
    for ( int r = 0; r < 4; r++ ) {
        M.mat[r]=new double[M.length];
        for ( int c = 0; c < pwm.getLength(); c++ ) {
            M.mat[r][c]=pwm.getPSSM() ( r,c );
        }
    }

    M.computesIntegerMatrix ( granularity );

    gScore sMin= ( range.isNA ( 0 ) ) ? ( pwm.getScoreRange().first ) :range[0];
    gScore sMax= ( range.isNA ( 1 ) ) ? ( pwm.getScoreRange().second ) :range[1];
    long long int _min=sMin*M.granularity+M.offset;
    long long int _max=sMax*M.granularity+M.offset;

    map<long long, double> *nbocc = M.calcDistribWithMapMinMax ( _min,_max );
    map<long long, double>::iterator iter;
    map<long long, double>::reverse_iterator riter = nbocc[M.length-1].rbegin();
    double sum = 0;

    while ( riter != nbocc[M.length-1].rend() ) {
        sum += riter->second;
        nbocc[M.length][riter->first] = sum;
        riter++;
    }
    iter = nbocc[M.length].begin();
    while ( iter != nbocc[M.length].end() && iter->first <= _max ) {
        insert ( pair<gScore,gScore> ( ( ( iter->first )-M.offset ) /M.granularity,iter->second ) );
        //cout << (((iter->first)-M.offset)/M.granularity) << " " << (iter->second) << " " << nbocc[M.length-1][iter->first] << endl;
        iter ++;
    }
    delete [] nbocc;
    insert ( pair<gScore,gScore> ( range[0],minAlpha ) );
}

gPWMScoreDistribution::gPWMScoreDistribution ( const gMatrix<gScore> & dmat )
{
    for ( gSize r=0; r<dmat.getRowsNum(); r++ ) {
        insert ( pair<gScore,gScore> ( dmat ( r,0 ),dmat ( r,1 ) ) );
    }
}

gMatrix<gScore> gPWMScoreDistribution::getDistributionAsMatrix() const
{
    gMatrix<gScore> M ( size(), ( gSize ) 2 );
    gSize count=0;
    for ( map<gScore,gScore>::const_iterator iter=begin(); iter!=end(); iter++ ) {
        M.setValue ( count,0,iter->first,false );
        M.setValue ( count++,1,iter->second,false );
    }
    return M;
}

gArray<gScore> gPWMScoreDistribution::pvalue ( gArray<gScore> score )
{
    gArray<gScore> ret ( 0,score.getSize(),true );
    gPWMScoreMap::iterator i1,i2;

    for ( gSize i=0; i<score.getSize(); i++ ) {
        i1=lower_bound ( score[i] );
        if ( i1==begin() ) {
            ret.setValue ( i,1,false );
        } else if ( i1==end() ) {
            ret.setValue ( i,1,false );
        } else {
            i2=i1;
            i2--;
            ret.setValue ( i, ( ( i2->second - i1->second ) * score[i]+ ( i1->second*i2->first - i2->second * i1->first ) ) / ( i2->first - i1->first ),false );
        }
    }
    return ret;
}

//END CLASS gPWMScoreDistribution


//BEGIN CLASS gPWMSet
typedef pair<gString,gPWM> namedPWM;

gPWMSet::gPWMSet()
{
}

void gPWMSet::addMatrix ( const gPWM & matrix,const gString & name )
{
    push_back ( namedPWM ( name,matrix ) );
    i_names_map.insert ( pair<gString,gSize> ( name,size()-1 ) );
}

const gPWM & gPWMSet::operator [] ( const gString & name ) const
{
    return at ( getIndex ( name ) ).second;
}

const gPWM & gPWMSet::operator [] ( gSize index ) const
{
    return at ( index ).second;
}

const gString & gPWMSet::getName ( gSize index ) const
{

    return at ( index ).first;
}

gPos gPWMSet::getIndex ( const gString & name ) const
{
    gSize index;
    try {
        index=i_names_map.at ( name );
    } catch ( out_of_range & r ) {
        throw gException ( r.what() );
    }
    return index;
}

gMatrix<gScore> gPWMSet::getScores ( const gSequence seq,gPWM::ScoringMode mode ) const
{
    gMatrix<gScore> ret ( size(),seq.getLength(),0,true );
    gSize pos=0;
    for ( const_iterator i=begin(); i!=end(); i++ ) {
        ret.setRow ( pos++,i->second.getScores ( seq,mode ) );
    }

//   for(gSize i=0;i<size();i++){
//     ret.setRow(i,at(i).second.getScores(seq,mode));
//   }

    return ret;
}

gArray<gScore> gPWMSet::getMaxScore ( const gSequence seq,bool doublestrand,gPWM::ScoringMode mode ) const
{
    gArray<gScore> ret ( size(),0,true );
    gSize pos=0;
    for ( const_iterator i=begin(); i!=end(); i++ ) {
        ret.setValue ( pos++,i->second.getMaxScore ( seq,doublestrand,mode ),false );
    }
//   for(gSize i=0;i<size();i++){
//     ret.setValue(i,at(i).second.getMaxScore(seq,doublestrand,mode),false);
//   }

    return ret;
}

gArray<gScore> gPWMSet::getAvgScore ( const gSequence seq,bool doublestrand,gPWM::ScoringMode mode ) const
{
    gArray<gScore> ret ( size(),0,true );
    gSize pos=0;
    for ( const_iterator i=begin(); i!=end(); i++ ) {
        ret.setValue ( pos++,i->second.getAvgScore ( seq,doublestrand,mode ),false );
    }
    return ret;
}




gMatrix<gScore> gPWMSet::getMaxScore(const vector<string> seqs,bool doublestrand,gPWM::ScoringMode mode) const{
    vector<gArray<gScore>> tret(seqs.size());
    gMatrix<gScore> ret(seqs.size(),size(),0,false);
#pragma omp parallel for
    for(gSize i=0;i<seqs.size();i++){
        tret[i]=getMaxScore(seqs[i],doublestrand,mode);
    }
    for(gSize i=0;i<seqs.size();i++){
        ret.setRow(i,tret[i]);
    }
    return ret;
}

gMatrix<gScore> gPWMSet::getAvgScore(const vector<string> seqs,bool doublestrand,gPWM::ScoringMode mode) const{
    vector<gArray<gScore>> tret(seqs.size());
    gMatrix<gScore> ret(seqs.size(),size(),0,false);
#pragma omp parallel for
    for(gSize i=0;i<seqs.size();i++){
        tret[i]=getAvgScore(seqs[i],doublestrand,mode);
    }
    for(gSize i=0;i<seqs.size();i++){
        ret.setRow(i,tret[i]);
    }
    return ret;
}




gMatrix<gScore> gPWMSet::getMaxScore(const vector<gSequence> seqs,bool doublestrand,gPWM::ScoringMode mode) const{
    vector<gArray<gScore>> tret(seqs.size());
    gMatrix<gScore> ret(seqs.size(),size(),0,false);
#pragma omp parallel for
    for(gSize i=0;i<seqs.size();i++){
        tret[i]=getMaxScore(seqs[i],doublestrand,mode);
    }
    for(gSize i=0;i<seqs.size();i++){
        ret.setRow(i,tret[i]);
    }
    return ret;
}

gMatrix<gScore> gPWMSet::getAvgScore(const vector<gSequence> seqs,bool doublestrand,gPWM::ScoringMode mode) const{
    vector<gArray<gScore>> tret(seqs.size());
    gMatrix<gScore> ret(seqs.size(),size(),0,false);
#pragma omp parallel for
    for(gSize i=0;i<seqs.size();i++){
        tret[i]=getAvgScore(seqs[i],doublestrand,mode);
    }
    for(gSize i=0;i<seqs.size();i++){
        ret.setRow(i,tret[i]);
    }
    return ret;
}


gMatrix<gScore> gPWMSet::getMaxFeatures ( const gElementCollection & elementCollection,const gElementRegions & regions,bool doublestrand,gPWM::ScoringMode mode ) const
{

    gScore * results=new gScore[elementCollection.getElementCount() *size()];
    #pragma omp parallel for
    for ( gSize i=0; i<elementCollection.getElementCount(); i++ ) {
        gArray<gScore> a=getMaxScore ( regions.at ( 0 ).getSequence ( elementCollection.getElement ( i ) ),doublestrand,mode );
        memcpy ( results+i*size(),a.getConstData(),size() *sizeof ( gScore ) );
    }
    gArray<gScore> b;
    b.ownData ( elementCollection.getElementCount() *size(),results );
    return gMatrix<gScore> ( elementCollection.getElementCount(),size(),b );
}

gMatrix<gScore> gPWMSet::getAvgFeatures ( const gElementCollection & elementCollection,const gElementRegions & regions,bool doublestrand,gPWM::ScoringMode mode ) const
{
    gMatrix<gScore> ret ( elementCollection.getElementCount(),size(),0,true );
    for ( gSize i=0; i<elementCollection.getElementCount(); i++ ) {
        ret.setRow ( i,getAvgScore ( regions.at ( 0 ).getSequence ( elementCollection.getElement ( i ) ),doublestrand,mode ) );
    }
    return ret;
}
//END CLASS gPWMSet


//BEGIN CLASS gPWMSetScoreDistributions
gPWMSetScoreDistributions::gPWMSetScoreDistributions ( const gPWMSet & set,gScore minAlpha,gScore granularity ) :vector<gPWMScoreDistribution> ( set.size() )
{
    #pragma omp parallel for
    for ( gSize i=0; i<set.size(); i++ ) {
        cerr << setw ( 25 ) << set.getName ( i ) << "\t" << set[i].getScoreRange().first << " - " <<  set[i].getScoreRange().second;
        at ( i ) =gPWMScoreDistribution ( set[i],minAlpha,granularity );
        cerr << "\t(" << at ( i ).size() << ")"  << endl;
    }
}
//END CLASS gPWMSetScoreDistributions


gPWMSet geco::kmers::readFromMemeFile ( const gString & memeFileName, const gArray<gScore> & background,gScore pcounts,gScore tolerance )
{
    gPWMSet res;
    bool readBackground=true;
    ifstream fin ( ( ( string ) memeFileName ).c_str() );
    if ( !fin.is_open() ) {
        throw gException ( "readFromMemeFile: cannot open meme file" );
    }
    string line;
    gArray<gScore> actbackground ( 0.25,4,false );
    if ( background.getSize() >0 ) {
        actbackground=background;
        readBackground=false;
    }

    gSize nmotif=0;
    vector< gMatrix<gScore> > PWMData;

    while ( !fin.eof() ) {
        getline ( fin,line );
        if ( fin.good() ) {
            if ( line=="Background letter frequencies (from uniform background):" && readBackground ) {
                getline ( fin,line );
                if ( !fin.good() ) {
                    throw gException ( "readFromMemeFile: bad format, missing background" );
                }
                gString ss;
                vector<gString> parts;
                ss=line;
                ss.split ( ' ',parts );
                actbackground.setValue ( 0,atof ( ( ( string ) parts[1] ).c_str() ),false );
                actbackground.setValue ( 1,atof ( ( ( string ) parts[3] ).c_str() ),false );
                actbackground.setValue ( 2,atof ( ( ( string ) parts[5] ).c_str() ),false );
                actbackground.setValue ( 3,atof ( ( ( string ) parts[7] ).c_str() ),false );
            } else if ( line.substr ( 0,5 ) =="MOTIF" ) {
                gString name=line.substr ( 6,1000 );
                while ( !fin.eof() ) {
                    getline ( fin,line );
                    gString ss1=line;
                    vector<gString> parts1;
                    ss1.split ( ' ',parts1 );
                    if ( ( parts1.size() >0 ) && ( parts1[0]=="letter-probability" ) ) {
                        gSize alength=0,w=0,nsites=1000000;
                        gDoubleScore E=0;

                        for ( gSize p=1; p<parts1.size(); p++ ) {
                            if ( parts1[p]=="alength=" ) {
                                alength=atoi ( ( ( string ) parts1[p+1] ).c_str() );
                            }
                            if ( parts1[p]=="w=" ) {
                                w=atoi ( ( ( string ) parts1[p+1] ).c_str() );
                            }
                            if ( parts1[p]=="nsites=" ) {
                                nsites=atoi ( ( ( string ) parts1[p+1] ).c_str() );
                            }
                            if ( parts1[p]=="E=" ) {
                                E=atof ( ( ( string ) parts1[p+1] ).c_str() );
                            }
                        }
                        if ( alength==4 ) {
                            gMatrix<gScore> temp ( 4,w,0,false );
                            for ( gSize i=0; i<w; i++ ) {
                                gScore A,C,G,T;
                                fin >> A >> C >> G >> T;
                                if ( !fin.good() ) {
                                    throw gException ( "readFromMemeFile: bad format, missing matrix row(s)" );
                                }
                                temp.setValue ( 0,i,A*nsites,false );
                                temp.setValue ( 1,i,C*nsites,false );
                                temp.setValue ( 2,i,G*nsites,false );
                                temp.setValue ( 3,i,T*nsites,false );
                                //temp.setCol(i,temp.getCol(i)/temp.getCol(i).getSum());
                            }

                            res.addMatrix ( gPWM ( temp,nsites,actbackground,pcounts ),name );
                            break;
                        } else {
                            throw gException ( "readFromMemeFile: bad alphabet length" );
                        }
                    } else if ( parts1.size() >0 ) {
                        throw gException ( "readFromMemeFile: bad format, missing motif" );
                    }
                }
            }
        }
    }
    fin.close();
    return res;
}

/*
gPWMSet geco::motifs::readFromHomerFile(const gString & memeFileName, const gArray<gScore> & background){
  gPWMSet res;

  ifstream fin(((string)memeFileName).c_str());
  if(!fin.is_open()) throw gException("readFromHomerFile: cannot open homer file");
  string line;
  gArray<gScore> actbackground=background();

  gSize nmotif=0;
  vector< gMatrix<gScore> > PWMData;

  while(!fin.eof()){
    while(line[0]!='>' && fin.good()){
      getline(fin,line);
    }
    if(fin.good()){
      gString ss;
      vector<gString> parts;
      ss.split('\t',parts);
      name
      strsplit
      if(line=="Background letter frequencies (from uniform background):" && readBackground){
	getline(fin,line);
	if(!fin.good()) throw gException("readFromMemeFile: bad format, missing background");
	gString ss;
	vector<gString> parts;
	ss=line;
	ss.split(' ',parts);
	actbackground.setValue(0,atof(((string)parts[1]).c_str()),false);
	actbackground.setValue(1,atof(((string)parts[3]).c_str()),false);
	actbackground.setValue(2,atof(((string)parts[5]).c_str()),false);
	actbackground.setValue(3,atof(((string)parts[7]).c_str()),false);
      }else if(line.substr(0,5)=="MOTIF"){
	gString name=line.substr(6,1000);
	while(!fin.eof()){
	  getline(fin,line);
	  gString ss1=line;
	  vector<gString> parts1;
	  ss1.split(' ',parts1);
	  if((parts1.size()>0) && (parts1[0]=="letter-probability")){
	    gSize alength=0,w=0,nsites=0;
	    gDoubleScore E=0;

	    for(gSize p=1;p<parts1.size();p++){
	      if(parts1[p]=="alength=") alength=atoi(((string)parts1[p+1]).c_str());
	      if(parts1[p]=="w=") w=atoi(((string)parts1[p+1]).c_str());
	      if(parts1[p]=="nsites=") nsites=atoi(((string)parts1[p+1]).c_str());
	      if(parts1[p]=="E=") E=atof(((string)parts1[p+1]).c_str());
	    }
	    if(alength==4){
	      gMatrix<gScore> temp(4,w,0,false);
	      for(gSize i=0;i<w;i++){
		gScore A,C,G,T;
		fin >> A >> C >> G >> T;
		if(!fin.good()) throw gException("readFromMemeFile: bad format, missing matrix row(s)");
		temp.setValue(0,i,A,false);
		temp.setValue(1,i,C,false);
		temp.setValue(2,i,G,false);
		temp.setValue(3,i,T,false);
	      }
	      //res.push_back(namedPWM(name,gPWM(temp,actbackground)));
	      res.insert(pair<string,gPWM>(name,gPWM(temp,actbackground)));
	      break;
	    }else throw gException("readFromMemeFile: bad alphabet length");
	  }else if(parts1.size()>0) throw gException("readFromMemeFile: bad format, missing motif");
	}
      }
    }
  }
  fin.close();
  return res;
}
*/
