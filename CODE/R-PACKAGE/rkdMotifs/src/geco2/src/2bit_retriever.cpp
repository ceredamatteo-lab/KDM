/***************************************************************************
 *   Copyright (C) 2009 by Uberto Pozzoli   *
 *   uberto@botia.bp.lnf.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "2bit_retriever.h"
#include<fstream>
#include<string>

using namespace std;
using namespace geco;
using namespace geco::sequence_retrievers;

//--------------------------------------------------------------------------
/// class gLocal2bitSequenceRetriever definition
//--------------------------------------------------------------------------

static const gChar * tetramers[256]={
 "TTTT", "TTTC", "TTTA", "TTTG", "TTCT", "TTCC", "TTCA", "TTCG", "TTAT", "TTAC", "TTAA", "TTAG", "TTGT", "TTGC", "TTGA", "TTGG",
 "TCTT", "TCTC", "TCTA", "TCTG", "TCCT", "TCCC", "TCCA", "TCCG", "TCAT", "TCAC", "TCAA", "TCAG", "TCGT", "TCGC", "TCGA", "TCGG",
 "TATT", "TATC", "TATA", "TATG", "TACT", "TACC", "TACA", "TACG", "TAAT", "TAAC", "TAAA", "TAAG", "TAGT", "TAGC", "TAGA", "TAGG",
 "TGTT", "TGTC", "TGTA", "TGTG", "TGCT", "TGCC", "TGCA", "TGCG", "TGAT", "TGAC", "TGAA", "TGAG", "TGGT", "TGGC", "TGGA", "TGGG",
 "CTTT", "CTTC", "CTTA", "CTTG", "CTCT", "CTCC", "CTCA", "CTCG", "CTAT", "CTAC", "CTAA", "CTAG", "CTGT", "CTGC", "CTGA", "CTGG",
 "CCTT", "CCTC", "CCTA", "CCTG", "CCCT", "CCCC", "CCCA", "CCCG", "CCAT", "CCAC", "CCAA", "CCAG", "CCGT", "CCGC", "CCGA", "CCGG",
 "CATT", "CATC", "CATA", "CATG", "CACT", "CACC", "CACA", "CACG", "CAAT", "CAAC", "CAAA", "CAAG", "CAGT", "CAGC", "CAGA", "CAGG",
 "CGTT", "CGTC", "CGTA", "CGTG", "CGCT", "CGCC", "CGCA", "CGCG", "CGAT", "CGAC", "CGAA", "CGAG", "CGGT", "CGGC", "CGGA", "CGGG",
 "ATTT", "ATTC", "ATTA", "ATTG", "ATCT", "ATCC", "ATCA", "ATCG", "ATAT", "ATAC", "ATAA", "ATAG", "ATGT", "ATGC", "ATGA", "ATGG",
 "ACTT", "ACTC", "ACTA", "ACTG", "ACCT", "ACCC", "ACCA", "ACCG", "ACAT", "ACAC", "ACAA", "ACAG", "ACGT", "ACGC", "ACGA", "ACGG",
 "AATT", "AATC", "AATA", "AATG", "AACT", "AACC", "AACA", "AACG", "AAAT", "AAAC", "AAAA", "AAAG", "AAGT", "AAGC", "AAGA", "AAGG",
 "AGTT", "AGTC", "AGTA", "AGTG", "AGCT", "AGCC", "AGCA", "AGCG", "AGAT", "AGAC", "AGAA", "AGAG", "AGGT", "AGGC", "AGGA", "AGGG",
 "GTTT", "GTTC", "GTTA", "GTTG", "GTCT", "GTCC", "GTCA", "GTCG", "GTAT", "GTAC", "GTAA", "GTAG", "GTGT", "GTGC", "GTGA", "GTGG",
 "GCTT", "GCTC", "GCTA", "GCTG", "GCCT", "GCCC", "GCCA", "GCCG", "GCAT", "GCAC", "GCAA", "GCAG", "GCGT", "GCGC", "GCGA", "GCGG",
 "GATT", "GATC", "GATA", "GATG", "GACT", "GACC", "GACA", "GACG", "GAAT", "GAAC", "GAAA", "GAAG", "GAGT", "GAGC", "GAGA", "GAGG",
 "GGTT", "GGTC", "GGTA", "GGTG", "GGCT", "GGCC", "GGCA", "GGCG", "GGAT", "GGAC", "GGAA", "GGAG", "GGGT", "GGGC", "GGGA", "GGGG"
};

typedef struct{
 unsigned int signature;     //- the number 0x1A412743 in the architecture of the machine that created the file
 unsigned int version;       //- zero for now. Readers should abort if they see a version number higher than 0.
 unsigned int sequenceCount; //- the number of sequences in the file.
 unsigned int reserved;      //- always zero for now
} em_2bit_header;

gLocal2bitSequenceRetriever::gLocal2bitSequenceRetriever(const gString & filename,bool useMask):gSequenceRetrieverImplementation(),i_filename(filename){
 //ifstream fin(((string) filename).c_str(),ios::binary);
 ifstream fin(((string) filename).c_str());
 em_2bit_header header;
 if(fin.is_open()){
  fin.read((char *) &header,16);
  if((header.signature==440477507)&&(header.version==0)){
   for(unsigned int i=0;i<header.sequenceCount;i++){
    char namelen;
    unsigned int offset;
    char *buffer;
    fin.read(&namelen,1);
    buffer=new char[(gSize)namelen+1];
    fin.read(buffer,(gSize)namelen);
    buffer[(gSize)namelen]=0;
    fin.read((char*)&offset,4);
    i_references.push_back(gString(buffer));
    i_offset.setValue((gPos) i,(gPos) offset,false);
    delete []buffer;
   }
   for(gPos i=0;i<i_offset.getSize();i++){
    unsigned int dnaSize,nBlockCount,maskBlockCount,reserved;
    unsigned int *nBlockStarts=NULL,*nBlockSizes=NULL,*maskBlockStarts=NULL,*maskBlockSizes=NULL;
    fin.seekg(i_offset[i],ios_base::beg);
    fin.read((char*)&dnaSize,4);
    i_dnaSize.setValue(i,dnaSize);
    fin.read((char*)&nBlockCount,4);
    if(nBlockCount>0){
     nBlockStarts=new unsigned int[nBlockCount];
     nBlockSizes=new unsigned int[nBlockCount];
     fin.read((char*) nBlockStarts,4*nBlockCount);
     fin.read((char*) nBlockSizes,4*nBlockCount);
     i_nBlockStarts.push_back(gArray<gPos>(0,nBlockCount,true));
     i_nBlockEnds.push_back(gArray<gPos>(0,nBlockCount,true));
     for(gPos j=0; j< nBlockCount;j++){
      i_nBlockStarts[i].setValue(j,(gPos) nBlockStarts[j],false);
      i_nBlockEnds[i].setValue(j,(gPos) nBlockStarts[j]+nBlockSizes[j],false);
     }
     delete [] nBlockStarts;
     delete [] nBlockSizes;
    }else{
     i_nBlockStarts.push_back(gArray<gPos>());
     i_nBlockEnds.push_back(gArray<gPos>());
    }

    fin.read((char*)&maskBlockCount,4);
    if(maskBlockCount>0){
     if(useMask){
      maskBlockStarts=new unsigned int[maskBlockCount];
      maskBlockSizes=new unsigned int[maskBlockCount];
      fin.read((char*) maskBlockStarts,4*maskBlockCount);
      fin.read((char*) maskBlockSizes,4*maskBlockCount);
      i_maskBlockStarts.push_back(gArray<gPos>(0,maskBlockCount,true));
      i_maskBlockEnds.push_back(gArray<gPos>(0,maskBlockCount,true));
      for(gPos j=0; j< maskBlockCount;j++){
       i_maskBlockStarts[i].setValue(j,(gPos) maskBlockStarts[j],false);
       i_maskBlockEnds[i].setValue(j,(gPos) maskBlockStarts[j]+maskBlockSizes[j],false);
      }
      delete [] maskBlockStarts;
      delete [] maskBlockSizes;
     }else{
      fin.seekg(8*maskBlockCount,ios_base::cur);
      i_maskBlockStarts.push_back(gArray<gPos>());
      i_maskBlockEnds.push_back(gArray<gPos>());
     }
    }else{
     i_maskBlockStarts.push_back(gArray<gPos>());
     i_maskBlockEnds.push_back(gArray<gPos>());
    }
    fin.read((char*)&reserved,4);
    if(reserved!=0) throw gException("gLocal2bitSequenceRetriever: Invalid 2bit file format");
    i_offset.setValue(i,i_offset[i]+16+(nBlockCount+maskBlockCount)*8,false);
   }
   fin.close();
  }else throw gException("gLocal2bitSequenceRetriever: Invalid sequence file");
 }else throw gException("gLocal2bitSequenceRetriever: Cannot open sequence file");
}

gLocal2bitSequenceRetriever::gLocal2bitSequenceRetriever(const gLocal2bitSequenceRetriever & retriever):gSequenceRetrieverImplementation(){
  i_filename=retriever.i_filename;
  i_references=retriever.i_references;
  i_dnaSize=retriever.i_dnaSize;
  i_offset=retriever.i_offset;
  i_nBlockStarts=retriever.i_nBlockStarts;
  i_nBlockEnds=retriever.i_nBlockEnds;
  i_maskBlockStarts=retriever.i_maskBlockStarts;
  i_maskBlockEnds=retriever.i_maskBlockEnds;
}

gLocal2bitSequenceRetriever::~gLocal2bitSequenceRetriever(){
}

gArrayRetrieverImplementation<gChar> * gLocal2bitSequenceRetriever::clone() const{
  return new gLocal2bitSequenceRetriever(*this);
}

gSequence gLocal2bitSequenceRetriever::getSequence_Internal(const gString & reference,gPos start,gPos end) const{
 size_t ref=0;
 while(ref<i_references.size()){
  if(i_references[ref]==reference){
   break;
  }else ref++;
 }
 if(ref<i_references.size()){
  if((start>=i_dnaSize[ref])||(end>i_dnaSize[ref])) throw gException("2bitlocalretreiver::Sequence positions out of boundaries");
  ifstream fin(((string) i_filename).c_str(),ios::binary);
  if(fin.is_open()){
   fin.seekg(i_offset[ref],ios_base::beg);
   size_t startbyte=start/4;
   size_t endbyte=end/4 + (size_t) (end % 4 !=0);
   unsigned char *bytes=new unsigned char[endbyte-startbyte];
   fin.seekg(startbyte,ios_base::cur);
   fin.read((char*) bytes,endbyte-startbyte);
   fin.close();   
   gString ret;
   gChar * Buffer=new char[4*(endbyte-startbyte)+1];
   for(size_t i=0;i<endbyte-startbyte;i++){
     memcpy(Buffer+4*i,tetramers[(size_t)(bytes[i])],4);
   }
   ret.ownData(4*(endbyte-startbyte),Buffer);
   delete [] bytes;
   gPos leading=(gPos)(start % 4);
   gPos trailing=ret.getLength()-(gPos)(3-((end-1) % 4));
   ret=gString(ret,leading,trailing);

   if(i_nBlockStarts[ref].getSize()>0){
    gArray<gPos> intv=which((i_nBlockStarts[ref]<end)&&(i_nBlockEnds[ref]>start));
    for(gPos i=0;i<intv.getSize();i++){
     gPos a=max(i_nBlockStarts[ref][intv[i]],start);
     gPos b=min(i_nBlockEnds[ref][intv[i]],end);
     ret.replace(a-start,gString(b-a,'N'));
    }
   }
   if(i_maskBlockStarts[ref].getSize()>0){
    gArray<gPos> intv=which((i_maskBlockStarts[ref]<end)&&(i_maskBlockEnds[ref]>start));
    for(gPos i=0;i<intv.getSize();i++){
     gPos a=max(i_maskBlockStarts[ref][intv[i]],start)-start;
     gPos b=min(i_maskBlockEnds[ref][intv[i]],end)-start;
     ret.replace(a,gString(ret,a,b).lower());
    }
   }
   return ret;
  }else throw gException("Invalid file");
 }else throw gException("Invalid sequence reference");
}

gSequence gLocal2bitSequenceRetriever::getRandomSequence(gSize length) const{
 gSequence ret;
 gArray<gPos> goodat=which(i_dnaSize>length);
 do{
  gPos choice=goodat[( rand() % goodat.getSize() )];
  gPos start=( rand() % (i_dnaSize[choice]-length));
  ret=getSequence_Internal(i_references[choice],start,start+length);
 }while(ret.getCounts('N').getSum()[0]!=0);
 return ret;
}

gSequence gLocal2bitSequenceRetriever::getRandomSequence(const gString & reference,gSize length) const{
 gSequence ret;
 vector<gString>::const_iterator id=find(i_references.begin(),i_references.end(),reference);
 if(id!=i_references.end()){
   gPos choice=id-i_references.begin();
   do{
     gPos start=( rand() % (i_dnaSize[choice]-length));
     ret=getSequence_Internal(i_references[choice],start,start+length);
   }while(ret.getCounts('N').getSum()[0]!=0);
 }else throw gException("gLocal2bitSequenceRetriever::getRandomSequence: invalid reference name");
 return ret;
}

gReferenceInterval gLocal2bitSequenceRetriever::getRandomInterval(gSize length) const{
  gArray<gPos> goodat=which(i_dnaSize>length);
  gPos choice=goodat[( rand() % goodat.getSize() )];
  gPos start=( rand() % (i_dnaSize[choice]-length));
  return gReferenceInterval(start,start+length);
}

gReferenceInterval gLocal2bitSequenceRetriever::getRandomInterval(const gString & reference,gSize length) const{
  gReferenceInterval ret;
  vector<gString>::const_iterator id=find(i_references.begin(),i_references.end(),reference);
  if(id!=i_references.end()){
    gPos choice=id-i_references.begin();
    gPos start=( rand() % (i_dnaSize[choice]-length));
    ret=gReferenceInterval(start,start+length);
  }else throw gException("gLocal2bitSequenceRetriever::getRandomInterval: invalid reference name");
  return ret;
}


gSize gLocal2bitSequenceRetriever::getReferenceLength(const gString & reference) const {
 gSize ret=0;
 size_t ref = 0;
 while(ref<i_references.size()){
  if(i_references[ref]==reference){
   break;
  }else ref++;
 }
 if(ref<i_references.size()){
  ret=i_dnaSize[ref];
 }
 return ret;
}

gSize gLocal2bitSequenceRetriever::getReferenceLength(gSize refIndex) const{
 gSize ret=0;
 if(refIndex<i_references.size()){
  ret=i_dnaSize[refIndex];
 }
 return ret;
}

gSize gLocal2bitSequenceRetriever::getReferenceCount() const{
  return i_references.size();
}

gString gLocal2bitSequenceRetriever::getReferenceName(gSize refIndex) const{
  gString ret;
  if(refIndex<i_references.size()){
    ret=i_references[refIndex];
  }else throw gException("gLocal2bitSequenceRetriever: invalid reference index");
  return ret;
}
//--------------------------------------------------------------------------


