/***************************************************************************
 *   Copyright (C) 2010 by Uberto Pozzoli and Matteo Cereda                *
 *   uberto.pozzoli@bp.lnf.it                                              *
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
#include "transcript.h"

using namespace std;
using namespace geco;

//--------------------------------------------------------------------------
// gElement implementation
//--------------------------------------------------------------------------
bool gTranscript::initInternal(const gArray<gPos> & refExonStarts,const gArray<gPos> & refExonEnds,bool forward,gPos cdsStart,gPos cdsEnd){
    i_exons=insertReferenceConnections(refExonStarts,refExonEnds,forward);
    gSize nexons=i_exons.getSize();
    gArray<gArrayIndex> e_startidx=getConnectionsStartSiteIndex(i_exons);
    gArray<gArrayIndex> e_endidx=getConnectionsEndSiteIndex(i_exons);
    i_transcript=connectSites(e_startidx[0],e_endidx[nexons-1]);
    i_sitetype.setValue(e_startidx[0],tss,false);
    i_sitetype.setValue(e_endidx[nexons-1],pas,false);
    if(nexons>1){
      gArray<gArrayIndex> i_startidx=e_endidx[getArray<gArrayIndex>(0,nexons-2,1)];
      gArray<gArrayIndex> i_endidx=e_startidx[getArray<gArrayIndex>(1,nexons-1,1)];
      i_sitetype.setValues(i_startidx,acceptor);
      i_sitetype.setValues(i_endidx,donor);
      i_introns=connectSites(i_startidx,i_endidx);
    }
    if (cdsEnd>cdsStart) {
      i_cds=insertReferenceConnections(cdsStart,cdsEnd,forward);
    } else {
      i_cds=gArray<gArrayIndex>(0,1,true,true);
    }
    return true;
}

/** @brief Default constructor
 *
 * Instatiates an empy gTranscript object
 */
gTranscript::gTranscript() {
}

/** @brief Constructor from reference
    *
    * Instatiates a gTranscript initializing it with the provided reference informations
    * @param refName name of the reference
    * @param refExonStarts array containing the starting reference position of exons (0 based)
    * @param refExonEnds array containing the ending reference position of exons (1 based)
    * @param forward strand of the transcript referred to the reference
    * @param startOffset number of BPs to be incuded at 5' end (non negative, defaults to 0)
    * @param endOffset  number of BPs to be incuded at 3' end (non negative, defaults to 0)
    * @param cdsStart the reference position of the cds start (0 based, defaults to 0).
    * @param cdsEnd the reference position of the cds end (1 based, defaults to 0, meaning no cds information).
    */
gTranscript::gTranscript(const gString & refName,const gArray<gPos> & refExonStarts,const gArray<gPos> & refExonEnds,bool forward,gPos startOffset,gPos endOffset,gPos cdsStart,gPos cdsEnd):gElement(refName,refExonStarts.getMin()[0]+((forward)?(-startOffset):(-endOffset)),refExonEnds.getMax()[0]+((forward)?(endOffset):(startOffset)),forward) {
 initInternal(refExonStarts,refExonEnds,forward,cdsStart,cdsEnd);
}

/** @brief Constructor from reference with sequence retriever
 *
 * Instatiates a gTranscript initializing it with the provided reference informations
 * A Sequence retriever must be provided.
 * @param refName name of the reference
 * @param refExonStarts array containing the starting reference position of exons (0 based)
 * @param refExonEnds array containing the ending reference position of exons (1 based)
 * @param forward strand of the transcript referred to the reference
 * @param retriever the sequence retriever object
 * @param mode a gElementSequenceMode value to spcify the way the sequence must be allocated
 * @param startOffset number of BPs to be incuded at 5' end (non negative, defaults to 0)
 * @param endOffset  number of BPs to be incuded at 3' end (non negative, defaults to 0)
 * @param cdsStart the reference position of the cds start (0 based, defaults to 0).
 * @param cdsEnd the reference position of the cds end (1 based, defaults to 0, meaning no cds information).
 */
gTranscript::gTranscript(const gString & refName,const gArray<gPos> & refExonStarts,const gArray<gPos> & refExonEnds,bool forward,const gSequenceRetriever & retriever,gElementSequenceMode mode,gPos startOffset,gPos endOffset,gPos cdsStart,gPos cdsEnd): gElement(refName,refExonStarts.getMin()[0]+((forward)?(-startOffset):(-endOffset)),refExonEnds.getMax()[0]+((forward)?(endOffset):(startOffset)),forward,retriever,mode) {
 initInternal(refExonStarts,refExonEnds,forward,cdsStart,cdsEnd);
}

/** @brief Constructor from reference with sequence
 *
 * Instatiates a gTranscript initializing it with the provided reference informations
 * A Sequence must be provided.
 * @param refName name of the reference
 * @param refExonStarts array containing the starting reference position of exons (0 based)
 * @param refExonEnds array containing the ending reference position of exons (1 based)
 * @param forward strand of the transcript referred to the reference
 * @param sequence the sequence
 * @param mode a gElementSequenceMode value to spcify the way the sequence must be allocated
 * @param startOffset number of BPs to be incuded at 5' end (non negative, defaults to 0)
 * @param endOffset  number of BPs to be incuded at 3' end (non negative, defaults to 0)
 * @param cdsStart the reference position of the cds start (0 based, defaults to 0).
 * @param cdsEnd the reference position of the cds end (1 based, defaults to 0, meaning no cds information).
 */
gTranscript::gTranscript(const gString & refName,const gArray<gPos> & refExonStarts,const gArray<gPos> & refExonEnds,bool forward,const gSequence & sequence,gElementSequenceMode mode,gPos startOffset,gPos endOffset,gPos cdsStart,gPos cdsEnd):gElement(refName,refExonStarts.getMin()[0]+((forward)?(-startOffset):(-endOffset)),refExonEnds.getMax()[0]+((forward)?(endOffset):(startOffset)),forward,sequence,mode) {
 initInternal(refExonStarts,refExonEnds,forward,cdsStart,cdsEnd);
}

/** @brief Copy constructor
 *
 * Instatiate a copy of the passed object optionally adding variations.
 * When variations are present all the transcript sites positions are updated.
 * if some of these positions are deleted their resulting position will be marked as NA.
 * @param transcript the obejct to copy from
 * @param variations the variation to add (defaults to an empty set);
 */
gTranscript::gTranscript(const gTranscript & transcript,const gVariations & variations):gElement(transcript,0,0,true,variations) {
    i_cds=transcript.i_cds;
    i_transcript=transcript.i_transcript;
    i_exons=transcript.i_exons;
    i_introns=transcript.i_introns;
    i_sitetype=transcript.i_sitetype;
}

/** @brief assignment operator
 *
 * Makes this object equal to the passed one
 * @param transcript the object to copy from
 * @return a reference to this object
 */
gTranscript & gTranscript::operator = (const gTranscript & transcript) {
    *(gElement *)this=transcript;
    i_cds=transcript.i_cds;
    i_transcript=transcript.i_transcript;
    i_exons=transcript.i_exons;
    i_introns=transcript.i_introns;
    i_sitetype=transcript.i_sitetype;
    return *this;
}

/** @brief Exon count
 *
 * Returns the number of exons in this transcript
 * @return number of exons
 */
gSize gTranscript::getExonCount() const {
    return i_exons.getSize();
}

/** @brief Intron count
 *
 * Return the number of introjns in this transcript
 * @return a gSize
 */
gSize gTranscript::getIntronCount() const{
    return i_exons.getSize()-1;
}

/** @brief pre-mRNA length
 *
 * Returns the length of the pre-mRNA corresponding to this transcript
 * @return gSize value
 */
gSize gTranscript::getPremRNALength() const{
  return getPremRNAInterval().getLength();
}

/** @brief Exon length
 *
 * Returns the length of the exon referred by exonNum.
 * If one (or both) boundaries have been deleted ('cause variations have been applied)
 * the resulting length will be returned.
 * @param exonNum exon number
 * @return gSize value
 */
gSize gTranscript::getExonLength(gArrayIndex exonNum) const {
    return getExonInterval(exonNum).getLength();
}

/** @brief Intron length
 *
 * Returns the length of the intron referred by intronNum
 * If one (or both) boundaries have been deleted ('cause variations have been applied)
 * the resulting length will be returned.
 * @param intronNum intron number
 * @return the intron length
 */
gSize gTranscript::getIntronLength(gArrayIndex intronNum) const {
    return getIntronInterval(intronNum).getLength();
}

/** @brief pre-mRNA interval
 *
 * Returns the element boundaries of the pre-mRNA of this transcript.
 * This is not trivial since gTranscripts can be constructed with boundaries offsets and
 * variations can modify or invalidate these positions.
 * If one (or both) boundaries have been deleted (variations have been applied)
 * the start and/or end positons will be marked as NA and will match to the last valid
 * position before the deletion and to the first valid position after the deletion respectively.
 * @return a two values gArray
 */
gElementInterval gTranscript::getPremRNAInterval() const{
  return getConnectionElementInterval(i_transcript[0]);
}

/** @brief Exon interval
 *
 * Returns the boundaries of one exon in element coordinates.
 * If one (or both) boundaries have been deleted (variations have been applied)
 * the start and/or end positons will be marked as NA and will match to the last valid
 * position before the deletion and to the first valid position after the deletion respectively.
 * @param exonNum exon number
 * @return a two valued gArray with exon boundaries information
 */
gElementInterval gTranscript::getExonInterval(gArrayIndex exonNum) const {
    return getConnectionElementInterval(i_exons[exonNum]);
}

/** @brief Intron interval
 *
 * Returns the boundaries of one intron in element coordinates.
 * If one (or both) boundaries have been deleted (variations have been applied)
 * the start and/or end positons will be marked as NA and will match to the last valid
 * position before the deletion and to the first valid position after the deletion respectively.
 * @param intronNum intron number
 * @return a two valued gArray with intron boundaries information
 */
gElementInterval gTranscript::getIntronInterval(gArrayIndex intronNum) const {
  gElementInterval tmp=getConnectionElementInterval(i_introns[intronNum]);
  return gElementInterval(tmp.getStart()+1,tmp.getEnd()-1,tmp.validStart(),tmp.validEnd());
}

/** @brief 5' UTR interval
 *
 * Returns the genomic boundaries  in element coordiantes of the 5' UTR.
 * If one (or both) boundaries have been deleted ('cause variations have been applied)
 * the corresponding positon(s) will be marked as NA and corresponding to the first
 * valid position after the deletion and to the last valid position before the deletion
 * for start and end respectively.
 * @return a two valued gArray with cds boundaries information
 */
gElementInterval gTranscript::get5UTRInterval() const{
  gElementInterval intv=getCdsInterval();
  gElementInterval mrna=getPremRNAInterval();
  return gElementInterval(mrna.getStart(),intv.getStart(),mrna.validStart(),intv.validStart());
}

/** @brief Coding sequence interval
 *
 * Returns the genomic boundaries of the coding sequence.
 * If one (or both) boundaries have been deleted ('cause variations have been applied)
 * the corresponding positon(s) will be marked as NA and corresponding to the first
 * valid position after the deletion and to the last valid position before the deletion
 * for start and end respectively.
 * @return a two valued gArray with cds boundaries information
 */
gElementInterval gTranscript::getCdsInterval() const{
 gElementInterval ret;
 if(!i_cds.isNA(0)){
   ret=getConnectionElementInterval(i_cds[0]);
 }
 return ret;
}

/** @brief 3' UTR interval
 *
 * Returns the genomic boundaries  in element coordiantes of the 3' UTR.
 * If one (or both) boundaries have been deleted ('cause variations have been applied)
 * the corresponding positon(s) will be marked as NA and corresponding to the first
 * valid position after the deletion and to the last valid position before the deletion
 * for start and end respectively.
 * @return a two valued gArray with cds boundaries information
 */
gElementInterval gTranscript::get3UTRInterval() const{
  gElementInterval intv=getCdsInterval();
  gElementInterval mrna=getPremRNAInterval();
  return gElementInterval(intv.getEnd(),mrna.getEnd(),intv.validStart(),mrna.validStart());

}

/** @brief mRNA length
 *
 * Returns the length of the mRNA associated to this transcript
 * considered as the sum of all exon lengths.
 * For this reason the resulting length will be markled as NA when
 * one or more exon boundaries have been deleted.
 * @return associated mRNA length information
 */
gArray<gSize> gTranscript::getmRNALength() const{
  gArray<gSize> ret(0,1,false);
  for(gSize i=0;i<i_exons.getSize();i++) ret+=getExonLength(i);
  return ret;
}

/** @brief 5' UTR length
 *
 * Returns the length of the 5' UTR portion of the mRNA associated to this transcript
 * @return associated mRNA length information
 */
gArray<gSize> gTranscript::get5UTRLength() const{
  gArray<gSize> ret(0,1,false);
  gArray<gRelativePos> mapped=mapElementPositionsOnmRNA(getCdsInterval().getStart());
  if(mapped.NACount()==0){
    ret=mapElementPositionsOnmRNA(getCdsInterval().getStart())[0];
  }
  return ret;
}

/** @brief Coding sequence length
 *
 * Returns the length of the coding portion of the mRNA associated to this transcript
 * @return a gArray<gSize>
 */
gArray<gSize> gTranscript::getCdsLength() const{
  gArray<gSize> len(0,1,true);
  gElementInterval intv=getCdsInterval();
  if(intv.getLength()>0){
    gArray<gPos> ps=mapElementPositionsOnmRNA(intv.getPositions());
    if(ps.NACount()==0){
      len=gArray<gSize>(ps[1]-ps[0]+1,1,false);
    }
  }
  return len;
}

/** @brief 3' UTR length
 *
 * Returns the length of the 3' UTR portion of the mRNA associated to this transcript
 * @return associated mRNA length information
 */
gArray<gSize> gTranscript::get3UTRLength() const{
  gArray<gSize> ret(0,1,true);
  gArray<gRelativePos> mapped=mapElementPositionsOnmRNA(getCdsInterval().getEnd()-1);
  gArray<gSize> mrnal=getmRNALength();
  if(!mapped.isNA(0) && (mrnal.NACount()==0)){
    ret=mrnal-mapped[0]-1;
  }
  return ret;
}

/** @brief Exon phases
 *
 * Returns a the starting and ending phases of the specidfied exon
 * @param exonNum the exon number (0 based)
 * @return a gArray<gShortUnsigned> object
 */
gArray<gShortUnsigned> gTranscript::getExonPhase(gArrayIndex exonNum) const{
  gArray<gShortUnsigned> ret(0,2,true);
  gArray<gPos> mm=mapElementPositionsOnCds(getExonInterval(exonNum).getPositions());
  ret.setValue(0,mm[0] % 3,mm.isNA(0));
  ret.setValue(1,(mm[1]+1) % 3,mm.isNA(1));
  return ret;
}

/** @brief mRNA sequence retriever
 *
 * Returns a gSequence containing the mRNA
 * @param start the starting postion of the portion of mRNA sequence to retrieve (defaults 0)
 * @param end the ending postion of the portion of mRNA sequence to retrieve (defaults 0, interpreted as the length of the mRNA)
 * @return a gSequence object
 */
gSequence gTranscript::getmRNASequence(gPos start,gPos end) const{
  gSequence ret;
  for(gArrayIndex i=0;i<i_exons.getSize();i++){
    gElementInterval intv=getExonInterval(i);
    if(intv.validStart() && intv.validEnd()){
      ret+=getSequence(getExonInterval(i));
    }else{
      break;
    }
  }
  gPos aend=(end==0)?(ret.getLength()):(end);
  return gSequence(ret,start,aend);
}

/** @brief mRNA 5' UTR sequence retriever
 *
 * Returns a gSequence containing the mRNA 5' UTR
 * @param start the starting postion of the portion of mRNA sequence to retrieve (defaults 0)
 * @param end the ending postion of the portion of mRNA sequence to retrieve (defaults 0, interpreted as the length of the 5' UTR)
 * @return a gSequence object
 */
gSequence gTranscript::get5UTRSequence(gPos start,gPos end) const{
  gSequence mRNA=getmRNASequence();
  gArray<gSize> utr5Len=get5UTRLength();
  gPos aend=start;
  if(!utr5Len.isNA(0)){
    aend=(end==0)?(utr5Len[0]):(end);
  }
  return gSequence(mRNA,start,aend);
}

/** @brief mRNA Coding sequence retriever
 *
 * Returns a gSequence containing the mRNA coding sequence
 * @param start the starting postion of the portion of mRNA sequence to retrieve (defaults 0)
 * @param end the ending postion of the portion of mRNA sequence to retrieve (defaults 0, interpreted as the length of the coding sequence)
 * @return a gSequence object
 */
gSequence gTranscript::getCdsSequence(gPos start,gPos end) const{
  gSequence ret;
  gSequence mRNA=getmRNASequence();
  gArray<gSize> utr5Len=get5UTRLength();
  gArray<gSize> cdsLen=getCdsLength();
  gArray<gPos> astart=utr5Len+start;
  gArray<gPos> aend=(end==0)?(cdsLen+utr5Len):(utr5Len+end);
  if(!astart.isNA(0)){
    if(!aend.isNA(0)){
      ret=gSequence(mRNA,astart[0],aend[0]);
    }
  }
  return ret;
}

/** @brief mRNA 3' UTR retriever
 *
 * Returns a gSequence containing the mRNA 3' UTR
 * @param start the starting postion of the portion of mRNA sequence to retrieve (defaults 0)
 * @param end the ending postion of the portion of mRNA sequence to retrieve (defaults 0, interpreted as the length of the 3' UTR)
 * @return a gSequence object
 */
gSequence gTranscript::get3UTRSequence(gPos start,gPos end) const{
  gSequence ret;
  gSequence mRNA=getmRNASequence();
  gArray<gSize> utr5len=get5UTRLength();
  gArray<gSize> utr3len=get3UTRLength();
  gArray<gSize> cdslen=getCdsLength();
  if(!utr5len.isNA(0) && !cdslen.isNA(0) && !utr5len.isNA(0)) {
    gPos aend=(end==0)?(utr3len[0]):(end);
    ret=gSequence(mRNA,utr5len[0]+cdslen[0]+start,utr5len[0]+cdslen[0]+aend);
  }
  return ret;
}

/** @brief Mapping of element positions to the mRNA
 *
 * Maps the positions to the mRAN returning a gArray containing positions from the beginning of the mRNA.
 * If a position is outside the mRNA it returns a NA.
 * @param elementPositions array of element positions
 * @return a garray<gPos> object
 */
gArray<gRelativePos> gTranscript::mapElementPositionsOnmRNA(const gArray<gRelativePos> & elementPositions) const{
 gArray<gRelativePos> ret(0,elementPositions.getSize(),true);
 gPos currentpos=0;
  for(gArrayIndex i=0;i<getExonCount();i++){
    gElementInterval ei=getExonInterval(i);
    if(ei.validStart() && ei.validEnd()){
      gArray<gPos> good=which(ei.contains(elementPositions));
      if(good.getSize()>0){
	ret.setValues(good,elementPositions[good]-ei.getStart()+currentpos);
      }
      currentpos+=getExonLength(i);
    }else{
      break;
    }
  }
 return ret;
}

/** @brief Mapping of element positions to the mRNAm 5'UTR
 *
 * Maps the positions to the mRNA 5'UTR returning a gArray containing positions from the beginning of the mRNA 5'UTR.
 * If a position is outside the mRNA 5'UTR it returns a NA.
 * @param elementPositions array of element positions
 * @return a garray<gPos> object
 */
gArray<gRelativePos> gTranscript::mapElementPositionsOn5UTR(const gArray<gRelativePos> & elementPositions) const{
  gArray<gPos> sd=mapElementPositionsOnmRNA(getCdsInterval().getPositions());
  gArray<gPos> tut=mapElementPositionsOnmRNA(elementPositions);
  return tut[which(tut<sd[gArray<gPos>(0)])];
}

/** @brief Mapping of element positions to the coding sequence
 *
 * Maps the positions to the mRNA coding sequence returning a gArray containing positions from the beginning of the coding sequence.
 * If a position is outside the mRNA 5'UTR it returns a NA.
 * @param elementPositions array of element positions
 * @return a garray<gPos> object
 */
gArray<gRelativePos> gTranscript::mapElementPositionsOnCds(const gArray<gRelativePos> & elementPositions) const {
  gArray<gRelativePos> sd=mapElementPositionsOnmRNA(getCdsInterval().getPositions());
  gArray<gRelativePos> tut=mapElementPositionsOnmRNA(elementPositions);
  return tut-sd[gArray<gPos>(0)];
}

/** @brief Mapping of element positions to the mRNA 3'UTR
 *
 * Maps the positions to the mRNA 3'UTR returning a gArray containing positions from the beginning of the mRNA 3'UTR.
 * If a position is outside the mRNA 5'UTR it returns a NA.
 * @param elementPositions array of element positions
 * @return a garray<gPos> object
 */
gArray<gRelativePos> gTranscript::mapElementPositionsOn3UTR(const gArray<gRelativePos> & elementPositions) const{
  gArray<gRelativePos> sd=mapElementPositionsOnmRNA(getCdsInterval().getPositions());
  gArray<gRelativePos> tut=mapElementPositionsOnmRNA(elementPositions);
  return tut[which(tut>=sd[gArray<gPos>(1)])]-sd[gArray<gPos>(1)];
}

/** @brief Mapping of cds positions onto the coding sequence
 *
 * Maps mRNA coding sequence positions to the corresponding element positions
 * If a position is outside the mRNA 5'UTR it returns a NA.
 * @param cdsPositions array of cds positions
 * @return a garray<gRelativePos> object
 */
gArray<gRelativePos> gTranscript::mapCdsPositionsOnElement(const gArray<gRelativePos> & cdsPositions) const{
  gArray<gRelativePos> ret(0,cdsPositions.getSize(),true);
  for(gPos i=0;i<getExonCount();i++){
   gElementInterval exon=getExonInterval(i);
   gElementInterval exonmap=mapElementIntervalOnCds(exon);
   gArray<gPos> gp=which(exonmap.contains(cdsPositions));
   ret.setValues(gp,cdsPositions[gp]-exonmap.getStart()+exon.getStart());
  }
  return ret;
}

/** @brief Mapping of element interval to the mRNA
 *
 * Maps the interval to the mRNA returning a gArray containing the interval mappend onto the mRNA.
 * If a position is outside the mRNA it returns a NA.
 * @param elementInterval two values array with starting and ending position of the interval
 * @return a garray<gPos> object
 */
gElementInterval  gTranscript::mapElementIntervalOnmRNA(const gElementInterval & elementInterval) const{
  gArray<gPos> rr=mapElementPositionsOnmRNA(elementInterval.getPositions());
  return gElementInterval(rr[0],rr[1]+1,!rr.isNA(0),!rr.isNA(1));
}

/** @brief Mapping of element interval to the mRNA 5'UTR
 *
 * Maps the interval to the mRNA 5'UTR returning a gArray containing the interval mappend onto the mRNA 5'UTR.
 * If a position is outside the mRNA 5'UTR it returns a NA.
 * @param elementInterval two values array with starting and ending position of the interval
 * @return a garray<gPos> object
 */
gElementInterval  gTranscript::mapElementIntervalOn5UTR(const gElementInterval & elementInterval) const{
  gArray<gPos> rr=mapElementPositionsOn5UTR(elementInterval.getPositions());
  return gElementInterval(rr[0],rr[1]+1,!rr.isNA(0),!rr.isNA(1));
}

/** @brief Mapping of element interval to the coding sequence
 *
 * Maps the interval to the coding sequence returning a gArray containing the interval mappend onto the coding sequence.
 * If a position is outside the coding sequence it returns a NA.
 * @param elementInterval two values array with starting and ending position of the interval
 * @return a garray<gPos> object
 */
gElementInterval  gTranscript::mapElementIntervalOnCds(const gElementInterval & elementInterval) const{
  gArray<gPos> rr=mapElementPositionsOnCds(elementInterval.getPositions());
  return gElementInterval(rr[0],rr[1]+1,!rr.isNA(0),!rr.isNA(1));
}

/** @brief Mapping of element interval to the mRNA 3'UTR
 *
 * Maps the interval to the mRNA 3'UTR returning a gArray containing the interval mappend onto the mRNA 3'UTR.
 * If a position is outside the mRNA 3'UTR it returns a NA.
 * @param elementInterval two values array with starting and ending position of the interval
 * @return a garray<gPos> object
 */
gElementInterval  gTranscript::mapElementIntervalOn3UTR(const gElementInterval & elementInterval) const{
  gArray<gPos> rr=mapElementPositionsOn3UTR(elementInterval.getPositions());
  return gElementInterval(rr[0],rr[1]+1,!rr.isNA(0),!rr.isNA(1));
}

/** @brief If element positions belong to the mRNA
 *
 * Checks if the position passed are inside the mRNA
 *
 * @param elementPositions array of element positions
 * @return a gArray<gBool> object
 */
gArray<gBool> gTranscript::ismRNA(gArray<gPos> elementPositions) const{
  gArray<gBool> ret(0,elementPositions.getSize(),false);
  for(gArrayIndex i=0;i<getExonCount();i++){
    ret=ret||getExonInterval(i).contains(elementPositions);
  }
  return ret;
}

/** @brief If element positions belong to the mRNA 5'UTR
 *
 * Checks if the position passed are inside the mRNA 5'UTR
 *
 * @param elementPositions array of element positions
 * @return a gArray<gBool> object
 */
gArray<gBool> gTranscript::is5UTR(gArray< gPos > elementPositions) const{
  gArray<gPos> sd=mapElementPositionsOnmRNA(getCdsInterval().getPositions());
  gArray<gPos> tut=mapElementPositionsOnmRNA(elementPositions);
  return tut<sd[gArray<gPos>(0)];
}

/** @brief If element positions belong to the coding sequence
 *
 * Checks if the position passed are inside the coding sequence
 *
 * @param elementPositions array of element positions
 * @return a gArray<gBool> object
 */
gArray<gBool> gTranscript::isCds(gArray<gPos> elementPositions) const{
 gArray<gPos> sd=mapElementPositionsOnmRNA(getCdsInterval().getPositions());
 gArray<gPos> tut=mapElementPositionsOnmRNA(elementPositions);
 return (tut>=sd[gArray<gPos>(0)])&&(tut<sd[gArray<gPos>(1)]);
}

/** @brief If element positions belong to the mRNA 3'UTR
 *
 * Checks if the position passed are inside the mRNA 3'UTR
 *
 * @param elementPositions array of element positions
 * @return a gArray<gBool> object
 */
gArray<gBool> gTranscript::is3UTR(gArray<gPos> elementPositions) const{
 gArray<gPos> sd=mapElementPositionsOnmRNA(getCdsInterval().getPositions());
 gArray<gPos> tut=mapElementPositionsOnmRNA(elementPositions);
 return tut>=sd[gArray<gPos>(1)];
}

/** @brief If element positions belong to an exon
 *
 * Checks if the position passed are inside an exon
 *
 * @param elementPositions array of element positions
 * @return a gArray<gBool> object
 */
gArray<gBool> gTranscript::isExon(gArray<gPos> elementPositions) const{
  gArray<gBool> ret;
  throw gException("gTranscript::isExon: not yet implemented");
  return ret;
}

/** @brief If element positions belong to an intron
 *
 * Checks if the position passed are inside an intron
 *
 * @param elementPositions array of element positions
 * @return a gArray<gBool> object
 */
gArray<gBool> gTranscript::isIntron(gArray<gPos> elementPositions) const{
  gArray<gBool> ret;
  throw gException("gTranscript::isIntron: not yet implemented");
  return ret;
}


/** @brief Protein length
 *
 * Returns the length (in aa) of the protein associated to this transcript
 *
 * @return a gArray<gSize> object
 */
gArray<gSize> gTranscript::getProteinLength() const{
  //should we check for stop codons possibly intruduced by variations????
  return getCdsLength()/3;
}

/** @brief Protein sequence retriever
 *
 * Returns a gString object containing a sub interval of the protein sequence
 * @param interval the interval
 * @return a gSequence object
 */
gString gTranscript::getProteinSequence(const gElementInterval & interval) const{
  gSequence seq=getCdsSequence().getTranslated(0);
  return gString(seq,interval.getStart(),min((gSize)interval.getEnd(),(gSize)seq.getLength()));
}

/** @brief Element positions mapping onto the protein sequence
 *
 * Maps the positions to the protein sequence returning a gArray containing positions from the beginning of the protein sequence.
 * If a position is outside the cds it returns a NA.
 * @param elementPositions array of element positions
 * @return a garray<gPos> object
 */
gArray<gPos> gTranscript::mapElementPositionsOnProtein(const gArray<gPos> & elementPositions) const{
  return mapElementPositionsOnCds(elementPositions)/3;
}

/** @brief Element positions mapping onto the protein sequence
 *
 * Maps the interval to the protein sequence returning a gArray containing the interval mapped onto the protein sequence.
 * If a position is outside the cds it returns a NA.
 * @param elementInterval two values array with starting and ending position of the interval
 * @return a garray<gPos> object
 */
gElementInterval gTranscript::mapElementIntervalOnProtein(const geco::gElementInterval& elementInterval) const{
  gArray<gPos> rr=mapElementPositionsOnProtein(elementInterval.getPositions());
  return gElementInterval(rr[0],rr[1]+1,!rr.isNA(0),!rr.isNA(1));
}

/** @brief Codon intrvals onto element
 *
 * Maps protein positions to the corresponding element positions (codon start)
 * @param proteinPositions array containing positions to be mapped
 * @return a garray<gPos> object
 */
gArray<gRelativePos> gTranscript::mapProteinPositionOnElement(const gArray<gPos> & proteinPositions) const{
  return mapCdsPositionsOnElement(proteinPositions*3);
}
//--------------------------------------------------------------------------
