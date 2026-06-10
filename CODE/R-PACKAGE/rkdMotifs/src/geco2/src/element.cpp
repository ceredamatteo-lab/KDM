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
#include "element.h"

using namespace std;
using namespace geco;

//--------------------------------------------------------------------------
/// gVariations implementation
//--------------------------------------------------------------------------
/** @brief Empty constructor
 *
 * Instatiate an enpty gVariations object
 */
gVariations::gVariations() {
}

/** @brief From reference constructor
 *
 * Instantiate a gVariations relative to a reference.
 * @param referenceName the reference name
 */
gVariations::gVariations(const gString & referenceName) {
    i_referenceName=referenceName;
}

/** @brief From element constructor
 *
 * Instantiate a gVariations realtive to an element reference.
 * @param element the element whos reference must be used
 */
gVariations::gVariations(const gElement & element) {
    i_referenceName=element.getReferenceName();
}

/** @brief Copy contructor
 *
 * Instatiate a copy of the provided object.
 * @param variations the object to copy from
 */
gVariations::gVariations(const gVariations & variations) {
    *this=variations;
}

/** @brief Assignment operator
 *
 * Copies everything from a provided object
 * @param variations the object to be copied into this one
 * @return a reference to this
 */
gVariations & gVariations::operator = (const gVariations & variations) {
    i_referenceName=variations.i_referenceName;
    i_varStart=variations.i_varStart;
    i_varEnd=variations.i_varEnd;
    i_varType=variations.i_varType;
    i_varSeq=variations.i_varSeq;
    i_pos=variations.i_pos;
    return *this;
}

/** @brief Variations count
 *
 * Returns the number of variations contained in this object
 * @return a gSize value
 */
gSize gVariations::getCount() const {
    return i_varStart.getSize();
}

/** @brief Reference name
 *
 * Returns the name of the reference to which variations are relative
 * @return a gString
 */
gString gVariations::getReferenceName() const {
    return i_referenceName;
}

/** @brief Subsetting operartor
 *
 * Returns a new object ogf type gVariations containing the variations
 * at the positions specified.
 * @param varNum an array containing the position in this object of the var to be retrieved
 * @return a newly instatiated object containing the specified variations
 */
gVariations gVariations::operator [] (const gArray<gArrayIndex> & varNum) const {
    gArray<gArrayIndex> avarNum=i_pos[varNum];
    gVariations ret(i_referenceName);
    ret.i_pos=getArray<gPos>(0,varNum.getSize()-1,1);
    ret.i_varStart=i_varStart[avarNum];
    ret.i_varEnd=i_varEnd[avarNum];
    ret.i_varType=i_varType[avarNum];
    for (gPos i=0;i<avarNum.getSize();i++) ret.i_varSeq.push_back(i_varSeq[avarNum[i]]);
    return ret;
}

/** @brief Add an insertion relative to an element
 *
 * Add an the insertion of seq at the given position referred to an element
 * @param element the element to which position is to be referred
 * @param elementPosition th position after which seq should be inserted
 * @param seq the insertion sequence
 * @return the index of the added insertion
 */
gArrayIndex gVariations::addInsertion(const gElement & element,gRelativePos elementPosition,const gSequence &seq) {
    gPos p=(element.isForward())?(elementPosition):(elementPosition+1);//QUESTO????
    return addVariation(element,gElementInterval(p,p+1),gInsertion,seq);
}

/** @brief Add a deletion relative to an element
 *
 * Add a deletion referred to an element
 * @param element the element to which positions are referred
 * @param elementInterval element interval correspnding to the deletion
 * @return the index of the added deletion
 */
gArrayIndex gVariations::addDeletion(const gElement & element,const gElementInterval & elementInterval) {
    return addVariation(element,elementInterval,gDeletion,gSequence());
}

/** @brief Add a substitution relative to an element
 *
 * Add a Substitution referred to an element
 * @param element the element to which positions are to be referred
 * @param elementInterval element interval corresponding to the substitution
 * @param seq the substituting sequence
 * @return the index of the added substitution
 */
gArrayIndex gVariations::addSubstitution(const gElement & element,const gElementInterval & elementInterval,const gSequence & seq) {
    return addVariation(element,elementInterval,gSubstitution,seq);
}

/** @brief Add an insertion relative to a rference
 *
 * Add an insertion of a sequence at the given reference position
 * @param referenceName the reference to which position is to be referred
 * @param referencePosition th position at which seq should be inserted
 * @param seq the insertion sequence
 * @return the index of the added insertion
 */
gArrayIndex gVariations::addInsertion(const gString & referenceName,gPos referencePosition,const gSequence &seq) {
    return addVariation(referenceName,gReferenceInterval(referencePosition,referencePosition+1),gInsertion,seq);
}

/** @brief Add a deletion relative to a reference
 *
 * Add a deletion of a given reference interval
 * @param referenceName the reference to which positions are to be referred
 * @param referenceInterval the reference interval to delete
 * @return the index of the added deletion
 */
gArrayIndex gVariations::addDeletion(const gString & referenceName,const gReferenceInterval & referenceInterval) {
    return addVariation(referenceName,referenceInterval,gDeletion,gSequence());
}

/** @brief Add a substitution relative to a reference
 *
 * Add a Substitution of the sequence of a given reference interval
 * @param referenceName the reference name to which positions are to be referred
 * @param referenceInterval the reference interval to substitute
 * @param seq the substituting sequence
 * @return the index of the added substitution
 */
gArrayIndex gVariations::addSubstitution(const gString & referenceName,const gReferenceInterval & referenceInterval,const gSequence &seq) {
    return addVariation(referenceName,referenceInterval,gSubstitution,seq);
}

/** @brief Type of a variation
 *
 * Returns the type of variation: gVariationType
 * @param varNum the number of the variation
 * @return type of variation
 */
gVariationType gVariations::getVariationType(gArrayIndex varNum) const {
    return (gVariationType) i_varType[i_pos[varNum]];
}

/** @brief Variation start
 *
 * Returns the reference starting position of a variation
 * @param varNum the number of the variation
 * @return the starting position (0 based)
 */
gPos gVariations::getVariationStart(gArrayIndex varNum) const {
    return i_varStart[i_pos[varNum]];
}

/** @brief Variation end
 *
 * Returns the reference ending position of a variation
 * @param varNum the number of the variation
 * @return the ending position (1 based)
 */
gPos gVariations::getVariationEnd(gArrayIndex varNum) const {
    return i_varEnd[i_pos[varNum]];
}

/** @brief Variation Interval
 *
 * Returns the reference interval corresponding to the variation
 * @param varNum the number of the variation
 * @return a reference interval
 */
gReferenceInterval gVariations::getVariationInterval(gArrayIndex varNum) const {
    return gReferenceInterval(i_varStart[i_pos[varNum]],i_varEnd[i_pos[varNum]]);
}

/** @brief Sequence of a variation
 *
 * Returns the variating (deleted for deletions) sequence corresponding to a given variation
 * @param varNum the number of the variation
 * @return a gSequence
 */
gSequence gVariations::getVariationSequence(gArrayIndex varNum) const {
    return i_varSeq[i_pos[varNum]];
}

//private THIS FUNCTION SHOULD BE OPTIMIZED FOR MULTIPLE INSERTIONS (EACH TIME A VAIATION IS INSERTED
// EVERYTHING IS SORTED WHILE KEEPING INDEXES)
gArrayIndex gVariations::addVariation(const gString & referenceName,const gReferenceInterval & referenceInterval,gVariationType type,const gSequence & seq) {
    gArrayIndex np=i_varStart.getSize();
    if (referenceName==i_referenceName){
        i_varType.setValue(np,type,false);
        i_varStart.setValue(np,referenceInterval.getStart(),!referenceInterval.validStart());
        i_varEnd.setValue(np,referenceInterval.getEnd(),!referenceInterval.validEnd());
        i_varSeq.push_back(seq);
        i_pos.setValue(np,np,false);
        if(referenceInterval.getStart()<i_varStart[np-1]){
          gArray<gPos> index=i_varStart.sort(true);
          i_varType=i_varType[index];
          i_varEnd=i_varEnd[index];
          vector<gSequence> newSeqs;
          for (gSize i=0;i<index.getSize();i++) newSeqs.push_back(i_varSeq[index[i]]);
          i_varSeq=newSeqs;
          i_pos=i_pos[index].sort(true);
        }
    } else throw gException("gVariations::addVariation, reference mismatch");
    return np;
}

gArrayIndex gVariations::addVariation(const gElement & element,const gElementInterval & elementInterval,gVariationType type,const gSequence &seq) {
    if (!(element.i_region->getCount()>0)) {
        gSequence tseq=seq;
        if (!element.isForward()) tseq.reverseComplement();
        return addVariation(element.getReferenceName(),element.getReferenceIntervalFromElement(elementInterval),type,tseq);
    } else throw gException("gVariations::addVariation, variations cannot be refferred to varied elements");
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// gSequenceRetrieverInternal implementation
//--------------------------------------------------------------------------
//This class is used internally by gRegion when no sequence retriever is passed
//it is initialized with a sequence (that can be zero length) and pad with "N" when
//required to return a sequence extendinf
class gSequenceRetrieverInternalImplementation:public gSequenceRetrieverImplementation,public gReferenceInterval {
private:
    gSequence i_sequence;
    virtual gSequence getSequence_Internal(const gString & reference,gPos start,gPos end) const;
protected:
public:
    gSequenceRetrieverInternalImplementation(const gReferenceInterval & referenceInterval,const gSequence & sequence);
    ~gSequenceRetrieverInternalImplementation();
    virtual gArrayRetrieverImplementation<gChar> * clone() const;
};

gSequenceRetrieverInternalImplementation::gSequenceRetrieverInternalImplementation(const gReferenceInterval & referenceInterval,const gSequence & sequence):gSequenceRetrieverImplementation(),gReferenceInterval(referenceInterval) {
    if (referenceInterval.getLength()>sequence.getLength()) {
        i_sequence=sequence+gSequence(geco::gArray<gChar>('N',referenceInterval.getLength()-sequence.getLength(),true));
    } else {
        i_sequence=gSequence(sequence,0,referenceInterval.getLength());
    }
}

gSequenceRetrieverInternalImplementation::~gSequenceRetrieverInternalImplementation(){
}

gArrayRetrieverImplementation<gChar> * gSequenceRetrieverInternalImplementation::clone() const {
    return new gSequenceRetrieverInternalImplementation(*this,i_sequence);
}

gSequence gSequenceRetrieverInternalImplementation::getSequence_Internal(const gString & reference,gPos start,gPos end) const {
    gReferenceInterval b(start,end);
    gReferenceInterval good=getIntersection(b);
    gString goodseq;
    if (good.getLength()>0) goodseq=gSequence(i_sequence,good.getStart()-getStart(),good.getEnd()-getStart());
    return gSequence(geco::gArray<gChar>('N',getLeftDiff(b).getLength(),true))+goodseq+gSequence(geco::gArray<gChar>('N',getRightDiff(b).getLength(),true));
}
//--------------------------------------------------------------------------



//--------------------------------------------------------------------------
/// gRegion implementation
//--------------------------------------------------------------------------
gRegion::gRegion(const gString & referenceName,const gReferenceInterval & referenceInterval,const gVariations & variations):gReferenceInterval(referenceInterval),gVariations(referenceName),i_retriever(gSequenceRetrieverInternalImplementation(referenceInterval,gSequence())){
    if (variations.getCount()>0) {
        if (getReferenceName()==variations.getReferenceName()) *((gVariations *) this)=variations;
        else throw gException("gRegion::gRegion: reference mismatch");
    }
    i_mode=gRet;
    //i_retriever= new gSequenceRetrieverInternalImplementation(referenceInterval,gSequence());
    i_refCount=0;
}

gRegion::gRegion(const gString & referenceName,const gReferenceInterval & referenceInterval,const gSequence & sequence,gElementSequenceMode mode,const gVariations & variations):gReferenceInterval(referenceInterval),gVariations(referenceName),i_retriever(gSequenceRetrieverInternalImplementation(referenceInterval,sequence)) {
    if (variations.getCount()>0) {
        if (getReferenceName()==variations.getReferenceName()) *((gVariations *) this)=variations;
        else throw gException("gRegion::gRegion: reference mismatch");
    }
    //i_retriever= new gSequenceRetrieverInternal(referenceInterval,sequence);
    i_mode=mode;
    //gArray<gPos> regionPositions=regionInterval.getPositions();
    geco::gArray<gPos> regionDisplacements(0,2,false);
    switch (mode) {
    case gRet:
        break;
    case gRef:
        i_sequence=i_retriever.getSequence(referenceName,referenceInterval.getStart(),referenceInterval.getEnd());
        break;
    case gReg:
        i_sequence=variateSequence(i_retriever.getSequence(referenceName,referenceInterval.getStart(),referenceInterval.getEnd()),referenceInterval,regionDisplacements);
        break;
    case gElm:
        i_mode=gRet;
        break;
    }
    i_refCount=0;
}

gRegion::gRegion(const gString & referenceName,const gReferenceInterval & referenceInterval,const gSequenceRetriever & sequenceRetriever,gElementSequenceMode mode,const gVariations & variations):gReferenceInterval(referenceInterval),gVariations(referenceName),i_retriever(sequenceRetriever) {
    if (variations.getCount()>0) {
        if (getReferenceName()==variations.getReferenceName()) *((gVariations *) this)=variations;
        else throw gException("gRegion::gRegion: reference mismatch");
    }
    //i_retriever = sequenceRetriever.clone();
    i_mode=mode;
    geco::gArray<gPos> regionDisplacements(0,2,false);
    switch (i_mode) {
    case gRet:
        break;
    case gRef:
        i_sequence=i_retriever.getSequence(referenceName,referenceInterval.getStart(),referenceInterval.getEnd());
        break;
    case gReg:
        i_sequence=variateSequence(i_retriever.getSequence(referenceName,referenceInterval.getStart(),referenceInterval.getEnd()),referenceInterval,regionDisplacements);
        break;
    case gElm:
        i_mode=gRet;
        break;
    }
    i_refCount=0;
}

gRegion::gRegion(const gRegion & region,const gVariations & variations):gReferenceInterval(region),gVariations(region),i_retriever(region.i_retriever){
    if (variations.getCount()>0) {
        if (getCount()==0) {
            if (getReferenceName()==variations.getReferenceName()) {
                *((gVariations *) this)=variations;
            } else throw gException("gRegion::gRegion: reference mismatch");
        } else throw gException("gRegion::gRegion: cannot add variations to a varied object");
    }
    //i_retriever=region.i_retriever->clone();
    i_mode=region.i_mode;
    geco::gArray<gPos> regionDisplacements(0,2,false);
    switch (i_mode) {
    case gRet:
        break;
    case gRef:
        i_sequence=region.i_sequence;
        break;
    case gReg:
        if (variations.getCount()==0) {
            i_sequence==region.i_sequence;
        } else {
            geco::gArray<gPos> bpos=getPositions();
            bpos.setNoNa();
            geco::gArray<gPos> regPos=getRegionPositionsFromReference(bpos);
            if (regPos.isNA(0)) regPos.setValue(0,regPos[0]+1,true);
            gReferenceInterval refIntv(getReferencePositionsFromRegion(regPos,regionDisplacements));
            //gReferenceInterval refIntv=getReferencePositionsFromRegion(region.getRegionPositionsFromReference(getPositions()),regionDisplacements);
            i_sequence=variateSequence(gSequence(region.i_sequence,refIntv.getStart()-region.getStart(),refIntv.getEnd()-region.getStart()),refIntv,regionDisplacements);
        }
        break;
    case gElm: //this is managed by gElement
        break;
    }
    i_refCount=0;
}

gRegion::~gRegion() {
}

geco::gArray<gPos> gRegion::getReferencePositionsFromRegion(const geco::gArray<gPos> & regionPositions,geco::gArray<gPos> & regionDisplacements) const {
    geco::gArray<gPos> ret(regionPositions.getSize(),0,true);
    gArrayIndex nv=getCount();
    bool invalid=false;
    gPos vstart,vend,vlen;
    gVariationType vtype;
    gRelativePos diff;
    gPos disp;
    gPos position;
    for (gSize pos=0; pos< regionPositions.getSize(); pos++) {
        invalid=false;
        diff=0;
        disp=0;
        position=regionPositions[pos];
        for (gArrayIndex i=0;i<nv;i++) {
            vstart=i_varStart[i];
            vend=i_varEnd[i];
            vtype=i_varType[i];
            if (vtype==gInsertion) {
                vlen=i_varSeq[i].getLength();
                //if insertion happened before
                if (vstart+vlen+diff<position) {
                    diff+=vlen;
                } else if (vstart+diff<position) {
                    disp=position-diff-vstart;
                    diff=position-vstart;
                    invalid=true;
                    break;
                } else break;
            } else if (vtype==gDeletion) {
                if (vstart+diff<=position) {
                    diff-=vend-vstart;
                }
            }
        }
        if (diff>(gRelativePos)position) throw gException("gRegion:getReferencePositions->Unhandled Error");
        ret.setValue(pos,position-diff,invalid||regionPositions.isNA(pos));
        regionDisplacements.setValue(pos,disp,false);
    }
    return ret;
}

geco::gArray<gPos> gRegion::getRegionPositionsFromReference(const geco::gArray<gPos> & referencePositions) const {
    geco::gArray<gPos> ret(referencePositions.getSize(),0,true);
    for (gSize pos=0; pos<referencePositions.getSize(); pos++) {
        gPos position=referencePositions[pos];
        bool invalid=false;
        gRelativePos diff=0;
		gArrayIndex nv=getCount();
        gArrayIndex i=0;
        while (i<nv) {
            if (i_varType[i]==gInsertion) {
                if (i_varEnd[i]<=position) {
                    diff+=i_varSeq[i].getLength();
                } else break;
            } else if (i_varType[i]==gDeletion) {
                if (i_varEnd[i]<=position) {
                    diff-=i_varEnd[i]-i_varStart[i];
                } else if (i_varStart[i]<=position) {
                    diff-=position-i_varStart[i]+1;
                    invalid=true;
                    break;
                } else break;
            }
            i++;
        }
        if (-diff>(gRelativePos) position) throw gException("Unhandled Error");
        ret.setValue(pos,position+diff,invalid||referencePositions.isNA(pos));
    }
    return ret;
}

gSequence gRegion::getSequenceFromRetriever(const gReferenceInterval & regionInterval) const {
    geco::gArray<gPos> regionDisplacements(0,2,false);
//  gArray<gPos> bpos=regionInterval.getPositions();
//  bpos.setNoNa();
//  gArray<gPos> refPos=getReferencePositionsFromRegion(bpos,regionDisplacements);
//  if(refPos.isNA(0)) refPos.setValue(0,refPos[0]+1,true);
    gReferenceInterval refIntv=gReferenceInterval(getReferencePositionsFromRegion(regionInterval.getPositions(),regionDisplacements));
    return variateSequence(i_retriever.getSequence(getReferenceName(),refIntv.getStart(),refIntv.getEnd()),refIntv,regionDisplacements);
}

gSequence gRegion::getSequenceFromReference(const gReferenceInterval & regionInterval) const {
    gString goodseq,leftseq,rightseq;
    //gReferenceInterval thisRegionInterval=gReferenceInterval(getRegionPositionsFromReference(getPositions()));
    geco::gArray<gPos> bpos=getPositions();
    bpos.setNoNa();
    geco::gArray<gPos> ipos=getRegionPositionsFromReference(bpos);
    if (ipos.isNA(0)) ipos.setValue(0,ipos[0]+1,true);
    gReferenceInterval thisRegionInterval=gReferenceInterval(ipos);


    gReferenceInterval good=thisRegionInterval.getIntersection(regionInterval);
    gReferenceInterval left=thisRegionInterval.getLeftDiff(regionInterval);
    gReferenceInterval right=thisRegionInterval.getRightDiff(regionInterval);
    if (good.getLength()>0) {
        geco::gArray<gPos> disp(0,2,false);
        gReferenceInterval refIntv=gReferenceInterval(getReferencePositionsFromRegion(good.getPositions(),disp));
        gSequence seq=gSequence(i_sequence,refIntv.getStart()-getStart(),refIntv.getEnd()-getStart());
        goodseq=variateSequence(seq,refIntv,disp);
    }
    if (left.getLength()>0) {
        leftseq=getSequenceFromRetriever(left);
    }
    if (right.getLength()>0) {
        rightseq=getSequenceFromRetriever(right);
    }
    gSequence ret=leftseq+goodseq+rightseq;
    return ret;
}

gSequence gRegion::getSequenceFromRegion(const gReferenceInterval & regionInterval) const {
    gString goodseq,leftseq,rightseq;

    //gReferenceInterval thisRegionInterval=gReferenceInterval(getRegionPositionsFromReference(getPositions()));
    geco::gArray<gPos> bpos=getPositions();
    bpos.setNoNa();
	geco::gArray<gPos> ipos=getRegionPositionsFromReference(bpos);
    if (ipos.isNA(0)) ipos.setValue(0,ipos[0]+1,true);
	gReferenceInterval thisRegionInterval=gReferenceInterval(ipos);

	gReferenceInterval good=thisRegionInterval.getIntersection(regionInterval);
	gReferenceInterval left=thisRegionInterval.getLeftDiff(regionInterval);
	gReferenceInterval right=thisRegionInterval.getRightDiff(regionInterval);

	if (good.getLength()>0) {
		goodseq=gSequence(i_sequence,good.getStart()-thisRegionInterval.getStart(),good.getEnd()-thisRegionInterval.getStart());
	}
	if (left.getLength()>0) {
		leftseq=getSequenceFromRetriever(left);
	}
	if (right.getLength()>0) {
		rightseq=getSequenceFromRetriever(right);
	}
	return leftseq+goodseq+rightseq;
}

gSequence gRegion::getSequence(const gReferenceInterval & regionInterval) const {
	gSequence ret;
	switch (i_mode) {
	case gRet:
		ret=getSequenceFromRetriever(regionInterval);
		break;
	case gRef:
		ret=getSequenceFromReference(regionInterval);
		break;
	case gReg:
		ret=getSequenceFromRegion(regionInterval);
		break;
	case gElm://this is managed by gElement
		break;		
	}


	return ret;
}

gSequence gRegion::variateSequence(const gSequence & seq,const gReferenceInterval & referenceInterval,const geco::gArray<gPos> & regionDisplacements) const {
	gSequence ret=seq;
	if (getCount()>0) {
		geco::gArray<gPos> overlapping=which((i_varEnd > geco::gArray<gPos>(referenceInterval.getStart(),1,false))&&(i_varStart < geco::gArray<gPos>(referenceInterval.getEnd(),1,false)));
		gRelativePos diff=0;
		for (gPos j=0; j<overlapping.getSize();j++) {
            gPos vnum=overlapping[j];
            gPos vstart=(i_varStart[vnum]>referenceInterval.getStart())?(i_varStart[vnum]):(referenceInterval.getStart());
            gPos vend=(i_varEnd[vnum]<referenceInterval.getEnd())?(i_varEnd[vnum]):(referenceInterval.getEnd());
            gSequence vseq=i_varSeq[vnum];
            gPos np=vstart - referenceInterval.getStart() + diff;
            switch (i_varType[vnum]) {
            case gInsertion:
                if ((vstart==referenceInterval.getStart())&&(vend!=referenceInterval.getEnd())) {
                    if ((np==0)&&(!referenceInterval.validStart())) {
                        if(regionDisplacements[0]<vseq.getLength()){
                         ret=gSequence(vseq,regionDisplacements[0]-1,vseq.getLength())+gSequence(ret,1,ret.getLength());
                        }else{
                         ret=gSequence(vseq,regionDisplacements[0]-1,vseq.getLength())+gSequence(ret,1,ret.getLength());
                        }
                    } else {
                      if(regionDisplacements[0]<vseq.getLength()){
                        ret.insert(np,gSequence(vseq,regionDisplacements[0],vseq.getLength()));
                      }else{
                        ret.insert(np,gSequence(vseq,regionDisplacements[0]-1,vseq.getLength()));
                      }
                    }
                    diff+=vseq.getLength()-regionDisplacements[0];
                } else if ((vend==referenceInterval.getEnd())&&(vstart!=referenceInterval.getStart())) {
                    ret.insert(np,gSequence(vseq,0,regionDisplacements[1]));
                    diff+=regionDisplacements[1];
                } else if ((vstart==referenceInterval.getStart())&&(vend==referenceInterval.getEnd())) {
                    //ret=gSequence(vseq,regionDisplacements[0],regionDisplacements[1]+1);
                    ret=gSequence(vseq,regionDisplacements[0],regionDisplacements[1]);
                    diff+=regionDisplacements[1]-regionDisplacements[0];
                } else {
                    ret.insert(np,vseq);
                    diff+=vseq.getLength();
                }
                break;
            case gDeletion:
                ret.remove(np,vend-vstart);
                diff-=vend-vstart;
                break;
            case gSubstitution:
                ret.replace(np,vseq);
                break;
            }
        }
    }
    return ret;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/// gBaseFeature class implementation
//--------------------------------------------------------------------------
gFeatureBase::gFeatureBase(const gString & name){
  i_name=name;
}

gFeatureBase::gFeatureBase(const gFeatureBase & feature){
  i_name=feature.i_name;
}

gFeatureBase::~gFeatureBase(){
}

const gString & gFeatureBase::getName() const{
  return i_name;
}
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
/// gElement implementation
//--------------------------------------------------------------------------
/** @brief Empty contructor
 *
 * Construct an empty gElement object
 */
gElement::gElement() {
    i_region=NULL;
}

/** @brief Constructor from reference
 *
 * Instatiate an element corresponding to a gReferenceInterval
 * @param referenceName name of the reference
 * @param referenceInterval reference interval
 * @param forward strand of the element referred to the reference
 * @param variations
 */
gElement::gElement(const gString & referenceName,const gReferenceInterval & referenceInterval,bool forward,const gVariations & variations) {
    i_region=new gRegion(referenceName,referenceInterval,variations);
    i_region->i_refCount++;
    i_forward=forward;
    i_boundaries=referenceInterval;
    i_boundaries_reg_disp=gArray<gPos>(0,2,false);
}

/** @brief Constructor from reference
 *
 * Instatiate an element corresponding to a starting and ending reference postion
 * @param referenceName name of the reference
 * @param referenceStart reference starting position of the element (0 based)
 * @param referenceEnd referrence ending position of the element (1 based)
 * @param forward strand of the element referred to the reference
 * @param variations
 */
gElement::gElement(const gString & referenceName, gPos referenceStart, gPos referenceEnd, bool forward,const gVariations & variations) {
    gReferenceInterval referenceInterval(referenceStart,referenceEnd);
    i_region=new gRegion(referenceName,referenceInterval,variations);
    i_region->i_refCount++;
    i_forward=forward;
    i_boundaries=referenceInterval;
    i_boundaries_reg_disp=gArray<gPos>(0,2,false);
}

/** @brief Contructor from reference with sequence
 *
 * Instatiate an element with respect to a reference range.
 * A Sequence must be provided.
 * @param referenceName name of the reference
 * @param referenceInterval reference interval
 * @param forward strand of the element referred to the reference
 * @param referenceSequence the sequence
 * @param mode a gElementSequenceMode value to spcify the way the sequence must be allocated
 * @param variations
 */
gElement::gElement(const gString & referenceName,const gReferenceInterval & referenceInterval,bool forward,const gSequence & referenceSequence,gElementSequenceMode mode,const gVariations & variations) {
    if (referenceInterval.getLength()==referenceSequence.getLength()) {
        i_region=new gRegion(referenceName,referenceInterval,referenceSequence,mode,variations);
        i_region->i_refCount++;
        i_forward=forward;
        i_boundaries=referenceInterval;
        i_boundaries_reg_disp=gArray<gPos>(0,2,false);
        if (mode==gElm) i_sequence=getSequenceFromRegion(getElementInterval());
    } else if (referenceSequence.getLength()==0) {
        i_region=new gRegion(referenceName,referenceInterval,variations);
        i_region->i_refCount++;
        i_forward=forward;
        i_boundaries=referenceInterval;
        i_boundaries_reg_disp=gArray<gPos>(0,2,false);
    } else throw gException("gElement: sequence length mismatch in constructor");
}

/** @brief Contructor from reference with sequence
 *
 * Instatiate an element with respect to a reference range.
 * A Sequence must be provided.
 * @param referenceName name of the reference
 * @param referenceStart reference starting position of the element (0 based)
 * @param referenceEnd referrence ending position of the element (1 based)
 * @param forward strand of the element referred to the reference
 * @param referenceSequence the sequence
 * @param mode a gElementSequenceMode value to spcify the way the sequence must be allocated
 * @param variations
 */
gElement::gElement(const gString & referenceName, gPos referenceStart, gPos referenceEnd, bool forward,const gSequence & referenceSequence,gElementSequenceMode mode,const gVariations & variations) {
    gReferenceInterval referenceInterval(referenceStart,referenceEnd);
    if (referenceInterval.getLength()==referenceSequence.getLength()) {
        i_region=new gRegion(referenceName,referenceInterval,referenceSequence,mode,variations);
        i_region->i_refCount++;
        i_forward=forward;
        i_boundaries=referenceInterval;
        i_boundaries_reg_disp=gArray<gPos>(0,2,false);
        if (mode==gElm) i_sequence=getSequenceFromRegion(getElementInterval());
    } else if (referenceSequence.getLength()==0) {
        i_region=new gRegion(referenceName,referenceInterval,variations);
        i_region->i_refCount++;
        i_forward=forward;
        i_boundaries=referenceInterval;
        i_boundaries_reg_disp=gArray<gPos>(0,2,false);
    } else throw gException("gElement: sequence length mismatch in constructor");
}

/** @brief Contructor from reference with sequence retriever
 *
 * Instatiate an element with respect to a reference range.
 * A Sequence Retriever object must be provided.
 * @param referenceName name of the reference
 * @param referenceInterval reference interval
 * @param forward strand of the element referred to the reference
 * @param sequenceRetriever the sequence retriever object
 * @param mode a gElementSequenceMode value to spcify the way the sequence must be allocated
 * @param variations
 */
gElement::gElement(const gString & referenceName,const gReferenceInterval & referenceInterval,bool forward,const gSequenceRetriever & sequenceRetriever,gElementSequenceMode mode,const gVariations & variations) {
    i_region=new gRegion(referenceName,referenceInterval,sequenceRetriever,mode,variations);
    i_region->i_refCount++;
    i_forward=forward;
    i_boundaries=referenceInterval;
    i_boundaries_reg_disp=gArray<gPos>(0,2,false);
    if (mode==gElm) i_sequence=getSequenceFromRegion(getElementInterval());
}

/** @brief Contructor from reference with sequence retriever
 *
 * Instatiate an element with respect to a reference range.
 * A Sequence Retriever object must be provided.
 * @param referenceName name of the reference
 * @param referenceStart reference starting position of the element (0 based)
 * @param referenceEnd referrence ending position of the element (1 based)
 * @param forward strand of the element referred to the reference
 * @param sequenceRetriever the sequence retriever object
 * @param mode a gElementSequenceMode value to spcify the way the sequence must be allocated
 * @param variations
 */
gElement::gElement(const gString & referenceName, gPos referenceStart, gPos referenceEnd, bool forward, const gSequenceRetriever & sequenceRetriever, gElementSequenceMode mode, const gVariations & variations) {
    gReferenceInterval referenceInterval(referenceStart,referenceEnd);
    i_region=new gRegion(referenceName,referenceInterval,sequenceRetriever,mode,variations);
    i_region->i_refCount++;
    i_forward=forward;
    i_boundaries=referenceInterval;
    i_boundaries_reg_disp=gArray<gPos>(0,2,false);
    if (mode==gElm) i_sequence=getSequenceFromRegion(getElementInterval());
}

/** @brief Copy constructor
 *
 * Instatiate an element as a copy of a range from the provided element.
 * Strand modification and variations can be optionally applyed.
 * If provided, the variations object must have the same reference as the passed element
 * In this case the provide element cannot be already variated.
 * @param element the element from which copy informations.
 * @param elementInterval the elment interval
 * @param forward the strand respect the passed element (defaults to true=forward)
 * @param variations the variation to be applyed to the new object (defaults to an empty set)
 */
gElement::gElement(const gElement & element,const gElementInterval & elementInterval,bool forward,const gVariations & variations) {
    //strand & interval
    i_forward=(element.i_forward)?(forward):(!forward);
    i_boundaries_reg_disp=gArray<gPos>(0,2,false);
    i_boundaries=element.getReferenceIntervalFromElement(elementInterval,i_boundaries_reg_disp);

    //variations & region
    if (variations.getCount()>0) {
        if (!(element.i_region->getCount()>0)) {
            i_region = new gRegion(*(element.i_region),variations);
            i_region->i_refCount++;
        } else throw gException("gElement::gElement: cannot add variations to a vried object");
    } else {
        i_region=element.i_region;
        i_region->i_refCount++;
    }

    //sequence
    if (element.i_sequence.getLength()>0) {
        if (variations.getCount()>0) {
            gArray<gPos> regionDisplacements(0,2,false);
            gReferenceInterval refIntv=element.getReferenceIntervalFromElement(elementInterval,regionDisplacements);
            gSequence seq=element.getSequence(elementInterval);

            if (!element.i_forward) seq.reverseComplement();
            i_sequence=i_region->variateSequence(seq,refIntv,regionDisplacements);
            if (!i_forward) i_sequence.reverseComplement();
        } else {
            i_sequence=element.getSequence(elementInterval);
            if (!forward) i_sequence.reverseComplement();
        }
    }
    //Sites & connetcions
    i_site_list=element.i_site_list;
    i_site_list_reg_disp=element.i_site_list_reg_disp;
    if (forward) {
        i_start_site_index=element.i_start_site_index;
        i_end_site_index=element.i_end_site_index;
    } else {
        i_start_site_index=element.i_end_site_index;
        i_end_site_index=element.i_start_site_index;
    }

    //Features
    for(size_t i=0;i<element.i_features.size();i++){
      i_features.push_back(element.i_features[i]->clone(*this,element,elementInterval,variations));
    }
}

/** @brief Copy constructor
 *
 * Instatiate an element as a copy of a range from the provided element.
 * Strand modification and variations can be optionally applyed.
 * If provided, the variations object must have the same reference as the passed element
 * and this latter cannot be already variated.
 * @param element the element from which copy informations.
 * @param elementStart the element start position
 * @param elementEnd the element end position
 * @param forward the strand respect the passed element (defaults to true=forward)
 * @param variations the variation to be applyed to the new object (defaults to an empty set)
 */
gElement::gElement(const gElement & element,gRelativePos elementStart, gRelativePos elementEnd,bool forward,const gVariations & variations) {
  if(element.i_region){
    gElementInterval elementInterval(elementStart,(elementEnd==0)?(element.getLength()):(elementEnd));

    //strand & interval
    i_forward=(element.i_forward)?(forward):(!forward);
    i_boundaries_reg_disp=gArray<gPos>(0,2,false);
    i_boundaries=element.getReferenceIntervalFromElement(elementInterval,i_boundaries_reg_disp);

    //variations & region
    if (variations.getCount()>0) {
        if (!(element.i_region->getCount()>0)) {
            i_region = new gRegion(*(element.i_region),variations);
            i_region->i_refCount++;
        } else throw gException("gElement::gElement: cannot add variations to a vried object");
    } else {
        i_region=element.i_region;
        i_region->i_refCount++;
    }

    //sequence
    if (element.i_sequence.getLength()>0) {
        if (variations.getCount()>0) {
            gArray<gPos> regionDisplacements(0,2,false);
            gReferenceInterval refIntv=element.getReferenceIntervalFromElement(elementInterval,regionDisplacements);
            gSequence seq=element.getSequence(elementInterval);

            if (!element.i_forward) seq.reverseComplement();
            i_sequence=i_region->variateSequence(seq,refIntv,regionDisplacements);
            if (!i_forward) i_sequence.reverseComplement();
        } else {
            i_sequence=element.getSequence(elementInterval);
            if (!forward) i_sequence.reverseComplement();
        }
    }

    i_site_list=element.i_site_list;
    i_site_list_reg_disp=element.i_site_list_reg_disp;
    if (true) {
        i_start_site_index=element.i_start_site_index;
        i_end_site_index=element.i_end_site_index;
    } else {
        i_start_site_index=element.i_end_site_index;
        i_end_site_index=element.i_start_site_index;
    }

    for(size_t i=0;i<element.i_features.size();i++){
      i_features.push_back(element.i_features[i]->clone(*this,element,elementInterval,variations));
    }
  }else{
    i_region=NULL;
  }
}

/** @brief Destructor */
gElement::~gElement() {
    for(size_t i=0; i < i_features.size(); i++){
      delete i_features[i];
    }
    if (i_region) {
        i_region->i_refCount--;
        if (i_region->i_refCount==0) delete i_region;
    }
}

/** @brief Assignment operator
 *
 * Make this object an exact copy of the provided one.
 * @param elm The element to copy
 * @return a reference to this object
 */
gElement & gElement::operator=(const gElement & elm) {
  if(elm.i_region){
   i_region=elm.i_region;
   i_region->i_refCount++;
  }else{
    i_region=NULL;
  }
  i_forward=elm.i_forward;
  i_sequence=elm.i_sequence;
  i_boundaries=elm.i_boundaries;
  i_boundaries_reg_disp=elm.i_boundaries_reg_disp;
  i_site_list=elm.i_site_list;
  i_site_list_reg_disp=elm.i_site_list_reg_disp;
  i_start_site_index=elm.i_start_site_index;
  i_end_site_index=elm.i_end_site_index;
  for(size_t i=0;i<elm.i_features.size();i++){
    i_features.push_back(elm.i_features[i]->clone(*this,elm,elm.getElementInterval(),gVariations()));
  }
  return *this;
}

/** @brief Strand information
 *
 * Returns true if this element is forward (relative to the reference)
 * @return a bool value: true if it's forward
 */
bool gElement::isForward() const {
    return i_forward;
}

/** @brief Element genomic length
 *
 * Returns the acutal element genomic size (considering variations)
 * @return actual genomic size
 */
gSize gElement::getLength() const {
  if(i_region){
    return getElementInterval().getLength();
  }else{
    return 0;
  }
}

/** @brief Reference name
 *
 * Returns the name of the reference this element refers to
 * @return a string containing the reference name
 */
gString gElement::getReferenceName() const {
    return i_region->i_referenceName;
}

/** @brief Reference interval
 *
 * Returns then reference interval corresponding to the element
 * @return gElementInterval
 */
gReferenceInterval gElement::getReferenceInterval() const {
    bool vstart=i_boundaries.validStart();
    return gReferenceInterval((vstart)?(i_boundaries.getStart()):(i_boundaries.getStart()+1),i_boundaries.getEnd(),vstart,i_boundaries.validEnd());
}

/** @brief Element interval
 *
 * Returns the element interval corresponding to the element
 * @return gReferenceInterval
 */
gElementInterval gElement::getElementInterval() const {
    return getElementIntervalFromReference(i_boundaries,i_boundaries_reg_disp);
}

/** @brief Sequence
 *
 * Returns the sequence corresponding to the provided element interval
 * with optional start and end offsets
 * @param elementInterval the element interval
 * @param startOffset optional offset relative to the interval start (defaults to 0)
 * @param endOffset optional offset relative to the interval end (defaults to 0)
 * @return a gSequence object
 */
gSequence gElement::getSequence(const gElementInterval & elementInterval, gRelativePos startOffset, gRelativePos endOffset) const {
    gSequence ret;
    gElementInterval b(elementInterval,startOffset,endOffset);
    if (i_sequence.getLength()>0) { //Sequence is in gElm mode
        gElementInterval a=getElementInterval();
        gElementInterval good=a.getIntersection(b);
        gString goodseq;
        if (good.getLength()>0) goodseq=gSequence(i_sequence,good.getStart(),good.getEnd());
        ret=getSequenceFromRegion(a.getLeftDiff(b))+goodseq+getSequenceFromRegion(a.getRightDiff(b));
    } else {
        //gReferenceInterval intv=gReferenceInterval(getRegionPositionsFromElement(elementInterval));
        ret=getSequenceFromRegion(b);
    }
    return ret;
}

gSequence gElement::getSequence(gRelativePos startOffset, gRelativePos endOffset) const{
  return getSequence(getElementInterval(),startOffset,endOffset);
}

/** @brief Sequence fragments merge
 *
 * Returns a gSequence object containing the merged sequences from the intervals defined by
 * starts and end postionis
 * @param elementStarts a gArray contining element positions corrsponding to the desired sequence fragments starts
 * @param elementEnds a gArray contining element positions corrsponding to the desired sequence fragments ends
 * @return a gSequence object
 */
gSequence gElement::getSequence(const gArray<gRelativePos> & elementStarts,const gArray<gRelativePos> & elementEnds) const{
  gSequence ret;
  for(gPos i=0;i<elementStarts.getSize();i++){
   ret+=getSequence(gElementInterval(elementStarts[i],elementEnds[i],!elementStarts.isNA(i),!elementEnds.isNA(i)));
  }
  return ret;
}

/** @brief Reference positions from element
 *
 * Returns a gArray contining reference positions mapping to the element positions given
 * @param elementPositions a gArray contining element positions
 * @return a gArray<gPos> object
 */
gArray<gPos> gElement::getReferencePositionsFromElement(const gArray<gRelativePos> & elementPositions) const {
    gArray<gPos> displacements(0,elementPositions.getSize(),false);
    return getReferencePositionsFromElement(elementPositions,displacements);
}

/** @brief Element positions from reference
 *
 * Returns a gArray containing the element positions corresponding to the given reference positions
 * @param referencePositions a gArray<gPos> contining reference positions
 * @return a gArray<gRelativePos> object
 */
gArray<gRelativePos> gElement::getElementPositionsFromReference(const gArray<gPos> & referencePositions) const {
    return getElementPositionsFromReference(referencePositions,gArray<gPos>(0,referencePositions.getSize(),false));
}

/** @brief Element to element positions mapping
 *
 * Returns a gArray<gRelativePos> containing the element positions corresponding to the given element's element positions
 * @param element a gElement to map elementPositions
 * @param elementPositions a gArray<gPos> contining reference positions
 * @return a gArray<gRelativePos> object
 */
gArray<gRelativePos> gElement::mapElementPositions(const gElement & element,const gArray<gRelativePos> & elementPositions) const {
    gArray<gRelativePos> ret(0,elementPositions.getSize(),true);
    if (i_region->getReferenceName()==element.i_region->getReferenceName()) {
        gArray<gPos> displacements(0,elementPositions.getSize(),false);
        gArray<gPos> refPos=getReferencePositionsFromElement(elementPositions,displacements);
        gArray<gPos> actualDisplacements(0,elementPositions.getSize(),false);
        for(gSize i=0;i<refPos.getSize();i++){
          if(displacements[i]>0){
            gArray<gPos> match1=element.i_region->i_varStart.find(refPos[i]);
            gArray<gPos> match2=i_region->i_varStart.find(refPos[i]);
            if(match1.getSize()>0){
              if((element.i_region->i_varType[match1[0]]==gInsertion)&&(element.i_region->i_varSeq[match1[0]]==i_region->i_varSeq[match2[0]])){
               actualDisplacements.setValue(i,displacements[i],false);
              }
            }
          }
        }
        ret=element.getElementPositionsFromReference(refPos,actualDisplacements);
        gArray<gPos> nona=which(actualDisplacements!=0);
        for(gSize i=0;i<nona.getSize();i++) ret.setValue(nona[i],ret[nona[i]],false);
    }
 return ret;
}

/** @brief Reference interval from element
 *
 * Returns a gReferenecInterval contining the reference interval corresponding to the given element interval
 * with optional offsets
 * @param elementInterval a gElementInterval defining the element interval
 * @param startOffset an optional offset for the element start position
 * @param endOffset an optional offset for the element end position
 * @return a gReferenceInterval
 */
gReferenceInterval gElement::getReferenceIntervalFromElement(const gElementInterval & elementInterval,gRelativePos startOffset,gRelativePos endOffset) const {
  gReferenceInterval ret;
  if(elementInterval.validStart()&&elementInterval.validEnd()){
    gArray<gPos> displacements(0,2,false);
    gElementInterval nint(elementInterval.getStart()+startOffset,elementInterval.getEnd()+endOffset,elementInterval.validStart(),elementInterval.validEnd());
    ret=getReferenceIntervalFromElement(nint,displacements);
  }
  return ret;
}

/** @brief Element interval from reference
 *
 * Returns a gElementInterval contining the element interval corresponding to the given reference interval
 * with optional offsets
 * @param referenceInterval a gReferenceInterval defining the element interval
 * @param startOffset an optional offset for the element start position
 * @param endOffset an optional offset for the element end position
 * @return a gReferenceInterval
 */
gElementInterval gElement::getElementIntervalFromReference(const gReferenceInterval & referenceInterval,gPos startOffset,gPos endOffset) const {
    gElementInterval ret;
    if(referenceInterval.validStart() && referenceInterval.validEnd()){
      gArray<gPos> displacements(0,2,false);
      gReferenceInterval nint(referenceInterval.getStart()+startOffset,referenceInterval.getEnd()+endOffset,referenceInterval.validStart(),referenceInterval.validEnd());
      ret=getElementIntervalFromReference(nint,displacements);
    }
    return ret;
}

/** @brief Element to element interval mapping
 *
 * Returns a gElementInterval containing the interval mappied to the element provided
 * @param element a gElement to map elementInterval
 * @param elementInterval a gElementInterval to map
 * @return a gElementInterval object
 */
gElementInterval gElement::mapElementInterval(const gElement & element,const gElementInterval & elementInterval) const {
  gArray<gRelativePos> pos=mapElementPositions(element,elementInterval.getPositions());
  bool sameStrand= (i_forward==element.i_forward);
  if(sameStrand){
    if(pos.isNA(1)) pos.setValue(1,pos[1]+1,true);
  }else{
    if(pos.isNA(0)) pos.setValue(0,pos[0]+0,true);
  }
  return gElementInterval(pos[ getArray<gPos>((sameStrand)?(0):(1),(sameStrand)?(1):(0),1)]);
}

/** @brief Sites insertion from element positions
 *
 * Insert the sites identified by element positions contined in g gArray<gRelativePos>.
 * Double occurring or already inserted positions are not duplicated.
 * Returns an gArray<gArrayIndex> with the indexes of the inserted sites
 * @param elementPositions the sites element positions
 * @return a gArray<gArrayIndex> object
 */
gArray<gArrayIndex> gElement::insertElementSites(const gArray<gRelativePos> & elementPositions) {
    gArray<gRelativePos> tmpelm=elementPositions.getUnique();
    gPos nsites=tmpelm.getSize();
    gArray<gPos> tmpregdisp(0,nsites,true);
    gArray<gPos> a=getRegionPositionsFromElement(tmpelm);
    gArray<gPos> tmpsites=i_region->getReferencePositionsFromRegion(a,tmpregdisp);
    gArray<gArrayIndex> npos=i_site_list.match(tmpsites);
    for (gPos i=0;i<tmpsites.getSize();i++){
      if(npos.isNA(i)){
        npos.setValue(i,i_site_list.getSize(),false);
        i_site_list.setValue(i_site_list.getSize(),tmpsites[i],tmpsites.isNA(i));
        i_site_list_reg_disp.setValue(i_site_list_reg_disp.getSize(),tmpregdisp[i],false);
      }else if(i_site_list_reg_disp[npos[i]]!=tmpregdisp[i]){
        npos.setValue(i,i_site_list.getSize(),false);
        i_site_list.setValue(i_site_list.getSize(),tmpsites[i],tmpsites.isNA(i));
        i_site_list_reg_disp.setValue(i_site_list_reg_disp.getSize(),tmpregdisp[i],false);
      }
    }
    if(elementPositions.getSize()!=nsites){
      tmpregdisp=gArray<gPos>(0,elementPositions.getSize(),true);
      a=getRegionPositionsFromElement(elementPositions);
      tmpsites=i_region->getReferencePositionsFromRegion(a,tmpregdisp);
      npos=i_site_list.match(tmpsites);
    }
    return npos;
}
/** @brief Sites insertion from reference positions
 *
 * Insert the sites identified by reference positions contained in g gArray<gPos>.
 * Double occurring or already inserted positions are not duplicated.
 * Returns an gArray<gArrayIndex> with the indexes of the inserted sites
 * @param referencePositions the sites reference positions
 * @return a gArray<gArrayIndex> object
 */
gArray<gArrayIndex> gElement::insertReferenceSites(const gArray<gPos> & referencePositions) {
    gArray<gPos> tmpref=referencePositions.getUnique();  
    gArray<gArrayIndex> npos=i_site_list.match(tmpref);
    for (gPos i=0;i<tmpref.getSize();i++){
      if(npos.isNA(i)){
        npos.setValue(i,i_site_list.getSize(),false);
        i_site_list.setValue(i_site_list.getSize(),tmpref[i],tmpref.isNA(i));
        i_site_list_reg_disp.setValue(i_site_list_reg_disp.getSize(),0,false);
      }
    }
    npos=i_site_list.match(referencePositions);
    return npos;
}

/** @brief Sites connection from indexs
 *
 * Creates sites connections returning indexes of the connection created
 * Duplicated connections are not re-inserted
 * Returns an gArray<gArrayIndex> with the indexes corresponding to the connections
 * @param startIndexes the start site indexes
 * @param endIndexes the end site indexes
 * @return a gArray<gArrayIndex> object
 */
gArray<gArrayIndex> gElement::connectSites(const gArray<gArrayIndex> & startIndexes,const gArray<gArrayIndex> & endIndexes) {
    if (startIndexes.getSize()!=endIndexes.getSize()) throw gException("gElement::connectSites: siteIndexes1 and siteIndexes2 must have the same number of elements");
    gArray<gArrayIndex> npos(0,startIndexes.getSize(),true,true);
    for (gPos i=0;i<startIndexes.getSize();i++) {

        gArray<gArrayIndex> pos=i_start_site_index.find(startIndexes[i]);
        for(gPos j=0;j<pos.getSize();j++){
          if (endIndexes[i]==i_end_site_index[pos[j]]) {
            npos.setValue(i,pos[j],false);
          }
        }
        if(npos.isNA(i)){
          npos.setValue(i,i_start_site_index.getSize(),false);
          i_start_site_index.setValue(i_start_site_index.getSize(),startIndexes[i],false);
          i_end_site_index.setValue(i_end_site_index.getSize(),endIndexes[i],false);
        }
/*
        gArray<gArrayIndex> pos=i_start_site_index.find(startIndexes[i]);
        if ((pos.getSize()>0) && (endIndexes[i]==i_end_site_index[pos[0]])) {
            npos.setValue(i,pos[0],false);
        } else {
            npos.setValue(i,i_start_site_index.getSize(),false);
            i_start_site_index.setValue(i_start_site_index.getSize(),startIndexes[i],false);
            i_end_site_index.setValue(i_end_site_index.getSize(),endIndexes[i],false);
        }
*/
    }
    return npos;
}

/** @brief Site connection from element Interval
 *
 * Inserts and connects sites corresponding to the boundaries of an element interval
 * Returns the index of the connection created
 * @param elementInterval the element interval from whixch connection should be inserted
 * @param forward the connection orientation realtive to the element
 * @return a gArrayIndex value
 */
gArrayIndex gElement::insertElementConnection(const gElementInterval & elementInterval,bool forward) {
    gArray<gArrayIndex> a=insertElementSites(elementInterval.getPositions());
    if (forward) {
        return connectSites(a[0],a[1])[0];
    } else {
        return connectSites(a[1],a[0])[0];
    }
}

/** @brief Multiple site connections from element positions
 *
 * Inserts sites from start and end element positions arrays and connect them
 * Returns gArray<gArrayIndex> containing the indexes of the connections created
 * @param elementStarts a gArray<gRelativePos> containing starting element positions
 * @param elementEnds a gArray<gRelativePos> containing ending element positions
 * @param forward
 * @return a gArray<gArrayIndex> value
 */
gArray<gArrayIndex> gElement::insertElementConnections(const gArray<gRelativePos> & elementStarts,const gArray<gRelativePos> & elementEnds,bool forward){
 gArray<gArrayIndex> ret(0,elementStarts.getSize(),false,false);
 for(gPos i=0;i<elementStarts.getSize();i++){
   gArrayIndex idx=insertElementConnection(gElementInterval(elementStarts[i],elementEnds[i],!elementStarts.isNA(i),!elementEnds.isNA(i)),forward);
   ret.setValue(i,idx,false);
 }
 return ret;
}

/** @brief Site connection from reference Interval
 *
 * Inserts and connects sites corresponding to the boundaries of a reference interval
 * Returns the index of the connection created
 * @param referenceInterval the reference interval from which connection should be inserted
 * @param forward
 * @return a gArrayIndex value
 */
gArrayIndex gElement::insertReferenceConnection(const gReferenceInterval & referenceInterval,bool forward) {
    gArray<gArrayIndex> a=insertReferenceSites(referenceInterval.getPositions());
    if (forward) {
        return connectSites(a[0],a[1])[0];
    } else {
        return connectSites(a[1],a[0])[0];
    }
}

/** @brief Multiple site connections from reference positions
 *
 * Inserts sites from start and end reference positions arrays and connect them
 * Returns gArray<gArrayIndex> containing the indexes of the connections created
 * @param referenceStarts a gArray<gPos> containing starting reference positions
 * @param referenceEnds a gArray<gPos> containing ending reference positions
 * @param forward
 * @return a gArray<gArrayIndex> value
 */
gArray<gArrayIndex> gElement::insertReferenceConnections(const gArray<gPos> & referenceStarts,const gArray<gPos> & referenceEnds,bool forward){
 gArray<gArrayIndex> ret(0,referenceStarts.getSize(),false,false);
 if (forward) {
   for(gPos i=0;i<referenceStarts.getSize();i++){
     gArrayIndex idx=insertReferenceConnection(gReferenceInterval(referenceStarts[i],referenceEnds[i],!referenceStarts.isNA(i),!referenceEnds.isNA(i)),forward);
     ret.setValue(i,idx,false);
   }
 }else{
   for(gPos i=0;i<referenceStarts.getSize();i++){
     gArrayIndex idx=insertReferenceConnection(gReferenceInterval(referenceStarts[i],referenceEnds[i],!referenceStarts.isNA(i),!referenceEnds.isNA(i)),forward);
     ret.setValue(referenceStarts.getSize()-i-1,idx,false);
   }
 }
 return ret;
}

/** @brief Number of sites
 *
 * returns the number of inserted sites
 * @return a gSize value
 */
gSize gElement::getSitesCount() const {
    return i_site_list.getSize();
}

/** @brief Number of connections
 *
 * returns the number of inserted connections
 * @return a gSize value
 */
gSize gElement::getConnectionsCount() const {
    return i_start_site_index.getSize();
}


/** @brief Get sites indexes from element positions
 *
 * Returns the indexes of the sites corresponding to element positions
 * @param elementPositions a gArray containing the element positions to check for
 * @return a gArray<gArrayindex> object with NA corresponding to non site positions
 */
gArray<gArrayIndex> gElement::getElementSiteIdxs(const gArray<gRelativePos> & elementPositions) const{
  gPos nsites=elementPositions.getSize();
  gArray<gPos> tmpregdisp(0,nsites,true);
  gArray<gPos> a=getRegionPositionsFromElement(elementPositions);
  gArray<gPos> tmpsites=i_region->getReferencePositionsFromRegion(a,tmpregdisp);
  gArray<gArrayIndex> npos=i_site_list.match(tmpsites);
  for(gPos i=0;i<npos.getSize();i++){
    if(!npos.isNA(i)){
      if(tmpregdisp[i]!=i_site_list_reg_disp[npos[i]]) npos.setValue(i,npos[i],true);
    }
  }
  return npos;
}

/** @brief Get sites indexes from reference positions
 *
 * Returns the indexes of the sites corresponding to reference positions
 * @param referencePositions a gArray containing the reference positions to check for
 * @return a gArray<gArrayindex> object with NA corresponding to non site positions
 */
gArray<gArrayIndex> gElement::getReferenceSiteIdxs(const gArray<gPos> & referencePositions) const{
  return i_site_list.match(referencePositions);
}

/** @brief Get the connection index corresponding to an element interval
 *
 * Returns the index of the connection corresponding to an element interval
 * @param elementInterval the element interval to check for
 * @return gArray<gArrayIndex> containin g the index marked as NA if there is no matching connectio
 */
gArray<gArrayIndex> gElement::getElementConnectionIdx(const gElementInterval & elementInterval) const{
 gArray<gArrayIndex> ret(0,1,true,true);
 gArray<gArrayIndex> sid=getElementSiteIdxs(elementInterval.getPositions());
 if(!sid.isNA(0) && !sid.isNA(1)){
   gArray<gArrayIndex> cid=i_start_site_index.find(sid[0]);
   for(gPos i=0;i<cid.getSize();i++){
     if(i_end_site_index[cid[i]]==sid[1]){
       ret.setValue(0,cid[i],false);
       break;
    }
  }
 }
 return ret;
}

/** @brief Get the connection index corresponding to a reference interval
 *
 * Returns the index of the connection corresponding to a reference interval
 * @param referenceInterval the reference interval to check for
 * @return gArray<gArrayIndex> containin g the index marked as NA if there is no matching connectio
 */
gArray<gArrayIndex> gElement::getReferenceConnectionIdx(const gReferenceInterval & referenceInterval,bool forward) const{
 gArray<gArrayIndex> ret(0,1,true,true);
 gArray<gArrayIndex> sid=getReferenceSiteIdxs(referenceInterval.getPositions());
 if(!sid.isNA(0) && !sid.isNA(1)){
   if(forward){
    gArray<gArrayIndex> cid=i_start_site_index.find(sid[0]);
    for(gPos i=0;i<cid.getSize();i++){
      if(i_end_site_index[cid[i]]==sid[1]){
        ret.setValue(0,cid[i],false);
        break;
      }
    }
   }else{
    gArray<gArrayIndex> cid=i_start_site_index.find(sid[1]);
    for(gPos i=0;i<cid.getSize();i++){
      if(i_end_site_index[cid[i]]==sid[0]){
        ret.setValue(0,cid[i],false);
        break;
      }
    }
   }
 }
 return ret;
}

/** @brief Get connections indexes between sites
 *
 * Returns the indexes of the connections connecting startIdxs and endIdxs if any (NA otherwise)
 * @param startIdxs index of the starting site sites
 * @param endIdxs index of the ending site sites
 * @return a gArray<gArrayindex> object
 */
gArray<gArrayIndex> gElement::getConnections(const gArray<gArrayIndex> &startIdxs,const gArray<gArrayIndex> &endIdxs) const{
 gArray<gArrayIndex> ret(0,startIdxs.getSize(),true);
 gArray<gPos> connected=which( !(i_start_site_index.match(startIdxs).isNA() || i_end_site_index.match(endIdxs).isNA()) );
 for(gPos i=0;i<connected.getSize();i++){
   gPos siteNum=connected[i];
   gArray<gPos> k=which((i_start_site_index==startIdxs[siteNum]) && (i_end_site_index==endIdxs[siteNum]));
   if(k.getSize()>0){
     ret.setValue(siteNum,k[0],false);
   }
 }
 return ret;
}

/** @brief Sites connecting to
 *
 * Returns the indexes of the sites connecting to the provided one
 * @param site the index of the site
 * @return a gArray<gArrayindex> object
 */
gArray<gArrayIndex> gElement::connectsTo(gArrayIndex site) const {
    return i_end_site_index[i_start_site_index.find(site)];
}

/** @brief Sites connecting from
 *
 * Returns the indexes of the sites connecting from the provided one
 * @param site the index of the site
 * @return a gArray<gArrayindex> object
 */
gArray<gArrayIndex> gElement::connectsFrom(gArrayIndex site) const {
    return i_start_site_index[i_end_site_index.find(site)];
}

/** @brief Sites element positions
 *
 * Returns a gArray<gRealtivePos> contining the element positions of the sites
 * corresponding to the indexes provided.
 * @param sites a gArray<gArrayindex> containing the indexes of the sites
 * @return a gArray<gRelativePos> object
 */
gArray<gRelativePos> gElement::getSitesElementPosition(const gArray<gArrayIndex> & sites) const {
  gArray<gRelativePos> ret;
  if(sites.getSize()>0){
    ret=getElementPositionsFromReference(i_site_list[sites],i_site_list_reg_disp[sites]);
  }else{
    ret=gArray<gRelativePos>();
  }
  return ret;
}

/** @brief Sites reference positions
 *
 * Returns a gArray<gPos> contining the reference positions of the sites
 * corresponding to the indexes provided.
 * @param sites a gArray<gArrayindex> containing the indexes of the sites
 * @return a gArray<gPos> object
 */
gArray<gPos> gElement::getSitesReferencePosition(const gArray<gArrayIndex> & sites) const{
  if(sites.getSize()>0){
    return i_site_list[sites];
  }else{
    return gArray<gPos>();
  }
}

/** @brief Element interval between two sites
 *
 * Returns a gElementInterval corresponding to the element interval defined by two sites
 * @param startSiteIndex the index of the starting site
 * @param endSiteIndex the index of the ending site
 * @return a gElementInterval object
 */
gElementInterval gElement::getConnectionElementInterval(gArrayIndex startSiteIndex,gArrayIndex endSiteIndex) const {
    gArray<gArrayIndex> sites(0,2,false);
    if (i_site_list[startSiteIndex]<=i_site_list[endSiteIndex]) {
        sites.setValue(0,startSiteIndex,false);
        sites.setValue(1,endSiteIndex,false);
    } else {
        sites.setValue(0,endSiteIndex,false);
        sites.setValue(1,startSiteIndex,false);
    }
    return getElementIntervalFromReference(gReferenceInterval(i_site_list[sites]),i_site_list_reg_disp[sites]);
}

/** @brief Connection strand relative to the element
 *
 * Returns a true if the corresponding connection has the same orientation than the element
 * @param connectionIndex the index of the connection
 * @return a bool value
 */
bool gElement::getConnectionElementStrand(gArrayIndex connectionIndex) const {
  if(i_forward){
   return  getConnectionReferenceStrand(connectionIndex);
  }else{
   return  !getConnectionReferenceStrand(connectionIndex);
  }
}

/** @brief Element interval from connection
 *
 * Returns a gElementInterval corresponding to the element interval defined by a connection
 * @param connectionIndex the index of the connection
 * @return a gElementInterval object
 */
gElementInterval gElement::getConnectionElementInterval(gArrayIndex connectionIndex) const {
    return getConnectionElementInterval(i_start_site_index[connectionIndex],i_end_site_index[connectionIndex]);
}

/** @brief Reference interval between two sites
 *
 * Returns a gREferernceInterval corresponding to the reference interval defined by two sites
 * @param startSiteIndex the index of the starting site
 * @param endSiteIndex the index of the ending site
 * @return a gReferenceInterval object
 */
gReferenceInterval gElement::getConnectionReferenceInterval(gArrayIndex startSiteIndex,gArrayIndex endSiteIndex) const{
    gArray<gArrayIndex> sites(0,2,false);
    if (i_site_list[startSiteIndex]<=i_site_list[endSiteIndex]) {
        sites.setValue(0,startSiteIndex,false);
        sites.setValue(1,endSiteIndex,false);
    } else {
        sites.setValue(0,endSiteIndex,false);
        sites.setValue(1,startSiteIndex,false);
    }
  return gReferenceInterval(i_site_list[sites]);
}


/** @brief Connection strand relative to the reference
 *
 * Returns a true if the corresponding connection has the same orientation than the reference
 * @param connectionIndex the index of the connection
 * @return a bool value
 */
bool gElement::getConnectionReferenceStrand(gArrayIndex connectionIndex) const {
  gPos sp=i_site_list[i_start_site_index[connectionIndex]]+i_site_list_reg_disp[i_start_site_index[connectionIndex]];
  gPos ep=i_site_list[i_end_site_index[connectionIndex]]+i_site_list_reg_disp[i_end_site_index[connectionIndex]];
  return sp<ep;
}

/** @brief Reference interval from connection
 *
 * Returns a gReferenceInterval corresponding to the reference interval defined by a connection
 * @param connectionIndex the index of the connection
 * @return a gReferenceInterval object
 */
gReferenceInterval gElement::getConnectionReferenceInterval(gArrayIndex connectionIndex) const{
  return getConnectionReferenceInterval(i_start_site_index[connectionIndex],i_end_site_index[connectionIndex]);
}

/** @brief Connection start sites
 *
 * Returns the start site indexes corrsponding to the connections identified by an array of indexes
 * @param connectionIndexes the connection indexes
 * @return a gArray<gArrayIndex> object
 */
gArray<gArrayIndex> gElement::getConnectionsStartSiteIndex(const gArray<gArrayIndex> & connectionIndexes) const{
  if(connectionIndexes.getSize()==0){
    return gArray<gArrayIndex>();
  }else{
    return i_start_site_index[connectionIndexes];
  }
}

/** @brief Connection end sites
 *
 * Returns the end site indexes corrsponding to the connections identified by an array of indexes
 * @param connectionIndexes the connection indexes
 * @return a gArray<gArrayIndex> object
 */
gArray<gArrayIndex> gElement::getConnectionsEndSiteIndex(const gArray<gArrayIndex> & connectionIndexes) const{
  if(connectionIndexes.getSize()==0){
    return gArray<gArrayIndex>();
  }else{
    return i_end_site_index[connectionIndexes];
  }
}

/** @brief Connections mapping
 *
 * Maps passed connections indexes from the element called upon to those of the passed one.
 * Returns the corresponding index if the connection do exists in the passed element; NA otherwise:
 * @param element the element to which map connections
 * @param connectionIdxs array containing the connections Idxs to map
 * @return a gArray<gArrayIndex> object
 */
gArray<gArrayIndex> gElement::mapElementConnections(const gElement & element, const gArray<gArrayIndex> & connectionIdxs){
  gArray<gArrayIndex> ret(0,connectionIdxs.getSize(),true);
  for(gPos i=0;i<connectionIdxs.getSize();i++){
    //gArray<gArrayIndex> anIdx=element.getReferenceConnectionIdx(getConnectionReferenceInterval(connectionIdxs[i]));
    gElementInterval eintv=getConnectionElementInterval(connectionIdxs[i]);
    gElementInterval owr=element.mapElementInterval(element,eintv);
    gArray<gArrayIndex> anIdx=element.getElementConnectionIdx(owr);
    ret.setValue(i,anIdx[0],anIdx.isNA(0));
  }
  return ret;
}


/** @brief Element from interval
 *
 * Returns a new gElement corresponding to the passed element interval and strand
 * with optional start and end offsets. Offsets are applied to the interval then the
 * strand is considered
 *
 * @param elementInterval the Interval
 * @param forward the strand
 * @param startOffset start offset relative to the interval (defaults to 0)
 * @param endOffset end offset (defaults to 0)
 * @return a gElement object
 */
gElement gElement::getElement(const gElementInterval & elementInterval,bool forward, gRelativePos startOffset, gRelativePos endOffset) const {
    //gElementInterval intv(elementInterval.getStart()+startOffset,elementInterval.getEnd()+endOffset,elementInterval.validStart()||(startOffset!=0),elementInterval.validEnd()||(endOffset!=0));
    gElementInterval intv(elementInterval.getStart()+startOffset,elementInterval.getEnd()+endOffset);
    return gElement(*this,intv,forward);
}

/** @brief Element from sites
 *
 * Returns a new gElement corresponding to the interval defined by two sites
 * with optional start and end offsets
 * @param startSiteIndex the starting site index
 * @param endSiteIndex the ending site index
 * @param startOffset oprion start offset (defaults to 0)
 * @param endOffset oprion end offset (defaults to 0)
 * @return a gElement object
 */
gElement gElement::getElement(gArrayIndex startSiteIndex, gArrayIndex endSiteIndex, gRelativePos startOffset, gRelativePos endOffset) const {
    gElementInterval interval=getConnectionElementInterval(startSiteIndex,endSiteIndex);
    if (i_site_list[startSiteIndex]<=i_site_list[endSiteIndex]) {
        return getElement(interval,i_forward,startOffset,endOffset);
    } else {
        return getElement(interval,!i_forward,startOffset,endOffset);
    }
}

/** @brief Element from connection
 *
 * Returns a new gElement corresponding to the interval defined by a connection
 * with optional start and end offsets. Offsets are relative to the connection
 * @param connection the connection
 * @param startOffset oprion start offset (defaults to 0)
 * @param endOffset oprion end offset (defaults to 0)
 * @return a gElement object
 */
gElement gElement::getElement(gArrayIndex connection, gRelativePos startOffset, gRelativePos endOffset) const {
  gElementInterval tintv=getConnectionElementInterval(connection);
  gElementInterval intv;
  bool fw=getConnectionElementStrand(connection);
  if(fw){
    intv=gElementInterval(tintv.getStart()+startOffset,tintv.getEnd()+endOffset);
  }else{
    intv=gElementInterval(tintv.getStart()+endOffset,tintv.getEnd()+startOffset);
  }
  return getElement(intv,fw);
}

//protected & private:
gRegion * gElement::getRegion() const {
 return i_region;  
}

gArray<gPos> gElement::getRegionPositionsFromElement(const gArray<gRelativePos> & elementPositions) const {
    //gReferenceInterval iv=gReferenceInterval(i_region->getRegionPositionsFromReference(i_boundaries.getPositions())+i_boundaries_reg_disp);
    gArray<gPos> bpos=i_boundaries.getPositions();
    bpos.setNoNa();
    gArray<gPos> ipos=i_region->getRegionPositionsFromReference(bpos)+i_boundaries_reg_disp;
    if (ipos.isNA(0)) ipos.setValue(0,ipos[0]+1,true);
    gReferenceInterval iv=gReferenceInterval(ipos);
    gArray<gRelativePos> tmp;
    if (i_forward) {
        tmp=gArray<gRelativePos>(iv.getStart(),elementPositions.getSize(),false)+elementPositions;
    } else {
        tmp=gArray<gRelativePos>(iv.getEnd()-1,elementPositions.getSize(),false)-elementPositions;
    }
    if ( (tmp<0).getSum()[0]>0 ) throw gException("gElement::getRegionPositionsFromElement: invalid positions");
    return (gArray<gPos>) tmp;
}

gArray<gRelativePos> gElement::getElementPositionsFromRegion(const gArray<gPos> & regionPositions) const {
    gArray<gPos> bpos=i_boundaries.getPositions();
    bpos.setNoNa();
    gArray<gPos> ipos=i_region->getRegionPositionsFromReference(bpos)+i_boundaries_reg_disp;
    if (ipos.isNA(0)) ipos.setValue(0,ipos[0]+1,true);
    gReferenceInterval iv=gReferenceInterval(ipos);
    gArray<gRelativePos> tmp;
    if (i_forward) {
        tmp=(gArray<gRelativePos>) regionPositions - gArray<gRelativePos>(iv.getStart(),regionPositions.getSize(),false);
    } else {
        tmp=gArray<gRelativePos>(iv.getEnd()-1,regionPositions.getSize(),false) - (gArray<gRelativePos>) regionPositions;
        for(gSize i=0;i<tmp.getSize();i++) if(tmp.isNA(i)) tmp.setValue(i,tmp[i]-1,true);
    }
    return tmp;
}

gArray<gPos> gElement::getReferencePositionsFromElement(const gArray<gRelativePos> & elementPositions,gArray<gPos> & displacements) const {
    gArray<gPos> a=getRegionPositionsFromElement(elementPositions);
    return i_region->getReferencePositionsFromRegion(a,displacements);
}

gReferenceInterval gElement::getReferenceIntervalFromElement(const gElementInterval & elementInterval,gArray<gPos> & displacements) const {
    gArray<gPos> a=getRegionPositionsFromElement(elementInterval.getPositions());
    gArray<gPos> b=i_region->getReferencePositionsFromRegion(a[getArray<gPos>((i_forward)?(0):(1),(i_forward)?(1):(0),1)],displacements);
    if(i_forward){
      if (b.isNA(0)) b.setValue(0,b[0]+1,true);
    }else{
      if (b.isNA(0)) b.setValue(0,b[0]+1,true);
    }
    return gReferenceInterval(b);
}

gArray<gRelativePos> gElement::getElementPositionsFromReference(const gArray<gPos> & referencePositions,const gArray<gPos> & displacements) const {
    gArray<gPos> positions=referencePositions;
    gArray<gPos> a=i_region->getRegionPositionsFromReference(positions)+displacements;
    return getElementPositionsFromRegion(a);
}

gElementInterval gElement::getElementIntervalFromReference(const gReferenceInterval & referenceInterval,const gArray<gPos> & displacements) const {
    gArray<gPos> positions=referenceInterval.getPositions();
    positions.setNoNa();
    gArray<gPos> a=i_region->getRegionPositionsFromReference(positions)+displacements;
    if(i_forward){
      if (a.isNA(0)) a.setValue(0,a[0]+1,true);
    }else{
      if (a.isNA(1)) a.setValue(1,a[1]-1,true);
    }
    gArray<gRelativePos> tmp=getElementPositionsFromRegion(a[getArray<gPos>((i_forward)?(0):(1),(i_forward)?(1):(0),1)]);

    return gElementInterval(tmp);
}

gSequence gElement::getSequenceFromRegion(const gElementInterval & interval) const {
    gSequence ret;
    if (interval.getLength()>0) {
      gArray<gPos> positions=getRegionPositionsFromElement(interval.getPositions());
      gReferenceInterval intv(positions[getArray<gPos>((i_forward)?(0):(1),(i_forward)?(1):(0),1)]);
      if (intv.getLength()>0) {
        if (i_forward) {
          ret=i_region->getSequence(intv);
        } else {
          ret=i_region->getSequence(intv).getReverseComplement();
        }
      }
    }
    return ret;
}

gMatrix<gRelativePos> gElement::getVariationElementIntervals() const {
  gMatrix<gRelativePos> ret;
  throw gException("gElement::getVariationElementIntervals: not yet implemented");
  return ret;
}

void gElement::setReferenceSites(gArray<gPos> referencePositions){
  i_site_list=referencePositions;
  i_site_list_reg_disp=gArray<gPos>(0,referencePositions.getSize(),false);
}

void gElement::setElementSites(gArray<gRelativePos> elementPositions){
  gArray<gPos> tmpregdisp(0,elementPositions.getSize(),true);
  gArray<gPos> a=getRegionPositionsFromElement(elementPositions);
  i_site_list=i_region->getReferencePositionsFromRegion(a,i_site_list_reg_disp);
}

void gElement::setConnectionIdxs(const gArray<gArrayIndex> &startIdxs,const gArray<gArrayIndex> &endIdxs){
 if(startIdxs.getSize()!=endIdxs.getSize()) throw gException("gElement::setConnectionIdxs: startIdxs and endIdxs must have equal length");
 i_start_site_index=startIdxs;
 i_end_site_index=endIdxs;
}
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
/// gRegionFeaturesParameters implementation
//--------------------------------------------------------------------------
gRegionFeaturesParameters::gRegionFeaturesParameters(){
}

gRegionFeaturesParameters::gRegionFeaturesParameters(const gArray<gPos> & sites,const std::vector< gArray<gRelativePos> > & siteOffsets,overlapControl control,gRelativePos shift):i_sites(sites),i_siteOffsets(siteOffsets),i_control(control),i_shift(shift){
  if(sites.getSize()!=siteOffsets.size()) throw gException("gRegionsFeatureExtractor: incompatible sites and offsets arrays");
  i_nwindows=0;
  for(gSize i=0;i<i_sites.getSize();i++){
    i_nwindows+=i_siteOffsets[i].getSize();
  }
}

gRegionFeaturesParameters::gRegionFeaturesParameters(const gRegionFeaturesParameters & parameters):i_sites(parameters.i_sites),i_siteOffsets(parameters.i_siteOffsets),i_control(parameters.i_control),i_shift(parameters.i_shift),i_nwindows(parameters.i_nwindows){
}

bool gRegionFeaturesParameters::operator == (const gRegionFeaturesParameters & parameters) const{
  bool ret=i_sites.equals(parameters.i_sites) && i_control==parameters.i_control && i_shift==parameters.i_shift && i_siteOffsets.size()==parameters.i_siteOffsets.size();
  for(gSize i=0;ret && (i<i_siteOffsets.size());i++){
    ret=ret && i_siteOffsets[i].equals(parameters.i_siteOffsets[i]);
  }
  return ret;
}

gArray<gPos> gRegionFeaturesParameters::getSites() const{
  return i_sites;
}

std::vector< gArray<gRelativePos> > gRegionFeaturesParameters::getOffsets() const{
  return i_siteOffsets;
}

gArray<gRelativePos> gRegionFeaturesParameters::getSiteOffsets(gPos siteNum) const{
  return i_siteOffsets[siteNum];
}

gRegionFeaturesParameters::overlapControl gRegionFeaturesParameters::getControl() const{
  return i_control;
}

gSize gRegionFeaturesParameters::getWindowsCount() const{
  return i_nwindows;
}
//--------------------------------------------------------------------------



//--------------------------------------------------------------------------
/// gElementCollection implementation
//--------------------------------------------------------------------------
gElementCollection::gElementCollection(){
}

void gElementCollection::insert(const gElement & e,gScore score, const gString & name){
  i_elements.push_back(e);
  i_scores.setValue(i_scores.getSize(),score,false);
  i_names.push_back(name);
}

gElementCollection & gElementCollection::append(const gElementCollection & elements){
  i_scores.concatenate(elements.i_scores);
  gSize cs=i_elements.size();
  i_elements.resize(cs+elements.i_elements.size());
  i_names.resize(cs+elements.i_names.size());
  for(gSize i=0;i<elements.i_elements.size();i++){
   i_elements[cs+i]=elements.i_elements[i];
   i_names[cs+i]=elements.i_names[i];
  }
  return *this;
}
 
gSize gElementCollection::getElementCount() const{
 return i_elements.size();
}

gArray<gScore> gElementCollection::getScores() const{
  return i_scores;
}

gScore gElementCollection::getScore(gPos id) const{
  return i_scores[id];
}

vector<gString> gElementCollection::getNames() const{
  return i_names;
}

gString gElementCollection::getName(gPos id) const{
  return i_names[id];
}

const gElement & gElementCollection::getElement(gPos id) const{
  return i_elements[id];
}

gElement & gElementCollection::getElement(gPos id){
  return i_elements[id];
}

gArray<gPos> gElementCollection::intersect(const gString & chrom,const gReferenceInterval & interval,gSize minOverlap){
  gBool * selector= new gBool[i_elements.size()];
  memset(selector,0,sizeof(gBool)*i_elements.size());
#pragma omp parallel for  
  for(gSize i=0;i<i_elements.size();i++){
    gElement & e = i_elements[i];
    if(e.getReferenceName()==chrom){
      gSize intersection=interval.getIntersection(e.getReferenceInterval()).getLength();
      if(minOverlap==0){
	selector[i]=(intersection==interval.getLength());
      }else{
	selector[i]=(intersection>=minOverlap);
      }
    }
  }
  gArray<gBool> sele2;
  sele2.ownData(i_elements.size(),selector);
  return which(sele2);
}

gElementCollection gElementCollection::select(const gArray<gPos> & positions) const{
  gElementCollection ret;
  gSize num=positions.getSize();
  ret.i_scores=i_scores[positions];  
  ret.i_elements=vector<gElement>(positions.getSize());
  ret.i_names=vector<gString>(positions.getSize());
  for(gSize i=0;i<num;i++){
    ret.i_elements[i]=i_elements[positions[i]];
    ret.i_names[i]=i_names[positions[i]];
  }
  return ret;
}


gElementCollection gElementCollection::selectN(gSize N) const{
  gElementCollection ret;
  gSize num=min(N,i_elements.size());
  ret.i_scores=gArray<gScore>(i_scores,0,num);
  ret.i_elements=i_elements;
  ret.i_elements.resize(num);
  ret.i_names=i_names;
  ret.i_names.resize(num);
  return ret;
}

gElementCollection gElementCollection::selectRange(gSize first,gSize last) const{
  gElementCollection ret;
  gSize num=last-first;
  ret.i_elements=vector<gElement>(num);
  ret.i_names=vector<gString>(num);
  ret.i_scores=gArray<gScore>(i_scores,first,last);
  for(gSize i=first;i<last;i++){
    ret.i_elements[i-first]=i_elements[i];
    ret.i_names[i-first]=i_names[i];
  }
  return ret;
}

gElementCollection gElementCollection::selectNRandom(gSize N) const{
  gElementCollection ret;
  gSize num=min(N,i_elements.size());
  gArray<gPos> order=randomSample(num,0,i_elements.size()-1);
  ret.i_elements=vector<gElement>(num);
  ret.i_names=vector<gString>(num);
  ret.i_scores=i_scores[order];
  for(gSize i=0;i<num;i++){
    ret.i_elements[i]=i_elements[order[i]];
    ret.i_names[i]=i_names[order[i]];
  }
  return ret;
}

gElementCollection gElementCollection::selectNBest(gSize N) const{
  gElementCollection ret;
  gSize num=min(N,i_elements.size());
  ret.i_scores=i_scores;
  gArray<gPos> order=ret.i_scores.sort(true);
  order.revert();
  ret.i_elements=vector<gElement>(num);
  ret.i_names=vector<gString>(num);
  ret.i_scores=i_scores[order];
  for(gSize i=0;i<num;i++){
    ret.i_elements[i]=i_elements[order[i]];
    ret.i_names[i]=i_names[order[i]];
  }
  return ret;
}

void gElementCollection::setScore(gPos i, gScore score){
  if(i<i_elements.size()){
    i_scores.setValue(i,score,false);
  }else throw gException("gElementCollection::setScore invalid index");
}

void gElementCollection::setName(gPos i, const gString & name){
  if(i<i_elements.size()){
    i_names[i]=name;
  }else throw gException("gElementCollection::setName invalid index");
}

void gElementCollection::setElement(gPos i,const gElement & e){
  if(i<i_elements.size()){
    i_elements[i]=e;
  }else throw gException("gElementCollection::setElement invalid index");
}


void gElementCollection::setScores(const gArray<gScore> & scores){
  if(scores.getSize()==i_scores.getSize()){
    i_scores=scores;
  }else throw gException("gElementCollection::setScores invalid number of scores");
  
}

void gElementCollection::setNames(const std::vector<gString> & names){
  if(names.size()==i_names.size()){
    i_names=names;
  }else throw gException("gElementCollection::setNames invalid number of names");
}
//--------------------------------------------------------------------------



//BEGIN CLASS gElementRegion definition
gElementRegion::gElementRegion(){
  i_connection=0;
  i_mode=0;
  i_offset=0;
  i_halfWindow=0;
}
      
gElementRegion::gElementRegion(gArrayIndex connection,gShortUnsigned mode,gRelativePos offset,gSize halfWindow){
  i_connection=connection;
  i_mode=mode;
  i_offset=offset;
  i_halfWindow=halfWindow;
}
      
gElementInterval gElementRegion::getInterval(const gElement & e) const{
  gElementInterval ret;
  gElementInterval intv=e.getConnectionElementInterval(i_connection);
  gRelativePos pos;
  switch(i_mode){
    case 0: 
      ret=intv;
      break;
    case 1:
      pos=intv.getStart();
      ret=gElementInterval(pos+i_offset-i_halfWindow,pos+i_offset+i_halfWindow);
      break;
    case 2:
      pos=intv.getEnd()-1;
      ret=gElementInterval(pos+i_offset-i_halfWindow,pos+i_offset+i_halfWindow);
      break;
    case 3:
      pos=(intv.getStart()+intv.getEnd())/2;
      ret=gElementInterval(pos+i_offset-i_halfWindow,pos+i_offset+i_halfWindow);
  }
  return ret;
}

gSequence gElementRegion::getSequence(const gElement & e) const{
  return e.getSequence(getInterval(e));
}
//END CLASS gElementRegion definition
