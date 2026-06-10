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
#include "sequence.h"

using namespace std;
using namespace geco;

//--------------------------------------------------------------------------
// static data
//--------------------------------------------------------------------------
static const gChar complement_map_data[256]={
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 45, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 84, 86, 71, 72, 78, 78, 67, 68, 78, 78, 77, 78, 75, 78, 78,
    78, 78, 89, 83, 65, 65, 66, 87, 78, 82, 78, 78, 78, 78, 78, 78,
    78,116,118,103,104,110,110, 99,100,110,110,109,110,107,110,110,
    110,110,121,115, 97, 97, 98,119,110,114,110, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78
};

static const gArray<gChar> complement_map ( complement_map_data,256 );
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// gSequence implementation
//--------------------------------------------------------------------------
/** @brief Empty contructor
*
* Instatiate an empty gSequence object
*/
gSequence::gSequence() :gString() {
}

/** @brief Copy (range) contructor
*
* Instatiate a new gSequence initialized with the values of a
* gArray<gChar> in the range specified.
* @param str the string to copy from
* @param start the beginning position in the provided string (0 based, defaults to 0)
* @param end the ending position in the provided string (1 based, defaults to the last character position)
*/
gSequence::gSequence ( const gArray<gChar> & str,gPos start,gPos end ) :gString ( str,start,end ) {
}

/** @brief c string range contructor
*
* Instatiate a new gSequence initialized with the values of a
* c string in the range specified.
* @param str the string to copy from
* @param start the beginning position in the provided string (0 based, defaults to 0)
* @param end the ending position in the provided string (1 based, defaults to the last character position)
*/
gSequence::gSequence ( const gChar * str,gPos start,gPos end) :gString ( str,start,end ) {
}

/** @brief c string range contructor
   *
   * Instatiate a new gSequence initialized with the values of a
   * std::string in the range specified.
   * @param str the string to copy from
   * @param start the beginning position in the provided string (0 based, defaults to 0)
   * @param end the ending position in the provided string (1 based, defaults to the last character position)
   */
gSequence::gSequence ( const std::string & str,gPos start,gPos end) :gString ( str,start,end ) {
}

/** @brief Given length contructor
 * Instatiate a gSequence object with a given length and initialied to the provided charatcer
 * @param length length of the sequence
 * @param fillvalue character used to initialize the sequence (defaults to 'N')
 */
gSequence::gSequence(gSize length, gChar fillvalue):gString(length,fillvalue){
}

/** @brief Reverese complement
 *
 * Reverese complement this sequence
 * @return a reference to this
 */
gSequence & gSequence::reverseComplement() {
  if(this->getLength()>0){
    gArray<gPos> napos=which(isNA());
    setNoNa();
    gSequence ret=complement_map[*this];
    ret.setNA(napos);
    *this=ret.getReverted();
  }
  return *this;
}

/** @brief Reverese complement
 *
 * Return a sequence equal to the reverse complement of this
 * @return a gSequence object
 */
gSequence gSequence::getReverseComplement() const{
    gSequence ret=*this;
    ret.reverseComplement();
    return ret;
}

/** @brief Sequence translation
 *
 * Translate the current sequence using the specified frame.
 * @param frame The reading frame relative to the beginning of the sequence
 * @return a gSequence object
 */
gString gSequence::getTranslated(gShortUnsigned frame){
 const string raw="TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG";
 const string aa="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
 gString seq(*this,frame);
 seq.upper();
 gSize slen=seq.getLength()/3;
 gString sout(slen,'X');
 for(gSize i=0; i < slen; i++){
  gString codon(seq,3*i,3*i+3);
  //codon.upper();
  string scod=(string)codon;
  size_t pos=raw.find(scod);
  if(pos!=string::npos){
   sout.setValue(i,aa[pos/4],false);
  }
 }
 return sout;
}


gString gSequence::getMtTranslated(gShortUnsigned frame){
 const string raw="TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG";
 const string aa="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
 gString seq(*this,frame);
 seq.upper();
 gSize slen=seq.getLength()/3;
 gString sout(slen,'X');
 for(gSize i=0; i < slen; i++){
  gString codon(seq,3*i,3*i+3);
  //codon.upper();
  string scod=(string)codon;
  size_t pos=raw.find(scod);
  if(pos!=string::npos){
   sout.setValue(i,aa[pos/4],false);
  }
 }
 return sout;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/// gSequenceGroup implementation
//--------------------------------------------------------------------------

//     std::vector<gString> i_sequenceNames;
//     std::vector<gString> i_featureNames;
//     std::vector<gSequence> i_sequences;
//     std::vector< vector< gMatrix<gDoubleScore> > > i_features;
//     gSize i_referenceId;
//     bool i_aligned;
//
//   public:
//     /** @name Constructors */
//     //@{
//     gSequenceGroup(bool aligned=false);
//     gSequenceGroup(const gSequenceGroup & group);
//     //@}
//
//     /** @name Sequence add functions*/
//     void addSequence(const gString & name,const gSequence & sequence,bool revert=false);
//     void addRefSequence(const gString & name,const gSequence & seq,bool revert=false);
//     void addFeature(const seqNum, const gString & name, const gMatrix<gDoubleScore> & feature);
//     //@}
//
//     /** @name Sequence info functions*/
//     gSize getSequenceCount() const;
//     gSize getFeaturesCount(gPos seqNum) const;
//     gString getSequenceName(gPos seqNum) const;
//     gSequence getSequence(gPos seqNum) const;
//     gMatrix<gDoubleScore> getFeature(gPos seqNum,gPos featNum) const;
//     gString getFeatureName(gPos seqNum,gPos featNum) const;
//     gSize getReferenceId() const;
//     bool isAligned() const;


gSequenceGroup::gSequenceGroup(bool aligned){
  i_aligned=aligned;
}

gSequenceGroup::gSequenceGroup(const gSequenceGroup & group){
  i_sequences=group.i_sequences;
  i_sequenceNames=group.i_sequenceNames;
  i_features=group.i_features;
  i_featureNames=group.i_featureNames;
  i_aligned=i_aligned;
}

void gSequenceGroup::addSequence(const gString & name,const gSequence & sequence,bool revert){
  i_sequenceNames.push_back(name);
  i_features.push_back(vector<gMatrix<gDoubleScore> >());
  i_featureNames.push_back(vector<gString>());
  if(revert){
    i_sequences.push_back(sequence.getReverseComplement());
  }else{
    i_sequences.push_back(sequence);
  }
}

void gSequenceGroup::setSequence(gSize seqNum, const gString & name,const gSequence & sequence,bool revert){
  if(i_sequences.size()-1<seqNum){
    i_sequences.resize(seqNum+1);
    i_sequenceNames.resize(seqNum+1);
    i_features.resize(seqNum+1);
    i_featureNames.resize(seqNum+1);
  }
  if(revert){
    i_sequences[seqNum]=sequence.getReverseComplement();
  }else{
    i_sequences[seqNum]=sequence;
  }
  i_sequenceNames[seqNum]=name;
}

void gSequenceGroup::addRefSequence(const gString & name,const gSequence & sequence,bool revert){
  this->addSequence(name,sequence,revert);
  i_referenceId=i_sequences.size()-1;
}

void gSequenceGroup::addFeature(const gSize seqNum, const gString & name, const gMatrix<gDoubleScore> & feature){
  i_features[seqNum].push_back(feature);
  i_featureNames[seqNum].push_back(name);
}

gSize gSequenceGroup::getSequenceCount() const{
 return i_sequences.size();
}

gSequence gSequenceGroup::getSequence(gPos seqNum) const{
  return i_sequences[seqNum];
}

gString gSequenceGroup::getSequenceName(gPos seqNum) const{
  return i_sequenceNames[seqNum];
}

gSize gSequenceGroup::getFeaturesCount(gPos seqNum) const{
  return i_features[seqNum].size();
}


gMatrix<gDoubleScore> gSequenceGroup::getFeature(gPos seqNum,gPos featNum) const{
  return i_features[seqNum][featNum];
}

gString gSequenceGroup::getFeatureName(gPos seqNum,gPos featNum) const{
  return i_featureNames[seqNum][featNum];
}

gSize gSequenceGroup::getReferenceId() const{
  return i_referenceId;
}

bool gSequenceGroup::isAligned() const{
  return i_aligned;
}
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
/// gSequenceRetrieverImplementation implementation
//--------------------------------------------------------------------------
/** @brief Empty constructor
 *
 * You must call this constructor in derived classes
 */
gSequenceRetrieverImplementation::gSequenceRetrieverImplementation():gArrayRetrieverImplementation<gChar>(1){
}

gSequenceRetrieverImplementation::~gSequenceRetrieverImplementation(){
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/// gSequenceRetriever implementation
//--------------------------------------------------------------------------
/** @brief Constructor
 *
 * You must call this constructor in derived classes
 */
gSequenceRetriever::gSequenceRetriever():gArrayRetriever<gChar>(){
}

/** @brief Constructor
 * Contructor from...
 * @param implementation: a reference to an oject of type gSequenceRetrieverImplementation,
 */
gSequenceRetriever::gSequenceRetriever(const gSequenceRetrieverImplementation & implementation):gArrayRetriever<gChar>(implementation){
}

/** @brief Virtual destructor
 *
 * virtual destructor (needed since gArrayRetreiver declares it as pure virtual)
 */
gSequenceRetriever::~gSequenceRetriever(){
}

/**
 * @brief Sequence retriever function member
 *
 * Call this function when you want to retreive a sequence. It will,
 * in turn call the getSequence_internal method and control for consistency of the obtained string.
 * @param reference a reference to the sequence you want to retrieve
 * @param start start position
 * @param end end position
 * @return a gSequence object containing the required sequence;
 */
gSequence gSequenceRetriever::getSequence ( const gString & reference,gPos start,gPos end ) const {
  gSequence ret;
  if(start>end) throw gException("gSequenceRetriever: start>end");
  if(start<end){
    try{
     ret=((const gSequenceRetrieverImplementation &) getImplementation()).getSequence_Internal( reference,start,end );
     gString tmpret=ret;
     tmpret.upper();
     gArray<gPos> inv=which(((gArray<gChar>)tmpret)=='N');
     if(inv.getSize()>0){
      ret.setValues(inv,gArray<gChar>('N',inv.getSize(),true));
     }
    }catch(gException & e){
      throw e;
    }
    if ( ret.getLength() !=end-start ) throw gException ( "gSequenceRetriever: incompatible string" );
  }
  return ret;
}
//--------------------------------------------------------------------------

