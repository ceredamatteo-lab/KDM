/***************************************************************************
 *   Copyright (C) 2010 by Uberto Pozzoli and Federico Marini              *
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
#include "genotype.h"
#include <iostream>
#include <algorithm>


using namespace std;
using namespace geco;

//------------------------------------------------------------------------------
// gGenotypeInterval implementation
//------------------------------------------------------------------------------
gGenotypeInterval::gGenotypeInterval():gReferenceInterval(0,0,false,false),gString(""){
}

gGenotypeInterval::gGenotypeInterval(const gString & name,gPos start,gPos end,bool validstart,bool validend):gReferenceInterval(start,end,validstart,validend),gString(name){
}

gGenotypeInterval::gGenotypeInterval(const gString & name,const gReferenceInterval & interval):gReferenceInterval(interval),gString(name){
}

gGenotypeInterval::gGenotypeInterval(const gString & name,const gArray<gPos> & positions):gReferenceInterval(positions),gString(name){
}

gGenotypeInterval::gGenotypeInterval(const gGenotypeInterval & interval):gReferenceInterval(interval),gString(interval){
}

bool gGenotypeInterval::isWholeReference() const{
  return !validStart() && !validEnd();
}

gArray<bool> gGenotypeInterval::contains(const gArray<gPos> & positions) const{
  gArray<bool> ret;
  if(validStart() && validEnd()){
    return ((gReferenceInterval *)this)->contains(positions);
  }else{
    return gArray<bool>(positions==0 || positions!=0);
  }
}

bool gGenotypeInterval::operator == (const gGenotypeInterval & interval) const{
  return ( ((gString)*this)==((gString) interval) ) && ( ((gReferenceInterval)*this)==((gReferenceInterval) interval) );
}

bool gGenotypeInterval::operator > (const gGenotypeInterval & interval) const{
  bool ret=false;
  if(((gString)*this) == ((gString) interval) ){
    if(!isWholeReference() and !interval.isWholeReference()){
      if( !(((gReferenceInterval)*this)==((gReferenceInterval) interval)) ){
        ret=((gReferenceInterval)*this)>((gReferenceInterval) interval);
      }else{
        ret = ( (gReferenceInterval) *this).getLength() < ((gReferenceInterval)interval).getLength();
      }
    }else{
      ret=!isWholeReference() && interval.isWholeReference();
    }
  }else{
    ret=!( (((gString)*this) < ((gString) interval)) || (((gString)*this) == ((gString) interval)) );
  }
  return ret;
}

bool gGenotypeInterval::operator < (const gGenotypeInterval & interval) const{
  return !((*this>interval)||(*this==interval));
}
//------------------------------------------------------------------------------

class idInterval:public gGenotypeInterval{
public:
  gPos i_id;
  gPos i_posInList;

  idInterval(gPos id,gPos posInList,gGenotypeInterval interval):gGenotypeInterval(interval){
    i_id=id;
    i_posInList=posInList;
  }

  idInterval(const idInterval & interval):gGenotypeInterval(interval){
    i_id=interval.i_id;
    i_posInList=interval.i_posInList;
  }

  bool operator < (const idInterval & interval) const{
    return ((gGenotypeInterval) *this) < ((gGenotypeInterval) interval );
  }


};

vector<gGenotypeInterval> mergeGenotypeIntervals(const vector<gGenotypeInterval> & i1,const vector<gGenotypeInterval> & i2,gArray<gPos> & map1,gArray<gPos> & map2){
 vector<gGenotypeInterval> ret;
 vector<idInterval> intervals;

 for(size_t i=0;i<i1.size();i++){
   intervals.push_back(idInterval(1,i,i1[i]));
 }
 for(size_t i=0;i<i2.size();i++){
   intervals.push_back(idInterval(2,i,i2[i]));
 }
 sort(intervals.begin(),intervals.end());
 ret.push_back(intervals[0]);

 size_t currintv=0;
 if(intervals[0].i_id==1) map1.setValue(intervals[0].i_posInList,currintv,false);
 if(intervals[0].i_id==2) map2.setValue(intervals[0].i_posInList,currintv,false);
 size_t testintv=1;

 while(testintv<intervals.size()){
   //cout << "currintv:\t" << ret[currintv] << " " << ret[currintv].getStart() << "-" << ret[currintv].getEnd() << endl;
   //cout << "currtest:\t" << intervals[testintv] << " " << intervals[testintv].getStart() << "-" << intervals[testintv].getEnd() << endl;

   if( ((gString) intervals[testintv])==((gString)ret[currintv])){
     if(ret[currintv].isWholeReference()||intervals[testintv].isWholeReference()){
       ret[currintv]=gGenotypeInterval(ret[currintv],0,0,false,false);
     }else{
       if(intervals[testintv].getStart()<ret[currintv].getEnd()){
         if(intervals[testintv].getEnd()>ret[currintv].getEnd()){
           ret[currintv]=gGenotypeInterval(ret[currintv],ret[currintv].getStart(),intervals[testintv].getEnd(),true,true);
         }
       }else{
         ret.push_back(intervals[testintv]);
         currintv++;
       }
     }
   }else{
      ret.push_back(intervals[testintv]);
      currintv++;
   }
   if(intervals[testintv].i_id==1) map1.setValue(intervals[testintv].i_posInList,currintv,false);
   if(intervals[testintv].i_id==2) map2.setValue(intervals[testintv].i_posInList,currintv,false);
   //cout << "newintv:\t" << ret[currintv] << " " << ret[currintv].getStart() << "-" << ret[currintv].getEnd() << endl << endl;
   testintv++;
 }
 return ret;
}

//------------------------------------------------------------------------------
// gGenotype implementation
//------------------------------------------------------------------------------
gGenotype::gGenotype(gExperimentType type){
  i_experimentType=type;
  i_isForward=true;
  i_hasReferenceAlleles=true;
}

gGenotype::gGenotype(gSize subjectNum,gSize variantNum,gExperimentType type){
  i_experimentType=type;
  i_isForward=true;
  i_hasReferenceAlleles=true;
  i_variantPositions=gArray<gPos>(0,variantNum,true);
  i_variantIntervalIds=gArray<gPos>(0,variantNum,true);
  i_variantAlleles=vector< vector<gSequence> >(variantNum);
  i_subjectNames=vector<gString>(subjectNum);
  i_allele1Ids=gMatrix<gShortUnsigned>(variantNum,subjectNum,0,true);
  i_allele2Ids=gMatrix<gShortUnsigned>(variantNum,subjectNum,0,true);
  i_isPhased=gMatrix<bool>(variantNum,subjectNum,0,true);
}

gGenotype::gGenotype(const vector<gString> & subjectsNames,gSize variantNum,gExperimentType type){
  i_experimentType = type;
  i_isForward = true;
  i_hasReferenceAlleles = true;
  i_variantPositions = gArray<gPos>(0,variantNum,true);
  i_variantIntervalIds=gArray<gPos>(0,variantNum,true);
  i_variantAlleles=vector< vector<gSequence> >(variantNum);
  i_subjectNames = subjectsNames;
  i_allele1Ids=gMatrix<gShortUnsigned>(variantNum,subjectsNames.size(),0,true);
  i_allele2Ids=gMatrix<gShortUnsigned>(variantNum,subjectsNames.size(),0,true);
  i_isPhased=gMatrix<bool>(variantNum,subjectsNames.size(),0,true);
}

gGenotype::gGenotype(const gGenotype & genotype){
  *this=genotype;
}

gGenotype::gGenotype(const gGenotype & genotype,gPos startSubject,gPos endSubject,bool removeInvariant){
  i_name=genotype.i_name;
  i_experimentType=genotype.i_experimentType;
  i_isForward=genotype.i_isForward;

  i_intervals=genotype.i_intervals;
  i_variantPositions=genotype.i_variantPositions;
  i_variantIntervalIds=genotype.i_variantIntervalIds;
  i_variantAlleles=genotype.i_variantAlleles;
  i_subjectNames.assign(genotype.i_subjectNames.begin()+startSubject,genotype.i_subjectNames.begin()+endSubject);
  gArray<gPos> subIds=getArray<gPos>(startSubject,endSubject-1,1);
  i_allele1Ids=genotype.i_allele1Ids.getColumns(subIds);
  i_allele2Ids=genotype.i_allele2Ids.getColumns(subIds);
  i_isPhased=genotype.i_isPhased.getColumns(subIds);
  if(removeInvariant){
    removeInvariantPositions();
  }
}

gGenotype::gGenotype(const gGenotype & genotype,gPos interval, gPos startPosition,gPos endPosition){
  i_name=genotype.i_name;
  i_experimentType=genotype.i_experimentType;
  i_isForward=genotype.i_isForward;
  gGenotypeInterval oldinterval=genotype.getInterval(interval);
  gArray<gPos> goodpos=which( (genotype.i_variantPositions>=startPosition) && (genotype.i_variantPositions<endPosition) && (genotype.i_variantIntervalIds==interval) );
  gPos newstart=(startPosition>=oldinterval.getStart())?(startPosition):(oldinterval.getStart());
  gPos newend=(endPosition<=oldinterval.getEnd())?(endPosition):(oldinterval.getEnd());
  gGenotypeInterval newinterval=gGenotypeInterval(oldinterval,newstart,newend,true,true);
  i_intervals=vector<gGenotypeInterval>(1,newinterval);
  i_subjectNames=genotype.i_subjectNames;
  if(goodpos.getSize()>0){
    i_variantPositions=genotype.i_variantPositions[goodpos];
    i_variantIntervalIds=gArray<gPos>(0,goodpos.getSize(),false);
    for(gPos i=0;i<goodpos.getSize();i++) i_variantAlleles.push_back(genotype.i_variantAlleles[goodpos[i]]);
    i_allele1Ids=genotype.i_allele1Ids.getRows(goodpos);
    i_allele2Ids=genotype.i_allele2Ids.getRows(goodpos);
    i_isPhased=genotype.i_isPhased.getRows(goodpos);
  }
}

gGenotype & gGenotype::operator =(const gGenotype & genotype){
  i_name=genotype.i_name;
  i_experimentType=genotype.i_experimentType;
  i_isForward=genotype.i_isForward;
  i_hasReferenceAlleles=genotype.i_hasReferenceAlleles;
  i_intervals=genotype.i_intervals;
  i_variantPositions=genotype.i_variantPositions;
  i_variantIntervalIds=genotype.i_variantIntervalIds;
  i_variantAlleles=genotype.i_variantAlleles;
  i_subjectNames=genotype.i_subjectNames;
  i_allele1Ids=genotype.i_allele1Ids;
  i_allele2Ids=genotype.i_allele2Ids;
  i_isPhased=genotype.i_isPhased;
  return *this;
}

void gGenotype::clear(){
  i_name=gString("");
  i_experimentType=gResequencing;
  i_isForward=true;
  i_hasReferenceAlleles=true;
  i_intervals.clear();
  i_variantPositions=gArray<gPos>();
  i_variantIntervalIds=gArray<gPos>();
  i_variantAlleles.clear();
  i_subjectNames.clear();
  i_allele1Ids=gMatrix<gShortUnsigned>();
  i_allele2Ids=gMatrix<gShortUnsigned>();
  i_isPhased=gMatrix<bool>();
}

void gGenotype::setName(const gString & name){
  i_name==name;
}

void gGenotype::setIntervals(const vector<gGenotypeInterval> intervals){
  i_intervals=intervals;
}

void gGenotype::setSubjects(const vector<gString> & subjectNames){
  if(subjectNames.size()==i_subjectNames.size()){
    i_subjectNames=subjectNames;
  }else throw gException("gGenotype::setSubject: incorrect number of subjects");
}

void gGenotype::setVariantPositions(const gArray<gPos> & intervalIds,const gArray<gPos> & positions,const vector< vector<gSequence> > & alleles){
  if(
     (positions.getSize()==i_variantPositions.getSize()) &&
     (intervalIds.getSize()==i_variantPositions.getSize()) &&
     (alleles.size()==i_variantPositions.getSize())
  ){
    i_variantAlleles=alleles;
    i_variantPositions=positions;
    i_variantIntervalIds=intervalIds;
  }else throw gException("gGenotype::setVariantPositions: incorrect number of variants");
}

void gGenotype::setGenotypes(const gMatrix<gShortUnsigned> & allele1Ids,const gMatrix<gShortUnsigned> & allele2Ids,const gMatrix<bool> & isPhased){
  i_allele1Ids=allele1Ids;
  i_allele2Ids=allele2Ids;
  i_isPhased=isPhased;
}

void gGenotype::setVariantGenotypes(gPos variantId, const gArray<gShortUnsigned> & allele1Id, const gArray<gShortUnsigned> & allele2Id,const gArray<bool> & isPhased){
  for(gPos subjectId=0;subjectId<i_subjectNames.size();subjectId++){
    i_allele1Ids.setValue(variantId,subjectId,allele1Id[subjectId],allele1Id.isNA(subjectId));
    i_allele2Ids.setValue(variantId,subjectId,allele2Id[subjectId],allele2Id.isNA(subjectId));
    i_isPhased.setValue(variantId,subjectId,isPhased[subjectId],isPhased.isNA(subjectId));
  }
}

void gGenotype::setSubjectGenotypes(gPos subjectId, const gArray<gShortUnsigned> & allele1Id, const gArray<gShortUnsigned> & allele2Id,const gArray<bool> & isPhased){
  for(gPos variantId=0;variantId<i_variantPositions.getSize();variantId++){
    i_allele1Ids.setValue(variantId,subjectId,allele1Id[variantId],allele1Id.isNA(variantId));
    i_allele2Ids.setValue(variantId,subjectId,allele2Id[variantId],allele2Id.isNA(variantId));
    i_isPhased.setValue(variantId,subjectId,isPhased[variantId],isPhased.isNA(variantId));
  }
}

void gGenotype::setGenotype(gPos subjectId,gPos variantId, const gArray<gShortUnsigned> & allele1Id,const gArray<gShortUnsigned> & allele2Id,const gArray<bool> & isPhased){
  i_allele1Ids.setValue(variantId,subjectId,allele1Id[0],allele1Id.isNA(0));
  i_allele2Ids.setValue(variantId,subjectId,allele2Id[0],allele2Id.isNA(0));
  i_isPhased.setValue(variantId,subjectId,isPhased[0],isPhased.isNA(0));
}

gPos gGenotype::addInterval(const gGenotypeInterval & interval){
  gPos nintv;
  vector<gGenotypeInterval>::iterator it= find(i_intervals.begin(),i_intervals.end(),interval);
  if(it==i_intervals.end()){
   nintv=i_intervals.size();
   i_intervals.push_back(interval);
  }else{
    nintv=it-i_intervals.begin();
  }
  return nintv;
}

gPos gGenotype::addSubject(const gString & subjectName){
  gPos col;
  vector<gString>::iterator it=find(i_subjectNames.begin(),i_subjectNames.end(),subjectName);
  if(it==i_subjectNames.end()){
    col=i_subjectNames.size();
    i_subjectNames.push_back(subjectName);
    if(i_variantPositions.getSize()>0){
      i_allele1Ids.setCol(col,gArray<gShortUnsigned>(0,i_variantPositions.getSize(),true));
      i_allele2Ids.setCol(col,gArray<gShortUnsigned>(0,i_variantPositions.getSize(),true));
      i_isPhased.setCol(col,gArray<gShortUnsigned>(false,i_variantPositions.getSize(),true));
    }
  }else{
    col=it-i_subjectNames.begin();
  }
  return col;
}

gPos gGenotype::addVariantPosition(gPos intervalId,gPos position, const vector<gSequence> & alleles){
  gPos row=i_variantPositions.getSize();
  if(intervalId<i_intervals.size()){
    i_variantPositions.setValue(row,position,false);
    i_variantIntervalIds.setValue(row,intervalId,false);
    i_variantAlleles.push_back(alleles);
    if(i_subjectNames.size()>0){
      i_allele1Ids.setRow(row,gArray<gShortUnsigned>(0,i_subjectNames.size(),true));
      i_allele2Ids.setRow(row,gArray<gShortUnsigned>(0,i_subjectNames.size(),true));
      i_isPhased.setRow(row,gArray<gShortUnsigned>(false,i_subjectNames.size(),true));
    }
  }else throw gException("gGenotype::addVariantPosition: invalid interval id");
  return row;
}

void gGenotype::setVariant(gPos intervalId,gPos variantId,gPos position,const vector<gSequence> & alleles, const gArray<gShortUnsigned> & allele1Id, const gArray<gShortUnsigned> & allele2Id,const gArray<bool> & isPhased){
  if(
      (intervalId<i_intervals.size()) &&
      (i_intervals[intervalId].contains(position)[0]) &&
      (variantId<i_variantPositions.getSize()) &&
      (allele1Id.getSize()==i_allele1Ids.getColsNum()) &&
      (allele2Id.getSize()==i_allele2Ids.getColsNum()) &&
      (isPhased.getSize()==i_isPhased.getColsNum())
    ){
    i_variantAlleles[variantId]=alleles;
    i_variantIntervalIds.setValue(variantId,intervalId,false);
    i_variantPositions.setValue(variantId,position,false);
    i_allele1Ids.setRow(variantId,allele1Id);
    i_allele2Ids.setRow(variantId,allele2Id);
    i_isPhased.setRow(variantId,isPhased);
  }else throw gException("gGenotype::setVariant: inconsistent data");
}

gExperimentType gGenotype::getExperimentType() const {
  return i_experimentType;
}

gSize gGenotype::getIntervalCount() const {
  return i_intervals.size();
}

gGenotypeInterval gGenotype::getInterval(gPos intervalNum) const {
  if(intervalNum<i_intervals.size()){
    return i_intervals[intervalNum];
  }else throw gException("gGenotype::getInterval: interval index out of bounds");
}

gSize gGenotype::getSubjectCount() const {
  return i_subjectNames.size();
}

gString gGenotype::getSubjectName(gPos subNum) const {
  if(subNum<i_subjectNames.size()){
    return i_subjectNames[subNum];
  }else throw gException("gGenotype::getSubjectName: subject index out of bounds");
}

gSize gGenotype::getVariantCount() const {
  return i_variantPositions.getSize();
}

gSize gGenotype::getVariantCount(gPos intervalNum) const{
  return getVariantPositions(intervalNum).getSize();
}

gPos gGenotype::getVariantIntervalNum(gPos variantNum) const {
  if(variantNum<i_variantIntervalIds.getSize()){
    return i_variantIntervalIds[variantNum];
  }else throw gException("gGenotype::getVariantIntervalNum: variant index out of bounds");
}

gArray<gPos> gGenotype::getVariantIntervals() const{
  return i_variantIntervalIds;
}

vector<gSequence> gGenotype::getVariantAlleles(gPos variantNum) const {
  if(variantNum<i_variantIntervalIds.getSize()){
    return i_variantAlleles[variantNum];
  }else throw gException("gGenotype::getVariantAlleles: variant index out of bounds");
}

vector< vector<gSequence> > gGenotype::getVariantAlleles() const{
  return i_variantAlleles;
}

gPos gGenotype::getVariantPosition(gPos variantNum) const {
  if(variantNum<i_variantIntervalIds.getSize()){
    return i_variantPositions[variantNum];
  }else throw gException("gGenotype::getVariantPosition: variant index out of bounds");
}

gMatrix<gShortUnsigned> gGenotype::getVariantGenotypes(gPos interval,bool asHaplotypes) const{
  gArray<gPos> nn=getVariantIds(interval);
  gMatrix<gShortUnsigned> ret(2*getSubjectCount(),nn.getSize(),0,true);
  if(asHaplotypes){
    for(gSize i=0;i<getSubjectCount();i++){
      gArray<gPos> nonhap=which(!i_isPhased.getCol(i)[nn]);
      gArray<gShortUnsigned> a1=i_allele1Ids.getCol(i)[nn];
      gArray<gShortUnsigned> a2=i_allele2Ids.getCol(i)[nn];
      a1.setNA(nonhap);
      a2.setNA(nonhap);
      ret.setRow(i*2,a1);
      ret.setRow(i*2+1,a2);
    }
  }else{
    for(gSize i=0;i<getSubjectCount();i++){
      ret.setRow(i*2,i_allele1Ids.getCol(i)[nn]);
      ret.setRow(i*2+1,i_allele2Ids.getCol(i)[nn]);
    }
  }
  return ret;
}

gArray<gShortUnsigned> gGenotype::getGenotype(gPos variantNum,gPos subjectNum) const{
  gArray<gShortUnsigned> gt(0,2,true);
  if( (variantNum<i_variantIntervalIds.getSize()) && (subjectNum<i_subjectNames.size()) ){
    gt.setValue(0,i_allele1Ids(variantNum,subjectNum),i_allele1Ids.isNA(variantNum,subjectNum));
    gt.setValue(1,i_allele2Ids(variantNum,subjectNum),i_allele2Ids.isNA(variantNum,subjectNum));
  }else throw gException("gGenotype::getGenotype: variant/subject index out of bounds");
  return gt;
}

gArray<bool> gGenotype::getPhaseStatus(gPos variantNum,gPos subjectNum) const{
  gArray<bool> ph(0,1,true);
  if( (variantNum<i_variantIntervalIds.getSize()) && (subjectNum<i_subjectNames.size()) ){
    ph.setValue(0,i_isPhased(variantNum,subjectNum),i_isPhased.isNA(variantNum,subjectNum));
  }else throw gException("gGenotype::getPhaseStatus: variant/subject index out of bounds");
  return ph;
}

gGenotype gGenotype::operator ()(gArray<gPos> variantIds,gArray<gPos> subjectIds) const {
  gGenotype res(subjectIds.getSize(),variantIds.getSize());
  if((subjectIds.getMax()[0]<=getSubjectCount()) && (variantIds.getMax()[0]<getVariantCount()) ){
    res.i_intervals=i_intervals; //Should exclude unused intervals (and change indexes)
    res.i_variantIntervalIds=i_variantIntervalIds[variantIds];
    res.i_variantPositions=i_variantPositions[variantIds];
    for(size_t i=0;i<variantIds.getSize();i++){
      res.i_variantAlleles[i]=i_variantAlleles[variantIds[i]];
    }
    for(size_t i=0;i<subjectIds.getSize();i++){
      res.i_subjectNames[i]=i_subjectNames[subjectIds[i]];
    }
    //there is a problem for this in geco::gArray CHECK!!!!
    //res.i_allele1Ids=gMatrix<gShortUnsigned>(variantIds.getSize(),subjectIds.getSize(),i_allele1Ids(variantIds,subjectIds));
    //res.i_allele2Ids=gMatrix<gShortUnsigned>(variantIds.getSize(),subjectIds.getSize(),i_allele2Ids(variantIds,subjectIds));
    //res.i_isPhased=gMatrix<bool>(variantIds.getSize(),subjectIds.getSize(),i_isPhased(variantIds,subjectIds));
    for(gPos v=0;v<variantIds.getSize();v++){
      for(gPos s=0;s<subjectIds.getSize();s++){
        res.i_allele1Ids.setValue(v,s,i_allele1Ids(variantIds[v],subjectIds[s]),i_allele1Ids.isNA(variantIds[v],subjectIds[s]));
        res.i_allele2Ids.setValue(v,s,i_allele2Ids(variantIds[v],subjectIds[s]),i_allele2Ids.isNA(variantIds[v],subjectIds[s]));
        res.i_isPhased.setValue(v,s,i_isPhased(variantIds[v],subjectIds[s]),i_isPhased.isNA(variantIds[v],subjectIds[s]));
      }
    }
  }else throw gException("gGenotype::(): invalid subjects and/or variants id specified");
 return res;
}

void gGenotype::addGenotypes(const gGenotype & genotype){

  //if(i_experimentType==gResequencing){
  if(i_experimentType==i_experimentType){
    if((i_variantPositions.getSize()>0)&&(i_subjectNames.size()>0)){
      gArray<gPos> imap1;
      gArray<gPos> imap2;
      vector<gGenotypeInterval> intervals=mergeGenotypeIntervals(i_intervals,genotype.i_intervals,imap1,imap2);

      vector<gString> subjects=i_subjectNames;
      gArray<gPos> smap1=getArray<gPos>(0,i_subjectNames.size()-1,1);
      gArray<gPos> smap2(0,i_subjectNames.size(),true);

      vector<gString>::iterator k;
      for(size_t i=0;i<genotype.i_subjectNames.size();i++){
	k=find(subjects.begin(),subjects.end(),genotype.i_subjectNames[i]);
	if(k!=subjects.end()){
	  smap2.setValue(k-subjects.begin(),i,false);
	}else{
	  smap2.setValue(subjects.size(),i,false);
	  smap1.setValue(subjects.size(),0,true);
	  subjects.push_back(genotype.i_subjectNames[i]);
	}
      }

      gArray<gPos> variantPositions=i_variantPositions;
      gArray<gPos> variantIntervalIds=imap1[i_variantIntervalIds];
      vector< vector<gSequence> > variantAlleles=i_variantAlleles;

      gGenotype tmpGenotype(genotype);
      gArray<gPos> vmap1=getArray<gPos>(0,i_variantPositions.getSize()-1,1);
      gArray<gPos> vmap2(0,i_variantPositions.getSize(),true);

      gArray<gPos> idxs=variantPositions.match(tmpGenotype.i_variantPositions);
      for(gPos j=0;j<idxs.getSize();j++){
	if(idxs.isNA(j)){
	  vmap1.setValue(variantPositions.getSize(),0,true);
	  vmap2.setValue(variantPositions.getSize(),j,false);
	  variantPositions.setValue(variantPositions.getSize(),tmpGenotype.i_variantPositions[j],false);
	  variantIntervalIds.setValue(variantIntervalIds.getSize(),imap2[tmpGenotype.i_variantIntervalIds[j]],false);
	  variantAlleles.push_back(tmpGenotype.i_variantAlleles[j]);

	}else{
	  if(variantAlleles[idxs[j]][0]==tmpGenotype.i_variantAlleles[j][0]){
	    vmap2.setValue(idxs[j],j,false);
	    for(size_t a=1;a<tmpGenotype.i_variantAlleles[j].size();a++){
	      vector<gSequence>::iterator np=find(variantAlleles[idxs[j]].begin()+1,variantAlleles[idxs[j]].end(),tmpGenotype.i_variantAlleles[j][a]);
	      if(np==variantAlleles[idxs[j]].end()){
		gPos newid=variantAlleles[idxs[j]].size();
		variantAlleles[idxs[j]].push_back(tmpGenotype.i_variantAlleles[j][a]);
		for(size_t s=0;s<tmpGenotype.i_subjectNames.size();s++){
		  if(tmpGenotype.i_allele1Ids(j,s)==a) tmpGenotype.i_allele1Ids.setValue(j,s,newid,tmpGenotype.i_allele1Ids.isNA(j,s));
		  if(tmpGenotype.i_allele2Ids(j,s)==a) tmpGenotype.i_allele2Ids.setValue(j,s,newid,tmpGenotype.i_allele2Ids.isNA(j,s));
		}
	      }else{
		gPos newid=np-variantAlleles[idxs[j]].begin();
		for(size_t s=0;s<tmpGenotype.i_subjectNames.size();s++){
		  if(tmpGenotype.i_allele1Ids(j,s)==a) tmpGenotype.i_allele1Ids.setValue(j,s,newid,tmpGenotype.i_allele1Ids.isNA(j,s));
		  if(tmpGenotype.i_allele2Ids(j,s)==a) tmpGenotype.i_allele2Ids.setValue(j,s,newid,tmpGenotype.i_allele2Ids.isNA(j,s));
		}
	      }
	    }
	  }else throw gException("gGenotype::addGenotype: inconsistent reference alleles");
	}
      }

      gMatrix<gShortUnsigned> a1(variantPositions.getSize(),subjects.size(),0,true);
      gMatrix<gShortUnsigned> a2(variantPositions.getSize(),subjects.size(),0,true);
      gMatrix<bool> aph(variantPositions.getSize(),subjects.size(),true,false);

      gArray<bool> mask1(false,variantPositions.getSize(),false);
      gArray<bool> mask2(false,variantPositions.getSize(),false);
      for(size_t iv=0;iv<intervals.size();iv++){
	gArray<gPos> iv1=which(imap1==iv);
	gArray<gPos> iv2=which(imap2==iv);
	gArray<gBool> civ=intervals[iv].contains(variantPositions);
	for(gPos l=0;l<iv1.getSize();l++){
	  mask1 = mask1 || ( i_intervals[iv1[l]].contains(variantPositions) && civ );
	}
	for(gPos l=0;l<iv2.getSize();l++){
	  mask2 = mask2 || ( tmpGenotype.i_intervals[iv2[l]].contains(variantPositions) && civ );
	}

      }
      gArray<gPos> gp1=which(mask1);
      gArray<gPos> gp2=which(mask2);

      for(size_t s=0;s<subjects.size();s++){
	if(!smap1.isNA(s)){
	  for(gPos l=0;l<gp1.getSize();l++){
	  a1.setValue(gp1[l],s,0,false);
	  a2.setValue(gp1[l],s,0,false);
	  aph.setValue(gp1[l],s,true,false);
	  }
	}
	if(!smap2.isNA(s)){
	  for(gPos l=0;l<gp2.getSize();l++){
	  a1.setValue(gp2[l],s,0,false);
	  a2.setValue(gp2[l],s,0,false);
	  aph.setValue(gp2[l],s,true,false);
	  }
	}
      }

      //genotypes;
      for(size_t s=0;s<subjects.size();s++){
	if((!smap1.isNA(s)) && smap2.isNA(s)){
	  for(gPos v=0;v<variantPositions.getSize();v++){
	    if(!vmap1.isNA(v)){
	      a1.setValue(v,s,i_allele1Ids(vmap1[v],smap1[s]),i_allele1Ids.isNA(vmap1[v],smap1[s]));
	      a2.setValue(v,s,i_allele2Ids(vmap1[v],smap1[s]),i_allele2Ids.isNA(vmap1[v],smap1[s]));
	      aph.setValue(v,s,i_isPhased(vmap1[v],smap1[s]),i_isPhased.isNA(vmap1[v],smap1[s]));
	    }
	  }
	}else if((!smap2.isNA(s)) && smap1.isNA(s)){
	  for(gPos v=0;v<variantPositions.getSize();v++){
	    if(!vmap2.isNA(v)){
	      a1.setValue(v,s,tmpGenotype.i_allele1Ids(vmap2[v],smap2[s]),tmpGenotype.i_allele1Ids.isNA(vmap2[v],smap2[s]));
	      a2.setValue(v,s,tmpGenotype.i_allele2Ids(vmap2[v],smap2[s]),tmpGenotype.i_allele2Ids.isNA(vmap2[v],smap2[s]));
	      aph.setValue(v,s,tmpGenotype.i_isPhased(vmap2[v],smap2[s]),tmpGenotype.i_isPhased.isNA(vmap2[v],smap2[s]));
	    }
	  }
	}else if((!smap1.isNA(s)) && (!smap2.isNA(s))){
	  for(gPos v=0;v<variantPositions.getSize();v++){
	    if((!vmap1.isNA(v)) && (vmap2.isNA(v))){
	      a1.setValue(v,s,i_allele1Ids(vmap1[v],smap1[s]),i_allele1Ids.isNA(vmap1[v],smap1[s]));
	      a2.setValue(v,s,i_allele2Ids(vmap1[v],smap1[s]),i_allele2Ids.isNA(vmap1[v],smap1[s]));
	      aph.setValue(v,s,i_isPhased(vmap1[v],smap1[s]),i_isPhased.isNA(vmap1[v],smap1[s]));
	    }else if((!vmap2.isNA(v)) && (vmap1.isNA(v))){
	      a1.setValue(v,s,tmpGenotype.i_allele1Ids(vmap2[v],smap2[s]),tmpGenotype.i_allele1Ids.isNA(vmap2[v],smap2[s]));
	      a2.setValue(v,s,tmpGenotype.i_allele2Ids(vmap2[v],smap2[s]),tmpGenotype.i_allele2Ids.isNA(vmap2[v],smap2[s]));
	      aph.setValue(v,s,tmpGenotype.i_isPhased(vmap2[v],smap2[s]),tmpGenotype.i_isPhased.isNA(vmap2[v],smap2[s]));
	    }else if((!vmap1.isNA(v)) && (!vmap2.isNA(v))){
	      if(
		  (i_isPhased(vmap1[v],smap1[s])==tmpGenotype.i_isPhased(vmap2[v],smap2[s])) &&
		  (i_isPhased.isNA(vmap1[v],smap1[s])==tmpGenotype.i_isPhased.isNA(vmap2[v],smap2[s])) && 
		  (
		  (
		    (
		      (i_allele1Ids(vmap1[v],smap1[s])==tmpGenotype.i_allele1Ids(vmap2[v],smap2[s])) &&
		      (i_allele1Ids.isNA(vmap1[v],smap1[s])==tmpGenotype.i_allele1Ids.isNA(vmap2[v],smap2[s]))
		    ) &&
		    (
		      (i_allele2Ids(vmap1[v],smap1[s])==tmpGenotype.i_allele2Ids(vmap2[v],smap2[s])) &&
		      (i_allele2Ids.isNA(vmap1[v],smap1[s])==tmpGenotype.i_allele2Ids.isNA(vmap2[v],smap2[s]))
		    )
		  ) ||
		  (
		    (
		      (i_allele1Ids(vmap1[v],smap1[s])==tmpGenotype.i_allele2Ids(vmap2[v],smap2[s])) &&
		      (i_allele1Ids.isNA(vmap1[v],smap1[s])==tmpGenotype.i_allele2Ids.isNA(vmap2[v],smap2[s]))
		    ) &&
		    (
		      (i_allele2Ids(vmap1[v],smap1[s])==tmpGenotype.i_allele1Ids(vmap2[v],smap2[s])) &&
		      (i_allele2Ids.isNA(vmap1[v],smap1[s])==tmpGenotype.i_allele1Ids.isNA(vmap2[v],smap2[s]))
		    )
		  )
		  )
	      ){
		a1.setValue(v,s,tmpGenotype.i_allele1Ids(vmap2[v],smap2[s]),tmpGenotype.i_allele1Ids.isNA(vmap2[v],smap2[s]));
		a2.setValue(v,s,tmpGenotype.i_allele2Ids(vmap2[v],smap2[s]),tmpGenotype.i_allele2Ids.isNA(vmap2[v],smap2[s]));
		aph.setValue(v,s,tmpGenotype.i_isPhased(vmap2[v],smap2[s]),tmpGenotype.i_isPhased.isNA(vmap2[v],smap2[s]));
	      }else throw gException("gGenotype::merge: inconsistent genotypes");
	    }
	  }
	}
      }

      i_intervals=intervals;
      i_variantPositions=variantPositions;
      i_variantIntervalIds=variantIntervalIds;
      i_variantAlleles=variantAlleles;
      i_subjectNames=subjects;
      i_allele1Ids=a1;
      i_allele2Ids=a2;
      i_isPhased=aph;
    }else{
      if(i_subjectNames.size()>0){
	if(genotype.getVariantCount()>0){
	  //*this=merge(*this,genotype);
	  *this=merge(genotype,*this);
	}else{
	  for(size_t i=0;i<genotype.i_subjectNames.size();i++){
	    this->addSubject(genotype.i_subjectNames[i]);
	  }
	}
      }else{
	*this=genotype;
      }
    }
    sortVariants();
  }else{

  }
  if(i_experimentType==gGenotyping){
    vector<gGenotypeInterval>::iterator i,j;
    gArray<gPos> valid(1,i_intervals.size(),false);
    vector<gGenotypeInterval> newIntervals;
    for(gPos i=0;i<i_intervals.size();i++){
      if(valid[i]){
	newIntervals.push_back(i_intervals[i]);
	gArray<gPos> tr=which(i_variantIntervalIds==i);
	i_variantIntervalIds.setValues(tr,newIntervals.size()-1);
	valid.setValue(i,0,false);
	gString contig1=i_intervals[i];
	for(gPos j=0;j<i_intervals.size();j++){
	  if(valid[j]){
	    gString contig2=i_intervals[j];
	    if(contig2==contig1){
	      valid.setValue(j,0,false);
	      gArray<gPos> pstart,pend;
	      pstart.setValue(0,i_intervals[i].getStart(),!i_intervals[i].validStart());
	      pstart.setValue(1,i_intervals[j].getStart(),!i_intervals[j].validStart());
	      pend.setValue(0,i_intervals[i].getEnd(),!i_intervals[i].validEnd());
	      pend.setValue(1,i_intervals[j].getEnd(),!i_intervals[j].validEnd());
	      gArray<gPos> minStart=pstart.getMin();
	      gArray<gPos> maxEnd=pend.getMax();
	      newIntervals[newIntervals.size()-1]=gGenotypeInterval(contig1,minStart[0],maxEnd[0],!(minStart.isNA(0)),!(maxEnd.isNA(0)));
	      gArray<gPos> tr=which(i_variantIntervalIds==j);
	      i_variantIntervalIds.setValues(tr,newIntervals.size()-1);
	    }
	  }
	}
      }
    }
    i_intervals=newIntervals;
  }
}

void gGenotype::sortVariants(){
  if(i_variantPositions.getSize()>0){
    gArray<gPos> opos=i_variantPositions.sort(true);
    if(!i_isForward){
      i_variantPositions.revert();
      opos.revert();
    }
    gArray<gPos> tmpIntervalIds=i_variantIntervalIds;
    i_variantIntervalIds=tmpIntervalIds[opos];
    vector< vector<gSequence> > tmpAlleles=i_variantAlleles;
    gMatrix<gShortUnsigned> a1=i_allele1Ids;
    gMatrix<gShortUnsigned> a2=i_allele2Ids;
    gMatrix<bool> ph=i_isPhased;
    for(gPos i=0;i<opos.getSize();i++){
      i_variantAlleles[i]=tmpAlleles[opos[i]];
      i_allele1Ids.setRow(i,a1.getRow(opos[i]));
      i_allele2Ids.setRow(i,a2.getRow(opos[i]));
      i_isPhased.setRow(i,ph.getRow(opos[i]));
    }
  }
}


void gGenotype::sortSubjects(const vector<gString> & names){
  if(names.size()==i_subjectNames.size()){
    gArray<gSize> cols(0,i_subjectNames.size(),true);
    for(gSize i=0;i<i_subjectNames.size();i++){
      vector<gString>::const_iterator it=find(names.begin(),names.end(),i_subjectNames[i]);
      if(it!=names.end()){
	cols.setValue(it-names.begin(),i,false);
      }else throw gException("gGenotype::sortSubjects: invalid sort name");
    }
    i_allele1Ids=i_allele1Ids.getColumns(cols);
    i_allele2Ids=i_allele2Ids.getColumns(cols);   
    i_isPhased=i_isPhased.getColumns(cols);
    i_subjectNames=names;
  }else throw gException("gGenotype::sortSubjects: invalid sort names length");
}


gArray<gPos> gGenotype::getVariantIds(gPos intervalNum) const{
  return which(i_variantIntervalIds==intervalNum);
}

gArray<gPos> gGenotype::getVariantPositions(gPos intervalNum) const{
  return i_variantPositions[getVariantIds(intervalNum)];
}

void gGenotype::removeVariants(const gArray<bool> & removeFlags,bool splitIntervals){
  gArray<bool> keepFlags=(!removeFlags);
  gArray<gPos> toKeep=which(!removeFlags);
  if(toKeep.getSize()!=i_variantPositions.getSize()){
    if(toKeep.getSize()>0){
      gArray<gPos>                  keepIntervals=i_variantIntervalIds[toKeep];
      gArray<gPos>                  keepPositions=i_variantPositions[toKeep];
      vector< vector<gSequence> >   keepAlleles(toKeep.getSize(),vector<gSequence>());
      gMatrix<gShortUnsigned>       keepAllele1Ids(toKeep.getSize(),i_subjectNames.size(),0,true);
      gMatrix<gShortUnsigned>       keepAllele2Ids(toKeep.getSize(),i_subjectNames.size(),0,true);
      gMatrix<bool>                 keepIsPhased(toKeep.getSize(),i_subjectNames.size(),0,true);

      if(splitIntervals){
        gPos currInterval=keepIntervals[0];
        for(gPos i=0;i<toKeep.getSize()-1;i++){

          if(keepIntervals[i]==keepIntervals[i+1]){
            if(toKeep[i+1]-toKeep[i]>1){
              cout << "intv1:\t" << i_intervals[currInterval] << ":" << i_intervals[currInterval].getStart() << "-" << i_variantPositions[toKeep[i]+1] << endl;
              i_intervals.push_back(gGenotypeInterval(i_intervals[currInterval],i_variantPositions[toKeep[i+1]-1]+1,i_intervals[currInterval].getEnd(),true,true));
              gPos start=i_intervals[currInterval].getStart();
              gPos end=i_variantPositions[toKeep[i]+1];
              i_intervals[currInterval]=gGenotypeInterval(i_intervals[currInterval],start,end,true,true);
              for(gPos p=i+1;p<toKeep.getSize();p++){
                if(keepIntervals[p]==currInterval){
                  keepIntervals.setValue(p,i_intervals.size()-1,false);
                }else{
                  break;
                }
              }
              currInterval=i_intervals.size()-1;
              cout << "intv2:\t" << i_intervals[currInterval] << ":" << i_variantPositions[toKeep[i+1]-1]+1 << "-" << i_intervals[currInterval].getEnd() << endl;
            }
          }else currInterval=keepIntervals[i+1];
        }
      }

      for(gPos i=0;i<toKeep.getSize();i++){
        keepAlleles[i]=i_variantAlleles[toKeep[i]];
        keepAllele1Ids.setRow(i,i_allele1Ids.getRow(toKeep[i]));
        keepAllele2Ids.setRow(i,i_allele2Ids.getRow(toKeep[i]));
        keepIsPhased.setRow(i,i_isPhased.getRow(toKeep[i]));
      }

      i_variantIntervalIds=keepIntervals;
      i_variantPositions=keepPositions;
      i_variantAlleles=keepAlleles;
      i_allele1Ids=keepAllele1Ids;
      i_allele2Ids=keepAllele2Ids;
      i_isPhased=keepIsPhased;
    }else{
      i_intervals.clear();
      i_variantPositions=gArray<gPos>();
      i_variantIntervalIds=gArray<gPos>();
      i_variantAlleles.clear();
      i_allele1Ids=gMatrix<gShortUnsigned>();
      i_allele2Ids=gMatrix<gShortUnsigned>();
      i_isPhased=gMatrix<bool>();
    }
  }
}

void gGenotype::removeVariantsByNumbers(const gArray<gPos> & variantNums,bool splitIntervals){
  gArray<bool> removeFlags(false,i_variantPositions.getSize(),false);
  removeFlags.setValues(variantNums,true);
  removeVariants(removeFlags,splitIntervals);
}

void gGenotype::removeVariantsByPositions(const gArray<gPos> & variantPositions,bool splitIntervals){
  //removeVariantsByNumbers(which(!i_variantPositions.match(variantPositions).isNA()),splitIntervals);
  gArray<gPos> a=i_variantPositions.match(variantPositions);
  gArray<gPos> b=which(!a.isNA());
  gArray<gPos> c=a[b];
  removeVariantsByNumbers(c,splitIntervals);
}

void gGenotype::removeVariantsByIntervals(const vector<gGenotypeInterval> & intervals,bool splitIntervals){
  gArray<bool> removeFlags(false,i_variantPositions.getSize(),false);
  for(size_t i=0;i<intervals.size();i++){
    for(size_t j=0;j<i_intervals.size();j++){
      if( ((gString) intervals[i])==((gString) i_intervals[j]) ){
        gArray<gPos> candnum=which(i_variantIntervalIds==j);
        gArray<gPos> candpos=i_variantPositions[candnum];
        gReferenceInterval intersect=intervals[i].getIntersection(i_intervals[j]);
        gArray<gPos> ptrm=which(intersect.contains(candpos));
        if(ptrm.getSize()>0){
          removeFlags.setValues(candnum[ptrm],true);
        }
      }
    }
  }
  removeVariants(removeFlags,splitIntervals);
}

void gGenotype::removeInvariantPositions(){
  gArray<bool> removeFlags(false,i_variantPositions.getSize(),false);
  for(gSize i=0;i<i_variantPositions.getSize();i++){
    gArray<bool> cond=(   i_allele1Ids.getRow(i)!=i_allele2Ids.getRow(i)).getSum() || i_allele1Ids.getRow(i).getUnique().getSize()>1;
    if(cond[0]==0){
     removeFlags.setValue(i,true,false);
    }
  }
  removeVariants(removeFlags,false);
}

vector<gSequence> gGenotype::getHaplotypes(gPos subjectNum,gPos intervalNum) const{
  gArray<gPos> varIds=getVariantIds(intervalNum);
  if(varIds.getSize()>0){
    vector<gSequence> ret(2,gSequence(varIds.getSize()));
    gArray<gShortUnsigned>  a1=i_allele1Ids.getCol(subjectNum)[varIds];
    gArray<gShortUnsigned>  a2=i_allele2Ids.getCol(subjectNum)[varIds];
    for(gPos i=0;i<varIds.getSize();i++){
      if(!a1.isNA(i)){
	ret[0].setValue(i,i_variantAlleles[varIds[i]][a1[i]],false);
      }else{
	ret[0].setValue(i,'N',true);
      }
      if(!a2.isNA(i)){
	ret[1].setValue(i,i_variantAlleles[varIds[i]][a2[i]],false);
      }else{
	ret[1].setValue(i,'N',true);
      }
    }
    return ret;
  }else{
    return vector<gSequence>(2);
  }
}

gArray<gDoubleScore> gGenotype::getMaf(gPos interval){
  gArray<gPos> selvar=getVariantIds(interval);
  gArray<gDoubleScore> ret(0,selvar.getSize(),true);
  for(gSize i=0;i<selvar.getSize();i++){
    if(i_variantAlleles[selvar[i]].size()==2){
      gArray<gDoubleScore> a1=i_allele1Ids.getRow(selvar[i])==0;
      gArray<gDoubleScore> a2=i_allele2Ids.getRow(selvar[i])==0;
      gArray<gDoubleScore> s1=a1.getSum();
      gArray<gDoubleScore> s2=a2.getSum();
      
//      gArray<gDoubleScore> s1=a1.getSum(0,0,true);
//      gArray<gDoubleScore> s2=a2.getSum(0,0,true);
      //gArray<gDoubleScore> n=gArray<gDoubleScore>(a1.getSize()-a1.NACount()+a2.getSize()-a2.NACount(),1,false);
      gArray<gDoubleScore> n=gArray<gDoubleScore>(2*a1.getSize(),1,false);
      gArray<gDoubleScore> ss=(s1+s2)/n;
      if(ss[0]<0.5){
	ret.setValue(i,ss[0],ss.isNA(0));
      }else{
	ret.setValue(i,(gDoubleScore)1-ss[0],ss.isNA(0));
      }
    }
  }
  return ret;
}

void gGenotype::reverse(){
  i_isForward=!i_isForward;
  for(size_t i=0;i<i_variantAlleles.size();i++){
    for (size_t a=0;a<i_variantAlleles[i].size();a++){
      i_variantAlleles[i][a].reverseComplement();
    }
  }
  sortVariants();
}

bool gGenotype::isForward() const{
  return i_isForward;
}

gArray<gPos> gGenotype::checkDiscordancies(bool markNA){
  gArray<gPos> errPos;
  for(size_t i=0;i<i_intervals.size();i++){
    gArray<gPos> variantIds=getVariantIds(i);
    gArray<gPos> pos=i_variantPositions[variantIds];
    gArray<gPos> ids=variantIds[pos.sort(true)];
    gArray<gPos> a(pos,0,pos.getSize()-1);
    gArray<gPos> b(pos,1,pos.getSize());
    gArray<gPos> err=which(a==b);
    if(err.getSize()>0){
      errPos=b[err];
      //genotype columns at ids[err] should be marked as NA
      for(gSize j=0;j<err.getSize();j++){
	for(gSize s=0;s<getSubjectCount();s++){
          i_allele1Ids.setValue(ids[err[j]],s,0,true);
          i_allele2Ids.setValue(ids[err[j]],s,0,true);
          i_isPhased.setValue(ids[err[j]],s,0,true);
	}
      }
      removeVariantsByNumbers(ids[err+1],false);
    }
  }
  return errPos;
}

void gGenotype::setExperimentType(gExperimentType type){
  if(type!=i_experimentType){
    i_experimentType=type;
    switch(type){
      case gResequencing:

	break;
      case gGenotyping:
	vector<gGenotypeInterval>::iterator i,j;
	gArray<gPos> valid(1,i_intervals.size(),false);
	vector<gGenotypeInterval> newIntervals;
	for(gPos i=0;i<i_intervals.size();i++){
	  if(valid[i]){
	    newIntervals.push_back(i_intervals[i]);
	    gArray<gPos> tr=which(i_variantIntervalIds==i);
	    i_variantIntervalIds.setValues(tr,newIntervals.size()-1);
	    valid.setValue(i,0,false);
	    gString contig1=i_intervals[i];
	    for(gPos j=0;j<i_intervals.size();j++){
	      if(valid[j]){
		gString contig2=i_intervals[j];
		if(contig2==contig1){
		  valid.setValue(j,0,false);
		  gArray<gPos> pstart,pend;
		  pstart.setValue(0,i_intervals[i].getStart(),!i_intervals[i].validStart());
		  pstart.setValue(1,i_intervals[j].getStart(),!i_intervals[j].validStart());
		  pend.setValue(0,i_intervals[i].getEnd(),!i_intervals[i].validEnd());
		  pend.setValue(1,i_intervals[j].getEnd(),!i_intervals[j].validEnd());
		  gArray<gPos> minStart=pstart.getMin();
		  gArray<gPos> maxEnd=pend.getMax();
		  newIntervals[newIntervals.size()-1]=gGenotypeInterval(contig1,minStart[0],maxEnd[0],!(minStart.isNA(0)),!(maxEnd.isNA(0)));
		  gArray<gPos> tr=which(i_variantIntervalIds==j);
		  i_variantIntervalIds.setValues(tr,newIntervals.size()-1);
		}
	      }
	    }
	  }
	}
	i_intervals=newIntervals;
	break;
    }
  }
}

ostream & gGenotype::print(ostream & out,bool useAlleleSeqs) const{
  for(size_t i=0;i < i_intervals.size();i++){
    out << "intervals:\t" << i_intervals[i] << ":" << i_intervals[i].getStart() << "-" << i_intervals[i].getEnd() << endl;
    out << "________________________________________________________________________" << endl;
    out << "\tREF\tALT\t";
    for(size_t j=0;j<i_subjectNames.size();j++){
      out << i_subjectNames[j] << "\t";
    }
    out << endl;
    gArray<gPos> selvariants=which(i_variantIntervalIds==i);
    for(gPos j=0;j<selvariants.getSize();j++){
      out << i_intervals[i_variantIntervalIds[selvariants[j]]] << ":" << i_variantPositions[selvariants[j]] << "\t";
      out << i_variantAlleles[selvariants[j]][0] << "\t";
      for (size_t k=1;k<i_variantAlleles[selvariants[j]].size();k++){
        out << i_variantAlleles[selvariants[j]][k] << ", ";
      }
      out << "\t";
      for(size_t k=0;k<i_subjectNames.size();k++){

        if(i_allele1Ids.isNA(selvariants[j],k)){
          out << "*";
        }else{
          if(useAlleleSeqs){
            out << i_variantAlleles[selvariants[j]][i_allele1Ids(selvariants[j],k)];
          }else{
            out << i_allele1Ids(selvariants[j],k);
          }
        }

        if(i_isPhased.isNA(selvariants[j],k)){
          out << "?";
        }else{
          if(i_isPhased(selvariants[j],k)){
            out << "|";
          }else{
            out << "/";
          }
        }

        if(i_allele2Ids.isNA(selvariants[j],k)){
          out << "*";
        }else{
          if(useAlleleSeqs){
            out << i_variantAlleles[selvariants[j]][i_allele2Ids(selvariants[j],k)];
          }else{
            out << i_allele2Ids(selvariants[j],k);
          }
        }
        out << "\t";
      }
      out << endl;
    }
    out << "________________________________________________________________________" << endl << endl;
  }
  return out;
}
//------------------------------------------------------------------------------


//--------------------------------------------------------------------------
/// class gGenotypeGroups definition
//--------------------------------------------------------------------------
gGenotypeGroups::gGenotypeGroups():gGenotype(){

}

gGenotypeGroups::gGenotypeGroups(const gGenotype & genotypes,const gString groupName):gGenotype(){
  addGroup(genotypes,groupName);
}

gGenotypeGroups::gGenotypeGroups(const vector<gGenotype> & genotypeGroups,const vector<gString> groupNames):gGenotype(){
  if(genotypeGroups.size()==groupNames.size()){
    for(size_t i=0;i<genotypeGroups.size();i++){
      addGroup(genotypeGroups[i],groupNames[i]);
    }
  }else throw(gException("gGenotypeGroups: size of gentypes does'nt match groupNames"));
}

gGenotypeGroups::gGenotypeGroups(const gGenotypeGroups & groups):gGenotype(groups){
  i_groupNames=groups.i_groupNames;
  i_groupStart=groups.i_groupStart;
  i_groupEnd=groups.i_groupEnd;
}

void gGenotypeGroups::addGroup(const gGenotype & groupGenotype,const gString groupName){
  i_groupStart.setValue(i_groupNames.size(),getSubjectCount(),false);
  addGenotypes(groupGenotype);
  i_groupEnd.setValue(i_groupNames.size(),getSubjectCount(),false);
  i_groupNames.push_back(groupName);
}

gGenotypeGroups gGenotypeGroups::getGroups(const gArray<gPos> & groupsIndex,bool removeInvariant) const{
 gGenotypeGroups ret;
 for(gPos i=0;i<groupsIndex.getSize();i++){
   ret.addGroup(gGenotype(*this,i_groupStart[groupsIndex[i]],i_groupEnd[groupsIndex[i]],false),i_groupNames[groupsIndex[i]]);
 }
 if(removeInvariant){
   ret.removeInvariantPositions();
 }
 return ret;
}

gArray<gSize> gGenotypeGroups::getGroupsSize(const gArray<gSize> & groupIndex) const{
  if(groupIndex.getSize()==0){
    return i_groupEnd-i_groupStart;
  }else{
    return i_groupEnd[groupIndex]-i_groupStart[groupIndex];
  }
}

gString gGenotypeGroups::getGroupName(gSize groupIndex) const{
 return i_groupNames[groupIndex];
}

gSize gGenotypeGroups::getGroupCount() const{
 return i_groupNames.size();
}
//--------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
gGenotype geco::merge(const gGenotype & genotype1,const gGenotype & genotype2){
  gGenotype ret=genotype1;
  ret.addGenotypes(genotype2);
  return ret;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
vector<gVariations> geco::getVariations(const gGenotype & genotypes,const vector<gString> & subjects,const gArray<gPos> & variantIds){
  vector<gVariations> ret;
  if(subjects.size()==0){
    for(size_t intv=0;intv<genotypes.getIntervalCount();intv++){
      gGenotypeInterval interval=genotypes.getInterval(intv);
      for(size_t sub=0;sub<genotypes.getSubjectCount();sub++){
        gVariations h1(interval),h2(interval);
	if(variantIds.getSize()==0){
	  for(gSize var=0;var<genotypes.getVariantCount();var++){
	    gPos pos=genotypes.getVariantPosition(var);
	    vector<gSequence> alleles=genotypes.getVariantAlleles(var);
	    gArray<gShortUnsigned> gt=genotypes.getGenotype(var,sub);
	    if(gt[0]!=0){
	      if(alleles[0].getLength()>0){
		if(alleles[gt[0]].getLength()==alleles[0].getLength()){
		  h1.addSubstitution(interval,gReferenceInterval(pos,pos+alleles[gt[0]].getLength()),alleles[gt[0]]);
		}
	      }
	    }
	    if(gt[1]!=0){
	      if(alleles[0].getLength()>0){
		if(alleles[gt[1]].getLength()==alleles[0].getLength()){
		  h2.addSubstitution(interval,gReferenceInterval(pos,pos+alleles[gt[1]].getLength()),alleles[gt[1]]);
		}
	      }
	    }
	  }
	}else{
	  for(gPos i=0;i<variantIds.getSize();i++){
	    gPos pos=genotypes.getVariantPosition(variantIds[i]);
	    vector<gSequence> alleles=genotypes.getVariantAlleles(variantIds[i]);
	    gArray<gShortUnsigned> gt=genotypes.getGenotype(variantIds[i],sub);
	    if(gt[0]!=0){
	      if(alleles[0].getLength()>0){
		if(alleles[gt[0]].getLength()==alleles[0].getLength()){
		  h1.addSubstitution(interval,gReferenceInterval(pos,pos+alleles[gt[0]].getLength()),alleles[gt[0]]);
		}
	      }
	    }
	    if(gt[1]!=0){
	      if(alleles[0].getLength()>0){
		if(alleles[gt[1]].getLength()==alleles[0].getLength()){
		  h2.addSubstitution(interval,gReferenceInterval(pos,pos+alleles[gt[1]].getLength()),alleles[gt[1]]);
		}
	      }
	    }
	  }
	}
        ret.push_back(h1);
        ret.push_back(h2);
      }
    }
  }else{


  }
  return ret;
}
//------------------------------------------------------------------------------
/*
void outputPedMapfiles(const gGenodata & object, gString& pedmapFile){
  ofstream ped(((string)(pedmapFile+".ped")).c_str());
  gReferenceInterval region=object.getCoveredRegion();
  ped << "# ped file generated from gGenodata object related to region " << object.getChrom() << ":" << region.getStart() <<"-"<< region.getEnd() << endl;
  ped << "# familyId, dadId, mumId, sex and phenotype are missing: use options --no-fid --no-parents --no-sex --no-pheno to read this file properly! " << endl;

  for (gSize i=0;i<object.getSubjectsCount();i++){
    ped << object.getSubjectName(i) << " ";
    for (gSize j=0;j<object.getVariantsCount();j++){
      gString h1=object.getHaplotype1Sequence(i,j);
      gString h2=object.getHaplotype1Sequence(i,j);
      if(h1.getLength()==1) ped << h1 << " ";
      else ped << "0 ";
      if(h2.getLength()==1) ped << h2 << " ";
      else ped << "0 ";
    }
    ped << endl;
  }
  ped.close();

  ofstream map(((string)(pedmapFile+".map")).c_str());
  gString chromId(object.getChrom(),3,object.getChrom().getLength());
  for (gSize i=0;i<object.getVariantsCount();i++) map << chromId << " snp" << i+1 << object.getVariantPosition(i) << " 0" << endl;
  map.close();
}
*/

