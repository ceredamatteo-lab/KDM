#pragma once
/**
 * @file geco_genotype.h
 * @author  Uberto Pozzoli  <uberto.pozzoli@bp.lnf.it>
 * @author  Federico Marini <federico.marini@bp.lnf.it>
 * @version 0.1
 *
 * @section LICENSE
 * Copyright (C) 2011 by Uberto Pozzoli and Federico Marini
 * uberto.pozzoli@bp.lnf.it
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * @section DESCRIPTION
 * ...
 */
#include "element.h"

namespace geco {

  typedef enum {gResequencing,gGenotyping} gExperimentType;

  class gGenotypeInterval:public gReferenceInterval,public gString{
  public:
    gGenotypeInterval();
    gGenotypeInterval(const gString & name,gPos start,gPos end,bool validstart,bool validend);
    gGenotypeInterval(const gString & name,const gReferenceInterval & interval);
    gGenotypeInterval(const gString & name,const gArray<gPos> & positions);
    gGenotypeInterval(const gGenotypeInterval & interval);

    bool isWholeReference() const;

    gArray<bool> contains(const gArray<gPos> & positions) const;
    //gArray<bool> contains(const gArray<gPos> & positions);
    bool operator == (const gGenotypeInterval & interval) const;
    bool operator < (const gGenotypeInterval & interval) const;
    bool operator > (const gGenotypeInterval & interval) const;
  };

  class gGenotype{
  private:
    gString                       i_name;
    gExperimentType               i_experimentType;
    bool                          i_isForward;
    bool                          i_hasReferenceAlleles;

    std::vector<gGenotypeInterval>     i_intervals;
    std::vector<gString>               i_subjectNames;

    gArray<gPos>                  i_variantPositions;
    gArray<gPos>                  i_variantIntervalIds;
    std::vector< std::vector<gSequence> >   i_variantAlleles;
    gMatrix<gShortUnsigned>       i_allele1Ids;
    gMatrix<gShortUnsigned>       i_allele2Ids;
    gMatrix<bool>                 i_isPhased;

  protected:
    void removeVariants(const gArray<bool> & keepFlags,bool splitIntervals);

  public:
    gGenotype(gExperimentType type=gResequencing);
    gGenotype(gSize subjectNum,gSize variantNum,gExperimentType type=gResequencing);

    gGenotype(const std::vector<gString> & subjectsNames, gSize variantNum,gExperimentType type=gResequencing);

    gGenotype(const gGenotype & genotype);
    gGenotype(const gGenotype & genotype,gPos startSubject,gPos endSubject,bool removeInvariant=false);
    gGenotype(const gGenotype & genotype,gPos interval, gPos startPosition,gPos endPosition);

    gGenotype & operator =(const gGenotype & genotype);
    void clear();
    void setName(const gString & name);
    void setIntervals(const std::vector<gGenotypeInterval> intervals);
    void setSubjects(const std::vector<gString> & subjectsNames);
    void setVariantPositions(const gArray<gPos> & intervalIds,const gArray<gPos> & positions,const std::vector< std::vector<gSequence> > & alleles);
    void setGenotypes(const gMatrix<gShortUnsigned> & allele1Ids,const gMatrix<gShortUnsigned> & allele2Ids,const gMatrix<bool> & isPhased);
    void setVariantGenotypes(gPos variantId, const gArray<gShortUnsigned> & allele1Id, const gArray<gShortUnsigned> & allele2Id,const gArray<bool> & isPhased);
    void setSubjectGenotypes(gPos subjectId, const gArray<gShortUnsigned> & allele1Id, const gArray<gShortUnsigned> & allele2Id,const gArray<bool> & isPhased);
    void setGenotype(gPos subjectId,gPos variantId, const gArray<gShortUnsigned> & allele1Id, const gArray<gShortUnsigned> & allele2Id,const gArray<bool> & isPhased);

    gPos addInterval(const gGenotypeInterval & interval);
    gPos addSubject(const gString & subjectName);
    gPos addVariantPosition(gPos intervalId,gPos position, const std::vector<gSequence> & alleles);
    void setVariant(gPos intervalId,gPos variantId,gPos position, const std::vector<gSequence> & alleles, const gArray<gShortUnsigned> & allele1Id, const gArray<gShortUnsigned> & allele2Id,const gArray<bool> & isPhased);

    gExperimentType getExperimentType() const;
    void setExperimentType(gExperimentType type);

    gSize getIntervalCount() const;
    gGenotypeInterval getInterval(gPos intervalNum) const;

    gSize getSubjectCount() const;
    gString getSubjectName(gPos subNum) const;

    gSize getVariantCount() const;
    gSize getVariantCount(gPos intervalNum) const;
    gPos getVariantIntervalNum(gPos variantNum) const;
    gArray<gPos> getVariantIntervals() const;
    gArray<gPos> getVariantIds(gPos intervalNum) const;

    std::vector<gSequence> getVariantAlleles(gPos variantNum) const;
    std::vector< std::vector<gSequence> > getVariantAlleles() const;
    gPos getVariantPosition(gPos variantNum) const;
    gArray<gPos> getVariantPositions(gPos intervalNum) const;

    gArray<gShortUnsigned> getGenotype(gPos variantNum,gPos subjectNum) const;
    gArray<bool> getPhaseStatus(gPos variantNum,gPos subjectNum) const;

    gGenotype operator () (gArray<gPos> variantIds,gArray<gPos> subjectIds) const;

    void addGenotypes(const gGenotype & genotype);

    void sortVariants();
    void sortSubjects(const std::vector<gString> & names);

    gArray<gPos> getInvariantPositions();
    gArray<gPos> getInvariantIndices();

    void removeVariantsByIntervals(const std::vector<gGenotypeInterval> & intervals,bool splitIntervals=true);
    void removeVariantsByNumbers(const gArray<gPos> & variantNums,bool splitIntervals=true);
    void removeVariantsByPositions(const gArray<gPos> & variantPositions,bool splitIntervals=true);
    void removeInvariantPositions();

    std::vector<gSequence> getHaplotypes(gPos subjectNum,gPos intervalNum=0) const;
    gMatrix<gShortUnsigned> getVariantGenotypes(gPos interval,bool asHaplotypes=false) const;
    gArray<gDoubleScore> getMaf(gPos interval);
    void reverse();
    bool isForward() const;
    gArray<gPos> checkDiscordancies(bool markNA);

    std::ostream & print(std::ostream & out,bool useAlleleSeq=false) const;
  };

  class gGenotypeGroups:public gGenotype{
  private:
    gArray<gPos> i_groupStart;
    gArray<gPos> i_groupEnd;
    std::vector<gString> i_groupNames;

  protected:
  public:
    gGenotypeGroups();
    gGenotypeGroups(const gGenotype & groupGenotype,const gString groupName);
    gGenotypeGroups(const std::vector<gGenotype> & groupsGenotype,const std::vector<gString> groupNames);
    gGenotypeGroups(const gGenotypeGroups & groups);
    void addGroup(const gGenotype & groupGenotype,const gString groupName);
    gGenotypeGroups getGroups(const gArray<gPos> & groupsIndex,bool removeInvariant=true) const;
    gArray<gSize> getGroupsSize(const gArray<gSize> & groupIndex=gArray<gSize>()) const;
    gString getGroupName(gSize groupIndex) const;
    gSize getGroupCount() const;

  };

  gGenotype merge(const gGenotype & genotype1,const gGenotype & genotype2);

  std::vector<gVariations> getVariations(const gGenotype & genotypes,const std::vector<gString> & subjects=std::vector<gString>(),const gArray<gPos> & variantIds=gArray<gPos>());

}
