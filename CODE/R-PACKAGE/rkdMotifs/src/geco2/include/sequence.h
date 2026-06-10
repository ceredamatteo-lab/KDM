#pragma once
/**
 * @file sequence.h
 * @author  Uberto Pozzoli  <uberto.pozzoli@bp.lnf.it>
  * @version 0.1
 *
 * @section LICENSE
 * Copyright (C) 2010 by Uberto Pozzoli and Matteo Cereda
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
 * This file contains the declarations of teh gElement class as well as
 * those of other classes used with gElement objects
 */
#include "array.h"

namespace geco {

  /** @brief Sequence Class
  *
  * String Class specialization to manage genomic Sequences.
  */
  class gSequence:public gString {
  public:
      /** @name Constructors */
      //@{
      gSequence();
      gSequence(const gArray<gChar> & string,gPos start=0,gPos end=0);
      gSequence(const gChar * str,gPos start=0, gPos end=0);
      gSequence(const std::string & str,gPos start=0, gPos end=0);
      gSequence(gSize length,gChar fillvalue='N');
      //@}

      /** @name Sequence manipulation function members (self)*/
      //@{
      gSequence & reverseComplement();
      //@}

      /** @name Sequence manipulation function members (instatiate new objects)*/
      //@{
      gSequence getReverseComplement() const;
      gString getTranslated(gShortUnsigned frame);
	  gString getMtTranslated(gShortUnsigned frame);
      //@}

  };


  /** @brief Sequence Group Class
  *
  * Represent a group of sequences possibly aligned and with quality scores
  */
  class gSequenceGroup{
  private:
    std::vector<gSequence> i_sequences;
    std::vector<gString> i_sequenceNames;
    std::vector< std::vector< gMatrix<gDoubleScore> > > i_features;
    std::vector< std::vector<gString> > i_featureNames;

    gSize i_referenceId;
    bool i_aligned;

  public:
    /** @name Constructors */
    //@{
    gSequenceGroup(bool aligned=false);
    gSequenceGroup(const gSequenceGroup & group);
    //@}

    /** @name Sequence add functions*/
    void addSequence(const gString & name,const gSequence & sequence,bool revert=false);
    void setSequence(gSize seqNum, const gString & name,const gSequence & sequence,bool revert=false);
    void addRefSequence(const gString & name,const gSequence & seq,bool revert=false);

    void addFeature(gSize seqNum, const gString & name, const gMatrix<gDoubleScore> & feature);
    //@}

    /** @name Sequence info functions*/
    gSize getSequenceCount() const;
    gString getSequenceName(gPos seqNum) const;
    gSequence getSequence(gPos seqNum) const;

    gSize getFeaturesCount(gPos seqNum) const;
    gMatrix<gDoubleScore> getFeature(gPos seqNum,gPos featNum) const;
    gString getFeatureName(gPos seqNum,gPos featNum) const;

    gSize getReferenceId() const;
    bool isAligned() const;
    //@}
  };


  /** @brief Sequence Retriever Implementation class
  *
  * Abstract base class from which users can derive specific objects to fetch sequence data from any source
  */
  class gSequenceRetrieverImplementation:public gArrayRetrieverImplementation<gChar>{
    friend class gSequenceRetriever;
  private:
    /** @brief Sequence retriever pure virtual function method
    *
    * You must implement this function in derive objects
    * @param reference a reference to the sequence you must retrieve
    * @param start start position
    * @param end end position
    * @return a gSequence object containing the required sequence;
    */
    virtual gSequence getSequence_Internal(const gString & reference,gPos start,gPos end) const = 0;

    /** @brief Object cloning pure virtual function method
    *
    * You must implement this function in derive objectsn to return a copy of the cyrrent object
    * @return a pointer to a copy of the current object.
    */
    virtual gRetrieverImplementation * clone() const = 0;
  public:
      gSequenceRetrieverImplementation();
      virtual ~gSequenceRetrieverImplementation();

  };

  /** @brief Sequence Retriever Class
  *
  * Class used to retrieve sequence information using a provided implementation
  * it provides gElement with sequence retieving capapbilities
  */
  class gSequenceRetriever:public gArrayRetriever<gChar>{
  private:
  public:
      gSequenceRetriever();
      gSequenceRetriever(const gSequenceRetrieverImplementation & implementation);
      virtual ~gSequenceRetriever();
      gSequence getSequence(const gString & reference,gPos start,gPos end) const;
  };

}

