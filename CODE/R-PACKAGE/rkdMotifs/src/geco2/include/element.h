#pragma once
/**
 * @file geco_element.h
 * @author  Uberto Pozzoli  <uberto.pozzoli@bp.lnf.it>
 * @author  Matteo Cereda   <matteo.cereda@bp.lnf.it>
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
#include "sequence.h"

#include <vector>
#include <map>


namespace geco {

  class gElement;

  /** @brief Genomic Variations
  *
  * A class to easily manage groups of genomic variations
  */
  class gVariations {
      friend class gRegion;
  protected:
      gString i_referenceName;
      gArray<gPos> i_pos;
      gArray<gPos> i_varStart;
      gArray<gPos> i_varEnd;
      gArray<gVariationType> i_varType;
      std::vector<gSequence> i_varSeq;

      gArrayIndex addVariation(const gString & referenceName,const gReferenceInterval & referenceInterval,gVariationType type,const gSequence & seq);
      gArrayIndex addVariation(const gElement & element,const gElementInterval & elementInterval,gVariationType type,const gSequence &seq);
  public:

      /** @name Constructors */
      //@{
      gVariations();
      gVariations(const gString & referenceName);
      gVariations(const gElement & element);
      gVariations(const gVariations & variations);
      //@}

      /** @name Copy operators */
      //@{
      gVariations & operator = (const gVariations & variations);
      //@}

      /** @name Object information operators */
      //@{
      gSize getCount() const;
      gString getReferenceName() const;
      //@}

      /** @name Subsetting function members and operators */
      //@{
      gVariations operator [] (const gArray<gArrayIndex> & varNum) const;
      //@}

      /** @name Variation adding function members */
      //@{
      gArrayIndex addInsertion(const gElement & element,gRelativePos elementPosition,const gSequence &seq);
      gArrayIndex addDeletion(const gElement & element,const gElementInterval & elementInterval);
      gArrayIndex addSubstitution(const gElement & element,const gElementInterval & elementInterval,const gSequence &seq);
      gArrayIndex addInsertion(const gString & referenceName,gPos referencePosition,const gSequence &seq);
      gArrayIndex addDeletion(const gString & referenceName,const gReferenceInterval & referenceInterval);
      gArrayIndex addSubstitution(const gString & referenceName,const gReferenceInterval & referenceInterval,const gSequence &seq);
      //@}

      /** @name Variation information function members */
      //@{
      gVariationType getVariationType(gArrayIndex varNum) const;
      gPos getVariationStart(gArrayIndex varNum) const;
      gPos getVariationEnd(gArrayIndex varNum) const;
      gReferenceInterval getVariationInterval(gArrayIndex varNum) const;
      gSequence getVariationSequence(gArrayIndex varNum) const;
      //@}
  };

  //A class to represent Genomics region as well as variations.
  //Used internally by the gElement class, you will never instantiate it.
  class gRegion:private gReferenceInterval,private gVariations {
      friend class gVariations;
      friend class gElement;
  private:

      //gString              i_referenceName;
      gSequenceRetriever   i_retriever;
      gElementSequenceMode i_mode;
      gSequence            i_sequence;
      gSize                i_refCount;

      gRegion(const gString & referenceName,const gReferenceInterval & referenceInterval,const gVariations & variations);
      gRegion(const gString & referenceName,const gReferenceInterval & referenceInterval,const gSequence & sequence,gElementSequenceMode mode,const gVariations & variations);
      gRegion(const gString & referenceName,const gReferenceInterval & referenceInterval,const gSequenceRetriever & sequenceRetriever,gElementSequenceMode mode,const gVariations & variations);
      gRegion(const gRegion & region,const gVariations & variations);
      ~gRegion();

      geco::gArray<gPos> getReferencePositionsFromRegion(const geco::gArray<gPos> & regionPositions, geco::gArray<gPos> & regionDisplacements) const;
      geco::gArray<gPos> getRegionPositionsFromReference(const geco::gArray<gPos> & referencePositions) const;

      gSequence getSequenceFromRetriever(const gReferenceInterval & regionInterval) const;
      gSequence getSequenceFromReference(const gReferenceInterval & regionInterval) const;
      gSequence getSequenceFromRegion(const gReferenceInterval & regionInterval) const;
      gSequence getSequence(const gReferenceInterval & regionInterval) const;

      gSequence variateSequence(const gSequence & seq, const gReferenceInterval & referenceInterval, const geco::gArray<gPos> & regionDisplacements) const;
  };

  template <class T> class gFeature;
  template <class T> class gFeatureRetriever;

  /** @brief Feature Retriever Implementation class
  *
  * Feature retriever implementation: this is an abstract template class from which
  * you can derive feature retrievers implementation by overloading teh virtual destructor
  * and the calculateFeature and clone memeber functions
  * implementation
  */
  template <class T> class gFeatureRetrieverImplementation:public gArrayRetrieverImplementation<T>{
    friend class gFeature<T>;
    friend class gFeatureRetriever<T>;
  protected:
      gSize i_startOffset;
      gSize i_endOffset;
      bool i_doubleStrand;


  public:
      /** @brief Constructor
      *
      * Instatiate an object
      * @param startOffset the number of bp needed for feature value calculation before the current position
      * @param endOffset the number of bp needed for feature value calculation after the current position
      * @param nFeatures the number of features this retriever is calculating.
      * @param doubleStrand whether this implementation calculates values that are valid ob both strands or not.
      */
      gFeatureRetrieverImplementation(gSize startOffset=0,gSize endOffset=0,gSize dimCount=1,bool doubleStrand=false);
      
      /** @brief Copy Constructor
      *
      * Instatiate an object
      * @param implementation tho object to copy from
      */
      gFeatureRetrieverImplementation(const gFeatureRetrieverImplementation<T> & implementation);
      
      /** @brief Virtaul destructor
      *
      * Deallocate the object
      */
      virtual ~gFeatureRetrieverImplementation();

      /** @brief Feature calculation in an interval
      *
      * This is a pure virtual function: it must be implemented when you derive a new retriever.
      * @param values the array that has to be recalculated
      * @param element the element
      * @param interval the interval
      * @return true if succesfull
      */
      virtual bool calculateFeature(gMatrix<T> & values,const gElement & element,const gElementInterval & interval) const = 0;

      /** @brief Start offset
      *
      * Returns the amount of positions before the current one needed to calculate the feature
      * @return a gSize value
      */
      gSize  getStartOffset() const;
      /** @brief End offset
      *
      * Returns the amount of positions before the current one needed to calculate the feature
      * @return a gSize value
      */
      gSize  getEndOffset() const;

      /** @brief Window length
      *
      * Returns the window comprising the start offset, the current position and the end offset
      * @return a gSize value
      */
      gSize  getWindowLenght() const;


      /** @brief Double Strandness
      *
      * Returns true if the feature doesn't need recalculation when reverting
      * @return a boolean value
      */
      bool isDoubleStrand() const;
  };


  /** @brief Feature Retriever Class
  *
  * Feature retriever: you can instantiate this template object providing a type and a suitable
  * implementation
  */
  template <class T> class gFeatureRetriever:public gArrayRetriever<T>{
    friend class gFeature<T>;
  public:
      /** @brief Default constructor
      *
      * Instatiante an empty gFEatureRetreiver
      */
      gFeatureRetriever();

      /** @brief Constructor
      *
      * gFEatureRetreiver constructor, you must specify the implementation
      * @param implementation the implementation fro this retriever
      */
      gFeatureRetriever(const gFeatureRetrieverImplementation<T> & implementation);

      /** @brief Destructor
      *
      * To be overloaded in derived classes if needed
      */
      virtual ~gFeatureRetriever();

      /** @brief Feautre calculation
      *
      * This member function is called internally by gElements whebn features needs to be calculate in an interval
      * it calls the appropriate implementation (prvided in the costructor of this object.
      * @param element The element upon which calculation should be made
      * @param interval The interval for which feature must be calculated
      * @param shift offset of the results
      * @return a gArray<T> object wih the results
      */
      gMatrix<T> getFeature(const gElement & element, const gElementInterval & interval, gRelativePos shift);

  };

  //This is used as a base class for gFeature
  class gFeatureBase{
    protected:
      gString i_name;
    public:
    gFeatureBase(const gString & name);
    gFeatureBase(const gFeatureBase & feature);
    virtual ~gFeatureBase();
    virtual gFeatureBase * clone(const gElement & newElement,const gElement & oldElement,const gElementInterval & interval,const gVariations & variations) const = 0;
    virtual bool calculate(const gElement & element,const gElementInterval &interval) = 0;
    const gString & getName() const;
  };

  //This is used internally in gElement
  template <class T> class gFeature: public gFeatureBase, public gMatrix<T>{
    friend class gElement;
  private:
    gFeatureRetriever<T>     i_retriever;
    gRelativePos             i_shift;
    gFeature(const gString & name, const gElement &element, const gFeatureRetriever<T> & retriever,gRelativePos shift);
    gFeature(const gFeature<T> & feature,const gElement & newElement,const gElement & oldElement,const gElementInterval & interval, const gVariations & variations);
    virtual ~gFeature();
    virtual bool calculate(const gElement & element,const gElementInterval & interval);
    virtual bool calculate(gMatrix<T> & values,const gElement & element,const gElementInterval & interval);
    virtual gFeatureBase * clone(const gElement & newElement,const gElement & oldElement,const gElementInterval & interval, const gVariations & variations) const;
  };

  /** @brief Genomic Element
  *
  * A class to easily manage genomic elements
  */
  class gElement {
      friend class gVariations;
  private:
      gRegion *                     i_region;
      bool                          i_forward;
      gSequence                     i_sequence;

      gReferenceInterval            i_boundaries;
      gArray<gPos>                  i_boundaries_reg_disp;

      gArray<gPos>                  i_site_list;
      gArray<gPos>                  i_site_list_reg_disp;

      gArray<gArrayIndex>           i_start_site_index;
      gArray<gArrayIndex>           i_end_site_index;

      std::vector< gFeatureBase * > i_features;

  protected:
      gRegion * getRegion() const;

      gArray<gPos> getRegionPositionsFromElement(const gArray<gRelativePos> & elementPositions) const;
      gArray<gRelativePos> getElementPositionsFromRegion(const gArray<gPos> & regionPositions) const;

      gArray<gPos> getReferencePositionsFromElement(const gArray<gRelativePos> & elementPositions,gArray<gPos> & displacements) const;
      gArray<gRelativePos> getElementPositionsFromReference( const gArray<gPos> & referencePositions,const gArray<gPos> & displacements) const;

      gReferenceInterval getReferenceIntervalFromElement(const gElementInterval & elementInterval,gArray<gPos> & displacements) const;
      gElementInterval getElementIntervalFromReference(const gReferenceInterval & referenceInterval, const gArray<gPos> & displacements) const;

      gSequence getSequenceFromRegion(const gElementInterval & interval) const;

      gMatrix<gRelativePos> getVariationElementIntervals() const;

      void setReferenceSites(gArray<gPos> referencePositions);
      void setElementSites(gArray<gRelativePos> elementPositions);

      void setConnectionIdxs(const gArray<gArrayIndex> &startIdxs,const gArray<gArrayIndex> &endIdxs);

  public:
      /** @name Constructors and destructor */
      //@{
      gElement();
      gElement(const gString & referenceName, const gReferenceInterval & referenceInterval, bool forward=true,const gVariations & variations=gVariations());
      gElement(const gString & referenceName, gPos referenceStart, gPos referenceEnd, bool forward=true,const gVariations & variations=gVariations());
      gElement(const gString & referenceName, const gReferenceInterval & referenceInterval, bool forward,const gSequence & referenceSequence,gElementSequenceMode mode=gElm,const gVariations & variations=gVariations());
      gElement(const gString & referenceName, gPos referenceStart, gPos referenceEnd, bool forward,const gSequence & referenceSequence,gElementSequenceMode mode=gElm,const gVariations & variations=gVariations());
      gElement(const gString & referenceName, const gReferenceInterval & referenceInterval, bool forward, const gSequenceRetriever & sequenceRetriever, gElementSequenceMode mode=gElm, const gVariations & variations = gVariations());
      gElement(const gString & referenceName, gPos referenceStart, gPos referenceEnd, bool forward, const gSequenceRetriever & sequenceRetriever, gElementSequenceMode mode=gElm, const gVariations & variations = gVariations());
      gElement(const gElement & element,const gElementInterval & elementInterval,bool forward=true,const gVariations & variations=gVariations());
      gElement(const gElement & element, gRelativePos elementStart = 0, gRelativePos elementEnd = 0, bool forward=true,const gVariations & variations = gVariations());
      ~gElement();
      //@}

      /** @name Assignment operators */
      //@{
      gElement & operator=(const gElement & elm);
      //@}

      /** @name Element information member functions */
      //@{
      bool isForward() const;
      gSize getLength() const;
      gString getReferenceName() const;
      gReferenceInterval getReferenceInterval() const;
      gElementInterval getElementInterval() const;
      //@}

      /** @name Element sequence member fuctions */
      //@{
      gSequence getSequence(const gElementInterval & elementInterval, gRelativePos startOffset=0, gRelativePos endOffset=0) const;
      gSequence getSequence(gRelativePos startOffset=0, gRelativePos endOffset=0) const;
      gSequence getSequence(const gArray<gRelativePos> & elementStarts,const gArray<gRelativePos> & elementEnds) const;
      //@}

      /** @name Positions & intervals conversion */
      //@{
      gArray<gPos> getReferencePositionsFromElement(const gArray<gRelativePos> & elementPositions) const;
      gArray<gRelativePos> getElementPositionsFromReference( const gArray<gPos> & referencePositions) const;
      gArray<gRelativePos> mapElementPositions(const gElement & element,const gArray<gRelativePos> & elementPositions) const;

      gReferenceInterval getReferenceIntervalFromElement(const gElementInterval & elementInterval, gRelativePos startOffset = 0, gRelativePos endOffset = 0) const;
      gElementInterval getElementIntervalFromReference(const gReferenceInterval & referenceInterval,gPos startOffset=0,gPos endOffset=0) const;
      gElementInterval mapElementInterval(const gElement & element,const gElementInterval & elementInterval) const;
      //gArray<gRelativePos> mapElementPositionsToInterval(const gElementInterval & elementInterval,const gArray<gRelativePos> & elementPositions) const;
      //@}

      /** @name Sites & connections */
      //@{
      gArray<gArrayIndex> insertElementSites(const gArray<gRelativePos> & elementPositions);
      gArray<gArrayIndex> insertReferenceSites(const gArray<gPos> & referencePositions);

      gArray<gArrayIndex> connectSites(const gArray<gArrayIndex> & siteIndex1,const gArray<gArrayIndex> & siteIndex2);

      gArrayIndex insertElementConnection(const gElementInterval & elementInterval,bool forward=true);
      gArray<gArrayIndex> insertElementConnections(const gArray<gRelativePos> & elementStarts,const gArray<gRelativePos> & elementEnds,bool forward=true);

      gArrayIndex insertReferenceConnection(const gReferenceInterval & referenceInterval,bool forward=true);
      gArray<gArrayIndex> insertReferenceConnections(const gArray<gPos> & referenceStarts,const gArray<gPos> & referenceEnds,bool forward=true);

      gSize getSitesCount() const;
      gSize getConnectionsCount() const;

      gArray<gArrayIndex> getElementSiteIdxs(const gArray<gRelativePos> & elementPositions) const;
      gArray<gArrayIndex> getReferenceSiteIdxs(const gArray<gPos> & referencePositions) const;

      gArray<gArrayIndex> getElementConnectionIdx(const gElementInterval & elementInterval) const;
      gArray<gArrayIndex> getReferenceConnectionIdx(const gReferenceInterval & referenceInterval,bool forward=true) const;

      gArray<gArrayIndex> getConnections(const gArray<gArrayIndex> &startIdxs,const gArray<gArrayIndex> &endIdxs) const;
      gArray<gArrayIndex> connectsTo(gArrayIndex site) const;
      gArray<gArrayIndex> connectsFrom(gArrayIndex site) const;

      gArray<gRelativePos> getSitesElementPosition(const gArray<gArrayIndex> & sites) const;
      gArray<gPos> getSitesReferencePosition(const gArray<gArrayIndex> & sites) const;

      gElementInterval getConnectionElementInterval(gArrayIndex startSiteIndex,gArrayIndex endSiteIndex) const;
      bool getConnectionElementStrand(gArrayIndex connectionIndex) const;
      gElementInterval getConnectionElementInterval(gArrayIndex connectionIndex) const;
      gReferenceInterval getConnectionReferenceInterval(gArrayIndex startSiteIndex,gArrayIndex endSiteIndex) const;
      bool getConnectionReferenceStrand(gArrayIndex connectionIndex) const;
      gReferenceInterval getConnectionReferenceInterval(gArrayIndex connectionIndex) const;
      gArray<gArrayIndex> getConnectionsStartSiteIndex(const gArray<gArrayIndex> & connectionIndexes) const;
      gArray<gArrayIndex> getConnectionsEndSiteIndex(const gArray<gArrayIndex> & connectionIndexes) const;
      gArray<gArrayIndex> mapElementConnections(const gElement & element, const gArray<gArrayIndex> & connectionIdxs);
      //@}


      /** @name Sub-elements*/
      //@{
      gElement getElement(const gElementInterval & elementInterval,bool forward=true, gRelativePos startOffset = 0, gRelativePos endOffset = 0) const;
      gElement getElement(gArrayIndex siteIndex1, gArrayIndex siteIndex2, gRelativePos startOffset, gRelativePos endOffset) const;
      gElement getElement(gArrayIndex connection, gRelativePos startOffset = 0, gRelativePos endOffset = 0) const;
      //gArray<gIndex> isContained(const gArray<gRelativePos> & elementPositions, const gElementInterval & elementInterval) const;
      //@}

      /** @name Features */
      //@{
      template <class T > gArrayIndex addFeature(const gString & name,const gFeatureRetriever<T> & retriever,gRelativePos shift=0);

      template <class T > gArray<T> getFeature(const gString & name, gSize dimNum, const gElementInterval & interval) const;
      template <class T > gArray<T> getFeature(gArrayIndex num, gSize dimNum,const gElementInterval & interval) const;
      template <class T > gMatrix<T> getFeature(const gString & name, const gElementInterval & interval, bool posbycolumns=true) const;
      template <class T > gMatrix<T> getFeature(gArrayIndex num, const gElementInterval & interval, bool posbycolumns=true) const;

      template <class T > gArray<T> getFeature(const gString & name, gSize dimNum, const gArray<gRelativePos> & positions) const;
      template <class T > gArray<T> getFeature(gArrayIndex num, gSize dimNum, const gArray<gRelativePos> & positions) const;
      template <class T > gMatrix<T> getFeature(const gString & name, const gArray<gRelativePos> & positions, bool posbyrow=true) const;
      template <class T > gMatrix<T> getFeature(gArrayIndex num,  const gArray<gRelativePos> & positions, bool posbyrow=true) const;

      template <class T > gArray<T> getFeatureVector(const gString & name, gPos position) const;
      template <class T > gArray<T> getFeatureVector(gArrayIndex num, gPos position) const;

      //@}
  };

  /** @brief Regions feature extractor parameters
  *
  * An helper class to contain parameters for gRegionFeaturesExtractor
  */
  class gRegionFeaturesParameters{
  public:
    typedef enum{none=0,otherSites=1,otherIntervals=2} overlapControl;      

  protected:
    gArray<gPos> i_sites;
    std::vector< gArray<gRelativePos> > i_siteOffsets;
    overlapControl i_control;
    gRelativePos i_shift;    
    gSize i_nwindows;
    
    gSize checkConnectionsAndOffset(const gArray<gArrayIndex> sites, const std::vector< gArray<gRelativePos> > & siteOffsets) const;    
  public:

    gRegionFeaturesParameters();
    gRegionFeaturesParameters(const gArray<gPos> & sites,const std::vector< gArray<gRelativePos> > & siteOffsets,overlapControl control=none,gRelativePos shift=0);
    gRegionFeaturesParameters(const gRegionFeaturesParameters & parameters);
    bool operator == (const gRegionFeaturesParameters & parameters) const;
    
    gArray<gPos> getSites() const;
    std::vector< gArray<gRelativePos> > getOffsets() const;      
    gArray<gRelativePos> getSiteOffsets(gPos siteNum) const;
    overlapControl getControl() const;
    gSize getWindowsCount() const;      
  };
  
  /** @brief Regions feature extractor
  *
  * An helper class to extract features from specific regions of an element
  */
  template<class T> class gRegionsFeaturesExtractor:public gRegionFeaturesParameters{
  public:
    
  private:
    gString i_featureName;
    gFeatureRetriever< T > i_retriever;
    
    gSize getIntervals( std::vector< std::pair< std::pair<gElementInterval,gRelativePos>, gArray<gRelativePos> > > & intervals,const gElement & element,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets) const;
  public:
    
    gRegionsFeaturesExtractor(const gFeatureRetriever<T> & retriever,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets, gRegionFeaturesParameters::overlapControl control=none, gSize shift=0);
    gRegionsFeaturesExtractor(const gFeatureRetriever<T> & retriever,const gRegionFeaturesParameters & parameters);
    gRegionsFeaturesExtractor(const gString & featureName,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets, gRegionFeaturesParameters::overlapControl control=none);
    gRegionsFeaturesExtractor(const gString & featureName,const gRegionFeaturesParameters & parameters);
    
    gMatrix<T> getRegionsFeature(const gElement & element,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets,bool posByColumns=true) const;
    gMatrix<T> getRegionsFeature(const gElement & element,bool posByColumns=true) const;
    gSequence getRegionsSequence(const gElement & element,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets) const;
    gSequence getRegionsSequence(const gElement & element) const;
  };

  /** @brief A class to manage gElement collections
  *
  * An helper class to manage gElements collections
  */
  class gElementCollection{
  private:
    std::vector< gElement > i_elements;
    gArray<gScore>          i_scores;
    std::vector<gString>    i_names;
    
  public:    
    gElementCollection();
    
    void insert(const gElement & e,gScore score, const gString & name);
    gElementCollection & append(const gElementCollection & elements);

    
    gSize getElementCount() const;
    gArray<gScore> getScores() const;
    gScore getScore(gPos id) const;
    std::vector<gString> getNames() const;
    gString getName(gPos id) const;
    
    const gElement & getElement(gPos id) const;
    gElement & getElement(gPos id);    
    

    gArray<gPos> intersect(const gString & chrom,const gReferenceInterval & interval,gSize minOverlap=0);
    gElementCollection select(const gArray<gPos> & positions) const;
    gElementCollection selectN(gSize N) const;
    gElementCollection selectNRandom(gSize N) const;
    gElementCollection selectNBest(gSize N) const;
    gElementCollection selectRange(gSize first,gSize last) const;
    
    void setScore(gPos i,gScore score);
    void setName(gPos i,const gString & name);
    void setElement(gPos i,const gElement & e);
    void setScores(const gArray<gScore> & scores);
    void setNames(const std::vector<gString> & names);
    
  };


  
  /** @brief Regions specification 
  *
  * An helper class to specify elements regions based on the connections
  */
  class gElementRegion{
  private:
    gArrayIndex i_connection;
    gShortUnsigned i_mode; //0=full,1=start,2=end,3=middle
    gRelativePos i_offset;
    gPos i_halfWindow;
  public:
    gElementRegion();      
    gElementRegion(gArrayIndex connection,gShortUnsigned mode,gRelativePos offset,gSize halfWindow);
    
    gElementInterval getInterval(const gElement & e) const;
    gSequence getSequence(const gElement & e) const;
  };
  
  typedef std::map<gSize,gElementRegion> gElementRegions;
  
//--------------------------------------------------------------------------
/// gFeature template class implementation
//--------------------------------------------------------------------------
template <class T> gFeature<T>::gFeature(const gString & name, const gElement & element, const gFeatureRetriever<T> & retriever,gRelativePos shift):gFeatureBase(name),i_retriever(retriever){
  const gFeatureRetrieverImplementation<T> & implementation = (const gFeatureRetrieverImplementation<T> & ) i_retriever.getImplementation();
  ((gMatrix<T> &) *this)=gMatrix<T>(element.getLength(),implementation.getDimCount(),0,true);
  i_shift=shift;
  calculate(element,element.getElementInterval());
}

template <class T> gFeature<T>::gFeature(const gFeature< T > & feature, const gElement & newElement, const gElement & oldElement, const gElementInterval & interval, const gVariations & variations):gFeatureBase(feature),i_retriever(feature.i_retriever){
  i_shift=feature.i_shift;
  const gFeatureRetrieverImplementation<T> & implementation = (gFeatureRetrieverImplementation<T> &) i_retriever.getImplementation();
  gElementInterval tint=newElement.getElementInterval();
  if(variations.getCount()>0){
    if(newElement.isForward()==oldElement.isForward()){
      //HERE WE SHOULD CALCULATE ONLY WHERE VARIATIONS INVALIDATE
      (gMatrix<T> &) *this=gMatrix<T>(tint.getLength(),implementation.getDimCount(),0,true);
      calculate(newElement,tint);
    }else{
      //HERE WE NEED To RECALCULATE EVERYTHING SO THIS IS THE RIGHT WAY
      (gMatrix<T> &) *this=gMatrix<T>(tint.getLength(),implementation.getDimCount(),0,true);
      calculate(newElement,tint);
    }
  }else{
    if(newElement.isForward()==oldElement.isForward()){
      (gMatrix<T> &) *this=gMatrix<T>(tint.getLength(),implementation.getDimCount(),0,true);
      gMatrix<T> & mat=*this;
      for(gSize i=0;i< implementation.getDimCount();i++){
       mat.setCol(i,gArray<T>(feature.getCol(i),interval.getStart(),interval.getEnd()));
      }
    }else{
      (gMatrix<T> &) *this=gMatrix<T>(tint.getLength(),implementation.getDimCount(),0,true);
      calculate(newElement,tint);
    }
  }
}

template <class T> gFeature<T>::~gFeature(){
}

template <class T> bool gFeature<T>::calculate(const gElement & element,const gElementInterval & interval){
  bool valid=true;
  if(interval.getLength()>0){
    const gFeatureRetrieverImplementation<T> & implementation = (const gFeatureRetrieverImplementation<T> & ) i_retriever.getImplementation();
    gSize so=implementation.getStartOffset();
    gSize eo=implementation.getEndOffset();
    gSize dsso=(implementation.isDoubleStrand())?(implementation.getEndOffset()):(0);
    gSize dseo=(implementation.isDoubleStrand())?(implementation.getStartOffset()):(0);

    gRelativePos totSOffset=so+dsso;
    gRelativePos totEOffset=eo+dseo;
    gRelativePos start=interval.getStart() - i_shift - totSOffset;
    gRelativePos end=interval.getEnd() - i_shift + totEOffset;
    //gRelativePos start=interval.getStart() - i_shift - implementation.getStartOffset();
    //gRelativePos end=interval.getEnd() - i_shift + implementation.getEndOffset();
    gElementInterval nintv(start,end);
    gMatrix<T> values(nintv.getLength(),implementation.getDimCount(),0,true);
    valid=implementation.calculateFeature(values,element,nintv);

    //IMPORTANT THIS SHOULD BE MADE EFFICIENT!
    if(valid){
      for(gSize i=0;i<implementation.getDimCount();i++){
        gArray<T> ocol=values.getCol(i);
        gArray<T> ncol(ocol,totSOffset,values.getRowsNum()-totEOffset);
        gMatrix< T >::setCol(i,ncol);
      }
    }

    //     This is a failed attempt check it better!      
//     if(valid){
//       gSize cc=values.getColsNum();
//       gSize rc=values.getRowsNum();
//       gMatrix<T> & m =*this;
//       m=gMatrix<T>(rc-totEOffset-totSOffset,cc,gArray<T>(values,cc*totSOffset,cc*(rc-totEOffset)));
//     }

    
  }
  return valid;
}

template <class T> bool gFeature<T>::calculate(gMatrix<T> & values,const gElement & element,const gElementInterval & interval){
  bool valid=true;
  if(interval.getLength()>0){
    const gFeatureRetrieverImplementation<T> & implementation = (const gFeatureRetrieverImplementation<T> & ) i_retriever.getImplementation();
    gSize so=implementation.getStartOffset();
    gSize eo=implementation.getEndOffset();
    gSize dsso=(implementation.isDoubleStrand())?(implementation.getEndOffset()):(0);
    gSize dseo=(implementation.isDoubleStrand())?(implementation.getStartOffset()):(0);

    gRelativePos totSOffset=so+dsso;
    gRelativePos totEOffset=eo+dseo;
    gRelativePos start=interval.getStart() - i_shift - totSOffset;
    gRelativePos end=interval.getEnd() - i_shift + totEOffset;
    //gRelativePos start=interval.getStart() - i_shift - implementation.getStartOffset();
    //gRelativePos end=interval.getEnd() - i_shift + implementation.getEndOffset();
    gElementInterval nintv(start,end);

    gMatrix<T> tcvalues(nintv.getLength(),implementation.getDimCount(),0,true);
    bool valid=implementation.calculateFeature(tcvalues,element,nintv);
    if(valid){
      for(gSize i=0;i<implementation.getDimCount();i++){
        //values.setCol(i,gArray<T>(tcvalues.getCol(i),implementation.getStartOffset(),tcvalues.getRowsNum()-implementation.getEndOffset()));
        values.setCol(i,gArray<T>(tcvalues.getCol(i),totSOffset,tcvalues.getRowsNum()-totEOffset));
      }
    }
    else throw gException("Unknown Error in feature calculation");
  }else values=gMatrix<T>();
  return valid;
}

template <class T> gFeatureBase * gFeature<T>::clone(const gElement & newElement, const gElement& oldElement,const gElementInterval & interval, const gVariations & variations) const {
  return new gFeature<T>(*this,newElement,oldElement,interval,variations);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/// gFeatureRetrieverImplementation template class definition
//--------------------------------------------------------------------------
template <class T> gFeatureRetrieverImplementation<T>::gFeatureRetrieverImplementation(gSize startOffset,gSize endOffset, gSize dimCount, bool doubleStrand):gArrayRetrieverImplementation<T>(dimCount){
  i_startOffset=startOffset;
  i_endOffset=endOffset;
  i_doubleStrand=doubleStrand;
}

template <class T> gFeatureRetrieverImplementation<T>::gFeatureRetrieverImplementation(const gFeatureRetrieverImplementation<T> & implementation):gArrayRetrieverImplementation<T>(implementation){
  i_startOffset=implementation.i_startOffset;
  i_endOffset=implementation.i_endOffset;
  i_doubleStrand=implementation.i_doubleStrand;
}


template <class T> gFeatureRetrieverImplementation<T>::~gFeatureRetrieverImplementation(){
}

template <class T> gSize gFeatureRetrieverImplementation<T>::getStartOffset() const{
  return i_startOffset;
}

template <class T> gSize gFeatureRetrieverImplementation<T>::getEndOffset() const{
  return i_endOffset;
}

template <class T> gSize gFeatureRetrieverImplementation<T>::getWindowLenght() const{
  return i_startOffset+i_endOffset+1;
}

template <class T> bool gFeatureRetrieverImplementation<T>::isDoubleStrand() const{
  return i_doubleStrand;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/// gFeatureRetriever template class definition
//--------------------------------------------------------------------------
template <class T> gFeatureRetriever<T>::gFeatureRetriever():gArrayRetriever<T>(){
}

template <class T> gFeatureRetriever<T>::gFeatureRetriever(const gFeatureRetrieverImplementation<T> & implementation):gArrayRetriever<T>(implementation){
}

template <class T> gFeatureRetriever<T>::~gFeatureRetriever(){
}

template <class T> gMatrix<T> gFeatureRetriever<T>::getFeature(const gElement & element, const gElementInterval & interval, gRelativePos shift){
  const gFeatureRetrieverImplementation<T> & implementation = (const gFeatureRetrieverImplementation<T> &) gRetriever::getImplementation();
  gRelativePos start=(gRelativePos) interval.getStart() - shift - (gRelativePos) implementation.getStartOffset();
  gRelativePos end=(gRelativePos) interval.getEnd() - shift + (gRelativePos) implementation.getEndOffset();
  gElementInterval nintv(start,end);
  gMatrix<T> values(nintv.getLength(),implementation.getDimCount(),0,true);
  bool valid=implementation.calculateFeature(values,element,nintv);
  if(valid){
    gMatrix<T> ret(nintv.getLength(),implementation.getDimCount(),0,true);
    for(gSize i=0;i<implementation.getDimCount();i++){
     ret.setCol(i,gArray<T>(values.getCol(i),implementation.getStartOffset(),values.getRowsNum()-implementation.getEndOffset()));
    }
    return ret;
  }else throw gException("getFeature: invalid result");
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/// gElement template function members implementation
//--------------------------------------------------------------------------
/** @brief Fature add template function member
 *
 * Add a feature to the current object.
 * The function is a template so you must specify the type of feature values.
 * You must also provide a name of the feature, a suitable feature retriver and an optional offset
 * The function returns the index of the newly inserted feature
 * @param name a gString containing a name fro this feature
 * @param retriever an instance of a suitable retriecver (i.e. same T)
 * @param shift an optional shift for the values positions with respect to the retriever calculations (defaults to 0)
 * @return a gArrayIndex value
 */
template <class T > gArrayIndex gElement::addFeature(const gString & name,const gFeatureRetriever<T> & retriever,gRelativePos shift){
  gArrayIndex ret=i_features.size();
  size_t i;
  for(i=0;i<i_features.size();i++) if(i_features[i]->getName()==name) break;
  if(i<i_features.size()){
    ret=i;
  }else{
    i_features.push_back(new gFeature<T>(name,*this,retriever,shift));
  }
  return ret;
}

/** @brief Fature get from name template function member
 *
 * Returns in a gArray<T> the values of a feature in a given element interval
 * The function is a template so you must specify the type for feature values.
 * If the type specified doesn't match the one of the feature
 * corresponding to the name provided then an exception is thrown.
 * @param name a gString containing a name fro this feature
 * @param dimNum the required dimension for the feature
 * @param interval the required element interval (i.e. same T)
 * @return a gArray<T> object
 *
*/
template <class T > gArray<T> gElement::getFeature(const gString & name,gSize dimNum, const gElementInterval & interval) const{
  size_t i;
  gArray<T> ret;
  for(i=0;i<i_features.size();i++) if(i_features[i]->getName()==name) break;
  if(i<i_features.size()){
    ret=getFeature<T>(i,dimNum,interval);
  }else throw gException("Required feature doesn't exist for this element");
  return ret;
}

/** @brief Fature get from index template function member
 *
 * Returns in a gArray<T> the values of a feature in a given element interval
 * The function is a template so you must specify the type for feature values.
 * If the type specified doesn't match the one of the feature
 * corresponding to the index provided then an exception is thrown.
 * @param num the index of the required feature
 * @param dimNum the required dimension
 * @param interval the required element interval (i.e. same T)
 * @return a gArray<T> object
 */
template <class T > gArray<T> gElement::getFeature(gArrayIndex num, gSize dimNum, const gElementInterval & interval) const{
  gFeature<T> * pFeat=dynamic_cast< gFeature<T> *> (i_features[num]);
  gArray<T> ret;
  if(interval.getLength()==0) return ret;
  if(pFeat){
    gElementInterval ei=getElementInterval();
    gElementInterval ii=ei.getIntersection(interval);
    gMatrix<T> ldvalues,rdvalues;
    gElementInterval ld=ei.getLeftDiff(interval);
    gElementInterval rd=ei.getRightDiff(interval);
    if((ld.getLength()>0) && (pFeat->calculate(ldvalues,*this,ld))){
      ret.concatenate(ldvalues.getCol(dimNum));
    }
    if(ii.getLength()>0){
      ret.concatenate(gArray<T>(((gMatrix<T> &)*pFeat).getCol(dimNum),(gPos) ii.getStart(),(gPos) ii.getEnd()));
    }
    if((rd.getLength()>0) &&(pFeat->calculate(rdvalues,*this,rd))){
      ret.concatenate(rdvalues.getCol(dimNum));
    }
  }else throw gException("getFeature: invalid template type");
  return ret;
}

template <class T > gMatrix<T> gElement::getFeature(const gString & name, const gElementInterval & interval,bool posbycolumns) const{
  size_t i;
  gMatrix<T> ret;
  for(i=0;i<i_features.size();i++) if(i_features[i]->getName()==name) break;
  if(i<i_features.size()){
    ret=getFeature<T>(i,interval,posbycolumns);
  }else throw gException("Required feature doesn't exist for this element");
  return ret;
}

template <class T > gMatrix<T> gElement::getFeature(gArrayIndex num, const gElementInterval & interval,bool posbycolumns) const{
  gFeature<T> * pFeat=dynamic_cast< gFeature<T> *> (i_features[num]);
  gMatrix<T> ret;
  if(interval.getLength()==0) return ret;
  gSize npos=interval.getLength();
  if(pFeat){
    gSize ndim=((gMatrix<T> &)*pFeat).getColsNum();
    gElementInterval ei=getElementInterval();
    gElementInterval ii=ei.getIntersection(interval);
    gMatrix<T> ldvalues,rdvalues;
    gElementInterval ld=ei.getLeftDiff(interval);
    gElementInterval rd=ei.getRightDiff(interval);
    if(!posbycolumns && ld.getLength()==0 && rd.getLength()==0){
      ret=((gMatrix<T> &)*pFeat).getRows(getArray<gPos>(ii.getStart(),ii.getEnd()-1,1));
    }else{
      gArray<T> arra(0,npos*ndim,false,false);
      if(posbycolumns){
	ret=gMatrix<T>(ndim,npos,arra);
      }else{
	ret=gMatrix<T>(npos,ndim,arra);
      }
      if(ld.getLength()>0) pFeat->calculate(ldvalues,*this,ld);
      if(rd.getLength()>0) pFeat->calculate(rdvalues,*this,rd);
      gArray<T> cdata(0,npos,false);
      for(gSize col=0;col<ndim;col++){
	//gArray<T> cdata;
	if(ld.getLength()>0) cdata.replace(0,ldvalues.getCol(col));//cdata.concatenate(ldvalues.getCol(col));
	if(ii.getLength()>0) cdata.replace(ld.getLength(),gArray<T>(((gMatrix<T> &)*pFeat).getCol(col),(gPos) ii.getStart(),(gPos) ii.getEnd()));//cdata.concatenate(gArray<T>(((gMatrix<T> &)*pFeat).getCol(col),(gPos) ii.getStart(),(gPos) ii.getEnd()));
	if(rd.getLength()>0) cdata.replace(ld.getLength()+ii.getLength(),rdvalues.getCol(col));//cdata.concatenate(rdvalues.getCol(col));
	if(posbycolumns){
	  ret.setRow(col,cdata);
	}else{
	  ret.setCol(col,cdata);
	}
      }
    }
  }else throw gException("getFeature: invalid template type");
  return ret;
}





template <class T > gArray<T> gElement::getFeature(const gString & name, gSize dimNum, const gArray<gRelativePos> & positions) const{
  gRelativePos minpos=positions.getMin()[0];
  gElementInterval elementInterval(minpos,positions.getMax()[0]+1);
  gArray<T> all=getFeature<T>(name,dimNum,elementInterval);
  gArray<T> ret=all[gArray<gPos>(positions-minpos)];
  return ret;
}

template <class T > gArray<T> gElement::getFeature(gArrayIndex num, gSize dimNum, const gArray<gRelativePos> & positions) const{
  gRelativePos minpos=positions.getMin()[0];
  gElementInterval intv=getElementInterval();
  gElementInterval elementInterval(minpos,positions.getMax()[0]+1);
  gArray<T> all=getFeature<T>(num,dimNum,elementInterval);
  gArray<T> ret=all[gArray<gPos>(positions-minpos)];
  return ret;
}

template <class T > gMatrix<T> gElement::getFeature(const gString & name, const gArray<gRelativePos> & positions, bool posbyrow) const{
  throw gException("Still to be implemented -> geco_element.h:804");
}

template <class T > gMatrix<T> gElement::getFeature(gArrayIndex num,  const gArray<gRelativePos> & positions, bool posbyrow) const{
  throw gException("Still to be implemented -> geco_element.h:804");
}

template <class T > gArray<T> gElement::getFeatureVector(const gString & name, gPos position) const{
  size_t i;
  gArray<T> ret;
  for(i=0;i<i_features.size();i++) if(i_features[i]->getName()==name) break;
  if(i<i_features.size()){
    ret=getFeatureVector<T>(i,position);
  }else throw gException("Required feature doesn't exist for this element");
  return ret;
}

template <class T > gArray<T> gElement::getFeatureVector(gArrayIndex num, gPos position) const{
  gFeature<T> * pFeat=dynamic_cast< gFeature<T> *> (i_features[num]);
  gArray<T> ret;
  if(pFeat){
    ret=pFeat->getRow(position);
  }else throw gException("getFeature: invalid template type");
  return ret;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/// gRegionsFeaturesExtractor2 template class definition
//--------------------------------------------------------------------------
template <class T> gRegionsFeaturesExtractor<T>::gRegionsFeaturesExtractor(const gFeatureRetriever<T> & retriever,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets, gRegionFeaturesParameters::overlapControl control, gSize shift):gRegionFeaturesParameters(sites,offsets,control,shift),i_retriever(retriever){
}

template <class T> gRegionsFeaturesExtractor<T>::gRegionsFeaturesExtractor(const gFeatureRetriever<T> & retriever,const gRegionFeaturesParameters & parameters):gRegionFeaturesParameters(parameters),i_retriever(retriever){
}

template <class T> gRegionsFeaturesExtractor<T>::gRegionsFeaturesExtractor(const gString & featureName,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets, gRegionFeaturesParameters::overlapControl control):gRegionFeaturesParameters(sites,offsets,control,0),i_featureName(featureName){
}

template <class T> gRegionsFeaturesExtractor<T>::gRegionsFeaturesExtractor(const gString & featureName,const gRegionFeaturesParameters & parameters):gRegionFeaturesParameters(parameters),i_featureName(featureName){
}

template <class T> gSize gRegionsFeaturesExtractor<T>::getIntervals( std::vector< std::pair< std::pair<gElementInterval,gRelativePos>, gArray<gRelativePos> > > & intervals,const gElement & element,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets) const{
  gSize nsites=sites.getSize();
  intervals=std::vector< std::pair< std::pair<gElementInterval,gRelativePos> , gArray<gRelativePos> > >(nsites);
  gArray<gRelativePos> sortedSitePositions=element.getSitesElementPosition(sites);
  gArray<gSize> sitesOrder=sortedSitePositions.sort(true);
  gSize matrixSize=0;
  
  for(gSize i=0;i<nsites;i++){
    intervals[sitesOrder[i]].second=offsets[sitesOrder[i]]+sortedSitePositions[i];
    intervals[sitesOrder[i]].first.first=gElementInterval(intervals[sitesOrder[i]].second.getMin()[0],intervals[sitesOrder[i]].second.getMax()[0]+1,true,true);
    intervals[sitesOrder[i]].first.second=matrixSize;
    matrixSize+=offsets[sitesOrder[i]].getSize();    
  }
  
//   std::cout << std::endl << "Expected positions" << std::endl;
//   for(gSize i=0;i<nsites;i++){
//     std::cout << intervals[sitesOrder[i]].first.first.getStart() << "-" << intervals[sitesOrder[i]].first.first.getEnd() << " nvar=" << intervals[sitesOrder[i]].second.getSize() << " offset: " << intervals[sitesOrder[i]].first.second << std::endl;
//   }
  
  
  //update intervals considering control over sites positions
  //std::cout << std::endl << "Site corrected Intervals" << std::endl;
  if(i_control&otherSites){
    for(gSize i=0;i<nsites;i++){
      gElementInterval & intv=intervals[sitesOrder[i]].first.first;
      
      if(intv.getLength()>0){
	if(i>0){
	  if(intv.getEnd()>sortedSitePositions[i-1]+1){
	    if(intv.getStart()<=sortedSitePositions[i-1]){
	      gSize a=which(intervals[sitesOrder[i]].first.first.contains(intervals[sitesOrder[i]].second)).getSize();
	      intervals[sitesOrder[i]].first.first=gElementInterval(sortedSitePositions[i-1]+1,intv.getEnd());
	      gSize b=which(intervals[sitesOrder[i]].first.first.contains(intervals[sitesOrder[i]].second)).getSize();
	      intervals[sitesOrder[i]].first.second+=a-b;
	    }
	  }else{
	    intervals[sitesOrder[i]].first.first=gElementInterval();
	  }
	}
      
        if(i<nsites-1){
	  if(intv.getStart()<sortedSitePositions[i+1]){
	    if(intv.getEnd()>sortedSitePositions[i+1]){
	      gSize a=which(intervals[sitesOrder[i]].first.first.contains(intervals[sitesOrder[i]].second)).getSize();
	      intervals[sitesOrder[i]].first.first=gElementInterval(intv.getStart(),sortedSitePositions[i+1]);
	      gSize b=which(intervals[sitesOrder[i]].first.first.contains(intervals[sitesOrder[i]].second)).getSize();
	      intervals[sitesOrder[i]].first.second+=a-b;
	    }
	  }else{
	    intervals[sitesOrder[i]].first.first=gElementInterval();
	  }
	}
      }
      intervals[sitesOrder[i]].second=intervals[sitesOrder[i]].second[which(intervals[sitesOrder[i]].first.first.contains(intervals[sitesOrder[i]].second))];
      //std::cout << intervals[sitesOrder[i]].first.first.getStart() << "-" << intervals[sitesOrder[i]].first.first.getEnd() << " nvar=" << intervals[sitesOrder[i]].second.getSize() << " offset: " << intervals[sitesOrder[i]].first.second << std::endl;
    }
    
  }
  
  //update intervals considering control over intervals positions and overlapping
  //std::cout << std::endl << "interval corrected Intervals" << std::endl;
  if(i_control&otherIntervals){
    for(gSize i=0;i<nsites-1;i++){
      gElementInterval & intv=intervals[sitesOrder[i]].first.first;
      gElementInterval & nextintv=intervals[sitesOrder[i+1]].first.first;
      if((intv.getLength()>0)&&(nextintv.getLength()>0)){
	if(intv.getStart()<nextintv.getStart()){
	  if(intv.getEnd()>nextintv.getStart()){
	    gRelativePos border=nextintv.getStart()+(intv.getEnd()-nextintv.getStart())/2;
	    intervals[sitesOrder[i]].first.first=gElementInterval(intv.getStart(),border);
	    gSize a=which(intervals[sitesOrder[i+1]].first.first.contains(intervals[sitesOrder[i+1]].second)).getSize();
	    intervals[sitesOrder[i+1]].first.first=gElementInterval(border,nextintv.getEnd());
	    gSize b=which(intervals[sitesOrder[i+1]].first.first.contains(intervals[sitesOrder[i+1]].second)).getSize();
	    intervals[sitesOrder[i+1]].first.second+=a-b;
	  }
	}else{
	  intervals[sitesOrder[i]].first.first=gElementInterval();
	  intervals[sitesOrder[i+1]].first.first=gElementInterval();
	}
      }
    }
    for(gSize i=0;i<nsites;i++){
      intervals[sitesOrder[i]].second=intervals[sitesOrder[i]].second[which(intervals[sitesOrder[i]].first.first.contains(intervals[sitesOrder[i]].second))];
      //std::cout << intervals[sitesOrder[i]].first.first.getStart() << "-" << intervals[sitesOrder[i]].first.first.getEnd() << " nvar=" << intervals[sitesOrder[i]].second.getSize() << " offset: " << intervals[sitesOrder[i]].first.second << std::endl;
    }
  }
  
  return matrixSize;
}

template <class T> gMatrix<T> gRegionsFeaturesExtractor<T>::getRegionsFeature(const gElement & element,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets,bool posByColumns) const{
  //checkConnectionsAndOffset(sites,offsets);
  std::vector< std::pair< std::pair<gElementInterval,gRelativePos> ,gArray<gRelativePos> > > intervals;
  gSize matsize=getIntervals(intervals,element,sites,offsets);
  gMatrix<T> result;
  if((intervals.size()==1) && (intervals[0].second.getSize()==1)){
    if(intervals[0].first.first.getLength()>0){
      if(i_featureName.getLength()>0){
	//gArray<gRelativePos> pos=intervals[i].second;
	//gMatrix<T> im=element.getFeature<T>(i_featureName,pos,true);
	gArray<gRelativePos> pos=intervals[0].second-intervals[0].first.first.getStart();
	if(posByColumns){
	  result=gMatrix<T>(i_retriever.getDimCount(),matsize,element.getFeature<T>(i_featureName,intervals[0].first.first,false));
	}else{
	  result=element.getFeature<T>(i_featureName,intervals[0].first.first,false);
	}
      }else{
	//std::cout <<  intervals[0].first.first.getStart() << "\t" << intervals[0].first.first.getEnd() << "\t" << intervals[0].first.first.getEnd()-intervals[0].first.first.getStart() << std::endl;
	gElement ne=element.getElement(intervals[0].first.first);
	gArray<gRelativePos> pos=intervals[0].second-intervals[0].first.first.getStart();
	gArrayIndex indx=ne.addFeature<T>("f",i_retriever,i_shift);
	if(posByColumns){
	  result=gMatrix<T>(i_retriever.getDimCount(),matsize,ne.getFeature<T>(indx,ne.getElementInterval(),false));
	}else{
	  result=ne.getFeature<T>(indx,ne.getElementInterval(),false);
	}
      }
    }
  }else{

    if(posByColumns){
      result=gMatrix<T>(i_retriever.getDimCount(),matsize,0,true);
    }else{
      result=gMatrix<T>(matsize,i_retriever.getDimCount(),0,true);
    }
    if(i_featureName.getLength()>0){
      for(gSize i=0;i<intervals.size();i++){
	if(intervals[i].first.first.getLength()>0){
	  gArray<gRelativePos> pos=intervals[i].second-intervals[i].first.first.getStart();
	  gMatrix<T> im=element.getFeature<T>(i_featureName,intervals[i].first.first,false);
	  if(posByColumns){
	    for(gSize c=0;c<pos.getSize();c++){
	      result.setCol(intervals[i].first.second+c,im.getRow(pos[c]));
	    }
	  }else{
	    for(gSize c=0;c<pos.getSize();c++){
	      result.setRow(intervals[i].first.second+c,im.getRow(pos[c]));
	    }
	  }
	}
      }   
      
    }else if(i_retriever.isImplemented()){
      for(gSize i=0;i<intervals.size();i++){
	if(intervals[i].first.first.getLength()>0){
	  gElement ne=element.getElement(intervals[i].first.first);
	  gArray<gRelativePos> pos=intervals[i].second-intervals[i].first.first.getStart();
	  gArrayIndex indx=ne.addFeature<T>("f",i_retriever,i_shift);
	  gMatrix<T> im=ne.getFeature<T>(indx,ne.getElementInterval(),false);
	  if(posByColumns){
	    for(gSize c=0;c<pos.getSize();c++){
	      result.setCol(intervals[i].first.second+c,im.getRow(pos[c]));
	    }
	  }else{
	    for(gSize c=0;c<pos.getSize();c++){
	      result.setRow(intervals[i].first.second+c,im.getRow(pos[c]));
	    }
	  }
	}
      }   
    }else throw gException("gRegionsFeaturesExtractor::getRegionsFeature: no suitable feature found");
  }
  
  return result;
}

template <class T> gMatrix<T> gRegionsFeaturesExtractor<T>::getRegionsFeature(const gElement & element,bool posByColumns) const{
  gMatrix<T> ret;
  if(i_sites.getSize()>0){
    ret=getRegionsFeature(element,i_sites,i_siteOffsets,posByColumns);
  }
  return ret;
}

template <class T> gSequence gRegionsFeaturesExtractor<T>::getRegionsSequence(const gElement & element,const gArray<gArrayIndex> & sites,const std::vector< gArray<gRelativePos> > & offsets) const{
  //checkConnectionsAndOffset(sites,offsets);  
  std::vector< std::pair< std::pair<gElementInterval,gRelativePos> ,gArray<gRelativePos> > > intervals;
  gSize matsize=getIntervals(intervals,element,sites,offsets);
  gSequence result(gArray<gChar>('N',matsize,true));
  for(gSize i=0;i<intervals.size();i++){
    if(intervals[i].first.first.getLength()>0){
      gSequence seq=element.getSequence(intervals[i].first.first);
      gArray<gSize> pos=intervals[i].second-intervals[i].first.first.getStart();
      result.setValues(getArray<gSize>(intervals[i].first.second,intervals[i].first.second+pos.getSize()-1,1),seq[pos]);
    }
  }
  return result;
}

template <class T> gSequence gRegionsFeaturesExtractor<T>::getRegionsSequence(const gElement & element) const{
  gSequence ret;
  if(i_sites.getSize()>0){
    ret=getRegionsSequence(element,i_sites,i_siteOffsets);
  }
  return ret;
}
//--------------------------------------------------------------------------
}

