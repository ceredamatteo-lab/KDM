#pragma once

/**
 * @file gPWM.h
 * @author  Uberto Pozzoli  <uberto.pozzoli@bp.lnf.it>*
 * @version 0.1
 *
 * @section LICENSE
 * Copyright (C) 2016 by Uberto Pozzoli
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
 * ....
 */
#include "element.h"

#include <vector>
#include <map>


namespace geco{
  namespace kmers{
    
    class gPWM{
    private:  
      gMatrix<gScore> i_PWM;
      gMatrix<gScore> i_PSSM;
      gArray<gScore> i_background;
    protected:
      
    public:
      typedef enum {PWMScore,PSSMScore} ScoringMode;
      gPWM();
      gPWM(const gMatrix<gScore> & PWMcounts, gScore N, const gArray<gScore> & background=gArray<gScore>(0.25,4,false),gScore pcounts=1e-9);
      gPWM(const gMatrix<gScore> & PWMdata, const gArray<gScore> & background=gArray<gScore>(0.25,4,false),gScore pcounts=1e-9,gScore tolerance=1e-5);

      void setData(const gMatrix<gScore> & PWMcounts, gScore N,const gArray<gScore> & background,gScore pcounts);      
      void setData(const gMatrix<gScore> & PWMdata, const gArray<gScore> & background,gScore pcounts,gScore tolerance);

      gArray<gScore> getScores(const gSequence seq,ScoringMode mode) const;
      gArray<gScore> getScores(const gElement & element,const gElementRegion & region ,bool doublestrand,ScoringMode mode) const;
      gArray<gScore> getScores(const gElement & element,bool doublestrand,ScoringMode mode) const;
      std::pair< gArray<gScore>,gArray<bool> > getDsScores(const gElement & element,ScoringMode mode) const;
      
      gScore getMaxScore(const gSequence seq,bool doublestrand,ScoringMode mode) const;
      gArray<gScore> getMaxScore(const gElement & element,const gElementRegions & regions ,bool doublestrand,ScoringMode mode) const;
      gArray<gScore> getMaxScore(const gElementCollection & elementCollection,const gElementRegions & regions ,bool doublestrand,ScoringMode mode) const;
      
      gScore getAvgScore(const gSequence seq,bool doublestrand,ScoringMode mode) const;
      gArray<gScore> getAvgScore(const gElement & element,const gElementRegions & regions ,bool doublestrand,ScoringMode mode) const;
      gArray<gScore> getAvgScore(const gElementCollection & elementCollection,const gElementRegions & regions ,bool doublestrand,ScoringMode mode) const;
      
      std::pair<gScore,gScore> getScoreRange() const;
      gArray<gScore> getPosIC() const;
      gScore getIC() const;
      
      gArray<gScore> computePvalues(const gArray<gScore> & scores,gScore initialGranularity=1e-4,bool forcedGranularity=true,gScore maxGranularity=1e-9) const;
      std::map<gScore,gScore> computeDistribution(const gArray<gScore> & range=gArray<gScore>(0,2,true),gScore granularity=1e-4) const;
      gScore computeThreshold(gScore alpha,gScore initialGranularity=1e-4,bool forcedGranularity=true,gScore maxGranularity=1e-10) const;
      
      gSize getLength() const;
      gMatrix<gScore> getPWM() const;
      gMatrix<gScore> getPSSM() const;
      gArray<gScore> getBackground() const;
      
    };

    typedef std::map<gScore,gScore> gPWMScoreMap;

    class gPWMScoreDistribution: public gPWMScoreMap{
    private:
    protected:
    public:
      gPWMScoreDistribution();
      gPWMScoreDistribution(const gPWM & pwm,gScore minAlpha,gScore granularity);
      gPWMScoreDistribution(const gMatrix<gScore> & dmat);
      gArray<gScore> pvalue(gArray<gScore> score);
      //gScore score(gScore pvalue);
      
      gMatrix<gScore> getDistributionAsMatrix() const;
    };

    class gPWMSet: public std::vector<std::pair<gString,gPWM> >{
    private:
      std::map<gString,gSize> i_names_map;
    public:
      gPWMSet();
      void addMatrix(const gPWM & matrix,const gString & name);
      const gPWM & operator [] (const gString & name) const;
      const gPWM & operator [] (gSize index) const;
      const gString & getName(gSize index) const;
      gPos getIndex(const gString & name) const;
      
      gMatrix<gScore> getScores(const gSequence seq,gPWM::ScoringMode mode) const;
      
      gArray<gScore> getMaxScore(const gSequence seq,bool doublestrand,gPWM::ScoringMode mode) const;
      gArray<gScore> getAvgScore(const gSequence seq,bool doublestrand,gPWM::ScoringMode mode) const;
      
      gMatrix<gScore> getMaxScore(const std::vector<std::string> seqs,bool doublestrand,gPWM::ScoringMode mode) const;
      gMatrix<gScore> getAvgScore(const std::vector<std::string> seqs,bool doublestrand,gPWM::ScoringMode mode) const;
      
      gMatrix<gScore> getMaxScore(const std::vector<gSequence> seqs,bool doublestrand,gPWM::ScoringMode mode) const;
      gMatrix<gScore> getAvgScore(const std::vector<gSequence> seqs,bool doublestrand,gPWM::ScoringMode mode) const;
      
      
      
      gMatrix<gScore> getMaxFeatures(const gElementCollection & elementCollection,const gElementRegions & regions,bool doublestrand,gPWM::ScoringMode mode) const;
      gMatrix<gScore> getAvgFeatures(const gElementCollection & elementCollection,const gElementRegions & regions,bool doublestrand,gPWM::ScoringMode mode) const;
      
      
    };

    class gPWMSetScoreDistributions:public std::vector<gPWMScoreDistribution>{
    private:
    protected:
    public:
      gPWMSetScoreDistributions(const gPWMSet & set,gScore minAlpha,gScore granularity);	
    };
    
    gPWMSet readFromMemeFile(const gString & memeFileName,const gArray<gScore> & background=gArray<gScore>(),gScore pcounts=1e-9,gScore tolerance=1e-5);
  }
}
