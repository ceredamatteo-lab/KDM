#pragma once
/**
 * @file geco_transcript.h
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
 * This file contains the declaration of the gTranscript class.
 */
#include "element.h"

namespace geco {


/** @brief Trascript Class
 *
 * A spcialization of gElement to manage transcripts
 */
class gTranscript:public gElement {
private:
    gArray<gArrayIndex>         i_cds;
    gArray<gArrayIndex>         i_transcript;
    gArray<gArrayIndex>         i_exons;
    gArray<gArrayIndex>         i_introns;
    gArray<gTranscriptSiteType> i_sitetype;

    bool initInternal(const gArray<gPos> & refExonStarts,const gArray<gPos> & refExonEnds,bool forward,gPos cdsStart,gPos cdsEnd);

protected:
public:
    /** @name Constructors and destructor */
    //@{
    gTranscript();
    gTranscript(const gString & refName,const gArray<gPos> & refExonStarts,const gArray<gPos> & refExonEnds,bool forward,gPos startOffset=0,gPos endOffset=0,gPos cdsStart=0,gPos cdsEnd=0);
    gTranscript(const gString & refName,const gArray<gPos> & refExonStarts,const gArray<gPos> & refExonEnds,bool forward,const gSequenceRetriever & retriever,gElementSequenceMode mode,gPos startOffset=0,gPos endOffset=0,gPos cdsStart=0,gPos cdsEnd=0);
    gTranscript(const gString & refName,const gArray<gPos> & refExonStarts,const gArray<gPos> & refExonEnds,bool forward,const gSequence & sequence,gElementSequenceMode mode,gPos startOffset=0,gPos endOffset=0,gPos cdsStart=0,gPos cdsEnd=0);
    gTranscript(const gTranscript & transcript,const gVariations & variations=gVariations());
    //@}

    /** @name Assignment operators */
    //@{
    gTranscript & operator = (const gTranscript & transcript);
    //@}

    /** @name pre-mRNA information */
    //@{
    gSize getExonCount() const;
    gSize getIntronCount() const;

    gSize getPremRNALength() const;
    gSize getExonLength(gArrayIndex exonNum) const;
    gSize getIntronLength(gArrayIndex intronNum) const;

    gElementInterval getPremRNAInterval() const;
    gElementInterval getExonInterval(gArrayIndex exonNum) const;
    gElementInterval getIntronInterval(gArrayIndex intronNum) const;
    gElementInterval get5UTRInterval() const;
    gElementInterval getCdsInterval() const;
    gElementInterval get3UTRInterval() const;
    //@}

    /** @name mRNA information */
    //@{
    gArray<gSize> getmRNALength() const;
    gArray<gSize> get5UTRLength() const;
    gArray<gSize> getCdsLength() const;
    gArray<gSize> get3UTRLength() const;

    gArray<gShortUnsigned> getExonPhase(gArrayIndex exonNum) const;

    gSequence getmRNASequence(gPos start=0,gPos end=0) const;
    gSequence get5UTRSequence(gPos start=0,gPos end=0) const;
    gSequence getCdsSequence(gPos start=0,gPos end=0) const;
    gSequence get3UTRSequence(gPos start=0,gPos end=0) const;

    gArray<gRelativePos> mapElementPositionsOnmRNA(const gArray<gRelativePos> & elementPositions) const;
    gArray<gRelativePos> mapElementPositionsOn5UTR(const gArray<gRelativePos> & elementPositions) const;
    gArray<gRelativePos> mapElementPositionsOnCds(const gArray<gRelativePos> & elementPositions) const;
    gArray<gRelativePos> mapElementPositionsOn3UTR(const gArray<gRelativePos> & elementPositions) const;

    gArray<gRelativePos> mapCdsPositionsOnElement(const gArray<gRelativePos> & cdsPositions) const;

    gElementInterval mapElementIntervalOnmRNA(const gElementInterval & elementInterval) const;
    gElementInterval mapElementIntervalOn5UTR(const gElementInterval & elementInterval) const;
    gElementInterval mapElementIntervalOnCds(const gElementInterval & elementInterval) const;
    gElementInterval mapElementIntervalOn3UTR(const gElementInterval & elementInterval) const;


    gArray<gBool> ismRNA(gArray<gPos> elementPositions) const;
    gArray<gBool> is5UTR(gArray<gPos> elementPositions) const;
    gArray<gBool> isCds(gArray<gPos> elementPositions) const;
    gArray<gBool> is3UTR(gArray<gPos> elementPositions) const;
    gArray<gBool> isIntron(gArray<gPos> elementPositions) const;
    gArray<gBool> isExon(gArray<gPos> elementPositions) const;
    //@}

    /** @name Protein information*/
    //@{
    gArray<gSize> getProteinLength() const;
    gString getProteinSequence(const gElementInterval & interval) const;
    gArray<gPos> mapElementPositionsOnProtein(const gArray<gPos> & elementPositions) const;
    gElementInterval mapElementIntervalOnProtein(const gElementInterval & elementInterval) const;
    gArray<gRelativePos> mapProteinPositionOnElement(const gArray<gPos> & proteinPositions) const;
    //@}
};





}

