 /**
 * @file 2bit.h
 * @author  Uberto Pozzoli  <uberto.pozzoli@bp.lnf.it>
 * @version 0.1
 *
 * @section LICENSE
 * Copyright (C) 2011 by Uberto Pozzoli
 * uberto.pozzoli@bp.lnf.it
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the* GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * @section DESCRIPTION
 * This file includes all the features of the dbseq module of teh geco library
 */
#ifndef __GECO_GENOMES_RETRIEVERS_2BIT__
#define __GECO_GENOMES_RETRIEVERS_2BIT__
#include "config.h"
#include "element.h"

namespace geco{
  namespace sequence_retrievers{

    class gLocal2bitSequenceRetriever:public gSequenceRetrieverImplementation{
    private:
      gString i_filename;
      std::vector<gString> i_references;
      gArray<gSize> i_dnaSize;
      gArray<gPos> i_offset;
      std::vector< gArray<gPos> > i_nBlockStarts;
      std::vector< gArray<gPos> > i_nBlockEnds;
      std::vector< gArray<gPos> > i_maskBlockStarts;
      std::vector< gArray<gPos> > i_maskBlockEnds;
      virtual gArrayRetrieverImplementation<gChar> * clone() const;
      virtual gSequence getSequence_Internal(const gString & reference,gPos start,gPos end) const;
    protected:
    public:
      gLocal2bitSequenceRetriever(const gString & filename,bool useMask);
      gLocal2bitSequenceRetriever(const gLocal2bitSequenceRetriever & retriever);
      virtual ~gLocal2bitSequenceRetriever();
      gSequence getRandomSequence(gSize length) const;
      gSequence getRandomSequence(const gString & reference,gSize length) const;
      gReferenceInterval getRandomInterval(gSize length) const;
      gReferenceInterval getRandomInterval(const gString & reference,gSize length) const;
      gSize getReferenceLength(const gString & reference) const;
      gSize getReferenceLength(gSize refIndex) const;
      gSize getReferenceCount() const;
      gString getReferenceName(gSize refIndex) const;
    };

  }
}

#endif
