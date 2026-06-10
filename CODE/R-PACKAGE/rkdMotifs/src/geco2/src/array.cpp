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
#include "array.h"
#include <string.h>


using namespace std;
using namespace geco;

//--------------------------------------------------------------------------
/// gString definition
//--------------------------------------------------------------------------
/** @brief Empty contructor
*
* Instatiate an empty gSequence object
*/
gString::gString() :gArray<gChar>() {
}

/** @brief Copy range contructor
*
* Instatiate a new gString initialized with the values from a
* gArray<gChar> in the range specified.
* @param str the string top copy from
* @param start the beginning position in the provided string (0 based, defaults to 0)
* @param end the ending position in the provided string (1 based, defaults to the last character position)
*/
gString::gString ( const geco::gArray< gChar >& str, gPos start, gPos end ) :gArray<gChar> ( str,start,end ) {
}

/** @brief c string range contructor
*
* Instatiate a gSequence object intializing with the specified range of the c string provided
* @param str a null terminated C string.
* @param start the beginning position in the provided string (0 based, defaults to 0)
* @param end the ending position in the provided string (1 based, defaults to the last character position)
*/
gString::gString ( const gChar * str ,gPos start,gPos end ) :gArray<gChar> ( str+start,((end==0)?(strlen(str)-start):(end-start)),gArray<gPos>() ) {
}

/** @brief std::string range contructor
*
* Instatiate a gSequence object intializing with the specified range of std:string provided
* @param str a std:string object
* @param start the beginning position in the provided string (0 based, defaults to 0)
* @param end the ending position in the provided string (1 based, defaults to the last character position)
*/
gString::gString ( const std::string & str , gPos start, gPos end ) :gArray<gChar> ( str.substr(start,((end==0)?(str.size()):(end))).c_str(),((end==0)?(str.size()):(end))-start,gArray<gPos>() ) {
}

/** @brief Same charater initializing constructor
 *
 * Instatiate a gString object of the specified initialized with a repeated character
 * @param length the desired length
 * @param symbol the initializing char.
 */
gString::gString ( gSize length,gChar symbol ) :gArray<gChar> ( symbol,length,false ) {
}

/** @brief std::string cast operator
 *
 * returns a std:string object initialized with this string
 */
gString::operator const std::string () const {
    string ret;
    if ( getSize() >0 ) {
        ret=string ( getData(),getSize() );
    } else {
        ret=string ( "" );
    }
    return ret;
}

/** @brief Length getter
 *
 * return the length of this gString
 * @return length of this object
 */
gSize gString::getLength() const {
    return getSize();
}

/** @brief Lower case
 *
 * Transform all this string characters to lower case.
 * @return a reference to this
 */
gString & gString::lower() {
  gArray<gPos> tochange=which(((gArray<gChar> &)(*this))<97);
  if(tochange.getSize()>0){
   setValues(tochange,(*this)[tochange]+32);
  }
  return *this;
}

/** @brief Upper case
 *
 * Transforms all of this string characters to upper case.
 * @return a reference to this
 */
gString & gString::upper() {
  gArray<gPos> tochange=which(((gArray<gChar> &)(*this))>90);
  if(tochange.getSize()>0){
   setValues(tochange,(*this)[tochange]-32);
  }
    //( gArray<gChar>& ) ( *this )-=32;
    return *this;
}

/** @brief std::string Assignment operator
 *
 * Copies the std:string provided into this
 * @param str the string to copy from
 * @return a reference to this
 */
gString & gString::operator = ( const string & str ) {
    setData ( str.size(),str.c_str(),gArray<gPos>(),false );
    return *this;
}

/** @brief c string assignment operator
 *
 * Copies the provided c string into this
 * @param str the string to copy from
 * @return a reference to this
 */
gString & gString::operator = ( const gChar *str ) {
    setData ( strlen ( str ),str,gArray<gPos>(),false );
    return *this;
}

/** @brief String concatenation operator
 *
 * Add the character from the provided string to the end of this
 * @param str the string to copy from
 * @return a reference to this
 */
gString & gString::operator += ( const gString & str ) {
    concatenate ( str );
    return *this;
}

/** @brief String concatenation operator
 *
 * return a new gString object composed by the concatenation o0f this
 * and the provided one
 * @param str the string to concatenate to this
 * @return gString newly instatiated obejct
 */
gString gString::operator + ( const gString & str ) const {
    gString ret ( *this );
    ret+=str;
    return ret;
}

/** @brief String splitting
 *
 * Fill the provided std:vector<gString> with tokens from this
 * separated by the provided character.
 * The character itsef is not included
 * @param splitChar The char used to split this
 * @param parts the vector object to fill woth tokens
 * @return a reference to the passed vector
 */
std::vector<gString> & gString::split ( gChar splitChar, vector< gString >& parts, bool skipEmpty ) const{
    gSize pos=0,last=0;
    if(getLength()>0){
      do {
        if ( (*this)[ pos ] ==splitChar ) {
          if ( pos>last ) parts.push_back ( gString ( *this,last,pos ) );
          else if(!skipEmpty) parts.push_back ( gString ( "" ) );
          last=pos+1;
        }
        pos++;
      } while ( pos<getSize() );
      if ( pos>last ) parts.push_back ( gString ( *this,last,pos ) );
      else if(!skipEmpty) parts.push_back ( gString ( "" ) );
    }
    return parts;
}

/** @brief Lower case
 *
 * Create a copy of this with all character set to lower case
 * @return gString newly instatiated obejct
 */
gString gString::getLower() const {
    gString ret=*this;
    ret.lower();
    return ret;
}

/** @brief Upper case
 *
 * Create a copy of this with all character set to upper case
 * @return gString newly instatiated obejct
 */
gString gString::getUpper() const {
    gString ret=*this;
    ret.upper();
    return ret;
}

/** @brief Equal operator
 *
 * Return true if the provided string is identical to this.
 * @param str The string to compare (c-string, gString or string)
 * @return true or false
 */
bool gString::operator == ( const char * str ) const {
    return ( strncmp ( getData(),str,getSize() ) ==0 ) && ( strlen ( str ) ==getSize() );
}

bool gString::operator== ( const gString & str ) const {
    return ( strncmp ( getData(),str.getData(),min ( getSize(),str.getSize() ) ) ==0 ) && ( getSize() ==str.getSize() );
}

bool gString::operator== ( const string & str ) const {
    return *this==str.c_str();
}

/** @brief Lesser than operator
 *
 * Return true if the provided string is lesser than this.
 * @param str The string to compare (c-string, gString or string)
 * @return true or false
 */
bool gString::operator < (const char * str) const{
  gPos l2=strlen(str);
  int diff=strncmp(getData(),str,min(getSize(),l2));
  return (diff<0) || (diff==0 && (getSize()<l2));
}

bool gString::operator < (const gString & str) const{
  int diff=strncmp(getData(),str.getData(),min(getSize(),str.getSize()));
  return (diff<0) || (diff==0 && (getSize()<str.getSize()));
}

bool gString::operator < (const std::string & str) const{
  return *this<str.c_str();
}



/** @brief Output streaming operator
 *
 * Write this to the provided output stream
 * @param stream The output stream to write to
 * @return a reference to the stream
 */
ostream & gString::out ( ostream & stream ) const {

    gChar *tp=new gChar[getSize() +1];
    memcpy(tp,getData(),getSize()*sizeof(gChar));
    gArray<gBool> nab=isNA();
    gArray<gPos> nap=which(nab);
    for(gPos i=0;i<nap.getSize();i++){
      tp[nap[i]]='N';
    }
    tp[getSize()]=0;

    stream << string(tp) ;
    delete [] tp;
    //stream << ((string) *this);
    return stream;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// non-template helper functions definition
//--------------------------------------------------------------------------
/** @brief Condition function
 *
 * @param condition
 * @param includeNA
 * @return
 * @ingroup arrays
 */
gArray<gPos> geco::which ( const gArray<gBool> & condition,bool includeNA) {
  return condition.find ( (gBool)1,0,0,!includeNA);
}

/** @brief Samplig function
 * returns an array of n random values between min and max 
 * @param n
 * @param min
 * @param max
 * @return an array of n values between min and max
 * @ingroup arrays
 */
gArray<gPos> geco::randomSample(gSize n, gPos min, gPos max){
  //srand (time(NULL));
  gSize s=max-min+1;
  gSize * list=new gSize[s];
  memcpy(list,getArray<gSize>(min,max,1).getConstData(),sizeof(gSize)*s);
  gSize * vektor=new gSize[n];
  for (gSize i = 0; i < n; i++) {
    int j = i + rand() % (s - i);
    int temp = list[i];
    list[i] = list[j];
    list[j] = temp;
    vektor[i] = list[i];
  }  
  delete [] list;
  gArray<gPos> ret;
  ret.ownData(n,vektor);
  return ret;
}

/** @brief String output stream operator
 *
 * String Simple output to a stream
 * @param sout the output stream
 * @param str the object to be streamed out
 * @return a reference to the output stream (for chainig outputs)
 * @ingroup arrays
 */
ostream & geco::operator << ( ostream &sout, const gString & str ) {
    return str.out ( sout );
}

/** @brief String input from stream operator
 *
 * String Simple intput from a stream
 * @param sin the input stream
 * @param str the object to be streamed in
 * @return a reference to the input stream (for chainig outputs)
 * @ingroup arrays
 */
istream & geco::operator >> (istream & sin, gString & str) {
  string ss;
  sin >> ss;
  str=ss;
  return sin;
}


//--------------------------------------------------------------------------
