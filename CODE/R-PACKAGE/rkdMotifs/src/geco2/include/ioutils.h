#pragma once
#include "csv_iterator.hpp"

#include <array.h>
#include <map>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>




namespace geco{
	namespace utils{
		
		class gInputFileStream:public boost::iostreams::filtering_istream {
		public:
			gInputFileStream();
		};

		bool gzopen ( const gString & filename,gInputFileStream & stream );
		
		gString getFileName(const gString & filePath, bool withExtension = true);
	
		class gLineParser{
		public:
			typedef enum{Skip,Numeric,Integer,Boolean,String,Character} columnType;
		private:
			std::vector<columnType> i_colType;		
			std::vector<boost::any> i_values;
			std::string i_sep;		
			std::map<gString,gShortUnsigned> i_colMap;
			utils::gInputFileStream i_fstream;
			bool i_header;
			bool i_open;
			bool i_lineOk;
			gSize i_nLine;
			gString i_message;
			gString i_fileName;
			std::string i_currentLine;
			boost::empty_token_policy i_empty_tokens_policy;

		public:
			gLineParser(const std::vector<columnType> & colTypes, const std::vector<gString> & colNames,const gString & separators=gString(" \t"),bool dropEmptyTokens=true, bool header=false):i_colType(colTypes),i_values(colTypes.size()),i_sep(separators){
				i_header=header;
				i_lineOk=false;
				i_open=false;
				i_nLine=0;
				
				if(colNames.size()!=colTypes.size()) throw gException("");
				for(gSize i=0;i<i_colType.size();i++){
					i_colMap.insert(std::make_pair(colNames[i],i));
				}
				if(dropEmptyTokens){
					i_empty_tokens_policy=boost::drop_empty_tokens;
				}else{
					i_empty_tokens_policy=boost::keep_empty_tokens;
				}
			}

			
			inline const gString & getError(){
				return i_message;
			}
			
			bool open(const gString & fileName){
				i_fileName=fileName;
				i_fstream.reset();
				gzopen(fileName,i_fstream);
				i_nLine=0;
				i_open=true;
				if(i_header){
					std::getline(i_fstream,i_currentLine);
					//check for compatibility
			// 		if(i_fstream.good()){
			// 			char_separator<char> sep(((string)separators).c_str(),"", keep_empty_tokens);	
			// 			typedef tokenizer<char_separator<char> > strtokenizer;
			// 			strtokenizer tokens(line, sep);
			// 			if(tokens.size()<colNames.size()) throw gException(  ((string)( gString("gLineParser: file ")+fileName+gString(" header has less columns than required colNames"))).c_str());
			// 			strtokenizer::iterator tok_iter=tokens.begin();
			// 			for (gSize i=0;i<colNames.size();i++){
			// 				cout << *tok_iter << "\t" colNames[i] << endl;
			// 				i_colMap.insert(make_pair(colNames[i],i));
			// 				tok_iter++;
			// 			}			
			// 			
			// 		}
				}
				return i_open;
			}
			
			
			inline bool nextRow(){
				i_lineOk=false;
				if(i_fstream.is_complete() & !i_fstream.eof() & i_fstream.good()){
					getline(i_fstream,i_currentLine);
					if(i_fstream.good()){
						boost::tokenizer<boost::char_separator<char>>  tokens(i_currentLine, boost::char_separator<char>(i_sep.c_str(),"", i_empty_tokens_policy));
						boost::tokenizer<boost::char_separator<char>>::iterator tok_iter=tokens.begin();
						i_nLine++;
						for (gSize i=0;i<i_colType.size();i++){
							if(tok_iter!=tokens.end()){
								i_lineOk=true;
								try{
									switch(i_colType[i]){
										case      Skip:	break;
										case   Numeric: i_values[i]=boost::lexical_cast<gScore>(*tok_iter);break;
										case   Integer: i_values[i]=(gSize) boost::lexical_cast<gScore>(*tok_iter);break;
										case   Boolean: i_values[i]=boost::lexical_cast<bool>(*tok_iter);break;
										case    String:	i_values[i]=gString(boost::lexical_cast<std::string>(*tok_iter));break;
										case Character: i_values[i]=boost::lexical_cast<char>(*tok_iter);break;
									}
								}catch(std::exception & e){
									i_message=gString(e.what()) + gString(" reading line ") + gString(std::to_string(i_nLine)) + gString(" column: ") + gString(std::to_string(i+1)) + gString(" file: ") +i_fileName;
									throw gException(((std::string)i_message).c_str());
								}catch(...){
									i_message=gString("Unknown error") + gString(" reading line ") + gString(std::to_string(i_nLine)) + gString(" column: ") + gString(std::to_string(i+1)) + gString(" file: ") +i_fileName;
									throw gException(((std::string)i_message).c_str());
								}
								tok_iter++;
							}else{
								
							}
						}
						
					}else{
						
					}
				}
				return i_lineOk;
			}
			
			template<class T>
			inline const T & get(const gString & colName) const{
				return boost::any_cast<const T &>(i_values[i_colMap.at(colName)]);
			}
			
		};

        template <typename entryType>
        class gCsvList:public std::list<entryType>{
        public:
            gCsvList(){
            }
            
            gCsvList(const std::string & fileName){
                typedef typename entryType::i_type myType;
                geco::utils::gInputFileStream in;
                geco::utils::gzopen(fileName,in);
                csv::iterator<myType> it(in);
                for_each(it,csv::iterator<myType>(),[this](const myType & rec){this->push_back(entryType(rec));});
            }
        };
		
	}
	
}
