#include "ioutils.h"
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/filesystem.hpp>



using namespace std;
using namespace geco;
using namespace geco::utils;

gInputFileStream::gInputFileStream() :boost::iostreams::filtering_istream() {
	
}

bool geco::utils::gzopen(const gString & filename,utils::gInputFileStream & stream){
	using namespace boost::iostreams;
	bool ret=true;
	string fname=filename;
	std::ifstream in(fname.c_str());
	if(!in.is_open()) throw gException( ((std::string)(gString("Error opening file ")+filename)).c_str());
	in.close();
	stream.reset();
	char c;
	filtering_istream gzin;
	gzin.push(gzip_decompressor());
	gzin.push(file_source(fname.c_str()),std::ios_base::in | std::ios_base::binary);
	try{
		gzin >> c;
		if(gzin.good()){
			stream.push(gzip_decompressor());
		}
	}catch(...){
	 //do nothing	
	}
    stream.push(file_source(fname.c_str()),std::ios_base::in | std::ios_base::binary);
	return in.is_open();
}

gString geco::utils::getFileName(const gString & filePath, bool withExtension){
	// Create a Path object from File Path
	boost::filesystem::path pathObj(filePath);
 
	// Check if file name is required without extension
	if(withExtension == false)
	{
		// Check if file has stem i.e. filename without extension
		if(pathObj.has_stem())
		{
			// return the stem (file name without extension) from path object
			return pathObj.stem().string();
		}
		return "";
	}
	else
	{
		// return the file name with extension from path object
		return pathObj.filename().string();
	}
 
}
