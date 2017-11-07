/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/

#include "fileparser.hpp"

namespace bayesopt
{   
  namespace utils 
  {
    namespace ublas = boost::numeric::ublas;

    FileParser::FileParser(std::string filename, int prec)
    : filename(filename), input(), output(){
        setPrecision(prec);
    }
    
    // For using to_value and to_string functions
    FileParser::FileParser(int prec){
        setPrecision(prec);
    }

    FileParser::~FileParser(){
        close();
    }
    
    void FileParser::open(bool readMode){
        if(readMode){
            openInput();
        }
        else{
            openOutput();
        }
    }
    void FileParser::openOutput(){
        close();
        output.open(filename.c_str());
    }
    void FileParser::openInput(){
        close();
        input.open(filename.c_str());
    }
    void FileParser::close(){
        output.close();
        input.close();
        currentLine = "";
    }
    /* Check fileparser mode */
    bool FileParser::isReading(){
        return input.is_open();
    }
    
    bool FileParser::isWriting(){
        return output.is_open();
    }
    
    /* Changes output precision of real numbers */
    void FileParser::setPrecision(int prec){
        if(prec > 0){
            precision = prec;
        }
        else{ // Default precision
            precision = 10;
        }
    }
    
    /* Checks if a file exists */
    bool FileParser::fileExists(){
        std::ifstream ifile(filename.c_str());
        bool result = ifile;
        ifile.close();
        return result;
    }
    
    /* Data write/read function */
    void FileParser::write(std::string name, std::string value){
        output << name << "=" << value << std::endl;
    }
    void FileParser::read(std::string name, std::string &value){
        //TODO: (Javier) use movePointer() returned bool, should be
        //used to avoid &value = "" when not found
        if(!movePointer(name,value))
	  {
            std::cerr << "Variable: " << name 
		      << " does not exist in file: " << filename 
		      << std::endl;
	  }
    }
    std::string FileParser::read(std::string name){
        std::string ret;
        read(name, ret);
        return ret;
    }
        
    /* Array write/read function */
    void FileParser::write(std::string name, 
			   const std::vector<std::string> &arr, 
			   const std::vector<int> &dims){
        // Write dimensions
        output << name << "=[";
        for(std::vector<int>::const_iterator it = dims.begin(); it != dims.end(); ++it) {
            if(it != dims.begin()){
                output << ",";
            }
            output << *it;
        }
        output << "]";
        
        // Write array
        output << "(";
        for(std::vector<std::string>::const_iterator it = arr.begin(); it != arr.end(); ++it) {
            if(it != arr.begin()){
                output << ",";
            }
            output << *it;
        }
        output << ")" << std::endl;
    }
    void FileParser::read(std::string name, std::vector<std::string> &arr, 
			  std::vector<int> &dims)
    {
        std::string contents;
        if(movePointer(name,contents)){
            parseArray(contents, arr, dims);
        }
        else{
            std::cerr << "Variable: " << name 
		      << " does not exist in file: " << filename 
		      << std::endl;
        }
    }
    
    /* 
     * Non-templated functions (types that requires special treatment) 
     */
    void FileParser::write_chars(std::string name, char* value){
        std::string str(value);
        write(name,str);
    }
    void FileParser::read_chars(std::string name, char* value){
        std::string str;
        read(name, str); 
        strcpy(value, str.c_str());
    }
    
    void FileParser::readOrWrite(std::string name, char* value){
        if(isReading()){
            read_chars(name, value);
        }
        else if(isWriting()){
            write_chars(name,value);
        }
    }
    
    void FileParser::write_ublas(std::string name, 
				 ublas::vector<double> &values)
    {
        std::vector<int> dims;
        dims.push_back(values.size());
        
        std::vector<std::string> arr;
        for(ublas::vector<double>::iterator it = values.begin(); 
	    it != values.end(); ++it)
	  {
            arr.push_back(to_string(*it));
	  }
        write(name, arr, dims);
    }
    void FileParser::read_ublas(std::string name, ublas::vector<double> &values){
        std::vector<std::string> arr;
        std::vector<int> dims;
        read(name, arr, dims);
        
        std::vector<double> doubles_arr;
        for(std::vector<std::string>::iterator it = arr.begin(); it != arr.end(); ++it) {
            doubles_arr.push_back(to_value<double>(*it));
        }
        
        values.resize(arr.size(), false);
        std::copy(doubles_arr.begin(), doubles_arr.end(), values.begin());        
    }
    void FileParser::readOrWrite(std::string name, ublas::vector<double> &values){
        if(isReading()){
            read_ublas(name, values);
        }
        else if(isWriting()){
            write_ublas(name,values);
        }
    }
    
    void FileParser::write_vecOfvec(std::string name, 
				    std::vector<ublas::vector<double> > &values){
        std::vector<int> dims;
        dims.push_back(values.size());
        dims.push_back(values.at(0).size());
        
        std::vector<std::string> arr;
        for(size_t i=0; i<values.size(); i++){
            ublas::vector<double> current = values.at(i);
            for(ublas::vector<double>::iterator it = current.begin(); 
		it != current.end(); ++it) {
                arr.push_back(to_string(*it));
            }
        }
        write(name, arr, dims);
    }
    void FileParser::read_vecOfvec(std::string name, std::vector<ublas::vector<double> > &values){
        std::vector<int> dims;
        std::vector<std::string> arr;
        read(name, arr, dims);
        
        size_t sample_dim = dims.at(1);
        
        values.resize(dims.at(0));
        for(size_t i=0; i<dims.at(0); i++){
            values.at(i).resize(sample_dim);
            for(size_t j=0; j<sample_dim; j++){
                values.at(i)[j] = to_value<double>(arr.at(i*sample_dim + j));
            }
        }
    }
    void FileParser::readOrWrite(std::string name, std::vector<ublas::vector<double> > &values){
        if(isReading()){
            read_vecOfvec(name, values);
        }
        else if(isWriting()){
            write_vecOfvec(name,values);
        }
    }
    
    void FileParser::write_double_array(std::string name, double values[], size_t length){
        std::vector<std::string> arr;
        std::vector<int> dims;
        dims.push_back(length);
        for(size_t i=0; i<length; i++){
            arr.push_back(to_string(values[i]));
        }
        write(name, arr, dims);
    }
    void FileParser::read_double_array(std::string name, double values[], size_t length){
        std::vector<std::string> arr;
        std::vector<int> dims;
        read(name, arr, dims);
        
        for(size_t i=0; i<length; i++){
            values[i] = to_value<double>(arr.at(i));
        }
    }
    void FileParser::readOrWrite(std::string name, double values[], size_t length){
        if(isReading()){
            read_double_array(name, values, length);
        }
        else if(isWriting()){
            write_double_array(name,values, length);
        }
    }
    
    /* Search variables in file */
    bool FileParser::movePointer(std::string name, std::string &contents){
        if(currentLine.length() > 0 && startsWith(currentLine, name+"=")){
            contents = currentLine.substr(name.length()+1);
            return true;
        }
        
        // TODO (Javier): detect when all the lines were read instead of 2 iterations over file
        // Wrap the file around in the search of a variable
        for(int i=0; i<2; i++){
            while (getline( input, currentLine)){
                if(currentLine.length() > 0 && startsWith(currentLine, name+"=")){
                    contents = currentLine.substr(name.length()+1);
                    return true;
                }
            }
            input.clear();
            input.seekg(0, std::ios::beg);
        }
        contents = "";
        return false;
        
    }
    
    bool FileParser::startsWith(std::string all, std::string sub){
        return all.rfind(sub, 0) == 0;
    }
    
    void FileParser::parseArray(std::string contents, 
				std::vector<std::string> &arr, 
				std::vector<int> &dims){
        // Parse data
        size_t init_arr = contents.find("(")+1;
        size_t end_arr = contents.find(")")-1;
        std::string input = contents.substr(init_arr, end_arr - init_arr +1);
        utils::split(input, ',', arr);        
        
        // Parse dimensions
        std::vector<std::string> dims_string;
        size_t init_dims = contents.find("[")+1;
        size_t end_dims = contents.find("]")-1;
        input = contents.substr(init_dims, end_dims - init_dims +1);
        utils::split(input, ',', dims_string);
        
        for(std::vector<std::string>::iterator it = dims_string.begin(); 
	    it != dims_string.end(); ++it) {
            dims.push_back(to_value<size_t>(*it));
        }
    }

  } //namespace utils
} //namespace bayesopt

