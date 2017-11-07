#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include "fileparser.hpp"
#include "ublas_extra.hpp"


void printDoubleArray(double* arr, size_t length, bayesopt::utils::FileParser &fp);
void printUblasVec(boost::numeric::ublas::vector<double> &ublas_vec, bayesopt::utils::FileParser &fp);

int main()  
{
    std::string filename = "test_fileparser.dat";
    int prec = 9; // max required precision in this test is 8, so one more digit is required
    
    size_t var_size_t = 24;
    int var_int_signed = -32;   // This 2 variables are necessary to check that FileParser
    int var_int = 25;           // does not retrieve "int_signed" when searching for only "int"
    double var_double = 1.2345678;

    size_t char_length = 5;
    char var_char_arr[5] = {'c', 'h', 'a', 'r'};
    std::string var_string = "string_value";
    
    size_t arr_length = 5;
    double var_double_arr[5] = {0,1.1,3.333,5.55555,8.88888888};
    boost::numeric::ublas::vector<double> var_ublas_vec = bayesopt::utils::array2vector(var_double_arr, arr_length);
    
    bayesopt::utils::FileParser fp(filename);
    fp.setPrecision(prec);
    // Write values
    fp.openOutput();
    
    // Note that variable names are the type, this is to easily check which type generated that value
    fp.write("size_t", var_size_t);
    fp.write("int_signed", var_int_signed);
    fp.write("int", var_int);
    fp.write("double", var_double);
    fp.write_chars("char_arr", var_char_arr);
    fp.write("string", var_string);
    fp.write_double_array("double_arr", var_double_arr,arr_length);    
    fp.write_ublas("ublas_vec", var_ublas_vec);
    
    fp.close();
    
    // Create new variables
    size_t cp_size_t;
    int cp_int_signed;
    int cp_int;
    double cp_double;
    char cp_char_arr[char_length];
    std::string cp_string;
    double cp_double_arr[arr_length];
    boost::numeric::ublas::vector<double> cp_ublas_vec;
    
    // Read values
    fp.openInput();
 
    // Note that variable names are the type, this is to easily check which type generated that value   
    fp.read("size_t", cp_size_t);
    fp.read("int_signed", cp_int_signed);
    fp.read("int", cp_int);
    fp.read("double", cp_double);
    fp.read_chars("char_arr", cp_char_arr);
    fp.read("string", cp_string);
    fp.read_double_array("double_arr", cp_double_arr,arr_length);
    fp.read_ublas("ublas_vec", cp_ublas_vec);
    
    fp.close();
    
    // Display .dat file
    std::cout << "Displaying saved file contents:" << std::endl;
    std::ifstream input;
    input.open(filename.c_str());
    std::string line;
    while(getline(input,line)){
        std::cout << line << std::endl;
    }
    input.close();
    std::cout << std::endl;
    
    // Display variable contents
    std::cout.precision(prec);
    std::cout << "Displaying variable contents:" << std::endl;
    std::cout << "size_t=" << cp_size_t << std::endl;
    std::cout << "int_signed=" << cp_int_signed << std::endl;
    std::cout << "int=" << cp_int << std::endl;
    std::cout << "double=" << cp_double << std::endl;
    std::cout << "char_arr=" << std::string(cp_char_arr) << std::endl;
    std::cout << "string=" << cp_string << std::endl;
    std::cout << "double_arr="; printDoubleArray(cp_double_arr, arr_length, fp);
    std::cout << "ublas_vec="; printUblasVec(cp_ublas_vec, fp);
    std::cout << std::endl;
    
    // Try to remove used .dat file
    if( remove( filename.c_str() ) == 0 ){
        std::cout << "File \"" << filename << "\" successfully removed" << std::endl;
    }
    else{
        std::cout << "Error: cannot remove \"" << filename << "\" file" << std::endl; 
    }
    
    /*
     * Compare restored values with original values
     */
    int returnValue = 0;
    // Testing size_t
    if(var_size_t != cp_size_t){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_size_t << " to be equals to " << cp_size_t << std::endl;
    }
    // Testing signed integer
    if(var_int_signed != cp_int_signed){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_int_signed << " to be equals to " << cp_int_signed << std::endl;
    }
    // Testing integer
    if(var_int != cp_int){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_int << " to be equals to " << cp_int << std::endl;
    }
    // Testing double
    if(var_double != cp_double){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_double << " to be equals to " << cp_double << std::endl;
            std::cout << '\t' << "(Note: check values, could be a rounding precision error)" << std::endl;
    }
    // Testing string
    if(var_string != cp_string){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_string << " to be equals to " << cp_string << std::endl;
    }
    // Testing double array and ublas vector
    for(int i=0; i<arr_length; i++){
        if(var_double_arr[i] != cp_double_arr[i]){
            returnValue = -1;
            std::cout << "ERROR: double array at " << i << " expected" << var_double_arr[i]
                << " to be equals to " << cp_double_arr[i] << std::endl;
            std::cout << '\t' << "(Note: check values, could be a rounding precision error)" << std::endl;
        }
        
        if(var_ublas_vec[i] != cp_ublas_vec[i]){
            returnValue = -1;
            std::cout << "ERROR: ublas vector at " << i << " expected" << var_ublas_vec[i]
                << " to be equals to " << cp_ublas_vec[i] << std::endl;
            std::cout << '\t' << "(Note: check values, could be a rounding precision error)" << std::endl;
        }
    }
    
    if(returnValue == 0){
        std::cout << "Tests completed without errors" << std::endl;
    }
    return returnValue;
}

void printDoubleArray(double* arr, size_t length, bayesopt::utils::FileParser &fp){
    std::cout << "[" << length << "](";
    for(size_t i=0; i<length; i++){
        if(i>0){
            std::cout << ",";
        }
        std::cout << fp.to_string(arr[i]);
    }
    std::cout << ")" << std::endl;
}

void printUblasVec(boost::numeric::ublas::vector<double> &ublas_vec, bayesopt::utils::FileParser &fp){
    std::cout << "[" << ublas_vec.size() << "](";
    for(boost::numeric::ublas::vector<double>::iterator it = ublas_vec.begin(); it != ublas_vec.end(); ++it){
        if(it != ublas_vec.begin()){
            std::cout << ",";
        }
        std::cout << fp.to_string(*it);
    }
    std::cout << ")" << std::endl;
}
