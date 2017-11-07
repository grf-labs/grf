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

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <boost/numeric/ublas/assignment.hpp>
#include "bayesopt/bayesopt.hpp"
#include "param_loader.hpp"
#include "fileparser.hpp"

#ifndef M_PI
/* It shouldn't be necessary, but Windows is completely nuts about
   math constants and I HATE when include order matters. */
    #define M_PI       3.14159265358979323846
#endif

class TemplateWritter{
    private:
    std::ifstream input;
    std::string prefix;
    
    public:
    TemplateWritter(std::string &template_filename, std::string &template_prefix)
    :input(), prefix(template_prefix)
    {
        input.open(template_filename.c_str());
    }
    
    ~TemplateWritter(){
        input.close();
    }
    
    // Create a file overwriting parameters from the template file
    void createFile(std::string &filename, const vectord& xin){
        // Move input pointer to the start of the file and create the output stream
        input.clear();
        input.seekg(0, std::ios::beg);
        std::ofstream output(filename.c_str());
        
        bayesopt::utils::FileParser fp; // to_value() and to_string() functions access
        
        std::string currentLine;
        while (getline( input, currentLine)){
            if(currentLine.length() > 0){
                // Search for prefixes
                size_t prefix_index = currentLine.find(prefix);
                while(prefix_index != std::string::npos){
                    size_t number_index = prefix_index + prefix.length();
                    
                    // Get the number after the prefix
                    std::string number_string = "";
                    while(number_index < currentLine.size() && std::isdigit(currentLine.at(number_index))){
                        number_string = number_string + currentLine.at(number_index);
                        number_index++;
                    }

                    // Index of the parameter
                    size_t x_index = fp.to_value<size_t>(number_string);
                    
                    // Check xin bounds before attempting to read data
                    if(x_index >= xin.size()){
                        // Show error and print "OUT_OF_BOUNDS" on the position of the template
                        std::cout << "ERROR: parameter " << x_index << " is out of bounds" << std::endl;
                        output << currentLine.substr(0, prefix_index) << "OUT_OF_BOUNDS";
                    }
                    else{
                        // And finally, write the value on the correct position of the template
                        double x = xin(x_index);
                        output << currentLine.substr(0, prefix_index) << fp.to_string(x);
                    }
                    
                    // Change currentLine to hold unwritten data and perform a new search
                    currentLine = currentLine.substr(number_index);
                    prefix_index = currentLine.find(prefix);
                }
                // Write unwritten data
                output << currentLine << std::endl;    
            }
        }
        
        output.close();
    }
};

class XMLCallsBranin: public bayesopt::ContinuousModel
{
public:
  XMLCallsBranin(bayesopt::Parameters par):
    ContinuousModel(2,par) {}

  double evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }

    float y = -1;
    double x1 = xin(0);
    double x2 = xin(1);
        
    // Create XML from template
    std::string template_filename("../examples/standalone_calls/problem_template.xml");
    std::string template_prefix("XXXX_");
    TemplateWritter tw(template_filename, template_prefix); 
    std::string created_filename("../examples/standalone_calls/created_file.xml");
    tw.createFile(created_filename, xin);
    
    // Results filename
    std::string results_filename("../examples/standalone_calls/results.txt");

    // Call python script
    std::string call = "python ../examples/standalone_calls/eval_branin_xml.py " + created_filename + " " + results_filename;
    system(call.c_str());
    
    // TODO (Javier): Change results to XML format
    bayesopt::utils::FileParser fp(results_filename.c_str());
    fp.openInput();
    fp.read("y",y);
    
    return y;
  }

  bool checkReachability(const vectord &query)
  {return true;};

  inline double sqr( double x ){ return x*x; };

  void printOptimal()
  {
    vectord sv(2);  
    sv(0) = 0.1238938; sv(1) = 0.818333;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.5427728; sv(1) = 0.151667;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.961652; sv(1) = 0.1650;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
  }

};

int main(int nargs, char *args[])
{  
  bayesopt::Parameters par;
  if(nargs > 1){
    if(!bayesopt::utils::ParamLoader::load(args[1], par)){
        std::cout << "ERROR: provided file \"" << args[1] << "\" does not exist" << std::endl;
        return -1;
    }
  }
  else{
    par = initialize_parameters_to_default();
    par.n_iterations = 190;
    par.random_seed = 0;
    par.verbose_level = 1;
    par.noise = 1e-10;
    //bayesopt::utils::ParamLoader::save("system_opt.txt", par);
  }

    
  XMLCallsBranin branin(par);
  vectord result(2);

  branin.optimize(result);
  std::cout << "Result: " << result << "->" 
	    << branin.evaluateSample(result) << std::endl;
  branin.printOptimal();
  
    // Remove results.txt file
  std::string filename("../examples/standalone_calls/created_file.xml");
  if( remove( filename.c_str() ) == 0 ){
    std::cout << "File \"" << filename << "\" successfully removed" << std::endl;
  }
  else{
    std::cout << "Error: cannot remove \"" << filename << "\" file" << std::endl; 
  }
  
    // Remove results.txt file
  filename = "../examples/standalone_calls/results.txt";
  if( remove( filename.c_str() ) == 0 ){
    std::cout << "File \"" << filename << "\" successfully removed" << std::endl;
  }
  else{
    std::cout << "Error: cannot remove \"" << filename << "\" file" << std::endl; 
  }
  return 0;
}


