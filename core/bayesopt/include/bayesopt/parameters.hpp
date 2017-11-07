/**  \file parameters.hpp \brief Parameter definitions. */
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


#ifndef __BOPT_PARAMETERS_HPP__
#define __BOPT_PARAMETERS_HPP__

#include <string>
#include <boost/numeric/ublas/vector.hpp>

#include "bayesopt/parameters.h"     // learning_type, score_type

typedef boost::numeric::ublas::vector<double>                   vectord;

/**
 * Namespace of the library interface
 */
namespace bayesopt {    
    class KernelParameters{
    public:
        /* Class member variables */
        std::string  name;          /**< Name of the kernel function */
        vectord hp_mean;            /**< Kernel hyperparameters prior (mean, log space) */
        vectord hp_std;             /**< Kernel hyperparameters prior (st dev, log space) */
    
        /* Class member functions */
        KernelParameters();
    };
    
    class MeanParameters{
    public:
        /* Class member variables */
        std::string  name;          /**< Name of the mean function */
        vectord coef_mean;          /**< Basis function coefficients (mean) */
        vectord coef_std;           /**< Basis function coefficients (std) */
        
        /* Class member functions */
        MeanParameters();
    };
    
    class Parameters{
    public:
        /* 
         * Class members variables
         */
        size_t n_iterations;        /**< Maximum BayesOpt evaluations (budget) */
        size_t n_inner_iterations;  /**< Maximum inner optimizer evaluations */
        size_t n_init_samples;      /**< Number of samples before optimization */
        size_t n_iter_relearn;      /**< Number of samples before relearn kernel */

        /** Sampling method for initial set 1-LHS, 2-Sobol (if available),
         *  other value-uniformly distributed */
        size_t init_method;          
        int random_seed;            /**< >=0 -> Fixed seed, <0 -> Time based (variable). */    

        int verbose_level;          /**< Neg-Error,0-Warning,1-Info,2-Debug -> stdout
                                        3-Error,4-Warning,5-Info,>5-Debug -> logfile */
        std::string log_filename;   /**< Log file path (if applicable) */

        size_t load_save_flag;      /**< 1-Load data,2-Save data,
                                             3-Load and save data. */
        std::string load_filename;  /**< Init data file path (if applicable) */
        std::string save_filename;  /**< Sava data file path (if applicable) */

        std::string surr_name;      /**< Name of the surrogate function */
        double sigma_s;             /**< Signal variance (if known). 
                                        Used in GaussianProcess and GaussianProcessNormal */
        double noise;               /**< Variance of observation noise (and nugget) */

        double alpha;               /**< Inverse Gamma prior for signal var. 
                                        Used in StudentTProcessNIG */
        double beta;                /**< Inverse Gamma prior for signal var. 
                                        Used in StudentTProcessNIG */

        score_type sc_type;         /**< Score type for kernel hyperparameters (ML,MAP,etc) */
        learning_type l_type;       /**< Type of learning for the kernel params */
        bool l_all;                 /**< Learn all hyperparameters or only kernel */

        double epsilon;             /**< For epsilon-greedy exploration */
        size_t force_jump;          /**< If >0, and the difference between two 
                                        consecutive observations is pure noise, 
                                        for n consecutive steps, force a random 
                                        jump. Avoid getting stuck if model is bad 
                                        and there is few data, however, it might 
                                        reduce the accuracy. */

        KernelParameters kernel;    /**< Kernel parameters */
        MeanParameters mean;        /**< Mean (parametric function) parameters */  

        std::string crit_name;      /**< Name of the criterion */
        vectord crit_params;        /**< Criterion hyperparameters (if needed) */
        
        /*
         * Class member functions 
         */
        Parameters();
        
        /* Constructor to get values from bopt_params */
        Parameters(bopt_params c_params);
        
        /* Generates a bopt_params struct from a Parameters instance */
        bopt_params generate_bopt_params();
        
        /* Allows to change the l_type with a string */
        void set_learning(std::string name);
        /* Returns the string that corresponds to the l_type value */
        std::string get_learning();
        
        /* Allows to change the sc_type with a string */
        void set_score(std::string name);
        /* Returns the string that corresponds to the sc_type value */
        std::string get_score();
        
    private:
        /* Encapsulated default values assigment operations */
        void init_default();

      void bostrdup (char* d, const char *s);
    }; //class Parameters
} //namespace bayesopt


#endif
