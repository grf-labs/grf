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

#include "param_loader.hpp"

namespace bayesopt
{       
  namespace utils{      
    bool ParamLoader::load(std::string filename, Parameters &par){
        utils::FileParser fp(filename);
        if(!fp.fileExists()){
            return false;
        }
        
        par = initialize_parameters_to_default();
        fp.openInput();
        
        loadOrSave(fp, par);
        return true;
    }
    
    void ParamLoader::save(std::string filename, Parameters &par){
        utils::FileParser fp(filename);
        fp.openOutput();
        
        loadOrSave(fp, par);
    }
    
    void ParamLoader::loadOrSave(utils::FileParser &fp, Parameters &par){
        fp.readOrWrite("n_iterations", par.n_iterations);
        fp.readOrWrite("n_inner_iterations", par.n_inner_iterations);
        fp.readOrWrite("n_init_samples", par.n_init_samples);
        fp.readOrWrite("n_iter_relearn", par.n_iter_relearn);
        fp.readOrWrite("init_method", par.init_method);
        fp.readOrWrite("random_seed", par.random_seed);
        fp.readOrWrite("verbose_level", par.verbose_level);
        fp.readOrWrite("log_filename", par.log_filename);
        fp.readOrWrite("load_save_flag", par.load_save_flag);
        fp.readOrWrite("load_filename", par.load_filename);
        fp.readOrWrite("save_filename", par.save_filename);
        fp.readOrWrite("surr_name", par.surr_name);
        fp.readOrWrite("sigma_s", par.sigma_s);
        fp.readOrWrite("noise", par.noise);
        fp.readOrWrite("alpha", par.alpha);
        fp.readOrWrite("beta", par.beta);
        
        // Enums
        if(fp.isReading()){
            par.sc_type = str2score(fp.read("sc_type").c_str());
            par.l_type = str2learn(fp.read("l_type").c_str());
        }
        else if(fp.isWriting()){
            fp.write("sc_type",score2str(par.sc_type));
            fp.write("l_type", learn2str(par.l_type));
        }
        
        fp.readOrWrite("l_all", par.l_all);
        fp.readOrWrite("epsilon", par.epsilon);
        fp.readOrWrite("force_jump", par.force_jump);
        
        fp.readOrWrite("kernel.name", par.kernel.name);
        fp.readOrWrite("kernel.hp_mean", par.kernel.hp_mean);
        fp.readOrWrite("kernel.hp_std", par.kernel.hp_std);
        
        fp.readOrWrite("mean.name", par.mean.name);
        fp.readOrWrite("mean.coef_mean", par.mean.coef_mean);
        fp.readOrWrite("mean.coef_std", par.mean.coef_std);
        
        fp.readOrWrite("crit_name", par.crit_name);
        fp.readOrWrite("crit_params", par.crit_params);
        
    }
  } //namespace utils
} //namespace bayesopt

