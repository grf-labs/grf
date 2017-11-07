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

#include "testfunctions.hpp"
#include "displaygp.hpp"



// Unfortunately OpenGL functions require no parameters, so the object has to be global.
bayesopt::utils::DisplayProblem1D GLOBAL_MATPLOT;

void display( void ){ GLOBAL_MATPLOT.display(); }
void reshape( int w,int h ){ GLOBAL_MATPLOT.reshape(w,h); }
void idle( void ) { glutPostRedisplay(); } 

void mouse(int button, int state, int x, int y ){ GLOBAL_MATPLOT.mouse(button,state,x,y); }
void motion(int x, int y ){ GLOBAL_MATPLOT.motion(x,y); }
void passive(int x, int y ){ GLOBAL_MATPLOT.passivemotion(x,y); }

void keyboard(unsigned char key, int x, int y)
{
    GLOBAL_MATPLOT.keyboard(key, x, y); 
    if(key=='r')   //Toogle run/stop
      { 
	GLOBAL_MATPLOT.toogleRUN();
      }
    if(key=='s')   //Activate one step
      { 
	GLOBAL_MATPLOT.setSTEP();
      }
}

int menu()
{
  std::string input;
  int option = 0;
  while ((option < 1) || (option > 5))
    {
      std::cout << "Please select an option for the parameters:\n\n"
		<< "  1- Default parameters.\n"
		<< "  2- Student t process model.\n"
		<< "  3- Combined kernel.\n"
		<< "  4- Lower Confidence Bound.\n"
		<< "  5- A-optimality criteria.\n\n"
		<< "Select [1-5]>";

      std::cin >> input;
      std::istringstream is(input);
      is >> option;
    }
  return option;
}

int main(int nargs, char *args[])
{
  bopt_params parameters = initialize_parameters_to_default();
  parameters.n_init_samples = 7;
  parameters.n_iterations = 100;
  parameters.verbose_level = 2;

  switch( menu() )
    {
    case 1: break;
    case 2: 
      {
	set_surrogate(&parameters,"sStudentTProcessNIG"); 
	parameters.n_iter_relearn = 5;
	break;
      }
    case 3:   
      { 
	set_kernel(&parameters,"kSum(kPoly3,kRQISO)");
	double mean[128] = {1, 1, 1, 1};
	double std[128] = {5, 5, 5, 5};
	size_t nhp = 4;
	memcpy(parameters.kernel.hp_mean, mean, nhp * sizeof(double));
	memcpy(parameters.kernel.hp_std,std, nhp * sizeof(double));
	parameters.kernel.n_hp = nhp;
	break;
      }
    case 4:
      set_criteria(&parameters,"cLCB");
      parameters.crit_params[0] = 5;
      parameters.n_crit_params = 1;
      break;      
    case 5:
      set_criteria(&parameters,"cAopt");
      parameters.n_crit_params = 0;
      break;
    default:
      break;
    };

  bayesopt::Parameters parameters_class(parameters);
  boost::scoped_ptr<ExampleOneD> opt(new ExampleOneD(parameters_class));
  GLOBAL_MATPLOT.init(opt.get(),1);

  glutInit(&nargs, args);
  glutCreateWindow(50,50,800,650);
  glutDisplayFunc( display );
  glutReshapeFunc( reshape );
  glutIdleFunc( idle );
  glutMotionFunc( motion );
  glutMouseFunc( mouse );
  glutPassiveMotionFunc(passive);    
  glutKeyboardFunc( keyboard );        
  glutMainLoop();    

  return 0;
}
