#ifndef INPUT_H_
#define INPUT_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <limits>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

extern int N;
extern double CFL;

/*! 
 * \brief This class is used for reading and storing the data from input file 
 */
class input
{
    private:
        /// PRIVATE VARIABLES
        string  title;      // title of the document
	string	wdir;	    // working directory
	int	N;	    // Number of control volumes
	double  t;	    // Time for simulation
	double  CFLn;	    // CFL Number
	string  FA;	    // Flux Approximation
	double  DP_x;	    // Diaphragm position
	double  rho_L;	    // Density on the left of diaphragm
	double  u_L;	    // Velocity on the left of diaphragm
	double  E_L;	    // Total energy on the left of diaphragm
	double  rho_R;	    // Density on the right of diaphragm
	double  u_R;	    // Velocity on the right of diaphragm
	double  E_R;	    // Total energy on the right of diaphragm
        
    protected:

    public:
        /// DEFAULT CONSTRUCTOR ///
        input();

        /// DESTRUCTOR
        ~input(){};

        /// GETTERS ///  
        string          getTitle()      {return title;};
        string          getWdir()       {return wdir;};
	int		getN()		{return N;};
	double		getTime()	{return t;};
	double		getCFL()	{return CFLn;};
	string		getFA()		{return FA;};
	double		getDPx()	{return DP_x;};
        double          getRhoL()       {return rho_L;};
        double          getUL()         {return u_L;};
        double          getEL()         {return E_L;};
        double          getRhoR()       {return rho_R;};
        double          getUR()         {return u_R;};
        double          getER()         {return E_R;};

        /// PUBLIC INTERFACE METHOD
        void readInputFile();
};

#endif
