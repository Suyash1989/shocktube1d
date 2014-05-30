
#include "Input.H"
#include "Functions.H"

//Global variables
int N;
double CFL;

int main(int argc, char* argv[]){

	// Input for initial condition 
	input* 	inputData = new input;

	// Read input file
	inputData->readInputFile();
	
	N = inputData->getN();				// Number of CVs
	CFL = inputData->getCFL();			// CFL Number
	double h = 1.0/N;			 	// delta x	
	double curr_time = 0;				// Keep track of current time
	int count=0;					// Count for no. of time steps
	int M = 3*N/CFL;				// Estimated number of time-levels
	
	// Allocate the arrays 
	double** w = AllocateDynamicArray(3,N);		// State vector
	double** f = AllocateDynamicArray(3,N);		// Flux vector
	double** rho = AllocateDynamicArray(M,N);	// Density
	double** E = AllocateDynamicArray(M,N);		// Total energy 
	double** e = AllocateDynamicArray(M,N);		// Specific internal energy
	double** p = AllocateDynamicArray(M,N);		// Pressure
	double** Mach = AllocateDynamicArray(M,N);	// Mach Number
	double** u = AllocateDynamicArray(M,N);		// Velocity

	//Initial condition
	Initialize_w(inputData,w,h);
	Initialize_f(inputData,f,h);

	cout << "\nInitial time step: "<<tau(w,h)<<endl;

	//Calculation of w at next time steps
	while(curr_time<=inputData->getTime()){

		//Current time
		curr_time = curr_time + tau(w,h);
		
		//Solve for vector w
		solve(inputData, w, f, h);
	
		//Update vector f with new values of vector w
		update_f(f,w);
		
		//Extract flow variables from vector w
		for(int j=0;j<N;j++){
			rho[count][j] = w[0][j];
			u[count][j] = w[1][j]/rho[count][j];
			E[count][j] = w[2][j];
			p[count][j] = 0.4*(E[count][j] - 0.5*rho[count][j]*u[count][j]*u[count][j]);
			e[count][j] = w[2][j]/rho[count][j] - 0.5*u[count][j]*u[count][j];
			Mach[count][j] = u[count][j]/SQRT(1.4*p[count][j]/rho[count][j]);
		}

		//Count number of time steps
		count++;	

		cout<<"\nTime step No. : "<<count<<"\tTime = "<<curr_time<<endl;
	}
	
	PrintFlowVariables(inputData, rho, u, E, e, p, Mach, count, h);

	cout<<"\nNumber of time steps = "<<count<<endl;
	cout<<"\nCurrent time = "<<curr_time<<endl;

	FreeDynamicArray(w);
	FreeDynamicArray(f);
	FreeDynamicArray(rho);
	FreeDynamicArray(E);
	FreeDynamicArray(p);
	FreeDynamicArray(u);
	FreeDynamicArray(Mach);
	delete inputData;

return 0;
}

