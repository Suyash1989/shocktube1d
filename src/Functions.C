
#include "Input.H"
#include "Functions.H"

//Allocate dynamic array
double** AllocateDynamicArray( int nRows, int nCols){
      double** dArray;
      dArray = new double*[nRows];
      for( int i = 0 ; i < nRows; i++)
	dArray[i] = new double [nCols];
return dArray;
}

//Free dynamic array
void FreeDynamicArray(double** dArray){
      delete [] *dArray;
      delete [] dArray;
return;
}

//Initialize state vector w
void Initialize_w(input* inputData, double** Array, double val){
	double rho, u, E;
	for(int j=0;j<N;j++){

		if(j*val<inputData->getDPx()){
			rho = inputData->getRhoL();
			u = inputData->getUL();
			E = inputData->getEL();
		}else{
			rho = inputData->getRhoR();
			u = inputData->getUR();
			E = inputData->getER();
		}
       		
		Array[0][j] = rho;
        	Array[1][j] = rho*u;
        	Array[2][j] = E;
	}
return;
}

// Initialize flux vector f
void Initialize_f(input* inputData, double** Array, double val){
	double rho, u, E, p;
	for(int j=0;j<N;j++){

		if(j*val<inputData->getDPx()){
			rho = inputData->getRhoL();
			u   = inputData->getUL();
			E   = inputData->getEL();
			p = 0.4 * (E - 0.5 * rho * u * u); 
		}else{
			rho = inputData->getRhoR();
			u   = inputData->getUR();
			E   = inputData->getER();
			p = 0.4 * (E - 0.5 * rho * u * u);
		}
		
		Array[0][j] = rho*u;
		Array[1][j] = rho*u*u + p;
		Array[2][j] = u*(E + p);
	}
return;
}

//Update vector f
void update_f(double**f, double** w){
	double rho, u, E, p;
	for(int j=0;j<N;j++){
	
		rho = w[0][j];
		u = w[1][j]/rho;
		E = w[2][j];
		p = 0.4 * (E - 0.5 * rho * u * u);

		f[0][j] = rho*u;
		f[1][j] = rho*u*u + p;
		f[2][j] = u*(E + p);
	}
return;
}

//Calculate alpha
double alpha(double** w, int Row, int Col, int sign){
	double rho, p, u_0, u_1, c_0, c_1;
	rho = w[0][Col];
	u_0 = w[1][Col]/rho;
	p = 0.4*(w[2][Col] - 0.5*rho*u_0*u_0);
	c_0 = SQRT(1.4*p/rho);
	rho = w[0][Col+sign];
	u_1 = w[1][Col+sign]/rho;
	p = 0.4*(w[2][Col+sign] - 0.5*rho*u_1*u_1);
	c_1 = SQRT(1.4*p/rho);

return C*MAXI(fabs(u_0) + c_0, fabs(u_1) + c_1);
}


double flux(double** f, double** w, int Row, int Col, int sign, string FA){
	if(!strcmp(FA.c_str(),"LF"))
		return flux_LF(f, w, Row, Col, sign);
	else if(!strcmp(FA.c_str(),"SIMPLE"))
		return flux_simple(f, w, Row, Col, sign);
	else if(!strcmp(FA.c_str(),"ROE"))
		return flux_Roe(f, w, Row, Col, sign);
	else if(!strcmp(FA.c_str(),"HLL"))
		return flux_HLL(f, w, Row, Col, sign);
	else if(!strcmp(FA.c_str(),"HLLC"))
		return flux_HLLC(f, w, Row, Col, sign);
	else{
		cout<<"Unknown type of flux approximation!"<<endl;
		cout<<"Implemented flux types:"<<endl;
		cout<<"LF : Lax-Friedrichs"<<endl;
		cout<<"ROE : Roe average"<<endl;
		cout<<"HLL : HLL"<<endl;
		cout<<"HLLC : HLLC"<<endl;
		exit(1);
	}
}

//Calculate flux: Lax-Friedrichs
double flux_LF(double** f, double** w, int Row, int Col, int sign){
double h = 1.0/N;
double dt = tau(w,h);
return (0.5*(f[Row][Col+sign] + f[Row][Col]) - sign*0.5*(w[Row][Col+sign] - w[Row][Col])*h/dt);	
}

//Calculate flux: Simple
double flux_simple(double** f, double** w, int Row, int Col, int sign){
return (0.5*(f[Row][Col+sign] + f[Row][Col]) - sign*0.5*alpha(w,Row,Col,sign)*(w[Row][Col+sign] - w[Row][Col]));
}

//Calculate flux by Roe average method
double flux_Roe(double** f, double** w, int Row, int Col, int sign){
	double rho_0, E_0, u_0, p_0, h_0, rho_1, E_1, u_1, p_1, h_1;
	
	//1.Roe averaging
	rho_0 = w[0][Col];
	E_0 = w[2][Col];
	u_0 = w[1][Col]/rho_0;
	p_0 = 0.4*(E_0 - 0.5*rho_0*u_0*u_0);
	h_0 = (E_0 + p_0)/rho_0;

	rho_1 = w[0][Col+sign];
	E_1 = w[2][Col+sign];
	u_1 = w[1][Col+sign]/rho_1;
	p_1 = 0.4*(E_1 - 0.5*rho_1*u_1*u_1);
	h_1 = (E_1 + p_1)/rho_1;
	
	//u_{i+/-1/2}
	double u = (SQRT(rho_0)*u_0 + SQRT(rho_1)*u_1)/(SQRT(rho_0) + SQRT(rho_1));

	//h_{i+/-1/2}
	double h = (SQRT(rho_0)*h_0 + SQRT(rho_1)*h_1)/(SQRT(rho_0) + SQRT(rho_1));

	//c_{i+/-1/2}
	double c = SQRT(0.4*(h - 0.5*u*u));
	
	//2.Calculate R
	double** R = AllocateDynamicArray(3,3);
	R[0][0] = 1.0;		R[0][1] = 1.0;		R[0][2] = 1.0;
	R[1][0] = (u-c);	R[1][1] = u;		R[1][2] = (u+c);
	R[2][0] = h - u*c;	R[2][1] = 0.5*u*u;	R[2][2] = h + u*c;

	//3.Calculate diagonal matrix
	double** Lambda = AllocateDynamicArray(3,3);
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			Lambda[i][j] = 0.0;
	double Al=1, A = Al*c;
	if(fabs(u-c)<A)	Lambda[0][0] = fabs(0.5*(A + (u-c)*(u-c)/A));		else Lambda[0][0] = fabs(u-c);
	if(fabs(u)<A)	Lambda[1][1] = fabs(0.5*(A + (u)*(u)/A));		else Lambda[1][1] = fabs(u);
	if(fabs(u+c)<A)	Lambda[2][2] = fabs(0.5*(A + (u+c)*(u+c)/A));		else Lambda[2][2] = fabs(u+c);

	//4.Calculate R inverse
	double** Rinv = AllocateDynamicArray(3,3);
	double G = 0.4/c/c;
	Rinv[0][0] = 0.5*(1.0+G*(u*u-h)+u/c);		Rinv[0][1] = -0.5*(G*u+(1.0/c));		Rinv[0][2] = 0.5*G;
	Rinv[1][0] = -G*(u*u-h);			Rinv[1][1] = G*u;				Rinv[1][2] = -G;
	Rinv[2][0] = 0.5*(1.0+G*(u*u-h)-u/c);		Rinv[2][1] = -0.5*(G*u-(1.0/c));		Rinv[2][2] = 0.5*G;
	
	//5. Matrix-Vector multiplication
	double W[3];
	for(int i=0;i<3;i++)
		W[i] = sign*(w[i][Col+sign] - w[i][Col]);
	
	MVMult(Rinv,W);
	MVMult(Lambda,W);
	MVMult(R,W);

	FreeDynamicArray(Rinv);
	FreeDynamicArray(Lambda);
	FreeDynamicArray(R);

return (0.5*(f[Row][Col+sign] + f[Row][Col]) - 0.5*C_Roe*W[Row]);
}

//Calculate flux: Harten, Lax, van Leer (HLL)
double flux_HLL(double** f, double** w, int Row, int Col, int sign){
	double rho_L = w[0][Col];
	double u_L = w[1][Col]/rho_L;
	double p_L = 0.4*(w[2][Col] - 0.5*rho_L*u_L*u_L);
	double a_L = SQRT(1.4*p_L/rho_L);

	double rho_R = w[0][Col+sign];
	double u_R = w[1][Col+sign]/rho_R;
	double p_R = 0.4*(w[2][Col+sign] - 0.5*rho_R*u_R*u_R);
	double a_R = SQRT(1.4*p_R/rho_R);


	double S_L, S_R, F_L, F_R;
	if(sign>0) S_L = u_L - a_L; else S_L = u_R - a_R;
	if(sign>0) S_R = u_R + a_R; else S_R = u_L + a_L;

//	S_L = MINI(u_L-a_L,u_R-a_R);
//	S_R = MAXI(u_L+a_L,u_R+a_R);

//	double u_star = 0.5*(u_L+u_R) + (a_L-a_R)/0.4;
//	double a_star = 0.5*(a_L+a_R) + 0.25*0.4*(u_L-u_R);
//	S_L = MINI(u_L-a_L,u_star-a_star);
//	S_R = MAXI(u_R+a_R,u_star+a_star);

/* Lax-Friedrichs
double h = 1.0/N;
double dt = tau(w,N,h);
S_R = h/dt;
S_L = -S_R;
*/
	if(sign>0){
		if(S_L >= 0.0)			return f[Row][Col];
		if( (S_L<=0.0) && (S_R>=0.0) )  return ( S_R*f[Row][Col] - S_L*f[Row][Col+sign] + S_L*S_R*( w[Row][Col+sign] - w[Row][Col]) )/(S_R - S_L);
		if(S_R <= 0.0)			return f[Row][Col+sign];
	}else{
		if(S_L >= 0.0)			return f[Row][Col+sign];
		if( (S_L<=0.0) && (S_R>=0.0) )  return ( S_R*f[Row][Col+sign] - S_L*f[Row][Col] + S_L*S_R*( w[Row][Col] - w[Row][Col+sign]) )/(S_R - S_L);
		if(S_R <= 0.0)			return f[Row][Col];
	}

}

//Calculate flux: HLLC
double flux_HLLC(double** f, double** w, int Row, int Col, int sign){

	double rho_L, u_L, p_L, rho_R, u_R, p_R, W_L, W_R, F_L, F_R;	
	if(sign>0){
		rho_L = w[0][Col];
		u_L = w[1][Col]/rho_L;
		p_L = 0.4*(w[2][Col] - 0.5*rho_L*u_L*u_L);

		rho_R = w[0][Col+sign];
		u_R = w[1][Col+sign]/rho_R;
		p_R = 0.4*(w[2][Col+sign] - 0.5*rho_R*u_R*u_R);

		W_L = w[Row][Col];
		W_R = w[Row][Col+sign];

		F_L = f[Row][Col];
		F_R = f[Row][Col+sign];
	}else{
		rho_L = w[0][Col+sign];
		u_L = w[1][Col+sign]/rho_L;
		p_L = 0.4*(w[2][Col+sign] - 0.5*rho_L*u_L*u_L);

		rho_R = w[0][Col];
		u_R = w[1][Col]/rho_R;
		p_R = 0.4*(w[2][Col] - 0.5*rho_R*u_R*u_R);

		W_L = w[Row][Col+sign];
		W_R = w[Row][Col];

		F_L = f[Row][Col+sign];
		F_R = f[Row][Col];
	}

	double a_L = SQRT(1.4*p_L/rho_L);
	double a_R = SQRT(1.4*p_R/rho_R);
	double S_L = u_L - a_L;
	double S_R = u_R + a_R;

//	double S_L = MINI(u_L-a_L,u_R-a_R);
//	double S_R = MAXI(u_L+a_L,u_R+a_R);

/*	double u_star = 0.5*(u_L+u_R) + (a_L-a_R)/0.4;
	double a_star = 0.5*(a_L+a_R) + 0.25*0.4*(u_L-u_R);
	double S_L = MINI(u_L-a_L,u_star-a_star);
	double S_R = MAXI(u_R+a_R,u_star+a_star);
*/
/* Lax-Friedrichs
double h = 1.0/N;
double dt = tau(w,N,h);
S_R = h/dt;
S_L = -S_R;
*/

	double S_star = ( p_R - p_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R) )/( rho_L*(S_L-u_L) - rho_R*(S_R-u_R) );

	double W_star_L, W_star_R;
	if(Row==0){
		W_star_L = rho_L*(S_L-u_L)/(S_L-S_star);
		W_star_R = rho_R*(S_R-u_R)/(S_R-S_star);
	}else if(Row==1){
		W_star_L = rho_L*(S_L-u_L)/(S_L-S_star)*S_star;
		W_star_R = rho_R*(S_R-u_R)/(S_R-S_star)*S_star;
	}else{
		W_star_L = rho_L*(S_L-u_L)/(S_L-S_star)*(W_L/rho_L + (S_star-u_L)*(S_star+p_L/rho_L/(S_L-u_L)));
		W_star_R = rho_R*(S_R-u_R)/(S_R-S_star)*(W_R/rho_R + (S_star-u_R)*(S_star+p_R/rho_R/(S_R-u_R)));
	}

        double F_star_L = F_L + S_L*(W_star_L-W_L);
        double F_star_R = F_R + S_R*(W_star_R-W_R);

	if(S_L >= 0.0)					return F_L;
	else if( (S_L <= 0.0) && (S_star >= 0.0) ) 	return F_star_L;
	else if( (S_star <= 0.0) && (S_R >= 0.0) ) 	return F_star_R;
	else if(S_R <= 0.0) 				return F_R;
}

//Matrix-Vector Multiplication
void MVMult(double** M, double* V){
	double sum[3]={0.0};
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			sum[i] = sum[i]  + M[i][j]*V[j];

	for(int i=0;i<3;i++)
		V[i] = sum[i];
}


//Calculate time step
double tau(double** w, double val){
	double lambda = 0.0, lambda_max = 0.0;
	double rho, p, u, c;
	for(int j=0;j<N;j++){
		rho = w[0][j];
		u = w[1][j]/rho;
		p = 0.4*(w[2][j] - 0.5*rho*u*u);
		c = SQRT(1.4*p/rho);
		lambda = fabs(u) + c;
		
		if(lambda>lambda_max)
			lambda_max = lambda;	
	}
return CFL*val/lambda_max;
}

//Solve for vector w
void solve(input* inputData, double** w, double** f, double h){
	double** wnew = AllocateDynamicArray(3,N);
	double dt=tau(w,h);

	string FA = inputData->getFA();

	for(int i=0;i<3;i++){
		for(int j=0;j<N;j++){
			if(j==0){
				//Reflection Boundary Conditions
				//w[0][j] = w[0][j+1];
				//w[1][j] = -w[1][j+1];
				//w[2][j] = w[2][j+1];							
				//Boundary conditions at x=0
				wnew[i][j] = w[i][j] - (1.0/h) * dt * (flux(f,w,i,j,+1,FA) - f[i][j]);
			}
			else if(j==N-1){
				//Reflection Boundary Conditions
				//w[0][j] = w[0][j-1];
				//w[1][j] = -w[1][j-1];
				//w[2][j] = w[2][j-1];						
				//Boundary conditions at x=1
				wnew[i][j] = w[i][j] - (1.0/h) * dt * (f[i][j] - flux(f,w,i,j,-1,FA));
			}
			else
				//Interior points
				wnew[i][j] = w[i][j] - (1.0/h) * dt * (flux(f,w,i,j,+1,FA) - flux(f,w,i,j,-1,FA));
			
		}
	}

	for(int i=0;i<3;i++)
		for(int j=0;j<N;j++)
			w[i][j] = wnew[i][j];

	FreeDynamicArray(wnew);
}

// Print flow variables
void PrintFlowVariables(input* inputData, double** rho, double** u, double** E, double** e, double** p, double** Mach, int count, double h){

	string FA		= inputData->getFA();
	string Density		= "Density_";
	string TEnergy		= "TEnergy_";
	string IEnergy		= "IEnergy_";
	string Pressure		= "Pressure_";
	string Velocity		= "Velocity_";
	string Mach_Number	= "Mach_Number_";

	//Print flow variables in respective files
	PrintArray(rho,count,N,h,Density.append(FA).append(".dat"));
	PrintArray(E,count,N,h,TEnergy.append(FA).append(".dat"));
	PrintArray(e,count,N,h,IEnergy.append(FA).append(".dat"));
	PrintArray(p,count,N,h,Pressure.append(FA).append(".dat"));
	PrintArray(u,count,N,h,Velocity.append(FA).append(".dat"));
	PrintArray(Mach,count,N,h,Mach_Number.append(FA).append(".dat"));
return;
}

//Print array
void PrintArray(double** Array, int nRows, int nCols, double h, string fileName){
	ofstream file(fileName.c_str());
	for(int j=0;j<nCols;j++){
			file<<(j*h+h/2)<<"\t";
		for(int i=0;i<nRows;i++){
			//cout<<Array[i][j]<<"\t";
			file<<Array[i][j]<<"\t";
		}
		//cout<<"\n";
		file<<"\n";
	}
}

// Antidiffusion for density
/*
void antidiffusionStep(double** w, double h){
 
    double* rho_old = new double[N]();
    double* rho = new double[N]();
    double* u = new double[N]();
    double* p = new double[N]();
 
    for(int j=0;j<N;j++){
        rho[j] = w[0][j];
        rho_old[j] = w[0][j];
        u[j] = w[1][j]/rho[j];
        p[j] = 0.4*(w[2][j] - 0.5*rho[j]*u[j]*u[j]);
    }
 
    // Calculation of gradients
    double* grad_rho = new double[N]();
    for(int j=0;j<N-1;j++)
        grad_rho[j] = (rho[j+1] - rho[j])/h;
 
    //Calculation of antidiffusion coefficient
    double a_L, a_R, S_L, S_R, eps, eps_max=0.0;
    for(int j=0;j<N-1;j++){
        a_L = SQRT(1.4*p[j]/rho[j]);
        a_R = SQRT(1.4*p[j+1]/rho[j+1]);
        S_L = MINI(fabs(u[j])-a_L, fabs(u[j+1])-a_R);
        S_R = MAXI(fabs(u[j])+a_L, fabs(u[j+1])+a_R);
        eps = fabs( S_L*S_R*(rho[j+1] - rho[j])/(S_R-S_L) );
         
        if(eps>eps_max) eps_max = eps;
    }
     
    //eps_max = eps_max*1e-5;
 
    // Calculation of fluxes
    double* F = new double[N]();
    for(int j=0;j<N-1;j++){
         
        if( fabs(grad_rho[j]) <= fabs(grad_rho[j+1]) )
            F[j] = -eps_max*grad_rho[j];
        else
            F[j] = -eps_max*grad_rho[j+1];
    }
 
    // Calculate the total antidiffusion
    double* AD = new double[N]();
    for(int j=0;j<N-1;j++){
        AD[j] = -(F[j] + F[j+1])/h;
    }
 
    // Calculate corrected rho
    double tau_p;
    for(int j=1;j<N;j++){
         
        // Calculate pseudo time step      
        tau_p = 0.25*h*h/eps_max*1e-2;
 
        rho[j] = rho_old[j] + tau_p*AD[j-1];
    }
 
    for(int j=0;j<N-1;j++)
            cout<<eps_max<<"\t"<<grad_rho[j]<<"\t"<<F[j]<<"\t"<<AD[j]<<"\t"<<rho[j]<<endl;
 
    std::memcpy(w[0], rho, N*sizeof(double));
 
    // Update other conservative variables
    for(int j=0;j<N;j++){
        w[1][j] = rho[j]*u[j];
        w[2][j] = 2.5*p[j] + 0.5*rho[j]*u[j]*u[j];
    }
     
    delete [] rho;
    delete [] rho_old;
    delete [] u;
    delete [] p;
    delete [] grad_rho;
    delete [] F;
    delete [] AD;
return;
}*/
