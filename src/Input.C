
#include "Input.H"

/*!
 *	Default constructor for input object
*/
input::input(){

    // Default values for the input parameters
    title = "Test_1";
    N	  = 100;
    t     = 0.2;
    CFLn  = 1.0;
    FA    = "LF";
    DP_x  = 0.5;
    rho_L = 1.0;
    u_L   = 0.0;
    E_L   = 2.5;
    rho_R = 0.125;
    u_R   = 0.0;
    E_R   = 0.25;
}

/*!
 *	Reads Input file
*/
void input::readInputFile(){

    string lineString;
    string dummyString;
    string separator;
    
    ifstream inputFile;
    inputFile.open("input",ios::in);

    if (inputFile.is_open()==false){
        cout << "Unable to open input file! Aborting... " << endl;
        exit(0);
    }

    cout.precision(16);
    cout << scientific;
    
    while (!inputFile.eof())
    {
        // Get a line and store in lineString
        getline(inputFile, lineString, '\n');

        // If the first character of the line is not a '#'
        if (lineString.c_str()[0] != '#')
        {
            istringstream iss(lineString);
            iss >> dummyString;
	    iss >> separator;
	    if(separator != ":"){
                cout << endl << "Please provide ':' as a separator between "<<dummyString<<" and "<<separator<<" !";
                cout << endl << "Aborting...";
		cout << endl;
                exit(0);    
            }else if(dummyString == "title")
                iss >> title;
            else if(dummyString == "wdir")
                iss >> wdir;
            else if(dummyString == "N")
                iss >> N;
            else if(dummyString == "End_time")
                iss >> t;
            else if(dummyString == "FA")
                iss >> FA;
            else if(dummyString == "CFL")
                iss >> CFL;
            else if(dummyString == "DP_x")
                iss >> DP_x;
            else if(dummyString == "rho_L")
                iss >> rho_L;
            else if(dummyString == "u_L")
                iss >> u_L;
            else if(dummyString == "E_L")
                iss >> E_L;
            else if(dummyString == "rho_R")
                iss >> rho_R;
            else if(dummyString == "u_R")
                iss >> u_R;
            else if(dummyString == "E_R")
                iss >> E_R;
            else{
                cout << endl << "Unknown keyword in the settings file : " << dummyString;
                cout << endl << "Aborting...";
                exit(0);    
            }
        }
        
    }

    // Report the settings read from the file.
    cout << endl;
    cout << "*-------------------------Input Parameters-------------------------*"<<endl;
    cout << endl;
    cout << "Title of the simulation			: " << title << endl;
    cout << "Working Directory			: " << wdir << endl;
    cout << "Number of control volumes		: " << N << endl;
    cout << "End time of simulation			: " << t << endl;
    cout << "Flux Approximation			: " << FA << endl;
    cout << "CFL Number				: " << CFL << endl;
    cout << "Diaphragm position			: " << DP_x << endl;
    cout << "Density on the left of diaphragm	: " << rho_L << endl;
    cout << "Velocity on the left of diaphragm	: " << u_L << endl;
    cout << "Total Energy on the left of diaphragm	: " << E_L << endl;
    cout << "Density on the right of diaphragm	: " << rho_R << endl;
    cout << "Velocity on the right of diaphragm	: " << u_R << endl;
    cout << "Total Energy on the right of diaphragm	: " << E_R << endl;
    cout << endl;
    cout << "*------------------------------------------------------------------*"<<endl;

    inputFile.close();
    
    return;
}


