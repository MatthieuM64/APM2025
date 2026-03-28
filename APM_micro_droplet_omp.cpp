/*C++ CODE - MANGEAT MATTHIEU - 2025*/
/*ACTIVE POTTS MODEL - MESTABILITITY OF THE LIQUID STATE - DROPLET INSERTION*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <omp.h>

using namespace std;

//Personal libraries.
#include "lib/random_OMP.cpp"
#include "lib/special_functions.cpp"

void modulo(int &x, const int &LX);

//////////////////////////
///// PARTICLE CLASS /////
//////////////////////////

//Position and spin defining the Potts particle
class particle
{
	public:
	
	int x,y;
	int spin;
	
	particle(const int &LX, const int &LY, const int &sigma);
	particle(const int &LX, const int &LY, const int &sigma, const int &Rd);
	void move(const double &r, const double &epsilon, const int &LX, const int &LY);
	void flip(const int &sigma);
};

//Creation of the spin
particle::particle(const int &LX, const int &LY, const int &sigma)
{
	spin=sigma;
	x=int(ran()*LX);
	y=int(ran()*LY);
}

//Creation of the spin inside a disk of radius Rd
particle::particle(const int &LX, const int &LY, const int &sigma, const int &Rd)
{
	spin=sigma;
	const double r=Rd*sqrt(ran()), phi=2*M_PI*ran();
	x=LX/2+int(r*cos(phi));
	y=LY/2+int(r*sin(phi));
	
	/*do
	{
		x=int((2*ran()-1)*Rd);
		y=int((2*ran()-1)*Rd);
	}while (x*x+y*y<Rd*Rd);
	x+=LX/2;
	y+=LY/2;*/
}

//Update the particle position.
void particle::move(const double &r, const double &epsilon, const int &LX, const int &LY)
{
	int p=spin; //Probability to hop in the favored direction (r<epsilon/3).
	if (r>epsilon/3.) //Probability to hop to a random direction (r>epsilon/3).
	{
		p=int(4*ran()); //Random initial direction (0,1,2,3)
	}
	
	if (p==0) //Move right.
	{
		x++;
	}
	else if (p==1) //Move up.
	{
		y++;
	}
	else if (p==2) //Move left.
	{
		x--;
	}
	else //Move down.
	{
		y--;
	}
	modulo(x,LX);
	modulo(y,LY);
}

//Flip the spin.
void particle::flip(const int &sigma)
{
	spin=sigma;
}

///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Modulo function.
void modulo(int &x, const int &LX)
{
	if (x<0)
	{
		x+=LX;
	}
	else if (x>=LX)
	{
		x-=LX;
	}
}

//Maximum state
int max_state(const vector<int> &rhos)
{
	int rhomax=rhos[0];
	for (int j=1;j<rhos.size();j++)
	{
		if (rhos[j]>rhomax)
		{
			rhomax=rhos[j];
		}
	}
	
	if (rhomax==0)
	{
		cerr << "THE MAX DENSITY IS ZERO -> ERROR" << endl;
		abort();
	}
	
	vector<int> kmax;
	for (int j=0;j<rhos.size();j++)
	{
		if (rhos[j]==rhomax)
		{
			kmax.push_back(j);
		}
	}
	
	return kmax[int(kmax.size()*ran())];
}

//Total average on all space.
double average(const vector< vector<int> > &RHO, const int &LX, const int &LY)
{
	int rhoAv=0;
	for (int x0=0; x0<LX; x0++)
	{
		for (int y0=0; y0<LY; y0++)
		{
			rhoAv+=RHO[x0][y0];
		}
	}
	return double(rhoAv)/(LX*LY);
}

//Export density.
void exportDensity(const vector< vector<int> > &RHO, const double &beta, const double &D, const double &epsilon, const double &rho0, const int &Rd, const double &rhod, const int &LX, const int &LY, const int &init, const int &RAN, const int &t)
{
	//Creation of the file.
	int returnSystem=system("mkdir -p data_APM_dynamics2d/");
	stringstream ssDensity;
	ssDensity << "./data_APM_dynamics2d/APM_density_beta=" << beta << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_Rd=" << Rd << "_rhod=" << rhod << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameDensity = ssDensity.str();
	
	ofstream fileDensity(nameDensity.c_str(),ios::trunc);
	fileDensity.precision(6);
		
	//Write in the file.
	for (int y0=0;y0<LY;y0++)
	{
		for (int x0=0; x0<LX; x0++)
		{
			fileDensity << RHO[x0][y0] << "\t";
		}
		fileDensity << endl;
	}
	fileDensity.close();
}

//Export state.
void exportState(const vector<vector< vector<int> > > &RHOS, const vector< vector<int> > &RHO, const double &beta, const double &D, const double &epsilon, const double &rho0, const int &Rd, const double &rhod, const int &LX, const int &LY, const int &init, const int &RAN, const int &t)
{
	//Creation of the file.
	int returnSystem=system("mkdir -p data_APM_dynamics2d/");
	stringstream ssState;
	ssState << "./data_APM_dynamics2d/APM_state_beta=" << beta << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_Rd=" << Rd << "_rhod=" << rhod << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameState = ssState.str();
	
	ofstream fileState(nameState.c_str(),ios::trunc);
	fileState.precision(6);
		
	//Write in the file.
	for (int y0=0;y0<LY;y0++)
	{
		for (int x0=0; x0<LX; x0++)
		{
			if (RHO[x0][y0]==0)
			{
				fileState << "nan" << "\t";
			}
			else
			{
				vector<int> rhos={RHOS[0][x0][y0],RHOS[1][x0][y0],RHOS[2][x0][y0],RHOS[3][x0][y0]};
				fileState << max_state(rhos) << "\t";
			}
		}
		fileState << endl;
	}
	fileState.close();
}

//Introduce a droplet in the middle of the domain.
void introduceDroplet(vector<particle> &POTTS, vector<vector<vector<int> > > &RHOS, vector<vector<int> > &RHO, int &Npart, const int &Rd, const int &sigma, const int &DelN, const int &LX, const int &LY)
{
	for (int i=0; i<Npart; i++) //change the state of the particles inside the droplet.
	{
		const int x0=POTTS[i].x, y0=POTTS[i].y;
		if (square(x0-LX/2)+square(y0-LY/2)<square(Rd))
		{
			RHOS[POTTS[i].spin][x0][y0]--;
			POTTS[i].spin=sigma;
			RHOS[sigma][x0][y0]++;
		}
	}
	
	for (int i=0; i<DelN; i++) //Add new particles inside the droplet.
	{
		particle Potts(LX,LY,sigma,Rd);
		POTTS.push_back(Potts);
		RHOS[Potts.spin][Potts.x][Potts.y]++;
		RHO[Potts.x][Potts.y]++;
	}
	
	//Update the particle number.
	Npart+=DelN;
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &beta, double &D, double &epsilon, double &rho0, int &Rd, double &rhod, int &LX, int &LY, int &tmax, int &init, int &RAN, int &THREAD_NUM)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-beta=" ))
		{
			beta=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-D=" ))
		{
			D=atof(argv[i]+3);
		}
		else if (strstr(argv[i], "-epsilon=" ))
		{
			epsilon=atof(argv[i]+9);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-Rd=" ))
		{
			Rd=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-rhod=" ))
		{
			rhod=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-ran=" ))
		{
			RAN=atoi(argv[i]+5);
		}
		else if (strstr(argv[i], "-threads=" ))
		{
			THREAD_NUM=atoi(argv[i]+9);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
	cout << "-beta=" << beta << "_D=" << D << " -epsilon=" << epsilon << " -rho0=" << rho0 << " Rd=" << Rd << " rhod=" << rhod << " -LX=" << LX << " -LY=" << LY << " -tmax=" << tmax << " -init=" << init << " -ran=" << RAN << " -threads=" << THREAD_NUM << endl;
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: beta=inverse temperature, D=diffusion rate, epsilon=self-propulsion parameter, rho0=average density, LX*LY=size of the box.
	double beta=1., D=1., epsilon=2.7, rho0=10, rhod=12;
	int LX=200, LY=200, Rd=10;
	
	//Numerical parameters: init=initial condition, tmax=maximal time, RAN=index of RNG, THREAD_NUM=number of threads.
	int init=0, tmax=600000, RAN=0, THREAD_NUM=4;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,beta,D,epsilon,rho0,Rd,rhod,LX,LY,tmax,init,RAN,THREAD_NUM);

	//OpenMP.
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	cout << OMP_MAX_THREADS << " maximum threads on this node. " << THREAD_NUM << " threads will be used." << endl;

	//Start the random number generator.
	init_gsl_ran();
	for (int k=0; k<THREAD_NUM; k++)
	{
		gsl_rng_set(GSL_r[k],THREAD_NUM*RAN+k);
	}
	
	//Number of particles.
	int Npart=int(rho0*LX*LY);
	
	//Number of particles to add inside the droplet.
	int DelN=int((rhod-rho0)*M_PI*Rd*Rd);
	
	//Density and magnetization on each sites.
	vector<vector<vector<int> > > RHOS(4,vector<vector<int> >(LX,vector<int>(LY,0)));
	vector<vector<int> > RHO(LX,vector<int>(LY,0));
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_APM_averages/");
	stringstream ssAverages;
	ssAverages << "./data_APM_averages/APM_averages_beta=" << beta << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_Rd=" << Rd << "_rhod=" << rhod << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	//Creation of particles.
	vector<particle> POTTS;	
	for (int i=0;i<Npart;i++)
	{
		particle Potts(LX,LY,2);
		POTTS.push_back(Potts);
		RHOS[Potts.spin][Potts.x][Potts.y]++;
		RHO[Potts.x][Potts.y]++;
	}
	
	//Time increment.
	const double dt=1./(4*D+exp(4*beta));
	
	//Probability to hop.
	const double proba_hop=4*D*dt;
	
	//Number of particles per threads (most accurate version when Npart is not a multiple of THREAD_NUM).
	int nth=int(Npart/THREAD_NUM);
	vector<int> Nth(THREAD_NUM,nth);
	Nth[THREAD_NUM-1]=Npart-(THREAD_NUM-1)*nth;
	
	//Number of MCS between two snaps.
	int DT=100;
	
	//Lock for the update of magnetization
	vector< vector<omp_lock_t> > lock_site(LX,vector<omp_lock_t>(LY));
	
	for (int x=0; x<LX; x++)
	{
		for (int y=0; y<LY; y++)
		{
			omp_init_lock(&lock_site[x][y]);
		}
	}
	
	//Time evolution.
	for(int t=0;t<=tmax;t++)
	{
		//Introduce the droplet after equilibration of the liquid.
		if (t==6000)
		{
			introduceDroplet(POTTS,RHOS,RHO,Npart,Rd,init,DelN,LX,LY);
			nth=int(Npart/THREAD_NUM);
			for (int threads=0;threads<THREAD_NUM-1;threads++)
			{
				Nth[threads]=nth;
			}
			Nth[THREAD_NUM-1]=Npart-(THREAD_NUM-1)*nth;
		}
		
		//Increase the number of MCS between two snaps.
		if (t==100000)
		{
			DT=1000;
		}
		
		//Export data.
		if (t%150==0 or t==tmax)
		{
			exportDensity(RHO,beta,D,epsilon,rho0,Rd,rhod,LX,LY,init,RAN,t);
			exportState(RHOS,RHO,beta,D,epsilon,rho0,Rd,rhod,LX,LY,init,RAN,t);
			
			const double ntot=average(RHO,LX,LY), p1=average(RHOS[0],LX,LY)/ntot, p2=average(RHOS[1],LX,LY)/ntot, p3=average(RHOS[2],LX,LY)/ntot, p4=average(RHOS[3],LX,LY)/ntot;
			cout << "time=" << t << " -rho=" << ntot << " -p1=" << p1 << " -p2=" << p2 << " -p3=" << p3/ntot << " -p4=" << p4 << running_time.TimeRun(" ") << endl;
			fileAverages <<  t << " " << ntot << " " << p1 << " " << p2 << " " << p3 << " " << p4 << endl;
		}
		
		//At each time-step move all particles randomly.
		#pragma omp parallel
		{
			const int actual_thread=omp_get_thread_num();
			for (int i=0; i<Nth[actual_thread]; i++)
			{
				//Choose a particle randomly (j), at the site (x0,y0).
				const int j=actual_thread*nth+int(Nth[actual_thread]*ran());
				const int x0=POTTS[j].x, y0=POTTS[j].y, spin0=POTTS[j].spin;
				
				//New possible spin-state
				int spin1=int(3*ran());
				if (spin1>=spin0)
				{
					spin1++;
				}
				
				//Probability to flip.
				omp_set_lock(&lock_site[x0][y0]);
				const int rhos=RHOS[spin0][x0][y0], rhosp=RHOS[spin1][x0][y0], rho=RHO[x0][y0];
				omp_unset_lock(&lock_site[x0][y0]);
				const double proba_flip=exp(-4*beta*double(rhos-rhosp-1)/double(rho))*dt;
				
				if (proba_hop+proba_flip>1)
				{
					cerr << "THE PROBABILITY TO WAIT IS NEGATIVE: " << 1-(proba_hop+proba_flip) << endl;
					abort();
				}
				
				const double random_number=ran();
				//The particle hops: perform the hopping on the particle and modify the density/magneization (const spin).
				if (random_number<proba_hop)
				{
					POTTS[j].move(random_number/proba_hop,epsilon,LX,LY);
					
					//Update density and magnetization on old site.
					omp_set_lock(&lock_site[x0][y0]);
					RHOS[spin0][x0][y0]--;
					RHO[x0][y0]--;
					omp_unset_lock(&lock_site[x0][y0]);
					
					//Update density and magnetization on new site.
					const int xN=POTTS[j].x, yN=POTTS[j].y;
					omp_set_lock(&lock_site[xN][yN]);
					RHOS[spin0][xN][yN]++;
					RHO[xN][yN]++;
					omp_unset_lock(&lock_site[xN][yN]);
				}
				//The particle flips: perform the flipping on the particle and modify the magnetization (on-site, ind. of epsilon).
				else if (random_number<proba_hop+proba_flip)
				{
					POTTS[j].flip(spin1);
					
					//Update the magnetization.
					omp_set_lock(&lock_site[x0][y0]);
					RHOS[spin0][x0][y0]--;
					RHOS[spin1][x0][y0]++;
					omp_unset_lock(&lock_site[x0][y0]);
				}
				//Else do nothing (proba_wait)
			}			
		}		
	}
	return 0;
}
