/*function to calculate filling for non-interacting fg in lattice at given temp t and mu*/
/*for some reason math.h header is not automatically inked with the rest of std library, therefore
to compile, need to run
$gcc nfg.c -lm
to link to it.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 50.0 /*number of sites before periodic*/

double dispersion(double kx, double ky)
{
	return (-2*cos(kx) - 2*cos(ky));
}

double fermi(double kx, double ky, double mu, double beta)
{
	return 1.0/(exp(beta*(dispersion(kx,ky)-mu))+1);
}

int main(int argc, char *argv[])
{
	double mu,beta;
	if(argc != 3){
		printf("Wrong number of command line arguments. Expected mu and beta.\n");
		return 0;
	}
	else{
		sscanf(argv[1],"%lf",&mu);
		sscanf(argv[2],"%lf",&beta);
	}

	double *kx,*ky;
	kx = malloc(sizeof(double)*(int)N);
	ky = malloc(sizeof(double)*(int)N);

	int ii;
	for(ii=0;ii<(int)N;ii++){
		kx[ii] = 2*M_PI*(double)ii/N;
		ky[ii] = 2*M_PI*(double)ii/N;
	}

	double *kxgrid = malloc(sizeof(double)*(int)N*(int)N);
	double *kygrid = malloc(sizeof(double)*(int)N*(int)N);

	for(ii=0;ii<(int)N;ii++){
		kxgrid[ii] = kx[(int)floor((double)ii/N)];
		kygrid[ii] = ky[ii %(int)N];
	}

	double n;
	for(ii=0;ii<(int)N*(int)N;ii++){
		n += fermi(kxgrid[ii],kygrid[ii],mu,beta);
	}
	n = n/(double)N/(double)N;
	printf("For Mu = %0.2f t, T = %0.2f t, n = %0.2f \n",mu,beta,n);

	free(kx);
	free(ky);
	free(kxgrid);
	free(kygrid);
	return 0;
}
