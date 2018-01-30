#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include "lapacke.h"
#include "lapacke_utils.h"

int _PRINT_WF=10;

double irand(unsigned int seed) {
	clock_t t1;
	double random;
	time_t rawtime;
	if (seed==0){
		t1 = clock();
		time(&rawtime);
		seed = (int)(rawtime / t1);
	}
	srand(seed);
	random = (double)rand() / (double)RAND_MAX;
	return(random);
}

double complex icrand(unsigned int seed1, unsigned int seed2) {
	clock_t t1;
	double random1;
	double random2;
	time_t rawtime;

	if (seed1==0){
		t1 = clock();
		time(&rawtime);
		seed1 = (int)(rawtime / t1);
	}
	if (seed2==0){
		t1 = clock();
		time(&rawtime);
		seed2 = (int)(rawtime / t1);
	}
	srand(seed1);
	random1 = (double)rand() / (double)RAND_MAX;
	srand(seed2);
	random2 = (double)rand() / (double)RAND_MAX;

	return(random1 + I*random2);
}


double powsymtri(double *alpha, double *beta, int msize, double *evector, double *tmpvector) {
	/*   alpha = digonals (msize)
	beta = diagonal+-1 (msize-1)
	msize = matrix size = size of alpha array */
	int i;
	int step = 0;
	int terminate = 1E6;
	double Normtmp = 0;
	double Norme = 0;
	double err = 100;
	double CONV = 1E-12;
	double Eigen, Eigen0;
	Eigen = 100;

	/*start with guess normalized vector x[i]
	START LOOP
	y = H.x
	y[i] = beta[i - 1]x[i - 1] + alpha[i]x[i] + beta[i]x[i + 1]
	Eigen = x.H.x = x.y
	x = y
	END LOOP
	finally y will converge to the true eigenvector */

	for (i = 0;i < msize;i++) {
		tmpvector[i] = irand(0);
		Normtmp += tmpvector[i] * tmpvector[i];
	} //generate random initial guess vector
	Normtmp = sqrt(Normtmp);

	//	while (err > CONV && step < terminate) {
	while (err > CONV) {
		Norme = 0;
		err = 0;
		Eigen0 = Eigen;
		Eigen = 0;
		for (i = 0;i < msize;i++) {

			if (step != 0) tmpvector[i] = evector[i];
			tmpvector[i] /= Normtmp;

			evector[i] = alpha[i] * tmpvector[i];

			if (i != 0) evector[i] += beta[i - 1] * tmpvector[i - 1];

			if (i != msize - 1){
				tmpvector[i+1] = evector[i+1];
				tmpvector[i+1] /= Normtmp;
				evector[i] += beta[i] * tmpvector[i + 1];
			} 

			Eigen += evector[i] * tmpvector[i];
			Norme += evector[i] * evector[i];
		} //generate new unnormalized vector

		Normtmp = sqrt(Norme);
		if (Eigen < 0) Normtmp *= -1.0000; //to ensure converge of negative eigens. otherwise the vectors get messed up
		err = fabs(Eigen - Eigen0);
		step++; //distinguish the first run
	}

	if (step >= terminate) return(-999); //NOT CONVERGED!

	for (i = 0;i < msize;i++) {
		evector[i] /= Normtmp; //normalize the last vector
	}
	return(Eigen);
}

double multiStateLsz(double complex *val, unsigned long long int *valloc, unsigned long long int nzsize, unsigned long long int index, int lszsteps, double *eigenvalues, char *name) {
	//val contains upper triangle nonzeros
	//valloc points to their location
	//nzsize is the number of nonzeros
	//index is the size of the hamiltonian
	//lszsteps in the size of the tridiagonal matrix to be calculated
	//name[] is the extension of the filename for wf. if NULL, then do not print wf's
	//returns the gs eigenvalue in case of error returns -999

	if (lszsteps>index) {
		printf("ERROR, cannot do lanczos iterations greater than the size of my hamiltonian\n");
		return(-999);
	}

	int k;
	unsigned long long int i, j, p;
	int Nit;
	FILE *gsvector;
	char string[1000];
        double complex inrprd;
	int info;
        double complex *gsH;
        double *z;



	//FIND EIGMAX _ THE GROUND STATE APPROXIMATION TO SHIFT THE VALUES //
	double eigen;
	double EIGMAX = 0;
	if (lszsteps != 1) {
		eigen = gslsz(val, valloc, nzsize, index, 1, "0");
		EIGMAX = abs(eigen) + 1000.0;
		if (eigen == -999) { //making sure everything is ok
			printf("ERROR in finding EIGMAX\nEXITING...\n");
			return(-999);
		}
		printf("EIGMAX=%lf\n", EIGMAX);
	}



	//ALLOCATE ARRAYS
	double complex *r;
	double complex *w;
	double complex *gramSchmidtShift;
	int saidItAlready = 0;
	r = (double complex*)calloc(index, sizeof(double complex));
	gramSchmidtShift = (double complex*)calloc(index, sizeof(double complex));
	z=(double *)calloc(lszsteps*lszsteps,sizeof(double));
	w = (double complex *)calloc((unsigned long long int)(lszsteps)*index, sizeof(double complex)); //IF DOING CONVERGE, IS THIS NECESSARY LIKE THIS?
	if (r == NULL) {
		printf("cannot allocate %llu to r\n", (unsigned long long int)index);
		return(-999);
	}
	if (w == NULL) {
		printf("cannot allocate %llu to w\n", (unsigned long long int)(lszsteps)*index);
		return(-999);
	}
	if (z == NULL) {
		printf("cannot allocate %d to z\n", (lszsteps)*(lszsteps));
		return(-999);
	}

	double alpha[lszsteps];
	double beta[lszsteps - 1];
	double evector[lszsteps];
	double tmpvector[lszsteps];

	double alphakeep[lszsteps];
        double betakeep[lszsteps-1];

	//CREATING RANDOM INITIAL VECTOR W
	double Norm = 0;
	for (j = 0; j<index; j++) {
		w[j] = icrand(0,0);
		Norm += (pow(creal(w[j]), 2.0) + pow(cimag(w[j]), 2.0));
		if (j < lszsteps) {  //init the td hamiltonian vectors
			alpha[j] = 0;
			if (j < lszsteps - 1) beta[j] = 0;
			evector[j] = 0; //?
			tmpvector[j] = 0; //?
		}
	}
	Norm = sqrt(Norm);
	for (j = 0; j<index; j++){ //NORMALIZE
		r[j] = 0.0;
		w[j] = w[j] / Norm;
	}
	printf("created random initial vector\n");

	//LANCZOS 

	// these are for gs convergence, maybe same model extend to excited states
	double eigen0 = -10000;
	double eigen1 = -10000;
	int kkeep=0;
	double complex store;
	double CONV = 1.0E-6;

	//LANCZOS CORE
	for (k = 0; k < lszsteps; k++) {
		alpha[k] = 0.0;
		//ALPHA
		for (p = 0; p < nzsize; p++) {    //Create A=H.w(n)-beta(n-1).w(n-1) and alpha(n)=w(n).H.w(n)
										  //determine the location of the element val[p]
			i = (int)((unsigned long long int)(valloc[p]) / ((unsigned long long int)(index)));
			j = (int)((unsigned long long int)(valloc[p]) - ((unsigned long long int)(index)*(unsigned long long int)(i)));

			//Error handling
			if (i < 0 || j < 0) {
				printf("valloc=%lld ==> i=%lld,j=%lld\nERROR with element locations\nEXITING LANCZOS", valloc[p], i, j);
				return(-999);
			}
			if (j < i) { //only stored upper triangle, requires j>i
				printf("ERROR...Nonzeros are filled wrong, j<i\t%lld>%lld  <===  %lld\n", i, j, valloc[p]);
				return(-999);
			}

			if (i == j && lszsteps != 1) { //diagonal elements, shift to make GS the greatest
				val[p] -= (EIGMAX);
			}
			r[i] += val[p] * w[k*index + j];
			alpha[k] += conj(w[k*index + i])*val[p] * w[k*index + j];
			if (i != j) { //lower triangle!
				r[j] += conj(val[p])*w[k*index + i];
				alpha[k] += conj(w[k*index + j])*conj(val[p])*w[k*index + i];
			}
			if (i == j && lszsteps != 1) { //revert the Hamiltonian to normal
				val[p] += (EIGMAX);
			}
		}//p END ALPHA
	 	// printf("alpha=%lf\t", alpha[k]);

		//BETA
		if (k < lszsteps - 1) {
			beta[k] = 0;
			for (i = 0; i < index; i++) {              // r=(H.w-(w-1).(beta-1))-w.alpha
				gramSchmidtShift[i]=0;
				r[i] = r[i] - w[k*index+i] * alpha[k];
				beta[k] += (pow(creal(r[i]), 2.0) + pow(cimag(r[i]), 2.0));   //beta=norm(r)
			}//i
			beta[k] = (double)sqrt(beta[k]);     // beta = sqrt(norm(r))  

			for (i = 0; i < index; i++) {              //(w+1)=r/norm(r), r=beta.w
				store = w[k*index+i];
				w[(k+1)*index+i] = r[i] / ((double complex)beta[k]);
				r[i] = -beta[k] * store;
			} //i
			// GRAM-SCHMIDT 
			for(j=0;j<=k;j++){ // for each previous state
				inrprd=0.0; 
				for(i=0;i<index;i++){    //<wj.wk+1>
	  				inrprd+=conj(w[j*index+i])*w[(k+1)*index+i];
				}	
				for(i=0;i<index;i++){   //<wj.wk+1>.wj
	  				gramSchmidtShift[i]+=inrprd*w[j*index+i];
				}
      		}	
		    Norm=0;
      		for(i=0;i<index;i++){   //wk+1-Sum_j(<wj.wk+1>.wj)
				w[(k+1)*index+i]-=gramSchmidtShift[i];
				Norm+=(pow(creal(w[(k+1)*index+i]),2.0)+pow(cimag(w[(k+1)*index+i]),2.0));
      		}
      		Norm=sqrt(Norm);
      		for(i=0;i<index;i++){ //normalize wk+1
				w[(k+1)*index+i]=w[(k+1)*index+i]/((double complex)(Norm));
      		}
		} //END BETA
	    // printf("beta=%lf\n", beta[k]);

		//EIGENVALUE
		if (lszsteps == 1 || k==lszsteps-1 || k%10 == 9) {
			// printf("Nit=%d, lszsteps=%d, k=%d, kkeep=%d, eigen0=%lf\n",Nit,lszsteps,k,kkeep,eigen0);
			eigen = powsymtri(alpha, beta, lszsteps, evector, tmpvector);
			if (k==lszsteps-1) {
				kkeep=k+1; // calculation is already over, must exit
			}
			if (eigen == -999) {
				printf("Error in diagonalization!\n");
				return(-999);
			}
			if (fabs(eigen - eigen0) < CONV) {
				kkeep = k+1; //keep the max lszsteps required
				if (saidItAlready==0){ 
					printf("ground state energy converged\n");
					saidItAlready=1;
				}
			}

			eigen0 = eigen;
			if (saidItAlready<1) printf("eigen = %lf, round = %d\n",eigen+EIGMAX,k);
			if (k%50==49){
				for(int copyCount=0;copyCount<lszsteps;copyCount++){
					alphakeep[copyCount] = alpha[copyCount];
					if (copyCount<lszsteps-1) betakeep[copyCount] = beta[copyCount];
				}
		  		info=LAPACKE_dsterf(lszsteps, alphakeep, betakeep);
		  		printf ("round %d: first 10 states:\n",k);
				for(int printCount=0;printCount<20;printCount++){
					printf("%d:%lf\n",printCount,alphakeep[printCount]+EIGMAX);
				}
				if (fabs(alphakeep[1] - eigen1) < CONV) {
					if (saidItAlready==1){
					    printf("first excited state energy converged\n");
					    saidItAlready=2;
				    }
				}
				eigen1 = alphakeep[1];
			}
		}
	} //k LANCZOS CORE END
	free(r);
  	free(gramSchmidtShift);
	// printf ("returning eigenvalue:%lf\n",eigen);
	if (lszsteps != 1){
		    //    printf("in GS printing\n");
    	FILE *sneakPeak;
    	FILE *wfs;
    	sneakPeak=fopen("sneak_peak.txt","w");
    	gsH=(double complex*)calloc(index,sizeof(double complex));
		for(i=0;i<lszsteps;i++){
      		z[i*lszsteps+i]=1.0;
    	}
        info=LAPACKE_dsteqr(101,'I',lszsteps,alpha,beta,z,lszsteps);
    	printf("done with diagonalization %d\n",info);
    	fprintf(sneakPeak,"%d %llu\n",lszsteps,index);
	    for(k=0;k<lszsteps;k++){
	    	if (k<_PRINT_WF){
		    	sprintf(string,"wf_st%d.txt",k);
	   	 		wfs = fopen(string,"w");
	    		fprintf(wfs,"%d %llu\n",lszsteps,index);
	        	fprintf(wfs,"Energy of %d: %lf\n",k,alpha[k]+EIGMAX);
	        }
	    	fprintf(sneakPeak,"Energy of %d: %lf\t\n",k,alpha[k]+EIGMAX);
	  	    for(i=0;i<index;i++){
				gsH[i]=0;
				for(j=0;j<lszsteps;j++){
			  		gsH[i]+=w[j*index+i]*((double complex)(z[j*lszsteps+k]));
				}
				if (k<_PRINT_WF) fprintf(wfs,"%llu\t%lf\t%lf\t%lf\n",i,creal(gsH[i]),cimag(gsH[i]),cabs(gsH[i]));      
		  		if(pow(cabs(gsH[i]),2)>=0.05){
		  			fprintf(sneakPeak,"%d: %llu\t{%lf %lf}\t%lf\n",k,i,creal(gsH[i]),cimag(gsH[i]),pow(cabs(gsH[i]),2));
		  		}
			}
			if (k<_PRINT_WF) fclose(wfs);
	    } // k
		fclose(sneakPeak);
		free(z);
	    free(gsH);
		free(w);
		for(int k=0;k<20;k++){
			printf("%d:%lf\n",k,alpha[k]+EIGMAX);
		}
		return(alpha[0]+EIGMAX);
	} 
	if (lszsteps == 1) return(eigen); //SHIFT
}



int main(void) {
	FILE *readfile;
	double reaval, comval;

	readfile = fopen("hamiltonian.0", "r");
	printf("file opened\n");
	int i;
	int tmp;
	unsigned long long int nzsize;
	unsigned long long int index;
	tmp = fscanf(readfile, "%llu", &nzsize);
	tmp = fscanf(readfile, "%llu", &index);
	printf("number of nonzeros: %llu, size of Hamiltonian: %llu\n", nzsize,index);

	double complex *val;
	unsigned long long int *valloc;
	val = (double complex*)calloc(nzsize, sizeof(double complex));
	valloc = (unsigned long long int *)calloc(nzsize, sizeof(unsigned long long int));
	int lszstep;
	tmp = fscanf(readfile, "%d", &lszstep);
	double eigen;
	printf("number of lanczos steps: %d\ninto the loop to read\n",lszstep);
	for (i = 0;i < nzsize;i++) {
		tmp = fscanf(readfile, "%llu %lf %lf", &valloc[i], &reaval, &comval);
		val[i] = reaval + comval*I;
	}
	fclose(readfile);
	printf("lets go to lanczos\n");
	clock_t t1;
	t1 = clock();
	double *eigenvalues;
	eigenvalues=(double*)calloc(lszstep, sizeof(double));
	eigen = multiStateLsz(val,valloc, nzsize, index, lszstep, eigenvalues, "test");
	printf("eigenvalue=%.8lf\n", eigen);

	free(val);
	free(valloc);
	t1 = clock() - t1;
	printf("diagonalized in %f secs.\n\n", ((float)t1) / CLOCKS_PER_SEC);
	free(eigenvalues);
	return(0);
}
