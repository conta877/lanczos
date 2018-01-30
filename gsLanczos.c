#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>


double irand(unsigned int seed) {
	double random;
	srand(seed);
	random = (double)rand() / (double)RAND_MAX;
	return(random);
}


double complex icrand(unsigned int seed1, unsigned int seed2) {
	double m_pi = 3.141592653589793238;
	double random1, random2;
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
	int terminate = 1E9;
	double Normtmp = 0;
	double Norme = 0;
	double err = 100;
	double CONV = 1E-10;
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
		tmpvector[i] = irand(123456);
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




double gslsz(double complex *val, unsigned long long int *valloc, unsigned long long int nzsize, unsigned long long int index, int lszsteps, char *name) {
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
	FILE *sneakpeak;
	char string[1000];



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
	r = (double complex*)calloc(index, sizeof(double complex));
	w = (double complex *)calloc(index, sizeof(double complex));
	if (r == NULL) {
		printf("cannot allocate %llu to r\n", (unsigned long long int)index);
		return(-999);
	}
	if (w == NULL) {
		printf("cannot allocate %llu to w\n", (unsigned long long int)index);
		return(-999);
	}

	double alpha[lszsteps];
	double beta[lszsteps - 1];
	double evector[lszsteps];
	double tmpvector[lszsteps];

	//CREATING RANDOM INITIAL VECTOR W
	double Norm = 0;
	for (j = 0; j<index; j++) {
		w[j] = icrand(12345,22222);
		// printf("%lf\t%lf\n",creal(w[j]),cimag(w[j]));
		Norm += (pow(creal(w[j]), 2.0) + pow(cimag(w[j]), 2.0));
		if (j < lszsteps) {  //init the td hamiltonian vectors
			alpha[j] = 0;
			if (j < lszsteps - 1) beta[j] = 0;
			evector[j] = 0;
			tmpvector[j] = 0;
		}
	}
	Norm = sqrt(Norm);
	// printf("created random initial vector");

	//STORE IF EIGENVECTOR NEEDED
	double complex *w0wf;
	int twice = 0;
	if (name != "0") {
		w0wf = (double complex *)calloc(index, sizeof(double complex));
		twice = 1;
	} //declare the array to keep the first random

	  //NORMALIZE THE RANDOM VECTOR
	for (j = 0; j<index; j++) { //normalize the random vector, also init the other
		r[j] = 0.0;
		w[j] = w[j] / Norm;
		if (name != "0") w0wf[j] = w[j]; //keep the first random to restart for the wf
	}

	// printf("normalized the random initial vector");


	//LANCZOS 

	double eigen0 = -1000;
	int kkeep;
	double complex store;
	double CONV = 1.0E-6;

	//Nit=0 calculate eigenvalue, Nit=1 calculate eigenvector
	for (Nit = 0;Nit <= twice;Nit++) {

		//IF SECOND ROUND, INITIALIZE
		if (Nit == 1) {
			k = 0;
			lszsteps = kkeep;
			for (j = 0;j < index;j++) {
				r[j] = 0.0;
				w[j] = w0wf[j]; //set w0 to stored random
				w0wf[j] = w[j] * ((double complex)evector[k]); //start calculating the gs vector
			}
		}

		//LANCZOS CORE
		// printf ("number lszsteps: %d\n",lszsteps);
		for (k = 0; k < lszsteps; k++) {
		   //printf("'k: %d\n",k);
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
				r[i] += val[p] * w[j]; //first round it was zero, after it is -beta_(k-1)v_(k-1)
				alpha[k] += conj(w[i])*val[p] * w[j];
				if (i != j) { //lower triangle!
					r[j] += conj(val[p])*w[i];
					alpha[k] += conj(w[j])*conj(val[p])*w[i];
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
					r[i] = r[i] - w[i] * alpha[k];
					beta[k] += (pow(creal(r[i]), 2.0) + pow(cimag(r[i]), 2.0));   //beta=norm(r)
				}//i
				beta[k] = (double)sqrt(beta[k]);     // beta = sqrt(norm(r))  

				for (i = 0; i < index; i++) {              //(w+1)=r/norm(r), r=beta.w
					store = w[i];
					w[i] = r[i] / ((double complex)beta[k]);
					r[i] = -beta[k] * store;
					if (Nit == 1) {
						w0wf[i] = w0wf[i] + w[i] * ((double complex)evector[k + 1]); //SECOND ROUND GS VECTOR
					}
				} //i
			} //END BETA
			// printf("beta=%lf\n", beta[k]);

			//EIGENVALUE
			if (Nit == 0 && (lszsteps == 1 || k==lszsteps-1 || k%10 == 9)) {
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
					k = lszsteps+10; //terminate the loop
				}
				eigen0 = eigen;
				printf("eigen = %lf, round = %d\n",eigen+EIGMAX,k);
			} //END EIGENVALUE - ONLY FIRST ROUND
		} //k LANCZOS CORE END

		// printf("eigen=%lf, round %d\n", eigen, kkeep);
		  //EIGENVECTOR PRINT
		if (Nit == 1) {
//			sprintf(string, "%s.evector", name);
			sprintf(string, "GS.evector");
			gsvector = fopen(string, "w");
			sprintf(string, "sneakpeak_GS.evector");
			sneakpeak = fopen(string, "w");

	    	fprintf(sneakpeak,"%d %llu\n",lszsteps,index);
	        fprintf(sneakpeak,"Energy of GS: %lf\n",eigen+EIGMAX);

	    	fprintf(gsvector,"%d %llu\n",lszsteps,index);
	        fprintf(gsvector,"Energy of GS: %lf\n",eigen+EIGMAX);
			for (i = 0;i < index;i++) {
				fprintf(gsvector, "%llu\t%lf\t%lf\t%lf\n", i, creal(w0wf[i]), cimag(w0wf[i]), cabs(w0wf[i]));
				if (cabs(w0wf[i]) >=0.02){
					fprintf(sneakpeak, "%llu\t%lf\t%lf\t%lf\n", i, creal(w0wf[i]), cimag(w0wf[i]), cabs(w0wf[i]));
				}
			}
			fclose(gsvector);
			fclose(sneakpeak);
		} //END EIGENVECTOR
	// printf("Nit=%d over\n",Nit);
	}//ROUNDS Nit

	free(r);
	free(w);
	if (name != "0") free(w0wf);
	// printf ("returning eigenvalue:%lf\n",eigen);
	if (lszsteps != 1) return(eigen + EIGMAX);
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
	printf("into the loop to read\n");
	for (i = 0;i < nzsize;i++) {
		tmp = fscanf(readfile, "%llu %lf %lf", &valloc[i], &reaval, &comval);
		val[i] = reaval + comval*I;
	}
	fclose(readfile);
	printf("lets go to lanczos\n");
	clock_t t1;
	t1 = clock();

	eigen = gslsz(val, valloc, nzsize, index, lszstep, "test");
	printf("eigenvalue=%.8lf\n", eigen);
	free(val);
	free(valloc);
	t1 = clock() - t1;
	printf("diagonalized in %f secs.\n\n", ((float)t1) / CLOCKS_PER_SEC);

	return(0);
}
