#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

//////////// BLAS functions for vectors ///////////////

// Computes X = alpha * X, alpha scalar and X vector
extern void dscal_(int* elems,double* alpha,double* X,int* xspace);
// Computes norm 2 of vector A
extern double dnrm2_(int* elems, double* A,int* aspace);
// Computes dot product of vectors A and B
extern double ddot_(int* elems, double* A, int* aspace, double* B, int* bspace);
// computes X = alpha*A + X with A and X vectors.
extern void daxpy_(int* elems, double* alpha, double* A, int* aspace,double* X,int* xspace);
///////////// Flags /////////////////////////

int DEBUG = 0;
int PRINTS = 0;

//////////// Global variables ///////////////

int M = 20;			// Mesh size
int sizeA = -1;		// A will be sizeA x sizeA
double L = 10;		// Length of the scheme (along the x axis)
double H = 5;		// Height of the scheme (along the y axis)

double deltaX, deltaY; // Distances between mesh nodes

int P = 86400; 				// Time period considered
int N = 20;			// Time discretization 

double Tmax = 20;	// Maximum temperature over the period
double Tmin = 8;	// Minimum temperature over the period

//////////// Thermal conductivity related variables ///////////////
//////////// For air :
double lambd_air = 20.0e-6;
double rho_air = 1;
double cp_air = 1;
//////////// For ground :
double lambd_ground = 0.25e-6;
double rho_ground = 1;
double cp_ground = 1;

////////////////// Data structures /////////////////////
/*
	Represents a vector of size N.
*/
typedef struct{
	int N;
	double* values;
}vector_t, *vector;
/*
	Represents a triple i,j,k. Used for matrix entries : a_{i,j} = k
*/
typedef struct ijk {  	// Contains one triplet (index, index, value)
	int i,j;			// Indices
	double k;			// Value
} ijk ;
/*
	Represents an MxN sparse matrix containing nnz elements. alloc is the number of entries in
	data. The matrix is stored using a COO representation, to do this we use the struct ijk to
	represent each entry.
*/
typedef struct COO {	// Compressed matrix
	int M,N;
	int nnz; 		// Number of non zeros in the matrix
	int alloc; 		// Number of entries that have been allocated
	int compressed; 	// Status : sorted/compressed or not
	ijk *data;			// Content of the matrix
} COO;

// Used for time measurements 
long timeval_diff(struct timeval *t2, struct timeval *t1)
{
  long diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
  return (diff);
}

// some utilities for the data structures

int compareijk (const void *_a,const void *_b) {
	ijk *a = (ijk*)_a;
	ijk *b = (ijk*)_b;
	if (a->i < b->i) return -1;
	if (a->i > b->i) return 1;
	if (a->j < b->j) return -1;
	if (a->j > b->j) return 1;
	return 0;
}

int compareijk2 (ijk* a, ijk* b) {
	if (a->i < b->i) return -1;
	if (a->i > b->i) return 1;
	if (a->j < b->j) return -1;
	if (a->j > b->j) return 1;
	return 0;
}

void allocateCOO ( COO** coo, size_t initial, size_t M,size_t N) {
	*coo = (COO*) malloc (sizeof(COO));
	(*coo)-> compressed = 0;
	(*coo)-> nnz = 0;
	(*coo)-> alloc = initial;
	(*coo)-> data = (ijk*) malloc ((*coo)->alloc*sizeof(ijk));
	(*coo)->M = M;
	(*coo)->N = N;
}

void freeCOO ( COO** coo) {
	if (*coo == NULL) return;
	free ((*coo)->data);
	free (*coo);
	*coo = NULL;;
}

void addToCoo (COO *coo, int I, int J, double K) {
	if(I >= coo->M || J>=coo->N){printf("WARNING : Inserting out of bound element in matrix : I=%d, J= %d, K=%f, MxN = %ux%u\n",I,J,K ,coo->M,coo->N);}
	// if max size is reached: we double the size of the memory allocated
	if (coo-> nnz == coo->alloc){
		ijk *temp = (ijk*) malloc (2*coo->alloc*sizeof(ijk));
		for (size_t i=0; i<coo->alloc;i++){
			temp[i].i = coo->data[i].i;
			temp[i].j = coo->data[i].j;
			temp[i].k = coo->data[i].k;
		}
		free (coo->data);
		coo->data = temp;
		coo->alloc *= 2;
	}

	coo-> compressed = 0;
	coo->data[coo->nnz].i =  I;
	coo->data[coo->nnz].j =  J;
	coo->data[coo->nnz++].k =  K;
}

// you should probably complete the following functions

// sort lexicographically then compress, unused
void sortAndCompress (COO *coo) {
    if (coo-> compressed == 1) return;

    qsort( (void*) coo->data , coo->nnz , sizeof(ijk) , compareijk); //sort first

    for(int i = 0 ; i < coo->nnz-1 ; i++){
    	if( compareijk2( (coo->data+i), (coo->data + i + 1)) == 0 ){
    		coo->data[i].k += coo->data[i+1].k ;
    		for( int j = i+1 ; j < coo->nnz-1 ; j++){
    			coo->data[j] = coo->data[j+1] ;
    		}
    		coo->nnz--;
    	}
    }

    coo-> compressed = 1;
}
// Binary search for COO matrix
ijk* findijk(COO* coo,int i,int j){
	ijk key;
	key.i = i;
	key.j = j;
	void* r = bsearch((void*) &key, (void*) coo->data, coo->nnz , sizeof(ijk), compareijk);
	if(r==NULL) { return NULL; }
	return (ijk*) r;
}

// All these functions are mandatory and implemented below the main !
void COOFullCopy(COO* src, COO* dest);
void COOICCopy(COO* src,COO* dest);
void printCOO(COO* coo);
void printVector(vector v);
void printVecor2(double* v,int M);
void COOToCSV(COO* coo,char* name);
void solToCSV(double* T,char* name);
void buildA(COO* coo);
double nrmA(COO* A);
double nrmILUL(COO* ILU);
double nrmILUU(COO* ILU);

/*
	Computes R := alpha * Av + beta * R, if beta = 0 R will be set to 0 elementwise.
	1-D Arrays are double* arguments
*/
void dcooemv(double alpha,COO* A, double* V, double beta, double* R){
	int ONE = 1;
	// If beta = 0 sets R to zero.
	if( beta == 0 ){  memset(R, 0, A->M * sizeof(double));  }
	// Otherwise applies dscal to R except if beta = 1.
	else if( beta != 1){  dscal_( &A->M, &beta, R, &ONE);  }

	int i = 0;
	for(; i < A->nnz ; i++ ){
		*(R + (A->data[i]).i) += alpha * ( (A->data[i]).k * *(V + (A->data[i]).j ) );
		//*(R->values + *(A->row + i)) += alpha * *(A->values + i) * *(V->values + *(A->col + i));
	}
}

void ILU(COO* A){

	ijk* diag = A->data;
	ijk* runnerdiag = diag;

	double coef = 0;

	int rowcurr;
	int rowdiag;

	ijk* curr = A->data; // First of second line
	// Go get the first :
	while(curr->i == (curr+1)->i){
		curr++;
	}
	curr++;
	ijk* runner = curr;

	for(int i = 1 ; i < A->M ; i++ ){ // Do M-1 times (almost all rows) ==> LE I

		while(curr->i > curr->j){	// tant que pas sur la diag ==> LE K
			// compute coef
			diag = findijk(A , curr->j , curr->j ); // diag always exists, binary search 
			if(diag==NULL){ printf("BIG ERROR ILU0, DIAG DOES NOT EXISTS WTF\n");exit(-1); }

			runnerdiag = diag+1;
			runner = curr+1;

			coef = curr->k / diag->k ;
			curr->k = coef;

			rowcurr = curr->i;
			rowdiag = diag->i;
			
			// faire le reste de la ligne ==> LE J
			while(runner->i == rowcurr && runnerdiag->i == rowdiag){ // Tant que les deux runner sont toujours sur leur ligne respective
				if( runnerdiag->j < runner->j){  runnerdiag++;  } // Si diag en retard
				else if(runner->j < runnerdiag->j){  runner++;  } // Si la ligne en retard
				// Si ils sont alignés alors update
				else {
					runner->k -= coef * runnerdiag->k ;
					runnerdiag++;
					runner++;
				}
			}
			curr++;
		}
		// METTRE CURR SUR LA LIGNE SUIVANTE
		while(curr->i == (curr+1)->i){
			curr++;
		}
		curr++;
	}
}

/*
	IC(0), adapted Cholesky from the book to work with COO and take advantage 
	of it.
*/
void IC(COO* A){
	ijk* currdiag = A->data; // First element, assume diagonal always exists
	ijk* runnerdiag = currdiag+1;
	ijk* rkj = runnerdiag;

	double sqrtdiag;

	ijk* currcol ;
	ijk* runnercol;

	// For each line
	for(int i = 0 ; i < A->M ; i++){

		sqrtdiag = sqrt(currdiag->k);

		// For j in the book :
		// We only have to edit a line if the R_k,j is non null
		//  so we iterate over all R_k,j thanks to the COO format.
		// For each non-null Rkj :
		while(rkj->i == currdiag->i){
			// Looking for a diag, as it always exists
			currcol = findijk(A,rkj->j,rkj->j);
			runnercol = currcol;

			runnerdiag = currdiag+1;
			// Edit the line according to Rkj and make the matching between 
			//  the line's and the diagonal's line nonzeros:
			while(runnercol->i == currcol->i && runnerdiag->i == currdiag->i){
				if(runnerdiag->j < runnercol->j) {  runnerdiag++;  }
				else if(runnercol->j < runnerdiag->j) {  runnercol++;  }
				else{
					runnercol->k -= runnerdiag->k * (rkj->k / currdiag->k );
					runnercol++;
					runnerdiag++;
				}
			}
			// Last line in the For k in the book:
			// R_k,k:m /= sqrt(R_k,k)
			// We can do this now instead of outside the loop as we won't use
			//  r_k,j afterwards.
			rkj->k = rkj->k / sqrtdiag;
			rkj++;
		}
		currdiag->k = currdiag->k / sqrtdiag;
		// We have to move currdiag to the next one :
		runnerdiag = currdiag+1;
		while(runnerdiag->i != runnerdiag->j){
			runnerdiag++;
		}
		currdiag = runnerdiag;
		runnerdiag++;
		rkj = runnerdiag;
	}
}

// solve LU x = b
void ILUSolve(COO *LU, double *x, double *b) {
    // 1 : solve L(Ux)=b for Ux
    double* Ux = calloc(LU->M , sizeof(double));
    if(Ux==NULL){ printf("Malloc error for Ux in ILUSolve.\n");exit(-1); }

    // Starting from the top
    ijk* runner = LU->data;

    // Variable to hold computations using discovered values of x
    double temp = 0;

    // For each line except the last one to avoid out of bounds
    for(int i = 0 ; i < LU->M - 1; i++){
    	// Up to the diag and make sure we don't go to the next line
    	while( runner->j < runner->i && runner->i == i ){
    		temp += runner->k * Ux[runner->j];
    		runner++;
    	}
    	// Diagonal element will match with the next unknow, in the case of L it's one
    	
    	Ux[runner->i] = b[runner->i]-temp;
    	temp = 0;

    	// Get the runner to the next line
    	// ATTENTION SEGFAULT LAST ITER
    	while(runner->i == i){
    		runner++;
    	}
    }
    // Up to the diag and make sure we don't go to the next line
	while( runner->j < runner->i && runner->i == LU->M - 1 ){
		temp += runner->k * Ux[runner->j];
		runner++;
	}
	// Diagonal element will match with the next unknow, in the case of L it's one
	Ux[runner->i] = b[runner->i]-temp;
	temp = 0;

    // 2 : Solve Ux = y for x, and y is the Ux variable.

    // Set the runner on the last element
    runner = LU->data + LU->nnz - 1;
    // Up to the second line to avoid out of bounds
    for(int i = LU->M - 1 ; i > 0 ; i--){
    	// Up to the diag and make sure we don't go to the next line
    	while( runner->j > runner->i && runner->i == i ){
    		temp += runner->k * x[runner->j];
    		runner--;
    	}

    	// Diagonal element will match with the next unknow
    	x[runner->i] = (Ux[runner->i]-temp)/runner->k;

    	temp = 0;

    	// Get the runner to the next line
    	while(runner->i == i){
    		runner--;
    	}
    }
    // First line threatment
	while( runner->j > runner->i ){
		temp += runner->k * x[runner->j];
		runner--;
	}
	// Diagonal element will match with the next unknow
	x[runner->i] = (Ux[runner->i]-temp)/runner->k;

	free(Ux);
}

void ICSolve(COO *IC, double *x, double *b){
	// TODO
	// NOTE : IC in upper part of matrix -> ignore lower part
	// Solve LL^T x = b, first solve L (L^Tx) = b :
	double* LTx = calloc(IC->M,sizeof(double));
	if(LTx==NULL){ printf("Calloc error for LTx.\n");exit(-1); }

	double* bcopy = calloc(IC->M,sizeof(double));
	if(bcopy==NULL){ printf("Calloc error for bcopy.\n");exit(-1); }
	memcpy(bcopy,b,IC->M * sizeof(double));

	ijk* runner = IC->data;

	for(int i = 0 ; i < IC->M - 1 ; i++){
		// Get runner on diagonal
		while(runner->i > runner->j){
			runner++;
		}

		// Solve for this diag elem
		*(LTx + runner->i) = *(bcopy + runner->i) / runner->k;
		runner++;
		// Update all values
		while(runner->i == i){
			*(bcopy + runner->j) -= runner->k * *(LTx + runner->i);
			runner++;
		}
	}
	// Last value is updated here 
	*(LTx + IC->M-1) = *(bcopy+IC->M-1) / IC->data[IC->nnz-1].k ;

	free(bcopy);

	// 2 : Solve L^T x = y where y = LTx variable
	double temp = 0;
	runner = IC->data + (IC->nnz - 1);
	for(int i = IC->M - 1 ; i > 0 ; i--){
		// Get older coefficients
		while(runner->j > runner->i && runner->i==i){
			temp+= runner->k * x[runner->j];
			runner--;
		}

		// Diagonal element so new unknown to discover
		x[runner->i] = (LTx[runner->i] - temp) / runner->k;
		temp = 0;

		// Go to the next line
		while(runner->i == i){
			runner--;
		}
	}
	// First line
	while(runner->j > runner->i){
		temp+= runner->k * x[runner->j];
		runner--;
	}
	// Diagonal element matches last unknown
	x[runner->i] = (LTx[runner->i]-temp)/runner->k;
	free(LTx);
}

double f(double t){
	return ((Tmax + Tmin)/2) + (((Tmax-Tmin)/2) * sin(2 * M_PI * t/P));
}

// Returns b_0 if we read the scheme line by line.
void createB(vector b){
	b->N = M * (2*M+1);
	b->values = malloc(b->N * sizeof(double));
	if(b->values == NULL){
		printf("Malloc error in createB.\n");
		exit(-1);
	}

	double dt = ((double)P)/N;
	double c = f(0)/(((double)P)/N);   // TODO : set back to f(0)

	for(int i = 0 ; i < b->N - 2*M - 1 ; i++){
		b->values[i] = c;
	}
	double Dground = lambd_ground / (rho_ground * cp_ground);
	double last = c + Dground * (f(dt)/pow(deltaY,2));
	for(int i = b->N - 2*M - 1 ; i < b->N ; i++){
		b->values[i] = last ;
	}
}

// preconditioned conjugate gradients  for solving Ax = b , PREC ~= A^{-1}
// PRECTYPE : 0 = none (so PREC=NULL) , 1 = ILU0 , 2 = IC0
// return 0 if OK
// return 1 if max iterations attained
// return 2 if division by 0
int  PCG (COO *A, COO *PREC, int PRECTYPE, double *b, double *x,short Xnull, double precision) {

	// Variables for BLAS
	int ONE = 1;
	double dONE = 1;

	double* r = malloc(A->M * sizeof(double));
	if(r==NULL){ printf("Malloc error for r in PCG.\n");exit(-1); }
	double* z = malloc(A->M * sizeof(double));
	if(z==NULL){ printf("Malloc error for z in PCG.\n");exit(-1); }
	double* p = malloc(A->M * sizeof(double));
	if(p==NULL){ printf("Malloc error for p in PCG.\n");exit(-1); }
	double alpha;
	double alphaneg;
	double beta;

	// r_0 = b ==> r_0 = b = b-Ax_0 where x_0 = 0
	memcpy((void*) r, (void*) b, A->M*sizeof(double));
	if(!Xnull){
		dcooemv(-1, A, x, 1, r);
	}
	

	// z_0 = M^-1 r_0   ==> solve(M z_0 = r_0)  NOTE : z_0 = r~_0
	switch (PRECTYPE){
		case 0:
			// No prec => M = I => M^-1 = I => z_0 = I r_0 = r_0
			memcpy((void*) z , (void*) r , A->M*sizeof(double));
			break;
		case 1:
			ILUSolve( PREC , z , r );
			break;
		case 2:
			ICSolve( PREC , z , r );
			break;
		default:
			// No prec => M = I => M^-1 = I => z_0 = I r_0 = r_0
			memcpy((void*) z , (void*) r , A->M*sizeof(double));
	}

	// p_0 = z_0
	memcpy((void*) p , (void*) z , A->M*sizeof(double));

	// Temporary vector for computing alpha
	double* Ap = malloc(A->M * sizeof(double));
	if(Ap==NULL){ printf("Malloc error for Ap in PCG.\n");exit(-1); }
	// Temporary value for p^T A p
	double pTAp;
	// Temporary value for r_k^T z_k
	double rTz;
	// Temporary value for r_k+1^T z_k+1
	double next_rTz = ddot_( &A->M , r , &ONE , z , &ONE);

	double nrmresidus = dnrm2_(&A->M , r , &ONE);
	double nrmr_0 = nrmresidus;

	int iter = 0;

	while(iter < A->M / 3){
		/////////////////////// alpha = <r,z> / <p,Ap> ///////////////////////
		// Ap = A*p
		dcooemv(1,A,p,0,Ap);
		// printf("ITER : %d : Ap[0] = %f\n",iter,*(Ap+100) );
		// dot product : <p,Ap>
		pTAp = ddot_( &A->M , p , &ONE , Ap , &ONE);
		// dot product : <r,z>
		rTz = next_rTz;
		// And finally :
		alpha = rTz/pTAp;

		///////////////////////// x = x + alpha * p //////////////////////////
		daxpy_( &A->M , &alpha , p , &ONE , x , &ONE);

		/////////////////////// r = r - alpha * A * p ////////////////////////
		alphaneg = -1*alpha;
		daxpy_( &A->M , &alphaneg , Ap , &ONE , r , &ONE);

		///////////////////// Stopping condition on r ////////////////////////
		nrmresidus = dnrm2_(&A->M , r , &ONE);

		if(nrmresidus/nrmr_0 < 1e-12){
			break;
		}

		/////////////////////// z = M^-1 r <==> Mz = r ///////////////////////
		switch (PRECTYPE){
			case 0:
				// No prec => M = I => M^-1 = I => z_0 = I r_0 = r_0
				memcpy((void*) z , (void*) r , A->M*sizeof(double));
				break;
			case 1:
				ILUSolve( PREC , z , r );
				break;
			case 2:
				ICSolve( PREC , z , r );
				break;
			default:
				// No prec => M = I => M^-1 = I => z_0 = I r_0 = r_0
				memcpy((void*) z , (void*) r , A->M*sizeof(double));
		}

		////////////// beta = < r_k+1 , z_k+1 > / < r_k , z_k > //////////////
		// Compute < r_k+1 , z_k+1 > 
		next_rTz = ddot_( &A->M , r , &ONE , z , &ONE);
		// < r_k , z_k > Has already been computed either before the loop, or 
		// 				 during last iteration.
		beta = next_rTz / rTz ;

		////////////////////////// p = z + beta * p //////////////////////////
		// First p = beta * p
		dscal_( &A->M , &beta , p , &ONE );
		// Then p = z + p
		daxpy_( &A->M , &dONE , z , &ONE , p , &ONE);

		iter++;

	}

	free(r);
	free(z);
	free(p);
	free(Ap);
	return iter;
}

void updateB(double* b,int iter){

	double invdt = ((double)N)/P;

	int ONE = 1;
	//////////////////////////// b = 1/dt * b ////////////////////////////
	dscal_( &sizeA , &invdt , b , &ONE);

	double Dground = lambd_ground / (rho_ground * cp_ground);

	double updator = Dground * (f(iter/invdt)/pow(deltaY,2));
	///////////////////////// Last line update ///////////////////////////
	for(int i = (M-1)*(2*M+1) ; i < sizeA ; i++){
		b[i] += updator ;
	}
}

int timeIter(COO *A, COO *PREC, int PRECTYPE, int repeat, double *b, double *x, double precision,char* name){

	int Xcentre = round(5.0/deltaX);
	int Ycentre = round(2.0/deltaY);
	int indexCentre = Ycentre*(2*M+1) + Xcentre;

	double dt = ((double)P)/N;

	double TminC = 1000;
	double TmaxC = -1000;

	double* Tc = malloc((repeat*N + 1)*sizeof(double));
	Tc[0] = f(0);

	short swap = 0;
	
	for(int i = 1 ; i <= repeat*N; i++){

		if(swap){

			PCG(A,PREC, PRECTYPE, x, b,0, precision);


			Tc[i] = b[indexCentre];
			if(i >= (repeat-1)*N){
				if(b[indexCentre] > TmaxC){  TmaxC = b[indexCentre];  }
				if(b[indexCentre] < TminC){  TminC = b[indexCentre];  }
			}

			updateB(b,i+1);

			swap = 0;
		} else {
			// Always first
			PCG(A,PREC, PRECTYPE, b, x,0, precision);

			Tc[i] = x[indexCentre];
			if(i >= (repeat-1)*N){
				if(x[indexCentre] > TmaxC){  TmaxC = x[indexCentre];  }
				if(x[indexCentre] < TminC){  TminC = x[indexCentre];  }
			}

			updateB(x,i+1);
			swap = 1;
		}
	}

	FILE* output = fopen(name,"w");
	if(output == NULL){printf("Failure with output file, name given : %s\n",name);exit(-1);}

	fprintf(output, "%d %.16f\n", repeat*N + 1 ,(TmaxC-TminC)/(Tmax-Tmin));
	for(int i = 0 ; i <=  repeat*N  ; i++){
		fprintf(output, "%.16f %.16f\n",i*dt,Tc[i] );
	}

	fclose(output);

	if(DEBUG){
		// printVecor2(iters,N);
		printf("Centre du tuyau : min = %.16f, max = %.16f\n",TminC,TmaxC);	
		printf("Efficacité = (TmaxC - TminC)/(Tmax - Tmin) = %.16f\n",(TmaxC-TminC)/(Tmax-Tmin));
	}
	
	free(Tc);

	return 0;
}

int main(int argc, char *argv[]) {

	int opt;
	int params = argc;						// Argument manager for -d & -v
	while((opt=getopt(argc,argv,"dve")) != -1){
		switch (opt){
			case 'd':
				DEBUG = 1;
				printf("============== Debug mode ==============\n");
				params--;
				break;
			case 'v':
				PRINTS = 1;
				printf("============== Print (verbose) mode ==============\n");
				params--;
				break;
			default:
				break;
		}
	}

	if(params != 7){printf("Not enough or too many arguments, use ./program DorY nb M N prec fichier.out [-d] [-v]\nWhere -d and -v are optional flags used for debugging.\n");exit(-1);}

	M = atoi(argv[optind+2]); 
	N = atoi(argv[optind+3]);
	int PRECTYPE = atoi(argv[optind+4]);  
	char* OUT = argv[optind+5];
	int repeat = atoi(argv[optind+1]);
	int day = atoi(argv[optind]);

	double precision = 1e-12;

	if(DEBUG){printf("ARGS : M = %d, N = %d, PRECTYPE = %d, OUT = %s, repeat = %d, day = %d\n",M,N,PRECTYPE,OUT,repeat,day );}

	// A nnz formula :
	// 2 lower corners : 6 , 2 upper corners : 6 , left & right side : 8*(M-2) , bottom & top : 8*(2*M -1) , center : 5*(M-2) * (2*M-1)
	// if M = 20 : nnz = 3978
	// M = 20 , A = 20x41 => 820 lines (note : from 1, not 0)
	int expectednnz = 12+ 8*(M-2) + 8*(2*M -1) + 5*(M-2) * (2*M-1);

	sizeA = (M)*(2*M+1);
	deltaX = ((double)L)/(2*M);
	deltaY = ((double)H)/M;

	if(!day){
		P = 31536000;     // secs in one year
		Tmin = -5;
		Tmax = 25;
	}

	COO* A;

	allocateCOO(&A,expectednnz,sizeA,sizeA);
	buildA(A);

	COO* PREC = malloc(sizeof(COO));
	if(PREC==NULL){ printf("Malloc error for PREC\n");exit(-1); }
	if(PRECTYPE == 1){    // ILU
		COOFullCopy(A,PREC);
		ILU(PREC);
	} else if(PRECTYPE==2){
		COOFullCopy(A,PREC);
		IC(PREC);
	}

	vector b_0 = malloc(sizeof(vector));
	if(b_0==NULL){ printf("Malloc error for b_0\n");exit(-1); }
	createB(b_0);


	double* T = calloc(sizeA,sizeof(double));
	if(T==NULL){ printf("Malloc error for T\n");exit(-1); }

	timeIter(A,PREC,PRECTYPE,repeat,b_0->values,T,precision,OUT);

	free(T);
	free(b_0->values);
	if(PREC){
		free(PREC);
	}
	freeCOO(&A);
	return 0;
}

//////////////// Helpers for the assignment and BuildA deep below /////////////
void COOToCSV(COO* coo,char* name){
	FILE* output = fopen(name,"w");
	if(output == NULL){printf("Failure with output file, name given : %s\n",name);exit(-1);}

	fprintf(output, "%d,%d\n",coo->M,coo->N );
	for(int i = 0 ; i < coo->nnz ; i++ ){
		fprintf(output, "%d,%d,%.16f\n", coo->data[i].i, coo->data[i].j, coo->data[i].k);
	}

	fclose(output);
}

void solToCSV(double* T,char* name){
	FILE* output = fopen(name,"w");
	if(output == NULL){printf("Failure with output file, name given : %s\n",name);exit(-1);}

	fprintf(output, "%d\n",M*(2*M+1) );
	for(int i = 0 ; i < M*(2*M+1)-1 ; i++ ){
		fprintf(output, "%.16f,", T[i]);
	}
	fprintf(output, "%.16f\n", T[M*(2*M+1)-1]);
	fclose(output);
}

void vectorToCSV(double* V, int len, char* name){

}

void printCOO(COO* coo){
	//int curr = coo->data[0].i;
	printf("i ,  j   ,   k\n");
	for(int index = 0 ; index < coo->nnz ; index++){
		if(index != coo->nnz - 1 && coo->data[index].i != coo->data[index+1].i){
			printf(" %u  ,  %u  ,  % 2.10f ===\n", coo->data[index].i , coo->data[index].j , coo->data[index].k );
		} else {
			printf(" %u  ,  %u  ,  % 2.10f\n", coo->data[index].i , coo->data[index].j , coo->data[index].k );
		}
	}
}



void printVector(vector v){
	for(int i = 0 ; i < v->N ; i++){
		printf("% 4d  : % 2.15f\n",i,v->values[i] );
	}
}

void printVecor2(double* v,int M){
	for(int i = 0 ; i < M ; i++){
		printf("% 4d  : % 2.18f\n",i,v[i] );
	}
}

double nrmILUL(COO* ILU){
	double tot = 0;
	tot += ILU->M;
	for(int i = 0 ; i < ILU->nnz ; i++){
		if(ILU->data[i].i > ILU->data[i].j){
			tot+=ILU->data[i].k;
		}
	}
	return sqrt(tot);
}

double nrmILUU(COO* ILU){
	double tot = 0;
	for(int i = 0 ; i < ILU->nnz ; i++){
		if(ILU->data[i].i <= ILU->data[i].j){
			tot+=ILU->data[i].k;
		}
	}
	return sqrt(tot);
}

double nrmA(COO* A){
	double tot = 0;
    for (int i = 0; i < A->nnz; ++i) {
        double k = A->data[i].k;
        tot += k*k;
    }
    return sqrt(tot);
}

void COOFullCopy(COO* src, COO* dest){
	dest->M = src->M;
	dest->N = src->N;
	dest->data = malloc(src->nnz*sizeof(ijk));
	if(dest->data == NULL){ printf("Malloc error in COOFullCopy.\n");exit(-1); }
	dest->alloc = src->nnz;
	dest->nnz = src->nnz;
	dest->compressed = src->compressed;
	for(int i = 0 ; i < dest->nnz ; i++){
		(dest->data + i)->i =(src->data + i)->i;
		(dest->data + i)->j =(src->data + i)->j;
		(dest->data + i)->k =(src->data + i)->k;
	}
}

// src assumed symmetric
void COOICCopy(COO* src,COO* dest){
	dest->M = src->M;
	dest->N = src->N;
	dest->data = malloc(((src->nnz + src->M)/2  )*sizeof(ijk));
	if(dest->data == NULL){ printf("Malloc error in COOFullCopy.\n");exit(-1); }
	dest->alloc = ((src->nnz + src->M)/2 );
	dest->nnz = ((src->nnz + src->M)/2 );
	dest->compressed = src->compressed;
	int inddest = 0;
	for(int i = 0 ; i < src->nnz ; i++){
		if((src->data + i)->i <= (src->data + i)->j){
			(dest->data + inddest)->i =(src->data + i)->i;
			(dest->data + inddest)->j =(src->data + i)->j;
			(dest->data + inddest)->k =(src->data + i)->k;
			inddest++;
		}
	}
}

/*
	Builds the system matrix A, yes it took a long time to write, it was not such
	a good idea to do it in order to avoid sorting and compressing .
*/
void buildA(COO* coo){

	double dx2 = pow( L/(2*M) , 2);
	double dy2 = pow( H/M , 2);

	double invdx2 = 1 / dx2;
	double invdy2 = 1 / dy2;

	int pipeInX = round(4.75 / deltaX);
	int pipeOutX = round(5.25 / deltaX);

	int pipeInY = round(1.75 / deltaY);
	int pipeOutY = round(2.25 / deltaY);

	double Dair = lambd_air / (rho_air * cp_air);
	double Dground = lambd_ground / (rho_ground * cp_ground);

	double Dmoy = (Dground + Dair) / 2 ;

	double dt = ((double)N)/P;

	double Ddx = -1 * (Dground/dx2);
	double Ddy = -1 * (Dground/dy2);

	double Dairdx = -1 * (Dair/dx2);
	double Dairdy = -1 * (Dair/dy2);

	double Dmoydx = -1 * (Dmoy/dx2);
	double Dmoydy = -1 * (Dmoy/dy2);

	double diag1 = dt + Dground*(1/dx2 + 1/dy2);
	double diag2 = dt + Dground*(1/dx2 + 2/dy2);
	double diag3 = dt + Dground*(2/dx2 + 1/dy2);
	double diag4 = dt + Dground*(2/dx2 + 2/dy2);

	int line = 0;
	// First line of matrix, lower left corner of scheme
	addToCoo(coo , line , line , diag1);
	addToCoo(coo , line , line+1 , Ddx);
	addToCoo(coo , line , line + 2 * M + 1 , Ddy );
	line++;
	for( ; line < 2*M ; line++){
		addToCoo(coo , line , line - 1 , Ddx);
		addToCoo(coo , line , line , diag3 );
		addToCoo(coo , line , line + 1 , Ddx);
		addToCoo(coo , line , line + 2 * M + 1 , Ddy );
	}
	// Line 2*M + 1, lower right corner of the scheme
	addToCoo(coo , line , line - 1 , Ddx );
	addToCoo(coo , line , line , diag1 );
	addToCoo(coo , line , line + 2*M + 1 , Ddy );
	line++;
	// Lines up to the line of the bottom of the pipe
	for( ; line < pipeInY * (2*M + 1); line++ ){
		// Left Neumann condition
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line , diag2 );
		addToCoo(coo , line , line + 1 , Ddx );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
		line++;
		// Center of the scheme
		for(int i = 0 ; i < 2*M - 1 ; i++){
			addToCoo(coo , line , line - 2*M - 1 , Ddy );
			addToCoo(coo , line , line - 1 , Ddx );
			addToCoo(coo , line , line , diag4 );
			addToCoo(coo , line , line + 1 , Ddx );
			addToCoo(coo , line , line + 2*M + 1 , Ddy );
			line++;
		}
		// Right Neumann condition
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line - 1 , Ddx );
		addToCoo(coo , line , line , diag2 );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
	}
	// First Node, Neumann condition
	addToCoo(coo , line , line - 2*M - 1 , Ddy );
	addToCoo(coo , line , line , diag2 );
	addToCoo(coo , line , line + 1 , Ddx );
	addToCoo(coo , line , line + 2*M + 1 , Ddy );
	line++;
	// Center until the pipe
	for(int i = 0 ; i < pipeInX - 1 ; i++){
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line - 1 , Ddx );
		addToCoo(coo , line , line , diag4 );
		addToCoo(coo , line , line + 1 , Ddx );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
		line++;
	}
	// Bottom of the pipe begins ============================
	// Lower left corner of the pipe
	addToCoo(coo , line , line - 2*M - 1 , Ddy);
	addToCoo(coo , line , line - 1 , Ddx );
	addToCoo(coo , line , line , dt + (Dmoy + Dground) * (invdx2 + invdy2) );
	addToCoo(coo , line , line + 1 , Dmoydx );
	addToCoo(coo , line , line + 2*M + 1 , Dmoydy);
	line++;
	// Bottom of the pipe excluding the lower right corner
	for(int i = 0 ; i < (pipeOutX - pipeInX) - 1 ; i++ ){
		addToCoo(coo , line , line - 2*M - 1 , Ddy);
		addToCoo(coo , line , line - 1 , Dmoydx);
		addToCoo(coo , line , line , dt - 2* Dmoydx - Dairdy - Ddy  );
		addToCoo(coo , line , line + 1 , Dmoydx);
		addToCoo(coo , line , line + 2*M + 1 , Dairdy);
		line++;
	}
	// Lower right corner of the pipe
	addToCoo(coo , line , line - 2*M - 1 , Ddy);
	addToCoo(coo , line , line - 1 , Dmoydx );
	addToCoo(coo , line , line , dt + (Dmoy + Dground) * (invdx2 + invdy2));
	addToCoo(coo , line , line + 1 , Ddx );
	addToCoo(coo , line , line + 2*M + 1 , Dmoydy);
	line++;
	// Gets to the right side border
	for( int i = 0 ; i < (2*M + 1) - pipeOutX - 2 ; i++ ){
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line - 1 , Ddx );
		addToCoo(coo , line , line , diag4 );
		addToCoo(coo , line , line + 1 , Ddx );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
		line++;
	}
	addToCoo(coo , line , line - 2*M - 1 , Ddy );
	addToCoo(coo , line , line - 1 , Ddx );
	addToCoo(coo , line , line , diag2 );
	addToCoo(coo , line , line + 2*M + 1 , Ddy );
	line++;
	// Up to the upper left corner of the pipe
	for( int i = 0 ; i < pipeOutY - pipeInY - 1 ; i++){
		// Left Neumann CL
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line , diag2 );
		addToCoo(coo , line , line + 1 , Ddx );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
		line++;
		// up until pipe
		for( int j = 0 ; j < pipeInX - 1 ; j++){
			addToCoo(coo , line , line - 2*M - 1 , Ddy );
			addToCoo(coo , line , line - 1 , Ddx );
			addToCoo(coo , line , line , diag4 );
			addToCoo(coo , line , line + 1 , Ddx );
			addToCoo(coo , line , line + 2*M + 1 , Ddy );
			line++;
		}
		// Pipe's left side
		addToCoo(coo , line , line - 2*M - 1 , Dmoydy);
		addToCoo(coo , line , line - 1 , Ddx );
		addToCoo(coo , line , line , dt - 2 * Dmoydy - Ddx - Dairdx);
		addToCoo(coo , line , line + 1 , Dairdx );
		addToCoo(coo , line , line + 2*M + 1 , Dmoydy);
		line++;
		// Inside the pipe
		for( int j = 0 ; j < pipeOutX - pipeInX - 1 ; j++){
			addToCoo(coo , line , line - 2*M - 1 , Dairdy );
			addToCoo(coo , line , line - 1 , Dairdx );
			addToCoo(coo , line , line , dt + 2*Dair * (invdx2 + invdy2) );
			addToCoo(coo , line , line + 1 , Dairdx );
			addToCoo(coo , line , line + 2*M + 1 , Dairdy );
			line++;
		}
		// Pipe's right side
		addToCoo(coo , line , line - 2*M - 1 , Dmoydy);
		addToCoo(coo , line , line - 1 , Dairdx );
		addToCoo(coo , line , line , dt - 2 * Dmoydy - Ddx - Dairdx);
		addToCoo(coo , line , line + 1 , Ddx );
		addToCoo(coo , line , line + 2*M + 1 , Dmoydy);
		line++;
		// Up until right side border
		for( int j = 0 ; j < 2*M - pipeOutX - 1 ; j++){
			addToCoo(coo , line , line - 2*M - 1 , Ddy );
			addToCoo(coo , line , line - 1 , Ddx );
			addToCoo(coo , line , line , diag4 );
			addToCoo(coo , line , line + 1 , Ddx );
			addToCoo(coo , line , line + 2*M + 1 , Ddy );
			line++;
		}
		// Right side, Neumann CL
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line - 1 , Ddx );
		addToCoo(coo , line , line , diag2 );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
		line++;
		//====================================================================
	}
	addToCoo(coo , line , line - 2*M - 1 , Ddy );
	addToCoo(coo , line , line , diag2 );
	addToCoo(coo , line , line + 1 , Ddx );
	addToCoo(coo , line , line + 2*M + 1 , Ddy );
	line++;
	for(int i = 0; i < pipeInX - 1 ; i++){
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line - 1 , Ddx );
		addToCoo(coo , line , line , diag4 );
		addToCoo(coo , line , line + 1 , Ddx );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
		line++;
	}
	// Upper left corner of the pipe
	addToCoo(coo , line , line - 2*M - 1 , Dmoydy);
	addToCoo(coo , line , line - 1 , Ddx );
	addToCoo(coo , line , line , dt + (Dmoy + Dground) * (invdx2 + invdy2) );
	addToCoo(coo , line , line + 1 , Dmoydx );
	addToCoo(coo , line , line + 2*M + 1 , Ddy);
	line++;
	// Top side of the pipe excluding upper right corner
	for(int i = 0 ; i < (pipeOutX - pipeInX) - 1 ; i++ ){
		addToCoo(coo , line , line - 2*M - 1 , Dairdy);
		addToCoo(coo , line , line - 1 , Dmoydx);
		addToCoo(coo , line , line , dt - 2* Dmoydx - Dairdy - Ddy );
		addToCoo(coo , line , line + 1 , Dmoydx);
		addToCoo(coo , line , line + 2*M + 1 , Ddy);
		line++;
	}
	// Upper right corner of the pipe
	addToCoo(coo , line , line - 2*M - 1 , Dmoydy);
	addToCoo(coo , line , line - 1 , Dmoydx );
	addToCoo(coo , line , line , dt + (Dmoy + Dground) * (invdx2 + invdy2) );
	addToCoo(coo , line , line + 1 , Ddx );
	addToCoo(coo , line , line + 2*M + 1 , Ddy);
	line++;
	// Up to the scheme's right side
	for( int j = 0 ; j < 2*M - pipeOutX - 1 ; j++){
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line - 1 , Ddx );
		addToCoo(coo , line , line , diag4 );
		addToCoo(coo , line , line + 1 , Ddx );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
		line++;
	}
	// Right side, Neumann CL
	addToCoo(coo , line , line - 2*M - 1 , Ddy );
	addToCoo(coo , line , line - 1 , Ddx );
	addToCoo(coo , line , line , diag2 );
	addToCoo(coo , line , line + 2*M + 1 , Ddy );
	line++;
	// Up until last line
	for( ; line < (M-1)*(2*M + 1) ; line++){
		// Left Neumann condition
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line , diag2 );
		addToCoo(coo , line , line + 1 , Ddx );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
		line++;
		// Center of the scheme
		for(int i = 0 ; i < 2*M - 1 ; i++){
			addToCoo(coo , line , line - 2*M - 1 , Ddy );
			addToCoo(coo , line , line - 1 , Ddx );
			addToCoo(coo , line , line , diag4 );
			addToCoo(coo , line , line + 1 , Ddx );
			addToCoo(coo , line , line + 2*M + 1 , Ddy );
			line++;
		}
		// Right Neumann condition
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line - 1 , Ddx );
		addToCoo(coo , line , line , diag2 );
		addToCoo(coo , line , line + 2*M + 1 , Ddy );
	}
	// Second to last line (not Dirichlet)
	// Last Left Neumann condition before Dirichlet
	addToCoo(coo , line , line - 2*M - 1 , Ddy );
	addToCoo(coo , line , line , diag2 );
	addToCoo(coo , line , line + 1 , Ddx );
	line++;
	// Center of the scheme but last line before Dirichlet
	for(int i = 0 ; i < 2*M - 1 ; i++){
		addToCoo(coo , line , line - 2*M - 1 , Ddy );
		addToCoo(coo , line , line - 1 , Ddx );
		addToCoo(coo , line , line , diag4 );
		addToCoo(coo , line , line + 1 , Ddx );
		line++;
	}
	//  LAST Right Neumann condition
	addToCoo(coo , line , line - 2*M - 1 , Ddy );
	addToCoo(coo , line , line - 1 , Ddx );
	addToCoo(coo , line , line , diag2 );

	coo->compressed = 1;

	if(PRINTS){ printCOO(coo); }

	if(DEBUG){
		printf("1/dx² = %.4f, 1/dy² = %.4f, 1/dt = %.8f, Dground = %e, Dair = %e\n", invdx2,invdy2,dt,Dground,Dair);
		printf("GROUND : diag1 = %.4f, diag2 = %.4f, diag3 = %.4f, diag4 = %.4f\n",diag1,diag2,diag3,diag4 );
		printf("M*(2M+1) = %d\n",(M)*(2*M+1) );
		printf("Index : X : pipeIn = %d, pipeOut = %d ; Y : pipeIn = %d, pipeOut = %d\n",pipeInX, pipeOutX ,pipeInY,pipeOutY);
		printf("Line number : In X : %d , Y : %d .   Out X : %d , Y : %d\n", pipeInY*(2*M+1) + pipeInX, 0 , pipeOutY*(2*M + 1) + pipeOutX , 0  );
	}
}
