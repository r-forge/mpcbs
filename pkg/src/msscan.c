#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Print.h>


//This function assumes that y[i,j] = i + y_dim[1] * j
//If this is not the case, the matrix probably just needs to be transposed before
//being passed in.  The return value will be a vector indexed in the same way.

SEXP ComputeZ(SEXP p_y, SEXP p_T, SEXP p_win, SEXP p_y_dim, SEXP p_alpha)
{
	//return value
	SEXP ret_dbl;

	//loop counters
	unsigned int i,j, k, k2;
	
	PROTECT(p_alpha = AS_NUMERIC(p_alpha));
	double alpha = NUMERIC_POINTER(p_alpha)[0];
	UNPROTECT(1);	

	PROTECT(p_T = AS_INTEGER(p_T));
	unsigned int T = INTEGER_POINTER(p_T)[0];
	UNPROTECT(1);

	PROTECT(p_win = AS_INTEGER(p_win));
	unsigned int win = INTEGER_POINTER(p_win)[0];
	UNPROTECT(1);

	PROTECT(p_y = AS_NUMERIC(p_y));
	double *y = NUMERIC_POINTER(p_y);
	//UNPROTECT(1);


	unsigned int y_dim[2];

	PROTECT(p_y_dim = AS_INTEGER(p_y_dim));
	y_dim[0] = INTEGER_POINTER(p_y_dim)[0];
	y_dim[1] = INTEGER_POINTER(p_y_dim)[1];
	UNPROTECT(1);



	//U=zeros(T,win);
    //Z=zeros(T,win);

	unsigned int Z_dim0 = T, Z_dim1 = win;
	unsigned int Z_sz = Z_dim0 * Z_dim1;
    double *U = malloc(sizeof(double)*Z_sz);

	PROTECT(ret_dbl = NEW_NUMERIC(Z_sz));  // Allocating storage space
    double *Z = NUMERIC_POINTER(ret_dbl);
    
	for(i = 0; i < T*win; i++) U[i] = Z[i] = 0.0;

    double *S = malloc(sizeof(double)*y_dim[1]);


	double log_alpha = log(alpha);
	double T_dbl = (double)T;
	double dfnum = 1.0;
	double dfden = T_dbl-2.0;

	for(i = 0; i < y_dim[0]; i++){

		// S=cumsum(y(i,:));
		//j = i*y_dim[1]; k = 0;

		S[0] = y[i];
		//for(k=1, j=j+1 ; k < y_dim[1]; k++, j++) S[k] = S[k-1] + y[j];
		for(k=1, j=(y_dim[0]+i) ; k < y_dim[1]; k++, j+= y_dim[0]) S[k] = S[k-1] + y[j];

		double S_last = S[T-1];

		//SST = sum((y(i,:)-S(T)/T).^2);
		double SST = 0.0;
		//j = i*y_dim[1];
		j = i;
		//for(k = 0; k < y_dim[1]; k++, j++){
		for(k = 0; k < y_dim[1]; k++, j+= y_dim[0]){
			double v1 = y[j] - S_last/T_dbl;

			SST += v1 * v1;
		}
		

		//for k=1:win
		double dbl_k;
		for(k = 0, dbl_k = 1.0; k < win; k++, dbl_k++){
			//un-vectorize the loop
			
			//To get the vector S(1:T-k) starting with k=1,2,...
			//the MATLAB loop would start k = 1, and the index 
			//set would be 1, 2, ..., T-1,  and then 1,2,...T-2, etc...

			//now we start k = 0, and index from 0,1,...,T-2
			//then 0,1,...,T-3, etc...

			//The vector S(k+1:T) would be indices 2,3,...T
			//and then 3,4,...,T, etc.
			//Now it will be 1,2,...,(T-2 + 1), and then 2,3,....(T-3+2), etc...
			for(k2 = 0; k2 < T-(k+1); k2++){

				//SSb = k*((S(k+1:T)-S(1:T-k))/k-S(T)/T).^2;
				double v1 = S[k2+k+1] - S[k2];
				double v2 = (v1/dbl_k - S_last/T_dbl);
				v2 *= dbl_k * v2;

				//SSb = SSb + (T-k)*((S(T)-(S(k+1:T)-S(1:T-k)))/(T-k)-S(T)/T).^2;
				double v3 = (S_last-v1)/(T_dbl-dbl_k) - S_last/T_dbl;
				double SSb = v2 + (T_dbl-dbl_k) * v3 * v3;

				//SSw = SST-SSb;
				double SSw = SST - SSb;

				//U(1:T-k,k) = (SSb/dfnum)./(SSw/dfden);

				//The convention is that  U[i,j] = i * u_dim[1] + j
				//U[ k2 * Z_dim1 + k ] = (SSb/dfnum)/(SSw/dfden);
				U[ k2 + Z_dim0 * k ] = (SSb/dfnum)/(SSw/dfden);
			}
		}

		double *Z2 = Z, *U2 = U;

		for(k = Z_sz; k; k--, Z2++, U2++){
			//g= @(u) u.*exp(u/2)./(ALPHA+exp(u/2));
			double exp_term = 1.0;
			if(*U2 < 10.0 + log_alpha){
				exp_term = exp((*U2) * 0.5);
				exp_term /= (exp_term + alpha);
			}

			(*Z2) += (*U2) * exp_term;
		};
	}

	free(S);
	free(U);

	UNPROTECT(2);

	//The last lines :
	//   Z = Z - N * psidot0
	//   Z=Z/sqrt(psidotdot0*N);
	//will be done in R code outside this function
	return( ret_dbl );
}


SEXP ChisqContrib(SEXP p_y, SEXP p_T, SEXP p_y_dim, SEXP p_leftpts, SEXP p_rightpts)
{
	int n_segs = length(p_leftpts);
	if(n_segs != length(p_rightpts)){
		;//return ( R_nilvalue);
	}

	//loop counters
	unsigned int i,j, k, k2;
	
	PROTECT(p_leftpts = AS_INTEGER(p_leftpts));
	int *leftpts = INTEGER_POINTER(p_leftpts);

	PROTECT(p_rightpts = AS_INTEGER(p_rightpts));
	int *rightpts = INTEGER_POINTER(p_rightpts);


	PROTECT(p_T = AS_INTEGER(p_T));
	unsigned int T = INTEGER_POINTER(p_T)[0];
	UNPROTECT(1);

	PROTECT(p_y = AS_NUMERIC(p_y));
	double *y = NUMERIC_POINTER(p_y);
	//UNPROTECT(1);

	unsigned int y_dim[2];

	PROTECT(p_y_dim = AS_INTEGER(p_y_dim));
	y_dim[0] = INTEGER_POINTER(p_y_dim)[0];
	y_dim[1] = INTEGER_POINTER(p_y_dim)[1];
	UNPROTECT(1);

	unsigned int N = y_dim[0];

	//return value
	SEXP ret_dbl;

	PROTECT(ret_dbl = NEW_NUMERIC(n_segs * y_dim[0]));  // Allocating storage space
    double *chisq_contr = NUMERIC_POINTER(ret_dbl);
    
	double T_dbl = (double)T;
	double dfnum = 1.0;
	double dfden = T_dbl-2.0;

	for(i = 0; i < n_segs; i++){

		int st=leftpts[i];
		int ed=rightpts[i];
		double w = ed-st+1;
    
		for(j = 0; j < N; j++){

			//SST = sum((y(j,:)-sum(y(j,:))/T).^2);

			double row_mean = 0.0;
			unsigned int yi = j;
			for(k = 0; k < y_dim[1]; k++, yi += y_dim[0]){
				row_mean += y[yi];
			}
			row_mean /= T_dbl;

			double SST = 0.0;
		
			for(yi = j, k = 0; k < y_dim[1]; k++, yi += y_dim[0]){
				double v2 = (y[yi] - row_mean);
				SST += v2 * v2;
			}
			

			//SSb=(sum(y(j,(st+1):ed))-w*sum(y(j,:))/T)^2/(w*(1-w/T));

			double SSb = 0.0;
			int ki;
			for(yi = j + st*y_dim[0], ki = st; ki < ed; ki++, yi += y_dim[0]){
				SSb += y[yi];
			}
			SSb = SSb - w * row_mean;
			SSb *= SSb;
			SSb /= (w*(1-w/T_dbl));
			
			//SSw = SST-SSb;
			double SSw = SST-SSb;


			//chisq(i,j) = (SSb/dfnum)/(SSw/dfden);

			chisq_contr[i + j*n_segs] = (SSb/dfnum)/(SSw/dfden);
		}
	}


	UNPROTECT(4);

	return( ret_dbl );
}





SEXP RemoveOverlap(SEXP p_leftpts, SEXP p_rightpts)
{
	int n_segs = length(p_leftpts);
	if(n_segs != length(p_rightpts)){
		;//return ( R_nilvalue);
	}

	//loop counters
	unsigned int i,j, k, k2;
	
	PROTECT(p_leftpts = AS_NUMERIC(p_leftpts));
	double *leftpts = NUMERIC_POINTER(p_leftpts);

	PROTECT(p_rightpts = AS_NUMERIC(p_rightpts));
	double *rightpts = NUMERIC_POINTER(p_rightpts);

	//return value
	SEXP ret_int;

	PROTECT(ret_int = NEW_INTEGER(n_segs));  // Allocating storage space
    int *segs = INTEGER_POINTER(ret_int);
    
	segs[0] = 1;
	for(i = 1; i < n_segs; i++) segs[i] = -1;

	unsigned int max_seg = 1;

	for(i = 1; i < n_segs; i++){
		unsigned int overlap=0;
		double st=leftpts[i];
		double ed=rightpts[i];

		for(j = 0; j < max_seg; j++){
			if(segs[j] == -1) continue;

			if(	(st <= leftpts[j] && ed >= rightpts[j]) || 
				(st >= leftpts[j] && st <= rightpts[j]) || 
				(ed >= leftpts[j] && ed <= rightpts[j]) )
			{
				overlap = 1;
				break;
			}
		}
		if(overlap == 0){
			segs[i] = 1;
			max_seg = i+1;
		}
	}


	UNPROTECT(3);

	return( ret_int );
}





