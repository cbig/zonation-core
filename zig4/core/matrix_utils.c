#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdint.h>

#define NR_END 1
#define FREE_ARG char*

void zig_legacy_alloc_error(char error_text[]);

float  **matrix( long nrl, long nrh, long ncl, long nch) {
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	float **m;

	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) {
	  zig_legacy_alloc_error("Allocation error 1 in matrix()");
	  free(m);
	  m=NULL;
	  return NULL;
	}

	//if (!m)
	//  Form1->Memo1->Lines->Add("Allocation error 1 in matrix");
	m += NR_END;
	m -= nrl;

	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl])
	{
		zig_legacy_alloc_error("Allocation error 2 in matrix()");
		free(m);
		m=0;
		return 0;
	}
	//if (!m[nrl])
	//  Form1->Memo1->Lines->Add("Allocation error 1 in matrix");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;
	return m;
}


char  **cmatrix( long nrl, long nrh, long ncl, long nch) {
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	char **m;

	m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char *)));
	if (!m) zig_legacy_alloc_error("Allocation error 1 in cmatrix()");
	m += NR_END;
	m -= nrl;

	m[nrl]=(char *) calloc(1, (size_t)((nrow*ncol+NR_END)*sizeof(char)));
	if (!m[nrl]) zig_legacy_alloc_error("Allocation error 2 in cmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;
	return m;
}

int    **imatrix(long nrl, long nrh, long ncl, long nch) {
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	int **m;

	m=(int **) calloc(1,(size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) zig_legacy_alloc_error("Allocation error 1 in imatrix()");
	m += NR_END;
	m -= nrl;

	m[nrl]=(int *) calloc(1, (size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) zig_legacy_alloc_error("Allocation error 2 in imatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;
	return m;
}

int64_t **int64_matrix(long nrl, long nrh, long ncl, long nch) {
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	int64_t **m;

	m=(int64_t **) malloc((size_t)((nrow+NR_END)*sizeof(int64_t*)));
	if (!m) zig_legacy_alloc_error("Allocation error 1 in imatrix()");
	m += NR_END;
	m -= nrl;

	m[nrl]=(int64_t *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int64_t)));
	if (!m[nrl]) zig_legacy_alloc_error("Allocation error 2 in int64_matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;
	return m;
}

float    **fmatrix(long nrl, long nrh, long ncl, long nch) {
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	float **m;

	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) zig_legacy_alloc_error("Allocation error 1 in fmatrix()");
	m += NR_END;
	m -= nrl;

	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) zig_legacy_alloc_error("Allocation error 2 in fmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;
	return m;
}

/************************** fpmatrix ***********************/

float  ***fpmatrix( long nrl, long nrh, long ncl, long nch) {
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	float ***m;

	m=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!m) {
	  zig_legacy_alloc_error("Allocation error 1 in matrix()");
	  free(m);
	  m=NULL;
	  return NULL;
	}
	//if (!m)
	//  Form1->Memo1->Lines->Add("Allocation error 1 in matrix");
	m += NR_END;
	m -= nrl;

	m[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float *)));
	if (!m[nrl]) {
	  zig_legacy_alloc_error("Allocation error 2 in fpmatrix()");
	  free(m);
	  m=NULL;
	  return NULL;
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;
	return m;
}

void   free_fpmatrix(float    ***m, long nrl, long nrh, long ncl, long nch) {
	if (m)
	{
		free((FREE_ARG) (m[nrl]+ncl-NR_END));
		free((FREE_ARG) (m+nrl-NR_END));
	}
	m = NULL;
}

void   free_matrix(float    **m, long nrl, long nrh, long ncl, long nch) {
	if (m)
	{
		free((FREE_ARG) (m[nrl]+ncl-NR_END));
		free((FREE_ARG) (m+nrl-NR_END));
	}
	m = NULL;
}


void   free_cmatrix(char    **m, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void   free_imatrix(int     **m, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void   free_fmatrix(int     **m, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
