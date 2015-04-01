#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <stdint.h> // for int64_t

extern "C" float  **matrix( long nrl, long nrh, long ncl, long nch);
extern "C" char   **cmatrix( long nrl, long nrh, long ncl, long nch);
extern "C" int    **imatrix(long nrl, long nrh, long ncl, long nch);
extern "C" int64_t **int64_matrix(long nrl, long nrh, long ncl, long nch);
extern "C" float  **fmatrix(long nrl, long nrh, long ncl, long nch);

extern "C" float  ***fpmatrix( long nrl, long nrh, long ncl, long nch);
extern "C" void   free_fpmatrix(float  ***m, long nrl, long nrh, long ncl, long nch);

extern "C" void   free_matrix(float    **m, long nrl, long nrh, long ncl, long nch);
extern "C" void   free_cmatrix(char    **m, long nrl, long nrh, long ncl, long nch);
extern "C" void   free_imatrix(int     **m, long nrl, long nrh, long ncl, long nch);
extern "C" void   free_fmatrix(int     **m, long nrl, long nrh, long ncl, long nch);

extern "C" void zig_legacy_alloc_error(char error_text[]);

#endif // MATRIX_UTILS_H
