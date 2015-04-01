//---------------------------------------------------------------------------

#ifndef amfftw_lnkH
#define amfftw_lnkH

extern int    IG_normalize, kernel_type;
extern float  cell_edge;

extern float  rat_min, rat_max, old_alpha;

// all of these are probably not worth exporting
// --------------- FFTW variables ----------------------
extern float **orig_LS, **conn_LS, **kernel;
//extern fftwf_complex *LS_fft, *kernel_fft;
//extern fftwf_plan plan_kernel_fwd, plan_LS_fwd, plan_LS_bwd;
extern int  dim_x, dim_y, LS_vec_len, fft_vec_len;
// -----------------------------------------------------

// mem required for FFTW smoothing kernels
float mem_required_smoothing_fftwf(int xd, int yd);

bool initialize_fftwf();
void uninitialize_fftwf();
void calc_fftwf_kernel(float alpha, int feat_num);
bool alloc_fftwf_matrixes();
void free_fftwf_matrixes();
void generate_fftwf_plans();
void free_fftwf_plans();
void copy_occ_to_fftwf_LS(int sp, int add_mode, float mult);
void copy_matrix_to_fftwf_LS(int add_mode, float **m);
void copy_matrix_to_fftwf_LS_edge(int add_mode, float **m);
void fftwf_connectivity();

//---------------------------------------------------------------------------
#endif
