#include "Unit1.h"
#include "GridMap.h"
#include "matrix_utils.h"
#include "fftw_lnk.h"
#include "VCL.h"
#include "bat_run.h"
#include <limits>
#include <ctime>
#include <fftw3.h>

int    IG_normalize=1, kernel_type=1;
float  cell_edge=1.0;

float  rat_min=100000, rat_max=-100000, old_alpha=-1;
int old_arb_kernel;

// --------------- FFTW variables ----------------------
float **orig_LS=0, **conn_LS=0, **kernel=0;
fftwf_complex *LS_fft=0, *kernel_fft=0;
fftwf_plan plan_kernel_fwd, plan_LS_fwd, plan_LS_bwd;
int  dim_x, dim_y, LS_vec_len, fft_vec_len;
// -----------------------------------------------------

// warning: this is too dependent on what is done in alloc_fftwf_matrixes() 
// but it has to be done before. - beware of possible future changes there!
float mem_required_smoothing_fftwf(int xd, int yd)
{
  float mem_req = 0;

  int d2_x = 2*xd;
  int d2_y = 2*yd;

  // 3 = orig_LS, conn_LS, and kernel
  mem_req += 3 * sizeof(typeof(**orig_LS))*d2_y*d2_x;

  // 2 = LS_fft, and kernel_fft  (both together should be ~equivalent to 1 of the above 3) - but x2 because fftwf_complex is made of 2 floats!
  fft_vec_len= max(d2_x,d2_y)*(min(d2_x,d2_y)/2+1);
  mem_req += 2 * fft_vec_len*sizeof(fftwf_complex);

  return mem_req;
}

bool initialize_fftwf()
{
#if !defined(__MINGW32__)
  return fftwf_init_threads() == 0 ? false : true;
#else
  return true;
#endif
}

void uninitialize_fftwf()
{
	//not needed, (there is actually some problems with releasing mutexes)
#if !defined(__MINGW32__)
	fftwf_cleanup_threads();
#endif
}

// kn: kernel number (assumed >= FIRST_NUM_ARB_KERNEL)
inline float
calc_arb_kernel_cell(Arb_Kernel& kernel, float dist, float alpha)
{
  float mult=1.0/(2.0*M_PI*alpha*alpha);
  return mult*exp(-alpha*dist);

  
  float extremum = 0.0;
  std::vector<float>& x = kernel.first;
  std::vector<float>& y = kernel.second;
  float slope = 0;
  // exception: 1 single point - 0,0 is implicit as initial point
  if (1 == x.size()) {
    if (0 == x[0]) {
      slope = 0;
    } else {
      slope = (y[0]-.0)/(x[0]-.0);
    }
    float val = slope * (dist-0);
    return val;
  }
  
  // More than 1 point:
  bool interv_found = false;
  float left = 0, right = 0;
  /* Sequential/linear search
  for (size_t i=1; i<x.size(); i++) {
    if (!interv_found && (i+1<x.size()) && (dist < x[i+1])) {
      left = y[i-1];
      right = y[i];
      interv_found = true;
      break;
    }
  }
  */
  // Binary search
  int idx_a = 0;
  int idx_b = x.size()-1;
  while (idx_b >= idx_a) {
    int mid = (idx_a + idx_b)/2;
    if (x[mid] <= dist && dist <= x[mid+1]) {
      interv_found = true;
      left = y[mid];
      right = y[mid+1];
      break;
    } else if (dist > x[mid+1]) {
      idx_a = mid+1;
    } else { // dist < x[mid]
      idx_b = mid-1;
    }
  }

  // 'dist' is higher than the highest x, extend last segment to infinity
  if (!interv_found) {
    right = y[y.size()-1];
    left = y[y.size()-2];
  }
  float val = left + slope * (right-left);
  // force >=0 values
  val = max(0.0f, val);
  return val;
}

// alpha is the alpha paremeter in the .spp file (features list)
// For individual species (distrib. smoothing), it is multiplied by the factor given 
// in the command line (second last parameter to the core)
void calc_fftwf_kernel(float alpha, int k_idx)
{
  // TODO: hard-wired tests
  //k_idx = 4;

  // If not particular kernel specified, use default:
  if (k_idx < 0)
    kernel_type = arb_kernels_default;

  if (alpha==old_alpha  
      && k_idx>0 && k_idx==old_arb_kernel) {
    Form1->Memo1->Lines->Add("Using same kernel as previous feature.");
    return;
  }

  Arb_Kernel* arb_kernel;
  if (k_idx >= FIRST_NUM_ARB_KERNEL) {
    Arbitrary_Kernels_Map::iterator it;
    it = arbitrary_kernels.find(k_idx);
    if (it == arbitrary_kernels.end()) {
      Form1->Memo1->Lines->Add("ERROR: arbitrary kernel not found (apparently it has not been loaded): "+IntToStr(k_idx));
      return;
    }
  }

  //Form1->Memo1->Lines->Add("USING KERNEL: "+IntToStr(k_idx));
  //Form1->Memo1->Lines->Add("calc_fftwf_kernel("+FloatToStr(alpha)+"): Kernel dim_x: "+IntToStr(dim_x)+", dim_y: "+IntToStr(dim_y));
  //Form1->Memo1->Lines->Add("calc_fftwf_kernel(), cell_edge: "+FloatToStr(cell_edge));

	cell_edge = obsmap[0].cell_size();
	//  Form1->Memo1->Lines->Add("CE="+FloatToStr(cell_edge));


	float mult=1.0/(2.0*M_PI*alpha*alpha);
	float gauss_mult = 1/(alpha*sqrt(2.0*M_PI));
	for(int y=0; y<dim_y; y++)
	{
		for(int x=0; x<dim_x; x++)
		{
			int dx = min(x, dim_x-x);
			int dy = min(y, dim_y-y);
			if ((dx>=(dim_x/2)) || (dy>=(dim_y/2)))
			{
				kernel[y][x] =0.0;
				continue;
			}
			float dist = sqrt(z_pow(cell_edge*dx, 2) + z_pow(cell_edge*dy ,2));
			if (k_idx<=1) {  // k_idx -1, 0, 1 are accepted as the traditional exponential decay
			  //         if ((alpha*dist)<35)
			  // traditional kernel
			  if (kernel_type==1)
			    kernel[y][x] = mult*exp(-alpha*dist);
			  else if (kernel_type==2)
			    kernel[y][x] =  mult*exp( -0.5*( z_pow(dx*cell_edge/alpha,2.0) +z_pow(dy*cell_edge/alpha,2.0) ) );
			  else if (kernel_type==3) {
			    if (dist<alpha)
			      kernel[y][x]=1.0;
			    else
			      kernel[y][x]=0.0;
			  }
			} else { // arbitrary kernels
			  if (k_idx >= FIRST_NUM_ARB_KERNEL) {
			    kernel[y][x] = calc_arb_kernel_cell(*arb_kernel, dist, alpha);
			  } else if (1 == k_idx) { // exponential
			    kernel[y][x] = mult*exp(-alpha*dist);
			  } else if (2 == k_idx) { // Gaussian
			    float dist_spread = (z_pow(cell_edge*dx, 2) + z_pow(cell_edge*dy ,2)) / (2.0*alpha*alpha);
			    kernel[y][x] = gauss_mult*exp(-dist_spread);
			  } else if (3 == k_idx) { // Threshold (square block)
			    if (dist<alpha)
			      kernel[y][x] = 1.0f;
			    else
			      kernel[y][x] = 0.0f;
			  } else if (4 == k_idx) { // Triangular
			    if (dist<alpha)
			      kernel[y][x] = 1.0f - dist/alpha;
			    else
			      kernel[y][x] = 0.0f;
			  } else { // default for unknown kernel: 0
			    kernel[y][x] = 0;
			  }
			}
		}
	}
	//  kernel[0][0]=0.0;
	//  Form1->Memo1->Lines->Add("Kernel dim(x+y)="+IntToStr(dim_x+dim_y));

	fftwf_execute(plan_kernel_fwd);

	old_alpha = alpha;
	old_arb_kernel = k_idx;
	Form1->Memo1->Lines->Add("DS connectivity calculations - Kernel calculated; alpha= "+FloatToStr(alpha));
}

void copy_occ_to_fftwf_LS(int sp, int add_mode, float mult)
{
	int x, y, loop;
	struct cell *c1;
	float *p, *v, sum, val;

	if (!add_mode)
	{
		for(y=0; y<dim_y; y++)
			for(x=0; x<dim_x; x++)
				orig_LS[y][x]=0.0;
	}

	sum=0.0;
	for(y=0; y<dim_y/2; y++)
		for(x=0; x<dim_x/2; x++)
		{
			if (vmat[y][x])
				val=vmat[y][x][sp]; //.matval(x,y);
			else
				val=0;
			if (val<0)
				val=0.0;
			if (add_mode)
				orig_LS[y][x] += mult*val;
			else
				orig_LS[y][x]  = val;

			sum +=val;
		}

	//  Form1->Memo1->Lines->Add("LSSum = "+FloatToStr(sum));
}


void copy_matrix_to_fftwf_LS(int add_mode, float **m)
{
	int x, y;

	if (!add_mode)
	{
		for(y=0; y<dim_y; y++)
			for(x=0; x<dim_x; x++)
				orig_LS[y][x]=0.0;
	}

	for(y=0; y<dim_y/2; y++)
		for(x=0; x<dim_x/2; x++)
		{
			if (m[y][x]>0.0)
				orig_LS[y][x] = m[y][x];
			else
				orig_LS[y][x] = 0.0;
		}
}

void copy_matrix_to_fftwf_LS_edge(int add_mode, float **m)
{
	int x, y;

	if (!add_mode)
	{
		for(y=0; y<dim_y; y++)
			for(x=0; x<dim_x; x++)
				orig_LS[y][x]=0.0;
	}

	for(y=0; y<dim_y/2; y++)
		for(x=0; x<dim_x/2; x++)
		{
			if (m[y][x]>=0.0)
				orig_LS[y][x] = m[y][x];
			else
				orig_LS[y][x] = 0.0; // missing not as edge "water"
		}
}


bool alloc_fftwf_matrixes()
{
	char txt[255];

	dim_x = 2*xd;  // fix here according to application
	dim_y = 2*yd;
	LS_vec_len = dim_x*dim_y;
	fft_vec_len= max(dim_x,dim_y)*(min(dim_x,dim_y)/2+1);
	//  fft_vec_len= dim_x*(dim_y/2+10); // xxx there is something wrong here; should be +1?

	//  sprintf(txt, "Xd=%i yd=%i", dim_x, dim_y);
	//  Form1->Memo1->Lines->Add(txt);

	if (!dim_x || !dim_y)
		return false;

	orig_LS = matrix(0,dim_y-1,0,dim_x-1);
	if (!orig_LS)
		return false;

	conn_LS = matrix(0,dim_y-1,0,dim_x-1);
	if (!conn_LS)
		return false;

	kernel = matrix(0,dim_y-1,0,dim_x-1);
	if (!kernel)
		return false;

	LS_fft     = (fftwf_complex*)fftwf_malloc(fft_vec_len*sizeof(fftwf_complex));
	if (!LS_fft)
		return false;

	//  kernel_fft     = (fftwf_complex*)fftwf_malloc((dim_x*(dim_y/2+1))*sizeof(fftwf_complex));
	kernel_fft     = (fftwf_complex*)fftwf_malloc(fft_vec_len*sizeof(fftwf_complex));
	if (!kernel_fft)
		return false;

	for(int loop=0; loop<fft_vec_len;loop++)
	{
		LS_fft[loop][0]=0.0;
		LS_fft[loop][1]=0.0;
		kernel_fft[loop][0]=0.0;
		kernel_fft[loop][1]=0.0;
	}
	return true;
}

void free_fftwf_matrixes()
{
	if (orig_LS)
		free_matrix(orig_LS, 0, dim_y-1, 0, dim_x-1);
	orig_LS=0;

	if (conn_LS)
		free_matrix(conn_LS, 0, dim_y-1, 0, dim_x-1);
	conn_LS=0;

	if (kernel)
		free_matrix(kernel, 0, dim_y-1, 0, dim_x-1);
	kernel=0;

	if (LS_fft)
		fftwf_free(LS_fft);
	LS_fft=0;

	if (kernel_fft)
		fftwf_free(kernel_fft);
	kernel_fft=0;
}

void  convolve(int len, fftwf_complex *LS_fft, fftwf_complex *kernel_fft)
{
	float ar, br, ai, bi;

	for(int loop=0; loop<len; loop++)
	{
		ar = LS_fft[loop][0];
		ai = LS_fft[loop][1];
		br = kernel_fft[loop][0];
		bi = kernel_fft[loop][1];
		LS_fft[loop][0]=ar*br-ai*bi;
		LS_fft[loop][1]=ai*br+bi*ar;
	}
}

void generate_fftwf_plans()
{
#if !defined(__MINGW32__)
        fftwf_plan_with_nthreads(nthreads);
#endif
	Form1->Memo1->Lines->Add("Using " + IntToStr(nthreads) + " thread(s) for preprocessing.");
	plan_kernel_fwd=fftwf_plan_dft_r2c_2d(dim_y, dim_x, &kernel[0][0],
					      kernel_fft, FFTW_ESTIMATE);
	plan_LS_fwd=fftwf_plan_dft_r2c_2d(dim_y, dim_x, &orig_LS[0][0], LS_fft, FFTW_ESTIMATE);
	plan_LS_bwd=fftwf_plan_dft_c2r_2d(dim_y, dim_x, LS_fft, &conn_LS[0][0], FFTW_ESTIMATE);
}

void free_fftwf_plans()
{
  if (plan_kernel_fwd)
    fftwf_destroy_plan(plan_kernel_fwd);
  if (plan_LS_fwd)
    fftwf_destroy_plan(plan_LS_fwd);
  if (plan_LS_bwd)
  fftwf_destroy_plan(plan_LS_bwd);
}

void fftwf_connectivity()
{
	int   loop, x, y;
	float maxv;

	fftwf_execute(plan_LS_fwd);
	convolve(fft_vec_len, LS_fft, kernel_fft);
	fftwf_execute(plan_LS_bwd);

	maxv=-std::numeric_limits<float>::max();

	for(y=0; y<dim_y/2; y++)
		for(x=0; x<dim_x/2; x++)
		{
			conn_LS[y][x]/=LS_vec_len;
			maxv = max(maxv, conn_LS[y][x]);
		}

	if (IG_normalize)
		for(y=0; y<dim_y/2; y++)
			for(x=0; x<dim_x/2; x++)
				conn_LS[y][x] = conn_LS[y][x]/maxv;

	// fedemp: 20120301 - removed msg too obscure for users
	//Form1->Memo1->Lines->Add("Maximum (smoothed) fraction of the distribution in any cell = "+FloatToStr(maxv));
}


