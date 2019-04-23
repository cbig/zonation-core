#include "bat_run.h"
#include "Unit1.h"
#include "ADMUs.h"
#include "marginal_loss.h"
#include "GridMap.h"
#include "defines.h"
#include "randz.h"
#include "LoadData.h"
#include "output.h"
#include "post_process.h"
#include "PLULA.h"
#include "LSIdent.h"
#include "ccsp.h"
#include "grid_utils.h"
#include "config.h"
#include "zig4lib/pod.h"
#include "zig4lib/zconf.h"
#include <cstdio>
#include <ctime>
#include <cstring>
#include <limits>
#include <queue>
#include <boost/thread.hpp>
#include <limits>

// more memory, faster search when removing cells from the edge list (warp blocks)
#define USE_MAT_EDGE_LIST_POS 1

//#define COST_DIV_DOUBLE
#undef COST_DIV_DOUBLE
#define CAZ_ALTERNATING_EPS 1

int nan_count = 0;

struct sort_xy *srtvec;
int    xd, yd, curve_store_interval, ADM_store_interval, ADM_per_admu_store_interval;
double *repr_lvls, *repr_after_rem;
int   rem_by_rarity, rem_step_by_rarity, c_pos, use_edge_rem, use_smoothing;
// double because they work as accumulator vars
double cost_used=0.0d, cost_in_area=0.0d, cost_out=0.0d, 
  cost_rat, cost_in_init_rem;
float alpha_mult;
float max_repr = -std::numeric_limits<float>::max(), min_repr = std::numeric_limits<float>::max();
String name_addition;
bool  bat_mode=false;
int   annotate_name = 0, all_zero_as_miss=0;

int    **nbm;

int   warp_factor, warp_ll, sliding_initial_min_mask_level=std::numeric_limits<int>::max(), current_min_mask_level=std::numeric_limits<int>::max();

// warp list
std::vector<int> wrx, wry, wnbc;  // (x,y) coord and number of neighbors
std::vector<float> wval;   // delta val
// only for PLU-mode (special warp list)
std::vector<int> plu_copy_wrx, plu_copy_wry;

int resamp_vec[MAX_SPP_COUNT];

int ties_cnt = 0;
//float wdirty[1000000];
#include "matrix_utils.h"
int** wdirty = 0;
bool use_8_connectivity = false;
//float** wdirty = 0;

int warp_dirty_count = 0, warp_nondirty_count = 0, warp_recalculations = 0;
int   BL;

// variables for solution loading
int    run_mode;
class  GridMap loaded1;
String load1fn="";
struct sort_xy *load_order;
int    load_max, loaded_cnt, logit;
struct IG_settings IG_set;

// variables for 
Grid_CCSP* cc_accounting = NULL;
float CCSP_prev_removal_threshold = 0;
float CCSP_prev_prev_removal_threshold = 0;
float CCSP_ma_removal_rate = 1;
float CCSP_last_failed_removals = 0;
FILE* wlist_file;
const bool CCSP_write_log_file = false;

const bool CCSP_slow_as_fast = true;
bool corr_ongoing = false;

// change of data structure
// species data array (matrix of arrays)
#ifdef COMPACT_VMAT
 Biodiv_Features_Occur_Container ** vmat = 0;
 BFOC_size_t Biodiv_Features_Occur_Container::priv_it = 0;
 float Biodiv_Features_Occur_Container::priv_empty_missing = -1;
 BFOC_size_t* Biodiv_Features_Occur_Container::plu_idx_vec = NULL;
#else
 // classic Z
 // This will not work anymore...
 float*** vmat = 0;
#endif

#ifdef COMPACT_VMAT_LOOKUP
 // Nothing else required
#elif defined NONCOMPACT_VMAT_FLOATP
 // needs global size of 'rowp' vectors
 size_t Biodiv_Features_Occur_Container::priv_size = 0;
#endif

// BQP variables
struct BQP_profile BQPs[100];
int    radii_cnt, radii[1000], prof_cnt, radii_sp[1000];
int    **nbms[1000], nbm_num[MAX_SPP_COUNT];
float  **nbms_orig[1000];
float  smult[MAX_SPP_COUNT], BQP_spp_losses[MAX_SPP_COUNT], sp_derivative[MAX_SPP_COUNT];
float  BQP_spp_loss_tmp[MAX_SPP_COUNT]; //this is being edited inside threads!
int    BQP_mode=2;

// SSI variables
struct SSI_sp  SSI[MAX_SPP_COUNT];
struct SSI_xydata   *SSI_list=0;
int    SSI_spp_cnt=0, SSI_loc_cnt, SSI_row_cnt=0, SSI_lpos=0, SSI_err_cnt;
int    **SSIxyPos=0, **SSIxyCnt=0;
//int    *SSI_sp_num;     // fedemp: old SSI stuff commented out, apparently not used
//float  *SSI_occ_size;   
double *SSI_repr_lvls;
float  SSI_curves[CURVES_LEN][3];
float **SSI_indiv_curves;

// Corridors
struct Corr_settings corr_set;

// Connected component variables
String output_ccsp_fn;
// >=1: periodic msgs, >=2: debug x,y msgs, >=3: debug even more details (get_max_..., cached, etc.)
int CCSP_verbose = 0;
size_t CCSP_info_period;
int CCSP_formula, CCSP_variant, CCSP_thickness;

bool glob_set_output_wrscr_map = true;
bool glob_set_output_prop_rank_map = false;

// only used when correcting the relative contributions of different cells to the total landscape size
double w_landscape_remaining = 1.0d;
// Must be initialized to the sum of all the cell occurrence weights
double total_landscape_occur_size_weight = .0d;

// assumes non-missing cell coordinates
float get_cell_area_correction_factor(int x, int y)
{
  if (!use_cost) {   // costmap.m has cell weights
    if (costmap.m[y][x] >=0)
      return costmap.m[y][x];
    else
      return 0;
  } else {  // need to keep in memory the cell weights map, and costmap.m is different
    if (occur_size_weights_map.m[y][x] >=0)
      return occur_size_weights_map.m[y][x];
    else
      return 0;
  }
}

void
graceful_exit()
{ 
  fflush(stdout);
  fflush(stderr);
  free_data();
  exit(1);
}


inline void
add_into_edge_list(int x, int y)
{
  if(edge[y][x]>0)
    return;

#ifdef USE_MAT_EDGE_LIST_POS
  mat_edge_list_pos[y][x] = ecnt;
  //  Form1->Memo1->Lines->Add("=e");
#endif
  edge[y][x] = 1;
  exl[ecnt]  = x;
  eyl[ecnt]  = y;
  ecnt++;
}

void Mark_PLULA_to_edge(int pln)
{
	int loop, end, i, j;

	//	  if (pln<0)
	//	    {
	//	      ShowMessage("Fooken PLN error");
	//	      return;
	//	    }

	if (PLvec[pln].at_edge || PLvec[pln].removed)   // pln=-1111321321  xxxPLULA error here
		return;                                       // xxxErrorHere

	end = PLvec[pln].start+PLvec[pln].el_cnt;
	for(loop=PLvec[pln].start;loop<end; loop++) {
	  i = PLX[loop];
	  j = PLY[loop];
	  if ((status[j][i]>0) && (edge[j][i]!=1)) {
	    add_into_edge_list(i, j);
	  }
	  //        nbm[j][i]--; // ok for BQP, only concerns mark to edge
	  // xxxPLULA, unlikely to work correct
	}
	PLvec[pln].at_edge=true;
}

void mark(int i, int j)
{
	if ((status[j][i]>0) && (edge[j][i]!=1))
	{
		if (use_PLULA) {
		  if (PLL[j][i]>=0)
		    Mark_PLULA_to_edge(PLL[j][i]);
		  //          else
			//            Form1->Memo1->Lines->Add("PLL Error at j="+IntToStr(j)+"  i="+IntToStr(i));
		} else {
		  add_into_edge_list(i, j);
		}
	}
	nbm[j][i]--;  // xxx ok for BQP_nbm, mark called for neighbours
	// xxx what the hell is this; add_edge_points ok???
}

void
Mark_neighbors_to_edge(int x, int y)
{
  int  sx, sy, ex, ey;

  sx = max(0,   x-1);
  ex = min(x+1, xd-1);
  sy = max(0,   y-1);
  ey = min(y+1, yd-1);

  /* TODO: test this. Seems to reduce comp time to 85%. BUT: problem with edge layer in ccsp
     /* See "alt imp below" */
  /**/
  //bool use_trick_minimal_edge_list = true;
  /*
    if (use_trick_minimal_edge_list) {
    if (x>0) {
    if ((status[y][x-1]>0) && (edge[y][x-1]!=1)) {
    add_into_edge_list(x-1, y);
    }
    nbm[y][x-1]--; // ok for BQP, only concerns mark to edge	  
    }
    if (x<xd) {
    if ((status[y][x+1]>0) && (edge[y][x+1]!=1)) {
    add_into_edge_list(x+1, y);
    }
    nbm[y][x+1]--; // ok for BQP, only concerns mark to edge	  
    }
    if (y>0) {
    if ((status[y-1][x]>0) && (edge[y-1][x]!=1)) {
    add_into_edge_list(x, y-1);
    }
    nbm[y-1][x]--; // ok for BQP, only concerns mark to edge	  
    }
    if (y<yd) {
    if ((status[y+1][x]>0) && (edge[y+1][x]!=1)) {
    add_into_edge_list(x, y+1);
    }
    nbm[y+1][x]--; // ok for BQP, only concerns mark to edge	  
    }
    return;
    }
  */
  /**/
  
  for(int j=sy; j<=ey; ++j) {
    for(int i=sx; i<=ex; ++i) {
      /* TODO: "alt impl below" */
      /*
	use_trick_minimal_edge_list = true;
	if (use_trick_minimal_edge_list && 
	!( (i==x) || (j==y) ) ) {
	//edge_ccsp[j][i] = 1;   // edge mark needed for ccsp
	continue;
	}
      */
      if ((status[j][i]>0) && (edge[j][i]!=1)) {
	if (use_PLULA) {
	  if (PLL[j][i]>=0)
	    Mark_PLULA_to_edge(PLL[j][i]);
	} else {
	  add_into_edge_list(i, j);
	}
      }
      nbm[j][i]--; // ok for BQP, only concerns mark to edge
    }
  }
}

// Marks as "edge" coordinates x,y if any neighbor cell has a higher level in the removal mask
// * Assumes x>0, x<xd, y>0, y<yd
inline void
mark_removal_mask_borders_to_edge(int x, int y)
{
  float level = maskmap.m[y][x];

  int minx = max(0,   x-1);
  int maxx = min(x+1, xd-1);
  int miny = max(0,   y-1);
  int maxy = min(y+1, yd-1);

  for(int j=miny; j<=maxy; j++) {
    for(int i=minx; i<=maxx; i++) {
      if (maskmap.m[j][i] > level && status[j][i]>0 && 1!=edge[y][x]) {
	if (use_PLULA) {
	  if (PLL[j][i]>=0)
	    Mark_PLULA_to_edge(PLL[j][i]);
	} else {
	  add_into_edge_list(x, y);
	}
	return;
      }
    }
  }
}

void  get_boundary_length(int prints, int prop_mode)
{
	int  x, y;
	const size_t LINE_LEN = 512;
	char txt[LINE_LEN];

	BL=0;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (status[y][x]<=0)
				continue;

			if (x==0)
				BL++;
			else if (status[y][x-1]<=0)
				BL++;

			if (y==0)
				BL++;
			else if (status[y-1][x]<=0)
				BL++;

			if (x==(xd-1))
				BL++;
			else if (status[y][x+1]<=0)
				BL++;

			if (y==(yd-1))
				BL++;
			else if (status[y+1][x]<=0)
				BL++;
		}
	}

	if (prints) {
	  snprintf(txt, LINE_LEN, "  BL=%i  Area=%i BL/A=%0.3f", BL, nonm1-removed, BL/(float)(nonm1-removed));
	  Form1->Memo1->Lines->Add(txt);
	}
}

void
init_corridors()
{
  // It is critical to initialize this after Calc_richness_et_al_matrixes() (depends on Rmax, vmat, etc. if using wrscr).
  // Here in Initialize_remove() it is safe.
  if (CCSP <= .0f) 
    return;

  corr_ongoing = true;

  String alloc_msg = " *****=====----- Allocating memory for corridor building and initializing...";
  if (CCSP_verbose <= 0)
    output_ccsp_fn = "";
  try {
    cc_accounting = new Grid_CCSP(CCSP, CCSP_variant, CCSP_formula, CCSP_thickness, use_8_connectivity, output_ccsp_fn, xd, yd, status);
    Form1->Memo1->Lines->Add(alloc_msg + " done!");
  } catch(std::bad_alloc const& ex) {
    Form1->Memo1->Lines->Add(alloc_msg + " error: " + String(ex.what()));
  }
}

// Important: transfers the -1s in Rmax[][] to status[][]
// sets edge[][] to 0s - but proper initialization happens in Initialize_remove()
void initialize_status()
{
  for(int y=0; y<yd; y++) {
    for(int x=0; x<xd; x++) {
      edge[y][x]=0;
      if ((x==0) || (y==0) || (x==(xd-1)) || (y==(yd-1))) { // xxx enabled 19.9.04
	status[y][x]=-1;
	continue;
      }
      
      if (use_PLULA && (PLL[y][x]<0)) { // if PLULA rem and plu num not exist
	status[y][x]=-1;
	continue;
      }
      
      if (Rmax[y][x]==-1) {
	status[y][x]=-1;
	//cnt++;
      } else {
	status[y][x]=1;
      }
    }
  }  
}

// Fills in edge[][] (should have been initialized to 0s in initialize_status()
void Initialize_remove()
{
	int   x,y, cnt, rnds;

	removed= 0;
	current= 1;

	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("******=====----- Preparing to start the ranking process... -----=====**********");

#ifdef USE_MAT_EDGE_LIST_POS
	mat_edge_list_pos = imatrix(0, yd, 0, xd);
	for(int j=0;j<yd;j++) {
	  for(int i=0;i<xd;i++) {
	    mat_edge_list_pos[j][i] = -1;
	  }
	}
#endif

	//  Form1->Memo1->Lines->Add("Cells with missing data count ="+IntToStr(cnt));
	if (!use_edge_rem) {
	  Form1->Memo1->Lines->Add("NOT using removal only from edges.");
	} else {
	  Form1->Memo1->Lines->Add("Note: using edge removal.");
	  if (add_edge_cnt>0)
	    Form1->Memo1->Lines->Add("   Count of additional edge points = "+IntToStr(add_edge_cnt));
	}

	if (!use_edge_rem) {
	  for(y=0; y<yd; y++)
	    for(x=0; x<xd; x++)
	      mark(x,y);
	} else {
	  for(y=0; y<yd; y++) {
	    for(x=0; x<xd; x++) {
	      if (status[y][x]==-1)
		Mark_neighbors_to_edge(x,y); // fixes also nmb counts
	    }
	  }
	  // Add borders between levels in the (hierarchical) removal mask
	  if (use_mask && glob_add_to_edge_borders_between_removal_mask_levels) {
	    Form1->Memo1->Lines->Add("Note: adding to edge list cells on the borders between different removal mask levels.");
	    for(y=1; y<yd-1; y++) {
	      for(x=1; x<xd-1; x++) {
		if (status[y][x] > 0)
		  mark_removal_mask_borders_to_edge(x,y);
	      }
	    }
	  }
	}

	if (BLP>0.0f)
		get_boundary_length(1, 0);

	if (100.0f <= corr_set.start_pc)
	  init_corridors();

	rnds=0;
	if ((add_edge_cnt>0) && use_edge_rem)
	{
		Form1->Memo1->Lines->Add("*********** Adding fake edge points; count = "+IntToStr(add_edge_cnt));
		Form1->Memo1->Lines->Add("");

		//cnt =0;
		do
		{
			x = static_cast<int>(xd*randz());
			y = static_cast<int>(yd*randz());
			if ((status[y][x]>0) && (1!=edge[y][x])) {
			  add_into_edge_list(x, y);
			}
			++rnds;
		}
		while ((cnt<add_edge_cnt) && (rnds<10*add_edge_cnt));
	}

	if (use_occur_size_weights_correct_landscape_fraction) {
	  total_landscape_occur_size_weight = .0d;
	  for (size_t y=0; y<yd; y++) {
	    for (size_t x=0; x<xd; x++) {  
	      if (status[y][x] > 0)  // count only effective cells! (after area mask, etc.)
		total_landscape_occur_size_weight += get_cell_area_correction_factor(x, y);
	    }
	  }
	  w_landscape_remaining = total_landscape_occur_size_weight;
	}
}

void Calc_richness_et_al_matrixes()
{
	int   x,y, s;
	float val;
	int   nok_cnt;

	nok_cnt=0;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			//          Rmin[y][x] = FLT_MAX;
			Rmax[y][x] =-std::numeric_limits<float>::max();
			if (!forget_Rsum && !mem_save_mode)
				Rsum[y][x]=0.0f; // Rave
		}
	}

	for(s=0;s<map_cnt;s++)
		spp[s].max_local_val = 0.0f;

	Rmax_max=Rsum_max=0.0f; //Rmin_max, Rave_max

	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (!vmat[y][x])
			{
				Rmax[y][x]=-1; // Rmin, Rave
				continue;
			}

			//for(s=0;s<map_cnt;s++)
			const Biodiv_Features_Occur_Container& rowp = vmat[y][x];
			for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s)) {
				val = vmat[y][x][s];
				if (val!=-1)
				{
					if (val>spp[s].max_local_val)
						spp[s].max_local_val = val;

					val *= spp[s].weight; // xxxW
					if (use_cost)
					{
						if (costmap.m[y][x]>0.0f)
							val /= costmap.m[y][x];
						else
						{
							costmap.m[y][x]=1.0f;
							++nok_cnt;
						}
					}
				}

				if (spp[s].weight<=0.0f)
					continue;  // xxxMCZ - ok, I guess 5/09

				//              Rmin[y][x]=min(Rmin[y][x], val);
				Rmax[y][x]=max(Rmax[y][x], val);
				if (val!=-1)
				{
					//                  Rave[y][x] += val;
					if (!forget_Rsum && !mem_save_mode)
						Rsum[y][x] += val;
				}
			}

		}
	}

	if (all_zero_as_miss)
	{
		for(y=0; y<yd; y++)
			for(x=0; x<xd; x++)
			{
				if (!vmat[y][x])
					continue;
				if (Rmax[y][x]<=0.0f)
				{
					Rmax[y][x]=-1;
					// free(vmat[y][x]); // COMPACT_VMAT
					// vmat[y][x]=0;
					vmat[y][x].clear();
				}
			}
	}

	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (Rmax[y][x]==-1)
			{
				Rmax[y][x]=-1;
				if (!forget_Rsum && !mem_save_mode)
					Rsum[y][x]=-1; //Rmin,  Rave
			}
		}
	}

	nonm1 = m1s = 0;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			//          Rmin_max = max(Rmin_max, Rmin[y][x]);
			Rmax_max = max(Rmax_max, Rmax[y][x]);
			//          Rave_max = max(Rave_max, Rave[y][x]);
			if (!forget_Rsum && !mem_save_mode)
				Rsum_max = max(Rsum_max, Rsum[y][x]);
			if (Rmax[y][x]==-1)
				m1s++;
			else
			{
				nonm1++;
				if (use_BQP && (BQP_mode==1))
					for(s=0;s<map_cnt;s++)  // xxx fixes BQP problem, but deletes missing data
					{
						if (spp[s].BQP_used)
							if (vmat[y][x][s]<0.0f)
								vmat[y][x][s]=0.0f;
					}
			}
		}
	}

	// here we have the total count of removable elements (nonm1) ready
	//printf("NONM1 READY\n");
	exl = new int[nonm1+1];
	eyl = new int[nonm1+1];
	srtvec= new struct sort_xy[nonm1+1];

	if (nok_cnt>0)
	{
		Form1->Memo1->Lines->Add("*****************************************************");
		Form1->Memo1->Lines->Add("ERROR: missing cost value where biological data exists; count of locations: "+IntToStr(nok_cnt));
		Form1->Memo1->Lines->Add("These have been replaced with cost=1.0");
		Form1->Memo1->Lines->Add("*****************************************************");
	}
}

// if use_hier_mask, a special sortfunc enforces the "hierarchical removal mask" (map in maskmap.m)
int  Get_sorted_cells_in_map(struct sort_xy *vec, float **R, bool use_hier_mask=false)
{
	int pos;

	pos=0;
	for(int y=0; y<yd; y++)
	{
		for(int x=0; x<xd; x++)
		{
			if (R[y][x]!=-1)
			{
				vec[pos].x  = x;
				vec[pos].y  = y;
				vec[pos].val= R[y][x];
				pos++;
			}
		}
	}

	if (!use_hier_mask)
	  qsort((void *)vec, nonm1, sizeof(struct sort_xy), sortfunc);
	else
	  qsort((void *)vec, nonm1, sizeof(struct sort_xy), sortfunc_with_hierarchy);

	return nonm1;
}

void Update_nbms(int x, int y)
{
	int  mat, rad, i, j, spnum;
	int  sx, sy, ex, ey, *row;

	// xxx BQP
	for(mat=0;mat<radii_cnt;mat++)
	{
		rad = radii[mat];
		if (BQP_mode==2) // BQPm2
		{
			spnum = radii_sp[mat];
			if (vmat[y][x][spnum]<0.0f)
				continue;
		}

		sx=max(0, x-rad);
		ex=min(x+rad, xd-1);
		sy=max(0, y-rad);
		ey=min(y+rad, yd-1);

		for(j=sy; j<=ey; ++j)
		{
			row = &nbms[mat][j][0];
			for(i=sx; i<=ex; ++i)
			{
				//              if (row[i]>0)
				--row[i]; // xxx BQP could check - should not get negative  even without if
			}
		}
	}
}

void Remove_site(int x, int y, int lpos)
{
	int   BL_nbc, pos;
	const size_t LINE_LEN = 512;
	char txt[LINE_LEN];

	//if (Rmax[y][x]==-1)
	if (-1 == status[y][x])	{
	  snprintf(txt, LINE_LEN, "Trying to remove x %i y %i lp=%i with no features present or already removed!", x, y, lpos);
	  Form1->Memo1->Lines->Add(txt);
	}
	if (use_BQP)
		get_delta_value(x,y, true);

	status[y][x] = 0;
	edge[y][x]   = 0;
	BL_nbc       = 0;
	//  Form1->Memo1->Lines->Add(FloatToStr(BLP));
	if (BLP>0.0f)
	{
                // this is get_nb_count(status, x, y)!
		if (status[y][x-1]>0)
			BL_nbc++;
		if (status[y][x+1]>0)
			BL_nbc++;
		if (status[y-1][x]>0)
			BL_nbc++;
		if (status[y+1][x]>0)
			BL_nbc++;

		switch(BL_nbc)
		{
		case 0: BL -= 4; break;
		case 1: BL -= 2; break;
		case 2: BL += 0; break;
		case 3: BL += 2; break;
		case 4: BL += 4; break;
		}
	}

	/*
	// update blob specific BL
	if (CCSP > .0f) {
	  // TODO: seems not to work!!!
	  //int bl_incs[] = {-4, -2, 0, 2, 4};
	  //cc_accounting->update_bl(x, y, bl_incs[BL_nbc]);
	}
	*/

	//  if (use_PLULA) // done in fillwarplist
	//    {
	//      int tmpnum;
	//      tmpnum = PLL[y][x];
	//      if (tmpnum>0)
	//        PLvec[tmpnum].removed = true;
	//    }

	Mark_neighbors_to_edge(x,y);

	if (use_BQP)
		Update_nbms(x,y); // updates nbhood cell counts for different buffers

	if (lpos!=-1) {
#ifdef USE_MAT_EDGE_LIST_POS
	  mat_edge_list_pos[y][x] = -1;
	  mat_edge_list_pos[eyl[ecnt-1]][exl[ecnt-1]] = lpos;
#endif
	  exl[lpos] = exl[ecnt-1];
	  eyl[lpos] = eyl[ecnt-1];
	  ecnt--;

	}

	removed++;
	if (use_cost) {
		cost_used += costmap.m[y][x];

		if(ADM_set.use_ADMUs && 
		   (use_groups || ADM_set.row_count_for_per_admu_curves > 0) ) {
		  int seq = ADMUs[y][x];
		  if (seq>0 && seq<ADM_set.count)
		    ADMU_cost_used[seq] += costmap.m[y][x];
		}
	} else {
	  if(ADM_set.use_ADMUs && 
	     (use_groups || ADM_set.row_count_for_per_admu_curves > 0) ) {
	    int seq = ADMUs[y][x];
	    if (seq>0 && seq<ADM_set.count)
	      ADMU_cost_used[seq]++;
	  }
	}

	//	  if (vmat[y][x]<=0)      // xxxPLULAError
	//	  {
	//	    Form1->Memo1->Lines->Add("working with location with missing data x="+IntToStr(x)+" y="+IntToStr(y));
	//		obsmap[0].set_color_for_pixel(Form1->Image1->Picture->Bitmap,
	//			x,y,1.0, sol[y][x], true);
	//		return;
	//	  }

	// rowp = &vmat[y][x][0];  // COMPACT_VMAT
	const Biodiv_Features_Occur_Container& rowp = vmat[y][x];
	if (!use_tree_conn && rowp) // if tree conn, reprs have been updated in fill warp list
	        // for(s=0; s<map_cnt; s++)
		for(int s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))
		{
			float rowp_s_val = rowp[s];
			if (rowp_s_val<0)  // BQPm2 changed
				continue;

			if (!use_BQP || (!spp[s].BQP_used))
			{
			        // repr_lvls[s] -= rowp[s];
			        repr_lvls[s] -= rowp_s_val; //rowp[s];
			}
			else
			{
				repr_lvls[s] -= BQP_spp_loss_tmp[s];
			}
		}

	if (removal_rule==3) // TBF
	{
		//    Form1->Memo1->Lines->Add("here ");
		for(int s=0; s<map_cnt; s++)
		{
			//        Form1->Memo1->Lines->Add("here "+FloatToStr(spp[s].sr_par)+"  "+FloatToStr(repr_lvls[s]) );
			if ( (spp[s].T_violation_fract<=-1.0f) &&  (repr_lvls[s]<spp[s].sr_par))
			{
				spp[s].T_violation_fract=removed/(float)(nonm1*(1.0f-rem_level));
				//          Form1->Memo1->Lines->Add("Violated");
			}
		}
	}

	if (use_SSI)
	{
		for(int s=0;s<SSIxyCnt[y][x];s++)
		{
			pos = SSIxyPos[y][x]+s;
			SSI_repr_lvls[SSI_list[pos].spnum] -= SSI_list[pos].val;  // xxx SSI
			if (SSI_repr_lvls[SSI_list[pos].spnum]<=0.0f && 
			    SSI[SSI_list[pos].spnum].lost_at_fraction==-1.0f)
			{
			  SSI[SSI_list[pos].spnum].lost_at_fraction = removed/(float)(nonm1*(1.0f-rem_level));
			}
			if (removal_rule==3) // TBF
			{
				if (SSI_repr_lvls[SSI_list[pos].spnum]<SSI[SSI_list[pos].spnum].sr_par &&
				    SSI[SSI_list[pos].spnum].T_violation_fract==-1.0f)
				{
				  SSI[SSI_list[pos].spnum].T_violation_fract = removed/(float)(nonm1*(1.0f-rem_level));
				}
			}
		}
	}

	// ADMxxx: if tree conn, reprs should be updated in fill warp list
	if (ADM_set.use_ADMUs && !use_tree_conn && rowp) {
	  int adu = ADMUs[y][x];
	  if (adu<0 || adu>=ADM_set.count) {
	    adu=0;
	  } else {
	    //for(s=0; s<map_cnt; s++) {
	    for(int si = rowp.first(); si != rowp.overflow(); si = rowp.next(si)) {
	      if (!use_BQP || (!spp[si].BQP_used)) {
		float rowp_s_val = rowp[si];
		if (rowp_s_val>0.0f)
		  ADMUxSP_repr[adu][si] -= rowp[si];
	      } else {
		ADMUxSP_repr[adu][si] -= BQP_spp_loss_tmp[si];
	      }
	    if (ADMUxSP_repr[adu][si] < 0)
	      ADMUxSP_repr[adu][si] = 0.0;
	    }
	  }
	}

	// Update GUI rank map and progress status
	obsmap[0].set_color_for_pixel(Form1->Image1->Picture->Bitmap,
				      x,y,1.0f, sol[y][x], true);

	// Will not be needed anymore. Free as soon as possible!!!
	// Beware: post-proc needs it!
	if (PPA_fname.isEmpty()) {
	  vmat[y][x].clear();
	}
}

// Prints: cells_removed, prop_landscape_remaining, edge_list_len, smallest_remaining_prop, time_seconds
void
print_rank_progress_line(size_t current_cell, float prop_land_remaining, float cost_remaining, size_t edge_list_len, float min_prop, time_t elapsed)
{
  String s; // only safe option for the %
  s.sprintf("%14lu %14g%s %14g %14d %14g %12lu", current_cell, 100.0f*prop_land_remaining, "%%", cost_remaining, edge_list_len, min_prop, elapsed);
  Form1->Memo1->Lines->Add(s);
  /*
  String cost_s;
  cost_s.sprintf("%14g", cost_remaining);
  Form1->Memo1->Lines->Add(IntToStrW(current_cell, 14)+" "+ 
			   FloatToStrF(100.0f*prop_land_remaining, ffFixed, 5, 12)+"%  "+ 
			   //cost_s + " " +
			   IntToStrW(edge_list_len, 14)+" "+
			   FloatToStrF(min_prop, ffFixed, 7, 20)+" "+
			   IntToStrW(elapsed, 12));
*/
}

bool time_to_print_info_line(int current_removed)
{
  int rem_period_print = 10000;
  // this is to print at least a few lines so you don't run out of patience
  const int min_lines = 25;

  if (nonm1 < min_lines*rem_period_print)
    rem_period_print = nonm1/min_lines;

  return (0==current_removed%rem_period_print);
}

void
Remove_site_and_record_curves(int rx, int ry, int rpos, float val, int rnd, time_t t_ini)
{
  float ave_rl, min_rl;
  int s, tmp_map_cnt, tmp_SSI_cnt;

  if (!use_occur_size_weights_correct_ranking)
    sol[ry][rx] = current/float(nonm1*(1.0f-rem_level)); // cell ranking
  else
    sol[ry][rx] = 1.0d - w_landscape_remaining/total_landscape_occur_size_weight; // ranking corrected by area

  sol_val[ry][rx] = 1.0f-val; // xxxMCZ error

  Remove_site(rx, ry, rpos);
  //      obsmap[0].set_color_for_pixel(Form1->Image1->Picture->Bitmap,rx,ry,
  //        1.0f, sol[ry][rx]);
  current++;

  float prop_lost;
  if (!use_occur_size_weights_correct_landscape_fraction) {
    prop_lost = rem_level + float(current)/float(nonm1+1.0f);
  } else {
    prop_lost = 1.0d - w_landscape_remaining/total_landscape_occur_size_weight;
  }
  float cost_remain;
  if (!use_cost) {
    cost_remain = static_cast<float>(nonm1-current+1);
  } else {
    cost_remain = static_cast<float>(cost_in_area-cost_used);
  }

  if (time_to_print_info_line(current)) {
    //              Form1->Memo1->Lines->Add("chk BL="+IntToStr(BL));
    time_t elapsed = time(NULL) - t_ini;
    print_rank_progress_line(current, 1.0f-prop_lost, cost_remain, ecnt, val, elapsed);
    if (BLP>0.0f)
      get_boundary_length(1,0);
  }

  if (use_occur_size_weights_correct_landscape_fraction) {
    w_landscape_remaining -= get_cell_area_correction_factor(rx, ry);
  }

  if ((rnd==0) || ((rnd%curve_store_interval)==0) || (rnd==(nonm1-1))) {
    curves[c_pos][0] = prop_lost;

#if 0
    if (!((curves[c_pos][0]>=0.0f) && (curves[c_pos][0]<=1.0f))) {
      ShowMessage("c0or "+FloatToStr(curves[c_pos][0])+" current="+IntToStr(current)+"  rl="+FloatToStr(rem_level)+"  nm="+IntToStr(nonm1)+" calc=");
      Form1->Memo1->Lines->Add("NAN error");
      Form1->Memo1->Lines->Add("curves0 ="+FloatToStr(curves[c_pos][0]));
      Form1->Memo1->Lines->Add("ftmp ="+FloatToStr(prop_lost));
      Form1->Memo1->Lines->Add("curent ="+IntToStr(current));
      Form1->Memo1->Lines->Add("nonm1 ="+IntToStr(nonm1));
      Form1->Memo1->Lines->Add("cpos ="+IntToStr(c_pos));
      Form1->Memo1->Lines->Add("rem level ="+FloatToStr(rem_level));
      Form1->Memo1->Lines->Add("divi ="+FloatToStr((float)((nonm1+1.0f)*(1.0f-rem_level))));
      Form1->Memo1->Lines->Add("num ="+FloatToStr((1.0f-rem_level)*(float)current));
      Form1->Memo1->Lines->Add("total ="+FloatToStr(rem_level+(1.0f-rem_level)*(float)current/(((float)nonm1+1.0f)*(1.0f-rem_level))));
    }
#endif

    curves[c_pos][1] = cost_remain;
				
    // per ADMU costs (only required for certain outputs
    if(ADM_set.use_ADMUs && 
       (use_groups || ADM_set.row_count_for_per_admu_curves > 0) ) {
      for(size_t adu=0; adu < ADM_set.count; adu++)
	ADMU_cost[adu] = ADMU_cost_in_area[adu] -
	  ADMU_cost_used[adu];
    }

    ave_rl = 0.0f;
    min_rl = std::numeric_limits<float>::infinity();
    tmp_map_cnt=0;
    float weighted_ave_rl = 0.0f;
    float tmp_weight_agg = 0.0f;
    for(s=0; s<map_cnt; s++) {
      const int FIRST_SPP_SHIFT = 5;
      curves[c_pos][s+FIRST_SPP_SHIFT]=(float)(repr_lvls[s]);
      // repr_lvls NaNs!!!
      if (spp[s].weight>=0.0f)
	{
	  ave_rl += (float)(repr_lvls[s]);
	  weighted_ave_rl += spp[s].weight * (float)repr_lvls[s];
	  if (repr_lvls[s]<min_rl)
	    min_rl = repr_lvls[s];
	  tmp_map_cnt++;
	  tmp_weight_agg += spp[s].weight;
	}
    }
    ave_rl /= tmp_map_cnt;
    if (0.0f!=tmp_weight_agg)
      weighted_ave_rl /= tmp_weight_agg;
    curves[c_pos][2]=min_rl;
    curves[c_pos][3]=ave_rl;
    curves[c_pos][4]=weighted_ave_rl;
    Form1->Series2->AddXY(1.0-curves[c_pos][0],     min_rl,   "", clRed);
    Form1->LineSeries2->AddXY(curves[c_pos][1], min_rl,   "", clRed);
    Form1->Series3->AddXY(1.0-curves[c_pos][0],     ave_rl,"", clBlue);
    Form1->Series4->AddXY(curves[c_pos][1],     ave_rl,"", clBlue);
    Form1->Series_23->AddXY(1.0-curves[c_pos][0],     weighted_ave_rl,   "", clBlack);
    Form1->LineSeries2_Series4->AddXY(curves[c_pos][1],     weighted_ave_rl,   "", clBlack);
    if (ave_rl>=0.0f)
      Form1->LineSeries3->AddXY(1.0-curves[c_pos][0], 1.0f-z_pow(ave_rl, SA_z), "", clBlue);
    else
      Form1->LineSeries3->AddXY(1.0-curves[c_pos][0], 1.0f, "", clBlue);

    if (use_SSI) {
      ave_rl = 0.0f;
      min_rl = 10.0f;
      tmp_SSI_cnt=0;
      weighted_ave_rl = 0.0f;
      tmp_weight_agg = 0.0f;
      for(s=0; s<SSI_spp_cnt; s++) {
	SSI_indiv_curves[c_pos][s]=(float)(SSI_repr_lvls[s]);

	if (SSI[s].weight<0.0f)
	  continue;
	ave_rl += (float)(SSI_repr_lvls[s]);
	weighted_ave_rl += SSI[s].weight * (float)(SSI_repr_lvls[s]);
	if (SSI_repr_lvls[s]<min_rl)
	  min_rl = SSI_repr_lvls[s];
	tmp_SSI_cnt++;
	tmp_weight_agg += SSI[s].weight;
      }
      ave_rl /= tmp_SSI_cnt;
      if (0.0f!=tmp_weight_agg)
	weighted_ave_rl /= tmp_weight_agg;
      SSI_curves[c_pos][0]=min_rl;
      SSI_curves[c_pos][1]=ave_rl;
      SSI_curves[c_pos][2]=weighted_ave_rl;
      Form1->LineSeries5->AddXY(1.0-curves[c_pos][0], SSI_curves[c_pos][0], "", clRed);
      Form1->LineSeries6->AddXY(1.0-curves[c_pos][0], SSI_curves[c_pos][1], "", clBlue);
      Form1->LineSeries_56->AddXY(1.0-curves[c_pos][0], SSI_curves[c_pos][2], "", clBlack);
    }

    c_pos++;
  }

  if (ADM_set.row_count_for_per_admu_curves > 0 && 
      (rnd%ADM_per_admu_store_interval)==0 ) {
    ADMUs_output_per_admu_curves_iter(prop_lost);
  }
  if (ADM_set.use_ADMUs && use_groups && (rnd%curve_store_interval)==0)
    ADMUs_output_per_admu_grp_curves_iter(prop_lost, glbl_groups_info);
  if (CCSP > .0f && corr_ongoing && (0 == rnd % CCSP_info_period)) {
    cc_accounting->output_info_line(prop_lost);
    time_t elapsed = time(NULL) - t_ini;
    if (CCSP > .0f && CCSP_verbose > 0)
      Form1->Memo1->Lines->Add("CCSP count: " + IntToStr(cc_accounting->n_penalty_calc) + " " +
			       IntToStr(elapsed) + " s" + "(CC calc:" + IntToStr(cc_accounting->calc_count()) + ")" +
			       " removed: " + IntToStr(current - 1));
  }

}

int Check_PLXY_vecs()
{
	int loop, cnt, pln, start, end;

	cnt=0;
	for(loop=0; loop<nme_cnt; loop++)
	{
		pln = PLL[PLY[loop]][PLX[loop]];
		if (pln<0)
			cnt++;
		else
		{
			start = PLvec[pln].start;
			end =  PLvec[pln].start+ PLvec[pln].el_cnt;
			if ((loop<start) || (loop>=end))
				cnt++;
		}
	}

	if (cnt>0)
		ShowMessage("Cnt ="+IntToStr(cnt));
	return cnt;
}

// Only for PLU
int fill_warp_list_plu(int plu_n, int pos=0)
{
	int pln, loop, loop2, w2, s;


	if ((plu_n<0) ||  (plu_n>PLcnt))
	{
		ShowMessage("ERROR: Unknown pl# at spp data location.");
		return 0;
	}
	else
	{
		//      Form1->Memo1->Lines->Add("Removed id= "+IntToStr(plu_n));
		if (PLvec[plu_n].removed)
		{
			Form1->Memo1->Lines->Add("   !!LIKELY ERROR: Trying to removing for 2nd time, plu = "+IntToStr(PLvec[plu_n].num));
			return 0;
		}
	}

	// xxx lpos is problem - truly!!!
	int start = PLvec[plu_n].start;
	int end = start+PLvec[plu_n].el_cnt;
	//pos = 0;
	for(loop=start;loop<end; loop++)
	{
		if (pos >= wrx.size()) {
		  if (use_PLULA) {
		    Form1->Memo1->Lines->Add("**** ERROR: maximum size for PLU exceeded ("
					     +IntToStr(wrx.size())+" grid cells). Planning unit number "+
					     IntToStr(plu_n)+" (id: "+IntToStr(PLvec[plu_n].num)+") has "+IntToStr(PLvec[plu_n].el_cnt)+
					     " elements (cells). min x: "+IntToStr(PLvec[plu_n].minx)+
					     ", min y: "+IntToStr(PLvec[plu_n].miny)+
					     ", max x: "+IntToStr(PLvec[plu_n].maxx)+
					     ", max y: "+IntToStr(PLvec[plu_n].maxy)+
					     " Anything can happen after this point ****");
		  } else {
		    Form1->Memo1->Lines->Add("**** ERROR: inconsistency in warp list, maximum size PLU exceeded. Anything can happen after this point ****");
		  }
		}

		int x = PLX[loop];
		int y = PLY[loop];
		if ((x<0)||(y<0))
			ShowMessage("coordinates out of range");
		if (PLL[PLY[loop]][PLX[loop]]!=plu_n)
			ShowMessage("plu number mismatch error");

		//        if ((x==597) && (y==510)) // xxxP3, clean these when all works ookie
		//			Form1->Memo1->Lines->Add("Doing 597 510 ");

		if (/*(Rmax[y][x] == -1) || */ (status[y][x]<=0)) // cells in PLU with no spp data
		{
			// status[y][x] = 0; // mark out just to be sure // xxxP3 xommented out
			// edge[y][x]   = 0;   // xxxP3, no apparent influence
			Form1->Memo1->Lines->Add("PLU coordinate has no spp data, skipping, PLU ID = "+IntToStr(PLvec[plu_n].num));
			//            if ((x==597) && (y==510))
			//	    		Form1->Memo1->Lines->Add("skipped 597 510 ");
			continue; // xxxPLULAERROR pot here. Could it be problem if first cell in plula evaluated is not data?
		}

		wrx[pos] = x;
		wry[pos] = y;
		pos++;

	}

	PLvec[plu_n].removed = true;
	if (use_PLULA && use_tree_conn)
	{
		for(s=0; s<map_cnt; s++)
			repr_lvls[s] -= PLULA_tree_loss(plu_n, s); // xxxP too much is removed due multiple removal of repr!
		traverse_tree_down(plu_n, 2, true, false, 0);
		traverse_tree_up(plu_n,   2, true, false, 0);
	}

	//  Form1->Memo1->Lines->Add("   Count of cells in pln"+IntToStr(PLvec[pln].el_cnt));
	//  Form1->Memo1->Lines->Add("   Count inserted "+IntToStr(pos));

	//  loop2=Check_PLXY_vecs();
	//  if (loop2>0)
	//     Form1->Memo1->Lines->Add("   PLXY vectors corrupt in #locations= "+IntToStr(loop2));

	return pos;
}

void copy_plu_warp_list(int cell_warp_ll)
{
  plu_copy_wrx.resize(cell_warp_ll);
  plu_copy_wry.resize(cell_warp_ll);

  std::copy(wrx.begin(), wrx.begin()+cell_warp_ll, plu_copy_wrx.begin());
  std::copy(wry.begin(), wry.begin()+cell_warp_ll, plu_copy_wry.begin());
}

int fill_as_much_as_poss_warp_list_plu(int warp_plu_fill_cnt, int plu_n, int cell_warp_ll)
{
  int last_plu_n = plu_n;
  // 0==cell_warp_ll can happen at the end of the removal process

  // plu_warp_pos will iterate from 0 through plu_copy_wrx.size() 
  // (potentially larger range than warp_plu_fill_cnt < cell_warp_ll
  int plu_warp_pos = warp_plu_fill_cnt;

  while (cell_warp_ll>0 && warp_plu_fill_cnt < cell_warp_ll && plu_warp_pos < plu_copy_wrx.size()) {
    int next_x = plu_copy_wrx[plu_warp_pos];
    int next_y = plu_copy_wry[plu_warp_pos];
    if (next_x <0 || next_y <0) {
      break;
    }
    
    int next_plu_n = PLL[next_y][next_x];
    if (status[next_y][next_x] <= 0)
      Form1->Memo1->Lines->Add(" ERROR, x: "+IntToStr(next_x)+", y: "+IntToStr(next_x)+
			       ", status: "+IntToStr(status[next_y][next_x]));
    
    if (next_plu_n < 0 || next_plu_n == last_plu_n ||
	PLvec[next_plu_n].removed) {
      /*
	Form1->Memo1->Lines->Add(" * break in PLU warp (<), prev: "+IntToStr(last_plu_n) +" ID: "+IntToStr(PLvec[last_plu_n].num)+
	", next: "+IntToStr(next_plu_n)+" ID: "+IntToStr(PLvec[next_plu_n].num) + ", REMOVED: " +
	IntToStr(PLvec[next_plu_n].removed));
	break;
      */
      plu_warp_pos++;
      continue;
    }
    
    if (PLvec[next_plu_n].el_cnt > (cell_warp_ll-warp_plu_fill_cnt)) {
      /*
	Form1->Memo1->Lines->Add(" * break in PLU warp (>), cell_warp_ll: "+IntToStr(cell_warp_ll) +
	", next el_cnt: " +IntToStr(PLvec[next_plu_n].el_cnt)+
	", next plu seq: "+IntToStr(next_plu_n)+" ID: "+IntToStr(PLvec[next_plu_n].num));
	break;
      */
      plu_warp_pos++;
      continue;
    }

    Form1->Memo1->Lines->Add(" Removing PLU, id: "+IntToStr(PLvec[next_plu_n].num)+", cells: "+
			     IntToStr(PLvec[next_plu_n].el_cnt)+", x: "+IntToStr(next_x)+
			     ", y: "+IntToStr(next_y));
    
    int next_pos = fill_warp_list_plu(next_plu_n, warp_plu_fill_cnt);
    if (next_pos>0) {
      warp_plu_fill_cnt += next_pos - warp_plu_fill_cnt;
    }
    last_plu_n = next_plu_n;

    plu_warp_pos++;
    
    // avoid inf loop if there are inconsistencies
    if (0 == warp_plu_fill_cnt)
      break;
  }
  return warp_plu_fill_cnt;
}


inline void
edge_list_remove(int x, int y)
{
#ifdef USE_MAT_EDGE_LIST_POS
  int pos = mat_edge_list_pos[y][x];
  exl[pos] = exl[ecnt-1];
  eyl[pos] = eyl[ecnt-1];
  ecnt--;

  mat_edge_list_pos[y][x] = -1;
  mat_edge_list_pos[eyl[pos]][exl[pos]] = pos;
#else
  for(int pos=0; pos<ecnt; pos++) {
    if ((exl[pos]==x) && (eyl[pos]==y)) {
      exl[pos] = exl[ecnt-1];
      eyl[pos] = eyl[ecnt-1];
      ecnt--;
      break;
    }
  }
#endif
}

void compact_list(int x, int y)
{
  edge_list_remove(x, y);
}

// Rmax -> sol[][] / sol_val[][]
void
init_sol_and_sol_val()
{
  for(int y=0; y<yd; y++) {
    for(int x=0; x<xd; x++) {
      if (Rmax[y][x] == -1) {
	sol[y][x]     = -1;
	sol_val[y][x] = -1;
      } else {
	sol[y][x]     = Rmax[y][x];
	sol_val[y][x] = Rmax[y][x];
      }
    }
  }
}

void Remove_bad(float &prop)
{
	int cnt, cur, x, y;

	cost_in_init_rem=0.0f;
	if (run_mode==1)
	{
		cnt = (int)(prop*nonm1);
		for(cur=0; cur<cnt; cur++)
		{
			x = srtvec[cur].x;
			y = srtvec[cur].y;
			sol[y][x]     = 0.0f;
			sol_val[y][x] = 0.0f;
			if (edge[y][x]==1)
				compact_list(x,y);
			Remove_site(x,y,-1);
			if (use_cost)
				cost_in_init_rem += costmap.m[y][x];
		}
		Form1->Memo1->Lines->Add("Initially removed cells count = "+IntToStr(cnt));
	}
	else  // loading sol
	{
		cnt=0;
		for(y=0; y<yd; y++)
		{
			for(x=0; x<xd; x++)
			{
				if (loaded1.m[y][x] == 0)
				{
					sol[y][x]      =  0.0f;
					sol_val[y][x]  =  0.0f;
					if (edge[y][x] == 1)
						compact_list(x,y);
					Remove_site(x,y, -1);
					cnt++;
					if (use_cost)
						cost_in_init_rem += costmap.m[y][x];
				}
			}
		}
		prop = cnt/(float)(nonm1);
		Form1->Memo1->Lines->Add("Loaded initially removed from file: count = "+IntToStr(cnt));
	}

	if (prop>0.0f)
	{
		Form1->Memo1->Lines->Add("Cost of initially removed cells = "+FloatToStr(cost_in_init_rem));
		Form1->Memo1->Lines->Add("***** Note: Using initial removal is not recommended. *****");
		Form1->Memo1->Lines->Add("*****       Use \"add edge points\" instead. *****");
	}
	else
		Form1->Memo1->Lines->Add("Note: No initial removal was used.");
}

int get_unsel_nbc(int x, int y)
{
	int sx, sy, ex, ey, cnt;

	sx=max(0, x-1);
	ex=min(x+1, xd-1);
	sy=max(0, y-1);
	ey=min(y+1, yd-1);
	cnt=0;

	for(int j=sy; j<=ey; j++)
		for(int i=sx; i<=ex; i++)
			if (status[j][i]<=0)
				++cnt;

	return cnt;
}

// ADMUxspp_loss: where per ADMUx spp repr losses are collected.
float BQP_dv_buf(int x, int y, int s)
{
	int   sx, sy, ex, ey, rad, adu;
	float dv, *dv_vec, *n_vec, *ovec, orig, fract_red, ratio;
	int   *rvec, remaining;

	if (!spp[s].BQP_used)
		return 0.0f;

	rad= spp[s].BQP_radius;
	sx = max(0, x-rad);
	ex = min(x+rad, xd-1);
	sy = max(0, y-rad);
	ey = min(y+rad, yd-1);

	n_vec  = spp[s].BQP_link;
	dv = 0.0f;
	for(int j=sy; j<=ey; j++)
	{
		ovec   = &nbms_orig[nbm_num[s]][j][0]; // orig # nbms
		rvec   = &nbms[nbm_num[s]][j][0]; // # nbs
		for(int i=sx; i<=ex; i++)
		{
			if ((i==x) && (j==y))
				continue;
			if (status[j][i]<=0)
				continue;
			if ((BQP_mode==2) && (vmat[j][i][s]<=0.0f))
				continue;
			orig      = ovec[i];
			remaining = rvec[i];
			//          ratio = 10000*remaining/orig;
			//              fract_red = n_vec[(int)(ratio)]-n_vec[(int)(ratio-10000/orig)];
			fract_red = n_vec[(int)(10000*remaining/orig)]-n_vec[(int)(10000*(remaining-1)/orig)];

			dv       += fract_red*vmat[j][i][s]; // slow indexing - outermost innermost

			if (ADM_set.use_ADMUs) // ADM addition 3/2010
			{
			  adu = ADMUs[y][x];
			  if (adu<0 || adu>=ADM_set.count) { // this may not really be needed
			    adu=0;
			  } else {
			    ADMUxspp_loss[adu][s] += dv;
			  }
			}

		}
	}

	return dv;  // returns the global BQP loss for spp s
}

//does not edit global variables
float PLULA_BQP_local_loss(int pln, int spnum)
{
	int   loop, end, x, y, s;
	float val, loss, remaining, orig;

	loss = 0.0f;
	end = PLvec[pln].start+PLvec[pln].el_cnt;
	for(loop=PLvec[pln].start;loop<end; loop++)
	{
		x = PLX[loop];
		y = PLY[loop];
		//if (Rmax[y][x]==-1)
		if (-1 == status[y][x])
		  continue;
		val = vmat[y][x][spnum];
		if (val<=0.0f)
		  continue;

		remaining = static_cast<float>(nbms[nbm_num[spnum]][y][x]); // number of nbs remaining
		orig      = nbms_orig[nbm_num[spnum]][y][x];
		if (orig>0)
			loss += val*spp[spnum].BQP_link[(int)(10000*remaining/orig)]; // matches initial #nbs and init prob
		else
			loss += val;
	}

	return loss;
}

//does not edit global variables
float PLULA_BQP_nbh_loss(int pln, int spnum) // xxx fix <0 vmat values also here
{
	int   loop, end, x, y, i, j, rad;
	float val, loss, remaining, orig;
	int   bsx, bex, bsy, bey;   // PL   BQP effects bounding box
	int   sx, ex, sy, ey;       // call BQP effects bounding box

	rad = spp[spnum].BQP_radius;
	bsx = max(0, PLvec[pln].minx-rad);
	bex = min(PLvec[pln].maxx+rad, xd-1);
	bsy = max(0, PLvec[pln].miny-rad);
	bey = min(PLvec[pln].maxy+rad, yd-1);

	for(j=bsy; j<=bey; j++)
	{
		for(i=bsx; i<=bex; i++)
			nwm[j][i] = 0;  // nwm is LSI matris, used in BQP as help during comp.
	}

	end  = PLvec[pln].start+PLvec[pln].el_cnt;
	for(loop=PLvec[pln].start;loop<end; loop++)
	{
		x = PLX[loop];
		y = PLY[loop];
		if (status[y][x]<=0)
			continue;
		if ((BQP_mode==2) && (vmat[y][x][spnum]<=0.0f))
			continue;
		sx = max(0,     x-rad);
		ex = min(x+rad, xd-1);
		sy = max(0,     y-rad);
		ey = min(y+rad, yd-1);
		for(j=sy; j<=ey; j++)  // loop box around xy and add loss
			for(i=sx; i<=ex; i++)
				nwm[j][i]++;
	}

	for(loop=PLvec[pln].start;loop<end; loop++) // null nb loss inside pln itself - all is lost there
	{
		x = PLX[loop];
		y = PLY[loop];
		nwm[y][x]=0;
	}

	loss = 0.0f;
	for(j=bsy; j<=bey; j++)  // loop box around PLN and tally loss
	{
		for(i=bsx; i<=bex; i++)
		{
			if (status[j][i]<=0) // xxx check if is necessary; idea is to skip missing and already removed locations
				continue;
			if (nwm[j][i]>0)
			{
				val = vmat[j][i][spnum];
				if (val<=0.0f)
					continue;

				remaining = static_cast<float>(nbms[nbm_num[spnum]][j][i]); // number of nbs remaining
				orig      = nbms_orig[nbm_num[spnum]][j][i];
				if (orig>0)
					loss += val*(  spp[spnum].BQP_link[(int)(10000*remaining/orig)]
						       - spp[spnum].BQP_link[(int)(10000*(remaining-nwm[j][i])/orig)] ); // matches initial #nbs and init prob
				else
				{
					loss += val;
					ShowMessage("PLULA-BQP error #1");
				}
			}
		}
	}

	return loss;
}

//does not edit global variables
float PLULA_BQP_loss(int pln, int spnum)
{
	float local, nbh;
	local = PLULA_BQP_local_loss(pln, spnum);
	nbh   = PLULA_BQP_nbh_loss(pln, spnum);

	return (local+nbh);
}

// ------------------ PLULA TREE LOSS ----------------------------
//possibly edits global variables
//PLvec[current].tree_conn_up;
//PLU_up[0];
//defined inside PLULA.h
float PLULA_tree_nbh_loss(int pln, int spnum)
{
	float loss;

	loss  =   traverse_tree_down(pln, 2, false, true, spnum);
	loss +=   traverse_tree_up(pln,   2, false, true, spnum);

	return loss;
}

///no global edits
float get_sp_loss_at_node_downriver(int spnum, int pln, float conn_loss)
{
	int   loop, end,  x, y, s;
	float val,  loss, remaining, orig, mult, mult2;

	orig      = PLvec[pln].orig_tree_conn_up;
	remaining = PLvec[pln].tree_conn_up; // loss downriver depends on BQPR up
	if (orig>0)                          // at downriver upconn changes
	{
		mult = ( spp[spnum].BQP_link[(int)(10000*remaining/orig)]-
			 spp[spnum].BQP_link[(int)(10000*(remaining-conn_loss)/orig)]  );
		if ((mult<0.0f) || (mult>1.0f))
			ShowMessage("Mult out of bounds in module gsland = "+FloatToStr(mult));
	}
	else
	{
		mult = 1;
		//     ShowMessage("here"); // never happened as it never should
	}

	orig      = PLvec[pln].orig_tree_conn_down;
	remaining = PLvec[pln].tree_conn_down;
	if (orig>0)
	{
		mult2 = spp[spnum].BQP_link2[(int)(10000*remaining/orig)]; // downriver downconn does not change
		if ((mult2<0.0f) || (mult2>1.0f))
			ShowMessage("Mult2 out of bounds in module gsland = "+FloatToStr(mult2));
	}
	else
	{
		mult2 = 1.0f;
		//     ShowMessage("here"); // never happened as it never should
	}

	if (PLvec[pln].allocated)  // if PLULA/spdata has been stored
	{
		//      Form1->Memo1->Lines->Add("allocated at loss");
		//      if ((PLvec[pln].datavec[spnum]<0.0f) || (PLvec[pln].datavec[spnum]>1.0f))
		//        ShowMessage("plndv out of bounds in module gsland = "+FloatToStr(PLvec[pln].datavec[spnum]));
		//      Form1->Memo1->Lines->Add("loss-A = "+FloatToStr(mult*PLvec[pln].datavec[spnum]));
		return mult2*mult*PLvec[pln].datavec[spnum];
	}

	//  ShowMessage("here");  // never happened
	loss = 0.0f;
	end = PLvec[pln].start+PLvec[pln].el_cnt;
	for(loop=PLvec[pln].start;loop<end; loop++)
	{
		x = PLX[loop];
		y = PLY[loop];
		//if (Rmax[y][x]==-1)
		if (-1 == status[y][x])
		  continue;
		val = vmat[y][x][spnum];
		if (val<=0.0f)
		  continue;

		loss += val;
	}

	//  loss *= mult;
	//     Form1->Memo1->Lines->Add("NOT allocated at loss");
	return mult2*mult*loss;
}

float get_sp_loss_at_node_upriver(int spnum, int pln, float conn_loss)
{
	int   loop, end,  x, y, s;
	float val,  loss, remaining, orig, mult, mult2;

	orig      = PLvec[pln].orig_tree_conn_down;
	remaining = PLvec[pln].tree_conn_down; // what happens upriver is relevant for
	if (orig>0)                            // spp with up response only
	{                                    // downconn changes for these spp
		mult = ( spp[spnum].BQP_link2[(int)(10000*remaining/orig)]-
			 spp[spnum].BQP_link2[(int)(10000*(remaining-conn_loss)/orig)]  );
		if ((mult<0.0f) || (mult>1.0f))
			ShowMessage("Mult out of bounds in module gsland-ur = "+FloatToStr(mult));
	}
	else
	{
		mult = 1;
		//     ShowMessage("here"); // never happened as it never should
	}

	orig      = PLvec[pln].orig_tree_conn_up;
	remaining = PLvec[pln].tree_conn_up;
	if (orig>0)
	{
		mult2 = spp[spnum].BQP_link[(int)(10000*remaining/orig)]; // upriver upconn does not change
		if ((mult2<0.0f) || (mult2>1.0f))
			ShowMessage("Mult2 out of bounds in module gsland-ur = "+FloatToStr(mult2));
	}
	else
	{
		mult2 = 1.0f;
		//     ShowMessage("here"); // never happened as it never should
	}

	if (PLvec[pln].allocated)  // if PLULA/spdata has been stored
	{
		//      Form1->Memo1->Lines->Add("allocated at loss");
		//      if ((PLvec[pln].datavec[spnum]<0.0f) || (PLvec[pln].datavec[spnum]>1.0f))
		//        ShowMessage("plndv out of bounds in module gsland = "+FloatToStr(PLvec[pln].datavec[spnum]));
		//      Form1->Memo1->Lines->Add("loss-A = "+FloatToStr(mult*PLvec[pln].datavec[spnum]));
		return mult2*mult*PLvec[pln].datavec[spnum];
	}

	//  ShowMessage("here");  // never happened
	loss = 0.0f;
	end = PLvec[pln].start+PLvec[pln].el_cnt;
	for(loop=PLvec[pln].start;loop<end; loop++)
	{
		x = PLX[loop];
		y = PLY[loop];
		//if (Rmax[y][x]==-1)
		if (-1 == status[y][x])
		  continue;
		val = vmat[y][x][spnum];
		if (val<=0.0f)
		  continue;

		loss += val;
	}

	//  loss *= mult;
	//     Form1->Memo1->Lines->Add("NOT allocated at loss");
	return mult2*mult*loss;
}

float PLULA_tree_local_loss(int pln, int spnum)
{
	int   s;
	float mult;

	float orig = PLvec[pln].orig_tree_conn_up;
	float remaining = PLvec[pln].tree_conn_up;
	if(orig > 0)
		mult = spp[spnum].BQP_link[static_cast<int>(10000.0f * remaining / orig)];
	else
		mult = 1.0f;

	orig = PLvec[pln].orig_tree_conn_down;
	remaining = PLvec[pln].tree_conn_down;
	if(orig > 0)
	{
		//      if (remaining<0)
		//        ShowMessage("Negative remaining");
		//      else if (!spp[spnum].BQP_link2)
		//        ShowMessage("No loink");
		//      else
		mult *= spp[spnum].BQP_link2[static_cast<int>(10000.0f * remaining / orig)]; //xxxTrror here, is remaining negative
	}
	else
		mult *= 1.0f;

	if (PLvec[pln].allocated)  // if was stored
	{
		return mult*PLvec[pln].datavec[spnum];
	}

	int end = PLvec[pln].start + PLvec[pln].el_cnt;
	float loss = 0.0f;
	float val = 0.0f;
	for(int loop = PLvec[pln].start; loop < end; ++loop)
	{
		int x = PLX[loop];
		int y = PLY[loop];
		//if (Rmax[y][x] == -1)
		if (-1 == status[y][x])
		  continue;
		val = vmat[y][x][spnum];
		if (val<=0.0f)
		  continue;

		loss += mult*val;
	}

	return loss;
}

float PLULA_tree_loss_only_local(int pln, int spnum)
{
	int   s;
	float mult;

	if (PLvec[pln].allocated)  // if was stored
	{
		return PLvec[pln].datavec[spnum];
	}

	int end = PLvec[pln].start + PLvec[pln].el_cnt;
	float loss = 0.0f;
	float val = 0.0f;
	for(int loop = PLvec[pln].start; loop < end; ++loop)
	{
		int x = PLX[loop];
		int y = PLY[loop];
		//if (Rmax[y][x] == -1)
		if (-1 == status[y][x])
		  continue;
		val = vmat[y][x][spnum];
		if (val<=0.0)
		  continue;

		loss += val;
	}

	return loss;
}

//possibly edits global variables
//PLvec[current].tree_conn_up;
//PLU_up[0];
//defined inside PLULA.h
float PLULA_tree_loss(int pln, int spnum)
{
	if (!spp[spnum].TREE_used)
	{
		return PLULA_tree_loss_only_local(pln, spnum);;
	}
	float local = PLULA_tree_local_loss(pln, spnum);
	float nbh   = PLULA_tree_nbh_loss(pln, spnum);
	return local + nbh;
}

/**	Get planning unit layer information
Does not edit global variables!
\param pln index of the planning unit
\param PLdv an array where delta values are stored
\param PLSSI an array where planning unit ssi indices are stored
\param PLSSIval an array where planning unit ssi values are stored
\param SSIcnt a reference where the number of ssi layers is stored
\return total cost of the planning unit
*/
float get_PLULA_info(int pln, float *PLdv, int *PLSSI, float *PLSSIval, int &SSIcnt)
{
	float total_cost; // total cost of the planning unit

	for(int s = 0; s < map_cnt; ++s) {
		PLdv[s] = 0.0f;
	}

	if(use_SSI) {
		for(int s = 0; s <= SSI_spp_cnt; ++s) {
			PLSSI[s]    = s;
			PLSSIval[s] = 0.0f;
		}
	}

	total_cost = 0.0f;
	SSIcnt = 0;

	int end = PLvec[pln].start + PLvec[pln].el_cnt; // index of the last index in the planning unit
	for(int loop = PLvec[pln].start; loop<end; ++loop) { //for each element
		int x = PLX[loop];
		int y = PLY[loop];
		if (!vmat[y][x])
			continue; // PLULA goes outside data area as in Alps
		const Biodiv_Features_Occur_Container& rowp = vmat[y][x]; // array of species representations
		//for(int s = 0; s < map_cnt; ++s) {
		for(int s = rowp.first(); s != rowp.overflow(); s = rowp.next(s)) {
                        float rowp_s_val = rowp[s];
		        //if(rowp[s] > 0.0f)
			//  PLdv[s] += rowp[s];
			if(rowp_s_val > 0.0f)
			  PLdv[s] += rowp_s_val;
		}

		if (use_SSI) {
			int SSIpos = SSIxyPos[y][x];
			for(int s = 0; s < SSIxyCnt[y][x]; ++s) {
				int SSIspnum = SSI_list[SSIpos+s].spnum;
				float SSIval = SSI_list[SSIpos+s].val;
				if((PLSSIval[SSIspnum] <= 0.0f) && (SSIval > 0.0f))
					++SSIcnt;
				PLSSIval[SSIspnum] += SSIval;
			}
		}

		if (use_cost)
			total_cost += costmap.m[y][x];
		else
			total_cost++;
	}

	if (!PLvec[pln].allocated) {
		PLvec[pln].datavec = new float[map_cnt];
		if (PLvec[pln].datavec != 0) {
			PLvec[pln].allocated = true;
			for(int s = 0; s < map_cnt; ++s)
				PLvec[pln].datavec[s] = PLdv[s];
			PLvec[pln].cost = total_cost;
		}
	}

	return total_cost;
}

void  check_PLULA_sp_data(int s)
{
	float sum;
	int   pln;
	const size_t LINE_LEN = 512; 
	char txt[LINE_LEN];

	sum=0.0f;
	for(pln=0; pln<PLcnt; pln++)
	{
		sum += PLvec[pln].datavec[s];
	}

	snprintf(txt, LINE_LEN, "feature %i, total PLULA distribution sum=%f",s, sum);
	Form1->Memo1->Lines->Add(txt);
}

/// Check whether planning unit layers are allocated
bool get_PLULA_sp_data()
{
  float  PLULA_data_vec[MAX_SPP_COUNT], PLULA_SSI_ps[MAX_SPP_COUNT], PLULA_cost;
  int    PLULA_SSI_spp[MAX_SPP_COUNT], PLULA_SSI_cnt, pln;
  bool   ok;

#ifdef COMPACT_VMAT
  Biodiv_Features_Occur_Container::init_plu(map_cnt);
#endif
	
  ok = true;
  for(pln=0; pln<PLcnt; pln++) {
    get_PLULA_info(pln, PLULA_data_vec,
		   PLULA_SSI_spp, PLULA_SSI_ps, PLULA_SSI_cnt);
    if (!PLvec[pln].allocated)
      ok=false;
  }

  // if (ok)
  //   for (int s=0; s<map_cnt; s++)
  //     check_PLULA_sp_data(s);
  
  return ok;
}

// GBF, generalized benefit function
float rr4(float fract, float w1, float Tj, float w2, float exp1, float exp2)
{
  if ((fract<=Tj) || (Tj>=1.0f))
    return w1*z_pow(fract/Tj, exp1);

  if (0!=w2)
    return (w1 + w2*z_pow((fract-Tj)/(1.0f-Tj), exp2));
  else
    return w1;
}

// TODO: in the end this is like grid_utils.h:get_nb8_count()
inline int
get_plus_minus_nb8_count(char**& sts, int x, int y)
{
  int cnt = 0;
  /* coastline killer;
  cnt = 
    sts[y-1][x-1] +
    sts[y-1][x] + 
    sts[y-1][x+1] +
    sts[y][x-1] +
    sts[y][x+1] +
    sts[y+1][x-1] +
    sts[y+1][x] +
    sts[y+1][x+1];
  */
  if (sts[y-1][x-1]>0)
    cnt++;
  if (sts[y-1][x]>0)
    cnt++;
  if (sts[y-1][x+1]>0)
    cnt++;
  if (sts[y][x-1]>0)
    cnt++;
  if (sts[y][x+1]>0)
    cnt++;
  if (sts[y+1][x-1]>0)
    cnt++;
  if (sts[y+1][x]>0)
    cnt++;
  if (sts[y+1][x+1]>0)
    cnt++;

  return cnt;
}

//possibly edits the following global variables
//BQP_spp_loss_tmp, defined here
//PLvec[current].tree_conn_up
//PLU_up[0]
//defined inside PLULA.h
float get_delta_value(int x, int y, bool single_cell_mode)
{
  float delta_value;
  if (ADM_set.use_ADMUs) {
    delta_value = delta_value_new(x, y, single_cell_mode);
  } else {
	/// planning unit index of the element
	int pln;

	/// species representation values in the element
	Biodiv_Features_Occur_Container rowp;

	float  orig, nbv_val, dv_val, tmpf, SSIval, oldval, newval;
	int    s, remaining, SSIspnum, SSIpos, end;
	float  min_delta_value; // 10.12.2009; added for CAZ MCZ

	struct sp *psp;
	// local xxxPLULA variables
	float  PLULA_data_vec[MAX_SPP_COUNT], PLULA_SSI_ps[MAX_SPP_COUNT];

	/// planning unit cost
	float PLULA_cost;
	int    PLULA_SSI_spp[MAX_SPP_COUNT], PLULA_SSI_cnt;
	bool   use_PLULA_mode;

	// ADMUs
	int   adu;
	float Wjadu;

	if (!single_cell_mode) // this is for removal, needs to evaluate single cell with BQP
		use_PLULA_mode = use_PLULA;
	else
		use_PLULA_mode = false;

	if (ADM_set.use_ADMUs){
	  adu = ADMUs[y][x];
	  if (adu<0 || adu>=ADM_set.count) // if wrong ADMU id, will use weights of the first admu!
	    adu=0;
	}
	else
		Wjadu = 1.0;

	if (use_PLULA_mode)
	{
		pln = PLL[y][x];
		if (pln<0)
		{
			ShowMessage("GDV PLULA number error");
			return 1.0f;
		}
		if (PLvec[pln].checked)  // all cells in the PLULA return same delta value
			return PLvec[pln].dv;
	}

	if (!use_PLULA_mode)  // xxx note: ADMUs work with PLUs because of rowp is fine and indep of ADM weights!
		rowp = vmat[y][x];
	else
	{
		if (use_SSI || (!PLvec[pln].allocated)) // xxx storing of SSI data with
		{                                     // plula has not been implemented
		  PLULA_cost = get_PLULA_info(pln, PLULA_data_vec,
					      PLULA_SSI_spp, PLULA_SSI_ps, PLULA_SSI_cnt);
		  // rowp = PLULA_data_vec;   // COMPACT_VMAT
		  rowp.assign_seq(PLULA_data_vec, map_cnt);			
		}
		else // spdata has been precomp + stored (but not yet SSI data)
		{
		  //Form1->Memo1->Lines->Add("heer");
		  // rowp = PLvec[pln].datavec; // COMPACT_VMAT
		  rowp.assign_seq(PLvec[pln].datavec, map_cnt);		  
		  PLULA_cost = PLvec[pln].cost;
		}
	}

	if (removal_rule==1) // original core area zonation  // xxxMCZ Q resolved 10.12.2009
	{
		delta_value     = -std::numeric_limits<float>::max();
		min_delta_value = std::numeric_limits<float>::max();

		//for(s=0; s<map_cnt; s++) // for all species layers
		for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))
		{
		  if (0==spp[s].weight)
		    continue;

                        float rowp_s_val = rowp[s];
		        //if (rowp[s]<0.0f)
			if (rowp_s_val<0.0f)
				continue;

			if (ADM_set.use_ADMUs)
				Wjadu = ADM_combined_weights[s][adu];
			psp = &spp[s];

			if ((use_PLULA_mode) && (use_tree_conn))
			{
				delta_value = max(delta_value, PLULA_tree_loss(pln, s)*smult[s]*Wjadu);
				min_delta_value = min(min_delta_value, PLULA_tree_loss(pln, s)*smult[s]*Wjadu);
			}
			else if (!use_BQP || !psp->BQP_used)  // no BQP, PLULA or not
			{
			        //delta_value = max(delta_value, rowp[s]*smult[s]*Wjadu);
				//min_delta_value = min(min_delta_value, rowp[s]*smult[s]*Wjadu);
			  float this_val = rowp_s_val*smult[s]*Wjadu;
#ifdef CAZ_ALTERNATING_EPS
			  {
			    //int nb_count = get_plus_minus_nb8_count(status, x, y);
			    int nb_count_shifted = get_nb8_count(status, x, y) - 3;   // nb8_count is in [0,7]
			    float round_bit = nb_count_shifted * this_val * std::numeric_limits<float>::epsilon();
			    this_val += round_bit;
			  }
			  /*else 
			    this_val -= round_bit;*/
			      // std::numeric_limits<float>::round_error(); is too big and dominating!
#endif
				delta_value = max(delta_value, this_val);
				min_delta_value = min(min_delta_value, this_val);
			}
			else
			{
				if (use_PLULA_mode)
				{
					delta_value = max(delta_value, PLULA_BQP_loss(pln, s)*smult[s]*Wjadu);
					min_delta_value = min(min_delta_value, PLULA_BQP_loss(pln, s)*smult[s]*Wjadu);
				}
				else
				{
					remaining = nbms[nbm_num[s]][y][x]; // number of nbs remaining
					orig      = nbms_orig[nbm_num[s]][y][x];
					if (orig>0)
					{
						nbv_val   = psp->BQP_link[(int)(10000*remaining/orig)]; // matches initial #nbs and init prob
						dv_val    = BQP_dv_buf(x,y,s);
					}
					else
					{
						nbv_val   = 1.0f;
						dv_val    = 0.0f;
					}
					dv_val    = BQP_dv_buf(x,y,s);
					//tmpf      = nbv_val*rowp[s]+dv_val;
					tmpf      = nbv_val*rowp_s_val+dv_val;
					delta_value     = max( delta_value, tmpf*smult[s]*Wjadu);  // xxxMCZ; how combines with wj<0???
					min_delta_value = min(min_delta_value, tmpf*smult[s]*Wjadu);  // xxxMCZ; how combines with wj<0???
					BQP_spp_loss_tmp[s] = tmpf; //global edit!
				}
			}			
			
		}		
		if (neg_weights_used)
		{
			if (min_delta_value<0.0)
				delta_value += min_delta_value; // MCZ 10.12.2009; mdv<0, thus addition ok.
		}										// missing 2nd variant; summing of neg effects
	}
	else if (removal_rule==2) // convex BF with par
	{
		delta_value=0.0f;
		// for(s=0; s<map_cnt; s++)
		for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))
		{
                        float rowp_s_val = rowp[s];
			//if (rowp[s]<0.0f)
			if (rowp_s_val<0.0f || 0==spp[s].weight)
				continue;
			if (ADM_set.use_ADMUs)
				Wjadu = ADM_combined_weights[s][adu];
			psp = &spp[s];
			// ax^b => D[]=a*b*x^(b-1)
			if ((use_PLULA_mode) && (use_tree_conn))
				delta_value += PLULA_tree_loss(pln, s)*sp_derivative[s]*Wjadu;
			else if (!use_BQP || !psp->BQP_used) // xxx error, not rowps
			        //delta_value += rowp[s]*sp_derivative[s]*Wjadu;
				delta_value += rowp_s_val*sp_derivative[s]*Wjadu;
			else
			{
				if (use_PLULA_mode)
					delta_value += PLULA_BQP_loss(pln, s)*sp_derivative[s]*Wjadu;
				else
				{
					remaining = nbms[nbm_num[s]][y][x]; // number of nbs remaining
					orig      = nbms_orig[nbm_num[s]][y][x];
					if (orig>0)
					{
						nbv_val   = psp->BQP_link[(int)(10000*remaining/orig)]; // matches initial #nbs and init prob
						dv_val    = BQP_dv_buf(x,y,s);
					}
					else
					{
						nbv_val   = 1.0f;
						dv_val    = 0.0f;
					}

					//tmpf      = nbv_val*rowp[s]+dv_val;
					tmpf      = nbv_val*rowp_s_val+dv_val;
					BQP_spp_loss_tmp[s] = tmpf; //global edit!
					delta_value += tmpf*sp_derivative[s]*Wjadu;
				}
			}
		}
	}
	else if (removal_rule==3) // target BF
	{
		delta_value=0.0f;
		// for(s=0; s<map_cnt; s++)
		for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))
		{
                        float rowp_s_val = rowp[s];
			//if (rowp[s]<0.0f)
			if (rowp_s_val<0.0f || 0==spp[s].weight)
				continue;
			psp = &spp[s];
			// ax^b => D[]=a*b*x^(b-1)
			if ((use_PLULA_mode) && (use_tree_conn))
			{
				tmpf = PLULA_tree_loss(pln, s)*sp_derivative[s];
				if ((repr_lvls[s]-tmpf)>psp->sr_par)
					delta_value += tmpf*sp_derivative[s];
				else
				{
					if (repr_lvls[s]<psp->sr_par)
						delta_value += 0.0f;
					else
						delta_value += (map_cnt+SSI_spp_cnt);
				}
			}
			else if (!use_BQP || !psp->BQP_used) // xxx error, not rowps
			{
			        //if ((repr_lvls[s]-rowp[s])>psp->sr_par)
				//	delta_value += rowp[s]*sp_derivative[s];
				if ((repr_lvls[s]-rowp_s_val)>psp->sr_par)
					delta_value += rowp_s_val*sp_derivative[s];
				else
				{
					if (repr_lvls[s]<psp->sr_par)
						delta_value += 0.0f;
					else
						delta_value += (map_cnt+SSI_spp_cnt);
				}
			}
			else
			{
				if (use_PLULA_mode)
					tmpf = PLULA_BQP_loss(pln, s);
				else
				{
					remaining = nbms[nbm_num[s]][y][x]; // number of nbs remaining
					orig      = nbms_orig[nbm_num[s]][y][x];
					if (orig>0)
					{
						nbv_val   = psp->BQP_link[(int)(10000*remaining/orig)]; // matches initial #nbs and init prob
						dv_val    = BQP_dv_buf(x,y,s);
					}
					else
					{
						nbv_val   = 1.0f;
						//                     dv_val    = 0.0f;
					}
					dv_val    = BQP_dv_buf(x,y,s);

					//tmpf      = nbv_val*rowp[s]+dv_val;
					tmpf      = nbv_val*rowp_s_val+dv_val;
					BQP_spp_loss_tmp[s] = tmpf; //global edit!
				}

				if ((repr_lvls[s]-tmpf)>psp->sr_par)
					delta_value += tmpf*sp_derivative[s];
				else
				{
					if (repr_lvls[s]<psp->sr_par)
						delta_value += 0.0f;
					else
						delta_value += (map_cnt+SSI_spp_cnt);
				}
			}
		}
	}
	else if (removal_rule==4) // xxxMCZ problem, really. (function inversion; (1-r) inversion?)
	{                         // 10.12.2009; no prob? Just invert power function parameters, as in UK MCZ paper.
		delta_value=0.0f;
		// for(s=0; s<map_cnt; s++)
		for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))
		{
		        float rowp_s_val = rowp[s];
			//if (rowp[s]<0.0f)
			if (rowp_s_val<0.0f)
				continue;
			if (ADM_set.use_ADMUs)
				Wjadu = ADM_combined_weights[s][adu];
			psp = &spp[s];
			// ax^b => D[]=a*b*x^(b-1)
			if ((use_PLULA_mode) && (use_tree_conn))
			{
				tmpf = PLULA_tree_loss(pln, s);
				//              delta_value += PLULA_tree_loss(pln, s)*sp_derivative[s];
			}
			else if (!use_BQP || !psp->BQP_used)
			{
			        //tmpf = rowp[s];
				tmpf = rowp_s_val;
				//              delta_value += rowp[s]*sp_derivative[s];
			}
			else
			{
				if (use_PLULA_mode)
				{
					// delta_value += PLULA_BQP_loss(pln, s)*sp_derivative[s];
					tmpf = PLULA_BQP_loss(pln, s);
				}
				else
				{
					remaining = nbms[nbm_num[s]][y][x]; // number of nbs remaining
					orig      = nbms_orig[nbm_num[s]][y][x];
					if (orig>0)
					{
						nbv_val   = psp->BQP_link[(int)(10000*remaining/orig)]; // matches initial #nbs and init prob
						dv_val    = BQP_dv_buf(x,y,s);
					}
					else
					{
						nbv_val   = 1.0f;
						dv_val    = 0.0f;
					}

					//tmpf      = nbv_val*rowp[s]+dv_val;
					tmpf      = nbv_val*rowp_s_val+dv_val;
					BQP_spp_loss_tmp[s] = tmpf;
					//                  delta_value += tmpf*sp_derivative[s];
				}
			}
			if (ADM_set.use_ADMUs)
				delta_value += (   rr4(repr_lvls[s],       Wjadu, spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2)
						   - rr4(repr_lvls[s]-tmpf,  spp[s].weight, spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2) );
			else
				delta_value += (   rr4(repr_lvls[s],       spp[s].weight, spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2)
						   - rr4(repr_lvls[s]-tmpf,  spp[s].weight, spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2) );
			// xxxMCZ ok? No sp_derivative?
			//float rr4(float fract, float w1, float Tj, float w2, float exp1, float exp2)
		}  // sp loop
	} // rr4

	if (use_SSI)  // xxxMCZ does not prevent wj<0 for SSI spp. SSI.wj<0 should work ok. 5/09
	{
		if (use_PLULA_mode)
			end = SSI_spp_cnt;
		else
		{
			end    = SSIxyCnt[y][x];
			SSIpos = SSIxyPos[y][x];
		}

		for(s=0; s<end; s++)
		{

			if (use_PLULA_mode)  // note end = SSI_spp_cnt
			{
				SSIspnum = s;
				SSIval   = PLULA_SSI_ps[s];
				if (SSIval<=0.0f)
					continue;
			}
			else
			{
				//          if ((SSIpos<0) || (SSIpos>SSI_lpos))
				//            ShowMessage("SSIerror");
				SSIspnum  = SSI_list[SSIpos+s].spnum;
				SSIval    = SSI_list[SSIpos+s].val;
			}

			if (removal_rule==1) // original core area zonation
			{
				//              Form1->Memo1->Lines->Add("predeltassival="+FloatToStr(delta_value));
				newval = SSIval*SSI[SSIspnum].weight/SSI_repr_lvls[SSIspnum];
				if (newval>delta_value)
					delta_value=newval;
				//              sprintf(txt, "(%i,%i) sp=%i, dv=%f %f", y,x,SSIspnum,delta_value,
				//                SSIval*SSI[SSIspnum].weight/SSI_repr_lvls[SSIspnum]);
				//              Form1->Memo1->Lines->Add(txt);
			}
			else if (removal_rule==2) // convex BF with par
			{
				//              if (SSI_repr_lvls[SSIspnum]>0.0f)
				delta_value += SSI[SSIspnum].weight*        // xxx problem with pow
						( z_pow(SSI_repr_lvls[SSIspnum], SSI[SSIspnum].sr_par) -
						  z_pow(SSI_repr_lvls[SSIspnum]-SSIval, SSI[SSIspnum].sr_par)   );
			}
			else if (removal_rule==3) // target BF
			{
				if ((SSI_repr_lvls[SSIspnum]-SSIval)>SSI[SSIspnum].sr_par)
				{ // here derivative not ok because can be large step change in repr
					oldval = z_pow( (SSI_repr_lvls[SSIspnum]-SSI[SSIspnum].sr_par)
						      /(1.0f-SSI[SSIspnum].sr_par), 0.25f);
					newval = z_pow( (SSI_repr_lvls[SSIspnum]-SSIval-SSI[SSIspnum].sr_par)
						      /(1.0f-SSI[SSIspnum].sr_par), 0.25f);
					delta_value += (oldval-newval);
				}
				else
				{
					if (SSI_repr_lvls[SSIspnum]<SSI[SSIspnum].sr_par)
						delta_value += 0.0f;
					else
						delta_value += (map_cnt+SSI_spp_cnt);
				}
			}
			else if (removal_rule==4)
			{
				delta_value += (   rr4(SSI_repr_lvls[SSIspnum],         SSI[SSIspnum].weight, SSI[SSIspnum].Tj, SSI[SSIspnum].w2, SSI[SSIspnum].exp1, SSI[SSIspnum].exp2)
						   - rr4(SSI_repr_lvls[SSIspnum]-SSIval,  SSI[SSIspnum].weight, SSI[SSIspnum].Tj, SSI[SSIspnum].w2, SSI[SSIspnum].exp1, SSI[SSIspnum].exp2) );
			}
		}  // SSI spp
	} // use SSI

	//float rr4(float fract, float w1, float Tj, float w2, float exp1, float exp2)

	if (use_PLULA_mode)
	{
		delta_value /= PLULA_cost; // or implicitly sum of area, if cost not used
	}
	else if (use_cost)
	{
#ifdef COST_DIV_DOUBLE
		delta_value = float(double(delta_value) / double(costmap.m[y][x]));
#else
		delta_value /= costmap.m[y][x];
#endif
	}

	if (use_PLULA_mode)
	{
		PLvec[pln].checked = true;
		PLvec[pln].dv      = delta_value;
	}

  }

  if (isnan(delta_value)) {
#ifndef ZCORE_DEBUG
    Form1->Memo1->Lines->Add("EXITING: ZIG has detected a NaN delta value, x: " + 
			     IntToStr(x) + ", y: " + IntToStr(y) + ", single cell mode: " +
			     IntToStr(single_cell_mode) + 
			     ", with fraction of landscape lost: " + 
			     FloatToStr(curves[c_pos][0]));
    graceful_exit();
#else
    if (0 == nan_count )
      Form1->Memo1->Lines->Add("First NaN, x: " + IntToStr(x) + ", y: " + 
			       IntToStr(y) + ", single cell mode: " +
			       IntToStr(single_cell_mode) + 
			       ", with fraction of landscape lost: " + 
			       FloatToStr(curves[c_pos][0]));
    nan_count++;
    if (0 == nan_count%100000)
      Form1->Memo1->Lines->Add("NaN count: " + IntToStr(nan_count));
#endif
  }

  return delta_value;
}

inline float
get_delta_value_plus_penalties_fast(int x, int y, bool single_cell_mode, float bl_inc, float removal_threshold)
{
  float tempval = get_delta_value(x, y, false);

  if (5 == removal_rule)
    tempval = randz();

  if (BLP > 0.0f) {
    tempval -= bl_inc;
  }

  tempval += cc_accounting->calc_ccsp_fast(x, y, removal_threshold);

  return tempval;
}

// values < removal_threshold will be removed
inline float
get_delta_value_plus_penalties(int x, int y, bool single_cell_mode, float bl_inc, bool calc_morpho_penalties, float removal_threshold)
{
  float tempval = get_delta_value(x, y, false);

  if (5 == removal_rule)
    tempval = randz();

  if (BLP > 0.0f) {
    tempval -= bl_inc;
  }

  if (CCSP > 0.0f && corr_ongoing) {
    if (calc_morpho_penalties) {
      // += (rather than penalty, bonus for connecting-cells)
      tempval +=  cc_accounting->calc_ccsp_as_if_removed(x, y, tempval, removal_threshold);
      
      /* TODO: wdirty penalization -> replaced by wdirty check in update_wlist
	 if (wdirty && 1 == wdirty[y][x]) tempval -= 1;
      */
      
    } else {
      if (!CCSP_slow_as_fast && 
	  Grid_CCSP::Slow == cc_accounting->calc_mode() && 
	  cc_accounting->is_cached_penalty(x, y)) {
	tempval += cc_accounting->calc_cached_penalty(x, y);
      }
    }
  }
  return tempval;
}

void update_wlist_sorted(int x, int y, float val,
			 float &max_val, int warp_rnds, int nbc)
{
  // use multiset::lower_bound??
  // priority_queue: 
  // NOP : sort_heap creates a data st. that only supports pop
  // It could be a priority queue if it wasn't because we need access to both front and back (max_warp_val)
  // We'll need to retrieve the warp elements / iterate through them, but no need for pop at all
  //
  // Best solution in theory: std::multimap
  // BUT since insertions happen only once... ???
  // SIMPLE SOL: std::vector, keep track of the max. all the time, and at the end, std::sort the vector.
  // but aren't there plenty of insertion attempts that fail???  <----- YES, but the failure is easy to detect in advance without any lookup, just by comparing against current_max_warp
  /*
struct Comp_Warp
{
public:
  
  bool operator() (int x, int y) const
  {
    return x > y;
    // compare wval, and in case of tie, wnbc
  }
};

 bool warp_comp_func() (int x, int y) const
 {
   return x > y;
   // compare wval, and in case of tie, wnbc
 }

  std::vector<int> w;
  Comp_Warp comp;
  //std::sort(w.begin(), w.end(), comp);
  std::sort(w.begin(), w.end(), warp_comp_func);
  */
}

// insert (x,y,val) in the warp list, with threshold max_val (I/O)
void update_wlist(int x, int y, float val, 
		  float &max_val, int warp_rnds, int nbc)
{
	int   ipos, w, max_nbc;
	float maxv;

		//wdirty[warp_ll] = false;
		if (!wdirty) {
		  wdirty = imatrix(0,yd,0,xd);
		  //wdirty = fmatrix(0,yd,0,xd);
		  for(int y=0;y<yd;y++) {
		    for(int x=0;x<xd;x++) {
		      wdirty[y][x] = 0;
		    }
		  }
		}

		/* Not used anymore...
		// Don't put dirty cells in the wlist
		if (wdirty[y][x] > 0) {
		  wdirty[y][x]--;
		  return;
		}
		*/

		// ******************** TODO: should be done before calling update_wlist!!! ************************** New version:
		// NOOOOOOOOOO, YOU CANNOT COMPARE AGAINT THE *LAST*/OUTDATED CCSP_last_removal_threshold (warp_max_val from Remove...)
		/*
		if (CCSP > 0.0f && cc_accounting->is_cached_penalty(x, y)) {
		  // a bit redundant, this also checks is_cached_penalty
		  //float tempval = get_delta_value_plus_penalties(x, y, false, BL_delta_val[nb_count], false, 0);
		  float tempval = get_delta_value_plus_penalties(x, y, false, 0, false, 0);
		  if (tempval > CCSP_last_removal_threshold )
		    if (CCSP_verbose)
		      std::cerr << " ********************  ******************** dirty to wlist!!! ********************  ********************  " << std::endl;
		    return;
		}
		*/
		
		/* HUGE RISK THAT NEAR THE END OF THE REMOVAL PROCESS WARP_LL=0
		if (cc_accounting->is_cached_penalty(x, y))
		  return;
		*/


		/*
		if (0!= wdirty[y][x] && wdirty[y][x] < ((float)removed/(float)nonm1)) {
		    return;
		  }
		wdirty[y][x]=0;
		*/

	if (warp_ll<warp_rnds)
	{
		wrx[warp_ll]  = x;
		wry[warp_ll]  = y;
		wval[warp_ll] = val;
		wnbc[warp_ll] = nbc;
		warp_ll++;
		return;
	}

	//  if ((x%100)==0)
	//    Form1->Memo1->Lines->Add("Here");
	ipos = -1;
	maxv = -std::numeric_limits<float>::max();
	for(w=0; w<warp_rnds; w++)
	{
		if (wval[w]>maxv)
		{
			maxv=wval[w];
			ipos=w;
		}
	}

	if ((val>maxv) || (maxv==-std::numeric_limits<float>::max()))
		return;

	//  if ((maxv>0.0f) && (val<maxv)) // xxxMCZ condition changed to allow delta<0.0f
	if (val<maxv)
	{
		wrx[ipos]  = x;
		wry[ipos]  = y;
		wval[ipos] = val;
		wnbc[ipos] = nbc;
		max_val=-std::numeric_limits<float>::max();
		for(w=0; w<warp_rnds; w++)
			if (wval[w]>max_val)
				max_val=wval[w];
		return;
	}

	ipos = -1;
	max_nbc=-1;
	for(w=0; w<warp_rnds; w++) // break ties with nbc
	{
		if (wval[w]==maxv)
		{
			if (wnbc[w]>max_nbc)
			{
				max_nbc=wnbc[w];
				ipos=w;
			}
		}
	}

	if (nbc>=max_nbc)
		return;

	if (ipos!=-1)
	{
		wrx[ipos]  = x;
		wry[ipos]  = y;
		wval[ipos] = val;
		wnbc[ipos] = nbc;
		max_val=-std::numeric_limits<float>::max();
		for(w=0; w<warp_rnds; w++)
			if (wval[w]>max_val)
				max_val=wval[w];
	}

	// else if ((x%100)==0)
	//    Form1->Memo1->Lines->Add(FloatToStr(maxv)+"   "+FloatToStr(val));
}

void fix_BL_delta(float *bldv)
{
	int remaining;

	remaining =nonm1-removed;
	if (remaining==1)
		remaining = 2;

	// the index in the bldv[] array is the nb - num
	bldv[0] = BLP*(BL/(float)remaining - (BL-4.0f)/(remaining-1.0f));
	bldv[1] = BLP*(BL/(float)remaining - (BL-2.0f)/(remaining-1.0f));
	bldv[2] = 0;
	bldv[3] = BLP*(BL/(float)remaining - (BL+2.0f)/(remaining-1.0f));
	bldv[4] = BLP*(BL/(float)remaining - (BL+4.0f)/(remaining-1.0f));
}

class SearchIteration
{

private:

	///global variables
	float &val;
	int &rpos;
	float &minval;
	int &rx;
	int &ry;
	int &max_nbc;
	float &min_repr;
	float &min_smult;
	float &ranval;
	float *BL_delta_val; //an array
	int &warp_rnds;
	
	///mutex for global references
	boost::mutex &global_resource_mutex;

public:

	SearchIteration(
			float &val,
			int &rpos,
			float &minval,
			int &rx,
			int &ry,
			int &max_nbc,
			float &min_repr,
			float &min_smult,
			float &ranval,
			float *BL_delta_val,
			int &warp_rnds,
			boost::mutex &global_resource_mutex
			):
		val(val),
		rpos(rpos),
		minval(minval),
		rx(rx),
		ry(ry),
		max_nbc(max_nbc),
		min_repr(min_repr),
		min_smult(min_smult),
		ranval(ranval),
		BL_delta_val(BL_delta_val),
		warp_rnds(warp_rnds),
		global_resource_mutex(global_resource_mutex)
	{
	}

	/// search
	void search(int start, int end) {
		for(int pos(start); pos < end; ++pos) {
			if (run_mode==2) { // solution reload
				rx=load_order[removed].x;
				ry=load_order[removed].y;
#ifdef USE_MAT_EDGE_LIST_POS
				rpos = mat_edge_list_pos[ry][rx];
#else
				for(int loop=0; loop<ecnt; loop++) {
					if ((exl[loop]==rx) && (eyl[loop]==ry)) {
						rpos=loop;
						break;
					}
				}
#endif
				if (use_BQP) {
					val=get_delta_value(rx, ry, false);
					for(int sloop = 0; sloop < map_cnt; ++sloop)
						BQP_spp_losses[sloop] = BQP_spp_loss_tmp[sloop];
				}

				break; // does not need to check all edge points
			}

			int x(exl[pos]);
			int y(eyl[pos]);
			int nb_count = 0;

			if (use_mask && 2!=run_mode) // in re-load mode, the edge list is already sorted according to mask (see Get_sorted_cells_in_map())
			{
				int mlvl = static_cast<int>(maskmap.m[eyl[pos]][exl[pos]]);
				if (mlvl>current_min_mask_level)
					continue; // higher mask cannot be removed
			}

			if (use_PLULA) {
				if (PLL[y][x]<0)
					continue;
				if (PLvec[PLL[y][x]].removed)
					continue;
			}
			//          if (status[y][x]==0)
			//            continue;

			if (BLP>0.0f)
			  nb_count = get_nb_count(status, x, y);

#if 1
			if ((!use_BQP) && (removal_rule==1) && (!neg_weights_used) && !ADM_set.use_ADMUs) { // xxx presently no shortcut for BF formulation
				float min_contr(min_smult*Rmax[y][x]);
				if (BLP>0.0f)
					min_contr -= BL_delta_val[nb_count]; // all old mask stuff removed from here after mask_lvls used
				if (min_contr>minval)  // this is fine even with PLULA, minval will just be large
					continue;  // xxxMCZ error, nok with wj<0.0!!! Thus, !nwu added to if condition.
			} // not problem for mask if does not work
#endif

			// Does not include morphological values/penalties
			float tempval = 0;
			//corr_ongoing = false;
			if (CCSP > .0f && !corr_ongoing /* && NULL == cc_accounting*/) {
			  float prop_remaining = 100.0f * (1.0f - (rem_level + float(current)/float(nonm1+1.0f)));
			  bool corr_start = (prop_remaining <= corr_set.start_pc);
			  if (corr_start && NULL == cc_accounting)
			    init_corridors();
			}
			//corr_ongoing = true;
			if ( CCSP > .0f && corr_ongoing) {    // TODO: simplify the ifs when done with the tests
			  if (Grid_CCSP::Fast_Warp == cc_accounting->calc_mode() || 
			      (CCSP_slow_as_fast && Grid_CCSP::Slow == cc_accounting->calc_mode()) ) {
			    tempval = get_delta_value_plus_penalties_fast(x, y, false, BL_delta_val[nb_count], 0);
			  } else if (Grid_CCSP::Slow == cc_accounting->calc_mode()) {
			    tempval = get_delta_value_plus_penalties(x, y, false, BL_delta_val[nb_count], false, 0);
			  }
			} else {
			  tempval = get_delta_value_plus_penalties(x, y, false, BL_delta_val[nb_count], false, 0);
			}
			// TODO: try (expensive)
			// float tempval = get_delta_value_plus_penalties(x, y, false, BL_delta_val[nb_count], true, warp_max_val);


			// xxx a load of old mask stuffitz removed from here
			//printf("warp: %d, pos: %d, val: %f, minval: %f\n", warp_rnds, pos, val, minval);
			//here we synchronize, because global references are changed
			boost::lock_guard<boost::mutex> guard(global_resource_mutex);

			if (tempval>minval)
				continue;

			if (warp_rnds>1) {
			  update_wlist(x, y, tempval, minval, warp_rnds, nbm[y][x]);
			  //               Form1->Memo1->Lines->Add("Warping");
			  continue;
			}

			bool save(false);
			int nbc(nbm[y][x]);
			if (tempval < minval) {
				save   = true;
				ranval = 1.0f;
			}
			else { // val must now be == minval
				if (nbc<max_nbc) { // #immediate neighbors breaks tie.
					save  = true;
					ranval= 1.0f;
				}
				else if (nbc==max_nbc) {
					float rtmp(randz()); // equal value surface randomization
					if (rtmp<ranval) {
						save  = true;
						ranval= rtmp;
					}
					ties_cnt++;
				}
			}

			if (save) { // did find new best cell for removal
				//             Form1->Memo1->Lines->Add(FloatToStr(val)+"  "+FloatToStr(minval)+"  "+FloatToStr(BLP));
				rpos = pos;
				minval = tempval;
				rx = x;
				ry = y;
				max_nbc = nbc;
				val = tempval; //synch addition - change global reference
				if (use_BQP)
					for(int sloop = 0; sloop < map_cnt; ++sloop)
						BQP_spp_losses[sloop]=BQP_spp_loss_tmp[sloop];
			}
		}
	}
};

struct IterationStepParameters
{
	int start_pos;
	int end_pos;
	IterationStepParameters() {}
	IterationStepParameters(int start_pos, int end_pos):
		start_pos(start_pos), end_pos(end_pos) {}
};

class ThreadQueue
{

private:

	typedef boost::unique_lock<boost::mutex> unique_lock;

	///thread pool
	unsigned int pool_size;
	unsigned int idle_count;
	boost::thread_group pool;
	bool running; ///< true during pool lifetime, set to false on destruction

	///pool queue
	std::queue<IterationStepParameters> pos_queue;

	///pool sync
	boost::mutex pool_mutex;
	boost::condition_variable empty_cond; ///< queue is not full
	boost::condition_variable full_cond; ///< queue is not empty

	///iteration specifics
	std::shared_ptr<SearchIteration> iteration;
	boost::mutex global_resource_mutex; ///< mutex for global references
	unsigned int iteration_round;

public:

	ThreadQueue(unsigned int pool_size):
		pool_size(pool_size),
		idle_count(0),
		running(true),
		iteration_round(0)
	{
		for(unsigned int i(0); i < pool_size; ++i) {
			// *this is callable, but can't be copied,
			// hence the boost::ref
			pool.add_thread(new boost::thread(boost::ref(*this)));
		}
	}

	/// notify all threads to stop immediately, join and terminate
	/// has to be called from the initiating thread
	~ThreadQueue() {
		{
			unique_lock lock(pool_mutex); // grab mutex
			running = false; // terminate processing
			empty_cond.notify_all(); // wake up worker threads
		}
		pool.join_all(); // wait until threads hev finished
	}

	/// has to be called from the initiating thread
	void nextIteration(
			float &val,
			int &rpos,
			float &minval,
			int &rx,
			int &ry,
			int &max_nbc,
			float &min_repr,
			float &min_smult,
			float &ranval,
			float *BL_delta_val,
			int &warp_rnds)
	{
		iteration.reset(new SearchIteration(
					val, rpos, minval, rx, ry, max_nbc, min_repr, min_smult, ranval, BL_delta_val, warp_rnds,
					global_resource_mutex));
		++iteration_round;
	}

	/// terminates iteration and waits until queue is empty AND all threas are idle
	/// has to be called from the initiating thread
	void completeIteration()
	{
		unique_lock lock(pool_mutex);
		empty_cond.notify_all(); //just in the unlikely case a thread is waiting while queue is not empty
		while(pos_queue.size() > 0 || idle_count < pool_size) {
			full_cond.wait(lock);
		}
	}

	/// create new execution thread
	/// return whether loop is still running
	void next(int start_pos, int end_pos) {
		unique_lock lock(pool_mutex);
		while(pos_queue.size() == pool_size) {
			full_cond.wait(lock);
		}
		pos_queue.push(IterationStepParameters(start_pos, end_pos));
		empty_cond.notify_one();
	}

	/// thread execution entry point
	void operator()() {
		while(true) {
			IterationStepParameters pos;
			{
				boost::unique_lock<boost::mutex> lock(pool_mutex);
				full_cond.notify_one();
				++idle_count;
				while(pos_queue.size() == 0 && running)
					empty_cond.wait(lock);
				// quit thread if loop has been terminated
				full_cond.notify_one();
				--idle_count;
				if(!running)
					break;
				pos = pos_queue.front();
				pos_queue.pop();
			}
			// this search is done in blocks created in ThreadQueue::next
			iteration->search(pos.start_pos, pos.end_pos);
		}
	}

	unsigned int size() {
		return pool_size;
	}
};


// here is where the grid/edge list may be split into blocks for different threads
// val: 1-rank
// val will be assigned the 'min representation across features' at the end
bool get_next_to_remove(float &val, ThreadQueue &queue)
{
	float minval, min_repr, min_smult, ranval;
	int   max_nbc;
	int   rx, ry, rpos;
	float BL_delta_val[5];

	if (BLP>0.0f)
		fix_BL_delta(BL_delta_val);

	if (use_PLULA) {
		for(int loop = 0; loop < PLcnt; ++loop)
		{
			PLvec[loop].checked = false;
			PLvec[loop].dv = -1.0f;
		}
	}

	max_nbc  = 8;
	rx=ry    = -1;
	max_repr = -1;
	min_repr = std::numeric_limits<float>::max();
	min_smult= std::numeric_limits<float>::max();
	minval   = std::numeric_limits<float>::max();

	if (2!=run_mode || time_to_print_info_line(current+1)) {  
	  // avoid if in re-load mode. This loop is not needed, and for many spp (like ~10^4) it takes too much
	  //  (exception: if the current min/etc. repr. value is needed to print info)
	  for(int s=0; s<map_cnt; s++) {
	    if (spp[s].weight>=0.0f) {
	      min_repr = min(min_repr, (float)repr_lvls[s]);
	      max_repr = max(max_repr, (float)repr_lvls[s]);
	    }
	    smult[s] = 0.0f;
	    sp_derivative[s] = 0.0f;
	    if (0.0f==spp[s].weight)
	      continue;
	    if (repr_lvls[s]>0.0f) {
	      if (!ADM_set.use_ADMUs)
		smult[s] = spp[s].weight/repr_lvls[s]; // xxxWj
	      else
		smult[s] = 1.0f/repr_lvls[s]; // xxxWj
	      
	      if (smult[s]<min_smult)
		min_smult=smult[s];
	      
	      if (removal_rule==2) { // MCZ directly in derivative
		if (!ADM_set.use_ADMUs) {
		  sp_derivative[s] = z_pow(repr_lvls[s], spp[s].sr_par-1.0f)
		    *spp[s].sr_par*spp[s].weight;
		} else {
		  sp_derivative[s] = z_pow(repr_lvls[s], spp[s].sr_par-1.0f)
		    *spp[s].sr_par;
		}
	      }
	      else if (removal_rule==3) { // xxxMCZ reverse functions? Done 10.12.2009
		//
		if (spp[s].weight>=0.0) {
		  if (repr_lvls[s]<=spp[s].sr_par)
		    sp_derivative[s]=0.0;
		  else
		    sp_derivative[s] = z_pow((repr_lvls[s]-spp[s].sr_par)/(1.0-spp[s].sr_par),-0.75)*0.25;
		} else { // invert for negatively weighted "get fast rid of almost fraction x", MCZ done 10.12.2009		  
		  if (repr_lvls[s]<=spp[s].sr_par) // enough has been removed
		    sp_derivative[s]=0.0;
		  else // function inversion to get rid of stuffitz fast
		    sp_derivative[s] = -4.0*z_pow((repr_lvls[s]-spp[s].sr_par)/(1.0-spp[s].sr_par),3.0);
		}
	      } else if (removal_rule==4) {
		// xxxrr4, MCZ derivative NOT used because rr4 uses difference between points
	      }	      
	    }
	  }
	}
	
	warp_ll=0;
	int warp_rnds(min(warp_factor, ecnt/100)); // xxxTBF fix
	/* ************************************************* TODO TODO *********************************************** */
	////////////////// warp_rnds += CCSP_last_failed_removals;

	if (warp_rnds<1)
		warp_rnds=1;

	{
		//create thread pool
		queue.nextIteration(
					val, rpos, minval, rx, ry, max_nbc, min_repr, min_smult, ranval, BL_delta_val, warp_rnds);
		int threads(queue.size());
		int task_size(ecnt / threads);
		if(ecnt % threads > 0)
			++task_size;
		for(int start_pos(0); start_pos < ecnt; start_pos += task_size) {
			int end_pos = min(start_pos + task_size,ecnt);
			queue.next(start_pos, end_pos);
		}
		queue.completeIteration();
	}

	if (1==warp_rnds) {  // this is also needed by the "re-load mode"
		wrx[warp_ll]  = rx; // what the fook is warp_ll, warp_rnds; warp_ll is 0.
		wry[warp_ll]  = ry;
		warp_ll++;
	}

	val = min_repr;

	if ((warp_rnds==1) && (rx==-1))
		return false;
	else
		return true;
}

int Get_smallest_mask_in_edge_list()
{
  int pos, mlvl;

  current_min_mask_level = std::numeric_limits<int>::max();
  for(pos=0; pos<ecnt; pos++) {
    mlvl = static_cast<int>(maskmap.m[eyl[pos]][exl[pos]]);
    if (mlvl<0)
      mlvl=0;
    if (mlvl<current_min_mask_level) {
      current_min_mask_level=mlvl;
      if (current_min_mask_level==sliding_initial_min_mask_level)  // no need to scan any more!
	break;
    }
  }

  // reset (increase) initial lvl once a lvl has been completely removed:
  if (current_min_mask_level > sliding_initial_min_mask_level)
    sliding_initial_min_mask_level = current_min_mask_level;
  
  return 0;
}

// @param t1 Start time of ranking process (after any initialization)
void Remove_and_rank_sites(time_t t1)
{
	int   rx, ry, rnd, rpos;
	int   w, wpos, wlpos;
	unsigned int poolsize;
	float val;

	if (run_mode==2)
		warp_factor=1;

	rnd = 0;
	c_pos = 0;
	loaded_cnt = 0;

	if (use_mask) {
	  bool ignore_hier_mask_in_reload_mode = false;
	  if (ignore_hier_mask_in_reload_mode && 2==run_mode) {
	    Form1->Memo1->Lines->Add("Running in re-load (-l) mode and because of that ignoring the (hierarchical) removal mask."
				     +IntToStr(current_min_mask_level));
	  } else {
	    Get_smallest_mask_in_edge_list();
	    Form1->Memo1->Lines->Add("Confirmed use of (hierarchical) removal mask. Initial smallest removal level: "
				     +IntToStr(current_min_mask_level));
	  }
	  sliding_initial_min_mask_level = current_min_mask_level;
	}


	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("The following lines provide information on how the ranking is done by iterative removal of cells, with 6 columns per line:");
	Form1->Memo1->Lines->Add("Cells removed, proportion of landscape remaining, cost of remaining landscape, edge list length, smallest remaining proportion in all biodiversity features (e.g. species), time elapsed (s)");
	Form1->Memo1->Lines->Add("----------------------------------------------------------------------------------------------------------------------------");


	print_rank_progress_line(0, 1.0f, cost_in_area, ecnt, min_repr, time(NULL)-t1);

	//create thread pool
	//	ThreadQueue queue(boost::thread::hardware_concurrency());
	if (use_tree_conn || ADM_set.use_ADMUs)
	  poolsize = 1;
	else
	  poolsize= VCLCommandLine::hardware_or_opt_concurrency();

	ThreadQueue queue(poolsize); // xxxP3

	bool push = false;    // Do remove cells that go over the removal_threshold (but not corridor critical cells)
	bool double_push = false;  // Do remove anything in the warp block


	float BL_delta_val[5];
	//if (BLP > .0f) { // && CCSP > .0f && Grid_CCSP::Fast_Warp == cc_accounting->calc_mode()) {
	if (BLP > .0f && CCSP > .0f  && corr_ongoing /*&& Grid_CCSP::Fast_Warp == cc_accounting->calc_mode()*/ ) {
	  fix_BL_delta(BL_delta_val);
	}

	time_t t_ini = time(NULL);	
	while(get_next_to_remove(val, queue))
	{

		if (terminate_flag==true)
		{
			if (Application->MessageBox("Terminate has been pressed, abort computations?", "Abort ZIG?", MB_YESNO)==ID_YES)
				return;
			terminate_flag=false;
		}

		if (use_PLULA) {

		  int cell_warp_ll = warp_ll;   // normal "cell removal" warp list length
		  int x = wrx[0];
		  int y = wry[0];
		  size_t plu_n = PLL[y][x];

		  int warp_plu_fill_cnt = 0;
		  if (cell_warp_ll>0) {
		    Form1->Memo1->Lines->Add(" Removing PLU, id: "+IntToStr(PLvec[plu_n].num)+", cells: "+
					     IntToStr(PLvec[plu_n].el_cnt)+", x: "+IntToStr(x)+
					     ", y: "+IntToStr(y));
		    // copy original list of cells to remove before it's too late
		    if (cell_warp_ll > PLvec[plu_n].el_cnt)
		      copy_plu_warp_list(cell_warp_ll);  // ->plu_copy_wrx / _wry
		    warp_plu_fill_cnt = fill_warp_list_plu(plu_n, 0);
		  }
		  if (cell_warp_ll > PLvec[plu_n].el_cnt)
		    warp_plu_fill_cnt = fill_as_much_as_poss_warp_list_plu(warp_plu_fill_cnt, plu_n, cell_warp_ll);

		  warp_ll = warp_plu_fill_cnt;
		  if (warp_factor>1)
		    Form1->Memo1->Lines->Add(" (PLU effective warp factor:  "+IntToStr(warp_ll)+")");

		  if (warp_ll==0) {
		    Form1->Memo1->Lines->Add("warp 0 error - while filling in warp list for PLU"); graceful_exit();}
		  if (warp_ll > wrx.size())
		    Form1->Memo1->Lines->Add("**** WARNING/ERROR: max grid cells for planning units exceeded, there are: "+IntToStr(wrx.size())+" ****");
		}


		int failed_removals = 0;
		int total_in_risk = 0, total_no_risk = 0;
		// TODO:
		float warp_max_val = -std::numeric_limits<float>::max();
		for(w=0; w < warp_ll; w++) {
		  if(wval[w]>warp_max_val) {
		    warp_max_val = wval[w];
		  }
		}
		// CCSP_last_removal_threshold = warp_max_val;  // <- for use outside this function (in update_wlist) = ABANDONED!!!
		//std::cerr << "for... ";

		if (CCSP > .0f && corr_ongoing &&
		    (Grid_CCSP::Fast_Warp == cc_accounting->calc_mode() ||
		     (CCSP_slow_as_fast && Grid_CCSP::Slow == cc_accounting->calc_mode()))
		    )
		  cc_accounting->init_next_warp_fast(warp_max_val);

		for(w=0; w < warp_ll; w++)
		{
			rx          = wrx[w];
			ry          = wry[w];
			rpos=-1;

#ifdef USE_MAT_EDGE_LIST_POS
			rpos = mat_edge_list_pos[ry][rx];
#else
			for(int w2=0; w2<ecnt; w2++)
			  if ((exl[w2]==rx) && (eyl[w2]==ry))
			    rpos = w2;
#endif

			if (rpos==-1)
			{
				Form1->Memo1->Lines->Add("Cell to remove not in the edge cell list, rx="+IntToStr(rx)+" ry="+IntToStr(ry)+
							 ", edge[ry][rx]: "+IntToStr(edge[ry][rx]));
				//if (Rmax[ry][rx]<0)
				if (-1 == status[ry][rx])
				  Form1->Memo1->Lines->Add("it has no data");
				else
				  Form1->Memo1->Lines->Add("but it has data");
				Form1->Memo1->Lines->Add("It has status = "+IntToStr(status[ry][rx]));
				Form1->Memo1->Lines->Add("It has edge = "+IntToStr(status[ry][rx]));
				if (NULL != PLL) {
				  Form1->Memo1->Lines->Add("Planning units in use...");
				  if (PLL[ry][rx]<0)
				    Form1->Memo1->Lines->Add("it has no valid planning unit ID: "+IntToStr(PLL[ry][rx]));
				  else
				    Form1->Memo1->Lines->Add("but it has valid planning unit ID: "+IntToStr(PLvec[PLL[rx][ry]].num));
				} else {
				  Form1->Memo1->Lines->Add("Planning units not in use...");
				}
				Form1->Memo1->Lines->Add("This implies some sort of inconsistency in the input raster layers");
				Form1->Memo1->Lines->Add("Normally this is  an unrecoverable error and I am probably going to stop or crash soon...");
				continue;
			}
			
			bool warp_fix = true;
			// safe to remove [rx,ry] ???
			//if (warp_fix && 0 == w)
			//  warp_nondirty_count++;

			if (CCSP > .0f && corr_ongoing && CCSP_write_log_file && 1 == current) {
			  String fn = String(outgridfn + ".debug_wlist.txt");
			  std::cerr << "Opening: " << fn.toStdString() << std::endl;
			  wlist_file = fopen(fn.toStdString().c_str(), "w");
			  fprintf(wlist_file, "# Round, w, warp_ll, rx, ry, lbl, area[lbl], wrscr[lbl], last_lbln, last_last_lbnln, uneroded_last_lbln, wval[w], tmpval, warp_max_val, \n");
			}


			//std::cerr << "a";
			if (CCSP > .0f && corr_ongoing && 
			    (!CCSP_slow_as_fast && 
			     Grid_CCSP::Slow == cc_accounting->calc_mode())
			    && warp_fix /*&& get_nb_count(status, rx, ry) >1*/ /*&& w > 0*/) {
			  // this is a minimal version of the value calculation done in SearchIteration.search()

			  
			  // INIT is strictly needed before any calculation of penalties (with explicit calc_ccsp_as_if_removed())
			  cc_accounting->slow_init_removal(rx, ry);

			  // (Removed redundancy - see ccscp.cpp)
			  // bool risk = cc_accounting->estimate_split_risk(status, rx, ry);
			  bool risk = true;

			  int nb_count = 0;
			  int bl_inc = 0;
			  if (BLP>0.0f) {
			    nb_count = get_nb_count(status, rx, ry);
			    // This is done in Remove_site() !!!
			    //int bl_inc = 0; // = BL_delta_val[nb_count];
			    bl_inc = BL_delta_val[nb_count];
			  }

			  /*
			  bool fix_bl = false;
			  if (fix_bl) {
			    float BL_delta_val[5];
			    if (BLP>0.0f) {
			      fix_BL_delta(BL_delta_val);
			      int nb_count = get_nb_count(status, rx, ry);
			      bl_inc = BL_delta_val[nb_count];
			    }
			  }			  */
			  /*
			  if (wdirty[ry][rx]>0) {
			    failed_removals++;
			    continue;
			  }			  */

			  float tmpval;
			  if (risk) {
			    total_in_risk++;
			    // This does apply the CCSP
			    tmpval = get_delta_value_plus_penalties(rx, ry, false, bl_inc, true, warp_max_val);
			    if (CCSP_write_log_file)
			      fprintf(wlist_file, "# *** CCL RECALC:\n");
			    warp_recalculations++;
			  } else {
			    total_no_risk++;
			    tmpval = get_delta_value_plus_penalties(rx, ry, false, bl_inc, false, warp_max_val);
			    if (CCSP_write_log_file)
			      fprintf(wlist_file, "# NO RISK:\n");
			    //wval[w] = tmpval;
			  }
			  
			  //bool split_happened = risk && cc_accounting->last_lbln > cc_accounting->last_last_lbln;
			  //bool split_risk = risk && cc_accounting->is_cached_penalty(rx, ry);
			  bool split_risk = risk; // || cc_accounting->is_cached_penalty(rx, ry);

			  if (CCSP_write_log_file && push)
			    fprintf(wlist_file, "# *******************^ PUSHING **************************************\n");
			  if (CCSP_write_log_file && double_push)
			    fprintf(wlist_file, "# *******************^ DOUBLE PUSHING **************************************\n");

			  //warp_max_val = wval[warp_ll-1];
			  bool blp_shift = (BLP>0.0f && BL_delta_val[nb_count]>0);
			  if (!(blp_shift) &&  /* exclude BLP-caused shift of delta */
			  //if (false && 
			      ( /* (wdirty[ry][rx] > 0) || */   // this first line should make no difference
			       (tmpval >  warp_max_val
				/* && (risk && cc_accounting->last_lbln > cc_accounting->last_last_lbln) */
				&& warp_ll > 5
				)
			       && 
			       !(push && !(split_risk /*split_happened*/) )   // push applies only to non-corridor cells? ---> may not finish!
			       &&
			       !double_push
				)
			      ) {

			    //fprintf(wlist_file, "C BEF %7d, %4d, %4d, %4d, %4d, %5d, %7d, %10.10f, %5d, %5d, %10.10f, %10.10f, %10.10f, \n", rnd, w, warp_ll, rx, ry, lbl, cc_accounting->areas[lbl], cc_accounting->wrscr[lbl], cc_accounting->last_lbln, cc_accounting->last_last_lbln, wval[w], tmpval, warp_max_val);

			    cc_accounting->slow_cancel_removal(rx, ry);
			    // *****TODO:
			    //wdirty[ry][rx] += 20; //20;
			    //wdirty[ry][rx] = max(1.10 * ((float)removed/(float)nonm1), 0.99);
			    warp_dirty_count++;
			    //Form1->Memo1->Lines->Add("Warp element dirty in round " + IntToStr(rnd) + ": " + IntToStr(w+1) + "/" + IntToStr(warp_ll) + " - " + IntToStr(rx) + "," + IntToStr(ry));			    
			    
			    failed_removals++;
			    
			    
			    //fprintf(wlist_file, "%7d, %4d, %4d, %4d, %4d, %5d, %7d, %7d, %5d, %5d, %10.10f, %5.5f, %10.10f, %10.10f, \n", rnd, w, warp_ll, rx, ry, cc_accounting->get_cc_label(rx, ry), cc_accounting->last_size, cc_accounting->safe_size, cc_accounting->last_lbln, cc_accounting->last_last_lbln, cc_accounting->wrscr[last_lbln], wval[w], tmpval, warp_max_val);
			    if (CCSP_thickness <=1) {
			      int lbl = cc_accounting->get_cc_label(rx, ry);
			      if (CCSP_write_log_file)
				fprintf(wlist_file, "C AFT %7d, %4d, %4d, %4d, %4d, %5d, %7d, %10.10f, %5d, %5d, %10.10f, %10.10f, %10.10f, \n", rnd, w, warp_ll, rx, ry, lbl, cc_accounting->areas[lbl], cc_accounting->wrscr[lbl], cc_accounting->last_lbln, cc_accounting->last_last_lbln, wval[w], tmpval, warp_max_val);
			    } else {
			      int lbl = cc_accounting->uneroded_lbl_buffer[ry][rx];
			      if (CCSP_write_log_file)
				fprintf(wlist_file, "C AFT %7d, %4d, %4d, %4d, %4d, %5d, %7d, %10.10f, %5d, %5d, %10.10f, %10.10f, %10.10f, \n", rnd, w, warp_ll, rx, ry, lbl, cc_accounting->uneroded_areas[lbl], cc_accounting->uneroded_wrscr[lbl], cc_accounting->last_lbln, cc_accounting->last_last_lbln, wval[w], tmpval, warp_max_val);			    
			    }
			    
			    if (CCSP_write_log_file)
			      fprintf(wlist_file, "# -------------------^ REMOVAL CANCELED\n");
			    continue;
			  } else {
			    
			    cc_accounting->slow_accept_removal(rx, ry);
			    warp_nondirty_count++;			    

			    if (CCSP_thickness <=1) {
			      int lbl = cc_accounting->get_cc_label(rx, ry);
			      if (CCSP_write_log_file)
				fprintf(wlist_file, "AFT %7d, %4d, %4d, %4d, %4d, %5d, %7d, %10.10f, %5d, %5d, %10.10f, %10.10f, %10.10f, \n", rnd, w, warp_ll, rx, ry, lbl, cc_accounting->areas[lbl], cc_accounting->wrscr[lbl], cc_accounting->last_lbln, cc_accounting->last_last_lbln, wval[w], tmpval, warp_max_val);
			    } else {
			      int lbl = cc_accounting->uneroded_lbl_buffer[ry][rx];
			      if (CCSP_write_log_file)
				fprintf(wlist_file, "AFT %7d, %4d, %4d, %4d, %4d, %5d, %7d, %10.10f, %5d, %5d, %10.10f, %10.10f, %10.10f, \n", rnd, w, warp_ll, rx, ry, lbl, cc_accounting->uneroded_areas[lbl], cc_accounting->uneroded_wrscr[lbl], cc_accounting->last_lbln, cc_accounting->last_last_lbln, wval[w], tmpval, warp_max_val);
			    }

			  }

			} else if ( CCSP > .0f && corr_ongoing && 
				    (Grid_CCSP::Fast_Warp == cc_accounting->calc_mode() ||
				     (CCSP_slow_as_fast && Grid_CCSP::Slow == cc_accounting->calc_mode()))
				    ) {

			  int nb_count = 0;
			  if (BLP>0.0f)
			    nb_count = get_nb_count(status, rx, ry);

			  float updated_val = get_delta_value_plus_penalties_fast(rx, ry, false, BL_delta_val[nb_count], 0);

			  if (CCSP_thickness <=1) {			    
			    int lbl = cc_accounting->get_cc_label(rx, ry);
			    if (CCSP_write_log_file)
			      fprintf(wlist_file, "%7d, %4d, %4d, %4d, %4d, %5d, %7d, %10.10f, %5d, %5d, %5d, %10.10f, %10.10f, %10.10f, \n", rnd, w, warp_ll, rx, ry, lbl, cc_accounting->areas[lbl], cc_accounting->wrscr[lbl], cc_accounting->last_lbln, cc_accounting->last_last_lbln, cc_accounting->uneroded_last_lbln, wval[w], updated_val, warp_max_val);
			  } else {
			    int lbl = cc_accounting->uneroded_lbl_buffer[ry][rx];
			    if (CCSP_write_log_file)
			      fprintf(wlist_file, "%7d, %4d, %4d, %4d, %4d, %5d, %7d, %10.10f, %5d, %5d, %5d, %10.10f, %10.10f, %10.10f, \n", rnd, w, warp_ll, rx, ry, lbl, cc_accounting->uneroded_areas[lbl], cc_accounting->uneroded_wrscr[lbl], cc_accounting->last_lbln, cc_accounting->last_last_lbln, cc_accounting->uneroded_last_lbln, wval[w], updated_val, warp_max_val);
			  }
			  
			  // detect (potentially big) shift of marginall loss because of BLP - and ignore!
			  bool blp_shift = (BLP>0.0f && BL_delta_val[nb_count]>0);
			  if ( !(blp_shift) &&  // not BLP-penalized (increased delta)
			      updated_val > warp_max_val && warp_ll > 2) { // implies updated_val > wval[w] 
			    failed_removals++;
			    continue;
			  } else {
			    cc_accounting->accept_fast_removal(rx, ry);
			  }
			}

			/*
			if (wdirty[ry][rx]) {
			  failed_removals--;
			  sol[ry][rx]     = 1; //current/float(nonm1*(1.0f-rem_level));
			  sol_val[ry][rx] = 1.0f-val; // xxxMCZ error
			} else {
			*/

			// --> Remove_site(rx, ry)
			Remove_site_and_record_curves(rx, ry, rpos, val, rnd, t_ini);
			rnd++;

			/*
			// warp list dirtyness report
			std::cerr << "warp list dirtyness : ";
			for(w=0; w<warp_ll; w++) {
			  //if (wdirty[w]) std::cerr << w << " ";
			  //wdirty[w] = false;			  
			} std::cerr << std::endl;
			//Form1->Memo1->Lines->Add("Warp dirtiness ratio in round " + IntToStr(rnd) + ": " + FloatToStr(float(dirty_count)/float(warp_ll)) + "! ");
			*/
		}

		CCSP_prev_prev_removal_threshold = CCSP_prev_removal_threshold;
		CCSP_prev_removal_threshold = warp_max_val;

		// note warp_max_val could be == CCSP_prev_removal_threshold => CCSP_ma_removal_rate = 0
		CCSP_ma_removal_rate = 0.95 * CCSP_ma_removal_rate + 0.05 * (warp_max_val-CCSP_prev_removal_threshold);
		/* This is handled in the penalty calculations, rather than fixing it here
		if (CCSP_ma_removal_rate < 0.000001f)
		  CCSP_ma_removal_rate = 0.000001f;
		*/

		// *****
		CCSP_verbose = 0;
		if (CCSP > .0f && corr_ongoing && CCSP_verbose > 0)
		  Form1->Memo1->Lines->Add("     >>> rnd: " + IntToStr(rnd) + ", warp_ll: " + IntToStr(warp_ll) + 
					 ", failed_removals: " + IntToStr(failed_removals) + 
					 ", RISK: " + IntToStr(total_in_risk) + ",NORSK: " + IntToStr(total_no_risk) + 
					 ", cached: " + IntToStr(cc_accounting->num_cached) + ", ever cached: " + IntToStr(cc_accounting->num_ever_cached) 
					 /* + ",  current-1: " + IntToStr(current - 1)*/ );
		CCSP_last_failed_removals = failed_removals;

		/*
		if (failed_removals > warp_ll / 0.65f)
		  cc_accounting->increase_cache_factor();
		else if (failed_removals < warp_ll / 0.15f)
		  cc_accounting->decrease_cache_factor();
		*/
		  
		if (failed_removals >= warp_ll){
		  if (push) {
		    double_push = true;
		    failed_removals = 0;
		    Form1->Memo1->Lines->Add("DOUBLE PUSH, rnd: " + IntToStr(rnd));
		  } else {
		    push = true;
		    double_push = false;
		    Form1->Memo1->Lines->Add("PUSH, rnd: " + IntToStr(rnd));
		    /*
		      for(int y=0;y<yd;y++) {
		      for(int x=0;x<xd;x++) {
		      wdirty[y][x] = 0;
		      }
		      }
		    */
		  }
		} else {
		  // drop 2 pushness levels
		  push = double_push = false;
		}

		double_push = false;


		// end of warp list processing
		if (use_mask && 2!=run_mode)   // ignore if in reload mode -> see Get_sorted_cells_in_map()
		  Get_smallest_mask_in_edge_list();
	}

	// Make sure that a last line is printed for fraction remaining 0%
	time_t elapsed_s = time(NULL) - t_ini;
	print_rank_progress_line(current-1, 0.0f, 0.0f, ecnt, .0f, elapsed_s);
}

float
Get_current_reprs(double* repr, double* repr2, bool op)
{
	int   s,x,y, pos;
	char  *sr;
	float val, *mr;

	if (op)	{
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("Biodiversity features performance levels check. Proportions remaining:");
	  Form1->Memo1->Lines->Add("---------");
	}

	for(s=0; s<map_cnt; s++) {
	  repr[s] = 0.0f;
	}

	if (use_SSI) {
	  for(s=0; s<SSI_spp_cnt; s++) {
	    repr2[s] = 0.0f;
	  }
	}

	for(y=0; y<yd; y++) {
	  sr = &status[y][0];
	  //      mr = &obsmap[s].m[y][0];
	  for(x=0; x<xd; x++) {
	    if (sr[x]>0) {
	      // rowp = &vmat[y][x][0];   // COMPACT_VMAT
	      const Biodiv_Features_Occur_Container& rowp = vmat[y][x];
	      // for(s=0;s<map_cnt;s++)
	      for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s)) {
		float rowp_s_val = rowp[s];
		//if (rowp[s]>0.0f)
		//	repr[s]+=rowp[s];  // xxxW
		if (rowp_s_val>0.0f)
		  repr[s]+=rowp_s_val;  // xxxW
#ifdef ZCORE_DEBUG
		if (isnan(rowp_s_val)) // NaNs show up here!
		  Form1->Memo1->Lines->Add("rowp[s] NaN: " + IntToStr(s) + ", x: " + IntToStr(x) + ", y: " + IntToStr(y));
#endif
	      }

	      if (use_SSI) {
		for(s=0; s < SSIxyCnt[y][x]; s++) {
		  pos = SSIxyPos[y][x]+s;
		  repr2[SSI_list[pos].spnum] += SSI_list[pos].val;  // xxx SSI
		}
	      }
	    }
	  }
	}

	if (op)
	{
	  for(s=0;s<map_cnt;s++)
	    Form1->Memo1->Lines->Add("Feature "+IntToStr(s+1)+": "+FloatToStrF(repr[s], ffFixed, 7, 4));
	  if (use_SSI) {
	    Form1->Memo1->Lines->Add("-------------");
	    for(s=0;s<SSI_spp_cnt;s++)
	      Form1->Memo1->Lines->Add("SSI feature: "+IntToStr(s+1)+": "+FloatToStrF(repr2[s], ffFixed, 7, 4));
	  }
	}

	min_repr=std::numeric_limits<float>::max();
	max_repr=-1;
	float ave = 0.0f;
	float wave = 0.0f;
	float w_accum = 0.0f;
	for(s=0;s<map_cnt;s++) {
	  if (repr[s]<min_repr)
	    min_repr = repr[s];
	  if (repr[s]>max_repr)
	    max_repr = repr[s];
	  
	  ave += repr[s];
	  wave += spp[s].weight * repr[s];
	  w_accum += spp[s].weight;
	}
	wave /= w_accum;
	ave /= map_cnt;

	if (op) {
	  Form1->Memo1->Lines->Add("---------");
	  Form1->Memo1->Lines->Add("Minimum proportion remaining: "+
				   FloatToStrF(min_repr, ffFixed, 7, 4)+
				   ", average: "+FloatToStrF(ave, ffFixed, 7, 4)+
				   ", weighted average: "+FloatToStrF(wave, ffFixed, 7, 4)+
				   ", maximum: "+FloatToStrF(max_repr, ffFixed, 7, 4));
	  Form1->Memo1->Lines->Add("");
	}
	return min_repr;
}

// Assumes all inputs have been loaded (and that loaded1.m is ready)
bool
check_consistency_reload_mode()
{
  bool ok;
  int MAX_ERR_MSG = 10;

  // A) make sure that all the cells in the old solution have at least one feature occurrency (are not missing in the new vmat)
  // load_order contains the x,y coordinates of cells in the original solution
  int err1_cnt = 0;  
  for (size_t i=0; i<nonm1; i++) {
    int rx = load_order[i].x;
    int ry = load_order[i].y;
    if (!vmat[ry][rx] && err1_cnt < MAX_ERR_MSG) {
      err1_cnt++;
      if (err1_cnt < MAX_ERR_MSG) {
	Form1->Memo1->Lines->Add("Error in reload mode, cell at x:"+IntToStr(rx)+", y:"+IntToStr(ry)+
				 " is part of the original (reloaded) solution but it has no occurrency of any of the features in the new list of features");
      } else if (MAX_ERR_MSG==err1_cnt) {
	Form1->Memo1->Lines->Add("Omitting next errors to avoid a flood of messages...");
      }
    }
  }

  // B) make sure that all the cells with occurrencies (non-missing cells in the new vmat) are non-missing in the old solution
  // errors in the opposite direction to err1_cnt
  int err2_cnt = 0;
  for (int y=0; y<yd; y++){
    for (int x=0; x<xd; x++){      
      if (vmat[y][x] && loaded1.m[y][x]<0) {
	err2_cnt++;
	if (err2_cnt < MAX_ERR_MSG) { 
	  Form1->Memo1->Lines->Add("Error in reload mode, cell at x:"+IntToStr(x)+", y:"+IntToStr(y)+
				   " has occurrencies of features but it is missing in the original (reloaded) solution");
	} else if (MAX_ERR_MSG==err2_cnt) {
	  Form1->Memo1->Lines->Add("Omitting next errors to avoid a flood of messages...");
	}
      }
    }
  }

  if (err1_cnt>0 || err2_cnt>0) {
    Form1->Memo1->Lines->Add("ERROR: Inconsistencies have been detected in reload mode, please check and fix the inputs.");
    if (err1_cnt>0) {
      Form1->Memo1->Lines->Add("ERROR: "+IntToStr(err1_cnt)+ 
			       " cells of the original (reloaded) solution have no occurrencies of any of the new features (the first "+
			       IntToStr(MAX_ERR_MSG)+" ones were shown above.");
    }
    if (err2_cnt>0) {
      Form1->Memo1->Lines->Add("ERROR: "+IntToStr(err2_cnt)+ 
			       " cells have occurrencies of the new features but are missing in the original solution (the first "+
			       IntToStr(MAX_ERR_MSG)+" ones were shown above.");
    }
    // be liberal in what you accept?
    return false;
  }
  
  return true;
}

bool write_output_files_at_the_end()
{
  Form1->Memo1->Lines->Add("============================================================"); 
  Form1->Memo1->Lines->Add("* Writing final output files..."); 
  Form1->Memo1->Lines->Add("Writing file of feature information: " +
			   features_info_fn);
  Form1->Memo1->Lines->Add("Writing file of performance (representation, coverage, etc.) curves: " +
			   curvesfn);
  output_features_info_n_curves();

  if (show_images && !mem_save_mode) {
    Form1->Memo1->Lines->Add("Writing output ranking as an image...");
    save_jpg();
  } else {
    Form1->Memo1->Lines->Add(" Note: not writing output ranking as an image");
  }

  Form1->Memo1->Lines->Add("Writing ranking as a GIS raster map...");
  output_grid_rank(0, xd, yd);

  if (glob_set_output_prop_rank_map) {
    Form1->Memo1->Lines->Add("Writing 'proportional loss map' as a GIS raster map...");
    output_grid_prop_rank(0, xd, yd);
  } else {
    Form1->Memo1->Lines->Add(" Note: not writing 'proportional loss map'");
  }

  //  output_cost_grid(0);

  if (ADM_set.use_ADMUs) {
    Form1->Memo1->Lines->Add("Preparing redistributed rank map for administrative units...");
    calculate_ADMU_redistributed_rank(0);
    Form1->Memo1->Lines->Add("Writing redistributed rank map for administrative units: " +
			     admu_redist_outgridfn);
    output_ADMU_redistributed_rank(0);
  }
  
  if (ADM_set.row_count_for_per_admu_curves > 0) {
    ADMUs_output_per_admu_curves_iter(1.0);
  }
  
  if (use_SSI) {
    Form1->Memo1->Lines->Add("Writing file of SSI feature information and performance (representation, coverage, etc.) curves: " +
			     SSI_features_info_fn);
    output_SSI_features_info_n_curves();
  }
  
  if (use_groups) {
    output_grp_curves(glbl_groups_info);
    Form1->Memo1->Lines->Add("Writing group curves file: " + grp_curves_fn);
  }
  
  if (ADM_set.use_ADMUs && use_groups)
    ADMUs_output_per_admu_grp_curves_iter(1.0, glbl_groups_info);
  
  post_process(PPA_fname);
  Form1->Memo1->Lines->Add("============================================================"); 

  return true;
}

bool bat_run2(ZConf const& conf)
{
	time_t t1, t2;
	int show_scale;

	Fix_output_file_names(conf);

	DecimalSeparator='.';
	run_finished = false;

	Form1->Visible=true;

	Screen->Cursor=crHourGlass;

	for(int loop=0; loop<MAX_SPP_COUNT; loop++)
		IG_set.IG_spw[loop]=1.0f;

	Form1->Caption="ZIG: Loading data files, please wait";

	t1=time(NULL);

	if (!data_input_v2())
	{
		free_data();
		//MessageBox(0, "Could not successfully process input files.", "Data input error", MB_OK);
		Form1->Memo1->Lines->Add("Data input error. Could not successfully process input files.");
		//       Delay(5);
		Screen->Cursor=crDefault;
		return false;
	}

	//  obsmap[0].draw_on_bitmap_smooth(Form1->Image1->Picture->Bitmap, 1.0f, 0.0f); // xxx2
	//  if (yd>650)
	Form1->Height=1024;
	//  else
	//    Form1->Height=yd+150;

	//  Form1->Memo1->Lines->Add("yd="+IntToStr(yd));
	//  Form1->Memo1->Lines->Add("FH="+IntToStr(Form1->Height));
	Form1->Caption=spp[0].fname;
	Form1->Image1->Repaint();
	AnnotForm->Edit6->Text = IntToStr(Form1->Image1->Width -200);
	AnnotForm->Edit7->Text = IntToStr(Form1->Image1->Height- 130);
	AnnotForm->Edit9->Text = IntToStr(Form1->Image1->Width -210);
	AnnotForm->Edit10->Text= IntToStr(Form1->Image1->Height- 70);

	Form1->Memo1->Lines->Add("Calculating richness across biodiversity features...");
	time_t t1_wrich, t2_wrich;
	t1_wrich = time(NULL);
	Calc_richness_et_al_matrixes();
	t2_wrich = time(NULL);
	Form1->Memo1->Lines->Add("Done in " +IntToStr(t2_wrich-t1_wrich)+  " seconds. Time now: " + zCurrentTime());

	// nonm1 (site count) calculated, init interprocess map
	//interprocessInitMap(nonm1, tstmap.cols(), tstmap.rows()); <- NEVER use tstmap here
	interprocessInitMap(nonm1, xd, yd);

	Fix_nbm();

	initialize_status();
	
	// Now output wrscr asap (if needed), so vmat can be freed soon (if possible)

	if (glob_set_output_wrscr_map) {
	  Form1->Memo1->Lines->Add("Saving raster map of weighted range size corrected richness");	
	  output_grid_wrscr(0, xd, yd);
	} else {
	  Form1->Memo1->Lines->Add("Note: not saving raster map of weighted range size corrected richness");	
	}
	if (!corr_set.use_corr)
	  free_wrscr_Rsum();

	Initialize_remove(); // after this, edge[][] and status[][] are initd -> Rmax not needed anymore

	Fix_BQP_nbms();

	if (use_cost) {
	  get_cost_info();
	  Form1->Chart2->BottomAxis->SetMinMax(0, cost_in_area);		
	} else {
	  Form1->Chart2->BottomAxis->SetMinMax(0, nonm1);
	}

	// per ADMU costs (only required for certain outputs
	if(use_cost && ADM_set.use_ADMUs && 
	   (use_groups || ADM_set.row_count_for_per_admu_curves > 0) ) {
	  // need all the info
	  get_per_ADMU_area_and_cost_info_and_non_missing();
	} else if (ADM_set.use_ADMUs) {
	  // need just # non missing for the redistributed rank
	  get_per_ADMU_area_and_non_missing();
	}

	if (rem_by_rarity)
	  Get_sorted_cells_in_map(srtvec, Rmax);
	//else
	//  Get_sorted_cells_in_map(srtvec, Rsum);

//#if 0
	obsmap[0].draw_on_bitmap_smooth_from_mat(Form1->Image1->Picture->Bitmap,
						 Rmax, Rmax_max, 0.0f, true, 1);
	if (use_SSI)
		obsmap[0].show_spots(Form1->Image1->Picture->Bitmap,  SSIxyCnt, 2);
//#endif

	init_sol_and_sol_val(); // Transfer Rmax into sol[][] and sol_val[][] before freeing it
	// After Get_sorted_cells_in_map(), Rmax is not needed anymore, except if using a shorcut for CAZ (see SearchIteration::search(int, int))
	if ( !((!use_BQP) && (removal_rule==1) && (!neg_weights_used) && !ADM_set.use_ADMUs))
	  free_wrscr_Rmax();


	t2=time(NULL);

	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("Loaded data and initialized in "+IntToStr(t2-t1)+" seconds. Current time: " + zCurrentTime());
	Form1->Memo1->Lines->Add("Cells with data = "+IntToStr(nonm1)+"; locations with missing values = "+IntToStr(m1s));
	Form1->Edit38->Text=IntToStr(nonm1);
	Form1->Edit35->Text=IntToStr(m1s);
	Form1->Memo1->Lines->Add("");

	Form1->Memo1->Lines->Add("---------------------======********** RANKING STARTS HERE ***********======-----------------------");
	Form1->Caption="ZIG: Iterative cell removal in progress";
	Form1->PageControl1->ActivePageIndex=0;

	if (run_mode==2)
	{
		load_order = new struct sort_xy[nonm1+1];
		if (!load_order)
		{
			ShowMessage("Out of memory for allocating load order vector, quitting");
			graceful_exit();
		}
		
		Get_sorted_cells_in_map(load_order, loaded1.m, use_mask);
		Form1->PageControl1->ActivePageIndex=3;

		// check consistency of load_order/loaded1 with vmat
		if (!check_consistency_reload_mode()) {
		  Form1->Memo1->Lines->Add("FATAL ERROR: reload (-l) mode cannot be run because of inconsistencies in the input files. The extent of the original solution does not match the extent of the features in the current settings (see messages above).\n");
		  graceful_exit();
		  return false;
		}
	}

	if (autoclose)
		Form1->PageControl1->ActivePageIndex=3;

	if (2!=run_mode)
	  Remove_bad(rem_level);  // this is an old leftover related to the old "initial removal"

	// Important for large maps
	// srtvec is only used for Get_sorted_cells_in_map() (fills it in), and Remove_bad() (if using initial removal level)
	if (srtvec)
	  delete [] srtvec;
	srtvec = NULL;
	//
	
	min_repr = Get_current_reprs(repr_after_rem, SSI_repr_lvls, true);

	for(int loop=0; loop<map_cnt;loop++)
		repr_lvls[loop]=repr_after_rem[loop];

	curve_store_interval=(int)(nonm1*(1.0f-rem_level))/1000;
	if (curve_store_interval<=0)
		curve_store_interval=1;

	ADM_store_interval=(int)(0.01*nonm1*(1.0f-rem_level));
	if (ADM_store_interval<=0)
		ADM_store_interval=1;

	if (ADM_set.row_count_for_per_admu_curves > 0) {
	  ADM_per_admu_store_interval = (int)(1.0/ADM_set.row_count_for_per_admu_curves
					      *nonm1*(1.0f-rem_level));
	  if (ADM_store_interval<=0)
	    ADM_store_interval=1;
	}

	//Application->Run();
	//free_data();
	//return false;

	if (use_PLULA)
	{
		bool ok = check_PLULA();
		if (!ok)
		  graceful_exit();


		Form1->Memo1->Lines->Add("Checking distributions of features in planning units, and allocating vectors of distributions within planning units. This may take a while...");
		if (get_PLULA_sp_data()) {
		  Form1->Memo1->Lines->Add("PLULA data vectors allocated with success for all PLULAs.");
		} else {
		  Form1->Memo1->Lines->Add("MEMORY LOW WARNING: Failure to allocate PLULA data vectors.");
		  Form1->Memo1->Lines->Add("More available memory would be needed. Giving up before it is too late...");
		  graceful_exit();
		}
	}


	t1=time(NULL);
	interprocessNotifyHardWorkBegin(); // *** hard work begins ***
	Remove_and_rank_sites(t1);
	interprocessNotifyHardWorkEnd(); // *** hard work ends ***
	Form1->Memo1->Lines->Add("----------------------------------------------------------------------------------------------------------------------------");
	Form1->ProgressBar1->Position = (int)((current+rem_level*nonm1)*100.0f/nonm1);
	Form1->Memo1->Lines->Add("Total count of cells removed = "+IntToStr(current - 1));
	//  sprintf(txt, "xd%i yd%i el%i non-1=%i",xd,yd,xd*yd,nonm1);
	//  Form1->Memo1->Lines->Add(txt);
	Form1->Memo1->Lines->Add("");


	if (terminate_flag)
	{
		free_data();
		return false;
	}

	t2=time(NULL);
	Form1->Memo1->Lines->Add("Done in "+IntToStr(t2-t1)+" seconds.");
	Form1->Memo1->Lines->Add("Found "+IntToStr(ties_cnt)+" ties.");

	Screen->Cursor=crDefault;
	//obsmap[0].draw_on_bitmap_smooth_from_mat(Form1->Image1->Picture->Bitmap,
	//										 sol, 1.0f, 0.0f, true, false);
	
	min_repr = Get_current_reprs(repr_lvls, SSI_repr_lvls, true);

	// Free edge list, etc buffers asap
	free_opt_structures();
	if (PPA_fname.isEmpty())
	  free_vmat();
	// After this there is additional memory need for buffers for the output layers

	write_output_files_at_the_end();

#if 1
	bool rem_ok;
	const size_t LINE_LEN = 512;
	char line[512];
	rem_ok=true;
	if (use_PLULA)
	{
		for(int loop=0; loop<PLcnt; loop++)
		{
			if (!PLvec[loop].removed)
			{
				Form1->Memo1->Lines->Add("not removed PLU "+IntToStr(loop));
				rem_ok=false;
			}
			snprintf(line, LINE_LEN, "%i  odc=%f dc=%f  ouc=%f uc=%f", PLvec[loop].num,
				PLvec[loop].orig_tree_conn_down,
				PLvec[loop].tree_conn_down, PLvec[loop].orig_tree_conn_up,
				PLvec[loop].tree_conn_up);
			//        Form1->Memo1->Lines->Add(line);
		}
		if (!rem_ok)
			Form1->Memo1->Lines->Add("Did not remove all planning units. Something may be wrong!");
		else
			Form1->Memo1->Lines->Add("Removed all planning units - OK.");
	}
#endif

	Form1->Memo1->Lines->Add("Finished spatial prioritization process. Freeing data structures in memory...");
	//	Form1->Caption = "ZIG2: - Done.";
	run_finished   = true;
	Form1->Button7->Enabled = true;

	Screen->Cursor=crDefault;
	if (!autoclose)
	{
		Form1->Button8->Enabled = true;
		Application->Run();
	}

	free_data();
	Form1->Memo1->Lines->Add("Finished at " + zCurrentTime());
	int elapsed = zElapsedTime();
	const uint64_t ms_per_h = 3600000;
	const uint64_t ms_per_day = ms_per_h * 24;
	Form1->Memo1->Lines->Add("Elapsed time: " + IntToStr(elapsed) + 
				 " ms == "+ FloatToStrF(float(elapsed)/ms_per_h, ffFixed, 3, 3) + " hours == " +
				 FloatToStrF(float(elapsed)/ms_per_day, ffFixed, 3, 3) + " days");
	Form1->Close();

	return true;
}

bool bat_run()
{
  time_t t1,t2;
  bat_mode = true;
  //Form1->Visible   = true; // not yet initialized!
  DecimalSeparator = '.';

  Randomize();  // <-- TODO: This is an empty function in VCL.cpp. I guess it was meant to 
  // generate a proper seed or something. Shouldn't be needed anymore, now that we have...
  init_randz();
#ifdef USE_POW_APPROX_IEEE_754
  powf_fast_init();
#endif

  ZCommandLine cmd;
  // we'll initialize VCL here
  String vers =  String("Zonation ") + zonation_VERSION_MAJOR + "." +
    zonation_VERSION_MINOR + "." zonation_VERSION_PATCH + ", build: " + __DATE__ + " " __TIME__;
  if (Parse_cmd_line(cmd)) {
    Form1->Memo1->Lines->Add(vers);
    Form1->Memo1->Lines->Add("==========================================================================");
    Form1->Memo1->Lines->Add("The Zonation software is distributed in the  hope that it will be useful, ");
    Form1->Memo1->Lines->Add("but WITHOUT ANY WARRANTY; without even the implied warranty of ");
    Form1->Memo1->Lines->Add("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.");	  
    Form1->Memo1->Lines->Add("==========================================================================");
    Form1->Memo1->Lines->Add(" For more information and conditions of use of this software, see the");
    Form1->Memo1->Lines->Add(" disclaimer in Help->About Zonation, in the graphical user interface");
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("Starting Zonation run on '" + zHostname() + "' at " + zCurrentTime());
    Form1->Memo1->Lines->Add("Locale is: " + String(setlocale(LC_ALL,NULL)));
    zTimerStart();

#ifdef COST_DIV_DOUBLE
    //Form1->Memo1->Lines->Add("-----------==================*********** WARNING: THIS IS ZIG4-DD *************=================-------------" );
#endif

#ifdef CAZ_ALTERNATING_EPS
    //Form1->Memo1->Lines->Add("-----------==================*********** WARNING: THIS IS ZIG4-AE 3/b *************=================-------------" );
#endif

    ZDATFile dat;
    if (!Read_settings_file(dat)) {
      free_data();
      String txt =  "Could not successfully read settings file " + setfn;
      MessageBox(0, txt, "Data input error", MB_OK);
      //          Delay(5);
      return false;
    }
    
    Application->ProcessMessages();
    Form1->Button8->Enabled=false;
    bat_run2(ZConf(cmd, dat));
  } else {
    //Form1->PageControl1->ActivePageIndex=1;
    //Application->Run();
    // Nothing more to say, specific error messages are printed 
    // printf(vers + ", could not process command line / project file.\n");
  }
  
#ifdef USE_POW_APPROX_IEEE_754
  powf_fast_release();
#endif
  return true;
}

/* xxxPLULA Tree conn

- Multiplication of up and down responses. Has to be accounted for in what
remains in a cell.
* local loss is fixed
* tree_nbh loss needs fixing

- Remember forced directory change above. Could this be made into option.

PLULA caz problem with ties when full distributions go at once.  Chk.

*/




