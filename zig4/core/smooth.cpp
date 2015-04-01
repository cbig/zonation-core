#include "Unit1.h"
#include "defines.h"
#include "matrix_utils.h"
#include "GridMap.h"
#include "fftw_lnk.h"
#include "bat_run.h"
#include "LoadData.h"
#include "VCL.h"
#include <cstdio>
#include <cstring>

// fedemp: commented out, not used anywhere
//float         **src_mat;

bool smoothing_allocated = false;

float ***hmat;

class GridMap EdgeFixMap;
bool  use_edge_fix;
float **edge_corr_mat;

bool use_arb_kernels = false;
int arb_kernels_default = 1;
float arb_kernels_constant = 1.0f;
String arb_kernels_prefix = "";

const int FIRST_NUM_ARB_KERNEL = 100;
// map of arbitrary kernel functions (sequences of x-y pairs)
Arbitrary_Kernels_Map arbitrary_kernels;

float
mem_required_smoothing(int xd, int yd)
{
  return mem_required_smoothing_fftwf(xd, yd);
}

void retrieve_smoothed_distr(int sp)
{
	int    x,y;
	float  sum;
	float  maxv;

	//  Form1->Memo1->Lines->Add("xd= "+IntToStr(xd));

	sum=0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][sp]!=-1)
					sum += conn_LS[y][x];
		}
	}

	//  Form1->Memo1->Lines->Add("sum= "+FloatToStr(sum));
	maxv=0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (!vmat[y][x])
				continue;

			if (vmat[y][x][sp]!=-1)
				vmat[y][x][sp] = conn_LS[y][x]/sum;

			if (vmat[y][x][sp]>maxv)
				maxv=vmat[y][x][sp];
		}
	}

	obsmap[sp].max_val=maxv;
}

void reduce_by_smoothed_distr(int spnum, float ia_par1) // xxxiatype2
{
	int    x,y, lcnt;
	float  sum;
	float  maxv, conn_max, conn_limit;

	//  return;
	//  Form1->Memo1->Lines->Add("xd= "+IntToStr(xd));

	conn_max = 0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][spnum]!=-1)
					if (conn_LS[y][x]>conn_max)
						conn_max = conn_LS[y][x];
		}
	}

	conn_limit = ia_par1*conn_max;
	lcnt=0;

	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][spnum]!=-1)
				{
					if (conn_LS[y][x]>=conn_limit)
					{
						vmat[y][x][spnum] = 0.0f;
						lcnt++;
					}
					else
						vmat[y][x][spnum] *= (1.0f-conn_LS[y][x]/conn_limit);
				}
		}
	}

	sum=0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][spnum]!=-1)
					sum += vmat[y][x][spnum];
		}
	}

	//  Form1->Memo1->Lines->Add("sum= "+FloatToStr(sum));
	maxv=0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (!vmat[y][x])
				continue;

			if (vmat[y][x][spnum]!=-1)
				vmat[y][x][spnum] /= sum;

			if (vmat[y][x][spnum]>maxv)
				maxv=vmat[y][x][spnum];
		}
	}

	// fedemp: 20120301 - removed msg too obscure for users
	//Form1->Memo1->Lines->Add(IntToStr(lcnt) + " locations exceeded the connectivity threshold for maximum negative effect");

	obsmap[spnum].max_val=maxv;
}

void multiply_smoothed_distr_to(int spnum, float ia_par1)
{
	int    x, y, lcnt;
	float  sum;
	float  maxv, conn_max, conn_limit;

	//  return;
	//  Form1->Memo1->Lines->Add("xd= "+IntToStr(xd));

	conn_max = 0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][spnum]!=-1)
					if (conn_LS[y][x]>conn_max)
						conn_max = conn_LS[y][x];
		}
	}

	conn_limit = ia_par1*conn_max;
	lcnt=0;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][spnum]!=-1)
				{
					if (conn_LS[y][x]<=conn_limit)
					{
						if (conn_LS[y][x]>0.0f)
							vmat[y][x][spnum] *= (conn_LS[y][x]/conn_limit);
						else
							vmat[y][x][spnum]  = 0.0f;
					}
					else
						lcnt++;
				}
		}
	}

	sum=0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][spnum]!=-1)
					sum += vmat[y][x][spnum];
		}
	}

	//  Form1->Memo1->Lines->Add("sum= "+FloatToStr(sum));
	maxv=0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (!vmat[y][x])
				continue;

			if (vmat[y][x][spnum]!=-1)
				vmat[y][x][spnum] /= sum;

			if (vmat[y][x][spnum]>maxv)
				maxv=vmat[y][x][spnum];
		}
	}

	Form1->Memo1->Lines->Add(IntToStr(lcnt) + " locations exceeded the connectivity threshold for maximum resource use density.");

	obsmap[spnum].max_val=maxv;
}

void save_smoothed_distr_to_matrix(int spnum, float ia_par1,
				   bool apply_edge_corr=0,  float **ecm=0)//, float **m)
{
	int    x, y, lcnt;
	float  sum;
	float  maxv, conn_max, conn_limit;

	//  return;
	//  Form1->Memo1->Lines->Add("xd= "+IntToStr(xd));

	if (apply_edge_corr && ecm)
		for(y=0; y<yd; y++)
		{
			for(x=0; x<xd; x++)
			{
				if (vmat[y][x])
					if (vmat[y][x][spnum]!=-1)
						if ((1.0f-ecm[y][x])>0.0f)
							conn_LS[y][x]/=(1.0f-ecm[y][x]);
			}
		}



	conn_max = 0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][spnum]!=-1)
					if (conn_LS[y][x]>conn_max)
						conn_max = conn_LS[y][x];
		}
	}

	conn_limit = ia_par1*conn_max;
	lcnt=0;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][spnum]!=-1)
				{
					if (conn_LS[y][x]<=conn_limit)
					{
						if (conn_LS[y][x]>0.0f)
							hmat[y][x][spnum]= vmat[y][x][spnum]*(conn_LS[y][x]/conn_limit);
						else
							hmat[y][x][spnum] = 0.0f;
					}
					else
					{
						hmat[y][x][spnum]= vmat[y][x][spnum];
						lcnt++;
					}
				}
		}
	}

	sum=0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (vmat[y][x])
				if (vmat[y][x][spnum]>0)
					sum += hmat[y][x][spnum];
		}
	}

	//  Form1->Memo1->Lines->Add("sum= "+FloatToStr(sum));
	maxv=0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (!vmat[y][x])
				continue;

			if (vmat[y][x][spnum]>0.0f)
				hmat[y][x][spnum] /= sum;

			if (hmat[y][x][spnum]>maxv)
				maxv=hmat[y][x][spnum];
		}
	}

	Form1->Memo1->Lines->Add(IntToStr(lcnt) + " locations exceeded the connectivity threshold for maximum resource use density.");

	obsmap[spnum].max_val=maxv; // xxx is this ok? - Is maxv used before layers returned from store?
}

bool
smoothing_is_allocated()
{
  return smoothing_allocated;
}

bool
Do_smoothings_alloc()
{
  //Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("Doing connectivity smoothings: initializing");

	if(!initialize_fftwf()) {
		Form1->Memo1->Lines->Add("ERROR: Could not initialize smoothing, skipping.");
		return false;
	}

	if (!alloc_fftwf_matrixes())
	{
		Form1->Memo1->Lines->Add("ERROR: Could not allocate connectivity matrixes, skipping smoothing.");
		Form1->Memo1->Lines->Add("ERROR: This means that there is not enough free memory available in this computer.");
		Form1->Memo1->Lines->Add("ERROR: Possible reasons are: ");
		Form1->Memo1->Lines->Add("ERROR:     1) other processes are consuming memory and competing with Zonation");
		Form1->Memo1->Lines->Add("ERROR:     2) not enough memory in this computer (considering the size and number of input layers)");
		Form1->Memo1->Lines->Add("ERROR: Please check section 3.7 Data limitations & system requirements / Memory requirements");
		return false;
	}

	generate_fftwf_plans();

	smoothing_allocated = true;
	return true;
}

void Do_smoothing_sp(int spnum)
{
  float dist_alpha = spp[spnum].alpha*alpha_mult;
  // adjust distance by relative size of cells at the spp distribution center
  if (!glob_set_occur_size_weights_layer_fn.isEmpty()) {
    float initial_alpha = dist_alpha;
    float area_at_center = get_cell_area_correction_factor(distrib_centers_x[spnum], 
						     distrib_centers_y[spnum]);

    float max_area;
    if (use_cost)
      max_area = occur_size_weights_map.max_val;
    else
      max_area = costmap.max_val;

    // should never happen, but in case cell area is not specified or wrongly given as 0, 
    //   avoid alpha=0  (and distance alpha is not modified)
    if (0 >= area_at_center)
      area_at_center = max_area;

    float cell_len_ratio = sqrt(area_at_center) / sqrt(max_area);
    dist_alpha *= cell_len_ratio;
    Form1->Memo1->Lines->Add("   Note: distrib. center is x: " + FloatToStr(distrib_centers_x[spnum]) +
			     ", y: " + FloatToStr(distrib_centers_y[spnum]) +"; cell area at distrib. center: " +
			     FloatToStr(area_at_center)+" (max: "+FloatToStr(max_area)+"); alpha corrected from " + 
			     FloatToStr(initial_alpha) + " to " + FloatToStr(dist_alpha));
  }

  calc_fftwf_kernel(dist_alpha, spp[spnum].grp_arb_kernel_ds);
  copy_occ_to_fftwf_LS(spnum, 0, 1.0f);
  fftwf_connectivity();
  
  retrieve_smoothed_distr(spnum);
  Application->ProcessMessages();
}


void Do_smoothing_interact(int resource, int consumer, float alpha, int iatype, float ia_par1)
{
  String txt;

  if ((iatype<1)  || (iatype>2)) {
    Form1->Memo1->Lines->Add("PROBABLE ERROR WARNING: Only interaction types 1 and 2 are supported presently - pls check interactions file.");
    return;
  }

  calc_fftwf_kernel(alpha, spp[resource].grp_arb_kernel_ia);
  copy_occ_to_fftwf_LS(consumer-1, 0, 1.0f);
  fftwf_connectivity();

  if (iatype==1) {
    multiply_smoothed_distr_to(resource-1, ia_par1);
    txt = String("IAtype1: Layer %1 (resource, %2) transformed by connectivity to L%3 (consumer, %4), alpha=%5")
      .arg(resource).arg(spp[resource-1].fname).arg(consumer).arg(spp[consumer-1].fname)
      .arg(alpha, 0, 'g', 5);
  } else {
    reduce_by_smoothed_distr(resource-1, ia_par1);
    txt = String("IAtype2: Layer %1 (%2) transformed by negative impact of connectivity to L%3 (source, %4), alpha=%5")
      .arg(resource).arg(spp[resource-1].fname).arg(consumer).arg(spp[consumer-1].fname)
      .arg(alpha, 0, 'g', 5);
  }

  Form1->Memo1->Lines->Add(txt);
  Application->ProcessMessages();
}

// community similarity transform (CST)
float
mem_required_cst_smoothing(int xd, int yd)
{
  return (xd*yd*comm_set.comm_ct_cnt+1) * sizeof(float);
}

// matrix connectivity transform  (MCT)
float
mem_required_mct_smoothing(int xd, int yd)
{
  return (xd*yd*comm_set.ct_cnt+1) * sizeof(float);
}

bool Do_matrix_smoothings()
{
	int   loop, col, loop2, x, y;
	float prev_alpha, **m_stack[MAX_SPP_COUNT];
	char  txt[1024];

	hmat = fpmatrix(0,yd,0,xd); // compulsory
	if (!hmat)
		return false;
	for(y=0;y<=yd;y++)
		for(x=0; x<=xd; x++)
			hmat[y][x]=0;

	if (comm_set.edge_effect_fn!="")
	{
		Form1->Memo1->Lines->Add("Loading edge correction file from "+comm_set.edge_effect_fn);
		EdgeFixMap.normalize = false;     // xxxEdgeFix
		if (!EdgeFixMap.load_from_file(comm_set.edge_effect_fn, 0, 0))
		{
			use_edge_fix=false;
			Form1->Memo1->Lines->Add("WARNING: load unsuccessful - NOT using edge corrections.");
			Form1->Memo1->Lines->Add("");
		}
		else
		{
			use_edge_fix=true;
			Form1->Memo1->Lines->Add("      *** load successful");
		}
		if (use_edge_fix)
		{
			edge_corr_mat = matrix(0,yd,0,xd);
		}
	}

	for(y=0;y<yd;y++)
	{
		for(x=0;x<xd;x++)
		{
			if (vmat[y][x])
			{
				hmat[y][x]= (float *)calloc(comm_set.ct_cnt+1, sizeof(float));
				if (!hmat[y][x])
				{
					ShowMessage("out of memory allocating for matrix connectivity");
					return false; // dealloc???
				}
				for(loop=0; loop<comm_set.ct_cnt; loop++)
					hmat[y][x][loop]=0.0f;
			}
		}
	}
#if 0
	for(loop=0; loop<comm_set.ct_cnt; loop++)
	{
		m_stack[loop] = matrix(0,yd,0,xd);
		if (!m_stack[loop])
		{
			ShowMessage("Out of memory allocating matrix connectivity clalculation memory, matrix # "+IntToStr(loop));
			for(loop2=0; loop2<loop; loop2++)
				free_matrix(m_stack[loop2], 0,yd,0,xd);
			return;
		}
	}
#endif

	prev_alpha=-123456.1f;
	float ed_max=-100000;
	float ed_min=1000000;
	for(loop=0; loop<comm_set.ct_cnt; loop++)
	{
		Application->ProcessMessages();

		if (spp[loop].alpha!=prev_alpha)
		{
		  calc_fftwf_kernel(spp[loop].alpha, spp[loop].grp_arb_kernel_matrix);
		  prev_alpha=spp[loop].alpha;
		  if (use_edge_fix)
		    {
		      copy_matrix_to_fftwf_LS_edge(0, EdgeFixMap.m);
		      //              copy_matrix_to_fftwf_LS(0, EdgeFixMap.m);
		      fftwf_connectivity();
		      for(y=0; y<yd; y++)
			for(x=0; x<xd; x++)
			  {
			    edge_corr_mat[y][x] = conn_LS[y][x];
			    if (edge_corr_mat[y][x]<ed_min)
			      ed_min=edge_corr_mat[y][x];
			    if (edge_corr_mat[y][x]>ed_max)
			      ed_max=edge_corr_mat[y][x];
			  }
		      Form1->Memo1->Lines->Add("Edge correction matrix min - max values are " + FloatToStr(ed_min)+ "  " +FloatToStr(ed_max));
		    }
		}

		copy_occ_to_fftwf_LS(0, 0, sim_mat[loop][0]);
		for(col=1; col<comm_set.ct_cnt; col++)
			copy_occ_to_fftwf_LS(col, 1, sim_mat[loop][col]);
		Application->ProcessMessages();

		fftwf_connectivity();

		//      multiply_smoothed_distr_to(loop, 1.0f);
		save_smoothed_distr_to_matrix(loop, 1.0f, use_edge_fix, edge_corr_mat);//, m_stack[loop]);

		sprintf(txt, "Input feature %i transformed via matrix connectivity calculations using alpha %f.",
			loop, spp[loop].alpha);

		Form1->Memo1->Lines->Add(txt);
		Application->ProcessMessages();
	}

	for(y=0; y<yd; y++)
		for(x=0; x<xd; x++)
			if (vmat[y][x])
			{
				for(loop=0; loop<comm_set.ct_cnt; loop++)
					if (vmat[y][x][loop]!=-1)
						vmat[y][x][loop]=hmat[y][x][loop]; //m_stack[loop][y][x];
			}

	//  for(loop=0; loop<comm_set.ct_cnt; loop++)
	//    free_matrix(m_stack[loop], 0,yd,0,xd);

	if (hmat)
	{
		for(y=0; y<yd; y++)
			for(x=0; x<xd; x++)
				if (hmat[y][x])
				{
					free(hmat[y][x]);
					hmat[y][x]=0;
				}
		free_fpmatrix(hmat, 0,yd,0,xd);
		hmat=0;
	}

	if (use_edge_fix)
	{
		EdgeFixMap.free_matrix_m();
		free_matrix(edge_corr_mat, 0,yd,0,xd);
	}

	return true;
}

bool Expand_community_similarity() // 10.12.2009
{
	int   loop, col, loop2, x, y;
	char  txt[256];
	double new_psum[MAX_SPP_COUNT];

	Form1->Memo1->Lines->Add("Starting expansion of community similarity.");
	hmat = fpmatrix(0,yd,0,xd); // compulsory storage block vmat copy
	if (!hmat)
		return false;
	for(y=0;y<=yd;y++)
		for(x=0; x<=xd; x++)
			hmat[y][x]=0;

	for(y=0;y<yd;y++)
	{
		for(x=0;x<xd;x++)
		{
			if (vmat[y][x])
			{
				hmat[y][x]= (float *)calloc(comm_set.comm_ct_cnt+1, sizeof(float));
				if (!hmat[y][x])
				{
					ShowMessage("Out of memory allocating matrix for community similarity expansion");
					return false; // dealloc???
				}
				for(loop=0; loop<comm_set.comm_ct_cnt; loop++)
					hmat[y][x][loop]=0.0;
			}
		}
	}

#if 1
	for(y=0; y<yd; y++) // 19.1.2011 to fix missing data
	{
		for(x=0; x<xd; x++)
		{
			if (!vmat[y][x])
				continue;
			for(loop=0; loop<comm_set.comm_ct_cnt; loop++)
			{
				if (vmat[y][x][loop]<0.0)
					vmat[y][x][loop]=0.0;
			}
		}
	}
#endif

	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (!vmat[y][x])
				continue;
			for(loop=0; loop<comm_set.comm_ct_cnt; loop++)
			{
				if (vmat[y][x][loop]>0.0)
				{
					for(loop2=0;loop2<comm_set.comm_ct_cnt; loop2++)
						hmat[y][x][loop2] +=  comm_sim_mat[loop2][loop]*vmat[y][x][loop]*spp[loop].prob_sum; // xxx *prob_sum added 2/2010
				}
			}
		}
	}

	for(loop=0; loop<MAX_SPP_COUNT; loop++)
		new_psum[loop]=0.0;

	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (!vmat[y][x])
				continue;
			for(loop=0; loop<comm_set.comm_ct_cnt; loop++)
			{
				if (hmat[y][x][loop]>0.0)
				{
					new_psum[loop] += hmat[y][x][loop];
				}
			}
		}
	}

	for(y=0; y<yd; y++)
		for(x=0; x<xd; x++)
			if (vmat[y][x])
			{
				for(loop=0; loop<comm_set.comm_ct_cnt; loop++)
					if (vmat[y][x][loop]!=-1)
						vmat[y][x][loop]=hmat[y][x][loop]/new_psum[loop];  // psum added
			}

	Form1->Memo1->Lines->Add("Note: Community similarity expansion modifies effective original extents of features.");
	Form1->Memo1->Lines->Add("Community count = "+IntToStr(comm_set.comm_ct_cnt));
	Form1->Memo1->Lines->Add("Original  expanded");
	for(loop=0; loop<comm_set.comm_ct_cnt; loop++)
	{
		sprintf(txt, "%-5i %0.2f\t%0.2f", loop, spp[loop].prob_sum, new_psum[loop]);
		Form1->Memo1->Lines->Add(txt);
	}

	for(loop=0; loop<comm_set.comm_ct_cnt; loop++)
		spp[loop].prob_sum = new_psum[loop];

	if (hmat)
	{
		for(y=0; y<yd; y++)
			for(x=0; x<xd; x++)
				if (hmat[y][x])
				{
					free(hmat[y][x]);
					hmat[y][x]=0;
				}
		free_fpmatrix(hmat, 0,yd,0,xd);
		hmat=0;
	}

#if 1
	float tmpf=0.0f;
	for(y=0; y<yd; y++)
		for(x=0; x<xd; x++)
			if (vmat[y][x])
			{
				if (vmat[y][x][0]!=-1)
					tmpf += vmat[y][x][loop];
			}
	Form1->Memo1->Lines->Add("So distr check= "+FloatToStr(tmpf));
#endif

	Form1->Memo1->Lines->Add("Expansion of community similarity succesfully completed.");
	return true;
}

void Free_smoothing_matrixes()
{
  if (!smoothing_allocated)
    return; // ignore and don't print anything

  free_fftwf_plans();
  free_fftwf_matrixes();
  uninitialize_fftwf();
  smoothing_allocated = false;
}


