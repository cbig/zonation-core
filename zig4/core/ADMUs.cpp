#include "Unit1.h"
#include "defines.h"
#include "matrix_utils.h"
#include "ADMUs.h"
#include "LoadData.h"
#include "GridMap.h"
#include "bat_run.h"
#include <cstdio>

// Administrative units
struct ADM_settings  ADM_set;
float **ADM_weights, **ADM_combined_weights, **ADMUxSP_repr, **ADMUxSP_repr_orig, **ADMUxspp_loss;
bool  *loss_in_ADU;
// \sum ADMU_n_non_missing should be nonm1!
size_t* ADMU_n_non_missing;
// array of per-admu costs, updated every output interval:
// ADMU_cost = ADMU_cost_in_area - ADMU_cost_used
double* ADMU_cost, * ADMU_cost_in_area, * ADMU_cost_used, * ADMU_area_in;
class GridMap ADMU_map;
int **ADMUs;
std::map<int, int> ADMUs_id_to_seq;  // from id in descriptions (and ADMU raster) to sequence 0...n

const String ADMUWeights_append = ".ADMU_weights.txt";

static std::vector<float> ADMU_redistributed_sol;

int alloc_ADMU_data()
{
	ADM_weights = matrix(0, map_cnt, 0, ADM_set.count);
	if (!ADM_weights)
	{
		Form1->Memo1->Lines->Add("Out of memory allocating ADMU weights matrix");
		return 0;
	}

	ADM_combined_weights = matrix(0, map_cnt, 0, ADM_set.count);
	if (!ADM_combined_weights)
	{
		Form1->Memo1->Lines->Add("Out of memory allocating 2nd ADMU weights matrix");
		return 0;
	}

	ADMUxSP_repr = matrix(0, ADM_set.count, 0, map_cnt);
	if (!ADMUxSP_repr)
	{
		Form1->Memo1->Lines->Add("Out of memory allocating ADMU repr matrix");
		return 0;
	}

	ADMUxSP_repr_orig = matrix(0, ADM_set.count, 0, map_cnt);
	if (!ADMUxSP_repr_orig)
	{
		Form1->Memo1->Lines->Add("Out of memory allocating ADMU repr orig matrix");
		return 0;
	}

	ADMUxspp_loss = matrix(0, ADM_set.count, 0, map_cnt);
	if (!ADMUxspp_loss)
	{
		Form1->Memo1->Lines->Add("Out of memory allocating ADMU repr orig matrix");
		return 0;
	}

	loss_in_ADU = new bool[ADM_set.count]();
	for(int a=0; a<ADM_set.count; a++)
	  loss_in_ADU[a]=false;


	if(ADM_set.use_ADMUs && 
	   (use_groups || ADM_set.row_count_for_per_admu_curves > 0) ) {
	  // Note they all are init. to 0 - some loops depend on that
	  ADMU_cost = new double[ADM_set.count]();
	  ADMU_cost_in_area = new double[ADM_set.count]();
	  ADMU_cost_used = new double[ADM_set.count]();
	} else {
	  ADMU_cost = ADMU_cost_in_area = ADMU_cost_used = NULL;
	}
	// required for the redistributed rank
	ADMU_n_non_missing = new size_t[ADM_set.count]();
	ADMU_area_in = new double[ADM_set.count]();

	return 1;
} 

int free_ADMU_data()
{
	if (ADM_set.descriptions)
	{
		delete[] ADM_set.descriptions;
	}

	if (ADMUs)
	{
		free_imatrix(ADMUs, 0, yd, 0, xd);
		ADMUs=0;
	}

	if (ADM_weights)
	{
		free_matrix(ADM_weights, 0, map_cnt, 0, ADM_set.count);
		ADM_weights=0;
	}

	if (ADM_combined_weights)
	{
		free_matrix(ADM_combined_weights, 0, map_cnt, 0, ADM_set.count);
		ADM_combined_weights=0;
	}

	if (ADMUxSP_repr)
	{
		free_matrix(ADMUxSP_repr, 0, ADM_set.count, 0, map_cnt);
		ADMUxSP_repr=0;
	}

	if (ADMUxSP_repr_orig)
	{
		free_matrix(ADMUxSP_repr_orig, 0, ADM_set.count, 0, map_cnt);
		ADMUxSP_repr_orig=0;
	}

	if (ADMUxspp_loss)
	{
		free_matrix(ADMUxspp_loss, 0, ADM_set.count, 0, map_cnt);
		ADMUxspp_loss=0;
	}

	delete[] loss_in_ADU;
	
	delete[] ADMU_n_non_missing, ADMU_cost, ADMU_cost_in_area, ADMU_cost_used, ADMU_area_in;

	return 1;
}

int get_ADMU_data()
{
  int x, y, loop, adu;
  for(adu=0; adu<=ADM_set.count; adu++) {
    for(loop=0; loop<map_cnt; loop++)	{
      ADMUxSP_repr[adu][loop] = 0.0;
      ADMUxSP_repr_orig[adu][loop]= 0.0;
    }
  }

  for(y=0;y<yd;y++) {
    for(x=0;x<xd;x++) {
      if (vmat[y][x]) {
	adu = ADMUs[y][x];
	if (adu<0 || adu >= ADM_set.count)
	  //adu=0;
	  continue;
	
	//for(loop=0; loop<map_cnt; loop++)
	Biodiv_Features_Occur_Container rowp;
	rowp= vmat[y][x];
	for(loop = rowp.first(); loop != rowp.overflow(); loop = rowp.next(loop)) {
	  if (vmat[y][x][loop]>0.0f) {
	    ADMUxSP_repr[adu][loop] += vmat[y][x][loop];
	    ADMUxSP_repr_orig[adu][loop] += vmat[y][x][loop];
	  }
	}
      }
    }
  }
  return 1;
}


int ADMUs_output_1()
{
	int    adu, loop;
	String aw_fname;
	FILE   *f;

	aw_fname  = outgridfn;
	aw_fname  = ChangeFileExt(aw_fname, ".");
	Fix_fname(aw_fname);
	aw_fname  = ChangeFileExt(aw_fname, ADMUWeights_append);
	f = fopen(aw_fname.toUtf8().constData(), "a+t");

	if (f)
	{
		fprintf(f,"\n");
		if (!use_condition)
			fprintf(f,"feature x ADMUs representation levels (of total - sums to 1.0) (one row per feature, one column per ADMU)\n");
		else
			fprintf(f,"feature x ADMUs representation levels (of total after condition transformation - can sum to <=1.0) (one row per feature, one column per ADMU)\n");
	}

	if (f)
		for(loop=0; loop<map_cnt; loop++)
		{
			for(adu=0; adu<ADM_set.count; adu++)
			{
				fprintf(f,"%-7.5f ", ADMUxSP_repr_orig[adu][loop]);
			}
			fprintf(f,"\n");
		}


	if (f)
		fclose(f);

	return 1;
}

static String per_admu_out_dir_suffix = ".per_ADMU_outputs";
static String per_admu_curves_suffix = ".curves.txt";

//ADMU.XX.curves.txt files
int ADMUs_output_per_admu_curves_init()
{
  String projname = ChangeFileExt(outgridfn, "");
  String subdirname = createSubDirIfNeeded(projname, per_admu_out_dir_suffix);

  for (int adu=0; adu < ADM_set.count; adu++) {
    
    String fname  = subdirname + "/ADMU." + IntToStr(ADM_set.descriptions[adu].id_number) + per_admu_curves_suffix;
    Fix_fname(fname);
    FILE* f = fopen(fname.toUtf8().constData(), "w");
    
    if(!f) {
      Form1->Memo1->Lines->Add("* Error opening per administrative unit curves file: "+fname);
      return 0;
    }

    ADM_set.curves_files.push_back(f);

    fprintf(f,"ADMU id: %d: Prop_landscape_lost, cost_needed_for_top_fraction, minimum_prop_remaining, average_prop_remaining, Weighted_prop_remaining,  per_feature_remaining_prop_at_level_of_removal...\n", ADM_set.descriptions[adu].id_number);

    //fclose(f);

  }

  return 1;
}

int ADMUs_output_per_admu_curves_iter(float prop_lost)
{
  String projname = ChangeFileExt(outgridfn, "");
  String subdirname = createSubDirIfNeeded(projname, per_admu_out_dir_suffix);

  // 1 file per ADMU
  for (int adu=0; adu < ADM_set.count; adu++) {
    FILE* f = ADM_set.curves_files[adu];
    if (!f)
      return 0;
    
    // avoid '-' sign -- would take some bytes!
    if (1.0 == prop_lost) {
      for (size_t spp=0; spp < map_cnt; spp++) {
	if (ADMUxSP_repr[adu][spp] <= 0 ) {
	  ADMUxSP_repr[adu][spp] = 0.0f;
	}
      }
    }

    // float loss = 0;
    float avg_prop = 0.0f;
    float wavg_prop = 0.0f;
    float min_prop = std::numeric_limits<float>::infinity();
    float tmp_weight_agg = 0.0f;
    for (int s=0; s < map_cnt; s++) {
      // loss += ADMUxspp_loss[adu][s];
      float val = ADMUxSP_repr[adu][s];
      if (val < 0.0f) // beware numerical issues
	val = 0.0f;

      if ( val < min_prop)
	min_prop = val;
      avg_prop += val;
      tmp_weight_agg += spp[s].weight;
      wavg_prop += spp[s].weight * val;
    }
    if (map_cnt >0)
      avg_prop /= map_cnt;
    if (tmp_weight_agg >0)
      wavg_prop /= tmp_weight_agg;

    // 1 line, fraction_remaining, 4 more things and then all the species
    if (ADMU_cost[adu] < 0.0f)
      ADMU_cost[adu] = 0.0d;
    fprintf(f, "%-9.6f%-10.5g%-8.5f%-8.5f%-8.5f", prop_lost, ADMU_cost[adu], min_prop, 
	    avg_prop, wavg_prop);
    for (int spp=0; spp < map_cnt; spp++) {
      fprintf(f, "%-7.4f", ADMUxSP_repr[adu][spp]);
    }
    fprintf(f,"\n");
  }

  // very last iteration. Done with these files
  if (1.0 == prop_lost) {
    for(size_t i=0; i < ADM_set.curves_files.size(); i++) {
      fclose(ADM_set.curves_files[i]);
      ADM_set.curves_files[i] = NULL;
    }
  }

  return true;
}

static String per_admu_grp_curves_suffix = ".grp_curves.txt";

int ADMUs_output_per_admu_grp_curves_init(Tgroups_info& groups_info)
{
  String projname = ChangeFileExt(outgridfn, "");
  String subdirname = createSubDirIfNeeded(projname, per_admu_out_dir_suffix);

  for (int adu=0; adu < ADM_set.count; adu++) {
    String fname  = subdirname + "/ADMU." + IntToStr(ADM_set.descriptions[adu].id_number) + 
      per_admu_grp_curves_suffix;
    Fix_fname(fname);
    FILE* f = fopen(fname.toUtf8().constData(), "w");
    
    if(!f)
      return 0;

    ADM_set.grp_curves_files.push_back(f);

    // fedemp: for uniformity, I use this "F-lost, TF-cost ... " format (as in
    // output_grp_curves, rather than the general .curves.txt format
    fprintf(f,"ADMU id: %d: F-lost, TF-cost\t", ADM_set.descriptions[adu].id_number);
    for(size_t pos=0; pos < groups_info.max_grp; pos++) {
      int g = groups_info.grpv[pos];
      fprintf(f, "min-%i mean-%i wmean-%i max-%i\t", g, g, g, g);
    }
    fprintf(f, "\n");  
  }

  groups_info.spp_indices.resize(groups_info.max_grp);
  for(size_t pos=0; pos < groups_info.max_grp; pos++) {
    // init tables/lists of spp included in this group - reused many times in ADMUS_output_..._grp_curves_iter()
    int grp_num = groups_info.grpv[pos];
    for (size_t s=0; s < map_cnt; s++) {
      if (spp[s].grp_op1 == grp_num) {
	groups_info.spp_indices[pos].push_back(s);
      }
    }
  }

  return ADM_set.count;
}

int ADMUs_output_per_admu_grp_curves_iter(float prop_lost, Tgroups_info& groups_info)
{
  String projname = ChangeFileExt(outgridfn, "");
  String subdirname = createSubDirIfNeeded(projname, per_admu_out_dir_suffix);

  for (int adu=0; adu < ADM_set.count; adu++) {
    FILE* f = ADM_set.grp_curves_files[adu];
    if(!f)
      return 0;

    if (ADMU_cost[adu] < 0.0f)
      ADMU_cost[adu] = 0;
    fprintf(f, "%0.5f\t%-10.5g\t", prop_lost, ADMU_cost[adu]);
    
    for (size_t grp_i=0; grp_i<groups_info.max_grp; grp_i++) {
      size_t cnt = 0;
      float mean = 0.0f;
      float wmean = 0.0f;
      float wmean_accum = 0.0f;
      float minv = std::numeric_limits<float>::infinity();
      float maxv = -std::numeric_limits<float>::infinity();

      // this is map_cnt/number_of_groups iterations
      for (size_t s_idx=0; s_idx < groups_info.spp_indices[grp_i].size(); s_idx++) {
	int s = groups_info.spp_indices[grp_i][s_idx];
	float val = ADMUxSP_repr[adu][s];
	if (val < 0.0f) // beware numerical issues
	  val = 0.0f;
	mean += val;
	wmean += spp[s].weight * val;
	if (val < minv)
	  minv = val;
	if (val > maxv)
	  maxv = val;
	wmean_accum += spp[s].weight;
	cnt++;
      }

      if (cnt > 0) {
	mean /= cnt;
	if (wmean > .0f) {
	  wmean /= wmean_accum;
	}
      } else {
	minv = mean = wmean = maxv = 0.0f;
      }
      fprintf(f, "%-8.5f\t%-8.5f\t%-8.5f\t%-8.5f\t", minv, mean, wmean, maxv);

    }
    // all groups done
    fprintf(f, "\n");
  }

  // very last iteration. Done with these files
  if (1.0 == prop_lost) {
    for(size_t i=0; i < ADM_set.grp_curves_files.size(); i++) {
      fclose(ADM_set.grp_curves_files[i]);
      ADM_set.grp_curves_files[i] = NULL;
    }
  }

  return true;
}


// this function will crash (qsort) if vector size is 0 - check before!
int
sort_intra_admu_removals(std::vector<struct sort_xy>& cells)
{
  qsort((void *)&cells[0], cells.size(), sizeof(struct sort_xy), sortfunc);
  
  return cells.size();
}

// redistributed vals in [0,1], uniform spacing
int
redistribute_ranks(std::vector<struct sort_xy>& rank)
{
  int result = 1;
  // just to avoid div-by-0
  if (1 == rank.size()) {
    rank[0].val = 1.0f;
    return result;
  }

  if (!use_occur_size_weights_correct_ranking) {
    for (size_t i=1; i<=rank.size(); i++) {
      rank[i-1].val = float(i)/rank.size();
    }
  } else {

    // init (reverse) accumulators for the adu-specific rankings
    double* ADMU_area_remaining = new double[ADM_set.count]();
    for (size_t a=0; a<ADM_set.count; a++) {
      ADMU_area_remaining[a] = ADMU_area_in[a];
    }

    for (size_t i=0; i<rank.size(); i++) {
      int adu = ADMUs[rank[i].y][rank[i].x];
      float cell_area = get_cell_area_correction_factor(rank[i].x, rank[i].y);

      ADMU_area_remaining[adu] -= cell_area;
      rank[i].val = 1.0d - ADMU_area_remaining[adu]/ADMU_area_in[adu];
    }

    delete [] ADMU_area_remaining;
  }

  return result;
}

int
calculate_ADMU_redistributed_rank(int num)
{
  int result = 1;

  float nodataval = -1;
  std::vector< std::vector<struct sort_xy> > ADMU_ranks;
  // A) redistribute
  ADMU_ranks.resize(ADM_set.count);

  for (size_t admu=0; admu<ADM_set.count; admu++) {  
    // allocate admu-vector, sort it by global val, and redistribute vals
    ADMU_ranks[admu].resize(ADMU_n_non_missing[admu]);
  }
  // fill in admu vectors - split map into pieces
  std::vector<size_t> ADMU_ranks_count;
  ADMU_ranks_count.resize(ADM_set.count, 0);
  for(int y=0; y<yd; y++) {
    for(int x=0; x<xd; x++) {
      //if (-1 != Rmax[y][x]) {
      if (-1 != status[y][x]) {
	int admu_seq = ADMUs[y][x];
	if (admu_seq < 0 || admu_seq >= ADM_set.count)
	  continue;

	size_t i = ADMU_ranks_count[admu_seq];
	ADMU_ranks[admu_seq][i].x = x;
	ADMU_ranks[admu_seq][i].y = y;
	ADMU_ranks[admu_seq][i].val = sol[y][x];
	ADMU_ranks_count[admu_seq]++;
      }
    }
  }
  // sort and do redist
  for (size_t admu=0; admu<ADM_set.count; admu++) {
    // beware of empty ADMUs!
    if (ADMU_ranks[admu].size() > 0) {
      sort_intra_admu_removals(ADMU_ranks[admu]);
      redistribute_ranks(ADMU_ranks[admu]);
    }
  }

  // B) put into output matrix (buffer), ready for saving
  int xsize = xd, ysize = yd;
  ADMU_redistributed_sol.resize(ysize * xsize, nodataval); // note nodataval background  
  for(size_t admu=0; admu < ADM_set.count; admu++) {
    for (size_t i=0; i<ADMU_ranks[admu].size(); i++) {
      int x = ADMU_ranks[admu][i].x;
      int y = ADMU_ranks[admu][i].y;
      int cell_pos =  x + y * xsize;
      if(sol[y][x] < 0) {
	ADMU_redistributed_sol[cell_pos] = nodataval;
      } else {
	ADMU_redistributed_sol[cell_pos] = ADMU_ranks[admu][i].val;
	// This would be the exact same solution (could be useful for testing)
	//ADMU_redistributed_sol[cell_pos] = sol[y][x].val;
      }
    }
  }
  
  return result;
}

int
output_ADMU_redistributed_rank(int num)
{
  int result = 1;
  int xsize = xd, ysize = yd;

  float nodataval = -1.0f;  
  SaveToRaster<float>(ChangeFileExt(admu_redist_outgridfn), &ADMU_redistributed_sol[0], nodataval, xsize, ysize);
  ADMU_redistributed_sol.resize(0);
  return result;
}

int calculate_and_output_ADMU_feature_weights()
{
	int    adu, s, loop;
	float  wL, wG, glob_pref, local_pref;
	String aw_fname;
	FILE   *f;

	aw_fname  = outgridfn;
	aw_fname  = ChangeFileExt(aw_fname, ".");
	Fix_fname(aw_fname);
	aw_fname  = ChangeFileExt(aw_fname, ADMUWeights_append);
	f = fopen(aw_fname.toUtf8().constData(), "w+t");

	fprintf(f,"WGlobal    W(first ADMU)    W(second ADMU)    .... \n");
	if (ADM_set.ADMU_mode==1)
		fprintf(f, "Note: Mode 1 weights integrate both local and global components.\n");
	else
		fprintf(f, "Note: Mode 2 weights are the effective feature-specific local component.\n");

	for (s=0; s<map_cnt; s++)
	{
		fprintf(f,"%-7.4f", spp[s].weight);
		for(adu=0; adu<ADM_set.count; adu++)
		{
			wG = spp[s].weight;
			wL = ADM_weights[s][adu];
			glob_pref  = ADM_set.descriptions[adu].global_weight;
			local_pref = ADM_set.descriptions[adu].local_weight;
			if (ADM_set.ADMU_mode==1)
				ADM_combined_weights[s][adu] = glob_pref*(local_pref*wL+(1.0-local_pref)*wG);
			else
				ADM_combined_weights[s][adu] = glob_pref*wL; // ADMxxx mode 2 weight calculation Q
			if (f)
				fprintf(f,"%0.5f ", ADM_combined_weights[s][adu]);
		}
		if (f)
		{
			fprintf(f,"\n");
		}
	}


	if (f)
		fclose(f);

	return 1;
}
