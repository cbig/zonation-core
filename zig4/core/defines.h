#ifndef DEFINES_H
#define DEFINES_H

#include "typedefs.h"
#include <vector>
class GridMap;

#define MAX_SPP_COUNT 65500
//#define _CRT_SECURE_NO_WARNINGS

// Used for curves[][] and SSI_curves[][]
const int CURVES_LEN = 10000;

/// A struct representing a species distribution
struct sp
{
  bool drop; // whether to drop/ignore it because it has no
             // occurrencies within the analysis area

	String fname;
	float  weight;
	float  alpha;
	float  prob_sum;
	float  IG_fract;
	float  max_local_val;
	float  sr_par;
	float  T_violation_fract;
	int    BQP_radius;
	int    BQP_prof_num;
	int    BQP_prof_down;
	int    nbmat_num;
	bool   BQP_used;
	bool   TREE_used;
	float  *BQP_nbv;
	float  *BQP_dnbv;
	float  *BQP_link;
	float  *BQP_link2;

	float  w2;  // rr4 pars
	float  Tj;
	float  exp1;
	float  exp2;

	int    grp_op1;
	int    grp_op2;
	int    grp_cond;
	int    grp_ret;
	int    grp_ret_mode;
        //int    grp_LEC;  // deprecated
        int    grp_arb_kernel_ds;
        int    grp_arb_kernel_matrix;
        int    grp_arb_kernel_ia;

	bool   DSd;
	bool   IAd;
	bool   CONDd;
	bool   IGd;
};

struct SSI_sp
{
	String fname;
	float  weight;
	float  prob_sum;
	int    N;
	int    *x;
	int    *y;

	float  sr_par;
	int    BQP_radius;
	int    BQP_prof_num;
	int    nbmat_num;
	bool   BQP_used;
	float  *BQP_nbv;
	float  *BQP_dnbv;
	float  *BQP_link;
	float  lost_at_fraction;
	float  T_violation_fract;

	float  w2;  // rr4 pars
	float  Tj;
	float  exp1;
	float  exp2;
};

struct SSI_xydata
{
	int   spnum;
	int   x;
	int   y;
	float val;
};

// Info about groups (max number and array of values). 
// Calculated in loaddata.cpp, and then used in bat_run2 to generate the per 
// ADMU .grp_curves.txt and the global .grp_curves.txt files.
struct Tgroups_info
{
  int max_grp;
  int* grpv;
  std::vector< std::vector<int> > spp_indices; // spp included in every group
};

struct sort_xy
{
	int x;
	int y;
	float val;
};

struct spot
{
	int num;
	int area;
	int nwn;
	float *bdv;
	float rank;
	int linked_to;
	bool in_min_percent;

	float Ex;
	float Ey;
	float one_x;
	float one_y;

	int   min_gx;
	int   min_gy;
	int   max_gx;
	int   max_gy;
	float mean_gx;
	float mean_gy;
};

struct IG_settings
{
	float  IGa;
	int    use_IG;
	int    use_IGw;
	int    normalize_IGw;
	int    IG_proportional;
	float  IG_spw[MAX_SPP_COUNT];
	String fnames[MAX_SPP_COUNT];
	String IGwfn;
};

struct community_settings
{
	int    load_sim_matrix;
	String sim_matrix_name;
	String comm_sim_matrix_name;
	int    ct_cnt;
	int    comm_ct_cnt;
	int    apply_to_conn;
	int    apply_to_repr;
	String edge_effect_fn;
};

struct BQP_profile
{
	int   num;
	int   pnt_cnt;
	float prop[21];
	float rem[21];
	float loss10000[10001];
};

/// A struct representing a planning unit
struct PLU // size 12 ints
{
	int   num;
	int   basin;

	/// index of the first element in planning unit in PLX and PLY arrays
	int   start;

	/// number of elements in the planning unit
	int   el_cnt;

	int   minx;
	int   maxx;
	int   miny;
	int   maxy;
	bool  at_edge;
	bool  removed;
	bool  checked;

	/// delta value
	float dv;
	int   PLULA_BQPmode;
	bool  allocated;
	float *datavec;
	float cost;
	// tree conn data
	int   down_link;
	int   up_link_cnt;
	int   up_link_start;
	float orig_tree_conn_up;
	float tree_conn_up;
	float orig_tree_conn_down;
	float tree_conn_down;
};

struct ADM_desc
{
  int number; // in [0, ADM_set.count-1]
  int id_number; // ID used in ADMU raster map
  float global_weight;
  float local_weight;
  String name;
};

struct ADM_settings
{
  int use_ADMUs;
  int ADMU_mode;
  int calc_from_condition;
  int count;
  float mode2_global_weight;
  int row_count_for_per_admu_curves;
  String ADM_layer_file;
  String ADM_weights_file;
  String ADM_weight_matrix_file;
  struct ADM_desc *descriptions;  // from ADMU descriptions input file
  // per-ADMU 'curves' and 'group curves' files are written on the fly
  // these vectors contain the file pointers, opened with the output_...init()
  // functions and closed in the last write iteration
  std::vector<FILE*> curves_files;
  std::vector<FILE*> grp_curves_files;
};

struct Corr_settings
{
  int    use_corr;
  int    use_domain_layers;
  String layers_fn;
  int    IG_proportional;
  std::vector<float> weights;
  std::vector<String> file_names;
  std::vector<GridMap> layers;
  // Percentages (rank) at which the 'corridor boundaries layers' 
  // should be generated
  std::vector<float> out_boundaries_pcs;
  float start_pc;
};

// These structures seem to not include file names (historic reasons?...)
// transformed layers: [Transformed layers] in settings file
struct transf_layers_settings
{
  bool some_output;
  bool output_final_layers;
  bool output_ig_layers;
  bool output_ds_layers;
  bool output_cst_layers;
  bool output_ct_layers;
  bool output_rt_layers;
  bool output_mct_layers;
  bool output_ia_layers;  
};

#endif /* DEFINES_H */
