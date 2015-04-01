#ifndef BAT_RUN_H
#define BAT_RUN_H

#include "typedefs.h"

extern int    xd, yd, curve_store_interval, ADM_store_interval;
extern double *repr_lvls, *repr_after_rem;
extern int   rem_by_rarity, rem_step_by_rarity, c_pos, use_edge_rem, use_smoothing;
extern double cost_used, cost_in_area, cost_out, cost_rat, cost_in_init_rem;
extern float alpha_mult;
extern float max_repr;
extern String name_addition;
extern bool  bat_mode;
extern int   annotate_name, all_zero_as_miss;

extern int    **nbm;

extern int   warp_factor, warp_ll, min_mask_level;
extern int   resamp_vec[];

extern std::vector<int> wrx, wry, wnbc;
extern std::vector<float> wval;
extern const size_t MAX_SIZE_PLU;

extern int   BL;

// variables for solution loading
extern int    run_mode;
extern class  GridMap loaded1;
extern String load1fn;
extern struct sort_xy *load_order;
extern int    load_max, loaded_cnt, logit;
extern struct IG_settings IG_set;

// change of data structure
// species data array (matrix of arrays)
#ifdef COMPACT_VMAT
 #include "occur_container.h"
 extern Biodiv_Features_Occur_Container **vmat;
#else
 extern float ***vmat;
#endif

// BQP variables
extern struct BQP_profile BQPs[];
extern int    radii_cnt, radii[], prof_cnt, radii_sp[];
extern int    **nbms[], nbm_num[];
extern float  **nbms_orig[];
extern float  smult[], BQP_spp_losses[], sp_derivative[];
extern float  BQP_spp_loss_tmp[]; //this is being edited inside threads!
extern int    BQP_mode;

// SSI variables
extern struct SSI_sp  SSI[];
extern struct SSI_xydata   *SSI_list;
extern int    SSI_spp_cnt, SSI_loc_cnt, SSI_row_cnt, SSI_lpos, SSI_err_cnt;
extern int    **SSIxyPos, **SSIxyCnt;
//extern int    *SSI_sp_num;   // fedemp: old SSI stuffcommented out, apparently not used
//extern float  *SSI_occ_size;
extern double *SSI_repr_lvls;
extern float  SSI_curves[][3];

void graceful_exit();

bool bat_run2();
float get_sp_loss_at_node_downriver(int spnum, int pln, float conn_loss);
float get_sp_loss_at_node_upriver(int spnum, int pln, float conn_loss);

///	Calculates marginal loss value for a single element.
///	\param x x-coordinate of the element
///	\param y y-coordinate of the element
///	\param single_cell_mode Evalute only single element, regardless of PLULA
///	\return calculated marginal loss
float get_delta_value(int x, int y, bool single_cell_mode);
float PLULA_tree_loss_only_local(int pln, int spnum);
float PLULA_tree_loss(int pln, int spnum);
float get_PLULA_info(int pln, float *PLdv, int *PLSSI, float *PLSSIval, int &SSIcnt);
float rr4(float fract, float w1, float Tj, float w2, float exp1, float exp2);
float PLULA_BQP_loss(int pln, int spnum);
float BQP_dv_buf(int x, int y, int s);

#endif // BAT_RUN_H
