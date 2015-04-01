#ifndef LOADDATA_H
#define LOADDATA_H

#include "typedefs.h"
#include <map>

// forward declarations
struct ZDATFile;
struct ZCommandLine;
struct Tgroups_info;
class ZConf;

// rr4
extern int  new_spp_file_format;

// condition xxx2.1
extern class  GridMap condition; // also used for retention v3

// retention xxx3.0
// class  GridMap retention;

// post processing
extern String PPA_fname;

// community settings
extern struct  community_settings comm_set;
extern float   **sim_mat, **comm_sim_mat;
extern float   ret_layers_rel_weight;

// area mask
extern class GridMap area_mask;

// Administrative units
extern struct ADM_settings  ADM_set;
extern float **ADM_weights, **ADM_combined_weights, **ADMUxSP_repr, **ADMUxSP_repr_orig, **ADMUxspp_loss;
extern bool  *loss_in_ADU;
extern size_t* ADMU_n_non_missing;
extern double* ADMU_cost, * ADMU_cost_in_area, * ADMU_cost_used, * ADMU_area_in;
extern class GridMap ADMU_map;
extern int **ADMUs;
// maps: ADMU ids => sequence numbers (in the sequence of ADMUs given in the ADMU description file)
extern std::map<int, int> ADMUs_id_to_seq;

void Fix_fname(String &fname);
void Fix_output_file_names(ZConf const& conf);
void free_vmat();  // frees only vmat, included in free_data, but should be called asap to free the big thing
void free_wrscr_Rsum(); // frees matrices used for calculating wrscr (Rsum) - use only after status/edge have been init
void free_wrscr_Rmax(); // frees matrices used for calculating wrscr (Rsum) - use only after status/edge have been init
void free_opt_structures();  // matrices that are not needed anymore after priorit. opt. finishes
void free_data();
void Fix_nbm();
void Fix_BQP_nbms();
void get_cost_info();
void get_per_ADMU_area_and_cost_info_and_non_missing();
void get_per_ADMU_area_and_non_missing();
bool Parse_cmd_line(ZCommandLine& command);
bool Read_settings_file(ZDATFile& dat);
bool data_input_v2();
bool Read_IGw_file(int num, const String& fname, float **wtmpm);
float BQP_interpolate(int remains,  int max_nbc, int prof_num);

#endif /* LOADDATA_H */

