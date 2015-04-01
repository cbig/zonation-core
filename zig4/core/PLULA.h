#ifndef PLULA_H
#define PLULA_H

#include "typedefs.h"

// PLULA module variables
extern const size_t MAX_NUM_PLU;
// extern const size_t MAX_SIZE_PLU;
extern size_t glob_plu_max_cnt;    // number of cells of biggest plu
extern size_t glob_plu_min_cnt;    // number of cells of smallest plu

/// planning unit index layer
extern int **PLL;
/// x indices of elements ordered by planning unit index
extern int *PLX;
/// y indices of elements ordered by planning unit index
extern int *PLY;
/// number of planning units
extern int PLcnt;

/// an array of planning units
extern struct PLU *PLvec; //this is being edited inside threads
extern int nme_cnt;
extern int PLL_xdim;
extern int PLL_ydim;
extern int *PLULA_nums_tmp; //indices in plula raster
extern int *PLULA_cnt_tmp; //number of index pixels in plula raster
extern int *PLU_uplinks;
extern int *PLU_up; //this is being edited inside threads
extern bool use_PLULA;

bool  Load_and_analyse_tree_file();
bool Load_PLULA(String fname);
void Free_PLULA_data();
float traverse_tree_down(int node, int mode, bool remove, bool nbh_loss_mode, int spnum);
float traverse_tree_up(int node, int mode, bool remove, bool nbh_loss_mode, int spnum);
bool check_PLULA();

#endif // PLULA_H
