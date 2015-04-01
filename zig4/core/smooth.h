#ifndef SMOOTH_H
#define SMOOTH_H

// fedemp: commented these externs out from here, as nobody else than smooth.cpp uses them.
//
// extern float         **src_mat;

// extern float ***hmat;

// extern class GridMap EdgeFixMap;
// extern bool  use_edge_fix;
// extern float **edge_corr_mat;

// additional memory required for smoothing kernel transforms (approx)
float mem_required_smoothing(int xd, int yd);
// community similarity transform (CST)
float mem_required_cst_smoothing(int xd, int yd);
// matrix connectivity transform  (MCT)
float mem_required_mct_smoothing(int xd, int yd);

bool Do_smoothings_alloc();
bool smoothing_is_allocated();
void Do_smoothing_sp(int spnum);
bool Expand_community_similarity();
bool Do_matrix_smoothings();
void Free_smoothing_matrixes();
void Do_smoothing_interact(int resource, int consumer, float alpha, int iatype, float ia_par1);

#endif // SMOOTH_H
