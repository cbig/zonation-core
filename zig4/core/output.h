#ifndef OUTPUT_H
#define OUTPUT_H

// Requires vmat. Call before starting opt.
void init_output_get_distr_centers_and_sums(bool do_sums=false);

float get_cost_fract(float cval);
void output_features_info_n_curves();
void output_SSI_features_info_n_curves();
void output_grid_rank(int num, int xsize, int ysize);
void output_grid_prop_rank(int num, int xsize, int ysize);
void output_grid_wrscr(int num, int xsize, int ysize);
//void output_cost_grid(int num);
void save_jpg();
//void grouped_outputs();
void output_grp_curves(Tgroups_info& groups_info);
// example o_subdir: project_distrib_smooth_layers
// reuses the same buffer for all the layers. Must be allocated/freed outside
bool generate_output_transf_layers(float* buff, const String& o_subdir, String*& trans_layers_str_append, const String must_suffix = "");
bool generate_output_transf_single_layer(size_t bf, float* buff, int xsize, int ysize,
				    const String& odirname, const String& trans_layer_str_append);

#endif /* OUTPUT_H */
