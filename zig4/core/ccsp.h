#ifndef CCSP_H
#define CCSP_H

#include "VCL.h"
#include "conn_components.h"

typedef struct  {
  int x;
  int y;

  inline void set(int new_x, int new_y)
  { x = new_x; y = new_y; }
} XY_Coord;

struct Grid_CC_Removal_State
{
  int x, y;
  int lbl;
  int num_labels;
  int diff;
  std::vector<size_t> sizes;
  std::vector<size_t> boundary_lens;
  std::vector<float> wrscr;
  std::vector<size_t> split_sizes;
  std::vector<size_t> split_boundary_lens;
  std::vector<float> split_wrscr;
};

class Grid_CCSP: public Grid_Connected_Components
{
 public:

  typedef enum { Fast_Warp = 0, Slow = 1 } Calc_Mode;

  /*
   * @param ofn: output file name. If empty string, no output is generated
   * @param variant: "fast"/warp mode, or "slow" mode
   */
  Grid_CCSP(float CCSP, int variant, int formula,
	    int thickness, bool use_8_conn, String ofn, size_t xdim, size_t ydim, char**& status);
  
  ~Grid_CCSP();

  Calc_Mode calc_mode()
  { return mode; }

  bool is_cached_penalty(int x, int y);
  float calc_cached_penalty(int x, int y);
    
  /*
   * Writes output info for the last round of removal iterations
   */
  bool output_info_line(float prop_lost);

  // delta: marginal loss without morpho penalty
  float calc_ccsp_as_if_removed(size_t x, size_t y, float delta, float removal_threshold);


  void discount_from_areas_wrscr_etc(int x,int y);
  void slow_init_removal(size_t x, size_t y);

  bool slow_accept_removal(int x, int y);
  bool slow_cancel_removal(size_t x, size_t y);
  bool recache_corridor(size_t x, size_t y);
  bool exterminate_corridor_cache(size_t x, size_t y);


  // Penalty rule formulas. Fast/approximate calculation
  float calc_penalty_rule_fast_warp(size_t rx, size_t ry, float width_loss_factor);

  // Penalty rule formulat. Exact version
  float calc_penalty(int diff, size_t area, float richness);
  // Alternatives:
  float calc_penalty_old_area_prop(int diff, size_t area, float richness);
  float calc_penalty_old_richness_prop(int diff, size_t area, float richeness);
  float calc_penalty_min_area_over_old_area_prop(int diff, size_t area, float richness);
  float calc_penalty_min_area_over_blob_area(int diff, size_t area, float richness);
  float calc_penalty_min_richness_over_old_richness_prop(int diff, size_t area, float richness);
  float calc_penalty_min_richness_over_blob_richness(int diff, size_t area, float richness);
  float calc_penalty_second_max_area_over_old_area(int diff, float area, float richness);
  float calc_penalty_second_max_richness_over_old_richness(int diff, float area, float richness);


  Grid_CC_Removal_State last_confirmed_removal;
  Grid_CC_Removal_State last_unconfirmed_removal;
  bool confirmed_removal_state;

  size_t boundary_len;
  /* std::vector<size_t>* */ 
  void calc_boundary_lengths_n_areas_n_wrscr(char**& status, float**& wrscr_mat, int last_lbln,
					     int**& lbl_buffer,
					     std::vector<size_t>& boundary_lengths,
					     std::vector<size_t>& areas,
					     std::vector<float>& wrscr);


  // quick/smart labelling of connected components and aggregation of boundary/area/wrscr
  int count_smart_labels_n_bl_n_areas_n_wrscr(char**& status, int**& lbl_layer, float**& wrscr_mat,
					      std::vector<size_t>& boundary_lengths,
					      std::vector<size_t>& areas,
					      std::vector<float>& wrscr);

  int label_cc_boundary(int orig_x, int orig_y, int lbln, char**& edge, int**& lbl_layer);


  //std::vector<size_t>* calc_areas();

  // 0: empty map
  // -1: error
  int count_labels_as_if_removed(int x, int y, char**& sts, float delta, float removal_threshold);

  static size_t n_penalty_calc;
  
  int last_x, last_y /*, last_size*/;
  bool current_removal_required_ccl;
  int safe_xy_lbl;
  int safe_last_lbln;
  int safe_last_last_lbln;
  int safe_size;

  std::vector<size_t> safe_areas;  
  std::vector<float> safe_wrscr;
  std::vector<size_t> safe_boundary_lengths;

  // For Fast_Warp mode:
  std::vector<float> cog_x, cog_x_prev;
  std::vector<float> cog_y, cog_y_prev;
  std::vector<size_t> cog_n;
  std::vector<float> cog_sum_weights;
  std::vector<float> cog_std_x; //, cog_std_x_prev;
  std::vector<float> cog_std_y; //, cog_std_y_prev;

  int last_lbln;
  int last_last_lbln;

  // need 2 to keep state, used to be int** lbl_layer;
  std::vector<int**> lbl_buffers;
  // switch between lbl buffers
  size_t current_lbl_buffer;

  // current CC label of cell x,y
  inline int get_cc_label(int x, int y) {
    //return lbl_layer[y][x];
    return (lbl_buffers[current_lbl_buffer])[y][x];
  }

  void update_bl(int x, int y, int inc) {
    //boundary_lengths[lbl_layer[y][x]] += inc;
    boundary_lengths[(lbl_buffers[current_lbl_buffer])[y][x]] += inc;
  }

  std::vector<size_t> boundary_lengths;


  void undiscount_from_areas_wrscr_etc(int x, int y, int lbl);
  size_t cancel_removal_count;

 private:
  float calc_domain_coefficient(int rx, int ry);

  bool do_save_corridor_boundaries_layer();

  void debug_save_raster(char**& m, String oname, float nodatavalue=-1.0f);
  void debug_save_raster_int(int**& m, String oname, int nodatavalue=-1);

  // Updates status and edge (sts and edge_bkp_eroded)
  void update_erode_foreground(int x, int y, size_t thickness, char**& sts_eroded, char**& edge_eroded);
  void edgify_neighbors(int x, int y, char**& sts, char**& edge_bkp_eroded);
  void undo_update_erode_foreground(int x, int y, size_t thickness, char**& sts, char**& edge_eroded);

  // sts_orig: normal/original status matrix (in)
  // sts: status matrix to erode (in/out)
  // 
  void do_erode_foreground(char**& sts_orig, char**& sts, size_t thickness, char**& edge_bkp);

  float CCSP;
  Calc_Mode mode;
  int rule_formula, old_variant, thickness, thickness_inc;
  FILE* of;

  // TODO: only for debug
  FILE *vdf;

  size_t xdim, ydim;

 public:
  // size (in # of cell) of connected components (cc)
  std::vector<size_t> areas;
  std::vector<float> wrscr;

  // Used only when thickness > 1
  std::vector<size_t> uneroded_boundary_lengths;
  std::vector<size_t> uneroded_areas;
  std::vector<float> uneroded_wrscr;
  int** uneroded_lbl_buffer;
  int uneroded_last_lbln;


  char** erosion_mask;
  char** state_backup_sts_bkp;
  char** erosion_update_state_sts;
  char** erosion_update_state_edge;
  // TODO: not needed anymore
  //  std::multimap<int, int> equivalences;


  // ref to the Z global 'char** status'
  char**& status;

  // this one is eroded as needed (if thickness > 1)
  char** status_bkp;
  char** edge_bkp_eroded;

  char** cc_status;


 public:
  char** cache_local_bfs_diamond_split;


  bool last_cached;
  // difference in number of labels
  int** diff_cache;
  float** p_cache;
  int** cached;
  float** cached_factor;
  float global_cache_factor;
  size_t num_cached, num_ever_cached;

  void increase_cache_factor();
  void decrease_cache_factor();


  bool estimate_split_risk(char**& sts, int rx, int ry, bool try_local_search=false);
  size_t dim_bfs_minilayer;
  char** bfs_minilayer;

  char** layer_smart_lbl_visited;

  // Warp/fast mode:
 public:
  void init_next_warp_fast(float removal_threshold);
  bool accept_fast_removal(int rx, int ry);
  float calc_ccsp_fast(int rx, int ry, float removal_threshold);
  bool check_basic_split_risk(int rx, int ry, char**& sts);

  bool check_basic_split_risk_thickness_bigger1(int rx, int ry, char**& sts);
  bool check_basic_split_risk_diamond(int rx, int ry, char**& sts);
  bool check_split_local_bfs(int rx, int ry, char**& sts);
  bool check_split_local_bfs_diamond(int rx, int ry, char**& sts);
  bool find_path_local_search_bfs(int orig_x, int orig_y, int dest_x, int dest_y);

 private:
  //int** local_bfs_mat; <--- duplicated! bfs_minilayer
  // size_t max_bfs_dim;


  // exp.
 private:
  size_t n_registered_removals;


  int get_get_get_cc_label(int x, int y);

  void get_split_sizes_n_lens_n_wrscr(int x, int y, std::vector<size_t>& sizes, std::vector<size_t>& lens, 
				      std::vector<float>& wrscr);

  float** wrscr_mat;
  float remaining_wrscr;
  bool use_8_conn;
  int max_dist_around(size_t x, size_t y);

//const int bfs_x_radius = 50;
//const int bfs_y_radius = 50;
//static const int bfs_x_radius = 200;
//static const int bfs_y_radius = 200;
static const int bfs_x_radius = 250;
static const int bfs_y_radius = 250;
//  static const int bfs_x_radius = 120;
//  static const int bfs_y_radius = 120;


  // For the new Slow/Fast mode  (uses lbl_buffer[current_lbl_buffer==0))
 private:

  void
    sf_get_split_sizes_n_lens_n_wrscr(int x, int y, char**& sts, 
				      std::vector<size_t>& sizes, std::vector<size_t>& lens, std::vector<float>& wrscr);

  int sf_label_before;
  int sf_label_cnt_before;  

  // coordinates where there's risk of split
  int sf_split_detected_x;
  int sf_split_detected_y;

  int** sf_lbl_buffer;

  int sf_label_after;
  int sf_label_cnt_after;
  char** sf_status_minilayer;
  std::vector<size_t> sf_before_boundary_lengths, sf_after_boundary_lengths;
  std::vector<size_t> sf_before_areas, sf_after_areas;
  std::vector<float> sf_before_wrscr, sf_after_wrscr;
};

#endif //  CCSP
