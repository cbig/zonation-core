#include "Unit1.h"
#include "bat_run.h"
#include "grid_utils.h"
#include "ccsp.h"

const bool CONN_THICKNESS = true;
extern bool save_cached_layer;
extern bool use_smart_count_labels;
extern bool layer_written;

// verbose debug file 
bool write_vdf = false;

// require 2,3,4,... consecutive in all directions from edge. -> example, thickness 2, if only 2 are found in one direction => split risk
int use_fast_mode_simple_width = 0;


void
print_mini_mat(FILE* vdf, int xmin, int xmax, int ymin, int ymax, int rx, int ry, char**& nr_mat)
{
  for (size_t j=ymin; j<=ymax; j++) {
    for (size_t i=xmin; i<=xmax; i++) {
      // avoid annoying '-1'
      if (nr_mat[j][i]<0) 
	fprintf(vdf, "-_");
      else 
	if (j==ry && i==rx)
	  fprintf(vdf, "%d*", nr_mat[j][i]);
	else
	  fprintf(vdf, "%d_", nr_mat[j][i]);
    }
    fprintf(vdf, "\n");
  }
}




#include <queue>
#include <set>

class XY_Map_Comp_f
{
public:
  
  bool operator()(XY_Coord p1, XY_Coord p2) const {
    if (p1.y < p2.y)
      return true;
    else if (p1.y == p2.y && p1.x < p2.x)
      return true;
    else
      return false;
  }
};

XY_Map_Comp_f map_comp;


bool
Grid_CCSP::find_path_local_search_bfs(/*int x, int y, */ int orig_x, int orig_y, int dest_x, int dest_y)
{
  std::queue<XY_Coord> paths;
  //std::set<XY_Coord, XY_Map_Comp_f> visited(map_comp);

  //char**& sts = status_bkp;
  //char**& sts = edge;

  //char** sts = status_bkp;
  char** sts;
  if (thickness <= 1)
    sts = edge; 
  else 
    sts = edge_bkp_eroded;
  //sts = status_bkp; 

  /*
  if (use_fast_mode_simple_width>0)
    sts = status;
  */

  memset(&(bfs_minilayer[0][0]), 0, dim_bfs_minilayer*dim_bfs_minilayer* sizeof(char));
  int x_offset = - orig_x + bfs_x_radius;
  int y_offset = - orig_y + bfs_y_radius;

  int minx = max(0, orig_x - bfs_x_radius);
  int maxx = min((int)xdim-1, orig_x + bfs_x_radius);
  int miny = max(0, orig_y - bfs_y_radius);
  int maxy = min((int)ydim-1, orig_y + bfs_y_radius);

  // local_bf_mat
  XY_Coord orig = { orig_x, orig_y };
  paths.push(orig);
  //visited.insert(orig);
  //std::cerr << "ORIG paths: " << paths.size() << ", at " << orig_x << ", " << orig_y << std::endl;

  do {
    XY_Coord p = paths.front();
    paths.pop();

    if (dest_x == p.x && dest_y == p.y) {
      //std::cerr << "Found!: " << dest_x << ", " << dest_y << std::endl;
      return true;
    }
    XY_Coord c;


    bool extra_cond = true;

    if (use_fast_mode_simple_width>0)
      extra_cond = c.x>minx && c.y>miny && c.x<maxx && c.y<maxy &&
							   (status[c.y-1][c.x-1]>0 && status[c.y][c.x-1]>0 
					       || status[c.y][c.x-1]>0 && status[c.y+1][c.x-1]>0);

    // 4-conn
    c.set(p.x-1, p.y);
    if (c.x > minx && sts[c.y][c.x] > 0 && bfs_minilayer[c.y+y_offset][c.x-1+x_offset] != 1 && extra_cond /*visited.find(c)==visited.end()*/) {
      paths.push(c);
      //visited.insert(c);
      bfs_minilayer[c.y+y_offset][c.x-1+x_offset] = 1;
      //std::cerr << "a";
    }

    if (use_fast_mode_simple_width>0)
      extra_cond = c.x>minx && c.y>miny && c.x<maxx && c.y<maxy &&
							   (status[c.y-1][c.x-1]>0 && status[c.y-1][c.x]>0 || status[c.y-1][c.x]>0 && status[c.y-1][c.x+1]>0);

    c.set(p.x, p.y-1);
    if (c.y > miny && sts[c.y][c.x] > 0 && bfs_minilayer[c.y-1+y_offset][c.x+x_offset] != 1 && extra_cond  /*visited.find(c)==visited.end()*/) {
      paths.push(c);
      //visited.insert(c);
      bfs_minilayer[c.y-1+y_offset][c.x+x_offset] = 1;
      //std::cerr << "b";
    }

    if (use_fast_mode_simple_width>0)
      extra_cond = c.x>minx && c.y>miny && c.x<maxx && c.y<maxy &&
							   (status[c.y-1][c.x+1]>0 && status[c.y][c.x+1]>0 || status[c.y][c.x+1]>0 && status[c.y+1][c.x+1]>0);

    c.set(p.x+1, p.y);
    if (c.x < maxx && sts[c.y][c.x] > 0 && bfs_minilayer[c.y+y_offset][c.x+1+x_offset] != 1 && extra_cond  /*visited.find(c)==visited.end()*/) {
      paths.push(c);
      //visited.insert(c);
      bfs_minilayer[c.y+y_offset][c.x+1+x_offset] = 1;
      //std::cerr << "c";
    }

    if (use_fast_mode_simple_width>0)
      extra_cond = c.x>minx && c.y>miny && c.x<maxx && c.y<maxy &&
							   (status[c.y+1][c.x-1]>0 && status[c.y+1][c.x]>0 || status[c.y+1][c.x]>0 && status[c.y+1][c.x+1]>0);

    c.set(p.x, p.y+1);
    if (c.y < maxy && sts[c.y][c.x] > 0 && bfs_minilayer[c.y+1+y_offset][c.x+x_offset] != 1 && extra_cond  /*visited.find(c)==visited.end()*/) {
      paths.push(c);
      //visited.insert(c);
      bfs_minilayer[c.y+1+y_offset][c.x+x_offset] = 1;
      //std::cerr << "d";
    }

    //std::cerr << "   v.size: " << visited.size() << ", p.size: " << paths.size() << "   ";
  } while (!paths.empty());
  //std::cerr << "paths END: " << paths.size() << "   " << std::endl;
  return false;
}

inline bool
Grid_CCSP::check_split_local_bfs(int rx, int ry, char**& sts)
{
  if (use_fast_mode_simple_width>0) {
    int nb4 = get_nb_count(sts, rx, ry);
    int nb8 = nb4 + get_nb_corners_count(sts, rx, ry);
    //if (nb8 < 3)  // with 2 neighbors there cannot be a connection of thickness 2 // TODO <- remove

    /* This alone (local_bfs always false? <--- NOP) generates plenty of width-1 corridors
    if (nb8 < 2)
      return false;
    */
    /* this generates exactly the same pattern!!! WTF???
    if (nb4 < 2)
      return false;
    */
    /* Same again!!!
    if (nb4 < 3)
      return false;
    */

    //    if (nb8 < 2)
    //      return false;
  }

  // as soon as a "no_alt_path" condition is found, a split is guaranteed
  bool no_alt_path = false;

  bool horiz_split = ((1 == sts[ry][rx-1] && 1 == sts[ry][rx+1])  && 
		      !(1==sts[ry-1][rx-1] && 1==sts[ry-1][rx] && 1==sts[ry-1][rx+1]) && 
		      !(1==sts[ry+1][rx-1] && 1==sts[ry+1][rx] && 1==sts[ry+1][rx+1]) );
  if (!no_alt_path && horiz_split) {
    no_alt_path |= !find_path_local_search_bfs(rx-1, ry, rx+1, ry);
    if (write_vdf && no_alt_path)
      fprintf(vdf, " ====> check_split_local_bfs(): horiz_split: %d, %d\n", rx, ry);
  }

  bool vert_split = ((1 == sts[ry-1][rx] && 1 == sts[ry+1][rx]) && 
		     !(1==sts[ry-1][rx-1] && 1==sts[ry][rx-1] && 1==sts[ry+1][rx-1]) && 
		     !(1==sts[ry-1][rx+1] && 1==sts[ry][rx+1] && 1==sts[ry+1][rx+1]) );
  if (!no_alt_path && vert_split) {
    no_alt_path |= !find_path_local_search_bfs(rx, ry-1, rx, ry+1);
    if (write_vdf && no_alt_path)
      fprintf(vdf, " ====> check_split_local_bfs(): vert_split: %d, %d\n", rx, ry);
  }
  
  bool corner1_split = ((1 == sts[ry][rx-1] && 1 == sts[ry-1][rx]) && 
			!(1==sts[ry-1][rx-1]) );
  if (!no_alt_path && corner1_split) {
    no_alt_path |= !find_path_local_search_bfs(rx-1, ry, rx, ry-1);
    if (write_vdf && no_alt_path)
      fprintf(vdf, " ====> check_split_local_bfs(): corner1_split: %d, %d\n", rx, ry);
  }
  
  bool corner2_split = ((1 == sts[ry-1][rx] && 1 == sts[ry][rx+1]) && 
			!(1==sts[ry-1][rx+1]) );
  if (!no_alt_path && corner2_split) {
    no_alt_path |= !find_path_local_search_bfs(rx, ry-1, rx+1, ry);
    if (write_vdf && no_alt_path)
      fprintf(vdf, " ====> check_split_local_bfs(): corner2_split: %d, %d\n", rx, ry);
  }
  
  bool corner3_split = ((1 == sts[ry][rx+1] && 1 == sts[ry+1][rx]) && 
			!(1==sts[ry+1][rx+1]) );
  if (!no_alt_path && corner3_split) {
    no_alt_path |= !find_path_local_search_bfs(rx+1, ry, rx, ry+1);
    if (write_vdf && no_alt_path)
      fprintf(vdf, " ====> check_split_local_bfs(): corner3_split: %d, %d\n", rx, ry);
  }
  
  bool corner4_split = ((1 == sts[ry+1][rx] && 1 == sts[ry][rx-1]) && 
			!(1==sts[ry+1][rx-1]) );
  if (!no_alt_path && corner4_split) {
    no_alt_path |= !find_path_local_search_bfs(rx, ry+1, rx-1, ry);
    if (write_vdf && no_alt_path)
      fprintf(vdf, " ====> check_split_local_bfs(): corner4_split: %d, %d\n", rx, ry);
  }

  if (write_vdf && no_alt_path) {  
    const int window = 5;
    int xmin = max(0, rx-window);
    int xmax = min((int)xdim-1, rx+window);
    int ymin = max(0, ry-window);
    int ymax = min((int)ydim-1, ry+window);    
    fprintf(vdf, "\t IN: check_split_local_bfs (%d, %d). sts (status_bkp) (%d-%d, %d-%d):\n", rx, ry, xmin, xmax, ymin, ymax);
    print_mini_mat(vdf, xmin, xmax, ymin, ymax, rx, ry, sts);
    fprintf(vdf, "\t AND edge_eroded_bkp\n", rx, ry);
    print_mini_mat(vdf, xmin, xmax, ymin, ymax, rx, ry, edge_bkp_eroded);
  }

  return no_alt_path;
}

inline bool
Grid_CCSP::check_split_local_bfs_diamond(int rx, int ry, char**& sts)
{
  bool split_detected = false;
  /*
  int inc = thickness_inc +1;
  if (rx-inc > 0)
    split_detected |= check_split_local_bfs(rx-inc, ry, sts);
  if (!split_detected && rx+inc < xdim)
    split_detected |= check_split_local_bfs(rx+inc, ry, sts);	
  if (!split_detected && ry-inc > 0)
    split_detected |= check_split_local_bfs(rx, ry-inc, sts);
  if (!split_detected && ry+inc < ydim)
    split_detected |= check_split_local_bfs(rx, ry+inc, sts);     
  */


  // THIS IS SQUARE
  int inc = thickness_inc;
  // Note the additional +/-1 (exclude border frame) cause neighbors are needed inside check_split_local_bfs
  int xmin = max(1, rx-inc);
  int xmax = min((int)xdim-2, rx+inc);
  int ymin = max(1, ry-inc);
  int ymax = min((int)ydim-2, ry+inc);

  // cache_local_bfs_diamond_split is used because when using corridor thicknes > 1, the same cells may be explored multiple times 
  // (for example from the left size and the right side of a vertical corridor.

  // Square rather than diamond!!!   
  for(int j=ymin; j<=ymax; j++) {
    for(int i=xmin; i<=xmax; i++) {
      if((j==ymin || j==ymax || i==xmin || i==xmax) && edge_bkp_eroded[j][i]>0) {  // ??? edge_bkp_eroded /OR status_bkp

	if (cache_local_bfs_diamond_split && cache_local_bfs_diamond_split[j][i]>0) {
	  split_detected = true;
	} else {
	  status_bkp[j][i] = 0;
	  edge_bkp_eroded[j][i] = 0;
	  split_detected |= check_split_local_bfs(i, j, sts);
	  status_bkp[j][i] = 1;
	  edge_bkp_eroded[j][i] = 1;

	  if (write_vdf && split_detected) {
	    fprintf(vdf, " *** diamond - split detected: %d, %d - center: %d, %d ***\n", i, j, rx, ry);
	  }
	}

	if (split_detected) {
	  if (cache_local_bfs_diamond_split)
	    cache_local_bfs_diamond_split[j][i] = 1;

	  if (Grid_CCSP::Slow == this->mode) {
	    sf_split_detected_x = i;
	    sf_split_detected_y = j;
	  }

	  return true;	    
	}
      }
    }
  }
  return split_detected;


  // THIS IS FOR the UPDATE_ERODE with simple star!!!
  inc = thickness_inc;
  bool ignore_status_on = false;

  if (rx-inc > 0) {
    if (ignore_status_on || edge_bkp_eroded[ry][rx-inc]>0) {
      status_bkp[ry][rx-inc] = 0;
      edge_bkp_eroded[ry][rx-inc] = 0;
      split_detected |= check_split_local_bfs(rx-inc, ry, sts);
      status_bkp[ry][rx-inc] = 1;
      edge_bkp_eroded[ry][rx-inc] = 1;
    }

    if (split_detected) {
      fprintf(vdf, " *** diamond (rx-inc) - split detected - center: %d, %d ***\n", rx, ry);
      return true;
    }
  }
  if (ry-inc > 0) {
    if (ignore_status_on || edge_bkp_eroded[ry-inc][rx]>0) {
      status_bkp[ry-inc][rx] = 0;
      edge_bkp_eroded[ry-inc][rx] = 0;
      split_detected |= check_split_local_bfs(rx, ry-inc, sts);
      status_bkp[ry-inc][rx] = 1;
      edge_bkp_eroded[ry-inc][rx] = 1;
    }

    if (split_detected) {
      fprintf(vdf, " *** diamond (ry-inc) - split detected - center: %d, %d ***\n", rx, ry);
      return true;
    }
  }
  if (rx+inc < xdim-1) {
    if (ignore_status_on || edge_bkp_eroded[ry][rx+inc]>0) {
      status_bkp[ry][rx+inc] = 0;
      edge_bkp_eroded[ry][rx+inc] = 0;
      split_detected |= check_split_local_bfs(rx+inc, ry, sts);
      status_bkp[ry][rx+inc] = 1;
      edge_bkp_eroded[ry][rx+inc] = 1;
    }

    if (split_detected) {
      fprintf(vdf, " *** diamond (rx+inc) - split detected - center: %d, %d ***\n", rx, ry);
      return true;
    }
  }
  if (ry+inc < ydim-1) {
    if (ignore_status_on || edge_bkp_eroded[ry+inc][rx]>0) {
      status_bkp[ry+inc][rx] = 0;
      edge_bkp_eroded[ry+inc][rx] = 0;
      split_detected |= check_split_local_bfs(rx, ry+inc, sts);
      status_bkp[ry+inc][rx] = 1;
      edge_bkp_eroded[ry+inc][rx] = 1;
    }

    if (split_detected) {
      fprintf(vdf, " *** diamond (ry+inc) - split detected - center: %d, %d ***\n", rx, ry);
      return true;
    }
  }

  fprintf(vdf, " >>> diamond, returning: %d - center: %d, %d\n", split_detected, rx, ry);
  return split_detected;



  /*

  // Traditional diamond:

  // [W->S]
  for(int j=ry; j<=ymax; j++) {
    for(int i=xmin; i<=rx; i++) {
      if (sts[j][i] >= 0 && 1 == edge_bkp_eroded[j][i])
	split_detected |= check_split_local_bfs(i, j, sts);
    }
  }
  if (split_detected)
    return true;
  // (S->E]
  for(int j=ymax-1; j>=ry; j--) {
    for(int i=rx+1; i<=xmax; i++) {
      if (sts[j][i] >= 0 && 1 == edge_bkp_eroded[j][i])
	split_detected |= check_split_local_bfs(i, j, sts);
    }
  }
  if (split_detected)
    return true;
  // (E->N]
  for(int j=ry-1; j>=ymin; j--) {
    for(int i=xmax-1; i>=rx; i--) {
      if (sts[j][i] >= 0 && 1 == edge_bkp_eroded[j][i])
	split_detected |= check_split_local_bfs(i, j, sts);
    }
  }
  if (split_detected)
    return true;
  // (N->W)
  for(int j=ymin+1; j<=ry-1; j++) {
    for(int i=rx-1; i>=xmin+1; i--) {
      if (sts[j][i] >= 0 && 1 == edge_bkp_eroded[j][i])
	split_detected |= check_split_local_bfs(i, j, sts);
    }
  }
  */

  return split_detected;
}

void
Grid_CCSP::init_next_warp_fast(float removal_threshold)
{

  if (thickness > 1 && cache_local_bfs_diamond_split) {
    for (size_t y=0; y<ydim; y++){
      for (size_t x=0; x<xdim; x++){
	cache_local_bfs_diamond_split[y][x] = 0;
      }
    }
  }

  // ***** std::cerr << "BEG init_next_warp" << std::endl;
  /*
  last_lbln = count_labels(status_bkp, lbl_buffers[0]);
  calc_boundary_lengths_n_areas_n_wrscr(status_bkp, wrscr_mat);
  */
  int ln_uneroded =0;
  if (use_smart_count_labels) {
    last_lbln = count_smart_labels_n_bl_n_areas_n_wrscr(status_bkp, lbl_buffers[current_lbl_buffer], wrscr_mat,
							boundary_lengths, areas, wrscr);

    if (thickness > 1) {
      count_smart_labels_n_bl_n_areas_n_wrscr(status, uneroded_lbl_buffer, wrscr_mat,
					      uneroded_boundary_lengths, uneroded_areas, uneroded_wrscr);
    }
    
  } else {
    // 0? why 0 and not 'lbl_buffers[current_lbl_buffer]'  <------ cause in fast mode only 0 is ever used

    if (thickness <= 1) {
      last_lbln = count_labels(status_bkp, lbl_buffers[0]);
      calc_boundary_lengths_n_areas_n_wrscr(status_bkp, wrscr_mat, last_lbln, lbl_buffers[0],
					    boundary_lengths, areas, wrscr);
    } else /*if (thickness > 1)*/ {
      // here status rather than status_bkp (eroded)
      ln_uneroded = count_labels(status, uneroded_lbl_buffer);
      calc_boundary_lengths_n_areas_n_wrscr(status, wrscr_mat, ln_uneroded, uneroded_lbl_buffer, 
					    uneroded_boundary_lengths, uneroded_areas, uneroded_wrscr);
      // also updates the cog_x, etc. vectors
    }
  }


  if (Grid_CCSP::Slow == this->mode) {
    if (thickness <= 1) { 
      sf_label_cnt_before = last_lbln;
      // TODO
    } else /*if (thickness > 1)*/ {
      
      //std::cerr << "********* POS 407 7: " << status_bkp[7][407] << std::endl;
      sf_label_cnt_before = count_labels(status_bkp, sf_lbl_buffer);
    }
    
    calc_boundary_lengths_n_areas_n_wrscr(status_bkp, wrscr_mat, sf_label_cnt_before, sf_lbl_buffer,
					  sf_before_boundary_lengths, sf_before_areas, sf_before_wrscr);
    
    // sf_label_before will be init later, in calc_ccsp
  }


  bool save_all_here = false;
  int thresh = 80200;
  if (save_all_here && thickness > 1 && !layer_written && (thresh >= (nonm1 - removed))) {
    layer_written = true;
    bool save_labels = true;
    if (save_labels) {
      debug_save_raster_int(lbl_buffers[current_lbl_buffer], "saved_layer_labels_"+IntToStr(thresh), -1);
    }
    bool save_eroded_sts = true;
    if (save_eroded_sts) {
      debug_save_raster(status_bkp, "saved_layer_status_bkp_eroded_"+IntToStr(thresh), 0);
    }
    bool save_sts = true;
    if (save_sts) {
      debug_save_raster(status, "saved_layer_status_"+IntToStr(thresh), 0);
    }
    bool save_edge = true;
    if (save_edge) {
      debug_save_raster(edge, "saved_layer_edge_"+IntToStr(thresh), 0);
    }
    bool save_eroded_edge = true;
    if (save_eroded_edge) {
      debug_save_raster(edge_bkp_eroded, "saved_layer_edge_eroded_"+IntToStr(thresh), 0);
    }
  }
}

void
Grid_CCSP::sf_get_split_sizes_n_lens_n_wrscr(int x, int y, char**& sts, 
					     std::vector<size_t>& sizes, std::vector<size_t>& lens, 
					     std::vector<float>& wrscr)
{
  int inc = 1;
    std::vector<int> labels;
  // 4 basic ones for 4-connected version
  if (sts[y][x-inc] > 0 && sf_lbl_buffer[y][x-inc] >= 0)
    labels.push_back(sf_lbl_buffer[y][x-inc]);
  if (sts[y-inc][x] > 0 && sf_lbl_buffer[y-inc][x] >= 0)
    labels.push_back(sf_lbl_buffer[y-inc][x]);
  if (sts[y][x+inc] > 0 && sf_lbl_buffer[y][x+inc] >= 0)
    labels.push_back(sf_lbl_buffer[y][x+inc]);
  if (sts[y+inc][x] > 0 && sf_lbl_buffer[y+inc][x] >= 0)
    labels.push_back(sf_lbl_buffer[y+inc][x]);

  sizes.clear();
  lens.clear();
  wrscr.clear();
  
  // Needs to distinguish connected cells that actually belong to the same blob
  // -> identify the set of CC labels
  std::set<int> cc_labels;
  for(size_t l = 0; l < labels.size(); l++) {
    cc_labels.insert(labels[l]);
  }
  std::set<int>::iterator i;
  for(i = cc_labels.begin(); i != cc_labels.end(); i++) {
    sizes.push_back(this->areas[*i]);
    lens.push_back(this->boundary_lengths[*i]);
    wrscr.push_back(this->wrscr[*i]);
  }
}


static int rule_num = 0;

static int num_total = 0;
static int num_risk = 0;
static int num_risk_and_split = 0;

float
Grid_CCSP::calc_ccsp_fast(int rx, int ry, float removal_threshold)
{
  float domain_coeff = calc_domain_coefficient(rx, ry);
  if (.0f == domain_coeff)
    return .0f;


  rule_num++;

  // NB: status (sts) here is always it's (eroded if thickness>1) backup
  char**& sts = status_bkp;

  float width_loss_factor = 1.0f;

  if (!use_8_connectivity) {

    num_total++;

    // Thickness - erosion before 
    if (CONN_THICKNESS && thickness > 1 ) {
      update_erode_foreground(rx, ry, thickness, sts, edge_bkp_eroded); //do_thinning_foreground(thickness);
      //do_dilate_background(thickness); //do_thickening_background(thickness);
    }

    /*
    // check width loss/degradation
    if (SIMPLE_THICKNESS && thickness > 1 ) {  // simple test that doesn't work-

      // get max degradation from the individual degradations in all directions 
      int xmin = max(0, rx-thickness);
      int xmax = min((int)xdim-thickness, rx+thickness);
      int ymin = max(0, ry-thickness);
      int ymax = min((int)ydim-1, ry+thickness);

      int west_len = 0, east_len = 0, south_len = 0, north_len = 0;
      for (int k=rx-1; k>=xmin; k--)
	if (sts[ry][k] >= 1)
	  west_len++;

      for (int k=rx+1; k<=xmax; k++)
	if (sts[ry][k] >= 1)
	  east_len++;

      for (int k=ry-1; k>=ymin; k--)
	if (sts[k][rx] >= 1)
	  south_len++;

      for (int k=ry+1; k<=ymax; k++)
	if (sts[k][rx] >= 1)
	  north_len++;

      int max_rem = thickness;
      if (west_len < thickness && east_len < thickness) {
	max_rem = max(west_len, east_len);
      }
      if (north_len < thickness && south_len < thickness) {
	int ns_rem = max(west_len, east_len);
	max_rem = max(max_rem, ns_rem);
      }
      if (max_rem < thickness) {
	float loss_ratio = 1.0f / (float(max_rem)+1);
	width_loss_factor = 1.0f-loss_ratio;
      }

      if (true && 0==(rule_num++%50000))
	Form1->Memo1->Lines->Add("w loss factor: "+FloatToStr(width_loss_factor));
    }
    */

    bool risk = false;
    if (!CONN_THICKNESS || thickness <= 1) {
	risk = check_basic_split_risk(rx, ry, sts);
    } else {
      // Is this =0 needed at all?
      //sts[ry][rx] = 0;
      //risk = check_basic_split_risk_diamond(rx, ry, status_bkp); // same: &sts = status_bkp

      risk = check_basic_split_risk_diamond(rx, ry, sts);

      //sts[ry][rx] = 1;
      /*if (!risk)
	risk = check_basic_split_risk(rx, ry, status);
        */

    }

    if (!risk) {
      if (CONN_THICKNESS && thickness > 1) {
	undo_update_erode_foreground(rx, ry, thickness, status_bkp, edge_bkp_eroded);  // NOT NEEDED here (if in the end update_erode does not erode (just save state)!!!
      }

      if (save_cached_layer)
	cached[ry][rx] = 0; //nonm1-removed;

      return 0;
    }

    num_risk++;

    // delete
    sts[ry][rx] = 0;
    edge[ry][rx] = 0;    // TODO!!!
    if (CONN_THICKNESS && thickness > 1) {
      edge_bkp_eroded[ry][rx] = 0;   // TODO: is this 0 enough - or what about the diamond around? -> update_erode_foreground() has been done above, but that didn't remove the thickness_inc+1 corners!!!

      /*
      // TODO ----> BUT this is what update_erode_foreground() does!!!
      // TODO: this should be: do_temporary_remove() (or so) which in thickness==1 is as simple as setting sts[ry][rx] = 0 (as 5 lines above)
      int inc = thickness_inc;
      // exclude border frame because later we need get_nb_count, etc.
      int xmin = max(1, rx-inc);
      int xmax = min((int)xdim-2, rx+inc);
      int ymin = max(1, ry-inc);
      int ymax = min((int)ydim-2, ry+inc);
      // Save state
      for (int j=ymin; j<=ymax; j++) {
	for (int i=xmin; i<=xmax; i++) {    
	  state_backup_sts_bkp[j-ymin][i-xmin] = sts[j][i];
	}
      }
      // Square rather than diamond !!!
      for(int j=ymin; j<=ymax; j++) {
	for(int i=xmin; i<=xmax; i++) {
	  sts[j][i] = 0;
	}
      }
      */
    }


    bool split_detected = false;
    if (CONN_THICKNESS && thickness > 1) {
      split_detected = check_split_local_bfs_diamond(rx, ry, sts);

      // ***
      /* this doesn't make sense with thickness > 1
      if (!split_detected) {
	status[ry][rx] = 0;
	edge[ry][rx] = 0;
	split_detected = check_split_local_bfs(rx, ry, status);
	status[ry][rx] = 1;
	edge[ry][rx] = 1;
	}*/
      // ***

    } else {
      /*
      if (use_fast_mode_simple_width>0) {
	thickness_inc = 1;
	split_detected = check_split_local_bfs_diamond(rx, ry, sts);   // contains 
      } else
      */

	split_detected = check_split_local_bfs(rx, ry, sts);   // contains 
    }


    // undelete
    status_bkp[ry][rx] = 1;
    edge[ry][rx] = 1;    // TODO!!!
    if (CONN_THICKNESS && thickness > 1) {
      // WTF NOOO!
      //      edge_bkp_eroded[ry][rx] = 1;

      /*
      // BUT: this is what undo_update_erode_foreground() does!!!
      // TODO: and this should be undo_temporary_remove()
      // Restore state
      for (int j=ymin; j<=ymax; j++) {
	for (int i=xmin; i<=xmax; i++) {    
	  sts[j][i] = state_backup_sts_bkp[j-ymin][i-xmin];
	}
      }
      */

    }

    if (CONN_THICKNESS && thickness > 1) {
      undo_update_erode_foreground(rx, ry, thickness, sts, edge_bkp_eroded);
    }

    if (!split_detected) {
      if (save_cached_layer)
	cached[ry][rx] = 0; //nonm1-removed;
      
      return 0;
    }

    if (save_cached_layer)
      cached[ry][rx] = removed; //nonm1-removed;

    float p;
    if (false  /* for now, disable cache in fast mode */ /*is_cached_penalty(rx, ry)*/) {
      p = p_cache[ry][rx];
    } else {

      if (Grid_CCSP::Slow == this->mode) {
	

	if (!CONN_THICKNESS || thickness <= 1) {
	  sf_split_detected_x = rx;
	  sf_split_detected_y = ry;
	}

	sf_label_before = sf_lbl_buffer[sf_split_detected_y][sf_split_detected_x];
	status_bkp[sf_split_detected_y][sf_split_detected_x] = 0;
	sf_label_cnt_after = count_labels(status_bkp, sf_lbl_buffer);
	calc_boundary_lengths_n_areas_n_wrscr(status_bkp, wrscr_mat, sf_label_cnt_after, sf_lbl_buffer,
					      sf_after_boundary_lengths, sf_after_areas, sf_after_wrscr);

	status_bkp[sf_split_detected_y][sf_split_detected_x] = 1;

	//std::cerr << "x: " << rx << ", y: " << ry << ", sf_split_x: " << sf_split_detected_x << ", sf_split_y: " << sf_split_detected_y << 
	//  ", lbl: " << sf_label_before << " , size: " << sf_before_areas.size() << std::endl;
	// fill in last_unconfirmed_removal (for calc_penalty_second_max_richness_over_old_richness(), etc.)
	//sf_get_split_sizes_n_lens_n_wrscr(rx, ry, status_bkp, 
	sf_get_split_sizes_n_lens_n_wrscr(sf_split_detected_x, sf_split_detected_y, status_bkp, 
					  last_unconfirmed_removal.split_sizes, last_unconfirmed_removal.split_boundary_lens, 
					  last_unconfirmed_removal.split_wrscr);

	p = calc_penalty(1 /*diff*/, sf_before_areas[sf_label_before], sf_before_wrscr[sf_label_before]);
	//std::cerr << "x: " << rx << ", y: " << ry << ", p: " << p << std::endl;


	sf_before_boundary_lengths = sf_after_boundary_lengths;
	sf_before_areas = sf_after_areas;
	sf_before_wrscr = sf_after_wrscr;
	sf_label_cnt_before = sf_label_cnt_after;


	if (true && 0==(rule_num%25000))
	  Form1->Memo1->Lines->Add("num_total: "+IntToStr(num_total)+", num_risk: "+IntToStr(num_risk)+", num_risk_and_split: "+IntToStr(num_risk_and_split)
				   + ". Was in x: "+IntToStr(rx)+", y: "+IntToStr(ry)+ ", p: "+FloatToStr(p));


      } else 
	p = calc_penalty_rule_fast_warp(rx, ry, 1.0f);
      /***********************/
      /*
      size_t area = areas[lbl];
      size_t richness = wrscr[lbl];
      float remaining = nonm1 - removed;
      //p = CCSP * float(area) / float(remaining);
      p = CCSP * float(richness) / float(remaining_wrscr);
      */
    }

    bool report_num_total_risk = false;
    num_risk_and_split++;
    if (report_num_total_risk) {
      if (true && 0==(rule_num++%50000))
	Form1->Memo1->Lines->Add("num_total: "+IntToStr(num_total)+", num_risk: "+IntToStr(num_risk)+", num_risk_and_split: "+IntToStr(num_risk_and_split)
				 + ". Was in x: "+IntToStr(rx)+", y: "+IntToStr(ry)+ ", p: "+FloatToStr(p));
    }

    bool use_verbose_debug_file = false;
    if (use_verbose_debug_file) {
      // VERBOSE DEGUB FILE
      bool was_cached = false;
      /*
	if (!is_cached_penalty(rx, ry)) {
	cached[ry][rx] = n_registered_removals;
	cached_factor[ry][rx] = 20;
	} else {
	was_cached = true;
	}
      */
      if (write_vdf && !was_cached && ry <= 250) {
	fprintf(vdf, "\n\n# Removing: %d, %d. cached: %d. thick_inc: %d. Split detected: %d Pty: %f\n", 
		rx, ry, was_cached, thickness_inc, split_detected, p);
	const int window = 5;
	int xmin = max(0, rx-window);
	int xmax = min((int)xdim-1, rx+window);
	int ymin = max(0, ry-window);
	int ymax = min((int)ydim-1, ry+window);
	fprintf(vdf, "\tstatus: (%d-%d, %d-%d)\n", xmin, xmax, ymin, ymax);
	print_mini_mat(vdf, xmin, xmax, ymin, ymax, rx, ry, status);
	
	fprintf(vdf, "\tedge: (%d-%d, %d-%d)\n", xmin, xmax, ymin, ymax);
	print_mini_mat(vdf, xmin, xmax, ymin, ymax, rx, ry, edge);
	
	fprintf(vdf, "\tstatus_bkp (eroded): (%d-%d, %d-%d)\n", xmin, xmax, ymin, ymax);
	print_mini_mat(vdf, xmin, xmax, ymin, ymax, rx, ry, status_bkp);
	
	fprintf(vdf, "\tedge_bkp (eroded): \n");
	print_mini_mat(vdf, xmin, xmax, ymin, ymax, rx, ry, edge_bkp_eroded);
      }
    }
    
    
    return domain_coeff * p;

  } else {
    // 8-conn
    return 0;
  }
}


bool
Grid_CCSP::accept_fast_removal(int rx, int ry)
{
  status_bkp[ry][rx] = 0;
  status[ry][rx] = 0;
  // edge[ry][rx] = 0;

  if (CONN_THICKNESS && thickness > 1) {
    edge_bkp_eroded[ry][rx] = 0;
    // *********************** should be STARRRRRRR/SQUAREEEEEEEEEE! 
    // update_erode_foreground(rx, ry, thickness, status_bkp, edge_bkp_eroded); //do_thinning_foreground(thickness);


    // SQUARE
    int inc = thickness_inc;
    int xmin = max(1, rx-inc);
    int xmax = min((int)xdim-2, rx+inc);
    int ymin = max(1, ry-inc);
    int ymax = min((int)ydim-2, ry+inc);
    // TODO: This should be do/confirm_erosion_update() or similar
    for(int j=ymin; j<=ymax; j++) {
      for(int i=xmin; i<=xmax; i++) {
	if ((i==xmin || i==xmax || j==ymin || j==ymax) && edge_bkp_eroded[j][i]>0) {
	  status_bkp[j][i] = 0;
	  edge_bkp_eroded[j][i] = 0;
	  edgify_neighbors(i, j, status_bkp, edge_bkp_eroded);
	}
      }
    }	

  }
  
  if (Grid_CCSP::Slow == this->mode) {
    
    // sf_label_before = sf_lbl_buffer[sf_split_detected_y][sf_split_detected_x];
    int& lbl = sf_label_before;
    sf_before_areas[lbl]--;
    if (sf_before_areas[lbl] <= 0) {
      sf_before_boundary_lengths[lbl] = 0;
      sf_before_areas[lbl] = 0;
      sf_before_wrscr[lbl] = 0;
      sf_label_cnt_after--;
    } else {
      sf_before_wrscr[lbl] -= wrscr_mat[sf_split_detected_y][sf_split_detected_x];
    }
    remaining_wrscr -= wrscr_mat[sf_split_detected_y][sf_split_detected_x];
  }
  
    // OLD version, simple star ---------------- TODO: This should be do/confirm_erosion_update() or similar
    /*
      char**& sts = status_bkp;
    if (rx>0) {
      // (rx-inc)
      sts[ry][rx-thickness_inc] = 0;
      edge_bkp_eroded[ry][rx-thickness_inc] = 0;
      edgify_neighbors(rx-thickness_inc, ry, sts, edge_bkp_eroded);
    }
    if (ry>0){
      // (ry-inc)
      sts[ry-thickness_inc][rx] = 0;
      edge_bkp_eroded[ry-thickness_inc][rx] = 0;
      edgify_neighbors(rx, ry-thickness_inc, sts, edge_bkp_eroded);      
    }
    if (rx<xdim-1) {
      // (rx+inc)
      sts[ry][rx+thickness_inc] = 0;
      edge_bkp_eroded[ry][rx+thickness_inc] = 0;
      edgify_neighbors(rx+thickness_inc, ry, sts, edge_bkp_eroded);
    }
    if (ry<ydim-1){
      // (ry+inc)
      sts[ry+thickness_inc][rx] = 0;
      edge_bkp_eroded[ry+thickness_inc][rx] = 0;
      edgify_neighbors(rx, ry+thickness_inc, sts, edge_bkp_eroded);      
    }
    */

    //  }

  do_save_corridor_boundaries_layer();

  return true;
}

// (rx,ry): coordinates of cell to remove
float
Grid_CCSP::calc_penalty_rule_fast_warp(size_t rx, size_t ry, float width_loss_factor)
{
  float p;

  // This only works for thickness==1, otherwise it'll be -1 all the time!!!
  // int lbl = get_cc_label(rx, ry);
  int lbl;
  if (thickness <= 1)
    lbl = get_cc_label(rx, ry);
  else
    lbl = uneroded_lbl_buffer[ry][rx];

  if (lbl<0)
    	Form1->Memo1->Lines->Add("ERRRRRR===========================, calc_penalty_rule_fast_warp(), x: "+IntToStr(rx)+
				 ", y: "+IntToStr(ry)+", lbln: "+IntToStr(lbl)+
				 ", areas.size(): "+IntToStr(areas.size())+
			       ", edge: "+IntToStr(edge[ry][rx])+
			       ", edge_bkp_eroded: "+IntToStr(edge_bkp_eroded[ry][rx])+
			       ", status: "+IntToStr(status[ry][rx])+
			       ", status_bkp: "+IntToStr(status_bkp[ry][rx])
				 );


  //float radius = sqrt(cog_std_x[lbl]*cog_std_x[lbl] + cog_std_y[lbl]*cog_std_y[lbl]);
  // note cog_std_x/y are actually var!
  // MOVE THIS TO THE END OF calc_boundary_lengths_n_areas_n_wrscr() - needed only once per lbl
  float radius;
  float r = 1.0/cog_sum_weights[lbl] * (float(cog_n[lbl])/float(cog_n[lbl]-1));
  if (11 <= rule_formula && rule_formula <= /*14*/ 20) {   // area
    radius = sqrt(cog_std_x[lbl]/(cog_n[lbl]-1) + cog_std_y[lbl]/(cog_n[lbl]-1));
  } else { // (1 <= rule_formula && rule_formula <= /*4*/ 9)  // richness
    radius = sqrt(cog_std_x[lbl]*r + cog_std_y[lbl]*r);
  }
  float diff_x = (rx-cog_x[lbl]);
  float diff_y = (ry-cog_y[lbl]);
  float dist = sqrt(diff_x*diff_x + diff_y*diff_y);
  float dist_factor = radius/dist; /// TEST works!!!!

  //dist_factor = dist_factor/dist;
  dist_factor *= dist_factor;

  /*
  if (width_loss_factor < 1.0f)
    width_loss_factor /= 5000.0f;
  */

  if (false && 0==(rule_num++%50000)) {
    //Form1->Memo1->Lines->Add("D: "+FloatToStr(dist)+", R: "+FloatToStr(radius)+", cog_std_x: "+FloatToStr(cog_std_x[lbl]/(cog_n[lbl]-1))+", cog_std_y: "+FloatToStr(cog_std_y[lbl]/(cog_n[lbl]-1)));
    Form1->Memo1->Lines->Add("D: "+FloatToStr(dist)+", R: "+FloatToStr(radius)+"r: "+FloatToStr(r)+", cog_std_x: "+FloatToStr(cog_std_x[lbl]*r)+", cog_std_y: "+FloatToStr(cog_std_y[lbl]*r));
    //rule_num--;
  }

  if (thickness <= 1) {
    
    // richness variants
    if (1 <= rule_formula && rule_formula <= 4) {
      p = CCSP * (width_loss_factor) * dist_factor * float(wrscr[lbl]) / float(remaining_wrscr);
    } else { // if (11 <= rule_formula && rule_formula <= 14)
      // area variants
      float remaining = nonm1 - removed;
      p = CCSP * (width_loss_factor) * dist_factor * float(areas[lbl]) / float(remaining);
    }
  } else {

    /* ***************************************************************** uneroded agg ******************************* */
    if (1 <= rule_formula && rule_formula <= 4) {
      p = CCSP * (width_loss_factor) * dist_factor * float(uneroded_wrscr[lbl]) / float(remaining_wrscr);
    } else { // if (11 <= rule_formula && rule_formula <= 14)
      // area variants
      float remaining = nonm1 - removed;
      p = CCSP * (width_loss_factor) * dist_factor * float(uneroded_areas[lbl]) / float(remaining);
    }
  }

  return p;
}
