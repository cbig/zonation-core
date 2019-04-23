#include "ccsp.h"

#include "Unit1.h"
#include "bat_run.h"
#include "randz.h"
#include "grid_utils.h"

#include <iostream>
//bool use_smart_count_labels = true;
bool use_smart_count_labels = false;

// TODO: add the final option for this: ? "save corridor map"?
bool save_cached_layer = true;

bool layer_written = false;

// Works for 4- and 8-connectivity

/*
template<class T> void
get_max_et_al(std::vector<T>& v, T* max, T* min, T* second_max)
//void get_max_et_al(std::vector<size_t>& v, size_t* max, size_t* min, size_t* second_max)
{
  // Assumes v.size() >= 2
  if (v.size() < 2)
    return;
  //std::cerr << "i1" << "max: " << *max << ", min: " << *min << "second_max: " << *second_max << ", v.size(): " << v.size() << std::endl;
  if (CCSP_verbose > 2) {
    std::cerr << "> get_max_...(): v[] = ";
    for (size_t i=0; i < v.size(); i++) {
      std::cerr << v[i] << " ";
    }
    std::cerr << std::endl;
  }
  if (v[0] > v[1]) {
    //std::cerr << ">" << "v[0]: " << v[0] << std::endl;
    *max = v[0];
    *min = *second_max = v[1];
  } else {
    //std::cerr << "no >" << "v[1]: " << v[1] << std::endl;
    *max = v[1];
    *min = *second_max = v[0];
  }
  for (size_t i=2; i < v.size(); i++) {
    if (v[i] < *min)
      *min = v[i];

    if (v[i] >= *max) {
      // second_max may be == max, if the max value is repeated
      *second_max = *max;
      *max = v[i];
    } else if (v[i] > *second_max) {
      *second_max = v[i];
    }
  }
}
*/

template<class T> void
get_max_et_al(std::vector<T>& v, T& maxi, T& mini, T& second_maxi)
// void get_max_et_al(std::vector<size_t>& v, size_t* maxi, size_t* mini, size_t* second_maxi)
{
  // Assumes v.size() >= 2
  if (v.size() < 2)
    return;
  if (CCSP_verbose > 2) {
    //std::cerr << "i1" << "maxi: " << *maxi << ", mini: " << *mini << "second_maxi: " << *second_maxi << ", v.size(): " << v.size() << std::endl;
    std::cerr << "> get_max_...(): v[] = ";
    for (size_t i=0; i < v.size(); i++) {
      std::cerr << v[i] << " ";
    }
    std::cerr << std::endl;
  }
  if (v[0] > v[1]) {
    //std::cerr << ">" << "v[0]: " << v[0] << std::endl;
    maxi = v[0];
    mini = second_maxi = v[1];
  } else {
    //std::cerr << "no >" << "v[1]: " << v[1] << std::endl;
    maxi = v[1];
    mini = second_maxi = v[0];
  }
  for (size_t i=2; i < v.size(); i++) {
    if (v[i] < mini)
      mini = v[i];

    if (v[i] >= maxi) {
      // second_max may be == max, if the max value is repeated
      second_maxi = maxi;
      maxi = v[i];
    } else if (v[i] > second_maxi) {
      second_maxi = v[i];
    }
  }

  //  std::cerr << "******** i1, " << "maxi: " << maxi << ", mini: " << mini << "second_maxi: " << second_maxi << ", v.size(): " << v.size() << std::endl;
}

// used from get_delta_value_plus_penalties!
float
Grid_CCSP::calc_cached_penalty(int x, int y)
{
  int lbln;
  if (thickness <= 1)
    lbln = get_cc_label(x, y);
  else
    lbln = uneroded_lbl_buffer[y][x];

  if (lbln<0)
	Form1->Memo1->Lines->Add("ERRRR_____________________________________, calc_cached_penalty(), x: "+IntToStr(x)+
				 ", y: "+IntToStr(y)+", lbln: "+IntToStr(lbln)+
				 ", areas.size(): "+IntToStr(areas.size())+
			       ", edge: "+IntToStr(edge[y][x])+
			       ", edge_bkp_eroded: "+IntToStr(edge_bkp_eroded[y][x])+
			       ", status: "+IntToStr(status[y][x])+
			       ", status_bkp: "+IntToStr(status_bkp[y][x])
				 );


  size_t area = areas[lbln];
  float richness = wrscr[lbln];
  int diff = diff_cache[y][x];

  // p_cache!
  float p = p_cache[y][x];
  return p;

  return cc_accounting->calc_penalty(diff, area, richness);
}

template <typename T> int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

float
Grid_CCSP::calc_penalty_second_max_richness_over_old_richness(int diff, float area, float richness)
{
  if (diff > 0) {
    int remaining = nonm1 - removed;
    float max_rich, min_rich, second_max_rich;
    get_max_et_al(last_unconfirmed_removal.split_wrscr, max_rich, min_rich, second_max_rich);

    /*
    Form1->Memo1->Lines->Add("calc_penalty_2nd_max_richness, diff: "+IntToStr(diff)
			     +", elements: "+IntToStr(last_unconfirmed_removal.split_wrscr.size())
			     +", safe_xy_lbl: "+IntToStr(safe_xy_lbl)
			     +", area: "+IntToStr(area)
			     +", rich: "+FloatToStr(richness)
			     +", 2nd_max: "+FloatToStr(second_max_rich)
			     +", max_rich:"+FloatToStr(max_rich)
			     +", min_rich:"+FloatToStr(min_rich)
			     +", remaining_wrscr: "+FloatToStr(remaining_wrscr)
			     );
    */
    return CCSP * float(second_max_rich) /* / float(richness) * float(richness) */ / remaining_wrscr; // // *sqrt(area) *  remaining;
  } else {
    return 0;
  }
}

float
Grid_CCSP::calc_penalty_second_max_area_over_old_area(int diff, float area, float richness)
{
  if (diff > 0) {
    int remaining = nonm1 - removed;
    size_t max_size, min_size, second_max_size;
    get_max_et_al(last_unconfirmed_removal.split_sizes, max_size, min_size, second_max_size);
    // LAST CHANGE / get_max_et_al(last_unconfirmed_removal.split_sizes, &max_size, &min_size, &second_max_size);
    return CCSP * float(second_max_size) /* / float(area) * float(area) */ / remaining;   // // *sqrt(area) *  remaining;
  } else {
    return 0;
  }
}

float
Grid_CCSP::calc_penalty_min_richness_over_old_richness_prop(int diff, size_t area, float richness)
{
  if (diff > 0) {
    int remaining = nonm1 - removed;
    float max_wrscr, min_wrscr, second_max_wrscr;
    get_max_et_al(last_unconfirmed_removal.split_wrscr, max_wrscr, min_wrscr, second_max_wrscr);
    //std::cerr << "*** diff: " << diff << ", richness: " << richness << ", min_wrscr: " << min_wrscr << ", max_wrscr: " << max_wrscr << ", 2nd_max_wrscr: " << second_max_wrscr << std::endl;
    return CCSP * float(min_wrscr) /* / float(richness) * float(richness) */ / remaining_wrscr; // // *sqrt(area) *  remaining;
  } else {
    return 0;
  }
}

float
Grid_CCSP::calc_penalty_min_richness_over_blob_richness(int diff, size_t area, float richness)
{
  if (diff > 0) {
    int remaining = nonm1 - removed;
    float max_wrscr, min_wrscr, second_max_wrscr;
    get_max_et_al(last_unconfirmed_removal.split_wrscr, max_wrscr, min_wrscr, second_max_wrscr);
    return CCSP * float(min_wrscr) / float(richness);
  } else {
    return 0;
  }
}

float
Grid_CCSP::calc_penalty_min_area_over_old_area_prop(int diff, size_t area, float richness)
{
  if (diff > 0) {
    int remaining = nonm1 - removed;
    size_t max_size = 0, min_size=0, second_max_size=0;
    // LAST CHANGE: get_max_et_al(last_unconfirmed_removal.split_sizes, &max_size, &min_size, &second_max_size);
    get_max_et_al(last_unconfirmed_removal.split_sizes, max_size, min_size, second_max_size);

    //std::cerr << "*** diff: " << diff << ", area: " << area << ", min_size: " << min_size << ", max_size: " << max_size << ", 2nd_max_size: " << second_max_size << std::endl;
    return CCSP * float(min_size) /* / float(area) * float(area) */ / remaining; // // *sqrt(area) *  remaining;
    // WHY this one doesn't behave well??????????????????
    //return CCSP * float(min_size) * float(min_size) / float(area) / remaining; // // *sqrt(area) *  remaining;
  } else
    return 0;
}

float
Grid_CCSP::calc_penalty_min_area_over_blob_area(int diff, size_t area, float richness)
{
  if (diff > 0) {
    int remaining = nonm1 - removed;
    size_t max_size = 0, min_size=0, second_max_size=0;
    // LAST CHANGE: get_max_et_al(last_unconfirmed_removal.split_sizes, &max_size, &min_size, &second_max_size);
    get_max_et_al(last_unconfirmed_removal.split_sizes, max_size, min_size, second_max_size);
    return CCSP * float(min_size) / float(area);
  } else
    return 0;
}

float
Grid_CCSP::calc_penalty_old_richness_prop(int diff, size_t area, float richness)
{
  int remaining = nonm1 - removed;
  if (diff > 0)
    return CCSP * float(diff) * richness / remaining_wrscr; // // *sqrt(area) *  remaining;
  else
    return 0;
}

float
Grid_CCSP::calc_penalty_old_area_prop(int diff, size_t area, float richness)
{
  int remaining = nonm1 - removed;

  //return CCSP * float(diff) * float(area) / remaining; // // *sqrt(area) *  remaining;
  if (diff > 0)
    return CCSP * float(diff) * float(area) / remaining; // // *sqrt(area) *  remaining;
  else
    return 0;
}

float
Grid_CCSP::calc_penalty(int diff, size_t area, float richness)
{
  // richness variants
  if (1 == rule_formula)
    return calc_penalty_second_max_richness_over_old_richness(diff, area, richness);
  else if (2 == rule_formula)
    return calc_penalty_old_richness_prop(diff, area, richness);
  else if (3 == rule_formula)
    return calc_penalty_min_richness_over_blob_richness(diff, area, richness);
  else if (4 == rule_formula)
    return calc_penalty_min_richness_over_old_richness_prop(diff, area, richness);
  // area variants
  else if (11 == rule_formula) 
    return calc_penalty_second_max_area_over_old_area(diff, area, richness);
  else if (12 == rule_formula) 
    return calc_penalty_old_area_prop(diff, area, richness);
  else if (13 == rule_formula) 
    return calc_penalty_min_area_over_blob_area(diff, area, richness);
  else // if (14 == rule_formula)
    return calc_penalty_min_area_over_old_area_prop(diff, area, richness);
}


////////////////////////////////////////////////////////////////////////////////////////////
// (x,y): starting point on edge
// lbl_layer has been prev. init to negative (-1). Valid label numbers are >=0
#include <queue>
int
Grid_CCSP::label_cc_boundary(int orig_x, int orig_y, int lbln, char**& edge, int**& lbl_layer)
{
  //memset(&(layer_smart_lbl_visited), 0, xdim*ydim* sizeof(int));

  //printf("AAAAA\n");

  std::queue<XY_Coord> edge_paths;
  XY_Coord orig = { orig_x, orig_y };  
  edge_paths.push(orig);

  do {
    XY_Coord p = edge_paths.front();
    edge_paths.pop();

    if (lbl_layer[p.y][p.x]>=0)
      continue;

    lbl_layer[p.y][p.x] = lbln;

    XY_Coord next;
    // 4-neighbors
    // west
    next.set(p.x-1, p.y);
    if (next.x>0 && edge[next.y][next.x]>0 && lbl_layer[next.y][next.x]<0) {
      edge_paths.push(next);
    }
    // east
    next.set(p.x+1, p.y);
    if (next.x<xdim && edge[next.y][next.x]>0 && lbl_layer[next.y][next.x]<0) {
      edge_paths.push(next);
    }
    // south
    next.set(p.x, p.y+1);
    if (next.y<ydim && edge[next.y][next.x]>0 && lbl_layer[next.y][next.x]<0) {
      edge_paths.push(next);
    }
    // north
    next.set(p.x, p.y-1);
    if (next.y>0 && edge[next.y][next.x]>0 && lbl_layer[next.y][next.x]<0) {
      edge_paths.push(next);
    }
    
  } while (!edge_paths.empty());

  return 0;
}

bool save_smart_layers = false;
int smart_it = 0;
int
Grid_CCSP::count_smart_labels_n_bl_n_areas_n_wrscr(char**& sts, int**& lbl_layer, float**& wrscr_mat,
						 std::vector<size_t>& boundary_lengths,
						 std::vector<size_t>& areas,
						 std::vector<float>& wrscr)
{
  int thresh = 80000;
  if (!layer_written && (thresh >= (nonm1 - removed))) {
    layer_written = true;
    bool save_labels = true;
    if (save_labels) {
      debug_save_raster_int(lbl_buffers[current_lbl_buffer], "saved_layer_labels_"+IntToStr(thresh));
    }
  }
  if (thickness > 1 && !layer_written && (thresh >= (nonm1 - removed))) {
    layer_written = true;
    bool save_eroded_sts = true;
    if (save_eroded_sts) {
      debug_save_raster(status_bkp, "saved_layer_status_bkp_eroded_"+IntToStr(thresh));
    }
    bool save_sts = true;
    if (save_sts) {
      debug_save_raster(status, "saved_layer_status_"+IntToStr(thresh));
    }
    bool save_edge = true;
    if (save_edge) {
      debug_save_raster(edge, "saved_layer_edge_"+IntToStr(thresh));
    }
    bool save_eroded_edge = true;
    if (save_eroded_edge) {
      debug_save_raster(edge_bkp_eroded, "saved_layer_edge_eroded_"+IntToStr(thresh));
    }
  }



  bool use_test_impl = true;
  if (use_test_impl) {
    // test impl. (as in non smart / normal labeling mode)
    //clock_t t_begin = clock();

    // Init here to get rid of old labels outside of current status
  for(size_t i = 0; i < ydim; i++) {	    
    for(size_t j = 0; j < xdim; j++) {
      lbl_layer[i][j] = -1;
    }
  }

    int lbl_num = count_labels(sts, lbl_buffers[current_lbl_buffer]);
    calc_boundary_lengths_n_areas_n_wrscr(sts, wrscr_mat, lbl_num, lbl_buffers[current_lbl_buffer], 
					  boundary_lengths, areas, wrscr);
    //float time = float(clock()-t_begin)/CLOCKS_PER_SEC;
    //Form1->Memo1->Lines->Add("Time slow: "+FloatToStr(time));

  if (0== ((smart_it++)%100))
    Form1->Memo1->Lines->Add("lbln: "+IntToStr(lbl_num) +", boundary len: "+IntToStr(boundary_len)+
			     +", area[0]: "+IntToStr(areas[0])+", wrscr[0]: "+FloatToStr(wrscr[0])
			     +", area[end]: "+IntToStr(areas[areas.size()-1])+", wrscr[end]: "+FloatToStr(wrscr[wrscr.size()-1])
			     );  

    return lbl_num;
  }

  //t_begin = clock();

  // count len, area, wrsnr, etc. and fill in inner cells of the cc
  char**& lbl_edge = edge;

  //printf("A\n");
  // A) Label edge cells
  //  memset(&(lbl_layer[0][0]), -1, xdim*ydim* sizeof(int));
  // This whole-map init should be replaced by init of the edge list coordinates (only relevant ones for the next loop)
  for(size_t i = 0; i < ydim; i++) {	    
    for(size_t j = 0; j < xdim; j++) {
      lbl_layer[i][j] = -1;
    }
  }

  int lbln = 0;
  for(size_t y=1; y<ydim-1; y++) {
    for(size_t x=1; x<xdim-1; x++) {
      if (1==lbl_edge[y][x] /* implied: && sts[y][x]>0 */
	  && lbl_layer[y][x]<0) {
	//lbl_layer[y][x] = lbln;
	label_cc_boundary(x, y, lbln, lbl_edge, lbl_layer);
	lbln++;
      }
    }
  }

  if (save_smart_layers && 0==smart_it)
    debug_save_raster_int(lbl_layer, "saved_labels_first_smart_after_first_scan");

  // TODO: This is wrong, is makes too many labels, it needs a join of equivalent labels
  // (and the lengths, areas, wrscr are also wrong)



  // B) count boundary len, area, wrscr, etc.
  // Lengths will be overestimated after the number of labels is corrected (reduced). don't care
  boundary_lengths.resize(lbln);
  boundary_lengths.assign(lbln, 0);
  areas.resize(lbln);
  areas.assign(lbln, 0);
  wrscr.resize(lbln);
  wrscr.assign(lbln, 0.0f);
  // TODO: MISSING: cog_x, cog_y, etc. for warp_fast if this is also applied there

  boundary_len = 0;
  int current_lbl = 0;
  int last_equivalence = -1; // null equivalence   ---> TODO: this equivalence trick is insufficient (-what about connections that are below the lines scanned so far? (U shape blobs)????????????????????????????????????)
  //int equiv_replaces = -1;  // the (higher) label that is replaced by the (lower) equivalence label
  std::vector<int> equiv_replaces;
  equiv_replaces.assign(lbln, -1);

  for(size_t y=1; y<ydim-1; y++) {
    for(size_t x=1; x<xdim-1; x++) {
      if (sts[y][x] <= 0)
	continue;

      if (1==lbl_edge[y][x]) {
	// fix for inner islands:
	if (lbl_layer[y-1][x]>=0 && lbl_layer[y][x-1]>=0) {

	  /*
	  if (0==smart_it && x==154 && y==331)
	    Form1->Memo1->Lines->Add("******************** ================================= BEFORE (x, y): (154,331), current_lbl: "+IntToStr(current_lbl)+", lbl_layer[y][x]: "+IntToStr(lbl_layer[y][x]) +", last_equivalence: "+IntToStr(last_equivalence)+ ", equiv_replaces: "+IntToStr(equiv_replaces[lbl_layer[y][x]])+"===================================================, up: "+IntToStr(lbl_layer[y-1][x])+", left: "+IntToStr(lbl_layer[y-1][x]));
	  */

	  current_lbl = min(lbl_layer[y-1][x], lbl_layer[y][x-1]);	  

	  if (current_lbl < lbl_layer[y][x]) { last_equivalence = current_lbl; equiv_replaces[lbl_layer[y][x]] = current_lbl;}
	  else last_equivalence = -1;

	  /*
	  if (0==smart_it && x==154 && y==331)
	    Form1->Memo1->Lines->Add("******************** ================================= AFTER (x, y): (154,331), current_lbl: "+IntToStr(current_lbl)+", lbl_layer[y][x]: "+IntToStr(lbl_layer[y][x]) +", last_equivalence: "+IntToStr(last_equivalence)+ ", equiv_replaces: "+IntToStr(equiv_replaces[lbl_layer[y][x]])+"=================================================up: "+IntToStr(lbl_layer[y-1][x])+", left: "+IntToStr(lbl_layer[y-1][x]));
	  */

	} else if (lbl_layer[y-1][x]>=0) {
	  current_lbl = lbl_layer[y-1][x];

	  if (current_lbl < lbl_layer[y][x]) { last_equivalence = current_lbl; 	equiv_replaces[lbl_layer[y][x]] = current_lbl;}
	  else last_equivalence = -1;

	} else if (lbl_layer[y][x-1]>=0) {
	  current_lbl = lbl_layer[y][x-1];

	  if (current_lbl < lbl_layer[y][x]) { last_equivalence = current_lbl; equiv_replaces[lbl_layer[y][x]] = current_lbl;}
	  else last_equivalence = -1;

	}
	  // TODO: note down the merge
	else {

	  /*
	  if (0==smart_it && x==156 && y==331)
	    Form1->Memo1->Lines->Add("******************** ================================= (x, y): (156,331), lbl_layer[y][x]: "+IntToStr(lbl_layer[y][x]) +", last_equivalence: "+IntToStr(last_equivalence)+ ", equiv_replaces: "+IntToStr(equiv_replaces[lbl_layer[y][x]])+"===================================================");
	  */

	  //if (last_equivalence<0) {
	  if (equiv_replaces[lbl_layer[y][x]]<0) {
	    current_lbl = lbl_layer[y][x];
	  } else {
	    //if (equiv_replaces == lbl_layer[y][x])
	    if (equiv_replaces[lbl_layer[y][x]] >= 0)
	      //current_lbl = /*lbl_layer[y][x] =*/ last_equivalence;
	      current_lbl = equiv_replaces[lbl_layer[y][x]];
	    else
	      current_lbl = lbl_layer[y][x];
	  }
	}
      }
      //if (sts[y][x] > 0)  // writes inner (and also overwrite edge cells, don't care)
      lbl_layer[y][x] = current_lbl;

      areas[current_lbl]++;
      wrscr[current_lbl] += wrscr_mat[y][x];

      if (0==x)
	boundary_lengths[current_lbl]++;
      else if (0 >= sts[y][x-1])
	boundary_lengths[current_lbl]++;

      if (0==y)
	boundary_lengths[current_lbl]++;
      else if (0 >= sts[y-1][x])
	boundary_lengths[current_lbl]++;

      if(xdim-1 == x)
	boundary_lengths[current_lbl]++;
      else if (0 >= sts[y][x+1])
	boundary_lengths[current_lbl]++;

      if(ydim-1 == y)
	boundary_lengths[current_lbl]++;
      else if (0 >= sts[y+1][x])
	boundary_lengths[current_lbl]++;


      // TODO: MISSING:
      /*
      if (Fast_Warp == mode) {
	...
      }
      */
      
    }
  }

  // Now lbln has to be corrected by subtracting the number of equivalences:
  for(size_t i=0; i<equiv_replaces.size(); i++){
    if (equiv_replaces[i]>=0)
      lbln--;
  }


  if (save_smart_layers && 0==smart_it++)
    debug_save_raster_int(lbl_layer, "saved_labels_first_smart_after_second_scan");

  boundary_len = 0;
  for(std::vector<size_t>::iterator j=boundary_lengths.begin(); j != boundary_lengths.end(); j++)
    boundary_len += *j;

  if (0== ((smart_it++)%100))
    Form1->Memo1->Lines->Add("lbln: "+IntToStr(lbln) +", boundary len: "+IntToStr(boundary_len)+
			     +", area[0]: "+IntToStr(areas[0])+", wrscr[0]: "+FloatToStr(wrscr[0]));  

  /*
  time =float(clock()-t_begin)/CLOCKS_PER_SEC;
  Form1->Memo1->Lines->Add("Time fast: "+FloatToStr(time));
  */
  return lbln;
}
/////////////////////////////

int
Grid_CCSP::count_labels_as_if_removed(int x, int y, char**& sts, float delta, float removal_threshold)
{
  last_x = x;
  last_y = y;
  if (thickness <= 1)
    safe_xy_lbl = get_cc_label(x, y);
  else 
    safe_xy_lbl = uneroded_lbl_buffer[y][x];

  if (safe_xy_lbl <0)
    Form1->Memo1->Lines->Add("ERRRRRRRRR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, safe_xy_lbl, x: "+IntToStr(x)+
			     ", y: "+IntToStr(y)+", safe_xy_lbl: "+IntToStr(safe_xy_lbl)+
			     ", areas.size(): "+IntToStr(areas.size())+
			     ", edge: "+IntToStr(edge[y][x])+
			     ", edge_bkp_eroded: "+IntToStr(edge_bkp_eroded[y][x])+
			     ", status: "+IntToStr(status[y][x])+
			     ", status_bkp: "+IntToStr(status_bkp[y][x])
			     );


  safe_last_lbln = last_lbln;
  safe_last_last_lbln = last_last_lbln;
  if (thickness <= 1) {
    safe_areas = areas;
    safe_wrscr = wrscr;
    safe_boundary_lengths = boundary_lengths;
  } else {
    safe_areas = uneroded_areas;
    safe_wrscr = uneroded_wrscr;
    safe_boundary_lengths = uneroded_boundary_lengths;
  }

  // TODO - not used anywhere, kill it
  //safe_size = areas[get_cc_label(x, y)];
  if (thickness <= 1)
    safe_size = areas[safe_xy_lbl];
  else
    safe_size = uneroded_areas[safe_xy_lbl];
    

  int lbln = 0;
  if (true /*x>0 && y>0*/) {

    current_removal_required_ccl = true;

    // switch to alternative buffer
    size_t alt_buffer;
    if (0 == current_lbl_buffer)
      alt_buffer = 1;
    else if (1 == current_lbl_buffer)
      alt_buffer = 0;    
    current_lbl_buffer = alt_buffer;

    // Remove, count, un-remove (until accept_slow_removal()...)
    // sts[y][x] = 0; <<<----- this removal is now done in init_removal()
    if (thickness > 1 ) {


      // DO this before update_erode_foreground
      // here use status rather than sts (status_bkp)
      if (use_smart_count_labels) {
	count_smart_labels_n_bl_n_areas_n_wrscr(status, uneroded_lbl_buffer, wrscr_mat, 
						       uneroded_boundary_lengths, uneroded_areas, uneroded_wrscr);
      } else {
	int ln = uneroded_last_lbln = count_labels(status, uneroded_lbl_buffer);
	// why not status_bkp?
	calc_boundary_lengths_n_areas_n_wrscr(status, wrscr_mat, ln, uneroded_lbl_buffer, 
					      uneroded_boundary_lengths, uneroded_areas, uneroded_wrscr);
      }
      

      update_erode_foreground(x, y, thickness, sts, edge_bkp_eroded);
    }

    /*
    //lbln = count_labels(sts, lbl_layer);
    // Connected components labeling (expensive)
    lbln = count_labels(sts, lbl_buffers[current_lbl_buffer]);
    // needed below for last_unconfirmed_...
    calc_boundary_lengths_n_areas_n_wrscr(sts, wrscr_mat, lbln);
    // This will be done, if needed, by slow_cancel_removal: 
    // sts[y][x] = 1;
    */

    if (use_smart_count_labels) {
      //lbln = count_smart_labels_n_bl_n_areas_n_wrscr(status_bkp, lbl_buffers[current_lbl_buffer], wrscr_mat);
      lbln = count_smart_labels_n_bl_n_areas_n_wrscr(sts, lbl_buffers[current_lbl_buffer], wrscr_mat, 
						     boundary_lengths, areas, wrscr);
    } else {
      lbln = count_labels(sts, lbl_buffers[current_lbl_buffer]);
      // why not status_bkp?
      calc_boundary_lengths_n_areas_n_wrscr(sts, wrscr_mat, lbln, lbl_buffers[current_lbl_buffer],
					    boundary_lengths, areas, wrscr);
    }
    
    last_last_lbln = last_lbln;
    last_lbln = lbln;

    // This would get a (null) label from the new layer where the cell (x,y) has been removed!
    // int lbl = get_cc_label(x, y);
    int lbl = safe_xy_lbl;

    if (lbl <0)
      Form1->Memo1->Lines->Add("ERRRRRRRRR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, count_labels_as_if_removed (2), x: "+IntToStr(x)+
			     ", y: "+IntToStr(y)+", lbl: "+IntToStr(lbl)+
			     ", last_last_lbln: "+IntToStr(last_last_lbln)+"last_lbln: "+IntToStr(last_lbln)+
			     ", areas.size(): "+IntToStr(areas.size())+
			     ", edge: "+IntToStr(edge[y][x])+
			     ", edge_bkp_eroded: "+IntToStr(edge_bkp_eroded[y][x])+
			     ", status: "+IntToStr(status[y][x])+
			     ", status_bkp: "+IntToStr(status_bkp[y][x])
			     + ", n_registered_removals: "+IntToStr(n_registered_removals)
			     + ", safe_xy_lbl: "+IntToStr(safe_xy_lbl)
			     );
      /*Form1->Memo1->Lines->Add("AND 321, 31: "+IntToStr(get_cc_label(321, 31)));
    Form1->Memo1->Lines->Add("AND 321, 32: "+IntToStr(get_cc_label(321, 32)));
    Form1->Memo1->Lines->Add("AND 321, 25: "+IntToStr(get_cc_label(321, 25)));*/


    // /*  TODOOOOOOOOOO
    // Fill in cache
    int diff = last_lbln - last_last_lbln;
    //int diff = last_unconfirmed_removal.num_labels - lbln;
    if (cached && diff > 0) { // cache "critical" (split) values 
      //cached[y][x] = n_cc_lbl_calc;
      cached[y][x] = n_registered_removals;
      diff_cache[y][x] = diff;

      // p_cache!
      if (cached[y][x]>0)
	num_cached--;
      p_cache[y][x] = calc_penalty(diff, safe_areas[safe_xy_lbl], safe_wrscr[safe_xy_lbl]);

      num_cached++;
      num_ever_cached++;

      size_t max_size = 0, min_size=0, second_max_size=0;
      // LAST CHANGE: get_max_et_al(last_unconfirmed_removal.split_sizes, &max_size, &min_size, &second_max_size);
      get_max_et_al(last_unconfirmed_removal.split_sizes, max_size, min_size, second_max_size);

      //cached_factor[y][x] = 1.0f; //float(min_size) / float(max_size);
      // ********************************************************************************************************************************************************************************************************************************************************************
      //cached_factor[y][x] = (delta+p_cache[y][x] - removal_threshold) / (removal_threshold - CCSP_prev_prev_removal_threshold / 2.0f);	
      //cached_factor[y][x] = (delta+p_cache[y][x] - removal_threshold) / (removal_threshold - CCSP_prev_removal_threshold);	

      // TODO: clarify the misterious isolated points that get a huge cached_factor!!!!
      float tmp = (delta+p_cache[y][x] - removal_threshold) / (CCSP_ma_removal_rate);
      /*
      if (tmp>25.0f)
	Form1->Memo1->Lines->Add("SMALL MA, x: "+IntToStr(x)+", y: "+IntToStr(y)+", ma: "+FloatToStr(CCSP_ma_removal_rate));
      */
      cached_factor[y][x] = min(tmp, 25.0f);
      //cached_factor[y][x] = tmp; //(delta+p_cache[y][x] - removal_threshold) / (CCSP_ma_removal_rate);
    }
    // */

    // split happened
    if (diff > 0) {

      /* ************** IF FIRST TIME DIFF>=4 -> SAVE STATUS, STATUS_BKP, LBLS_STATUS, LBLS_STATUS_BKP ***************** */


      // Fill in unconfirmed removal state
      last_unconfirmed_removal.x = x;
      last_unconfirmed_removal.y = y;
      last_unconfirmed_removal.lbl = lbl;
      last_unconfirmed_removal.num_labels = lbln;
      last_unconfirmed_removal.diff = diff;

      if (thickness <= 1) {
	last_unconfirmed_removal.sizes = areas;
	last_unconfirmed_removal.boundary_lens = boundary_lengths;
	last_unconfirmed_removal.wrscr = wrscr;
      } else {
	last_unconfirmed_removal.sizes = uneroded_areas;
	last_unconfirmed_removal.boundary_lens = uneroded_boundary_lengths;
	last_unconfirmed_removal.wrscr = uneroded_wrscr;
      }
      //last_unconfirmed_removal.split_sizes.clear();
      //last_unconfirmed_removal.split_boundary_lens.clear();
      get_split_sizes_n_lens_n_wrscr(x, y, last_unconfirmed_removal.split_sizes, last_unconfirmed_removal.split_boundary_lens, 
				     last_unconfirmed_removal.split_wrscr);
      confirmed_removal_state = false;

    }

  } else {
    //last_size = 0;
    lbln = 0;
  }
  return lbln;
}

void
Grid_CCSP::discount_from_areas_wrscr_etc(int x,int y)
{
  //int lbln = get_cc_label(x, y);
  int lbln;
  if (thickness <= 1)
    lbln = get_cc_label(x, y);
  else 
    lbln = uneroded_lbl_buffer[y][x];


  areas[lbln]--;
  if (areas[lbln] <= 0) {   // only exactly == 0 should happen!
    // CRITICAL: otherwise, e.g., after removing 1 isolated cell, the
    // number of labels next time will be interpreted wrongly 
    wrscr[lbln] = 0;
    boundary_lengths[lbln] = 0;    
    last_lbln--;
    last_last_lbln = last_lbln;
    //last_removal_vanished_one_cc = true; ===> safe_last/last_last_lbln
	  /* 	  last_unconfirmed_removal.num_labels = last_lbln;	  */
  } else {
    //last_removal_vanished_one_cc = true;
    wrscr[lbln] -= wrscr_mat[y][x];
  }
  remaining_wrscr -= wrscr_mat[y][x];
}

void
Grid_CCSP::undiscount_from_areas_wrscr_etc(int x, int y, int lbl)
{
  // undo init removal!
  // normally the removal of an isolated cell should not need to be canceled, but who knows...
  if (areas[lbl] <= 0) {
    last_lbln++;
    last_last_lbln++;
  }
  areas[lbl]++;
  wrscr[lbl] += wrscr_mat[y][x];
  remaining_wrscr += wrscr_mat[y][x];
}

void
Grid_CCSP::slow_init_removal(size_t x, size_t y)
{
  if (thickness <= 1 ) {
    status_bkp[y][x] = 0;
    
    // MOVED FROM accept_slow_removal!!!
    
    // TODO: remove here? OR rather in accept_slow_removal - using safe_xy_lbl, safe_areas, etc.?   --- ADN THEN YOU DON'T NEED TO UNDO THIS REMOVAL IN CANCEL_REMOVAL
    //int lbln = get_cc_label(x, y);
    int lbln;
  if (thickness <= 1)
    lbln = get_cc_label(x, y);
  else 
    lbln = uneroded_lbl_buffer[y][x];


      if (lbln >= areas.size() || lbln < 0) {
	// would be WRONG: lbln = 0;
	return;
	Form1->Memo1->Lines->Add("ERRRRRRRRRRRRRRRRRRRRRRRRRRRROR, init_removal(), x: "+IntToStr(x)+
				 ", y: "+IntToStr(y)+", lbln: "+IntToStr(lbln)+
				 ", areas.size(): "+IntToStr(areas.size())+
			       ", edge: "+IntToStr(edge[y][x])+
			       ", edge_bkp_eroded: "+IntToStr(edge_bkp_eroded[y][x])+
			       ", status: "+IntToStr(status[y][x])+
			       ", status_bkp: "+IntToStr(status_bkp[y][x])
				 );
      } else {
	// Form1->Memo1->Lines->Add("OK, init_removal(), x: "+IntToStr(x)+
	// 			 ", y: "+IntToStr(y)+", lbln: "+IntToStr(lbln)+
	// 			 ", areas.size(): "+IntToStr(areas.size())+
	// 			 ", edge: "+IntToStr(edge[y][x])+
	// 		       ", edge_bkp_eroded: "+IntToStr(edge_bkp_eroded[y][x])+
	// 		       ", status: "+IntToStr(status[y][x])+
	// 		       ", status_bkp: "+IntToStr(status_bkp[y][x])
	// 			 );
      }

      discount_from_areas_wrscr_etc(x,y);


  } else { // thickness > 1

    /*
    TODO thickness > 1;

    decrease bl, area, wrscr in labels in diamond => to the diamond, apply 'discount_from_areas_wrscr_etc(x,y)'
    */
  }

  // for now...
  current_removal_required_ccl = false;
}

bool
Grid_CCSP::slow_cancel_removal(size_t x, size_t y)
{
  if (CCSP_verbose > 1)
    std::cerr << " ========== canceling, x: " << x << ", y: " << y << std::endl;

  // if TOO BIG penalty -> will not be removed
  // restore the last!
  last_lbln = safe_last_lbln;
  last_last_lbln = safe_last_last_lbln;
  areas = safe_areas;
  wrscr = safe_wrscr;
  boundary_lengths = safe_boundary_lengths;
  
  //status[y][x] = 1;
  status_bkp[y][x] = 1;
  
  undiscount_from_areas_wrscr_etc(x, y, safe_xy_lbl);
  //

  if (thickness <= 1) {

  } else { //  if (thickness > 1) {

    /*
    TODO thickness > 1;

    for every cell in the diamond: apply 'undiscount_from_areas_wrscr_etc(safe_xy_lbl);' --- arg, needs the safe_xy_lbl in the diamond ==> create a diamond_list_safe_xy_lbl!!!
    */

    /*    undo_update_erode_foreground(x, y, thickness, status_bkp, edge_bkp_eroded); */
    // This should be done more efficiently: last_last_unconfirmed...
    //calc_boundary_lengths_n_areas_n_wrscr(status_bkp, wrscr_mat);
  }


  if (current_removal_required_ccl) {
    // If canceling, and labels were recalculated for this tentative removal ->
    // go back to prev. labels. Beware, if this is the first cancellation, there's
    // no prev. labels layer to switch to.
    //if (!last_cached /*&& cancel_removal_count > 0*/) {
      if (1 == current_lbl_buffer)
	current_lbl_buffer = 0;
      else if (0 == current_lbl_buffer)
	current_lbl_buffer = 1;
      //}    
  }
  // wrong!
  //last_unconfirmed_removal = last_confirmed_removal;

  /*
  if (is_cached_penalty(x, y))
    recache_corridor(x, y);
  */
  cancel_removal_count++;
  return true;
}

//const bool CONN_THICKNESS = false;
const bool CONN_THICKNESS = true;

bool
Grid_CCSP::slow_accept_removal(int x, int y)
{
  confirmed_removal_state = true;
  n_registered_removals++;
  /*
  Form1->Memo1->Lines->Add("n_registered_removals++, x: "
			   +IntToStr(x)+", y: "+IntToStr(y));
  */
  
  if (CCSP_verbose>1)
    std::cerr << " ========== confirming, x: " << x << ", y: " << y << std::endl;
  status_bkp[y][x] = 0;  // shouldn't be needed  - is not needed for status[][] (removed in Remove_Site) but IT IS for status_bkp (ccsp internal stuff)
  
  /*
  // CRITICAL: otherwise, e.g., after removing 1 isolated cell, the
  // number of labels next time will be interpreted wrongly 
  int lbln = cc_label(x, y);
  // THIS IS WRONG. THIS CELL DOES NOT EXIST AFTER THE REMOVAL DONE IN COUNT_LABELS!!!
  areas[lbln]--;
  if (areas[lbln] <= 0) {
    wrscr[lbln] = 0;
    boundary_lengths[lbln] = 0;    
    last_lbln--;
    last_last_lbln = last_lbln;
    last_unconfirmed_removal.num_labels = last_lbln;
  } else {
    wrscr[lbln] -= wrscr_mat[y][x];
  }
  remaining_wrscr -= wrscr_mat[y][x];
  */


  last_confirmed_removal = last_unconfirmed_removal;
    

  if (is_cached_penalty(x, y))
    exterminate_corridor_cache(x, y);


  if (CONN_THICKNESS && thickness > 1) {
    int& rx = x;
    int& ry = y;

    edge_bkp_eroded[ry][rx] = 0;  // not needed in Slow mode?

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

  do_save_corridor_boundaries_layer();

  return true;
}

int before_a =0, after_a=0, after_b=0, after_c =0;

#include <iostream>
float
Grid_CCSP::calc_ccsp_as_if_removed(size_t x, size_t y, float delta, float removal_threshold)
{
  n_penalty_calc++;

  //
  if (0 == x || xdim-1 == x || 0 == y || ydim-1 == y) {
    return 0;
  }

  // this shouldn't happen!
  if (0 == status[y][x])
    return 0;


  int remaining = nonm1 - removed;

//#define XBLP_RANDOM
//#define XBLP_FAKE_BLP
//#define XBLP_SIMPLE_RULE
//#define XBLP_CCSP

//#ifdef XBLP_RANDOM

  if (0 == old_variant)
    return (randz()-0.5f)*boundary_len/remaining;
  //#elif defined(XBLP_FAKE_BLP)
  else if (1 == old_variant) {
    int diff[5] = {4, 2, 0, -2, -4};
    int nbc = diff[get_nb_count(status, x, y)];
    int BL = boundary_lengths[0]; //get_boundary_length();
    return CCSP*(BL/(float)remaining - (BL-nbc)/(remaining-1.0f));
  }
  //#elif defined(XBLP_SIMPLE_RULE)
  else if (2 == old_variant) {
    int maxdim = max(xdim, ydim);
    //return float(max_dist_around(x, y))/maxdim/remaining;
    return - CCSP * float(max_dist_around(x, y))/remaining;
    //return -0.03 * float(max_dist_around(x, y))/remaining;
  }
  //#elif defined(XBLP_CCSP)
  else if (3 == old_variant) {
    // proper CCSP
    
    float domain_coeff = calc_domain_coefficient(x, y);
    if (.0f == domain_coeff)
      return .0f;

    
    if (0==before_a%500)
      Form1->Memo1->Lines->Add("******before_a: "+IntToStr(before_a)+
			       ", after_a: "+IntToStr(after_a)+
			       ", after_b: "+IntToStr(after_b)+
			       ", after_c: "+IntToStr(after_c)
			       );


    char**& sts = status_bkp;
    if (CONN_THICKNESS && thickness > 1 ) {
      update_erode_foreground(x, y, thickness, sts, edge_bkp_eroded); 
    }


    before_a++;
    // A) Discard obvious cases (local alternative connection)
    // TODO: this is redundant with estimate_split_risk? (called from bat_run.cpp)
    bool risk = false;
    if (!CONN_THICKNESS || thickness <= 1) {
      risk = check_basic_split_risk(x, y, status_bkp);
    } else {
      risk = check_basic_split_risk_diamond(x, y, sts);
    }

    if (!risk) {
      if (CONN_THICKNESS && thickness > 1) {
	undo_update_erode_foreground(x, y, thickness, status_bkp, edge_bkp_eroded);
      }
      if (save_cached_layer)
	cached[y][x] = 0; //nonm1-removed;

      return 0;
    }

    after_a++;
    int lbln;
    size_t area;
    float richness;

    // B) check if cached
    if (thickness <= 1) {
      //lbln = get_cc_label(x, y);
      if (thickness <= 1)
	lbln = get_cc_label(x, y);
      else 
	lbln = uneroded_lbl_buffer[y][x];

      area = areas[lbln];
      richness = wrscr[lbln];
    } else {
      uneroded_last_lbln = lbln = uneroded_lbl_buffer[y][x];
      area = uneroded_areas[lbln];
      richness = uneroded_wrscr[lbln];
    }

  if (lbln<0)
    	Form1->Memo1->Lines->Add("ERRRRRR((((((((((((((((((((((((((((((((((((((((, calc_ccsp_as_if_removed(), x: "+IntToStr(x)+
				 ", y: "+IntToStr(y)+", lbln: "+IntToStr(lbln)+
				 ", areas.size(): "+IntToStr(areas.size())+
				 ", edge: "+IntToStr(edge[y][x])+
				 ", edge_bkp_eroded: "+IntToStr(edge_bkp_eroded[y][x])+
				 ", status: "+IntToStr(status[y][x])+
				 ", status_bkp: "+IntToStr(status_bkp[y][x])
				 );



    // cache:
    const float tolerance = .0f;  
    last_cached = false;
    //size_t cache_size = min(20*warp_ll, CACHE_LEN);  
    //if (cached && cached[y][x] > 0 && cached[y][x] > n_cc_lbl_calc - cache_size) {  //  /* && (cached[y][x] > n_cc_lbl_calc - CACHE_LEN) */) {
    //if (cached && cached[y][x] > 0 && (((cached[y][x] + cache_size) * cached_factor[y][x]) > n_cc_lbl_calc) ) {
    if (is_cached_penalty(x, y)) {
      int diff = diff_cache[y][x];
      float p = calc_penalty(diff, area, richness);

      // p_cache!
      p = p_cache[y][x];

      if (delta + p > removal_threshold /*+ tolerance*/) {
	// Will (should) not be removed -> cancel will not do anything
	last_cached = true;
	if (CCSP_verbose > 2)
	  std::cerr << " ++++++++++++++++++++++++++++ CACHE returning: " << p << std::endl;
	return p;
      }
    }

    after_b++;
    // /*  BIIIIIIIIIIIIIIIG TODO: why this doesn't work here
    // C) search alternative paths in the neighborhood  
    /*
    status_bkp[y][x] = 0;
    edge[y][x] = 0;
    bool no_alt_path = check_split_local_bfs(x, y, status_bkp);
    status_bkp[y][x] = 1;
    edge[y][x] = 1;
    if (no_alt_path)
      return 0.0f;
    */

    after_c++;
    // D) expensive map-wide operation
    // if (CONN_THICKNESS && thickness > 1) {  .. else {
    count_labels_as_if_removed(x, y, status_bkp, delta, removal_threshold);
    //cc_accounting->calc_boundary_lengths_n_areas_n_wrscr(status_bkp, wrscr_mat);

    //lbln = cc_label(x, y);
    // this will be fixed by slow_accept_removal(x,y) / slow_cancel_removal(x, y) - if called.
    //areas[lbln]++;

    int diff = last_lbln - last_last_lbln;
    //int diff = last_unconfirmed_removal.num_labels - last_lbln;
    
    if (CCSP_verbose > 2) {
      // the -1 incs here don't work if thickness>1
      std::cerr << ", AFTER SPLIT, x: " << x << ", y: " << y << ", DIFF: " << diff << std::endl;
      if (status[y][x-1] > 0 && get_cc_label(x-1,y) != 0)
	std::cerr << "   y,x-1: " << areas[get_cc_label(x-1,y)] << " (lbl: " << get_cc_label(x-1,y) << ")" << std::endl;
      if (status[y-1][x] > 0  && get_cc_label(x, y-1) != 0)
	std::cerr << "   y-1,x: " << areas[get_cc_label(x, y-1)] << " (lbl: " << get_cc_label(x, y-1) << ")" << std::endl;
      if (status[y][x+1] > 0  && get_cc_label(x+1, y) != 0)
	std::cerr << "   y,x+1: " << areas[get_cc_label(x+1, y)] << " (lbl: " << get_cc_label(x+1, y) << ")" << std::endl;
      if (status[y+1][x] > 0  && get_cc_label(x, y+1) != 0)
	std::cerr << "   y+1,x: " << areas[get_cc_label(x, y+1)] << " (lbl: " << get_cc_label(x, y+1) << ")" << std::endl;
    }


    size_t& rx = x;
    size_t& ry = y;

    if (CONN_THICKNESS && thickness > 1) {
      undo_update_erode_foreground(rx, ry, thickness, sts, edge_bkp_eroded);
    }

    if (!(diff >0)) {
      if (save_cached_layer)
	cached[ry][rx] = 0; //nonm1-removed;
      return 0;
    }
    if (save_cached_layer)
      cached[ry][rx] = removed; //nonm1-removed;


    if (diff > 0) {

      float p = calc_penalty(diff, area, richness);
      /*
      Form1->Memo1->Lines->Add("calc_penalty, x: "+IntToStr(x)
			       +", y: "+IntToStr(y)
			       +", diff: "+IntToStr(diff)
			       +", area: "+IntToStr(area)
			       +", rich: "+FloatToStr(richness)
			       +", p: "+FloatToStr(p)
			       );
      */

      if (CCSP_verbose > 2)
	std::cerr << " calc_ccsp_as_if_removed() - ++++++++++++++++++++++++++++ returning: " << p << std::endl;
      return domain_coeff * p;
    } else
      return 0;


  }
  //#elif defined(XBLP_MBLP)
  else if (4 == old_variant) {
    // multiple, component specific BL/A
    // This is done once per warp: calc_boundary_lengths_n_areas(status_bkp, wrscr_mat);
    int lbl; // = get_cc_label(x, y);
    if (thickness <= 1)
      lbl = get_cc_label(x, y);
    else 
      lbl = uneroded_lbl_buffer[y][x];


  if (lbl<0)
    	Form1->Memo1->Lines->Add("ERRRRRR===========================, 4==variant, x: "+IntToStr(x)+
				 ", y: "+IntToStr(y)+", lbln: "+IntToStr(lbl)+
				 ", areas.size(): "+IntToStr(areas.size())+
			       ", edge: "+IntToStr(edge[y][x])+
			       ", edge_bkp_eroded: "+IntToStr(edge_bkp_eroded[y][x])+
			       ", status: "+IntToStr(status[y][x])+
			       ", status_bkp: "+IntToStr(status_bkp[y][x])
				 );

    size_t bl = boundary_lengths[lbl];
    size_t area = areas[lbl];

    if (area >= 10) {
      
      int diff[5] = {4, 2, 0, -2, -4};
      //int diff[5] = {4, 2, 0, 0, 0};
      int nbc = diff[get_nb_count(status, x, y)];

      /*
      if (status[y-1][x-1]>0 && status[y+1][x+11])
	nbc -= 2;
      if (status[y+1][x-1]>0 && status[y-1][x+11])
	nbc -= 2;
      */
      
      nbc = diff_neighbors_8(x,y);

      if (0 == nbc)
	return 0;
      else {
	float p = CCSP*abs(sqrt(bl/(float)area - (bl-nbc)/(area-1.0f)));
	//float p = CCSP*sqrt(bl/(float)area - (bl-nbc)/(area-1.0f));
	// NOP float p = CCSP*z_pow(bl/(float)area - (bl-nbc)/(area-1.0f),2);

	p = CCSP*(bl/(float)area - (bl-nbc)/(area-1.0f));
	return p;
	if (p>0)
	  return p;
	else
	  return 0;
      }
    } else {
      return CCSP;
    }
  } else if (5 == old_variant) {
    int min_s = std::numeric_limits<int>::max();
    int max_s = std::numeric_limits<int>::min();
    for(size_t i=0; i<areas.size();i++) {
      if (areas[i] < min_s)
	min_s = areas[i];
      if (areas[i] > min_s)
	max_s = areas[i];
    }
    return CCSP * float(min_s)/float(max_s) / float(remaining);
  }

  //#endif
  return 0;
}
