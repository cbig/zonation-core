#include "ccsp.h"
// for status[][], etc.
#include "Unit1.h"
// for vmat[][]
#include "bat_run.h"
// for struct sp
#include "defines.h"
#include "GridMap.h"
#include "matrix_utils.h"
#include "randz.h"
#include "grid_utils.h"

#include <iostream>

// Speedup-URGENT TODO for "Slow" mode: fix the "TODO thickness" (remove and unremove with diamond)
//          in Grid_CCSP::cancel_removal and Grid_CCSP::init_removal
// THEN, for debugging: "IF FIRST TIME DIFF>=4 -> SAVE STATUS, STATUS_BKP, LBLS_STATUS, LBLS_STATUS_BKP"


extern bool use_smart_count_labels;
extern bool write_vdf;
extern bool save_cached_layer;

size_t Grid_CCSP::n_penalty_calc=0;

//const int CACHE_LEN = 5000;
//const int CACHE_LEN = 2000;
//const int CACHE_LEN = 1500;
int CACHE_LEN = 500;

void
Grid_CCSP::increase_cache_factor()
{
  global_cache_factor++;
}

void
Grid_CCSP::decrease_cache_factor()
{
  global_cache_factor--;
}

#include <boost/scoped_array.hpp>


void
Grid_CCSP::debug_save_raster(char**& m, String oname, float nodatavalue)
{
  // save raster
  //float nodatavalue(-1.0f);
  boost::scoped_array<float> plane(new float[xdim * ydim]);
  for(int y = 0; y < ydim; ++y) {
    for(int x = 0; x < xdim; ++x) {
      int i(x + y * xdim);
      if(m[y][x] < 0) {
	plane[i] = nodatavalue;
      } else {	
	plane[i] = m[y][x];
      }
    }
  }
  SaveToRaster<float>(oname, plane.get(), nodatavalue, xdim, ydim);  
}

void
Grid_CCSP::debug_save_raster_int(int**& m, String oname, int nodatavalue)
{
  // save raster
  //float nodatavalue(-1.0f);
  boost::scoped_array<float> plane(new float[xdim * ydim]);
  for(int y = 0; y < ydim; ++y) {
    for(int x = 0; x < xdim; ++x) {
      int i(x + y * xdim);
      if(m[y][x] < 0) {
	plane[i] = nodatavalue;
      } else {	
	plane[i] = m[y][x];
      }
    }
  }
  SaveToRaster<float>(oname, plane.get(), nodatavalue, xdim, ydim);  
}


Grid_CCSP::Grid_CCSP(float CCSP, int variant, int formula, int thickness, bool use_8_conn, String ofname, size_t xdim, size_t ydim, char**& status): 
  Grid_Connected_Components(xdim, ydim),
  CCSP(CCSP), rule_formula(formula), old_variant(3 /*variant*/), thickness(thickness), 
  use_8_conn(use_8_conn), thickness_inc(thickness / 2) /* round down */, xdim(xdim), ydim(ydim), 
  status(status), last_x(0), last_y(0), last_last_lbln(0), 
  uneroded_last_lbln(0),
  current_removal_required_ccl(false), safe_last_lbln(0), safe_xy_lbl(0), safe_last_last_lbln(0), last_lbln(0), /*last_size(0),*/ safe_size(0), 
  confirmed_removal_state(false), current_lbl_buffer(0), n_registered_removals(0), num_cached(0), num_ever_cached(0),
  // , boundary_lengths(0), areas(0), wrscr(0)
  cancel_removal_count(0), /*confirm_removal_count(0), */
  dim_bfs_minilayer(bfs_x_radius*2), global_cache_factor(5)
{
  /*
  if (thickness >1)
    this->thickness = 1;
  */

  //Form1->Memo1->Lines->Add("THICKNESS: "+IntToStr(thickness)+", INC: "+IntToStr(thickness_inc));

  if (1 == variant)
    this->mode = Grid_CCSP::Fast_Warp;
  else if (2 == variant)
    this->mode = Grid_CCSP::Slow;

  std::string path;
  if (!ofname.isEmpty()) {
    path = ofname.toUtf8().constData();
    path += ".corridors_log.txt";
    of = fopen(path.c_str(), "w");
    
    if(!of)
      return;
  }

  // TODO: only for verbose debug
  if (!ofname.isEmpty() && write_vdf) {
    vdf = fopen((path+"_zzz_corridors_verbose_log.txt").c_str(), "w");
    if (!vdf)
      return;
  }


  bfs_minilayer = cmatrix(0, dim_bfs_minilayer, 0, dim_bfs_minilayer);

  // Not needed, use lbl_layer (with int label numbers) as "visited" mask
  //if (use_smart_count_labels)
  //  layer_smart_lbl_visited = cmatrix(0, ydim, 0, xdim);
  
  //current_lbl_buffer = 0;
  lbl_buffers.resize(2);
  // lbl_buffers.assign(2, NULL);
  // TODO: memset!!!
  for(size_t buf=0; buf < lbl_buffers.size(); buf++) {
    lbl_buffers[buf] = imatrix(0, ydim, 0, xdim);
    for(size_t i = 0; i < ydim; i++) {	    
      for(size_t j = 0; j < xdim; j++) {
	(lbl_buffers[buf])[i][j] = -1; // critical to init lbl layers to whatever is expected as "undefined label" in count_labels
	//(lbl_buffers[buf])[i][j] = 0; // For Chang-C-L <0 / 0 / >0
      }
    }
    // for Chang-C-L
    //memset(&((lbl_buffers[buf])[0][0]), 0, xdim*ydim* sizeof(int));
    //memset(&((lbl_buffers[buf])[0][0]), -1, xdim*ydim* sizeof(int));
  }
  if (thickness > 1) {
    uneroded_lbl_buffer = imatrix(0, ydim, 0, xdim);
    for(size_t i = 0; i < ydim; i++) {	    
      for(size_t j = 0; j < xdim; j++) {
	uneroded_lbl_buffer[i][j] = -1;
      }
    }
  }

  if (Grid_CCSP::Slow == this->mode) {
    sf_lbl_buffer = imatrix(0, ydim, 0, xdim);
    for(size_t i = 0; i < ydim; i++) {	    
      for(size_t j = 0; j < xdim; j++) {
	sf_lbl_buffer[i][j] = -1;
      }
    }
  }


  status_bkp = cmatrix(0, ydim, 0, xdim);
  edge_bkp_eroded = cmatrix(0, ydim, 0, xdim);
  for(size_t i = 0; i < ydim; i++) {	    
    for(size_t j = 0; j < xdim; j++) {
      // Label status[][] or edge[][]?
      status_bkp[i][j] = status[i][j];
      //status_bkp[i][j] = edge[i][j];
    }
  }


  // Calculated as in output_grid_3
  wrscr_mat = fmatrix(0, ydim, 0, xdim);

  for(size_t y = 0; y < ydim; y++) {	    
    for(size_t x = 0; x < xdim; x++) {
      //if(-1 == Rmax[y][x]) {
      if(-1 == status[y][x]) {
	wrscr_mat[y][x] = -1;
      } else { 
	wrscr_mat[y][x] = 0;
	//Biodiv_Features_Occur_Container& rowp = vmat[y][x];
	//for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))
	for(size_t s = 0; s < map_cnt; s++) {
	  float v = vmat[y][x][s];   // beware of negative weights
	  if (v > .0f)
	    wrscr_mat[y][x] += v*spp[s].weight;
	}
      }
    }
  }

  remaining_wrscr = 0;
  for(size_t y = 0; y < ydim; y++) {	    
    for(size_t x = 0; x < xdim; x++) {
      if (wrscr_mat[y][x] > 0.0f)
	remaining_wrscr += wrscr_mat[y][x];
    }
  }

  if (thickness > 1) {

    // note +1: (thickness_inc+1)*2 <=== WRONG? - NOP
    state_backup_sts_bkp = cmatrix(0, thickness+1, 0, thickness+1);

    erosion_mask = cmatrix(0, ydim, 0, xdim);
    erosion_update_state_sts = cmatrix(0, thickness /*+1*/, 0, thickness /*+1*/);
    // eroded edge needs two more rows and cols (edge around removed cells)
    erosion_update_state_edge = cmatrix(0, thickness+2, 0, thickness+2);
    do_erode_foreground(status, status_bkp, thickness, edge_bkp_eroded); //do_thinning_foreground(thickness);

    if (Grid_CCSP::Slow == this->mode) {
      // Reconsider remaining wrscr (excluding corridor / boundary / buffer cells)
      remaining_wrscr = 0;
      for(size_t y = 0; y < ydim; y++) {	    
	for(size_t x = 0; x < xdim; x++) {
	  if (status_bkp[y][x] > 0 && wrscr_mat[y][x] > 0.0f)
	    remaining_wrscr += wrscr_mat[y][x];
	}
      }
    }


    bool save_all_here = false;
    bool save_eroded_sts = true;
    if (save_all_here && save_eroded_sts) {
      debug_save_raster(status_bkp, "saved_layer_status_bkp_eroded.1st", 0);
    }
    bool save_sts = true;
    if (save_all_here && save_sts) {
      debug_save_raster(status, "saved_layer_status.1st", 0);
    }
    bool save_edge = true;
    if (save_all_here && save_edge) {
      debug_save_raster(edge, "saved_layer_edge.1st", 0);
    }
    bool save_eroded_edge = true;
    if (save_all_here && save_eroded_edge) {
      debug_save_raster(edge_bkp_eroded, "saved_layer_edge_eroded.1st", 0);
    }
  }

  // Transfer current buffer to the alternative (or otherwise the first time it will switch to an empty (-1) labels layer):
  ///*
  size_t alt_buffer;
  if (0 == current_lbl_buffer)
    alt_buffer = 1;
  else if (1 == current_lbl_buffer)
    alt_buffer = 0;    
  for(size_t i = 0; i < ydim; i++) {	    
    for(size_t j = 0; j < xdim; j++) {
      if ((lbl_buffers[current_lbl_buffer])[i][j] >= 0) {
	(lbl_buffers[alt_buffer])[i][j] = (lbl_buffers[current_lbl_buffer])[i][j];
      }
    }
  }
  //*/

  //max_bfs_dim = 80;
  //local_bfs_mat = imatrix(0, max_bfs_dim, 0, max_bfs_dim);

  last_cached = false;
  cached = imatrix(0, ydim, 0, xdim);
  cached_factor = fmatrix(0, ydim, 0, xdim);
  for (size_t y=0; y<ydim; y++){
    for (size_t x=0; x<xdim; x++){
      cached[y][x] = -1;
      cached_factor[y][x] = .0f;
    }
  }

  diff_cache = imatrix(0, ydim, 0, xdim);
  p_cache = fmatrix(0, ydim, 0, xdim);
  // init to avoid use of unitialized values
  for (size_t y=0; y<ydim; y++){
    for (size_t x=0; x<xdim; x++){
      diff_cache[y][x] = 0;
      p_cache[y][x] = .0f;
    }
  }

  if (thickness > 1) {
    cache_local_bfs_diamond_split = cmatrix(0, ydim, 0, xdim);
    for (size_t y=0; y<ydim; y++){
      for (size_t x=0; x<xdim; x++){
	cache_local_bfs_diamond_split = 0;
      }
    }
  }

  //last_last_lbln = last_lbln = count_labels(-1, -1);
  //last_confirmed_removal.num_labels = last_last_lbln = last_lbln = count_labels(status_bkp, lbl_layer);
  if (use_smart_count_labels) {
    if (thickness > 1) {
	uneroded_last_lbln = count_smart_labels_n_bl_n_areas_n_wrscr(status, uneroded_lbl_buffer, wrscr_mat,
								     uneroded_boundary_lengths, uneroded_areas, uneroded_wrscr);
    }

    last_lbln = count_smart_labels_n_bl_n_areas_n_wrscr(status_bkp, lbl_buffers[current_lbl_buffer], wrscr_mat,
							boundary_lengths, areas, wrscr);
    debug_save_raster_int(lbl_buffers[current_lbl_buffer], "saved_labels_first_smart");

    int lll = count_labels(status_bkp, lbl_buffers[current_lbl_buffer]);
    debug_save_raster_int(lbl_buffers[current_lbl_buffer], "saved_labels_first_unionfind");
    Form1->Memo1->Lines->Add("******************** ================================= Consturctor, smart lbln: "+IntToStr(last_lbln) +", slow lbln: "+IntToStr(lll)+"===================================================");
      

  } else {
    if (thickness > 1) {
      // here status rather than status_bkp (eroded)
      int ln = uneroded_last_lbln = count_labels(status, uneroded_lbl_buffer);
      calc_boundary_lengths_n_areas_n_wrscr(status, wrscr_mat, ln, uneroded_lbl_buffer,
					    uneroded_boundary_lengths, uneroded_areas, uneroded_wrscr);
    }

    last_lbln = count_labels(status_bkp, lbl_buffers[current_lbl_buffer]);
    calc_boundary_lengths_n_areas_n_wrscr(status_bkp, wrscr_mat, last_lbln, lbl_buffers[current_lbl_buffer], 
					  boundary_lengths, areas, wrscr);

    bool save_layer_start = false;
    if (save_layer_start)
      debug_save_raster_int(lbl_buffers[current_lbl_buffer], "saved_labels_first_unionfind_traditional");

  }
  last_confirmed_removal.num_labels = safe_last_lbln = safe_last_last_lbln = 
    last_last_lbln = last_lbln;

  //std::sort(boundary_lengths.begin(), boundary_lengths.end(), std::greater<size_t>());
  //std::sort(areas.begin(), areas.end(), std::greater<size_t>());



  if (Grid_CCSP::Slow == this->mode) {
    if (thickness <= 1) { 
      sf_label_cnt_before = count_labels(status_bkp, sf_lbl_buffer);
      // TODO
      calc_boundary_lengths_n_areas_n_wrscr(status_bkp, wrscr_mat, sf_label_cnt_before, sf_lbl_buffer,
					    sf_before_boundary_lengths, sf_before_areas, sf_before_wrscr);

    } else /*if (thickness > 1)*/ {

      //std::cerr << "CONSTRUCTOR ********* POS 407 7: " << (int)status_bkp[7][407] << std::endl;
      sf_label_cnt_before = count_labels(status_bkp, sf_lbl_buffer);
      calc_boundary_lengths_n_areas_n_wrscr(status_bkp, wrscr_mat, sf_label_cnt_before, sf_lbl_buffer,
					    sf_before_boundary_lengths, sf_before_areas, sf_before_wrscr);
    }
    // sf_label_before will be init later, in calc_ccsp
  }



  if (false && of) {
    fprintf(of, "# CCSP debug/output file. Format:\n");
    fprintf(of, "# fract_lost, cells_removed, edge_count,  calculation_rounds_penalty, calculations_connected_components, lbl_count, lbl_alg_equivs, 6_biggest_boundary_lengths, 6_biggest_areas...\n" );
    output_info_line(0.0);
  /*
    fprintf(of, "%-9.6f, %d, %d, %d, %d, ", 0.0, 1, ecnt, last_lbln, last_lbl_equiv);
    for(size_t i=0; i<5 && i<boundary_lengths.size(); i++) {
    fprintf(of, "%d, ", boundary_lengths[i]);
    }
    fprintf(of,"\n");
  */
  }


  if (Grid_CCSP::Slow == this->mode) {
    sf_status_minilayer = cmatrix(0, 1+2*(thickness-1), 0, 1+2*(thickness-1));
  }
}

Grid_CCSP::~Grid_CCSP()
{
  //free_imatrix(status_bkp, 0, ydim, 0, xdim);
  //free_imatrix(lbl_layer, 0, ydim, 0, xdim);
  for (size_t b=0; b < lbl_buffers.size(); b++){
    //free_imatrix(lbl_buffers[b], 0, ydim, 0, xdim);
    lbl_buffers[b] = NULL;
  }
  if (thickness > 1) {
    free_imatrix(uneroded_lbl_buffer, 0, ydim, 0, xdim);
    free_cmatrix(cache_local_bfs_diamond_split, 0, ydim, 0, xdim);
  }

  free_cmatrix(erosion_mask, 0, ydim, 0, xdim);
  free_cmatrix(erosion_update_state_sts, 0, thickness /*+1*/, 0, thickness /*+1*/);
  free_cmatrix(state_backup_sts_bkp, 0, thickness+1, 0, thickness+1);
  free_cmatrix(erosion_update_state_edge, 0, thickness+2, 0, thickness+2);

  free_cmatrix(bfs_minilayer, 0, dim_bfs_minilayer, 0, dim_bfs_minilayer);
  //free_cmatrix(layer_smart_lbl_visited, 0, ydim, 0, xdim);

  if (Grid_CCSP::Slow == this->mode) {
    free_imatrix(sf_lbl_buffer, 0, ydim, 0, xdim);
    free_cmatrix(sf_status_minilayer, 0, 1+2*(thickness-1), 0, 1+2*(thickness-1));
  }

  fclose(of);
}

bool
Grid_CCSP::is_cached_penalty(int x, int y)
{
  // cache_size 4*warp_ll ---> CCSP=0.1 -> 300s
  // cache_size 4*warp_ll ---> CCSP=1 -> 423s
  // size_t cache_size = 10 * warp_ll; AFTER cached[][]=-1 corridor removal acceleration: CCSP=0.1 -> 207s
  // size_t cache_size = 10 * warp_ll; AFTER cached[][]=-1 corridor removal acceleration: CCSP=10 -> 233s
  // size_t cache_size = 10 * warp_ll; AFTER cached[][]=-1 corridor removal acceleration: CCSP=1000 -> x195s
  //size_t cache_size = 10 * warp_ll;

  bool adaptive_cache_rate = false;
  if (!adaptive_cache_rate) {
    global_cache_factor = 2;
    //int f = max((int)global_cache_factor, 1);
    //size_t cache_size = global_cache_factor * warp_ll;

    int warp_rnds(min(warp_factor, ecnt/100)); // xxxTBF fix

    /* TODO: global factor set to 1 here */
    global_cache_factor = 5;
    size_t cache_size = global_cache_factor * warp_rnds;

    //return (cached && cached[y][x] > 0 && cached[y][x] > n_cc_lbl_calc - cache_size);
    
    //if (cached && cached[y][x] > 0 && cached[y][x] > n_cc_lbl_calc - cache_size) {  //  /* && (cached[y][x] > n_cc_lbl_calc - CACHE_LEN) */) {
    //return (cached && cached[y][x] > 0 && (((cached[y][x] + cache_size) * cached_factor[y][x]) > n_cc_lbl_calc) );
    if (cached && cached[y][x] > 0 ) {
      //if (cached && cached[y][x] > 0 && (((cached[y][x] + cache_size) * cached_factor[y][x]) > n_registered_removals) ) {
      //if ((float(cached[y][x]) + float(warp_rnds) * cached_factor[y][x]) > float(n_registered_removals)) {
      if ((float(cached[y][x]) + cache_size) > float(n_registered_removals)) {
	return true;
      } else {
	cached[y][x] = -1;
	num_cached--;
      }
    }
  } else {
    

  }
  return false; 
}

//std::vector<size_t>*
void
Grid_CCSP::calc_boundary_lengths_n_areas_n_wrscr(char**& sts, float**& wrscr_mat, int last_lbln, 
						 int**& lbl_buffer,
						 std::vector<size_t>& boundary_lengths,
						 std::vector<size_t>& areas,
						 std::vector<float>& wrscr)
{
  // ***** std::cerr << "0!: " << xdim << ", " << ydim << ": " << last_lbln;
  if (last_lbln <= 0)
    return;

  //  Form1->Memo1->Lines->Add("CALC 0: "+IntToStr(last_lbln));
  //std::cerr << "0: (" << last_lbln << ") ";
  boundary_lengths.resize(last_lbln);
  //Form1->Memo1->Lines->Add("CALC 01: "+IntToStr(last_lbln));
  //std::cerr << "a ";
  boundary_lengths.assign(last_lbln, 0);
  //Form1->Memo1->Lines->Add("CALC 02: "+IntToStr(last_lbln));
  //std::cerr << "1 ";
  areas.resize(last_lbln);
  //Form1->Memo1->Lines->Add("CALC 03: "+IntToStr(last_lbln));
  //std::cerr << "a ";
  areas.assign(last_lbln, 0);
  //std::cerr << "2 ";
  wrscr.resize(last_lbln);
  //std::cerr << "a ";
  wrscr.assign(last_lbln, 0.0);
  //std::cerr << "3 ";

  //Form1->Memo1->Lines->Add("CALC A: "+IntToStr(last_lbln));

  cog_x.resize(last_lbln); cog_y.resize(last_lbln); cog_x_prev.resize(last_lbln); cog_y_prev.resize(last_lbln); cog_n.resize(last_lbln); cog_sum_weights.resize(last_lbln); cog_std_x.resize(last_lbln); cog_std_y.resize(last_lbln); //cog_std_x_prev.resize(last_lbln); cog_std_y_prev.resize(last_lbln);
  cog_x.assign(last_lbln, 0.0); cog_y.assign(last_lbln, 0.0); cog_x_prev.assign(last_lbln, 0.0); cog_y_prev.assign(last_lbln, 0.0); cog_n.assign(last_lbln, 0.0); cog_sum_weights.assign(last_lbln, 0.0); cog_std_x.assign(last_lbln, 0.0); cog_std_y.assign(last_lbln, 0.0); //cog_std_x_prev.assign(last_lbln, 0.0); cog_std_y_prev.assign(last_lbln, 0.0);

  //Form1->Memo1->Lines->Add("CALC B: "+IntToStr(last_lbln));

  // size_t bl = 0;  // global length!
  boundary_len = 0;
  for(size_t y=0; y<ydim; y++) {
    for(size_t x=0; x<xdim; x++) {
      if (0 >= sts[y][x])
	continue;

      int lbl = lbl_buffer[y][x];
      
      if (lbl < 0)
	Form1->Memo1->Lines->Add("ERRRR}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}, calc_boundary_lengths_n_...(), x: "+IntToStr(x)+
				 ", y: "+IntToStr(y)+", lbln: "+IntToStr(lbl)+
				 ", areas.size(): "+IntToStr(areas.size())+
				 ", edge: "+IntToStr(edge[y][x])+
				 ", edge_bkp_eroded: "+IntToStr(edge_bkp_eroded[y][x])+
				 ", status: "+IntToStr(status[y][x])+
				 ", status_bkp: "+IntToStr(status_bkp[y][x])
				 );


      /*
      if(lbl >= last_lbln || lbl<0)  
	Form1->Memo1->Lines->Add("lbl: "+IntToStr(lbl)+"(and lbln: "+IntToStr(last_lbln));
      */

      areas[lbl]++;
      wrscr[lbl] += wrscr_mat[y][x];

      if (0==x)
	boundary_lengths[lbl]++;
      else if (0 >= sts[y][x-1])
	boundary_lengths[lbl]++;

      if (0==y)
	boundary_lengths[lbl]++;
      else if (0 >= sts[y-1][x])
	boundary_lengths[lbl]++;

      if(xdim-1 == x)
	boundary_lengths[lbl]++;
      else if (0 >= sts[y][x+1])
	boundary_lengths[lbl]++;	

      if(ydim-1 == y)
	boundary_lengths[lbl]++;
      else if (0 >= sts[y+1][x])
	boundary_lengths[lbl]++;	

      if (Fast_Warp == mode || true) {
	// http://math.stackexchange.com/questions/102978/incremental-computation-of-standard-deviation
	// *** http://mathcentral.uregina.ca/QQ/database/QQ.09.02/carlos1.html
	// and the weighted incremental algorithm: http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
	float n = ++cog_n[lbl];
	if (1==n) {
	  cog_x[lbl] = x;
	  cog_y[lbl] = y;
	  cog_std_x[lbl] = 0;
	  cog_std_y[lbl] = 0;
	  if (1 <= rule_formula && rule_formula <= 4)
	    cog_sum_weights[lbl] = wrscr_mat[y][x];
	} else {
	  if (11 <= rule_formula && rule_formula <= 14) {
	    cog_x_prev[lbl] = cog_x[lbl];
	    cog_y_prev[lbl] = cog_y[lbl];
	    cog_x[lbl] = cog_x_prev[lbl] + (x - cog_x[lbl]) / n;
	    cog_y[lbl] = cog_y_prev[lbl] + (y - cog_y[lbl]) / n;
	    
	    //cog_std_x_prev[lbl] = cog_std_x[lbl];
	    //cog_std_y_prev[lbl] = cog_std_y[lbl];
	    //cog_std_x[lbl] = cog_std_x_prev[lbl] + (x - cog_x_prev[lbl])*(x-cog_x[lbl]);
	    //cog_std_y[lbl] = cog_std_y_prev[lbl] + (y - cog_y_prev[lbl])*(y-cog_y[lbl]);
	    cog_std_x[lbl] = cog_std_x[lbl] + (x - cog_x_prev[lbl])*(x-cog_x[lbl]);
	    cog_std_y[lbl] = cog_std_y[lbl] + (y - cog_y_prev[lbl])*(y-cog_y[lbl]);
	  } else { // (1 <= rule_formula && rule_formula <= 4)
	    float w = wrscr_mat[y][x];
	    float tmp = cog_sum_weights[lbl] + w;
	    cog_x_prev[lbl] = cog_x[lbl];
	    cog_y_prev[lbl] = cog_y[lbl];
	    cog_x[lbl] = cog_x_prev[lbl] + w*(x - cog_x[lbl]) / tmp;
	    cog_y[lbl] = cog_y_prev[lbl] + w*(y - cog_y[lbl]) / tmp;

	    //cog_std_x_prev[lbl] = cog_std_x[lbl];
	    //cog_std_y_prev[lbl] = cog_std_y[lbl];
	    //cog_std_x[lbl] = cog_std_x_prev[lbl] + w*(x - cog_x_prev[lbl])*(x-cog_x[lbl]);
	    //cog_std_y[lbl] = cog_std_y_prev[lbl] + w*(y - cog_y_prev[lbl])*(y-cog_y[lbl]);

	    cog_std_x[lbl] = cog_std_x[lbl] + w*(x - cog_x_prev[lbl])*(x-cog_x[lbl]);
	    //cog_std_x[lbl] = cog_std_x[lbl] + cog_sum_weights[lbl]*(x - cog_x_prev[lbl])* (x-cog_x_prev[lbl])*w/tmp;
	    cog_std_y[lbl] = cog_std_y[lbl] + w*(y - cog_y_prev[lbl])*(y-cog_y[lbl]);
	    cog_sum_weights[lbl] = tmp;
	  }
	}
      }

    }
  }
  //  Form1->Memo1->Lines->Add("CALC C: "+IntToStr(last_lbln));
  //std::cerr << "4 " << std::endl;
  // TODO: remove this from here, reuse old stuff!
  for(std::vector<size_t>::iterator j=boundary_lengths.begin(); j != boundary_lengths.end(); j++)
    boundary_len += *j;
  // Form1->Memo1->Lines->Add("CALC D: "+IntToStr(last_lbln));
}

bool
Grid_CCSP::output_info_line(float prop_lost)
{
  // gather status info and fprintf
  //fprintf(of, "%-9.6f %9d %9d\n", prop_lost, n_penalty_calc, ecnt);
  return false;
  if (!of)
    return false;

  fprintf(of, "%-9.6f, %9d, %9d,  %8d, %6d, %6d, %6d,   ", prop_lost, current-1, ecnt, n_penalty_calc, this->calc_count(), last_lbln, -1 /*last_lbl_equiv*/);
  for(size_t i=0; i<6 && i<boundary_lengths.size(); i++) {
    fprintf(of, "%d, ", boundary_lengths[i]);
  }
  for(size_t i=0; i<6 && i<areas.size(); i++) {
    fprintf(of, "%d, ", areas[i]);
  }
  fprintf(of,"\n");
}

int
Grid_CCSP::max_dist_around(size_t x, size_t y)
{
  int dx = 0;
  // E-W
  while (x-dx>0 && x+dx<xdim-1 ) {
    if (0 >= status[y][x-dx] || 0 >= status[y][x+dx])
      break;
    dx++;
  }

  // N-S
  int dy = 0;
  while (y-dy>0 && y+dy<ydim-1) {
    if (0 >= status[y-dy][x] || 0 >= status[y+dy][x])
      break;
    dy++;    
  }

  // NW-SE
  int dd1 = 0;
  while (x-dd1>0 && x+dd1<xdim-1 && y-dd1>0 && y+dd1<ydim-1) {
    if (0 >= status[y-dd1][x-dd1] || 0 >= status[y+dd1][x+dd1])
      break;
    dd1++;    
  }

  // NE-SW
  int dd2 = 0;
  while (x-dd2>0 && x+dd2<xdim-1 && y-dd2>0 && y+dd2<ydim-1) {
    if (0 >= status[y+dd2][x-dd2] || 0 >= status[y-dd2][x+dd2])
      break;
    dd2++;
  }
  return sqrt(dx*dx + dy*dy + dd1*dd1 + dd2*dd2);
}




bool
Grid_CCSP::exterminate_corridor_cache(size_t x, size_t y)
{
  // accelerate removal of neighbors which may have been part of a corridor but most likely are not anymore...
  // assumes 4-neighborhood
  // std::cerr << "cacheds: (" << x << "," << y << "): " << cached[y][x-1] << " " << cached[y-1][x] << " " << cached[y][x+1] << " " << cached[y+1][x] << "   ";
  if (cached[y][x-1]>0) {
    cached[y][x-1] = -1;
    //cached_factor[y][x-1] = 1.0f;
    num_cached -= 1;
    if (x-1 > 1)
      exterminate_corridor_cache(x-1, y);
  }
  if (cached[y-1][x]>0) {
    cached[y-1][x] = -1;
    //cached_factor[y-1][x] = 1.0f;
    num_cached -= 1;
    if (y-1 > 1)
      exterminate_corridor_cache(x, y-1);
  }
  if (cached[y][x+1]>0) {
    cached[y][x+1] = -1;
    //cached_factor[y][x+1] = 1.0f;
    num_cached -= 1;
    if (x+1 < xdim-1)
      exterminate_corridor_cache(x+1, y);
  }
  if (cached[y+1][x]>0) {
    cached[y+1][x] = -1;
    //cached_factor[y+1][x] = 1.0f;
    num_cached -= 1;
    if (y+1 < ydim-1)
      exterminate_corridor_cache(x, y+1);
  }
  return true;
}

bool
Grid_CCSP::recache_corridor(size_t x, size_t y)
{
  // If a cell removal has been canceled because it is still a valuable corridor, update the cache life of all the connected corridor cells
  if (CCSP_verbose > 2)
    std::cerr << "cacheds: (" << x << "," << y << "): " << cached[y][x-1] << " " << cached[y-1][x] << " " << cached[y][x+1] << " " << cached[y+1][x] << "   ";

  if (cached[y][x-1] > 0 && cached[y][x-1] < n_registered_removals) {
    cached[y][x-1] = n_registered_removals;
    //cached_factor[y][x-1] = 1.0f;
    if (x-1 > 1)
      recache_corridor(x-1, y);
  }
  if (cached[y-1][x] > 0 && cached[y-1][x] < n_registered_removals) {
    cached[y-1][x] = n_registered_removals;
    //cached_factor[y-1][x] = 1.0f;
    if (y-1 > 1)
      recache_corridor(x, y-1);
  }
  if (cached[y][x+1] > 0 && cached[y][x+1] < n_registered_removals) {
    cached[y][x+1] = n_registered_removals;
    //cached_factor[y][x+1] = 1.0f;
    if (x+1 < xdim-1)
      recache_corridor(x+1, y);
  }
  if (cached[y+1][x] > 0 && cached[y+1][x] < n_registered_removals) {
    cached[y+1][x] = n_registered_removals;
    //cached_factor[y+1][x] = 1.0f;
    if (y+1 < ydim-1)
      recache_corridor(x, y+1);
  }
  return true;
}








// just to discard disconnected/isolated cells
bool
Grid_CCSP::check_basic_split_risk(int rx, int ry, char**& sts)
{  
  //risk = (get_nb_count(sts, rx, ry) > 1 ) && get_nb8_count(sts, rx, ry) < 7;
  
  // /* would work only for thickness==1
  // In case "edge removal" if off
  /*  if (1 != edge[ry][rx])
      return false;*/
  // */

  int nb4 = get_nb_count(sts, rx, ry);
  int nb8 = nb4 + get_nb_corners_count(sts, rx, ry);

  bool risk = (nb4 > 1) && (nb8 < 7);
  return risk;
}

bool
Grid_CCSP::check_basic_split_risk_diamond(int rx, int ry, char**& sts)
{  
  bool risk = false;
  // 4-star
  //int inc = thickness_inc + 1; // would be wrong!
  int& inc = thickness_inc;

  /*
  if (rx-inc > 0 && sts[ry][rx-inc]>0)
    risk |= check_basic_split_risk(rx-inc, ry, sts);
  if (!risk && (rx+inc < xdim) && sts[ry][rx+inc] >0)
    risk |= check_basic_split_risk(rx+inc, ry, sts);
  if (!risk && (ry-inc > 0) && sts[ry-inc][rx]>0)
    risk |= check_basic_split_risk(rx, ry-inc, sts);
  if (!risk && (ry+inc < ydim) && sts[ry+inc][rx]>0)
  risk |= check_basic_split_risk(rx, ry+inc, sts);
  */


  /*
  // OLD SIMPLE Star
  if (rx>0)
    risk |= check_basic_split_risk(rx-inc, ry, sts);
  if (!risk && ry>0)
    risk |= check_basic_split_risk(rx, ry-inc, sts);
  if (!risk && rx<xdim-1)
    risk |= check_basic_split_risk(rx+inc, ry, sts);
  if (!risk && ry<ydim-1)
    risk |= check_basic_split_risk(rx, ry+inc, sts);
  
  return risk;
  */

  /* All this is not needed anymore with (eroded) status_bkp[][]... it would only apply to the normal status[][]
  // A) with thickness > 1, there must be a thickness-long row either vertically or horizontally
  // both horiz_len and vert_len count *contiguous* cells
  // upwards
  int vert_len = 0;
  for(int j=ry-1; j>=ymin; j--) {
    if (sts[j][rx])
      vert_len++;
    else
      break;
  }
  vert_len++;
  // downwards
  for(int j=ry+1; j<=ymax; j++) {
    if (sts[j][rx])
      vert_len++;
    else
      break;
  }
  if (vert_len < thickness) {
    return false;
  }
  int horiz_len = 0;
  // to left
  for (int i=rx-1; i>=xmin; i--) {
    if (sts[ry][i])
      horiz_len++;
    else
      break;
  }
  horiz_len++;
  // to right
  for (int i=rx+1; i<=xmax; i++) {
    if (sts[ry][i])
      horiz_len++;
    else
      break;
  }
  if (horiz_len < thickness) {
    return false;
    }
  */

  // Check (edge of) the square around
  int xmin = max(1, rx-inc);
  int xmax = min((int)xdim-2, rx+inc);
  int ymin = max(1, ry-inc);
  int ymax = min((int)ydim-2, ry+inc);
  // A bit faster (2-3% total run time save) below...
  // for(int j=ymin; j<=ymax; j++) {
  //   for(int i=xmin; i<=xmax; i++) {
  //     if((j==ymin || j==ymax || i==xmin || i==xmax) /* && edge_bkp_eroded[j][i]*/) // ??? edge_bkp_eroded /OR status_bkp
  // 	risk |= check_basic_split_risk(i, j, sts);
  //   }
  // }

  // top row
  for (int i=xmin; i<=xmax; i++) {
    risk |= check_basic_split_risk(i, ymin, sts);
    if (risk)
      return risk;
  }
  // bottom row
  for (int i=xmin; i<=xmax; i++) {
    risk |= check_basic_split_risk(i, ymax, sts);
    if (risk)
      return risk;
  }
  // left (excluding corners)
  for(int j=ymin+1; j<ymax; j++) {
    risk |= check_basic_split_risk(xmin, j, sts);
    if (risk)
      return risk;
  }
  // right (excluding corners)
  for(int j=ymin+1; j<ymax; j++) {
    risk |= check_basic_split_risk(xmax, j, sts);
    if (risk)
      return risk;
  }
  
  /*
    // old diamond

  // 0 all cells on the diamond  
  for(int j=ry; j<=ymax; j++) {
    for(int i=xmin; i<=rx; i++) {
      sts[j][i] = 0;
    }
  }
  for(int j=ymax-1; j>=ry; j--) {
    for(int i=rx+1; i<=xmax; i++) {
      sts[j][i] = 0;
    }
  }
  for(int j=ry-1; j>=ymin; j--) {
    for(int i=xmax-1; i>=rx; i--) {
      sts[j][i] = 0;
    }
  }
  for(int j=ymin+1; j<=ry-1; j++) {
    for(int i=rx-1; i>=xmin+1; i--) {
      sts[j][i] = 0;
    }
  }

  // [W->S]
  for(int j=ry; j<=ymax; j++) {
    for(int i=xmin; i<=rx; i++) {
      if (sts[j][i] >= 0 && edge_bkp_eroded[j][i] ) { // works if the bfs search is done on edge_bkp_eroded
	risk |= check_basic_split_risk(i, j, sts);
      }
    }
  }
  // (S->E]
  for(int j=ymax-1; j>=ry; j--) {
    for(int i=rx+1; i<=xmax; i++) {
      if (sts[j][i] >= 0 && edge_bkp_eroded[j][i]) {
	risk |= check_basic_split_risk(i, j, sts);
      }
    }
  }
  // (E->N]
  for(int j=ry-1; j>=ymin; j--) {
    for(int i=xmax-1; i>=rx; i--) {
      if (sts[j][i] >= 0 && edge_bkp_eroded[j][i]) {
	risk |= check_basic_split_risk(i, j, sts);
      }
    }
  }
  // (N->W)
  for(int j=ymin+1; j<=ry-1; j++) {
    for(int i=rx-1; i>=xmin+1; i--) {
      if (sts[j][i] >= 0 && edge_bkp_eroded[j][i]) {
	risk |= check_basic_split_risk(i, j, sts);
      }
    }
  }

  // Restore state
  for (int j=ymin; j<=ymax; j++) {
    for (int i=xmin; i<=xmax; i++) {    
      sts[j][i] = state_backup_sts_bkp[j-ymin][i-xmin];
    }
  }

  */

  return risk;
}






bool
Grid_CCSP::estimate_split_risk(char**& sts, int rx, int ry, bool try_local_search)
{
  // when using thickness 
  // TODO THICKNESS: This should be like check_basic_split_risk / diamond and  a bit of check_split_local_bfs_diamond
  if (thickness > 1)
    return true;


  bool risk = true;
  if (!use_8_connectivity) {
    risk = (get_nb_count(sts, rx, ry) > 1 ) && get_nb8_count(sts, rx, ry) < 7;
    //risk = (get_nb_count(sts, rx, ry) > 1 ) && get_nb8_count(sts, rx, ry) <=6 ;  ==> does not work in 4-conn, because:
    // 1 1 0
    // 1 X 1
    // 1 1 0
    bool horiz_split = ((1 == sts[ry][rx-1] && 1 == sts[ry][rx+1])  && !(1==sts[ry-1][rx-1] && 1==sts[ry-1][rx] && 1==sts[ry-1][rx+1]) && !(1==sts[ry+1][rx-1] && 1==sts[ry+1][rx] && 1==sts[ry+1][rx+1]) );
    bool vert_split = ((1 == sts[ry-1][rx] && 1 == sts[ry+1][rx]) && !(1==sts[ry-1][rx-1] && 1==sts[ry][rx-1] && 1==sts[ry+1][rx-1]) && !(1==sts[ry-1][rx+1] && 1==sts[ry][rx+1] && 1==sts[ry+1][rx+1]) );
    bool corner1_split = ((1 == sts[ry][rx-1] && 1 == sts[ry-1][rx]) && !(1==sts[ry-1][rx-1]) );
    bool corner2_split = ((1 == sts[ry-1][rx] && 1 == sts[ry][rx+1]) && !(1==sts[ry-1][rx+1]) );
    bool corner3_split = ((1 == sts[ry][rx+1] && 1 == sts[ry+1][rx]) && !(1==sts[ry+1][rx+1]) );
    bool corner4_split = ((1 == sts[ry+1][rx] && 1 == sts[ry][rx-1]) && !(1==sts[ry+1][rx-1]) );
    risk &= (horiz_split || vert_split || 
	     corner1_split || corner1_split || corner2_split || corner3_split || corner4_split);

    return risk;

  } else {
    risk = (get_nb_count(sts, rx, ry) > 1 ) && get_nb8_count(sts, rx, ry) < 7;

    risk &= (
	     // horizontal and vertical lines across center
	     ( (1 == sts[ry-1][rx] && 1 == sts[ry+1][rx]) && !(1 == sts[ry][rx-1]
								     ||
								     1 == sts[ry][rx+1]
								     ) )
	     ||
	     ( (1 == sts[ry][rx-1] && 1 == sts[ry][rx+1]) && !(1 == sts[ry-1][rx]
								     ||
								     1 == sts[ry+1][rx]
								     ) )
	     ||
	     // diagonal lines across center
	     ( (1 == sts[ry-1][rx-1] && 1 == sts[ry+1][rx+1]) && !(1 == sts[ry-1][rx] && 1 == sts[ry][rx+1]
									 ||
									 1 == sts[ry][rx-1] && 1 == sts[ry+1][rx]
									 ) )
	     ||
	     ( (1 == sts[ry-1][rx+1] && 1 == sts[ry+1][rx-1]) && !(1 == sts[ry-1][rx] && 1 == sts[ry][rx-1]
									 ||
									 1 == sts[ry][rx+1] && 1 == sts[ry+1][rx]
									 ) ) 
	     ||
	     // horizontal and vertical (edge/around) lines
	     ( (1 == sts[ry-1][rx-1] && 1 == sts[ry-1][rx+1]) && !(1 == sts[ry-1][rx]) ) // and the circle (3 cells around)....!!! TODOOOOO MISSING!!!
	     ||
	     ( (1 == sts[ry-1][rx-1] && 1 == sts[ry+1][rx-1]) && !(1 == sts[ry][rx-1]) ) // and the circle (3 cells around)....!!!
	     ||
	     ( (1 == sts[ry-1][rx+1] && 1 == sts[ry+1][rx+1]) && !(1 == sts[ry][rx+1]) ) // and the circle (3 cells around)....!!!
	     ||
	     ( (1 == sts[ry+1][rx+1] && 1 == sts[ry+1][rx-1]) && !(1 == sts[ry+1][rx]) ) // and the circle (3 cells around)....!!! 
	     // TODOOOOO: MISSING? edges around, with distance 2/3 between '1' cells

	       );
  }

  return risk;
}


// TODO: connectivity_check4()

// TODO: connectivity_check8()

#include <boost/scoped_array.hpp>

//#include "UnionFind.hpp"
//#include "DisjointSets.h"
//#include <boost/pending/disjoint_sets.hpp>

inline int Grid_CCSP::get_get_get_cc_label(int x, int y)
{
  if (thickness <= 1)
    return get_cc_label(x, y);
  else 
    return uneroded_lbl_buffer[y][x];
}

void
Grid_CCSP::get_split_sizes_n_lens_n_wrscr(int x, int y, std::vector<size_t>& sizes, 
					  std::vector<size_t>& lens, std::vector<float>& wrscr)
{
  int inc = thickness_inc + 1;  // when thickness==1, thickness_inc==0, inc==1
  // TODO: is this correct, or does it need to check the diamond around?

  std::vector<int> labels;
  // 4 basic ones for 4-connected version
  if (status[y][x-inc] > 0 && get_get_get_cc_label(x-inc, y) >= 0)
    labels.push_back(get_get_get_cc_label(x-inc, y));
  if (status[y-inc][x] > 0 && get_get_get_cc_label(x, y-inc) >= 0)
    labels.push_back(get_get_get_cc_label(x, y-inc));
  if (status[y][x+inc] > 0 && get_get_get_cc_label(x+inc, y) >= 0)
    labels.push_back(get_get_get_cc_label(x+inc, y));
  if (status[y+inc][x] > 0 && get_get_get_cc_label(x, y+inc) >= 0)
    labels.push_back(get_get_get_cc_label(x, y+inc));
  
  if (use_8_conn) {
    // additional ones for 8-connected version
    if (status[y-inc][x-inc] >= 0 && get_get_get_cc_label(x-inc, y-inc) >= 0) 
      labels.push_back(get_get_get_cc_label(x-inc, y-inc));
    if (status[y-inc][x+inc] >= 0 && get_get_get_cc_label(x+inc, y-inc) >= 0)
      labels.push_back(get_get_get_cc_label(x+inc, y-inc));
    if (status[y+inc][x+inc] >= 0 && get_get_get_cc_label(x+inc, y+inc) >= 0)
      labels.push_back(get_get_get_cc_label(x+inc, y+inc));
    if (status[y+inc][x-inc] >= 0 && get_get_get_cc_label(x-inc, y+inc) >= 0)
      labels.push_back(get_get_get_cc_label(x-inc, y+inc));
  }
  
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

  bool debug_lbls_sizes = false;
  if (CCSP_verbose > 2) {
    debug_lbls_sizes = true;
  }
  if (debug_lbls_sizes) {
    std::cerr << "get_split_sizes_n_lens_n_wrscr(), x: " << x << ", y: " << y << ", sizes-labels: ";
    for(i = cc_labels.begin(); i != cc_labels.end(); i++) {
      std::cerr << "lbl: " << *i << ", a: " << this->areas[*i] << ", r: " << this->wrscr[*i] << "   ";
    }
    std::cerr << std::endl;
  }
}

// expects status_bkp (eroded) in sts
// If there are no more changes -> TODO: rename to save_state_erode_foreground
void
Grid_CCSP::update_erode_foreground(int x, int y, size_t thickness, char**& sts, char**& edge_bkp_eroded)
{
  // erode a diamond around x,y
  int& inc = thickness_inc; // round down
  int xmin = max(0, x-inc);
  int xmax = min((int)xdim-1, x+inc);
  int ymin = max(0, y-inc);
  int ymax = min((int)ydim-1, y+inc);
  // save state
  for(int j=ymin; j<=ymax; j++) {
    for(int i=xmin; i<=xmax; i++) {
      erosion_update_state_sts[j-ymin][i-xmin] = sts[j][i];
    }
  }

  // for eroded edge within removed
  int inc_edge = thickness_inc+1; // beware of +1
  int xmin_edge = max(0, x-inc_edge);
  int xmax_edge = min((int)xdim-1, x+inc_edge);
  int ymin_edge = max(0, y-inc_edge);
  int ymax_edge = min((int)ydim-1, y+inc_edge);
  for(int j=ymin_edge; j<=ymax_edge; j++) {
    for(int i=xmin_edge; i<=xmax_edge; i++) {
      erosion_update_state_edge[j-ymin_edge][i-xmin_edge] = edge_bkp_eroded[j][i];
    }
  }

  // DONT ERODE MORE!!!!!!!!!!!!!!!!! LEAVE IT FOR THE TENTATIVE ERODE IN CHECK_SPLIT_LOCAL_BFS+DIAMOND
  return;

  // Star!
  // (x-inc)
  sts[y][xmin] = 0;
  edge_bkp_eroded[y][xmin] = 0;
  edgify_neighbors(xmin, y, sts, edge_bkp_eroded);
  // (y-inc)
  sts[ymin][x] = 0;
  edge_bkp_eroded[ymin][x] = 0;
  edgify_neighbors(x, ymin, sts, edge_bkp_eroded);
  // (x+inc)
  sts[y][xmax] = 0;
  edge_bkp_eroded[y][xmax] = 0;
  edgify_neighbors(xmax, y, sts, edge_bkp_eroded);
  // (y+inc)
  sts[ymax][x] = 0;
  edge_bkp_eroded[ymax][x] = 0;
  edgify_neighbors(x, ymax, sts, edge_bkp_eroded);
  return;


  /* ************************************************************************************************************************************ */
  // Square rather than diamond!
  for(int j=ymin; j<=ymax; j++) {
    for(int i=xmin; i<=xmax; i++) {
      if (sts[j][i] >= 0) {
	sts[j][i] = 0;
	// remove remainders (prev. erosion step)
	edge_bkp_eroded[j][i] = 0;
	// and make new edge cells
	edgify_neighbors(i, j, sts, edge_bkp_eroded);
      }
    }
  }
  return;




  // [W->S]
  for(int j=y; j<=ymax; j++) {
    for(int i=xmin; i<=x; i++) {
      if (sts[j][i] >= 0) {
	status_bkp[j][i] = 0;
	// remove remainders (prev. erosion step)
	edge_bkp_eroded[j][i] = 0;
	// and make new edge cells
	edgify_neighbors(i, j, sts, edge_bkp_eroded);
      }
    }
  }
  // (S->E]
  for(int j=ymax-1; j>=y; j--) {
    for(int i=x+1; i<=xmax; i++) {
      if (sts[j][i] >= 0) {
	status_bkp[j][i] = 0;
	edge_bkp_eroded[j][i] = 0;
	edgify_neighbors(i, j, sts, edge_bkp_eroded);
      }
    }
  }
  // (E->N]
  for(int j=y-1; j>=ymin; j--) {
    for(int i=xmax-1; i>=x; i--) {
      if (sts[j][i] >= 0) {
	status_bkp[j][i] = 0;
	edge_bkp_eroded[j][i] = 0;
	edgify_neighbors(i, j, sts, edge_bkp_eroded);
      }
    }
  }
  // (N->W)
  for(int j=ymin+1; j<=y-1; j++) {
    for(int i=x-1; i>=xmin+1; i--) {
      if (sts[j][i] >= 0) {
	status_bkp[j][i] = 0;
	edge_bkp_eroded[j][i] = 0;
	edgify_neighbors(i, j, sts, edge_bkp_eroded);	
      }
    }
  }

}

void 
Grid_CCSP::edgify_neighbors(int x, int y, char**& sts, char**& edge_bkp_eroded)
{
  // all neighbors to edge
  size_t xl = max(0, x-1);
  size_t xr = min(x+1, (int)xdim-1);
  size_t yl = max(0, y-1);
  size_t yr = min(y+1, (int)ydim-1);
  for(size_t i=xl; i<=xr; i++) {
    for(size_t j=yl; j<=yr; j++) {
      // 2nd term is critical: do not include what is normal/uneroded edge!
      if((j==yl || j==yr || i==xl || i==xr) &&
	 (sts[j][i] > 0) /*&& (1!=edge[j][i])*/ /*&& (1!=edge_bkp_eroded[j][i])*/)
	edge_bkp_eroded[j][i] = 1;
    }
  }
}

void
Grid_CCSP::undo_update_erode_foreground(int x, int y, size_t thickness, char**& sts, char**& edge_eroded)
{
  int& inc = thickness_inc; // round down
  int xmin = max(0, x-inc);
  int xmax = min((int)xdim-1, x+inc);
  int ymin = max(0, y-inc);
  int ymax = min((int)ydim-1, y+inc);

  // restore state before update_erode_foreground(x,y,thickness)
  for(size_t j=ymin; j<=ymax; j++) {
    for(size_t i=xmin; i<=xmax; i++) {
      sts[j][i] = erosion_update_state_sts[j-ymin][i-xmin];
      // edge_bkp_eroded[j][i] = erosion_update_state_edge[j-ymin][i-xmin]; -> below
    }
  }

  // for eroded edge within removed
  int inc_edge = thickness_inc+1; //  ARRRRRRRRRGGGG? - NOP???
  int xmin_edge = max(0, x-inc_edge);
  int xmax_edge = min((int)xdim-1, x+inc_edge);
  int ymin_edge = max(0, y-inc_edge);
  int ymax_edge = min((int)ydim-1, y+inc_edge);
  for(int j=ymin_edge; j<=ymax_edge; j++) {
    for(int i=xmin_edge; i<=xmax_edge; i++) {
      edge_bkp_eroded[j][i] = erosion_update_state_edge[j-ymin_edge][i-xmin_edge];
    }
  }
}

void
Grid_CCSP::do_erode_foreground(char**& sts_orig, char**& sts, size_t thickness, char**& edge_bkp)
{
  // this init wouldn't be needed if do_erode were run only once at the beginning => which is not the case -> TODO: remove
  for(size_t y=1; y < ydim-1; y++) {	    
    for(size_t x=1; x < xdim-1; x++) {
      erosion_mask[y][x] = 0;
    }
  }

  //for (int r=0; r<int(thickness/2); r++) {
  for (int r=0; r<thickness_inc; r++) {
    /*
    // erosion based on a 3x3 kernel/structuring element
    for (int r=0; r<thickness_inc; r++) {
      for(size_t y=1; y < ydim-1; y++) {	    
	for(size_t x=1; x < xdim-1; x++) {

	  erosion_mask[u][x] = 0;

	  if (0 >= sts[y][x])
	    continue;
	  
	  if (sts[y-1][x-1]<=0 || sts[y-1][x]<=0 || sts[y-1][x+1]<=0
	      ||
	      sts[y][x-1]<=0  ||  sts[y][x+1]<=0
	      ||
	      sts[y+1][x-1]<=0 || sts[y+1][x]<=0 || sts[y+1][x+1]<=0)
	    erosion_mask[y][x] = 1;
	}
      }
    }
*/
    // This would be the "negative" alternative
    // erode around 0/removed/background cells
    // CANNOT FORGET THE 1st/second row / column, and the second last/last!
    // the first and last are already 0 from the beginning, now go to the second and second last...
    if (0 == r) {
      for(size_t y=1; y < ydim-2; y++) {
	erosion_mask[y][1] = 1;
	erosion_mask[y][xdim-2] = 1;
      }
      for(size_t x=1; x < xdim-2; x++) {
	erosion_mask[1][x] = 1;
	erosion_mask[ydim-2][x] = 1;
      }
    }
    for(size_t y=1; y < ydim-1; y++) {	    
      for(size_t x=1; x < xdim-1; x++) {
	
	edge_bkp[y][x] = 0; // init for later

	if (1 == sts[y][x])
	  continue;

	// -1/0 == sts[y][x])

	// W-E-N-S
	if (1 == sts[y][x-1])
	  erosion_mask[y][x-1] = 1;
	if (1 == sts[y][x+1])
	  erosion_mask[y][x+1] = 1;
	if (1 == sts[y-1][x])
	  erosion_mask[y-1][x] = 1;
	if (1 == sts[y+1][x])
	  erosion_mask[y+1][x] = 1;

	// corners
	if (1 == sts[y+1][x+1])
	  erosion_mask[y+1][x+1] = 1;
	if (1 == sts[y-1][x-1])
	  erosion_mask[y-1][x-1] = 1;
	if (1 == sts[y+1][x-1])
	  erosion_mask[y+1][x-1] = 1;
	if (1 == sts[y-1][x+1])
	  erosion_mask[y-1][x+1] = 1;
      }
    }


    for(size_t y = 0; y < ydim; y++) {	    
      for(size_t x = 0; x < xdim; x++) {
	if (sts_orig[y][x]>0 && 1!=erosion_mask[y][x])
	  //if (/*1 == sts[y][x] && */ 1 == erosion_mask[y][x])
	  sts[y][x] = 1;
	else
	  sts[y][x] = 0;  // TODOOOOOOOOOOOOOOO: -1 or 0? => 0

      }
    }
  }

  /*
  // this init wouldn't be needed if do_erode were run only once at the beginning
  for(size_t y=1; y < ydim-1; y++) {	    
    for(size_t x=1; x < xdim-1; x++) {
      edge_bkp_eroded[y][x] = 0;
    }
    }*/


  // Generate edge for eroded map
  // the cells that have been eroded above are now 'removed'
  for(int y=0; y < ydim; y++) {	    
    for(int x=0; x < xdim; x++) {
      //if(-1==sts[y][x]) {  <- this is the condition for edge[][]
      if(1==erosion_mask[y][x]) {
	// all neighbors to edge
	int xl = max(0, x-1);
	int xr = min(x+1, (int)xdim-1);
	int yl = max(0, y-1);
	int yr = min(y+1, (int)ydim-1);
	for(size_t i=xl; i<=xr; i++) {
	  for(size_t j=yl; j<=yr; j++) {
	    // sts[j][i] > 0 (==1) is guaranteed by the loop above (in this init method)
	    if((sts[j][i] > 0) /*&& (1!=edge_bkp_eroded[j][i])*/  /* && (1!=edge[j][i]))   // never put an eroded-edge cell on top of an edge cell*/
	       )
	      edge_bkp[j][i] = 1;
	  }
	}
      }
    }
  }

}

bool
Grid_CCSP::do_save_corridor_boundaries_layer()
{
  // 20% const int remaining_thresh = 22056;
  // 12% const int remaining_thresh = 97048; remaining: 13234
  // 25%, remaining: 27570
  //int remaining_thresh = 27570;
  //int rem_level = nonm1 - remaining_thresh;
  //if (save_cached_layer && (remaining_thresh >= nonm1-removed)) {
  if (save_cached_layer && corr_set.out_boundaries_pcs.size()>0 && 
      (corr_set.out_boundaries_pcs[0]/100.0f*float(nonm1) >= nonm1-removed)) {
    //save_cached_layer = false;
      // save raster
    String name_fract_appendix = FloatToStr(corr_set.out_boundaries_pcs[0]);
    int remaining_thresh = corr_set.out_boundaries_pcs[0]/100.0f*float(nonm1);
    int rem_level = nonm1 - remaining_thresh;
    //corr_set.out_boundaries_pcs.pop_back();
    corr_set.out_boundaries_pcs.erase(corr_set.out_boundaries_pcs.begin());
    float nodatavalue(-1.0f);
    boost::scoped_array<float> plane(new float[xdim * ydim]);
    for(int y = 0; y < ydim; ++y) {
      for(int x = 0; x < xdim; ++x) {
	int i(x + y * xdim);
	//	  if(cached[y][x] <= -1) {
	//if (cached[y][x] >= 0)  // 
	//if (cached[y][x] >= rem_level-110)  // 110 bigger than warp factor
	if (cached[y][x] >= rem_level-warp_ll)  // 
	  plane[i] = float(cached[y][x])/float(nonm1);
	else
	  plane[i] = nodatavalue;
	/*
	  if(cached[y][x] <= 0) {
	  plane[i] = nodatavalue;
	  } else {	
	  if (cached[y][x] > 0) {
	  plane[i] = 1;//cached[y][x];
	  } else {
	  plane[i] = nodatavalue;
	  }
	  }*/
      }
    }

    String output_corr_suffix = "_corridor_boundaries";
    String projname = emffn;
    String projdir = projname.left(projname.lastIndexOf('/')+1);  // path of output (with trailing '/' !)
    String out_base = emffn.mid(emffn.lastIndexOf('/')+1);
    String output_corr_subdir = ChangeFileExt(out_base) + output_corr_suffix;
    String subdirname = createSubDirIfNeeded(projdir, output_corr_subdir); // full path

    String oname = subdirname + "/" + "corridor_boundaries_top_"+name_fract_appendix;
    SaveToRaster<float>(oname, plane.get(), nodatavalue, xdim, ydim);
    // SaveToRaster<float>(oname"z_cached_layer_remaining_"+IntToStr(remaining_thresh), plane.get(), nodatavalue, xdim, ydim);

    return true;
  } else
    return false;
}

float
Grid_CCSP::calc_domain_coefficient(int rx, int ry)
{
  float domain_coeff = 1.0f;
  /* OLD Single domain layer code:
  if (corr_set.use_domain_layers) {
    // Apply domain layers:
    if (corr_set.layers.size() > 0 && (corr_set.layers[0].m[ry][rx] <= 0 || isnan(corr_set.layers[0].m[ry][rx])))
      return 0;

    domain_coeff = corr_set.weights[0] * corr_set.layers[0].m[ry][rx];
  }
  */

  // Process the list of domain layers
  if (corr_set.use_domain_layers) {
    // Apply domain layers:
    domain_coeff = .0f;
    for(size_t lidx=0; lidx<corr_set.layers.size(); lidx++) {
      if (corr_set.layers[0].m[ry][rx] >= .0f /* implied: && isnan(corr_set.layers[0].m[ry][rx]))*/ ) {
	domain_coeff += corr_set.weights[lidx] * corr_set.layers[lidx].m[ry][rx];
      }
    }
    if (.0f == domain_coeff)
      return .0f;
  }

  return domain_coeff;
}
