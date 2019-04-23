#include "conn_components.h"

#include <algorithm>
#include <iostream>
#include <boost/scoped_array.hpp>

#include "matrix_utils.h"
// only for saveToRaster!  <- TODO: remove
#include "Unit1.h"

//size_t Grid_Connected_Components::n_cc_lbl_calc=0;

Grid_Connected_Components::Grid_Connected_Components(size_t xdim, size_t ydim):
  xdim(xdim), ydim(ydim), n_cc_lbl_calc(0)
{
}

Grid_Connected_Components::~Grid_Connected_Components()
{
}

int
Grid_Connected_Components::count_labels(char**& sts, int**& lbl_layer)
{
  if (!sts || !lbl_layer)
    return -99;

  n_cc_lbl_calc++;
  //int variant = 1;  // Chang-Chen-Lu
  int ccl_variant = 2;  // Chang-Chen-Lu
  //int ccl_variant = 2;  // 2 Pass with union-find
  if (1 == ccl_variant) {
    CC_Labeling_Chang_Chen_Lu ccl(xdim, ydim);
    int v = ccl.do_cc_labeling_Chang_Chen_Lu(sts, lbl_layer, false);
    std::cerr << "LABELS CHANG: " << v;
    return v;
  } else {
    CC_Labeling ccl(16000, xdim, ydim);
    return ccl.do_cc_labeling_2pass(sts, lbl_layer);
  }

}

int
CC_Labeling::do_cc_labeling_2pass(char**& sts, int**& lbl_layer)
{
  static size_t n_calc = 0;
  // A) first pass
  // label 0 is background
  last_lbl_equiv = 0;

  //RZ::UnionFind<int> uf;
  //DisjointSets dj;
  //std::vector<int> rank(1);
  //std::vector<int> parent(1);
  //boost::disjoint_sets<int*,int*> bdj(&rank[0], &parent[0]);
  for(size_t i = 0+1; i < ydim; i++) {	    
    for(size_t j = 0+1; j < xdim; j++) {
      // Missing/NaN cell
      // For "edge removal = 0"
      if (sts[i][j] <= 0)
	continue;

      // For "edge removal = 1"
      //if (sts[i][j] <= 0)   /* || 1 != edge[i][j] ===> requires changes in the labeling (state: last label used on the left, and discontinuities)*/
	//continue;

      // now, 0 == sts[i][j] are background cells
      //    , 1 == sts[i][j] are ON cells

      // TODO: use a decision tree / Wu-Otoo-Suzuki
    
      // Same as left (W), there is no usable (N)
      if (sts[i][j-1] == 1 
	  &&
	  sts[i-1][j] != 1
	  ) {
	lbl_layer[i][j] = lbl_layer[i][j-1];
	continue;
      }

      // Merge above (N) and left (W) 
      if (sts[i-1][j] == 1 && sts[i][j-1] == 1) {
	if ( lbl_layer[i-1][j] == lbl_layer[i][j-1] ) {
	  lbl_layer[i][j] = lbl_layer[i][j-1];	  
	} else {
	  int min_lbl = std::min(lbl_layer[i-1][j], lbl_layer[i][j-1]);
	  //int max_lbl = std::max(lbl_layer[i-1][j], lbl_layer[i][j-1]);
	  lbl_layer[i][j] = min_lbl;
	  if (min_lbl >= 0)
	    //uf.setEqual(min_lbl, max_lbl);
	    //dj.Union(min_lbl-1, max_lbl-1);
	    //bdj.link(min_lbl, max_lbl);
	    this->merge(lbl_layer[i-1][j], lbl_layer[i][j-1]);
	  //equivalences.insert(std::pair<int,int>(min_lbl, max_lbl));
	  //equivalences.insert(std::pair<int,int>(min_lbl, max_lbl));
	  //std::cerr << "EQUIV! ";
	  // RECORD EQUIVALENCE
	  last_lbl_equiv++;
	}
	continue;
      }

      // Take same label as above (N)
      //if (sts[i][j-1] != 1 /*sts[i][j]*/ || 
      if (sts[i-1][j] == 1 && sts[i][j-1] != 1) {
	lbl_layer[i][j] = lbl_layer[i-1][j];
	//std::cerr << "SAA ";
	continue;
      }
      
      if (sts[i][j-1] != 1 /*sts[i][j]*/ &&
	  sts[i-1][j] != 1 /*sts[i][j]*/) {
	//lbl_layer[i][j] = ++lbln;
	//uf.insert(lbln);
	//dj.AddElements(1);
	//bdj.make_set(lbln);
	lbl_layer[i][j] = this->new_label();
      }
	
    }
  }

  bool save_layers = false;//true;
  if (save_layers && n_calc < 1) {
    // save raster
    float nodatavalue(-1.0f);
    boost::scoped_array<float> plane(new float[xdim * ydim]);
    for(int y = 0; y < ydim; ++y) {
      for(int x = 0; x < xdim; ++x) {
	int i(x + y * xdim);
	if(sts[y][x] == -1) {
	  plane[i] = nodatavalue;
	} else {	
	  plane[i] = lbl_layer[y][x];
	}
      }
    }
    SaveToRaster<float>("lbl_layer.1st", plane.get(), nodatavalue, xdim, ydim);
  }


  //  std::cerr << "size of ccl.labels set before merge: " << ccl.labels.size(); << std::endl;

  // A) second pass
  //int key_total = uf.size();
  //int key_total = dj.NumElements();
  //int key_total = 0; //bdj.count_sets(0, std::numeric_limits<int>::max());
  //int key_total = this->.labels.size();
  // should be compress_labels()
  int newtag = 0;
  for(int id=0; id<this->labels.size(); ++id)
    if(this->is_root_label(id))
      this->labels[id].tag = newtag++;

  for(size_t i = 0+1; i < ydim; i++) {	    
    for(size_t j = 0+1; j < xdim; j++) {
      if (sts[i][j] <= 0)   /* || 1 != edge[i][j] */
	continue;

      //      std::cerr << " second pass: " << lbl_layer[i][j] << std::endl;  
      //int me = uf.findd(lbl_layer[i][j]);
      //int me = dj.FindSet(lbl_layer[i][j] - 1) + 1;
      //int me = bdj.find_set(lbl_layer[i][j]);
      lbl_layer[i][j] = this->labels[this->root_of(lbl_layer[i][j])].tag;

      //if (lbl_layer[i][j] > lbln ) {
      /*
      if (lbl_layer[i][j] > this->labels.size() ) {
	std::cerr << "E2: i: " << i << ", j: " << j << ", label: " << lbl_layer[i][j] << ", newtag: " << newtag << ", ccl.labels.size(): " << ccl.labels.size();
	lbl_layer[i][j] = 0;
      }
      */
    }
  }  

  if (save_layers && n_calc < 1) {
    // save raster
    float nodatavalue(-1.0f);
    boost::scoped_array<float> plane2(new float[xdim * ydim]);
    for(int y = 0; y < ydim; ++y) {
      for(int x = 0; x < xdim; ++x) {
	int i(x + y * xdim);
	if(sts[y][x] == -1) {
	plane2[i] = nodatavalue;
	} else {	
	  plane2[i] = lbl_layer[y][x];
	}
      }
    }
    SaveToRaster<float>("lbl_layer.2nd", plane2.get(), nodatavalue, xdim, ydim);
  }

  n_calc++;
  // corrected by subtracting redundant labels
  return newtag;
}

// Chang-Chen-Lu 1-pass method
/*
static int search_dir[8][2] = { {0,1},{1,1},{1,0},{1,-1},
				{0,-1},{-1,-1},{-1,0},{-1,1}
};
*/
static int search_dir[4][2] = { {0,1},{1,0},
				{0,-1},{-1,0}
};
// unused for now
static int **contourmap;


CC_Labeling_Chang_Chen_Lu::CC_Labeling_Chang_Chen_Lu(int xdim, int ydim): 
  xdim(xdim), ydim(ydim)
{
}

inline void
CC_Labeling_Chang_Chen_Lu::tracer(/*char**& contour_sts, int**& lbl_layer, */ int& cy, int& cx, int& tracingdirection)
{
  //for(int i = 0; i < 7; i++) { // 8-conn
  for(int i = 0; i < 4; i++) {
    int y = cy + search_dir[tracingdirection][0];
    int x = cx + search_dir[tracingdirection][1];

    if(contour_sts[y][x] <= 0) {
      contour_lbl_layer[y][x] = -1;
      //tracingdirection = (tracingdirection + 1) % 8; // 8-conn
      tracingdirection = (tracingdirection + 1) % 4;
    } else {
      cy = y;
      cx = x;
      break;
    }
  }
}

inline void
CC_Labeling_Chang_Chen_Lu::contour_tracing(/*char**& sts, int**& lbl_layer, */ int cy, int cx, int label_idx, 
					    int tracingdirection, bool gen_contour)
{
  char tracingstopflag = 0, SearchAgain = 1;
  int fx, fy, sx = cx, sy = cy;

  tracer(/*sts, lbl_layer, */ cy, cx, tracingdirection);
  if (gen_contour)
    contourmap[cy][cx] = label_idx;

  if (cx != sx || cy != sy) {
    fx = cx;
    fy = cy;	  
    while (SearchAgain) {
      //tracingdirection = (tracingdirection + 6) % 8; // 8-conn
      tracingdirection = (tracingdirection + 3) % 4;
      contour_lbl_layer[cy][cx] = label_idx;
      tracer(/* sts, lbl_layer, */ cy, cx, tracingdirection);
      if (gen_contour)
	contourmap[cy][cx] = label_idx;
      
      if (cx == sx && cy == sy) {
	tracingstopflag = 1;
      }
      else if(tracingstopflag) {
	if(cx == fx && cy == fy) {
	  SearchAgain = 0;
	} else {
	  tracingstopflag = 0;
	}
      }
    }
  }
}

// Labels are nonnegative integers, <0 for undefined cells
int
CC_Labeling_Chang_Chen_Lu::do_cc_labeling_Chang_Chen_Lu(char**& sts, int**& lbl_layer, bool gen_contour)
{
  size_t cy, cx;
  int ccl_num = 0, label_idx = -1, tracingdirection = 0;

  contour_sts = sts;
  contour_lbl_layer = lbl_layer;

  memset(&(contour_lbl_layer[0][0]), 0, xdim*ydim* sizeof(int));
  /*
  for(size_t i = 0; i < ydim; i++) {	    
    for(size_t j = 0; j < xdim; j++) {
      contour_lbl_layer[i][j] = 0;
    }
  }
  */
  for(cy = 1; cy < ydim - 1; cy++){
    for(cx = 1, label_idx = -99; cx < xdim - 1; cx++) {
      //std::cerr << cx << ","<< cy << " ";      
      if(sts[cy][cx] > 0) { // ON cell
	if(label_idx != 0) { // use pre-cell label
	  //std::cerr << "label_idx!=0: " << label_idx << " ";
	  lbl_layer[cy][cx] = label_idx;
	} else {
	  label_idx = lbl_layer[cy][cx];
	  //std::cerr << "label_idx==0: " << label_idx << " ";
	  if(label_idx == 0) {
	    // std::cerr << "label_idx==0: " << label_idx << " ";
	    label_idx = ++ccl_num;
	    tracingdirection = 0;
	    contour_tracing(/* sts, lbl_layer, */ cy, cx, label_idx, tracingdirection, false);// external contour
	    lbl_layer[cy][cx] = label_idx;
	  }
	}
      } else if(label_idx != 0) { // OFF cell & pre-cell has been labeled
	if(lbl_layer[cy][cx] == 0) {
	  //std::cerr << "lbl_layer[cy][cx]==0: " << label_idx << " ";
	  tracingdirection = 1;
	  contour_tracing(/* sts, lbl_layer, */ cy, cx - 1, label_idx, tracingdirection, false);// internal contour
	}
	label_idx = 0;
      }
    }
  }
  
  return ccl_num;
}
