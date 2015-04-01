#ifndef AM_GRIDMAP_H
#define AM_GRIDMAP_H
#include <VCL.h>

class GridMap
{
private:
	float  xc, yc; // corners of data matrix
	float  xl, yl;
	float  cs;     // cell size
	float rr, rp;
	int   xdim, ydim;   // dimensions of landscape
	int   ixdim, iydim; // dimensions of aligned image
	int   ok_cells;     // number of cells with data
	float   nodataval;

	//      int    *r, *g, *b;
	//      bool   *is_suitable;
	//      TColor *mcols;
	bool  flip;
	bool  loaded;
	int   matrix_passed;

public:
	size_t non_missing_cnt;
	double sum, sum_orig;  // sum_orig: before normalization

	double dxc, dyc;
	double dcs;
	//double sum;
	bool  normalize;
	bool  miss_as_zero;
	GridMap();
	~GridMap();
	float **m;
	float max_val;

	void dealloc();
	void pixelblock(TCanvas *c, int x, int y, TColor col, bool missing);
	bool load_from_file(const String& fname, int mask_data, float **area_mask, bool setGlobalProjection = false);
	bool load_from_file_IG(const String& fname,
			       struct IG_settings& IG_set, float sp_w, float **IGw,
			       int logit, float &psum, float &IG_sum,  int mask_data,
			       float **area_mask, float **tmpm, float** occur_w_map=NULL);

	// counts non-missing occurrencies into occur_count_matrix. 
	bool load_from_file_IG_just_to_count(const String& fname, int mask_data, float **area_mask, 
					     float **tmpm, int** occur_count_matrix);

	int  rows(){return ydim;}
	int  cols(){return xdim;}
	float  nodata(){return nodataval;}
	float  getxc(){return xc;}
	float  getyc(){return yc;}
	float  matval(int x, int y) {return m[y][x];}
	float  cell_size(){return cs;}
	void set_flip(bool fl){flip=fl;}
	float get_pixel_size(){return rr;}
	void get_dimensions(float &xcorner, float &ycorner,
			    float &cell_size, int &xd, int &yd);
	bool is_loaded(){return loaded;}

	float calc_r_etal(int ixd, int iyd);
	bool align_with_image(int xd, int yd, int &new_xd, int &new_yd);
#if 0
	// commented out to et rid of use of xc, yc
	void get_im_coordinates(float x, float y, int &new_xd, int &new_yd);
	//      bool set_color_for_type(int val, int R, int G, int B, bool suitable); not good for floats

	float get_value_at_point(float x, float y);
	float get_value_at_matrix_point(int x, int y, float &cx, float &cy);
	int   get_cell_center_etc_at_matrix(int X, int Y, float &cx, float &cy);
	void get_matrix_coord_at_point(float x, float y, int &lx, int &ly);
	TColor get_color_at_point(float x, float y);

	float get_value_at_image_point(int x, int y);
	TColor get_color_at_image_point(int x, int y);
	TColor get_color_for_val(float val);
	bool draw_on_bitmap(Graphics::TBitmap *bm);
#endif

//#if 0
	void Annotate(Graphics::TBitmap *bm);
	bool draw_on_bitmap_smooth(Graphics::TBitmap *bm, float maxv,
				   float cut_level, bool rank_mode);
	bool Is_on_border(float **mm, int y, int x);
	bool draw_on_bitmap_smooth_from_mat(Graphics::TBitmap *bm, float **m,
					    float maxv, float cut_level, bool rank_mode,
					    int green_scale=0);
	void set_color_for_pixel(Graphics::TBitmap *bm,
				 int x, int y, float maxv, float val, bool rank_mode);
	//      bool add_to_bitmap(Graphics::TBitmap *bm);     // xc issue
	//      float get_suitable_area(float x, float y, float a);  // xc issue
	bool show_spots(Graphics::TBitmap *bm, int **cm, int comp_mode);
//#endif

	void export_GIS_INT(int **nwm, const String& fname);
	void set_no_normalize(){normalize=false;}
	void average_normalize();
	// does not do a deep-free if the matrix had been passed from a previously allocated buffer
	void free_matrix_m();
	// frees buffer m regardless of anything else. Needed when mapgrid objects are loaded more than once
	void force_free_matrix_m();
	bool loadFromFile(const char *fileName, float **tmpm = 0, bool setGlobalProjection = false);
};

#endif
