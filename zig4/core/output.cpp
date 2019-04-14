#include "GridMap.h"
#include "defines.h"
#include "output.h"
#include "bat_run.h"
#include "Unit1.h"
#include "matrix_utils.h"
#include <boost/scoped_array.hpp>
#include <cstdio>
#include <cmath>

float get_cost_fract(float cval)
{
	float c1, c2, f1, f2, f, diff;
	int   y, ok;

	ok=0;
	for(y=0; y<(c_pos-1); y++)
	{
		c1=curves[y][1];   // c2<c1
		c2=curves[y+1][1];
		f1=curves[y][0];
		f2=curves[y+1][0]; // f2>f1
		if ((c1>=cval) && (c2<=cval))
		{
			ok=1;
			break;
		}
	}
	if (!ok)
		return -1.0f;

	f = f2 - (cval-c2)*(f2-f1)/(c1-c2);

	return (1.0f-f);
}

std::vector<float> distrib_centers_x;
std::vector<float> distrib_centers_y;

// This iterates through y,x,features -> takes time
// to speed it up it could be done together with richness calculations in Calc_richness_et_al_matrixes()
void
init_output_get_distr_centers_and_sums(bool do_sums)
{
  if (do_sums)
    Form1->Memo1->Lines->Add("Calculating centers and sums of the distributions of "+IntToStr(map_cnt)+" features...");
  else
    Form1->Memo1->Lines->Add("Calculating centers of the distributions of "+IntToStr(map_cnt)+" features...");

  distrib_centers_x.resize(map_cnt);
  distrib_centers_x.assign(distrib_centers_x.size(),0);
  distrib_centers_y.resize(map_cnt);
  distrib_centers_y.assign(distrib_centers_y.size(),0);
  // for(int s=0; s<map_cnt; s++) {
  //   get_distr_center(s, distrib_centers_x[s], distrib_centers_y[s]);
  // }
  float wxsum = 1.0f;
  float wysum = 1.0f;
  std::vector<float> sum(map_cnt, 1.0f);
  for (int y=0; y<yd; y++) {
    for (int x=0; x<xd; x++) {
      const Biodiv_Features_Occur_Container& rowp = vmat[y][x];
      for (int s = rowp.first(); s != rowp.overflow(); s = rowp.next(s)) {
	float val = rowp[s]; // == vmat[y][x][s];
	if (val != -1) {
	  distrib_centers_x[s] += val*x;
	  distrib_centers_y[s] += val*y;
	  sum[s] += val;
	}
      }
    }
  }
  for (int s=0; s<map_cnt; s++) {
    if (0 != sum[s]) {
      distrib_centers_x[s] /= sum[s];
      distrib_centers_y[s] /= sum[s];
    }

    if (do_sums) {
      spp[s].prob_sum = sum[s];
      spp[s].IG_fract = sum[s];
    }
  }
  
#ifdef ZCORE_DEBUG
  if ( isnan(dx) || isnan(dy) )
    Form1->Memo1->Lines->Add("NAN in get_distr_center!!!");
#endif
}

void
output_features_info_n_curves()
{
  FILE *finfo_of;
  int x,y, pos_w_cnt;
  double ext1, ext2, wsum, Tviol, dx, dy;

  finfo_of = fopen(features_info_fn.toUtf8().constData(), "w+t");	
  if (!finfo_of) {
    Form1->Memo1->Lines->Add("************** ERROR ***************");
    Form1->Memo1->Lines->Add("  Error opening feature info data file "+features_info_fn);
    return;
  }

  fprintf(finfo_of, "List of biodiversity features (e.g. species) and weights used in analysis in order of columns:\n");
  fprintf(finfo_of, "Weight  distribution-sum  IGRetained TviolationFractRem  Distr-mean-X   Distr-mean-Y  MapFileName\n");
  Tviol=100;
  for(x=0; x<map_cnt; x++) {
    // get_distr_center(x, dx, dy); // removed from here to avoid dependence on vmat -> now its init at the loading stage
    fprintf(finfo_of, "%-6.3f\t%-6.3f\t%0.2f\t%0.5f\t%-.2f\t%-.2f\t%s\n", spp[x].weight, spp[x].prob_sum,
	    spp[x].IG_fract, 1.0f-fabs(spp[x].T_violation_fract), distrib_centers_x[x], distrib_centers_y[x], spp[x].fname.toUtf8().constData());
    if (spp[x].T_violation_fract<Tviol)
      Tviol=spp[x].T_violation_fract;
  }

  if (removal_rule==3) {
    fprintf(finfo_of, "First biodiversity feature target violated at fraction %f of landscape removed, %f remaining\n", Tviol, 1.0f-Tviol);
    Form1->Memo1->Lines->Add("First biodiversity feature target violated at fraction of landscape remaining = "+FloatToStr(1.0f-Tviol));
  }

  fclose(finfo_of);


  FILE *curves_of;
  curves_of = fopen(curvesfn.toUtf8().constData(), "w+t");	
  if (!curves_of) {
    Form1->Memo1->Lines->Add("************** ERROR ***************");
    Form1->Memo1->Lines->Add("  Error opening curves data file "+curvesfn);
    return;
  }

  if (rem_level>0)
    fprintf(curves_of, "# Note: initial removal of %f was used\n",rem_level);
  fprintf(curves_of, "# Prop_landscape_lost  cost_needed_for_top_fraction  min_prop_rem  ave_prop_rem  W_prop_rem  ext-1  ext-2  prop_for_each_biodiversity_feature_(e.g._species)_remaining_at_level_of_removal...\n");
  pos_w_cnt=0;
  
  const int FIRST_SPP_SHIFT = 5;
  // avoid -0.000, which should only show up by the very end
  for(x=1; x<=(map_cnt+(FIRST_SPP_SHIFT-1)); x++) {
    for (int row=0; row<c_pos; row++) {
      if (curves[row][x] < 0.0f )
	curves[row][x] = 0.0f;
    }
  }
  float wproprem;
  for(y=0; y<c_pos; y++) {
    ext1=ext2=0.0f;
    wsum = 0.0f;
    wproprem=0.0f;
    for(x=FIRST_SPP_SHIFT; x<=(map_cnt+(FIRST_SPP_SHIFT-1)); x++) {
      if (spp[x-FIRST_SPP_SHIFT].weight<0.0f)
	continue;  // xxxMCZ; weighed means now calc for w>0 layers only
      
      wsum     += spp[x-FIRST_SPP_SHIFT].weight;
      wproprem += spp[x-FIRST_SPP_SHIFT].weight*curves[y][x];
      pos_w_cnt++;
      
      if (curves[y][x]>=0.0f) {
	ext1 += (1.0f-z_pow(curves[y][x], SA_z));
	ext2 += spp[x-FIRST_SPP_SHIFT].weight*(1.0f-z_pow(curves[y][x], SA_z));
      } else {
	ext1 += 1.0f;
	ext2 += spp[x-FIRST_SPP_SHIFT].weight;
      }
    }
    if (pos_w_cnt>0)
      ext1 /= pos_w_cnt;  // xxx MCZ mean for w>0 layers only
    
    if (wsum>0.0f) {
      ext2 /= wsum;
      wproprem /= wsum;
    }
    
    fprintf(curves_of, "%0.5f  %-10.5g  %-7.5f %-7.5f %-7.5f   %-7.5f %-7.5f   ",
	    curves[y][0], curves[y][1], curves[y][2], curves[y][3], wproprem,
	    ext1, ext2);
    for(x=FIRST_SPP_SHIFT; x<=(map_cnt+(FIRST_SPP_SHIFT-1)); x++) {
      fprintf(curves_of, "%0.3f ",curves[y][x]);
    }
    fprintf(curves_of, "\n");
  }
  
  fclose(curves_of);
}


void output_SSI_features_info_n_curves()
{
  // A) .SSI_features_info.txt
  FILE* ssi_finfo_of;
  ssi_finfo_of = fopen(SSI_features_info_fn.toUtf8().constData(), "w+t");
  if (!ssi_finfo_of) {
    Form1->Memo1->Lines->Add("************** ERROR ***************");
    Form1->Memo1->Lines->Add("  Error opening SSI feature info data file "+SSI_features_info_fn);
    return;
  }

  float first_Tviol = 1.0f;
  fprintf(ssi_finfo_of, "List of SSI biodiversity features (e.g. species) and weights used in analysis in order of columns:\n");
  fprintf(ssi_finfo_of, "Weight distribution_sum T_violation_fract_remaining map_file_name\n");
  for(int x=0; x<SSI_spp_cnt; x++) {
    fprintf(ssi_finfo_of, "%-6.3f %-6.3f %0.6f %s\n", SSI[x].weight, SSI[x].prob_sum, 
	    1.0f-fabs(SSI[x].T_violation_fract), SSI[x].fname.toUtf8().constData());
    if (SSI[x].T_violation_fract<first_Tviol)
      first_Tviol = SSI[x].T_violation_fract;    
  }

  if (removal_rule==3) {
    fprintf(ssi_finfo_of, "First SSI biodiversity feature target violated at fraction %f of landscape removed, %f remaining\n", 
	    first_Tviol, 1.0f-first_Tviol);
    Form1->Memo1->Lines->Add("First SSI biodiversity feature target violated at fraction of landscape remaining= "+FloatToStr(1.0f-first_Tviol));
  }

  fclose(ssi_finfo_of);


  // B) .SSI_curves.txt
  FILE  *of;
  int   x,y;
  of = fopen(SSI_curvesfn.toUtf8().constData(), "w+t");
  if (!of)
    {
      Form1->Memo1->Lines->Add("************** ERROR ***************");
      Form1->Memo1->Lines->Add("  Error opening SSI curves data file "+curvesfn);
      return;
    }

  if (rem_level>0)
    fprintf(of, "# Note: initial removal of %f was used\n",rem_level);

  fprintf(of, "Prop_landscape_lost  cost_needed_for_top_fraction  min_prop_rem  ave_prop_rem  weighted_ave_prop_rem   prop_remaining_at_level_of_removal_for_each_ssi_feature...\n");
  for(y=0; y<c_pos; y++) {
    // avoid -0.000, which should only show up by the very end, with large datasets
    for(x=0; x<3; x++) {
      if (SSI_curves[y][x]<0.0f)
	SSI_curves[y][x] = 0.0f;
    }
    for(x=1; x<SSI_spp_cnt; x++) {
      for (int row=0; row<c_pos; row++) {
	if (SSI_indiv_curves[row][x] < 0.0f )
	  SSI_indiv_curves[row][x] = 0.0f;
      }
    }
    
    fprintf(of, "%0.4f  %-10.5g  %-6.4f %-6.4f %-6.4f   ",
	    curves[y][0], curves[y][1], SSI_curves[y][0], SSI_curves[y][1], SSI_curves[y][2]);
    for(x=0; x<SSI_spp_cnt; x++) {
      fprintf(of, "%0.3f ", SSI_indiv_curves[y][x]);
    }
    fprintf(of, "\n");
  }

  fprintf(of, "\n\nFeature  lost_at_fraction_removed  file_name\n");
  for(x=0; x<SSI_spp_cnt; x++) {
    if (-1.0f==SSI[x].lost_at_fraction)
      SSI[x].lost_at_fraction = 1.0f;
    fprintf(of, "%-6i %-7.5f %s\n", x, SSI[x].lost_at_fraction, SSI[x].fname.toUtf8().constData());
  }
  fclose(of);
}

// rank
void output_grid_rank(int num, int xsize, int ysize)
{
  float nodatavalue(-1.0f);
  boost::scoped_array<float> plane(new float[xsize * ysize]);
  for(int y = 0; y < ysize; ++y) {
    for(int x = 0; x < xsize; ++x) {
      int i(x + y * xsize);
      //if(Rmax[y][x] == -1)
      if(-1 == status[y][x])
	plane[i] = nodatavalue;
      else
	plane[i] = sol[y][x];
    }
  }
  SaveToRaster<float>(ChangeFileExt(outgridfn), plane.get(), nodatavalue, xsize, ysize);
}

// prop
void output_grid_prop_rank(int num, int xsize, int ysize)
{
  float nodatavalue(-1.0f);
  boost::scoped_array<float> plane(new float[xsize * ysize]);
  for(int y = 0; y < ysize; ++y) {
    for(int x = 0; x < xsize; ++x) {
      int i(x + y * xsize);
      //if(Rmax[y][x] == -1)
      if(-1 == status[y][x])
	plane[i] = nodatavalue;
      else
	plane[i] = sol_val[y][x];
    }
  }
  SaveToRaster(ChangeFileExt(outgridfn2), plane.get(), nodatavalue, xsize, ysize);
}

// wrscr
// Assumes Rmax is still allocated
void output_grid_wrscr(int num, int xsize, int ysize)
{
  float nodatavalue(-1.0f);
  boost::scoped_array<float> plane(new float[xsize * ysize]);
  for(int y = 0; y < ysize; ++y) {
    for(int x = 0; x < xsize; ++x) {
      int i(x + y * xsize);
      if(Rmax[y][x] == -1) { 
	plane[i] = nodatavalue;
      } else {
	plane[i] = 0.0f;
	
	//for(int spn = 0; spn != map_cnt; ++spn) {
	const Biodiv_Features_Occur_Container& rowp = vmat[y][x];
	for(int spn = rowp.first(); spn != rowp.overflow(); spn = rowp.next(spn)) {
	    float val = vmat[y][x][spn];  // xxxMCZ, here includes neg weights.
	    if(val > 0.0f)
	      plane[i] += val * spp[spn].weight;
	}
      }
    }
  }
  SaveToRaster(ChangeFileExt(outgridfn3), plane.get(), nodatavalue, xsize, ysize);
}

// Save the rank using a color bar (traditional Zonation)
// It should be called save_image() or similar (format not necessarily jpg)
void save_jpg()
{
  // Before, the grid of colors (4 byte per cell) was kept in memory
  // from the beginning. 
  //Form1->Image1->Picture->Bitmap->SaveToFiles(ChangeFileExt(emffn, ""));

  // Now: avoid having potentially several GB around all the time
  // generate only here, after other big blocks have likely been freed.
  // like in gridmap.cpp:set_color_for_pixel(...) and
  // GridMap::draw_on_bitmap_smooth(...)
  int xsize = xd;
  int ysize = yd;
  float nodatavalue(-1.0f);
  int** buffer = imatrix(0, yd, 0, xd);

  //SaveToRaster(ChangeFileExt(outgridfn3), plane.get(), nodatavalue, xsize, ysize);
  // Like the old vcl.cpp:Graphics::TBitmap::SaveToFiles(ChangeFileExt(emffn, ""));
  for(int y=0; y<ysize; y++) {
    for(int x=0; x<xsize; x++) {
      int gray = 255.*sol[y][x]; // /1.0f;
      TColor col = (TColor)RGB(255-gray, 255-gray, 255-gray);
      if (sol[y][x] < .0f) {
	col = (TColor)RGB(255, 255, 255);
      } else {
	if (sol[y][x]<0.2f)
	  col =(TColor)RGB(0, 0, 0);
	else if (sol[y][x]<0.5f)
	  col =(TColor)RGB(0, 0, 100);
	else if (sol[y][x]<0.75f)
	  col =(TColor)RGB(0, 0, 255);
	else if (sol[y][x]<0.9f)
	  col =clYellow;
	else if (sol[y][x]<0.95f)
	  col =(TColor)RGB(255, 0, 255);
	else if (sol[y][x]<0.98f)
	  col =(TColor)RGB(150, 0, 0);
	else
	  col =(TColor)RGB(255, 0, 0);
      }

      buffer[y][x] = col;
    }
  }

  save_image_to_file(emffn, xsize, ysize, buffer);
  free_imatrix(buffer, 0, yd, 0, xd);
}

void output_grp_curves(Tgroups_info& groups_info)
{
	FILE  *f;
	int   loop, loop2, num, row, col, pos, cnt;
	float minv, maxv, mean, wmean, wmean_accum, sa_risk, val;

	int max_grp = groups_info.max_grp;  // this is the max grp index!
	int* grpv = groups_info.grpv;

	f=fopen(grp_curves_fn.toUtf8().constData(), "w+t");
	if (!f)
		return;
	fprintf(f, "F-lost\t   TF_cost\t");
	for(pos=0; pos<max_grp; pos++)
	  fprintf(f, "min-%i\tmean-%i\tmax-%i\tw-mean-%i\text2-%i\t", grpv[pos], grpv[pos], grpv[pos], grpv[pos], grpv[pos]);
	fprintf(f, "\n");

	const int FIRST_SPP_SHIFT = 5;
	for(row=0; row<c_pos; row++)
	{
		fprintf(f, "%0.4f\t%-10.5g\t", curves[row][0], curves[row][1]);
		for(pos=0; pos<max_grp; pos++)
		{
			num  = grpv[pos];
			mean = 0.0f;
			cnt  = 0;
			wmean = 0.0f;
			wmean_accum = 0.0f;
			sa_risk = 0.0f;
			maxv = -std::numeric_limits<float>::infinity();
			minv = std::numeric_limits<float>::infinity();
			for(col=FIRST_SPP_SHIFT; col<(map_cnt+FIRST_SPP_SHIFT); col++)
			{
				if (spp[col-FIRST_SPP_SHIFT].grp_op1==num)
				{
					val = curves[row][col];
					mean += val;
					wmean += spp[col-FIRST_SPP_SHIFT].weight * val;
					if (val<minv)
						minv=val;
					if (val>maxv)
						maxv=val;
					cnt++;
					wmean_accum += spp[col-FIRST_SPP_SHIFT].weight;
					sa_risk += spp[col-FIRST_SPP_SHIFT].weight*(1.0f-z_pow(curves[row][col], SA_z));
				}
			}
			if (cnt>0) {
			  mean/=cnt;
			  if (0!=wmean_accum) {
			    wmean/=wmean_accum;
			    sa_risk/=wmean_accum;
			  } else {
			    wmean = sa_risk = 0;
			  }
			} else {
			  minv = mean = wmean = maxv = sa_risk = 0.0f;
			}
			fprintf(f, "%-.4f\t%-.4f\t%-.4f\t%-.4f\t%-.4f\t", minv, mean, maxv, wmean, sa_risk);
		}

		fprintf(f, "\n");
	}

	fclose(f);
}

// Assumes Rmax is still allocated
bool
generate_output_transf_single_layer(size_t bf, float* buff, int xsize, int ysize,
				    const String& odirname, const String& trans_layer_str_append)
{
  float nodatavalue = -1.0f;
  int64_t pos = 0;
  for(int y = 0; y < ysize; ++y) {
    for(int x = 0; x < xsize; ++x) {
      if(Rmax[y][x] < 0 || !(vmat[y][x])) {
	buff[pos] = nodatavalue;
      } else {	
	buff[pos] = vmat[y][x][bf];
      }
      pos++;
    }
  }

  QString layer_name = spp[bf].fname;
  // beware of species like tut_input/species1.asc, I want just 'species1'
  layer_name = ChangeFileExt(layer_name.mid(layer_name.lastIndexOf('/')+1));
  String ofname = odirname + "/feat_" + IntToStr(bf+1) + "_" + layer_name + "_" + trans_layer_str_append;
  SaveToRaster<float>(ofname, buff, nodatavalue, xsize, ysize);

  return true;
}

bool
generate_output_transf_layers(float* buff, const String& o_subdir, String*& trans_layers_str_append, const String must_suffix)
{
  GridMap& map(obsmap[0]);
  int xsize = map.cols(), ysize = map.rows();

  String odirname = createSubDirIfNeeded(o_subdir, "");
  for(size_t bf = 0; bf < map_cnt;  bf++) {
    // skip layers for which no transform has been done
    if ( trans_layers_str_append[bf].isEmpty() || 
	 !must_suffix.isEmpty() && (0 == trans_layers_str_append[bf].endsWith(must_suffix)) )
      continue;

    generate_output_transf_single_layer(bf, buff, xsize, ysize, odirname, trans_layers_str_append[bf]);
  }
  return true;
}

// fedemp 20120213: call to get_grp_max() moved to
// loaddata.cpp:data_input_v2() (this is now needed for
// ADMUs_output_grp_curves_iter() which is called during the iterative
// removal.
// Also, bat_run.cpp now calls output_grp_curves() directly
/*
void grouped_outputs(int max_grp, int *grpv)
{
   int max_grp;
   int grps[MAX_SPP_COUNT];   // ---> dynamic with size map_cnt
   max_grp = get_grp_max(grps);

   output_grp_curves(max_grp, grps);
}
*/
