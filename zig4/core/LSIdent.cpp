#include "VCL.h"
#include "Unit1.h"
#include "GridMap.h"
#include "defines.h"
#include "matrix_utils.h"
#include "bat_run.h"
#include "LoadData.h"
#include <cstdio>
#include <ctime>
#include <cstring>
#include <cmath>

// special LSM mode with "-1 distance" = no aggregation of spots, just statistics per unit
bool LSM_minus1_mode = false;

float   top_percent, min_percent, max_dist, min_simil;
String  LSIRfn, LSIDfn, LSImask;

int    spot_cnt, nwc, ccnt;
int    **cm=0, **nwm, nwns_cnt;
std::vector<struct spot> spots;
std::vector<int> s_in_nw;
std::vector<bool> nw_in_min;
std::vector<float> Ebd;

class  GridMap LSI_maskmap;

void get_candidate_cells()
{
	float val_th;
	int   x, y;

	val_th = 1.0f-top_percent;

	ccnt=0;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (sol[y][x]>=val_th)
			{
				cm[y][x]=0;
				ccnt++;
			}
			else
				cm[y][x]=-1;
		}
	}

	Form1->Memo1->Lines->Add("Potential cells count = "+IntToStr(ccnt));
}

bool nb_in_s(int x, int y, int s)
{
	int minx, miny, maxx, maxy, gx,gy;

	minx=max(0, x-1);
	miny=max(0, y-1);
	maxx=min(xd,x+1);
	maxy=min(yd,y+1);

	for(gy=miny; gy<=maxy; gy++)
	{
		for(gx=minx; gx<=maxx; gx++)
			if (cm[gy][gx]==s)
				return true;
	}

	return false;
}

bool add_to_spot(int s)
{
	int  minx, maxx, miny, maxy, x, y, loop;
	bool added;
	float val, *rowp;

	minx=max(0, spots[s].min_gx-1);
	miny=max(0, spots[s].min_gy-1);
	maxx=min(xd,spots[s].max_gx+1);
	maxy=min(yd,spots[s].max_gy+1);

	added=false;
	for(y=miny; y<=maxy; y++)
	{
		for(x=minx; x<=maxx; x++)
		{
			//         if (cm[y][x]==0)
			if ((cm[y][x]!=s) && (cm[y][x]!=-1))
			{
				//             Form1->Memo1->Lines->Add("HERE");
				if (nb_in_s(x,y,s))
				{
					//               Form1->Memo1->Lines->Add("YES");
					cm[y][x]=s;
					spots[s].area++;
					spots[s].rank += sol[y][x];
					spots[s].mean_gx += x;
					spots[s].mean_gy += y;
					spots[s].min_gx=min(spots[s].min_gx, x);
					spots[s].min_gy=min(spots[s].min_gy, y);
					spots[s].max_gx=max(spots[s].max_gx, x);
					spots[s].max_gy=max(spots[s].max_gy, y);

					if (sol[y][x]>(1.0f-min_percent))
						spots[s].in_min_percent=true;
					// rowp=&vmat[y][x][0];  // COMPACT_VMAT
					Biodiv_Features_Occur_Container& rowp = vmat[y][x];
					for(loop=0; loop<map_cnt; loop++)
					{
					  //std::cerr << "rowp: " << rowp.size() << std::endl;
					  val = rowp[loop];
					  if (val!=-1)
					    spots[s].bdv[loop] += val;
					  // else
					  //  spots[s].bdv[loop] = 0.0f;
					}
					added=true;
				}
			}
		}
	}

	return added;
}

void expand_spot(int x, int y)
{
	bool added;
	int loop;
	
	spots[spot_cnt].bdv     = 0;
	spots[spot_cnt].min_gx  = x;
	spots[spot_cnt].min_gy  = y;
	spots[spot_cnt].max_gx  = x;
	spots[spot_cnt].max_gy  = y;
	spots[spot_cnt].mean_gx = static_cast<float>(x);
	spots[spot_cnt].mean_gy = static_cast<float>(y);
	spots[spot_cnt].area    = 1;
	spots[spot_cnt].rank    = sol[y][x];
	if (sol[y][x]>=(1.0f-min_percent))
		spots[spot_cnt].in_min_percent=true;
	else
		spots[spot_cnt].in_min_percent=false;
	spots[spot_cnt].bdv = new float[map_cnt];

	for(loop=0; loop<map_cnt; loop++)
		spots[spot_cnt].bdv[loop] =0.0f;

	cm[y][x]=spot_cnt;
	do {
		added = add_to_spot(spot_cnt);
		Application->ProcessMessages();
	} while(added);
#if 0
	char txt[255];
	sprintf(txt,"Spot %i A=%i Xmin=%i xmax=%i ymin=%i ymax=%i  mean-x=%0.3f  mean-y=%0.3f",
		spot_cnt, spots[spot_cnt].area, spots[spot_cnt].min_gx,
		spots[spot_cnt].max_gx, spots[spot_cnt].min_gy,
		spots[spot_cnt].max_gy,
		spots[spot_cnt].mean_gx, spots[spot_cnt].mean_gy);
	//  if ((spot_cnt%10)==0)
	Form1->Memo1->Lines->Add(txt);
#endif
}

void get_spots()
{
	float val_th;
	int   x, y, in_cnt;

	spot_cnt=1;
	const size_t DEF_MAX_SPOTS = 2048;
	spots.reserve(DEF_MAX_SPOTS);
	try {
	  spots.resize(spot_cnt+1);
	} catch(std::bad_alloc const& ex) {
	  Form1->Memo1->Lines->Add("Out of memory in landscape identification: "+String(ex.what()));
	}
	in_cnt  =0;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (cm[y][x]==0)
			{
				Application->ProcessMessages();
				expand_spot(x,y);
				if (spots[spot_cnt].in_min_percent)
					in_cnt++;
				//               Form1->Memo1->Lines->Add("New spot, area = "
				//                 +IntToStr(spots[spot_cnt].area));
				spot_cnt++;
				spots.resize(spots.size()+1);

				if ((spot_cnt%1000)==0)
					Form1->Memo1->Lines->Add("Spot count = "+IntToStr(spot_cnt-1));
			}
		}
	}
	Form1->Memo1->Lines->Add("Spot count = "+IntToStr(spot_cnt-1));
	Form1->Memo1->Lines->Add("Spots including best areas count = "+IntToStr(in_cnt));
}

float calc_dist(int s1, int s2)
{
	float dij, dx, dy, dm2;
	float minx1, minx2, maxx1, maxx2, miny1, miny2, maxy1, maxy2;
	int   x1, x2, y1, y2;

	dm2 = max_dist*max_dist;
	if (dm2==0)
		return (max_dist+1.0f); // with zero dist, separate spots cannot be joined

	minx1=static_cast<float>(spots[s1].min_gx);
	maxx1=static_cast<float>(spots[s1].max_gx);
	miny1=static_cast<float>(spots[s1].min_gy);
	maxy1=static_cast<float>(spots[s1].max_gy);

	minx2=static_cast<float>(spots[s2].min_gx);
	maxx2=static_cast<float>(spots[s2].max_gx);
	miny2=static_cast<float>(spots[s2].min_gy);
	maxy2=static_cast<float>(spots[s2].max_gy);

	//  Form1->Memo1->Lines->Add("corners");
	//  sprintf(txt, "minx1=%f  maxx1=%f  miny1=%f  maxy1=%f", minx1, maxx1, miny1, maxy1);
	//  Form1->Memo1->Lines->Add(txt);
	//  sprintf(txt, "minx2=%f  maxx2=%f  miny2=%f  maxy2=%f", minx2, maxx2, miny2, maxy2);
	//  Form1->Memo1->Lines->Add(txt);
	//            Form1->Memo1->Lines->Add("yxxxvc");

	if (minx1>(maxx2+max_dist))
		return (max_dist+1.0f);
	if (minx2>(maxx1+max_dist))
		return (max_dist+1.0f);
	if (miny1>(maxy2+max_dist))
		return (max_dist+1.0f);
	if (miny2>(maxy1+max_dist))
		return (max_dist+1.0f);

	//            Form1->Memo1->Lines->Add("here");
	for(y1=static_cast<int>(miny1); y1<=maxy1;y1++)
		for(x1=static_cast<int>(minx1); x1<=maxx1;x1++)
		{
			//            Form1->Memo1->Lines->Add("y1loop"+IntToStr(y1));
			if (cm[y1][x1]==s1)
			{
				for(y2=static_cast<int>(miny2); y2<=maxy2;y2++)
					for(x2=static_cast<int>(minx2); x2<=maxx2;x2++)    // xxx stuck in this loop.
					{
						if (cm[y2][x2]==s2)
						{
							dij = z_pow(x1-x2,2)+z_pow(y1-y2,2);
							if (dij<=dm2)
								return 0.0f;
						}
					}
			}
		}

	return (max_dist+1.0f);
}

float calc_sim(int s1, int s2)
{
	int loop, s, lvl1, lvl2;
	float diff, val;

	diff=0.0f;
	for(loop=0; loop<map_cnt; loop++)
	{
		val = spots[s1].bdv[loop];
		if (val<0.01f*Ebd[loop])
			lvl1=0;
		else if (val<0.1f*Ebd[loop])
			lvl1=1;
		else if (val<Ebd[loop])
			lvl1=2;
		else if (val<10*Ebd[loop])
			lvl1=3;
		else if (val<100*Ebd[loop])
			lvl1=4;
		else
			lvl1=5;

		val = spots[s2].bdv[loop];
		if (val<0.01f*Ebd[loop])
			lvl2=0;
		else if (val<0.1f*Ebd[loop])
			lvl2=1;
		else if (val<Ebd[loop])
			lvl2=2;
		else if (val<10*Ebd[loop])
			lvl2=3;
		else if (val<100*Ebd[loop])
			lvl2=4;
		else
			lvl2=5;

		diff += fabs(lvl1-lvl2);
	}
	diff/=map_cnt;

	return diff;
}

void add_nb_to_nw(int spot)
{
	int   loop;
	float dist, diff;

	for(loop=1; loop<spot_cnt; loop++)
	{
		//  Form1->Memo1->Lines->Add("spot = "+IntToStr(loop));
		if (spots[loop].nwn==-1)
		{
			//  Form1->Memo1->Lines->Add("dist next");
			dist = calc_dist(spot, loop);
			//  Form1->Memo1->Lines->Add("sim next");
			diff  = calc_sim(spot, loop);
			if ((dist<=max_dist) && (diff<min_simil))
			{
				spots[loop].nwn = nwc;
				nw_in_min[nwc] = (nw_in_min[nwc] || spots[loop].in_min_percent);
				s_in_nw[nwns_cnt] = loop;
				nwns_cnt++;
				//              Form1->Memo1->Lines->Add("joined");
			}
		}
	}
}

void expand_network(int s)
{
	int pos;

	spots[s].nwn = nwc;
	nwns_cnt     = 1;
	s_in_nw[0]   = s;
	nw_in_min[nwc] = (nw_in_min[nwc] || spots[s].in_min_percent);

	pos=0;
	while(pos<nwns_cnt)
	{
		//  Form1->Memo1->Lines->Add("pos="+IntToStr(pos));
		add_nb_to_nw(s_in_nw[pos]);
		pos++;
	}
}

void  Fix_bd_values()
{
	int  loop, s;

	for(loop=0; loop<map_cnt; loop++)
		Ebd[loop]=0.0f;

	for(s=1;s<spot_cnt;s++)
	{
		for(loop=0; loop<map_cnt; loop++)
		{
			Ebd[loop] += spots[s].bdv[loop];
			spots[s].bdv[loop] /= spots[s].area;
		}
	}

	for(loop=0; loop<map_cnt; loop++)
		Ebd[loop] /= ccnt; // Ebd[loop]=average over cells in cut
}

void  get_networks()
{
	int   loop, x, y, cilcnt;
	/*
	const size_t DEF_MAX_NW = 2048;
	nw_in_min.reserve(DEF_MAX_NW);
	nw_in_min.assign(2, false);
	*/
	// get_networks always called after get_spots() -> spot_cnt known
	try {
	  s_in_nw.resize(spot_cnt, 0);
	  Ebd.resize(map_cnt, 0);
	  nw_in_min.resize(spot_cnt+1, false);
	} catch(std::bad_alloc const& ex) {
	  Form1->Memo1->Lines->Add("Out of memory in landscape identification: "+String(ex.what()));
	}

	nwc = 1; // xxx
	for(loop=0; loop<spot_cnt; loop++) {
		spots[loop].nwn=-1;
	}

	// uses Ebd[]
	Fix_bd_values();

	for(loop=1; loop<spot_cnt; loop++)
	{
		if ((spots[loop].nwn==-1) && (spots[loop].in_min_percent))
		{
			//  Form1->Memo1->Lines->Add(IntToStr(loop));
			expand_network(loop);
			nwc++;
			nw_in_min.push_back(false);
		}
	}

	for(y=0; y<yd; y++) {
	  for(x=0; x<xd; x++) {
	    if (false && cm[y][x] >0 && spots[cm[y][x]].nwn >= nw_in_min.size()) {
	      Form1->Memo1->Lines->Add("SEGFAUUUUUUUUUUUUUUUUUULT, cm: "+ IntToStr(cm[y][x]));
	      Form1->Memo1->Lines->Add("SEGFAUUUUUUUUUUUUUUUUUULT, cm: "+ IntToStr(cm[y][x]) +  ", spots:"+IntToStr(spots[cm[y][x]].nwn));
	    }
	    //if (Rmax[y][x]==-1)
	    if (-1 == status[y][x])
	      nwm[y][x]=-1;
	    else if (cm[y][x]==0)
	      nwm[y][x]=-2;
	    else if (cm[y][x]==-1)
	      nwm[y][x]=-2;
	    // the -1== check is critical, or likely segfault!
	    else if ((cm[y][x]>0) && (-1==spots[cm[y][x]].nwn || (!nw_in_min[spots[cm[y][x]].nwn])))
	      nwm[y][x]=-2;
	    else
	      nwm[y][x]=spots[cm[y][x]].nwn; // xxx error
	    }
	}

	Form1->Memo1->Lines->Add("Found networks count = "+IntToStr(nwc-1));

	std::vector<float> spdat;
	try {
	  spdat.resize(map_cnt, 0.0);
	} catch(std::bad_alloc const& ex) {
	  Form1->Memo1->Lines->Add("Out of memory in landscape identification: "+String(ex.what()));
	}
	cilcnt=0;
	for(y=0; y<yd; y++) {
	  for(x=0; x<xd; x++) {
	    //if (Rmax[y][x]==-1)
	    if (-1 == status[y][x])
	      continue;
	    if (nwm[y][x]>0) {
	      // rowp=&vmat[y][x][0];  // COMPACT_VMAT
	      Biodiv_Features_Occur_Container& rowp = vmat[y][x];
	      //for(loop=0;loop<map_cnt;loop++)
	      for(loop=rowp.first(); loop!=rowp.overflow(); loop=rowp.next(loop)) {
		if (rowp[loop]>0.0f)
		  spdat[loop] += rowp[loop];
	      }
	      cilcnt++;
	    }
	  }
	}
	Form1->Memo1->Lines->Add("Cells in classified landscapes = "+IntToStr(cilcnt));

	const size_t MAX_STR_LEN = 512;
	char  txt[MAX_STR_LEN];
	for(loop=0;loop<map_cnt;loop++)
	{
		sprintf(txt, "%-6.3f %-5.3f %s\n", spp[loop].weight, spdat[loop], spp[loop].fname.toUtf8().constData());
		Form1->Memo1->Lines->Add(txt);
	}
}

// Gets statistics for units specified in the PPA mask, 
// so networks will actually be the units 
// all the nwarea, nwx, nwy, nwrank are aggregated (not normalized) and will be divided by nwarea[] later on
void
get_fake_networks_from_mask(std::vector<int>& nwarea, std::vector<float>& nwx, std::vector<float>& nwy, 
			    std::vector<float>& nwrank, float**& mat)
{
  // max_val pre-calculated in load_from_file...
  spot_cnt = (int)LSI_maskmap.max_val+1;
  try {
    spots.resize(spot_cnt);
  } catch(std::bad_alloc const& ex) {
    Form1->Memo1->Lines->Add("Out of memory in landscape identification: "+String(ex.what()));
  }

  /*
  // Init to 0. units have numbers >=1
  for (size_t i=0; i<spot_cnt; i++)
    spots[i].num = 0;

  // Find used unit numbers. Use the array spots[].num 
  // for the unit numbers
  for (size_t y=0; y<yd; y++) {
      for (size_t x=0; x<xd; x++) {
	// make sure the mask doesn't include "missing" cells
	if (sol[y][x] < 0.0f)
	  continue;
	int unit_idx = LSI_maskmap.m[y][x];
	if (unit_idx > 0 && unit_idx <= spots.size())
	  if (0==spots[unit_idx].num)
	  spots[unit_idx].num++;
      }
  }
  // unit numbers actually used in the LSI analysis mask
  std::vector<int> unit_nums;
  nwc = 0;
  for (size_t i=0; i<spot_cnt; i++) {
    if (0 < spots[i].num) {
      nwc++;     // nwc is global
      unit_nums.push_back(i);
    }
  }
  nwc++; // yes, it is number of networks/units +1
  */
  // bypass the 2 loops above. Use as many networks as the biggest number in the planning
  // units layer/mask. This avoids crashes if non-consecutive numbers are used.
  nwc = spot_cnt;

  try {
    nwarea.resize(nwc+1, 0);
    nwrank.resize(nwc+1, 0);
    nwrank.resize(nwc+1, 0);
    nwx.resize(nwc+1, 0);
    nwy.resize(nwc+1, 0);	  
  } catch(std::bad_alloc const& ex) {
    Form1->Memo1->Lines->Add("Out of memory in landscape identification: "+String(ex.what()));
  }
  mat = matrix(0, nwc+1, 0, map_cnt);
  if (!mat) {
    ShowMessage("Out of memory when doing LSIdent");
    return;
  }
  for(int nw=0; nw<=nwc; nw++) {
    nwarea[nw]=0;
    nwrank[nw]=0.0f;
    nwx[nw]=0.0f;
    nwy[nw]=0.0f;
    for(int sp=0; sp<map_cnt; sp++)
      mat[nw][sp]=0.0f;
  }

  for (size_t y=0; y<yd; y++) {
    for (size_t x=0; x<xd; x++) {
      // make sure the mask doesn't include "missing" cells
      if (sol[y][x] < 0.0f)
	continue;

      int unit_idx = LSI_maskmap.m[y][x];

      if (unit_idx <= 0 || unit_idx >= spot_cnt)
	continue;

      nwarea[unit_idx]++;
      nwx[unit_idx] += x;
      nwy[unit_idx] += y;
      nwrank[unit_idx] += sol[y][x];

      // float* rowp = &vmat[y][x][0];  // COMPACT_VMAT
      Biodiv_Features_Occur_Container& rowp = vmat[y][x];
      if (rowp)
	//for(size_t spp_idx=0; spp_idx<map_cnt; spp_idx++)
	for(size_t spp_idx=rowp.first(); spp_idx!=rowp.overflow(); spp_idx=rowp.next(spp_idx))
	  if (rowp[spp_idx] > .0f)
	    mat[unit_idx][spp_idx] += rowp[spp_idx];
    }
  }
  // And the nwout raster (variable nwm) is not generated (it's = LSImask)
}

void print_network_data()
{
	int   loop, nw, sp, num, c10, c1, c01, c001, c0001, sp_at_zero;
	float **mat, nwtot;
	FILE  *f;

	f=fopen(LSIDfn.toUtf8().constData(), "w+t");
	if (!f)
	{
		ShowMessage("Could not open output file " + LSIDfn);
		return;
	}


	std::vector<int> nwarea;
	std::vector<float> nwrank, nwx, nwy;
	
	if (LSM_minus1_mode) {
	  // All params by ref./output
	  get_fake_networks_from_mask(nwarea, nwx, nwy, nwrank, mat);
	} else {
	  try {
	    nwarea.resize(nwc+1, 0);
	    nwrank.resize(nwc+1, 0);
	    nwrank.resize(nwc+1, 0);
	    nwx.resize(nwc+1, 0);
	    nwy.resize(nwc+1, 0);	  
	  } catch(std::bad_alloc const& ex) {
	    Form1->Memo1->Lines->Add("Out of memory in landscape identification: "+String(ex.what()));
	  }	  
	  //  sprintf(txt, "nwc=%i spots=%i",nwc, spot_cnt);
	  //  ShowMessage(txt);
	  mat = matrix(0, nwc+1, 0, map_cnt);
	  if (!mat)
	    {
	      ShowMessage("Out of memory when doing LSIdent");
	      fclose(f);
	      return;
	    }
	  for(nw=0; nw<=nwc; nw++)
	    {
	      nwarea[nw]=0;
	      nwrank[nw]=0.0f;
	      nwx[nw]=0.0f;
	      nwy[nw]=0.0f;
	      for(sp=0; sp<map_cnt; sp++)
		mat[nw][sp]=0.0f;
	    }
	  for(loop=1; loop<spot_cnt; loop++)
	    {
	      num = spots[loop].nwn;
	      if (num == -1)
		continue;
	      
	      for(sp=0; sp<map_cnt; sp++)
		mat[num][sp] += spots[loop].bdv[sp]*spots[loop].area;
	      nwarea[num] += spots[loop].area;
	      nwrank[num] += spots[loop].rank;
	      nwx[num]    += spots[loop].mean_gx;
	      nwy[num]    += spots[loop].mean_gy;
	    }
	}
	
	
	std::string nets_or_units;
	if (LSM_minus1_mode)
	  nets_or_units = "units";
	else 
	  nets_or_units = "networks";
	fprintf(f, "Most important biodiversity features (e.g. species) in %s; those occurring at a 1%%+ level\n", nets_or_units.c_str());
	fprintf(f, "of original distribution\n");

	std::string net_or_unit;
	if (LSM_minus1_mode)
	  net_or_unit = "Unit";
	else 
	  net_or_unit = "Network";
	fprintf(f, "%s  Area  Mean-Rank  X   Y  Spp_distribution_sum  spp occurring at >10%%  >1%%  >0.1%%  >0.01%% >0.001%%\n", net_or_unit.c_str());

	std::vector<float> sptot;
	try {
	  sptot.resize(map_cnt, 0);
	} catch(std::bad_alloc const& ex) {
	  Form1->Memo1->Lines->Add("Out of memory in landscape identification: "+String(ex.what()));
	}
	float tottot=0.0f;
	for(nw=1; nw<nwc; nw++)
	{
	  // do not calc/output results for empty/missing unit numbers
	  if (LSM_minus1_mode && nwarea[nw]<=0)
	    continue;

		c10 = c1 = c01 = c001 = c0001 = 0;
		nwtot=0.0f;
		for(sp=0; sp<map_cnt; sp++)
		{
			nwtot += mat[nw][sp];
			sptot[sp] += mat[nw][sp];
			if (mat[nw][sp]>0.1f)
			{
				c10++;
				c1++;
				c01++;
				c001++;
				c0001++;
			}
			else if (mat[nw][sp]>0.01f)
			{
				c1++;
				c01++;
				c001++;
				c0001++;
			}
			else if (mat[nw][sp]>0.001f)
			{
				c01++;
				c001++;
				c0001++;
			}
			else if (mat[nw][sp]>0.0001f)
			  {
				c001++;
				c0001++;
			  }
			else if (mat[nw][sp]>0.00001f)
			  {
				c0001++;
			  } 
		}
		tottot += nwtot;

		//           nw   area  rnk    x      y      tot
		fprintf(f, "%-5i %-6i %-6.3f %-6.3f %-6.3f %-6.3f %-5i %-5i %-5i %-5i %-5i\n",
			nw, nwarea[nw], nwrank[nw]/nwarea[nw], nwx[nw]/nwarea[nw], nwy[nw]/nwarea[nw],
			nwtot, c10, c1, c01, c001, c0001);
		for(sp=0; sp<map_cnt; sp++)
		{
			if (mat[nw][sp]>0.01f)
				fprintf(f, "                    Feature %s, %-5.2f%% of full distribution\n",
					spp[sp].fname.toUtf8().constData(),100*mat[nw][sp]);
		}
		//      fprintf(f, "\n");
	}

	fprintf(f, "Repeat without spp info for easy import\n");
	fprintf(f,"%s  Area  Mean-Rank  X   Y  Spp_distribution_sum  spp occurring at >10%%  >1%%  >0.1%%  >0.01%% >0.001%%\n", net_or_unit.c_str());

	//tottot=0.0f;
	for(nw=1; nw<nwc; nw++)
	{
	  // do not cal/output results for empty/missing unit numbers
	  if (LSM_minus1_mode && nwarea[nw]<=0)
	    continue;

		c10 = c1 = c01 = c001 = c0001 = 0;
		nwtot=0.0f;
		for(sp=0; sp<map_cnt; sp++)
		{
			nwtot += mat[nw][sp];
			// would be the second time!
			//sptot[sp] += mat[nw][sp];
			if (mat[nw][sp]>0.1f)
			{
				c10++;
				c1++;
				c01++;
				c001++;
				c0001++;
			}
			else if (mat[nw][sp]>0.01f)
			{
				c1++;
				c01++;
				c001++;
				c0001++;
			}
			else if (mat[nw][sp]>0.001f)
			{
				c01++;
				c001++;
				c0001++;
			}
			else if (mat[nw][sp]>0.0001f)
			  {
				c001++;
				c0001++;
			  } else if (mat[nw][sp]>0.00001f)
			  {
				c0001++;
			  }
		}
		// this would make tottot 2x the true totot
		//tottot += nwtot;

		//           nw   area  rnk    x      y      tot
		fprintf(f, "%-5i %-6i %-6.3f %-6.3f %-6.3f %-6.3f %-5i %-5i %-5i %-5i %-5i\n",
			nw, nwarea[nw], nwrank[nw]/nwarea[nw], nwx[nw]/nwarea[nw], nwy[nw]/nwarea[nw],
			nwtot, c10, c1, c01, c001, c0001);
	}

	fprintf(f, "\n\nAverage proportion remaining over all spp in %s = %f\n",nets_or_units.c_str(), tottot/map_cnt);

	sp_at_zero=0;
	for(sp=0; sp<map_cnt; sp++)
	{
		if (sptot[sp]<=0.0f)
			sp_at_zero++;
	}
	fprintf(f, "Count of biodiversity features (e.g. species) with nothing remaining in the network = %i\n",sp_at_zero);

	fprintf(f, "Total proportion and sum remaining for biodiversity features\n");
	for(sp=0; sp<map_cnt; sp++)
	{
		fprintf(f, "%s  %-5.4f %0.4f\n",
			spp[sp].fname.toUtf8().constData(), sptot[sp], sptot[sp]*spp[sp].prob_sum);
	}

	if (LSM_minus1_mode)
	  fprintf(f, "\n\nBiological data of %i %s.\n",nwc-1, nets_or_units.c_str());
	else
	  fprintf(f, "\n\nBiological data of %i %s (spots=%i).\n",nwc-1, nets_or_units.c_str(), spot_cnt-1);
	fprintf(f, "%s x biodiversity features matrix\n", nets_or_units.c_str());

	if (LSM_minus1_mode)
	  fprintf(f, "Unit_number  area[cells]  sp_data .....\n");
	else
	  fprintf(f, "Nw_number  area[cells]  sp_data .....\n");
	for(nw=1; nw<nwc; nw++)
	{
	  // do not calc/output results for empty/missing unit numbers
	  if (LSM_minus1_mode && nwarea[nw]<=0)
	    continue;

		fprintf(f, "%-5i %-6i ", nw, nwarea[nw]);
		for(sp=0; sp<map_cnt; sp++)
			fprintf(f,"%-6.4f ", mat[nw][sp]);
		fprintf(f, "\n");
	}

	fclose(f);
	free_matrix(mat, 0, nwc+1, 0, map_cnt);
}

bool read_LSI_mask(int top_fraction_mode)
{
	int x, y;

	LSI_maskmap.normalize=false;
	if (!LSI_maskmap.load_from_file(LSImask, mask_data, area_mask.m))
	{
		Form1->Memo1->Lines->Add("************** ERROR ***************");
		Form1->Memo1->Lines->Add("  FAILURE attempting LSI mask map load.");
		return false;
	}
	Form1->Memo1->Lines->Add("LSI mask map loaded.");

	val_th = 1.0f-top_percent;
	ccnt = 0;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (LSI_maskmap.m[y][x]>=1)
			{
				if (top_fraction_mode)
				{
					if (sol[y][x]>=val_th)
					{
						cm[y][x]=0;
						ccnt++;
					}
					else
						cm[y][x]=-1;
				}
				else
				{
					cm[y][x]=0;
					ccnt++;
				}
			}
			else
				cm[y][x]=-1;
		}
	}

	Form1->Memo1->Lines->Add("Potential cells count = "+IntToStr(ccnt));

	return true;
}

int LSIdent(int LSI_mode)
{
	const size_t MAX_STRLEN = 2048;
	char txt[MAX_STRLEN];

	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("NEW LANDSCAPE IDENTIFICATION ANALYSIS");
	// This is done now in the visitor class in post_process.cpp
	//top_percent = StrToFloat(Form1->Edit4->Text)/100.0f;
	//min_percent = StrToFloat(Form1->Edit5->Text)/100.0f;
	//max_dist    = StrToFloat(Form1->Edit6->Text);
	//min_simil   = StrToFloat(Form1->Edit7->Text);
	//LSIRfn      = Form1->Edit8->Text;
	
	bool lsi_mask_ok = true;

	if (LSI_mode==0) // "LSB"
	{
	  sprintf(txt, "Running LSIdent with top%%=%0.3f min%%=%0.3f max-d=%0.3f min-s=%0.3f", 
		  top_percent*100, min_percent*100, max_dist, min_simil);
	  Form1->Memo1->Lines->Add(txt);
	  Form1->Memo1->Lines->Add("1. Getting candidate cells.");
	  get_candidate_cells();
	}
	else if (LSI_mode==1) // "LSM"
	{
	  if (0.0f > max_dist)
	    {
	      LSM_minus1_mode = true;
	      sprintf(txt, "Running LSIdent with mask file (%s). Note: LSM special case with max. distance -1, ignoring top%%=%0.3f and using whole landscape", 
		      LSImask.toStdString().c_str(), top_percent*100);
	      Form1->Memo1->Lines->Add(txt);
	      top_percent = 1.0f;
	    } 
	  else 
	    {
	      sprintf(txt, "Running LSIdent with mask file (%s) max-d=%0.3f min-s=%0.3f",
		      LSImask.toStdString().c_str(), max_dist, min_simil);
	      Form1->Memo1->Lines->Add(txt);
	      Form1->Memo1->Lines->Add("1. Reading relevant areas from mask file "+Form1->Edit42->Text);
	    }

	    lsi_mask_ok = read_LSI_mask(0);	      
	    if (!lsi_mask_ok) {
	      Form1->Memo1->Lines->Add("ERROR! failed to read LSM areas from mask file: " + LSImask);
	    }
	}
	else // 2==LSI_mode "LSB"
	{
	  sprintf(txt, "Running LSIdent for top fraction within masked area, mask file %s fract=%0.4f  max-d=%0.3f min-s=%0.3f",
		  LSImask.toStdString().c_str(), top_percent, max_dist, min_simil);
	  Form1->Memo1->Lines->Add(txt);
	  Form1->Memo1->Lines->Add("1. Reading relevant areas from mask file "+Form1->Edit42->Text);
	  lsi_mask_ok = read_LSI_mask(1);
	  if (!lsi_mask_ok) {
	    Form1->Memo1->Lines->Add("ERROR! failed to read LSB areas from mask file: " + LSImask);
	  }
	}

	if (!lsi_mask_ok) {
	  Form1->Memo1->Lines->Add("Please fix the mask file name. No results will be generated for this post-processing analysis.");
	  return false;
	}

	if (!LSM_minus1_mode) {
	  // free only if traditional modes
	  LSI_maskmap.free_matrix_m();

	  Form1->Memo1->Lines->Add("2. Identifying spots.");
	  get_spots(); // spots[s].bdv[sp] = 0.0f; sisaltaa prop of sp spotissa

	  Form1->Memo1->Lines->Add("3. Finding networks.");
	  get_networks();
	}

	print_network_data();

	if (LSM_minus1_mode) {
	  LSI_maskmap.free_matrix_m();
	} else {
#if 0
	  //  obsmap[0].show_spots(Form1->Image1->Picture->Bitmap,  cm);
	  obsmap[0].show_spots(Form1->Image1->Picture->Bitmap,  nwm, 0);
#endif

	  obsmap[0].export_GIS_INT(nwm, LSIRfn);

	  // the spots[].bdv are not allocated in LSM "-1 distance" mode
	  for(int loop=0; loop<spots.size(); loop++)
	    {
	      if (spots[loop].bdv)
		delete[] spots[loop].bdv;
	      
	      spots[loop].bdv=0;
	    }
	}

	Screen->Cursor=crDefault;

	return true;
}

void LSCAnalysis(float f1, float f2, const String& cfn, const String& comp_outfn)
{
	//float f1, f2;
	//String cfn, comp_outfn;
	float  f1cells, f2cells, bothcells, rodiff;
	class  GridMap cmpmap;
	int    x, y, z1, z2, rcnt;
	bool   f1ok, f2ok;

	DecimalSeparator='.';
	//cfn = Edit23->Text;

	cmpmap.set_no_normalize();

	if (!cmpmap.load_from_file(cfn, mask_data, area_mask.m))
	{
		ShowMessage("Could not load given comparison solution");
		return;
	}

	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("Solution comparison stats");
	//f1 = StrToFloat(Edit22->Text);
	//f2 = StrToFloat(Edit24->Text);
	Form1->Memo1->Lines->Add("S1 cut level = "+FloatToStrF(f1, ffFixed, 7, 4));
	Form1->Memo1->Lines->Add("S2 cut level = "+FloatToStrF(f2, ffFixed, 7, 4));

	f1cells=f2cells=bothcells=rodiff=0.0f;
	z1=z2=rcnt=0;
	for(y=0; y<yd; y++)
		for(x=0; x<xd; x++)
		{
			f1ok=f2ok=false;

			if (f1>0.0f)
			{
				if ((sol[y][x]!=-1) && (sol[y][x]>=(1.0f-f1)))
					f1ok=true;
			}
			else
			{
				if ((sol[y][x]!=-1) && (sol[y][x]<=(-f1)))
					f1ok=true;
			}

			if (f2>0.0f)
			{
				if ((cmpmap.m[y][x]!=-1) && (cmpmap.m[y][x]>=(1.0f-f2)))
					f2ok=true;
			}
			else
			{
				if ((cmpmap.m[y][x]>0.0f) && (cmpmap.m[y][x]<=(-f2)))
					f2ok=true;
			}

			if (f1ok)
				f1cells++;
			if (f2ok)
				f2cells++;
			if (f1ok && f2ok)
				bothcells++;
			if (sol[y][x]==0.0f)
				z1++;
			if (cmpmap.m[y][x]==0.0f)
				z2++;
			if ((sol[y][x]!=-1) && (cmpmap.m[y][x]!=0.0f))
			{
				++rcnt;
				rodiff+= fabs(sol[y][x]-cmpmap.m[y][x]);
			}

			nwm[y][x]  = 0;
			//if (Rmax[y][x]==-1)
			if (-1 == status[y][x])
				nwm[y][x] = -1;
			else if (f1ok && f2ok)
				nwm[y][x] = 1;
			else if (f1ok)
				nwm[y][x] = 2;
			else if (f2ok)
				nwm[y][x] = 3;
		}
	Form1->Memo1->Lines->Add("Cells in present solution fraction = "      +IntToStr((int)f1cells));
	Form1->Memo1->Lines->Add("Cells in comparison solution fraction = "   +IntToStr((int)f2cells));
	Form1->Memo1->Lines->Add("Cells included in both solutions = "        +IntToStr((int)bothcells));
	Form1->Memo1->Lines->Add("Initially removed in present solution = "   +IntToStr(z1));
	Form1->Memo1->Lines->Add("Initially removed in comparison solution = "+IntToStr(z2));
	Form1->Memo1->Lines->Add("Similarity f1 = "+FloatToStrF(bothcells/f1cells, ffFixed, 7, 4));
	Form1->Memo1->Lines->Add("Similarity f2 = "+FloatToStrF(bothcells/f2cells, ffFixed, 7, 4));
	Form1->Memo1->Lines->Add("Average difference in removal order = "+FloatToStrF(rodiff/rcnt, ffFixed, 7, 4));

	const size_t MAX_STR_LEN =  512;
	char   txt[MAX_STR_LEN];
	sprintf(txt, "Overlap f1 = %0.4f,  f1 = %0.4f. Average order diff=%0.4f. See also memo.",
		bothcells/f1cells,bothcells/f2cells, rodiff/rcnt);
	Form1->Memo1->Lines->Add(txt);
	if (!bat_mode)
		ShowMessage(txt);

	if (Form1->CheckBox7->Checked)
	{
#if 0
		//comp_outfn = Edit29->Text;
		obsmap[0].show_spots(Form1->Image1->Picture->Bitmap,  nwm, 1);
#endif
		obsmap[0].export_GIS_INT(nwm, comp_outfn);
	}
}
