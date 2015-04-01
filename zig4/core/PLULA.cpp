#include "defines.h"
#include "Unit1.h"
#include "matrix_utils.h"
#include "GridMap.h"
#include "bat_run.h"
#include "VCL.h"
#include "LoadData.h"
#include "zig4lib/raster.h"
#include <boost/scoped_array.hpp>
#include <cstdio>
#include <cstdlib>

// PLULA module variables
//const size_t MAX_SIZE_PLU =  5000000;
size_t glob_plu_max_cnt = 0;
size_t glob_plu_min_cnt = std::numeric_limits<size_t>::max();
const size_t MAX_NUM_PLU  = 10000000;

/// planning unit index layer
int **PLL = NULL;
/// x indices of elements ordered by planning unit index
int *PLX = NULL;
/// y indices of elements ordered by planning unit index
int *PLY = NULL;
/// number of planning units
int PLcnt = 0;

/// an array of planning units
struct PLU *PLvec; //this is being edited inside threads
int nme_cnt = 0;
int PLL_xdim = 0;
int PLL_ydim = 0;
int *PLULA_nums_tmp = NULL; //indices in plula raster
int *PLULA_cnt_tmp = NULL; //number of index pixels in plula raster
int *PLU_uplinks = NULL; 
int *PLU_up = NULL; //this is being edited inside threads
bool use_PLULA = false;

float traverse_tree_down(int node, int mode, bool remove,
			 bool nbh_loss_mode, int spnum)
{
	float sum, amount, loss;
	int   current, link_cnt;
	bool  prints = false;

	if (prints)
		Form1->Memo1->Lines->Add("NODE: "+IntToStr(PLvec[node].num));

	switch(mode)  // amount of decrease in connectivity
	{
	case 1:
		amount = 1.0f;
		break;
	case 2:
		amount = (float)PLvec[node].el_cnt;
		break;
	}

	loss = sum     = 0.0f;
	current = node;
	link_cnt=0;
	while(PLvec[current].down_link!=-1)
	{
		current = PLvec[current].down_link;
		if (prints)
			Form1->Memo1->Lines->Add("   D: "+IntToStr(PLvec[current].num));

		if (!nbh_loss_mode)
		{
			switch(mode)
			{
			case 1:
				sum++;
				break;
			case 2:
				sum += PLvec[current].el_cnt;
				break;
			}
		}
		else // hypothetical loss if removed
		{
			if (!PLvec[current].removed)
				loss += get_sp_loss_at_node_downriver(spnum, current, amount);  // amount has been calc earlier
		}

		if (remove)
			PLvec[current].tree_conn_up -= amount;  // note up - removing from upstream, global edit!
		link_cnt++;
		if (link_cnt>(2*PLcnt))
		{
			ShowMessage("EXITING: ZIG has identified circular PLU linkage in directed connectivity file, including PLU# "+IntToStr(PLvec[current].num));
			free_data();
			exit(1);
		}
	}

	if (prints)
		Form1->Memo1->Lines->Add("td-end");

	if (!nbh_loss_mode)
		return sum;
	else
		return loss;
}

//possibly edits global variables
//PLvec[current].tree_conn_up;
//PLU_up[0];
//defined inside PLULA.h
float traverse_tree_up(int node, int mode, bool remove, bool nbh_loss_mode, int spnum)
{
	float sum, amount, loss;
	int   current, link_cnt;
	bool  prints = false, found;
	int   loop, start, end, pos, tpos;

	if (prints)
		Form1->Memo1->Lines->Add("UP from node: "+IntToStr(PLvec[node].num));

	if (PLvec[node].up_link_cnt<=0)
		return 0.0f;

	switch(mode)
	{
	case 1:
		amount = 1;
		break;
	case 2:
		amount = static_cast<float>(PLvec[node].el_cnt);
		break;
	}

	loss = sum = 0.0f;
	PLU_up[0]=node; //global edit!
	pos =1;
	tpos=0;
	link_cnt=0;

	while(tpos<pos)
	{
		current = PLU_up[tpos];

		if (prints)
			Form1->Memo1->Lines->Add("   U: "+IntToStr(PLvec[current].num));

		if ((remove) && (tpos>0))
			PLvec[current].tree_conn_down -= amount;

		if (tpos>0)
		{
			if (!nbh_loss_mode)
			{
				switch(mode)
				{
				case 1:
					sum++;
					break;
				case 2:
					sum += PLvec[current].el_cnt;
					break;
				}
			}
			else // hypothetical loss if removed
			{
				if (!PLvec[current].removed)
					loss += get_sp_loss_at_node_upriver(spnum, current, amount);  // amount has been calc earlier                                                                    //
			}
		}

		start = PLvec[current].up_link_start;
		end   = PLvec[current].up_link_start+PLvec[current].up_link_cnt;
		if (start!=-1)
			for(loop=start; loop<end; loop++)
			{
				PLU_up[pos]= PLU_uplinks[loop];
				pos++;
			}
		tpos++;

		link_cnt++;
		if (link_cnt>(2*PLcnt))
		{
			ShowMessage("EXITING: ZIG has identified circular PLU linkage in directed connectivity file, including PLU# "+IntToStr(PLvec[current].num));
			free_data();
			exit(1);
		}
	}

	if (prints)
		Form1->Memo1->Lines->Add("t-up-end");

	if (!nbh_loss_mode)
		return sum;
	else
		return loss;
}

int mark_tree_up(int node, int basin_num)
{
	int   current;
	//  bool  prints = false;
	int   loop, start, end, pos, tpos, cnt;

	PLU_up[0]=node;
	pos =1;
	tpos=0;

	cnt=1;
	while(tpos<pos)
	{
		current = PLU_up[tpos];

		//     if (prints)
		//       Form1->Memo1->Lines->Add("   U: "+IntToStr(PLvec[current].num));

		PLvec[current].basin = basin_num;
		cnt++;

		start = PLvec[current].up_link_start;
		end   = PLvec[current].up_link_start+PLvec[current].up_link_cnt;
		if (start!=-1)
			for(loop=start; loop<end; loop++)
			{
				PLU_up[pos]= PLU_uplinks[loop];
				pos++;
			}
		tpos++;
	};

	return cnt;
}

void Add_to_PLULA_list(int elem)
{
	int  loop;

	for(loop=0; loop<PLcnt; loop++)
	{
		if (PLULA_nums_tmp[loop]==elem)
		{
			PLULA_cnt_tmp[loop]++;
			return;
		}
	}

	// implicitly not found, add new PLU to list
	PLULA_nums_tmp[PLcnt] = elem;
	PLULA_cnt_tmp[PLcnt]  = 1;
	PLcnt++;
}

int Get_PLULA_num(int elem)
{
	int  loop;

	for(loop=0; loop<PLcnt; loop++)
	{
		if (PLvec[loop].num==elem)
			return loop;
	}

	return -1;
}


bool  Load_and_analyse_tree_file()
{
	FILE *f;
	char line[256];
	int  num, down_num, id, down_id, load_cnt, rows, downlets, leafs,
			uplinks, pos, loop, start, end, basin_cnt, plu_in_basin, links_in_chain;

	if (PLcnt<=0)
		return false;

	f=fopen(tree_fname.toUtf8().constData(), "r+t");
	if (!f)
		return false;

	rows = 0;
	downlets = uplinks = 0;
	while(fgets(line, 256, f))
	{
		if (strlen(line)<2)
			continue;

		load_cnt=sscanf(line, "%i %i", &id, &down_id);
		if (load_cnt!=2)
		{
			Form1->Memo1->Lines->Add("Could not load two integers from row in planning units hierarchy file "+tree_fname);
			fclose(f);
			return false;
		}

		if (id!=-1)
			num      = Get_PLULA_num(id);
		else
			num = -2;  // distinguish between missing data and plula not found

		if (down_id!=-1)
			down_num = Get_PLULA_num(down_id);
		else
			down_num = -2;

		if (num==-1)
		{
			Form1->Memo1->Lines->Add("WARNING: Unknown planning unit number in directed connectivity description file. Skipping ID "+IntToStr(id));
			continue; // xxx changed to continut 4/2010
			//			fclose(f);
			//			return false;
		}

		if (down_num==-1)
		{
			Form1->Memo1->Lines->Add("WARNING: Unknown downriver planning unit number in directed connectivity description file, "+IntToStr(down_id)+" - was changed to sea connection.");
			down_num=-2;
			//          fclose(f);
			//          return false;
		}

		if (num==down_num)
		{
			Form1->Memo1->Lines->Add("WARNING: planning unit flows into itself "+IntToStr(down_id)+" - was changed to sea connection.");
			down_num=-2;
			//          fclose(f);
			//          return false;
		}


		if (down_num==-2)
		{
			downlets++;
			PLvec[num].down_link = -1;
		}
		else
		{
			PLvec[num].down_link = down_num;
			uplinks++;
			PLvec[down_num].up_link_cnt++;
		}

		rows++;
	}

	Form1->Memo1->Lines->Add("Successfully loaded 1st stage of planning units connectivity hierarchy from "+tree_fname);
	Form1->Memo1->Lines->Add("    Count of rows  "+IntToStr(rows));
	Form1->Memo1->Lines->Add("    Root PLU count "+IntToStr(downlets));
	Form1->Memo1->Lines->Add("    Uplinked PLUs  "+IntToStr(uplinks));

	pos = leafs = 0;
	for(loop=0; loop<PLcnt; loop++)
	{
		PLvec[loop].up_link_start=pos;
		if (PLvec[loop].up_link_cnt==0)
			leafs++;
		pos += PLvec[loop].up_link_cnt;
	}
	Form1->Memo1->Lines->Add("    Leaf PLU count "+IntToStr(leafs));

	rewind(f);
	while(fgets(line, 256, f))
	{
		if (strlen(line)<2)
			continue;

		load_cnt=sscanf(line, "%i %i", &id, &down_id);
		if (load_cnt!=2)
		{
			Form1->Memo1->Lines->Add("Could not load two integers from row in planning units hierarchy file "+tree_fname);
			fclose(f);
			return false;
		}

		if (id!=-1)
			num      = Get_PLULA_num(id);
		else
			num = -2;  // distinguish between missing data and plula not found

		if (id==-1)
			continue; // xxx changed to continue 4/2010

		if (down_id!=-1)
		{
			down_num = Get_PLULA_num(down_id);
			if (down_num==-1)
				down_num=-2; // "fix to sea connection"
			if (down_num==num)
				down_num=-2;
		}
		else
			down_num = -2;

		if (down_num!=-2)
		{
			start = PLvec[down_num].up_link_start;
			end   = PLvec[down_num].up_link_start+PLvec[down_num].up_link_cnt;
			for(loop=start; loop<end; loop++)
			{
				if (PLU_uplinks[loop]==-1)
				{
					PLU_uplinks[loop]=num;
					break;
				}
			}
		}
	}

	fclose(f);
	// At this stage data loaded but not linked up. xxx where to resolve circular reference?

	for(loop=0; loop<PLcnt; loop++)
	{
		PLvec[loop].tree_conn_down = PLvec[loop].orig_tree_conn_down = traverse_tree_down(loop, 2, false, false, 0);
		PLvec[loop].tree_conn_up   = PLvec[loop].orig_tree_conn_up   = traverse_tree_up(loop, 2, false, false, 0);
		if (PLvec[loop].up_link_cnt>0)
			sprintf(line, "%i up-link-cnt=%i (first=%i) dc=%f  uc=%f",PLvec[loop].num, PLvec[loop].up_link_cnt,
				PLvec[PLU_uplinks[PLvec[loop].up_link_start]].num, PLvec[loop].orig_tree_conn_down,
				PLvec[loop].orig_tree_conn_up);
		else
			sprintf(line, "%i  up-link-cnt=%i  dc=%f  uc=%f",PLvec[loop].num, PLvec[loop].up_link_cnt,
				PLvec[loop].orig_tree_conn_down, PLvec[loop].orig_tree_conn_up);
		//      Form1->Memo1->Lines->Add(line);
		//      Form1->Memo1->Lines->Add("");
	}

	basin_cnt=0;
	Form1->Memo1->Lines->Add("Count of PLUs in...");
	for(loop=0; loop<PLcnt; loop++)
	{
		if (PLvec[loop].down_link!=-1)
			continue;

		plu_in_basin=mark_tree_up(loop, basin_cnt);
		Form1->Memo1->Lines->Add("               basin #"+IntToStr(PLvec[loop].num)+" = "+IntToStr(plu_in_basin));

		basin_cnt++;
	}

	Form1->Memo1->Lines->Add("Count of basins = "+IntToStr(basin_cnt));

	return true;
}

// loads plula raster indices to PLULA_nums_tmp and PLULA_cnt_tmp
// does not load the actual raster data yet
bool load_PLULA_from_file_phase1(Raster<int>& plula, float **area_mask)
{
	PLL_xdim = plula.xsize();
	PLL_ydim = plula.ysize();
	Form1->Memo1->Lines->Add(" Checking planning units identifiers/numbers...");
	for(int y = 0; y != PLL_ydim; ++y) {
		for(int x = 0; x != PLL_xdim; ++x) {
			int itmp = plula.value(x, y);
			if(area_mask) {
				if(area_mask[y][x] <= 0)
					itmp = plula.nodatavalue();
			}
			if(isnan(itmp) || itmp == plula.nodatavalue()) {
				continue;
			}
			Add_to_PLULA_list(itmp);
			++nme_cnt;
		}
	}
	return true;
}

bool load_PLULA_from_file_phase2(Raster<int>& plula, float **area_mask)
{
	int missing = 0;
	boost::scoped_array<int> vec_pos(new int[PLcnt]);
	if(!vec_pos)
		return false;
	for(int loop = 0; loop < PLcnt; loop++)
		vec_pos[loop] = PLvec[loop].start;
	Form1->Memo1->Lines->Add(" Processing planning units layer, phase 2...");
	for(int y = 0; y != PLL_ydim; ++y) {
		for(int x = 0; x != PLL_xdim; ++x) {
			int itmp = plula.value(x, y);
			if(area_mask) {
				if(area_mask[y][x] <= 0)
					itmp = plula.nodatavalue();
			}
			if(isnan(itmp) || itmp == plula.nodatavalue()) {
				++missing;
				continue;
			}
			int pln = Get_PLULA_num(itmp);
			if(pln<0)
				return false;
			PLX[vec_pos[pln]] = x;
			PLY[vec_pos[pln]] = y;
			PLL[y][x] = pln;
			++vec_pos[pln];
			if(x > PLvec[pln].maxx)
				PLvec[pln].maxx = x;
			if(y > PLvec[pln].maxy)
				PLvec[pln].maxy = y;
			if(x < PLvec[pln].minx)
				PLvec[pln].minx = x;
			if(y < PLvec[pln].miny)
				PLvec[pln].miny = y;
		}
	}
	Form1->Memo1->Lines->Add("Missing elements in planning units layer "+IntToStr(missing));
	return true;
}

bool Alloc_PLULA_data()
{
	int loop, pos, x, y;

	PLL   = 0;
	PLX   = 0;
	PLY   = 0;
	PLvec = 0;

	PLvec = new struct PLU[PLcnt];
	if (!PLvec)
		return false;

	PLU_uplinks = new int[PLcnt];
	PLU_up      = new int[PLcnt];

	pos = 0;
	for(loop=0; loop<PLcnt; loop++)
	{
		PLvec[loop].num     = PLULA_nums_tmp[loop];
		PLvec[loop].start   = pos;
		PLvec[loop].el_cnt  = PLULA_cnt_tmp[loop];
		if (PLvec[loop].el_cnt > glob_plu_max_cnt)
		  glob_plu_max_cnt = PLvec[loop].el_cnt;
		if (PLvec[loop].el_cnt < glob_plu_min_cnt)
		  glob_plu_min_cnt = PLvec[loop].el_cnt;
		PLvec[loop].removed = false;
		PLvec[loop].at_edge = false;
		PLvec[loop].maxx    = 0;
		PLvec[loop].maxy    = 0;
		PLvec[loop].minx    = PLL_xdim;
		PLvec[loop].miny    = PLL_ydim;
		PLvec[loop].allocated = false;
		//      Form1->Memo1->Lines->Add("PLULA "+IntToStr(PLvec[loop].num)+ " start "+IntToStr(pos));
		PLvec[loop].up_link_cnt   = 0;
		PLvec[loop].up_link_start = -1;
		PLvec[loop].down_link     = -1;
		PLU_uplinks[loop]         = -1;
		pos += PLULA_cnt_tmp[loop];
	}

	if (PLULA_nums_tmp)
	{
		delete[] PLULA_nums_tmp;
		PLULA_nums_tmp = 0;
	}

	if (PLULA_cnt_tmp)
	{
		delete[] PLULA_cnt_tmp;
		PLULA_cnt_tmp = 0;
	}

	PLL = imatrix(0, PLL_ydim, 0, PLL_xdim);
	if (!PLL)
		return false;
	for(y=0; y<PLL_ydim; y++)
		for(x=0; x<PLL_xdim;x++)
			PLL[y][x]=-1;

	PLX   = new int[nme_cnt];
	PLY   = new int[nme_cnt];
	if (!PLX || !PLY)
		return false;
	for(x=0; x<nme_cnt;x++)
	{
		PLX[x]=PLY[x]=-1;
	}

	return true;
}

void Free_PLULA_data()
{
	int loop;

	if (PLvec)
		for(loop=0; loop<PLcnt; loop++)
			if (PLvec[loop].allocated)
				delete[] PLvec[loop].datavec;

	if (PLvec)
		delete[] PLvec;
	PLvec=0;

	if (PLU_uplinks)
		delete[] PLU_uplinks;
	PLU_uplinks=0;

	if (PLU_up)
		delete[] PLU_up;
	PLU_up=0;

	if (PLX)
		delete[] PLX;
	PLX=0;

	if (PLY)
		delete[] PLY;
	PLY=0;

	if (PLL)
		free_imatrix(PLL, 0, PLL_ydim, 0, PLL_xdim);
	PLL=0;

#ifdef COMPACT_VMAT
	Biodiv_Features_Occur_Container::free_plu();
#endif
}

bool check_PLULA()
{
  int err_cnt, ok_cnt, ok_data, miss_cnt, miss_data_cnt;

  Form1->Memo1->Lines->Add("Checking consistency of planning units (PLU) layer...");
  err_cnt  = 0;
  ok_cnt   = 0;
  ok_data  = 0;
  miss_cnt = 0;
  miss_data_cnt=0;
  for(size_t y=0; y<PLL_ydim; y++) {
    for(size_t x=0; x<PLL_xdim;x++) {
      //          if ((status[y][x]==-1) && (PLL[y][x]>=0))
      if (PLL[y][x]<0)
	miss_cnt++;

      if ((status[y][x]>=0) && (PLL[y][x]<0))
	err_cnt++;
      else if ((status[y][x]==-1) && (PLL[y][x]>=0))
	miss_data_cnt++;
      else {
	ok_cnt++;
	if (status[y][x]!=-1)
	  ok_data++;
      }
      
    }
  }

  if ((err_cnt+miss_data_cnt)>0) {
    Form1->Memo1->Lines->Add("WARNING: PLU check not ok - PLUs and data are not aligned:");
    Form1->Memo1->Lines->Add("         Count of locations with spp data but missing PLU info = "+IntToStr(err_cnt));
    Form1->Memo1->Lines->Add("         Count of locations with plu number but no spp data    = "+IntToStr(miss_data_cnt));
    Form1->Memo1->Lines->Add("         OK total = "+IntToStr(ok_cnt));
    Form1->Memo1->Lines->Add("         OK with data cnt = "+IntToStr(ok_data));
  } else {
    Form1->Memo1->Lines->Add("PLU check ok - PLUs and data are aligned, OK total         = "+IntToStr(ok_cnt));
    Form1->Memo1->Lines->Add("                                          OK with data cnt = "+IntToStr(ok_data));
  }
  Form1->Memo1->Lines->Add("         Missing data elements in PLU file = "+IntToStr(miss_cnt));

  if (err_cnt>0) {
    Form1->Memo1->Lines->Add("   PLU layer consistency check failed!");
    return false;
  } else {
    Form1->Memo1->Lines->Add("   PLU layer consistency check successful!");
    return true;
  }
}

void  Output_PLULA_data()
{
	char  txt[255];
	int   x, y, loop;
	FILE  *f;

	f=fopen("PLULA_converted_test.txt", "w+t");
	if (!f)
		return;
	for(y=0; y<PLL_ydim; y++)
	{
		for(x=0; x<PLL_xdim;x++)
		{
			fprintf(f, "%-2i ", PLL[y][x]);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
	fprintf(f, "\n");
	for(loop=0; loop<nme_cnt; loop++)
		fprintf(f, "%-3i  %-2i  %-2i\n", loop, PLX[loop], PLY[loop]);

	fprintf(f, "\n");
	fprintf(f, "\n");
	for(loop=0; loop<PLcnt; loop++)
		fprintf(f, "%-3i  n=%-2i  s=%-3i  c=%-3i x=%i-%i  y=%i-%i\n", loop, PLvec[loop].num,  PLvec[loop].start,
			PLvec[loop].el_cnt, PLvec[loop].minx, PLvec[loop].maxx,PLvec[loop].miny,PLvec[loop].maxy);

	fclose(f);
}

bool Load_PLULA(String fname)
{
	char txt[255];
	int  loop, sum;

	if (!use_PLULA)
		return false;

	//  ShowMessage("Loading PLULA data "+fname);

	Form1->Memo1->Lines->Add("---------");
	Form1->Memo1->Lines->Add("Loading planning units layer (PLULA) from "+fname);
	PLULA_nums_tmp = new int[MAX_NUM_PLU];
	PLULA_cnt_tmp  = new int[MAX_NUM_PLU];

	if (!PLULA_nums_tmp || !PLULA_cnt_tmp)
	{
		Form1->Memo1->Lines->Add("Out of memory analysing planning units layer.");
		return false;
	}

	for(int loop=0; loop<MAX_NUM_PLU; loop++)
	{
		PLULA_nums_tmp[loop] = -1;
		PLULA_cnt_tmp[loop]  = 0;
	}

	{ // load plula raster
		Raster<int> plula(fname.toUtf8().constData());

		// loads plula raster indices to PLULA_nums_tmp and PLULA_cnt_tmp
		// does not load the actual raster data yet
		if (!load_PLULA_from_file_phase1(plula, area_mask.m))
		{
			Form1->Memo1->Lines->Add("Failure loading planning units layer.");
			ShowMessage("Failure loading planning units layer from "+fname);
			delete[] PLULA_nums_tmp;
			delete[] PLULA_cnt_tmp;
			return false;
		}

		sprintf(txt, "Planning units layer has rows: %i,  cols: %i,  data elements: %i,  number of units: %i",
			PLL_ydim, PLL_xdim, nme_cnt, PLcnt);
		Form1->Memo1->Lines->Add(txt);

		if (!Alloc_PLULA_data())
		{
			ShowMessage("Failure allocating planning units data "+fname);
			Free_PLULA_data();
			return false;
		} else {
		  Form1->Memo1->Lines->Add("Planning unit structures allocated successfully!");
		}

		load_PLULA_from_file_phase2(plula, area_mask.m);
	}

	sum=0;
	for(loop=0; loop<=PLcnt; loop++)
		sum+=PLvec[loop].el_cnt;
	Form1->Memo1->Lines->Add("Sum of map/grid elements (cells) in planning units layer = "+IntToStr(sum));

	//  check_PLULA();

	if (false)
	  Output_PLULA_data(); // xxx only for checking

	//  Free_PLULA_data(); // xxx relocate this when calling from another prog.

	Form1->Memo1->Lines->Add("Number of planning units found: " + IntToStr(PLcnt));
	Form1->Memo1->Lines->Add("   minimum number of cells in any planning unit: " + IntToStr(glob_plu_min_cnt));
	Form1->Memo1->Lines->Add("   maximum number of cells in any planning unit: " + IntToStr(glob_plu_max_cnt));
	Form1->Memo1->Lines->Add("Planning units layer loaded and analysed successfully.");
	Form1->Memo1->Lines->Add("---------");
	Form1->Edit34->Text=IntToStr(PLcnt);

	return true;
}

