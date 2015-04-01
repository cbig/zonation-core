//---------------------------------------------------------------------------

#include "Unit1.h"
#include "GridMap.h"
#include "defines.h"
#include "matrix_utils.h"
#include "LSIdent.h"
#include "bat_run.h"
#include "PLULA.h"
#include "output.h"
#include "LoadData.h"
#include <cstring>
#include <cmath>
#include <cstdio>
#include <QFile>
#include <QTextStream>

TForm1 *Form1;

String setfn, sppfn, outgridfn, outgridfn2, outgridfn3, emffn, costfn, features_info_fn,
  curvesfn, maskfn, BQP_prof_fn, memofn, SSI_fname,
  SSI_features_info_fn, SSI_curvesfn, PLULAfn, tree_fname, ia_fname, condition_fname, retention_fname,
  groups_fname, grp_curves_fn, area_mask_fname, LEC_fname;

String admu_redist_outgridfn;
String output_transf_ig_layers_dn, output_transf_ds_layers_dn, output_transf_cst_layers_dn, 
  output_transf_ct_layers_dn, output_transf_rt_layers_dn, output_transf_mct_layers_dn, 
  output_transf_ia_layers_dn, 
  output_transf_final_layers_dn;

bool   autoclose, run_finished, terminate_flag=false, neg_weights_used=false;
//float alpha, weight;
float rem_level=0.0f, val_th, BLP, SA_z;
int   normalize, smooth, show_images, miss_as_zero,
use_cost=0, use_mask=0, use_BQP=0, add_edge_cnt=0, removal_rule,
use_SSI=0, use_tree_conn=0, use_interactions=0, mem_save_mode=0,
use_condition=0, use_retention=0, use_groups=0, mask_data=0, use_LEC=0;
bool glob_add_to_edge_borders_between_removal_mask_levels = true;

// corridors
int use_corridors = false;
String corridors_layers_fn = "";

class  GridMap obsmap[MAX_SPP_COUNT], wmap[MAX_SPP_COUNT], costmap, obs_map, maskmap, tstmap;
struct sp* spp;

struct Tgroups_info glbl_groups_info;

/// number of feature/species layers
int map_cnt;
int    IG_map_cnt, ia_cnt=0, cur_map, nonm1, m1s, removed, current;

bool forget_Rsum = true;
float  **Rmax=0, **Rsum=0, **sol=0, **sol_val=0, **curves;
char   **edge=0, **status=0;
float  Rmax_max, Rsum_max; // Rmin_max, **Rmin, **Rave, Rave_max
int    *exl = NULL, *eyl = NULL, ecnt = 0;
int **mat_edge_list_pos = NULL; // position/index of cell [y][x] in the edge list (exl, eyl, ecnt)

float CCSP = 0;
extern String output_ccsp_fn;

float z_pow(float x, float ex)
{
// some worst performance cases of powf() can be very frequent here (repr levels slightly <1.0!
// http://entropymine.com/imageworsener/slowpow/
#define USE_FAST_POW_ABF 1 
//#define USE_POW_APPROX_IEEE_754 1

#ifdef USE_FAST_POW_ABF
  if (0.25==ex)
    return sqrt(sqrt(x));
  else if (0.5==ex)
    return sqrt(x);
  else if (0.75==ex)
    return sqrt(x)*sqrt(sqrt(x));
  else if (0.125==ex) 
    return sqrt(sqrt(sqrt(x)));
  else if (1.0==ex)
    return x;
  else if (2.0==ex)
    return x*x;
  else if (3.0==ex)
    return x*x*x;
  else if (4.0==ex)
    return x*x*x*x;
  else
# ifdef USE_POW_APPROX_IEEE_754
    return powf_fast(x, ex);
# else
    return powf(x, ex);
# endif
#endif
}

static unsigned int* powf_fast_table = NULL;
void pow_fast_set_table(unsigned int* const pTable, const unsigned int precision);
float pow_fast_lookup(const float val, const float ilog2, unsigned int* const pTable, const unsigned int precision);

void
powf_fast_init(int prec)
{
  size_t powf_table_size = pow(2, powf_fast_precision) + 1;
  powf_fast_table = new unsigned int[powf_table_size];
  pow_fast_set_table(powf_fast_table, powf_fast_precision);
}

void
powf_fast_release()
{
  delete [] powf_fast_table;
}

// fast replacement for std powf()
float
powf_fast(const float x, const float y)
{
  float ilog2 = logf(x) * 1.44269504088896d; // log(r) / log(2)
  return pow_fast_lookup(y, ilog2, powf_fast_table, powf_fast_precision);
}


// Fast pow 
// Adapted from code by Harrison Ainsworth / HXA7241, Copyright (c) 2007, BSD license (new) 
const float _2p23 = 8388608.0f;

/**
 * Initialize powFast lookup table.
 *
 * @pTable     length must be 2 ^ precision
 * @precision  number of mantissa bits used, >= 0 and <= 18
 */
void
pow_fast_set_table(unsigned int* const pTable, const unsigned int precision)
{
  /* step along table elements and x-axis positions */
  float zeroToOne = 1.0f / ((float)(1 << precision) * 2.0f);
  int   i;
  for( i = 0;  i < (1 << precision);  ++i) {
    /* make y-axis value for table element */
    const float f = ((float)pow( 2.0f, zeroToOne ) - 1.0f) * _2p23;
    pTable[i] = (unsigned int)( f < _2p23 ? f : (_2p23 - 1.0f) );
    
    zeroToOne += 1.0f / (float)(1 << precision);
  }
}

/**
 * Get pow (fast!).
 *
 * @val        power to raise radix to
 * @ilog2      one over log, to required radix, of two
 * @pTable     length must be 2 ^ precision
 * @precision  number of mantissa bits used, >= 0 and <= 18
 */
float
pow_fast_lookup(const float val, const float ilog2, unsigned int* const pTable, const unsigned int precision)
{
  /* build float bits */
  const int i = (int)( (val * (_2p23 * ilog2)) + (127.0f * _2p23) );

  /* replace mantissa with lookup */
  const int it = (i & 0xFF800000) | pTable[(i & 0x7FFFFF) >>
					   (23 - precision)];
  /* convert bits to float */
  return *(const float*)( &it );
}


int sortfunc( const void *a, const void *b)
{
	struct sort_xy *s1, *s2;
	s1 = (struct sort_xy *)a;
	s2 = (struct sort_xy *)b;
	if (s1->val<s2->val)
		return -1;
	else if (s1->val>s2->val)
		return 1;
	else
		return 0;
}

int sortfunc_with_hierarchy(const void *a, const void *b)
{
  struct sort_xy *s1, *s2;
  s1 = (struct sort_xy *)a;
  s2 = (struct sort_xy *)b;

  // primary criterion: hierarchical removal level
  int lvl1 = maskmap.m[s1->y][s1->x];
  int lvl2 = maskmap.m[s2->y][s2->x];
  if (lvl1<lvl2)
    return -1;
  else if (lvl1>lvl2)
    return 1;

  // secondary criterion: cell value
  if (s1->val<s2->val)
    return -1;
  else if (s1->val>s2->val)
    return 1;
  else
    return 0;
}

int gcc_sortfunc( const void *a, const void *b)
{
	float *s1, *s2;
	s1 = (float *)a;
	s2 = (float *)b;
	if (*s1 > *s2)
		return -1;
	else if (*s1 < *s2)
		return 1;
	else
		return 0;
}

float get_cell_cnt(int s, float frac, int mode, float lvl)
{
	float *vec, sum, total_sum;
	int   x, y, pos, idx;

	//  vec  = new float[(int)(lvl*nonm1+10)]; alloc size was wrong
	vec  = new float[xd*yd];
	//  Form1->Memo1->Lines->Add("Vlen="+IntToStr((int)(lvl*nonm1+10)));

	pos  = 0;
	total_sum = 0.0f;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (status[y][x]<0)
				continue;

			if (mode==1)
			{
				if (sol[y][x]<lvl)
					continue;
			}
			else
			{
				if (sol_val[y][x]<lvl)
					continue;
			}
			vec[pos] = vmat[y][x][s];
			total_sum += vec[pos];
			++pos;
		}
	}

	qsort((void *)vec, pos, sizeof(float), gcc_sortfunc);

	sum=0.0f;
	for(idx=0; idx<pos; idx++)
	{
		sum += vec[idx];
		if (sum>(frac*total_sum))
			break;
	}

	delete[] vec;

	return (float)(idx+1);
}

void  get_solution_boundary_length(int mode, float lvl, int collect_hist, float *hvec)
{
	int   x, y, BL, AR, s, tmp_map_cnt;// sp_at_zero;
	char  txt[256];
	float min_solval, cost, cc50, cc90, local, ave_rem=0.0f, wprop_rem=0.0f, tmp_wsum;

	Screen->Cursor=crHourGlass;

	if (collect_hist!=-1)
		for(x=0; x<=20; x++)
			hvec[x]=0;

	min_solval  = 1.00f;
	BL   = 0;
	AR   = 0;
	cost = 0.0f;
	//  sp_at_zero = 0;
	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (status[y][x]<0)
				continue;

			if (mode==1)
			{
				if (sol[y][x]<lvl)
					continue;
			}
			else
			{
				if (sol_val[y][x]<lvl)
					continue;
			}

			++AR;
			if (use_cost)
				cost += costmap.m[y][x];
			else
				cost++;

			if (sol_val[y][x]<min_solval)
				min_solval = sol_val[y][x];

			if (x==0)
				BL++;
			else if (sol[y][x-1]<lvl)
				BL++;

			if (y==0)
				BL++;
			else if (sol[y-1][x]<lvl)
				BL++;

			if (x==(xd-1))
				BL++;
			else if (sol[y][x+1]<lvl)
				BL++;

			if (y==(yd-1))
				BL++;
			else if (sol[y+1][x]<lvl)
				BL++;

			for(s=0; s<map_cnt; s++)
			{
				if (vmat[y][x][s]>0)
				{
					if (spp[s].weight<0.0f)
						continue; // xxxMCZ  sol characteristics only comp for features with wj>=0
					ave_rem   += vmat[y][x][s];
					wprop_rem += fabs(spp[s].weight)*vmat[y][x][s];  // xxxMCZ check
				}
			}

			if (collect_hist!=-1)
			{
				local = vmat[y][x][collect_hist];
				//              Form1->Memo1->Lines->Add(FloatToStr(local)+" / "+FloatToStr(spp[collect_hist].max_local_val));
				if ((local>=0) && (spp[collect_hist].max_local_val>0))
					hvec[(int)(20*local/spp[collect_hist].max_local_val)]++;
			}
		}
	}

	cc50=cc90=0.0f;
	tmp_wsum=0.0f;
	tmp_map_cnt = 0;
	for(s=0; s<map_cnt; s++)
	{
		if (spp[s].weight<0.0f)
			continue;

		cc90 += get_cell_cnt(s, 0.9f, mode, lvl)/map_cnt;
		cc50 += get_cell_cnt(s, 0.5f, mode, lvl)/map_cnt;
		tmp_wsum += spp[s].weight;
		tmp_map_cnt++;
	}

	Screen->Cursor=crDefault;

	ave_rem   /= tmp_map_cnt; //xxxMCZ ok - includes only wj>0 layers?
	wprop_rem /= tmp_wsum;
	if (collect_hist==-1)
	{
		//      sprintf(txt, "  Selection data: BL=%i  Area=%i  BL/A=%0.3f  Cost=%f  min-prop-rem=%f  average_90%%_cell_count=%0.0f  50%%_cell_count=%0.0f",
		//         BL, AR, BL/(float)AR, cost, 1.0f-min_solval, cc90, cc50);
		sprintf(txt, " Data for selected top fraction: BL=%i  Area=%i  BL/A=%0.3f  Cost=%f  min-prop-rem=%f  mean-prop_rem=%f  weighted_mean_prop_rem=%f  average_90%%_cell_count=%0.0f  50%%_cell_count=%0.0f",
			BL, AR, BL/(float)AR, cost, 1.0f-min_solval, ave_rem, wprop_rem, cc90, cc50);
		Form1->Memo1->Lines->Add(txt);
		Application->ProcessMessages();
		ShowMessage(txt);
	}
}

//---------------------------------------------------------------------------
TForm1::TForm1(TComponent* Owner)
	: TForm(Owner)
{
}

void TForm1::FormClose(TObject *Sender, TCloseAction &Action)
{
  if (run_finished) {
    QFile file(memofn);
    if(file.open(QIODevice::Text | QIODevice::WriteOnly | QIODevice::Truncate)) {
      QTextStream out(&file);
      //out.setCodec("UTF-8");
      out << compileMsgList();
    } else {
      Form1->Memo1->Lines->Add("Note: error while closing the memo output text file!");
    }
  }
}
//---------------------------------------------------------------------------

#if 0
//---------------------------------------------------------------------------
void TForm1::Button3Click(TObject *Sender)
{
	Form1->Caption="Richness";
	obsmap[0].draw_on_bitmap_smooth_from_mat(Form1->Image1->Picture->Bitmap,
						 Rsum, Rsum_max, 0.0f, false, false);
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------
void TForm1::Button4Click(TObject *Sender)
{
	Form1->Caption="Rarity";
	obsmap[0].draw_on_bitmap_smooth_from_mat(Form1->Image1->Picture->Bitmap,
						 Rmax, Rmax_max, 0.0f, false, false);
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------
void TForm1::Button5Click(TObject *Sender)
{
	float prop, dummy[21];
	char txt[128];

	DecimalSeparator='.';

	val_th = StrToFloat(Form1->Edit3->Text);
	sprintf(txt,"Selection when %0.2f of the landscape remains", val_th);
	Form1->Caption=txt;

	obsmap[0].draw_on_bitmap_smooth_from_mat(Form1->Image1->Picture->Bitmap,
						 sol, 1.0f, 1.0f-val_th, true, false);

	if (use_SSI && AnnotForm->CheckBox6->Checked)
		obsmap[0].show_spots(Form1->Image1->Picture->Bitmap,  SSIxyCnt, 2);

	get_solution_boundary_length(1, 1.0f-val_th, -1, dummy);

	Form1->VertScrollBar->Position=0;
}

void TForm1::Button6Click(TObject *Sender)
{
	char txt[128];
	float dummy[21];

	DecimalSeparator='.';

	val_th = StrToFloat(Form1->Edit2->Text);
	sprintf(txt,"Area required for getting %0.2f of distributions for all features", val_th);
	Form1->Caption=txt;
	obsmap[0].draw_on_bitmap_smooth_from_mat(Form1->Image1->Picture->Bitmap,
						 sol_val, 1.0f, 1.0f-val_th, true, false);

	if (use_SSI && AnnotForm->CheckBox6->Checked)
		obsmap[0].show_spots(Form1->Image1->Picture->Bitmap,  SSIxyCnt, 2);

	get_solution_boundary_length(2, 1.0f-val_th, -1, dummy);

	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

void TForm1::ComboBox1Change(TObject *Sender)
{
	int   point, idx, loop;
	float lvl, histo[21];
	char  txt[256];

	Series1->Clear();
	LineSeries1->Clear();
	idx = ComboBox1->ItemIndex;
	Series5->Clear();

	val_th = lvl = 1.0f-StrToFloat(Edit31->Text);
	if ((lvl>0.0f) && (lvl<1.0f))
	{
		get_solution_boundary_length(2, val_th, idx, histo);
		Form1->Memo1->Lines->Add("");
		sprintf(txt, "Histogram data for species %s in %0.3f top fraction.", spp[idx].fname.toUtf8().constData(), 1.0f - lvl);
		Form1->Memo1->Lines->Add("Occurrence level is given as density relative to maximal local density.");
		Form1->Memo1->Lines->Add(txt);
		for(loop=0; loop<=20; loop++)
		{
			Series5->AddXY(0.05f*loop, histo[loop], "", clBlack);
			sprintf(txt, "%-.2f  %0.0f", 0.05f*loop, histo[loop]);
			Form1->Memo1->Lines->Add(txt);
		}
		Form1->Memo1->Lines->Add("");
		Form1->Memo1->Lines->Add("");
	}

	Series1->LinePen->Width=2;
	LineSeries1->LinePen->Width=2;
	for(point=0; point<c_pos; point++)
	{
		Series1->AddXY(1.0-curves[point][0], curves[point][idx+4],"", clBlack);
		LineSeries1->AddXY(curves[point][1], curves[point][idx+4],"", clBlack);
	}
	Form1->Memo1->Lines->Add("updated curve"+IntToStr(idx+1)+"  cpos="+IntToStr(c_pos) );
}
//---------------------------------------------------------------------------

void TForm1::Chart1DblClick(TObject *Sender)
{
	if(SaveDialog1->Execute())
	{
		String pnfilename = SaveDialog1->FileName;     //AnsiString
		Chart1->SaveToMetafileEnh(pnfilename);
	}
}
//---------------------------------------------------------------------------


void TForm1::Button7Click(TObject *Sender)
{
	DecimalSeparator='.';
	LSIdent(Form1->RadioGroup5->ItemIndex);
}
//---------------------------------------------------------------------------

void TForm1::Label13Click(TObject *Sender)
{
	AckForm->Show();
}
//---------------------------------------------------------------------------

void TForm1::ComboBox2Change(TObject *Sender)
{
	int  x,y;
	bool ok;

	if (mem_save_mode)
		return;

	cur_map    = ComboBox2->ItemIndex;

	if (CheckBox3->Checked)
	{
		if (!IG_set.use_IG)
		{
			ShowMessage("Error weight maps have not been loaded");
			return;
		}
		else
		{
			Form1->Caption="Error rates for "+spp[cur_map].fname;
			ok=Read_IGw_file(cur_map, IG_set.fnames[cur_map],0);
			if (ok)
			{
				if (IG_set.normalize_IGw)
					wmap[cur_map].average_normalize();
				wmap[cur_map].draw_on_bitmap_smooth(Form1->Image1->Picture->Bitmap,
								    wmap[cur_map].max_val, 0.0f, false);
				wmap[cur_map].free_matrix_m();
			}
		}
		return;
	}


	for(y=0;y<yd;y++)
		for(x=0;x<xd;x++)
		{
			if (vmat[y][x])
				display_mat[y][x] = vmat[y][x][cur_map];
			else
				display_mat[y][x] = -1;
		}
	obsmap[0].draw_on_bitmap_smooth_from_mat(Form1->Image1->Picture->Bitmap,
						 display_mat, obsmap[cur_map].max_val,
						 0.0f, false, false);

	if (use_SSI && AnnotForm->CheckBox6->Checked)
		obsmap[0].show_spots(Form1->Image1->Picture->Bitmap,  SSIxyCnt, 2);

	Form1->Caption=spp[cur_map].fname;
}
//---------------------------------------------------------------------------

void TForm1::UpDown1ChangingEx(TObject *Sender,
			       bool &AllowChange, short NewValue, TUpDownDirection Direction)
{
	bool show_wmap;
	int  x, y, ok;
	int  show_scale;

	if (mem_save_mode)
		return;

	if (CheckBox3->Checked && IG_set.use_IG)
		show_wmap=true;
	else
	{
		CheckBox3->Checked=false;
		show_wmap=false;
	}

	if (Direction==updUp)
	{
		if (cur_map<(map_cnt-1))
			cur_map++;
	}
	else
	{
		if (cur_map>0)
			cur_map--;
	}

	ComboBox2->ItemIndex=cur_map;

	if (!show_wmap)
	{
		//      obsmap[cur_map].draw_on_bitmap_smooth(Form1->Image1->Picture->Bitmap, 1.0f, 0.0f);
		for(y=0;y<yd;y++)
			for(x=0;x<xd;x++)
				if (vmat[y][x])
					display_mat[y][x] = vmat[y][x][cur_map];
				else
					display_mat[y][x] = -1;

		obsmap[0].draw_on_bitmap_smooth_from_mat(Form1->Image1->Picture->Bitmap,
							 display_mat, obsmap[cur_map].max_val, 0.0f, false, false);
		//  obsmap[0].draw_on_bitmap_smooth_from_mat(Form1->Image1->Picture->Bitmap,
		//                                          display_mat, obsmap[cur_map].max_val,
		//                                          0.0f, false, false);
		if (IG_set.use_IG)
			Form1->Caption="IG discounted species distribution for "+spp[cur_map].fname;
		else
			Form1->Caption="Species distribution for "+spp[cur_map].fname;
	}
	else
	{
		Form1->Caption="Error rates for "+spp[cur_map].fname;
		ok=Read_IGw_file(cur_map, IG_set.fnames[cur_map],0);
		if (ok)
		{
			if (IG_set.normalize_IGw)
				wmap[cur_map].average_normalize();
			wmap[cur_map].draw_on_bitmap_smooth(Form1->Image1->Picture->Bitmap,
							    wmap[cur_map].max_val, 0.0f, false);
			wmap[cur_map].free_matrix_m();
		}
		//      Form1->Caption="Error rates for "+spp[cur_map].fname;
		//      wmap[cur_map].draw_on_bitmap_smooth(Form1->Image1->Picture->Bitmap,
		//         wmap[cur_map].max_val, 0.0f, false);  // xxx2 fix
	}

}
//---------------------------------------------------------------------------


void TForm1::TabSheet4Show(TObject *Sender)
{
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

void TForm1::TabSheet4Enter(TObject *Sender)
{
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

void TForm1::PageControl1Enter(TObject *Sender)
{
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

void TForm1::TabSheet1Enter(TObject *Sender)
{
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

void TForm1::TabSheet2Enter(TObject *Sender)
{
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

void TForm1::TabSheet3Enter(TObject *Sender)
{
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

void TForm1::TabSheet5Enter(TObject *Sender)
{
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

void TForm1::TabSheet6Enter(TObject *Sender)
{
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

void TForm1::PageControl1Change(TObject *Sender)
{
	Form1->VertScrollBar->Position=0;
}
//---------------------------------------------------------------------------

float get_nbl(int x, int y, float **m, float lvl, int br)
{
	int sx, sy, ex, ey, cnt;

	sx = max(0, x-br);s
	ex = min(x+br, xd-br);
	sy = max(0, y-br);
	ey = min(y+br, yd-br);
	cnt= 0;

	for(int j=sy; j<=ey; j++)
		for(int i=sx; i<=ex; i++)
			if (m[j][i]>=lvl)
				++cnt;

	--cnt;
	return ( cnt/((2.0f*br+1.0f)*(2.0f*br+1.0f)-1.0f) );
}

float SPIG(float **m, float **nbm, float beta, float lvl)
{  // Biol Cons frag UCA
	int   x, y, spn;
	float spv[MAX_SPP_COUNT], mult, minv;
	Biodiv_Features_Occur_Container rowp;

	for(x=0;x<map_cnt;x++)
		spv[x]=0.0f;

	for(y=0; y<yd; y++)
		for(x=0; x<xd; x++)
		{
			if (m[y][x]<lvl)
				continue;

			mult = exp(-beta*(1.0f-nbm[y][x]));
			// rowp=&vmat[y][x][0];  // COMPACT_VMAT
			rowp = vmat[y][x];
			for(spn=0;spn<map_cnt;spn++)
			{
				spv[spn] += mult*rowp[spn];
			}
		}

	minv=10000.0f;
	for(x=0;x<map_cnt;x++)
		if (spv[x]<minv)
			minv=spv[x];

	return minv;
}

void TForm1::Button2Click(TObject *Sender)
{ // frag UCA analysis; as far as I can see, ok in ZIG2 without change
	String cfn;
	class  GridMap cmpmap;
	float  f2, beta, mp, **nbm;
	int    x, y, br;
	char   txt[128];

	DecimalSeparator='.';
	cfn = Edit21->Text;
	f2  = StrToFloat(Edit20->Text);
	br  = StrToInt(Edit26->Text);
	nbm = matrix(0,yd,0,xd);
	if (!nbm)
		return;

	cmpmap.set_no_normalize();

	if (!cmpmap.load_from_file(cfn, mask_data, area_mask.m))
	{
		ShowMessage("Could not load given comparison solution");
		return;
	}

	for(y=0; y<yd; y++)
		for(x=0; x<xd; x++)
			nbm[y][x] = get_nbl(x, y, cmpmap.m, 1.0f-f2, br);

	Memo1->Lines->Add("");
	Memo1->Lines->Add("Info-gapping given solution using presently loaded data");
	Memo1->Lines->Add("  Remaining-columns give the fraction of biological value lost");
	Memo1->Lines->Add("  at the respective IG-alpha when 87.5%, 50% and 0% of neighboring cells remain.");
	Memo1->Lines->Add("SPIG-IG-alpha  Remaining87.5%  Remaining50%   Remaining0%   min-prop-over-feature");

	beta=0.0f;
	mp = SPIG(cmpmap.m, nbm, beta, 1.0f-f2);
	sprintf(txt, "%-5.3f\t%-6.2f\t%-6.2f\t%-6.2f\t%-6.4f", beta, 1.0f-exp(-0.125f*beta), 1.0f-exp(-0.5f*beta), 1.0f-exp(-beta), mp);
	Memo1->Lines->Add(txt);

	for(beta=0.01f; beta<40; beta *= 1.5f)
	{
		mp = SPIG(cmpmap.m, nbm, beta, 1.0f-f2);
		sprintf(txt, "%-5.3f\t%-6.2f\t%-6.2f\t%-6.2f\t%-6.4f", beta, 1.0f-exp(-0.125f*beta), 1.0f-exp(-0.5f*beta), 1.0f-exp(-beta), mp);
		//      sprintf(txt, "%-5.3f  %-6.2f  %-6.2f  %-6.4f", beta, 1.0f-exp(-0.5f*beta), 1.0f-exp(-beta), mp);
		Memo1->Lines->Add(txt);
	}

	Memo1->Lines->Add("SPIG Done.");
	Memo1->Lines->Add("");
	free_matrix(nbm,0,yd,0,xd);
}
//---------------------------------------------------------------------------

void TForm1::Button9Click(TObject *Sender)
{
	terminate_flag = true;
}
//---------------------------------------------------------------------------

void TForm1::Button10Click(TObject *Sender)
{
	AnnotForm->Show();
}
//---------------------------------------------------------------------------

void TForm1::Label47Click(TObject *Sender)
{
	DisclaimerForm->Show();
}
//---------------------------------------------------------------------------

void TForm1::Label29Click(TObject *Sender)
{
	ReferenceForm->Show();
}
//---------------------------------------------------------------------------

void TForm1::Chart2DblClick(TObject *Sender)
{
	if(SaveDialog1->Execute())
	{
		String pnfilename = SaveDialog1->FileName;     //AnsiString
		Chart2->SaveToMetafileEnh(pnfilename);
	}
}
//---------------------------------------------------------------------------

void TForm1::Chart3DblClick(TObject *Sender)
{
	if(SaveDialog1->Execute())
	{
		String pnfilename = SaveDialog1->FileName;     //AnsiString
		Chart3->SaveToMetafileEnh(pnfilename);
	}
}
//---------------------------------------------------------------------------


void TForm1::Chart5DblClick(TObject *Sender)
{
	if(SaveDialog1->Execute())
	{
		String pnfilename = SaveDialog1->FileName;     //AnsiString
		Chart5->SaveToMetafileEnh(pnfilename);
	}
}
//---------------------------------------------------------------------------

void TForm1::Button11Click(TObject *Sender)
{
	SpecMap->Show();
}
//---------------------------------------------------------------------------

void TForm1::Button12Click(TObject *Sender)
{
	float fract, cval;
	char  txt[128];

	cval = StrToFloat(Edit36->Text);

	fract = get_cost_fract(cval);
	if (fract!=-1.0f)
		sprintf(txt, "Landscape top fraction corresponding to cost %f is approx. %f", cval, fract);
	else
		sprintf(txt, "Failure to identify equal cost fraction. Pls check input.");
	Form1->Memo1->Lines->Add(txt);
	ShowMessage(txt);
}
//---------------------------------------------------------------------------

void TForm1::Chart4DblClick(TObject *Sender)
{
	if(SaveDialog1->Execute())
	{
		String pnfilename = SaveDialog1->FileName;     //AnsiString
		Chart4->SaveToMetafileEnh(pnfilename);
	}
}
//---------------------------------------------------------------------------

/*
void TForm1::Button13Click(TObject *Sender)
{
 char txt[512];

 if (strlen(Edit10->Text)>1)
 {
  sprintf(txt, "notepad %s", Edit10->Text);
  // system(txt);
 }
 else
  ShowMessage("No name given for file to edit.");
}
*/
//---------------------------------------------------------------------------

/*
void TForm1::Button14Click(TObject *Sender)
{
 char txt[512];

 if (strlen(Edit16->Text)>1)
 {
  sprintf(txt, "notepad %s", Edit16->Text);
  //system(txt);
 }
 else
  ShowMessage("No name given for file to edit.");

}
*/
//---------------------------------------------------------------------------

/*
void TForm1::Button15Click(TObject *Sender)
{
 char txt[512];

 if (strlen(Edit1->Text)>1)
 {
  sprintf(txt, "notepad %s", Edit1->Text);
  //system(txt);
 }
 else
  ShowMessage("No name given for file to edit.");

}
*/
//---------------------------------------------------------------------------

/*
void TForm1::Button16Click(TObject *Sender)
{
 char txt[512];

 if (strlen(Edit30->Text)>1)
 {
  sprintf(txt, "notepad %s", Edit30->Text);
  //system(txt);
 }
 else
  ShowMessage("No name given for file to edit.");

}
*/
//---------------------------------------------------------------------------

/*
void TForm1::Button17Click(TObject *Sender)
{
 char txt[512];

 if (strlen(Edit40->Text)>1)
 {
  sprintf(txt, "notepad %s", Edit40->Text);
  //system(txt);
 }
 else
  ShowMessage("No name given for file to edit.");

}
*/
//---------------------------------------------------------------------------
/*
void TForm1::Button18Click(TObject *Sender)
{
 char txt[512];

 if (strlen(Edit9->Text)>1)
 {
  sprintf(txt, "notepad %s", Edit9->Text);
  //system(txt);
 }
 else
  ShowMessage("No name given for file to edit.");

}
*/
//---------------------------------------------------------------------------
/*
void TForm1::Button19Click(TObject *Sender)
{
 char txt[512];

 if (strlen(Edit41->Text)>1)
 {
  sprintf(txt, "notepad %s", Edit41->Text);
  //system(txt);
 }
 else
  ShowMessage("No name given for file to edit.");

}
*/
//---------------------------------------------------------------------------

void TForm1::Button20Click(TObject *Sender)
{
	FILE  *of;
	int   x,y, num;
	char  fname[512];
	float val;

	num = ComboBox1->ItemIndex;
	//sprintf(fname, "%s.#%i.exported.asc", spp[num].fname, num+1);
	of = fopen(fname, "w+t");
	if (!of)
	{
		Form1->Memo1->Lines->Add("************** ERROR ***************");
		Form1->Memo1->Lines->Add("  Error opening output file "+String(fname));
		return;
	}

	fprintf(of, "ncols\t%i\n", obsmap[num].cols());
	fprintf(of, "nrows\t%i\n", obsmap[num].rows());
	fprintf(of, "xllcorner\t%lf\n", obsmap[num].dxc);
	fprintf(of, "yllcorner\t%lf\n", obsmap[num].dyc);
	//  fprintf(of, "xllcorner\t%0.4f\n", obsmap[num].getxc());
	//  fprintf(of, "yllcorner\t%0.4f\n", obsmap[num].getyc());
	fprintf(of, "cellsize\t%0.8f\n",  obsmap[num].cell_size());
	fprintf(of, "NODATA_value\t%i\n", -1);

	for(y=0; y<obsmap[num].rows(); y++)
	{
		for(x=0; x<obsmap[num].cols(); x++)
		{
			if (!vmat[y][x])
				fprintf(of, "-1 ");
			else
			{
				val = vmat[y][x][num];

				if (val==-1)
					fprintf(of, "-1 ");
				else
					fprintf(of, "%0.5g ",val);
			}
		}
		fprintf(of, "\n");
	}

	fclose(of);
	ShowMessage("Exported data layer "+String(fname));
}
//---------------------------------------------------------------------------
#endif

