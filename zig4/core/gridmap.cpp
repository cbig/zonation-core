#include "GridMap.h"
#include "defines.h"
#include "matrix_utils.h"
#include "Unit1.h"
#include "bat_run.h"
#include "zig4lib/raster.h"
#include "VCL.h"
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <boost/scoped_array.hpp>

GridMap::GridMap()
{
	flip      = false;
	normalize = true;
	m = NULL;
	matrix_passed=0;
	//  r=g=b=0;
	//  is_suitable=0;
	//  mcols=0;
	// fedemp: I keep this initialization, but why not -1, if that's the
	// standard missing value later on?
	nodataval=-9999;
}

GridMap::~GridMap()
{
	dealloc();
}

void GridMap::dealloc()
{
  if (m && !matrix_passed) {
    free_matrix(m, 0, ydim, 0, xdim);
    m = NULL;
  }
#if 0
	if (r)
		delete[]r;
	if (g)
		delete[]g;
	if (b)
		delete[]b;
	if (is_suitable)
		delete[]is_suitable;
	if (mcols)
		delete[]mcols;
	m=0;
	r=g=b=0;
	mcols=0;
	is_suitable=0;
#endif
	loaded=false;
}

int tokenize_gridmap(char *dataitems[], char *line)
{
	int cnt, linelen, pos;

	linelen = static_cast<int>(strlen(line));

	cnt=pos=0;
	while(pos<linelen)
	{
		while ((iswspace(line[pos]) || ((char)line[pos]==(char)'\t')))
			pos++;
		if (pos<linelen)
		{
			dataitems[cnt]=&line[pos];
			cnt++;
		}

		while ( (isalnum(line[pos]) // || (line[pos]=='*') || (line[pos]=='/')
			 || (line[pos]==',') || (line[pos]=='.') || (line[pos]=='-') || (line[pos]=='+')
			 || (line[pos]=='_') || (line[pos]=='\\'))
				&& (pos<linelen) )  //jumps over alphabets, digits etc.
			pos++;
		line[pos]=0;

		pos++;
	}

	return cnt;
}
//------------------------------------------------------------------

void GridMap::pixelblock(TCanvas *c, int x, int y, TColor col, bool missing=true)
{
	c->Pixels[x][y]=col;
}

bool GridMap::loadFromFile(const char *fileName, float **tmpm, bool setGlobalProjection) {
	Raster<float> raster(fileName, 1, setGlobalProjection);

	if(!raster.rasterLoaded()) {
		Form1->Memo1->Lines->Add("Load failed: Could not load raster file "+String(fileName));
		return false;
	}

	xdim = raster.xsize();
	ydim = raster.ysize();

	// First check consistency of sizes. Too big -> crash/corruption when later trying to load into smaller buffer
	if (0!=xd && 0!=yd && 
	    (xdim !=xd || ydim !=yd)) {
	  Form1->Memo1->Lines->Add("*** ERROR in raster map dimensions: " +String(fileName)+ ". It has "+IntToStr(ydim)+ " rows x "+
				   IntToStr(xdim)+ " columns, but it should be: "+IntToStr(yd)+ " rows x "+
				   IntToStr(xd)+ " columns!");
	  return false;
	}

	xc = dxc = raster.transform(0) + ydim * raster.transform(2);
	yc = dyc = raster.transform(3) + ydim * raster.transform(5);
	cs = dcs = sqrt(raster.transform(1) * raster.transform(1) + raster.transform(2) * raster.transform(2));
	// This may be NaN, for example in ESRI EHFA .img files...
	// so it needs special handling!
	nodataval = raster.nodatavalue();
	max_val = raster.max();

	if (m != tmpm)  // for 2nd, 3rd... calls on the same gridmap object
	  dealloc();
	if (!tmpm) {
		m = matrix(0,ydim,0,xdim);
		loaded = false;
		if (!m)
		{
			Form1->Memo1->Lines->Add("Load failed: Could not allocate for map grid " + String(fileName));
			return false;
		}

	} else {
		m = tmpm;
		matrix_passed = 1;
		//      Form1->Memo1->Lines->Add("Using temporary matrix.");
	}
	if(raster.read(*m, xdim + 1) != CE_None) {
		Form1->Memo1->Lines->Add("Load failed: Could not read from raster file " + String(fileName));
		return false;
	}

	return true;
}

bool GridMap::load_from_file(const String &fname, int mask_data, float **area_mask, bool setGlobalProjection)
{
	int    row, loop, ndc=0;
	float  mv, ftmp;
	double sum;

	if(!loadFromFile(fname.toUtf8().constData(), NULL /*tmpm*/, setGlobalProjection)) {
		return false;
	}

	// reading ends, now we have m
	mv  = 0.0;
	sum = 0.0;
	for(row=0; row<ydim; row++) {
		for(loop=0;loop<xdim;loop++) {
			ftmp = m[row][loop];
			if((row==0) || (row==(ydim-1)) || (loop==0) || (loop==(xdim-1))) {
				ftmp = nodataval;
			} else if (mask_data && area_mask) {
				if (area_mask[row][loop]<=0)
					ftmp = nodataval;
			}
			// Note (NaN == NaN) is false!
			if (isnan(ftmp) || ftmp==nodataval)
			{
				ftmp=-1;
				ndc++;
				if (miss_as_zero)
					ftmp=0.0;
			}

			m[row][loop] = ftmp;
			if (ftmp!=-1)
			{
				sum += ftmp;
				//              if (ftmp<0.0)
				//                Form1->Memo1->Lines->Add("negative prob!");
			}
			mv = max(mv, ftmp);
		}
		//      if ((row%op_int)==0)
		//        Form1->Memo1->Lines->Add(IntToStr(row)+" rows read");
	}

	nodataval=-1;
	Form1->Memo1->Lines->Add(fname + ": " + IntToStr(row) + " rows, " + IntToStr(loop) + " columns read. "
				 "Nodata cells =  " + IntToStr(ndc)+ ", sum of elements = " + FloatToStr(sum));

	max_val=mv;

	if (normalize)
	{
		for(int y=0; y<ydim; y++)
		{
			for(int x=0; x<xdim; x++)
			{
				if (m[y][x]!=-1)
					m[y][x]/=sum;
			}
		}
		max_val=max_val/sum;
	}

	xl = xdim*cs;
	yl = ydim*cs;
#if 0
	if (r)
		delete[]r;
	r = new int[mv+1];
	if (g)
		delete[]g;
	g = new int[mv+1];
	if (b)
		delete[]b;
	b = new int[mv+1];
	if (mcols)
		delete[] mcols;
	mcols = new TColor[mv+1];
	if (!cols)
	{
		Form1->Memo1->Lines->Add("Load failed: could not allocate for colors in"+fname);
		return false;
	}
	if (is_suitable)
		delete[]is_suitable;
	is_suitable = new bool[mv+1];
#endif
	loaded=true;
	return true;
}

bool GridMap::load_from_file_IG(const String& fname,
				struct IG_settings& IG_set, float sp_w, float **IGw, int logit,
				float &psum, float &IG_sum, int mask_data, float **area_mask, float **tmpm, float** occur_weights_map)
{
	int    row, loop; // op_int;
	float  mv, ftmp, decr, IGa, logbt_sum;
	double val;

	if(!loadFromFile(fname.toUtf8().constData(), tmpm)) {
		return false;
	}

	non_missing_cnt = 0;
	mv  = 0.0;
	sum = 0.0;
	// Do frame outside of big loop: 
	//    ((row==0) || (row==(ydim-1)) || (loop==0) || (loop==(xdim-1)))
	for (int y=0; y<ydim; y++) {   // left, rigth
	  m[y][0] = m[y][xdim-1] = -9999;
	  if (miss_as_zero) {
	    float val = 0.0;
	    m[y][0] = m[y][xdim-1] = val;
	    non_missing_cnt +=2;
	    mv = max(mv, val);
	  }
	}
	for (int x=0; x<xdim; x++) {   // top, bottom
	  m[0][x] = m[ydim-1][x] = -9999;
	  if (miss_as_zero) {
	    float val = 0.0;
	    m[0][x] = m[ydim-1][x] = val;
	    non_missing_cnt +=2;
	    mv = max(mv, val);
	  }
	}
	for(row=1; row<ydim-1; row++)
	{
		for(loop=1;loop<xdim-1;loop++)
		{
			/* moved outside of this loop (up)
			if ((row==0) || (row==(ydim-1)) || (loop==0) || (loop==(xdim-1)))
			   ftmp = nodataval;
			else*/
			if (mask_data && area_mask && area_mask[row][loop]<=0) {
			  ftmp = nodataval;
			} else { 
			  ftmp = m[row][loop];
			}

			if (isnan(ftmp) || ftmp==nodataval)
			{
				ftmp=-9999;
				if (miss_as_zero)
					ftmp=0.0;
			} else{  
			  if ((NULL != occur_weights_map) /*&& !isnan(occur_weights_map[row][loop])*/)
			    ftmp *= occur_weights_map[row][loop];
			}

			m[row][loop] = ftmp;
			if (ftmp!=-9999) {
			  non_missing_cnt++;
			  sum += ftmp;
			}
			mv = max(mv, ftmp);
		}
	}

	nodataval=-1;
	/*
	Form1->Memo1->Lines->Add(fname + ": " + IntToStr(row) + " rows, " + IntToStr(loop)
				 + " columns read; non-missing cells: " 
				 + IntToStr(non_missing_cnt) + "; their sum: "
				 + FloatToStr(sum));
	*/
	max_val=mv;

#ifdef ZCORE_DEBUG
	for(row=0; row<ydim; row++) {
	  for(loop=0;loop<xdim;loop++) {
	    if ( isnan(m[row][loop]) ){
	      Form1->Memo1->Lines->Add("load_from_file_IG(): NaN: " + IntToStr(row) + ", " + IntToStr(row));
	      exit(1);
	    }
	  }
	}
#endif

	// xxxv3; all IG things happen between here and line 630
	// 10.12.2009; opportunity seems to operate on negative IGa on zig command line

	IG_sum = psum = sum;
	IGa=sp_w*IG_set.IGa;
	if (IGa!=0.0f) // xxx was >0.0 before opportunity
	{
		if (IGa>0.0f)
			Form1->Memo1->Lines->Add("Used Info-gap robustness(alpha>0)/discounting; alpha="+FloatToStr(IGa)+" feature weight = "+FloatToStr(sp_w));
		else
			Form1->Memo1->Lines->Add("Used Info-gap opportunity (alpha<0); alpha="+FloatToStr(IGa)+" feature weight = "+FloatToStr(sp_w));
	}
	if (IGw)
		Form1->Memo1->Lines->Add("Using Info-gap weight matrix");

	if ((sum>0) && (!logit))
	{
		decr = IGa/sum;
	}
	else if (logit)
		decr=IGa;
	else
		decr=0.0f;


	if (normalize) // xxx normalization not forced for wmaps, but yes for omaps
	{
		if (logit && (decr>0.0f)) // logit+ig
		{
			logbt_sum = 0.0f;
			for(int y=0; y<ydim; y++)
			{
				for(int x=0; x<xdim; x++)
				{
					if (m[y][x]!=-9999)
					{
						val = m[y][x];
						logbt_sum += exp(val)/(1.0f+exp(val));
					}
				}
			}
		}

		for(int y=0; y<ydim; y++)
		{
			for(int x=0; x<xdim; x++)
			{
				if (m[y][x]!=-9999)
				{
					if (!logit && sum>0)
						m[y][x]/=sum;
					if (decr==0.0f)
						continue;


					if (IGw)
					{
						if (IG_set.IG_proportional)
							m[y][x] *= (1.0f-IGw[y][x]*decr); // xxxv3 IG min removed 10.12.2009 to allow prop opportunity
						//                        m[y][x] *= (1.0f-min(1.0f,IGw[y][x]*decr));
						else
							m[y][x] -= IGw[y][x]*decr;
					}
					else
					{
						if (IG_set.IG_proportional)  // xxxv3 IG min removed 10.12.2009 to allow prop opportunity
							m[y][x] *= (1.0-decr);
						//                        m[y][x] *= (1.0-min(1.0, decr));
						else
							m[y][x] -= decr;

					}
					if (!logit)
					{
						if (m[y][x]<0.0) // occurrence level cannot become negative
							m[y][x]=0.0;
						//                      else if (m[y][x]>1.0) // xxx issue? Probability cannot exceed 1.0, but we don't really work on probs.
						//                        m[y][x]=1.0;  // commented out 10.12.2009

					}
				}
			}
		}
	}


	if (logit) // log back-transform
	{
		Form1->Memo1->Lines->Add("Backtransforming from logit space");
		for(int y=0; y<ydim; y++)
		{
			for(int x=0; x<xdim; x++)
			{
				if (m[y][x]!=-9999)
				{
					val = m[y][x];
					m[y][x] = exp(val)/(1.0f+exp(val));
				}
			}
		}
	} // xxx normalization moved out of this loop 15.8.06

	sum_orig = sum;
	sum = 0.0f;
	for(int y=0; y<ydim; y++)
		for(int x=0; x<xdim; x++)
			if (m[y][x]!=-9999)
				sum += m[y][x];

	if (IGa!=0.0f)
	{
		if (!logit)
		{
		  if (IGa<0.0)
		    Form1->Memo1->Lines->Add(fname+": fraction of original occurrences remaining after info-gap OPPORTUNITY analysis = "+FloatToStr(sum));
		  else
		    Form1->Memo1->Lines->Add(fname+": fraction of original occurrences remaining after info-gap robustness, DISTRIBUTION DISCOUNTING = "+FloatToStr(sum));
		  IG_sum=sum;
		}
		else if (logbt_sum>0.0f)
		{
			Form1->Memo1->Lines->Add(fname+": fraction of original occurrences remaining after info-gap robustness, DISTRIBUTION DISCOUNTING (info-gap) = "+FloatToStr(sum/logbt_sum));
			IG_sum=sum/logbt_sum;
		}

		if (sum<=0.0f)
		{
			Form1->Memo1->Lines->Add("Distribution discounting: alpha is too large, causing problem - all occurrences disappear for feature (e.g. species) "+fname);
			return false;
		}
	}

	// if sum==0, all of them should already be 0s
	if (normalize && sum>0)  // the IG error matrix is not normalized
		for(int y=0; y<ydim; y++)
			for(int x=0; x<xdim; x++)
			  if (m[y][x]!=-9999)
			    m[y][x]/=sum;


#ifdef ZCORE_DEBUG
	for(row=0; row<ydim; row++) {
	  for(loop=0;loop<xdim;loop++) {
	    if ( isnan(m[row][loop]) ){
	      Form1->Memo1->Lines->Add("At the end of IG: NaN: " + IntToStr(row) + ", " + IntToStr(row));
	      exit(1);
	    }
	  }
	}
#endif

	for(int y=0; y<ydim; y++)
	{
		for(int x=0; x<xdim; x++)
		{
			if (m[y][x]==-9999)
				m[y][x] = -1;
		}
	}

#if 0
	max_val=max_val/sum;
	max_val -= decr;
	// explicitly check NaN
	if (isnan(max_val) || max_val<0)
		max_val=0.0f;
	else if (max_val>1.0f)
		max_val=1.0f;
#endif

	max_val=0.0f;
	if (sum>0) {
	  for(int y=0; y<ydim; y++) {
	    for(int x=0; x<xdim; x++)
	      if (m[y][x]>max_val)
		max_val = m[y][x];
	  }
	}
	xl = xdim*cs;
	yl = ydim*cs;

	loaded=true;
	return true;
}

// just load values, don't care about normalizations, etc.
bool GridMap::load_from_file_IG_just_to_count(const String& fname,
					      int mask_data, float **area_mask, float **tmpm,
					      int** occur_count_matrix)
{
  if(!loadFromFile(fname.toUtf8().constData(), tmpm)) {
    return false;
  }

  float mv = 0, ftmp;
  double sum = 0; 
  // Just ignore/skip frame: left, bottom, right, up
  //(row==0) || (row==(ydim-1)) || (loop==0) || (loop==(xdim-1))
  bool mask_yeah = mask_data && area_mask;
  int row, loop;
  non_missing_cnt = 0;
  for(row=1; row<ydim-1; row++) {
    for(loop=1; loop<xdim-1; loop++) {
      ftmp = m[row][loop];
      if ((nodataval==ftmp || isnan(ftmp))
	  ||
	  (mask_yeah && area_mask[row][loop]<=0)
	  ) {
	continue;
      }

      non_missing_cnt++;
      // m[row][loop] = ftmp;
      //sum += ftmp; 
      // mv = max(mv, ftmp); // save this
      // this is what really matters
      occur_count_matrix[row][loop]++;
    }
  }

  nodataval = -1;
  //max_val = mv;
  loaded = true;
  return true;
}
void GridMap::get_dimensions(float &xcorner, float &ycorner,
			     float &cell_size, int &xd, int &yd)
{
	xcorner  = xc;
	ycorner  = yc;
	cell_size= cs;
	xd       = xdim;
	yd       = ydim;
}

#if 0
float GridMap::get_value_at_point(float x, float y)
{
	int lx, ly;

	lx = (int)((x-xc)/cs);
	ly = (int)((y-yc)/cs);

	if ((lx>=xdim) || (ly>=ydim) || (lx<0) || (ly<0))
		return -1;
	return m[ly][lx];
}

float GridMap::get_value_at_matrix_point(int x, int y, float &cx, float &cy)
{
	cx = xc+x*cs;
	cy = yc+y*cs;

	return m[y][x];
}

TColor GridMap::get_color_for_val(float val)
{
	float lv;
	int intens;

	if (val<0)
		return clWhite;
	else if (val>1.0f)
		return clGreen;

	if (val>0)
	{
		lv = log10(val);
		intens = (int)(-30.0f*lv);
	}
	else
	{
		intens = 0;
	}

	if (intens>255)
		intens=255;
	return (TColor)RGB(intens, intens, intens);
}

TColor GridMap::get_color_at_point(float x, float y)
{
	int   lx, ly;
	float val;

	lx = (int)((x-xc)/cs);
	ly = (int)((y-yc)/cs);

	if ((lx>=xdim) || (ly>=ydim))
		val=-1;
	else
		val= m[ly][lx];

	if ((val>=0)  && (val<=max_val))
		return get_color_for_val(val);  // return mcols[val];
	else
		return clWhite;
}

void GridMap::get_matrix_coord_at_point(float x, float y, int &lx, int &ly)
{
	lx = (int)((x-xc)/cs);
	ly = (int)((y-yc)/cs);
}

int GridMap::get_cell_center_etc_at_matrix(int X, int Y, float &cx, float &cy)
{
	float ht;

	cx = xc+(X+0.5f)*cs;
	cy = yc+(Y+0.5f)*cs;

	if ((Y<0) || (X<0) ||(Y>ydim) || (X>xdim))
		return -1;

	ht = m[Y][X];
	if (ht>0)
		return 1;
	else
		return 0;

	//  if (is_suitable[ht])
	//    return 1;
	//  else
	//    return 0;
}

bool GridMap::set_color_for_type(float val, int R, int G, int B, bool suitable)
{
	int col;

	if (!loaded)
		return false;
	if ((val<0) || (val>max_val))
	{
		Form1->Memo1->Lines->Add("Color error");
	}

	r[val]=R;
	g[val]=G;
	b[val]=B;
	is_suitable[val]=suitable;
	col = B*0x10000+G*0x100+R;
	return mcols[val]=(TColor)(col);

}
#endif

bool GridMap::align_with_image(int xd, int yd, int &new_xd, int &new_yd)
{
	char tmp[111];

	if (!loaded)
		return false;

	rr = min(xd/xl,yd/yl);

	if (rr>0)
		rp=1/rr;
	else
		return false;

	ixdim=xd;
	iydim=yd;
	//sprintf(tmp,"%f %f %f %i %i", xl, yl, rr, xd, yd); // CHECK THIS, xl, yl are rubbish
	//Form1->Memo1->Lines->Add(tmp);
	return true;
}

#if 0
void GridMap::get_im_coordinates(float x, float y, int &new_xd, int &new_yd)
{
	new_xd = (int)((x-xc)*rr);
	if (flip)
		new_yd = (int)(iydim-(y-yc)*rr);
	else
		new_yd = (int)((y-yc)*rr);
}

float GridMap::get_value_at_image_point(int x, int y)
{
	float vx, vy;
	vx = xc + rp*x;
	vy = yc + rp*y;

	return   get_value_at_point(vx,vy);
}

TColor GridMap::get_color_at_image_point(int x, int y)
{
	float vx, vy;
	vx = xc + rp*x;
	vy = yc + rp*y;

	return   get_color_at_point(vx,vy);
}

bool GridMap::draw_on_bitmap(Graphics::TBitmap *bm)
{
	int loop, loop2;
	char tmp[200];
	TCanvas *c;

	c = bm->Canvas;

	if ((bm->Width!=ixdim) || (bm->Height!=iydim))
	{
		sprintf(tmp,"%i %i %i %i", ixdim, iydim, bm->Width, bm->Height);
		Form1->Memo1->Lines->Add(tmp);
		return false;
	}


	Screen->Cursor=crHourGlass;
	for(loop=0;loop<ydim;loop++)
	{
		for(loop2=0;loop2<xdim;loop2++)
		{

			//          c->Pixels[loop][loop2]=(TColor)(loop*loop2);
			//          if (flip)
			//            c->Pixels[loop2][iydim-1-loop]=get_color_at_image_point(loop2, loop);
			//          else
			//            c->Pixels[loop2][loop]=get_color_at_image_point(loop2, loop);
			if (flip)
				pixelblock(c, loop2, ydim-1-loop, get_color_at_image_point(loop2, loop)); // xxx prblem
			else
				pixelblock(c, loop2, loop, get_color_at_image_point(loop2, loop)); // xxx prblem
		}
	}
	Screen->Cursor=crDefault;

	return true;
}
#endif

void GridMap::Annotate(Graphics::TBitmap *bm)
{
	int    x, y, s;
	String txt;

	bm->Canvas->Font->Style = TFontStyles()<<fsBold;
	bm->Canvas->Pen->Color   = clWhite;
	bm->Canvas->Brush->Color = clWhite;
	bm->Canvas->Font->Color  = clBlack;

	if (AnnotForm->CheckBox1->Checked)
	{
		x   = StrToInt(AnnotForm->Edit2->Text);
		y   = StrToInt(AnnotForm->Edit3->Text);
		s   = StrToInt(AnnotForm->Edit4->Text);
		txt = AnnotForm->Edit1->Text;
		bm->Canvas->Font->Name="Arial";
		bm->Canvas->Font->Size  = s;
		bm->Canvas->TextOutA(x,y,txt);
	}

	if (AnnotForm->CheckBox2->Checked)
	{
		x   = StrToInt(AnnotForm->Edit6->Text);
		y   = StrToInt(AnnotForm->Edit7->Text);
		s   = StrToInt(AnnotForm->Edit8->Text);
		txt = AnnotForm->Edit5->Text;
		bm->Canvas->Font->Name="Arial";
		bm->Canvas->Font->Size  = s;
		bm->Canvas->TextOutA(x,y,txt);
	}

	if (AnnotForm->CheckBox3->Checked)
	{
		x   = StrToInt(AnnotForm->Edit9->Text);
		y   = StrToInt(AnnotForm->Edit10->Text);
		s   = StrToInt(AnnotForm->Edit11->Text);
		bm->Canvas->Pen->Width = 3;
		bm->Canvas->Pen->Color = clBlack;
		bm->Canvas->MoveTo(x,y);
		bm->Canvas->LineTo(x+s,y);
		//      bm->Canvas->Font->Style = TFontStyles()<<fsBold;
		//      bm->Canvas->TextOutA(x,y,txt);
	}
}

bool GridMap::draw_on_bitmap_smooth(Graphics::TBitmap *bm,
				    float maxv, float cut_level, bool rank_mode)
{
	int loop, loop2, intens;
	char tmp[200];
	TCanvas *c;
	TColor  col;
	bool    print_cols, show_outline, missing;
	float   val;

	c = bm->Canvas;

	if ((bm->Width!=ixdim) || (bm->Height!=iydim))
	{
		sprintf(tmp,"%i %i %i %i", ixdim, iydim, bm->Width, bm->Height);
		Form1->Memo1->Lines->Add(tmp);
		return false;
	}

	show_outline=AnnotForm->CheckBox4->Checked;
	if (AnnotForm->RadioGroup1->ItemIndex==0)
		print_cols=false;  // = colours!!
	else
		print_cols=true;   // B&W!

	maxv=max_val;
	if (maxv<=0.0f)
		maxv=1.0f;

	Screen->Cursor=crHourGlass;
	for(loop=0;loop<ydim;loop++)
	{
		for(loop2=0;loop2<xdim;loop2++)
		{
			missing=false;
			intens = static_cast<int>(255*m[loop][loop2]/maxv);
			//          val = m[loop][loop2];
			if (print_cols) // B&W
			{
				col =(TColor)RGB(255-intens, 255-intens, 255-intens); // xxxcolfix

				if (intens <= cut_level*255)
					col =(TColor)RGB(0, 0, 0); // xxxcolfix

				if (m[loop][loop2]<0)
				{
					if (show_outline && Is_on_border(m, loop, loop2))
						col =(TColor)RGB(0, 0, 0);
					else
						col =(TColor)RGB(255, 255, 255);
					missing=true;
				}
			}
			else
			{
				if (!rank_mode)
				{
					col =(TColor)RGB(255-intens, 255-intens, 255-intens);

					if (m[loop][loop2]<0.0f)
					{
						col =(TColor)RGB(0, 0, 100);
						missing=true;
					}
					else if (m[loop][loop2]==0.0f)
						col =(TColor)RGB(0, 0, 255);
				}
				else
				{
					if (m[loop][loop2]<0.0f)
					{
						col =(TColor)RGB(255, 255, 255);
						missing=true;
					}
					else if (m[loop][loop2]<0.2f)
						col =(TColor)RGB(0, 0, 0);
					else if (m[loop][loop2]<0.5f)
						col =(TColor)RGB(0, 0, 100);
					else if (m[loop][loop2]<0.75f)
						col =(TColor)RGB(0, 0, 255);
					else if (m[loop][loop2]<0.9f)
						col =clYellow;
					else if (m[loop][loop2]<0.95f)
						col =(TColor)RGB(255, 0, 255);
					else if (m[loop][loop2]<0.98f)
						col =(TColor)RGB(150, 0, 0);
					else
						col =(TColor)RGB(255, 0, 0);

					if ((m[loop][loop2]<cut_level) && (m[loop][loop2]>0))
						col=(TColor)RGB(0, 0, 0);
				}
			}

			if (flip)
				pixelblock(c, loop2, ydim-1-loop, col, missing);
			else
				pixelblock(c, loop2, loop, col, missing);
		}
	}

	Annotate(bm);
	Screen->Cursor=crDefault;

	return true;
}

void GridMap::set_color_for_pixel(Graphics::TBitmap *bm,
				  int x, int y, float maxv, float val, bool rank_mode)
{  int intens;
	TCanvas *c;
	TColor col;
	bool missing;

	c = bm->Canvas;
	if (maxv<=0.0f)
		maxv=1.0f;
	intens = static_cast<int>(255*val/maxv);
	col =(TColor)RGB(intens, intens, intens);
	missing=false;

	if (!rank_mode)
	{
		if (intens <0)
		{
			col =(TColor)RGB(0, 0, 150);
			missing=true;
		}
		else if (intens==0)
			col =(TColor)RGB(0, 0, 250);
	}
	else
	{
		if (val<0.0f)
		{
			col =(TColor)RGB(255, 255, 255);
			missing=true;
		}
		else if (val<0.2f)
			col =(TColor)RGB(0, 0, 0);
		else if (val<0.5f)
			col =(TColor)RGB(0, 0, 100);
		else if (val<0.75f)
			col =(TColor)RGB(0, 0, 255);
		else if (val<0.9f)
			col =clYellow;
		else if (val<0.95f)
			col =(TColor)RGB(255, 0, 255);
		else if (val<0.98f)
			col =(TColor)RGB(150, 0, 0);
		else
			col =(TColor)RGB(255, 0, 0);

	}
	c->Pixels.Image.removeSite(x, y, col, val);
	/*
	  if (flip)
	  pixelblock(c, x, ydim-1-y, col, missing);
	  else
	  pixelblock(c, x, y, col, missing);
	*/
	//  if (flip)
	//     c->Pixels[x][ydim-1-y]=col;
	//  else
	//     c->Pixels[x][y]=col;
}

bool GridMap::Is_on_border(float **mm, int y, int x)
{
	int  xmin, ymin, xmax, ymax, w=1;
	int  loop, loop2;
	bool dfound;

	if (mm[y][x]>=0)
		return false;

	xmin = max(x-w,0);
	ymin = max(y-w,0);
	xmax = min(x+w,xdim-1);
	ymax = min(y+w,ydim-1);
	for(loop2=ymin; loop2<=ymax; loop2++)
		for(loop=xmin; loop<=xmax; loop++)
			if (mm[loop2][loop]>=0)
				return true;

	return false;
}

bool GridMap::draw_on_bitmap_smooth_from_mat(Graphics::TBitmap *bm, float **mm,
					     float maxv, float cut_level,
					     bool rank_mode, int green_mode)
{
	int loop, loop2, intens;
	char tmp[200];
	TCanvas *c;
	TColor  col;
	bool    print_cols, show_outline, missing;
	float   val;

	c = bm->Canvas;
	if (maxv<=0.0f)
		maxv=1.0f;

	//  Form1->Memo1->Lines->Add("SFM");
	if ((bm->Width!=ixdim) || (bm->Height!=iydim))
	{
		sprintf(tmp,"Error in layer / display map dimensions %i %i / %i %i", ixdim, iydim, bm->Width, bm->Height);
		Form1->Memo1->Lines->Add(tmp);
		return false;
	}

	show_outline=AnnotForm->CheckBox4->Checked;
	if (AnnotForm->RadioGroup1->ItemIndex==0)
		print_cols=false; // colour
	else
		print_cols=true;  // B&W
	Screen->Cursor=crHourGlass;

	for(loop=0;loop<ydim;loop++)
	{
		for(loop2=0;loop2<xdim;loop2++)
		{
			missing=false;
			if (!print_cols) // color
			{
				intens = static_cast<int>(255*mm[loop][loop2]/maxv);
				//              val    = mm[loop][loop2];
				if (!rank_mode)
				{
					if (!green_mode)
						col =(TColor)RGB(255-intens, 255-intens, 255-intens);  //xxxv2 changed from straight intens
					else
					{
						intens = static_cast<int>(220*mm[loop][loop2]/maxv);
						col = (TColor)RGB(0, 220-intens, 0);
					}

					if (mm[loop][loop2]<0.0f)
					{
						col =(TColor)RGB(0, 0, 100);
						missing=true;
					}
					else if (mm[loop][loop2]==0.0f)
						col =(TColor)RGB(0, 0, 255);
					else if (intens<(255*cut_level))
						col =(TColor)RGB(0, 255, 0);
				}
				else
				{
					if (green_mode)
					{
						intens = static_cast<int>(220*mm[loop][loop2]/maxv);
						if (mm[loop][loop2]<0)
						{
							col =(TColor)RGB(255, 255, 255);
							missing=true;
						}
						else
							col = (TColor)RGB(0, 220-intens, 0);
					}
					else
					{
						if (mm[loop][loop2]<0.0f)
						{
							col =(TColor)RGB(255, 255, 255);
							missing=true;
						}
						else if (mm[loop][loop2]<0.2f)
							col =(TColor)RGB(0, 0, 0);
						else if (mm[loop][loop2]<0.5f)
							col =(TColor)RGB(0, 0, 100);
						else if (mm[loop][loop2]<0.75f)
							col =(TColor)RGB(0, 0, 255);
						else if (mm[loop][loop2]<0.9f)
							col =clYellow;
						else if (mm[loop][loop2]<0.95f)
							col =(TColor)RGB(255, 0, 255);
						else if (mm[loop][loop2]<0.98f)
							col =(TColor)RGB(150, 0, 0);
						else
							col =(TColor)RGB(255, 0, 0);

						if ((mm[loop][loop2]<cut_level) && (mm[loop][loop2]>0))
							col=(TColor)RGB(0, 0, 0);
					}
				}
			}
			else // yes print_cols = B&W!
			{
				intens = static_cast<int>(255*mm[loop][loop2]/maxv);
				col =(TColor)RGB(255-intens, 255-intens, 255-intens);
				//              intens = 210*mm[loop][loop2]/maxv;
				//              col =(TColor)RGB(211-intens, 211-intens, 211-intens);

				if (intens <0)
				{
					col =(TColor)RGB(255, 255, 255);
					if (show_outline && Is_on_border(mm, loop, loop2))
						col =(TColor)RGB(0, 0, 0);
					missing=true;
				}
				else if (intens==0)
					col =(TColor)RGB(255, 255, 255);
				else if (intens<(255*cut_level))
					col =(TColor)RGB(225, 225, 225);
			}

			if (flip)
				pixelblock(c, loop2, ydim-1-loop, col, missing);
			else
				pixelblock(c, loop2, loop, col, missing);
#if 0
			if (flip)
				c->Pixels[loop2][ydim-1-loop]=col;
			else
				c->Pixels[loop2][loop]=col;
#endif
			//           c->Brush->Color=col;
			//           c->FillRect(TRect(loop2,loop, loop2+1,loop+1));
		}
	}
#if 0
	if (AnnotForm->CheckBox5->Checked) // thicken lines
	{
		for(loop=0;loop<(ydim-1);loop++)
		{
			for(loop2=0;loop2<(xdim-1);loop2++)
			{
				ok =false;
				col = c->Pixels[loop][loop2];
				if (!print_cols)
				{
					if (col==(TColor)RGB(0, 0, 100))
						ok =true;
				}
				else
				{
					if (col==(TColor)RGB(255, 255, 255))
						ok =true;
				}
				if (ok)
					c->Pixels[loop2][loop]=c->Pixels[loop2][loop+1];
				ok =false;
				col = c->Pixels[loop][loop2];
				if (!print_cols)
				{
					if (col==(TColor)RGB(0, 0, 100))
						ok =true;
				}
				else
				{
					if (col==(TColor)RGB(255, 255, 255))
						ok =true;
				}
				if (ok)
					c->Pixels[loop2][loop]=c->Pixels[loop2+1][loop];
			}
		}
	}
#endif
	Annotate(bm);
	Screen->Cursor=crDefault;

	return true;
}

#if 0
bool GridMap::add_to_bitmap(Graphics::TBitmap *bm)
{
	int loop, loop2;
	char tmp[200];
	TCanvas *c;
	int val, val2;
	TColor col;

	c = bm->Canvas;

	if ((bm->Width!=ixdim) || (bm->Height!=iydim))
	{
		sprintf(tmp,"%i %i %i %i", ixdim, iydim, bm->Width, bm->Height);
		Form1->Memo1->Lines->Add(tmp);
		return false;
	}


	Screen->Cursor=crHourGlass;
	for(loop=0;loop<ydim;loop++)
	{
		for(loop2=0;loop2<xdim;loop2++)
		{

			//          c->Pixels[loop][loop2]=(TColor)(loop*loop2);
			if (flip) // xxx
			{
				val = (int)c->Pixels[loop2][(ydim-1-loop)];
				val2 = (int)get_color_at_image_point(loop2, loop);
				pixelblock(c, loop2, loop, (TColor)((val+val2)/2) );
			}
			else
			{
				val = (int)c->Pixels[loop2][loop];
				val2 = (int)get_color_at_image_point(loop2, loop);
				pixelblock(c, loop2, loop, (TColor)((val+val2)/2) );
				//              c->Pixels[loop2][loop]=(TColor)((val+val2)/2);
			}
		}
	}
	Screen->Cursor=crDefault;

	return true;
}

float GridMap::get_suitable_area(float x, float y, float a)
{
	float L, prop, cx, cy, ok, nok;
	int   ML, mx, my, X, Y;
	int   hab_ok;
	char  txt[255];

	if (a<(cs*cs))
		return a;

	L  = sqrt(a)/2.0f;
	ML = (int)ceil(L/cs);

	get_matrix_coord_at_point(x, y, mx, my);

	ok=nok=0.0f;
	for(X=mx-ML;X<=mx+ML;X++)
	{
		for(Y=my-ML;Y<=my+ML;Y++)
		{
			hab_ok=get_cell_center_etc_at_matrix(X, Y, cx, cy);
			if ((cx>(x-L)) && (cx<(x+L)) &&
					(cy>(y-L)) && (cy<(y+L)))
			{
				if (hab_ok==1)
					ok++;
				else if (hab_ok==0)
					nok++;
			}
		}
	}
	//  sprintf(txt, "%f %0.0f %0.0f %0.5f (%f %f %f %f %i %i %i)",a,ok,nok, a*ok/(ok+nok),
	//    x,y,cx,cy,mx,my,ML);
	//  if (nok>0.0f)
	//    Form1->Memo1->Lines->Add(txt);
	return a*ok/(ok+nok);
}
#endif

TColor col_list[8]={clRed, clYellow,clFuchsia,
		    clPurple,clTeal,clMaroon, clOlive, clLime};

bool GridMap::show_spots(Graphics::TBitmap *bm, int **cm, int comp_mode)
{
	int loop, loop2, cmax=0, intens;
	TCanvas *c;
	TColor  *cols, col;
	//  bool    missing;

	c = bm->Canvas;
	//  missing = true;
	//  add = time(NULL)%8;
	for(loop=0;loop<ydim;loop++)
		for(loop2=0;loop2<xdim;loop2++)
		{
			cmax=max(cmax, cm[loop][loop2]);
		}

	cols = new TColor[cmax+1];
	if (comp_mode==0)
	{
		for(loop=0; loop<=cmax;loop++)
		{
			if (AnnotForm->RadioGroup1->ItemIndex==0)
				cols[loop] = (TColor)(RGB(random(255),random(255),random(255)));
			else
			{
				intens     = random(12)*20;
				cols[loop] = (TColor)(RGB(intens, intens, intens));
			}
		}
	}
	else if (comp_mode==1)
	{
		if (AnnotForm->RadioGroup1->ItemIndex==0)
		{
			cols[1] = clYellow;
			cols[2] = clLime;
			cols[3] = clGreen;
		}
		else
		{
			cols[1] = (TColor)RGB(0,0,0);
			cols[2] = (TColor)RGB(70, 80, 80);
			cols[3] = (TColor)RGB(140, 140, 160);
		}
	}
	else  // show_SSI
	{
		for(loop=0; loop<=cmax;loop++)
			cols[loop] = (TColor)(RGB(0,255,0));
	}

	for(loop=0;loop<ydim;loop++)
	{
		for(loop2=0;loop2<xdim;loop2++)
			if (cm[loop][loop2]>0)
			{
				//            col = col_list[(cm[loop][loop2]+add)%8];
				col = cols[cm[loop][loop2]]; //(TColor)(RGB(random(255),random(255),random(255)));
				pixelblock(c, loop2, loop, col);
				if (comp_mode==2) // the SSI
				{
					pixelblock(c, loop2, loop-1, col);
					pixelblock(c, loop2-1, loop, col);
					pixelblock(c, loop2-1, loop-1, col);
				}
			}
			else if (cm[loop][loop2]==-2)
			{
				if (AnnotForm->RadioGroup1->ItemIndex==1)
					col = (TColor)RGB(255, 255, 255);
				else
					col = clBlue;
				pixelblock(c, loop2, loop, col, true);
			}
			else if ((cm[loop][loop2]==0) && (comp_mode==1))
			{
				if (AnnotForm->RadioGroup1->ItemIndex==1)
					col = (TColor)RGB(210, 210, 210);
				else
					col = clBlue;
				pixelblock(c, loop2, loop, col);
				//            c->Pixels[loop2][loop]=col;
			}
	}

	Annotate(bm);

	if (cols)
		delete[] cols;

	return true;
}

// essentially the same as output.cpp -> void output_grid(int num)
void GridMap::export_GIS_INT(int **nwm, const String& fname)
{
	String out(fname);
	int xsize = cols(), ysize = rows();
	float nodatavalue(-1.0f);
	boost::scoped_array<int> plane(new int[xsize * ysize]);
	for(int y = 0; y < ysize; ++y) {
		for(int x = 0; x < xsize; ++x) {
			int i(x + y * xsize);
			plane[i] = nwm[y][x];
		}
	}
	SaveToRaster<int>(ChangeFileExt(out), plane.get(), nodatavalue, xsize, ysize);
	/*
 FILE *of;
 int x,y;

 if (!nwm)
  return;

 of = fopen(fname.toUtf8().constData(), "w+t");
 if (!of)
    {
  Form1->Memo1->Lines->Add("Error opening output file "+fname);
  return;
    }

 fprintf(of, "ncols\t%i\n", cols());
 fprintf(of, "nrows\t%i\n", rows());
 fprintf(of, "xllcorner\t%0.4f\n", getxc());
 fprintf(of, "yllcorner\t%0.4f\n", getyc());
 fprintf(of, "cellsize\t%0.8f\n",  cell_size());
 fprintf(of, "NODATA_value\t%i\n", -1);

 for(y=0; y<rows(); y++)
    {
  for(x=0; x<cols(); x++)
	{
   fprintf(of, "%i ",nwm[y][x]);
	}
  fprintf(of, "\n");
    }

 fclose(of);
 */
}

void GridMap::average_normalize()
{
	float ave_val, ave_cnt;

	ave_val=0.0f;
	ave_cnt=0;

	for(int y=0; y<ydim; y++)
	{
		for(int x=0; x<xdim; x++)
		{
			if (m[y][x]>=0)
			{
				ave_cnt++;
				ave_val += m[y][x];
			}
		}
	}

	ave_val /= ave_cnt;
	if (!ave_val)
		return;

	for(int y=0; y<ydim; y++)
	{
		for(int x=0; x<xdim; x++)
		{
			if (m[y][x]>=0)
			{
				m[y][x] /= ave_val;
			}
		}
	}

	max_val=-std::numeric_limits<float>::max();
	for(int y=0; y<ydim; y++)
	{
		for(int x=0; x<xdim; x++)
		{
			if (m[y][x]>max_val)
				max_val=m[y][x];
		}
	}
}

void GridMap::free_matrix_m()
{
  if (!matrix_passed) {
    free_matrix(m,0,ydim,0,xdim);
  }
  m = NULL;
}

void GridMap::force_free_matrix_m()
{
  if (m) {
    free_matrix(m,0,ydim,0,xdim);
    m = NULL;
  }
}
