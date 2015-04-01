#include <fstream>
#include <sys/stat.h>
#include "mem_dump.h"
#include "Unit1.h"
#include "bat_run.h"  // vmat declared here
#include <ctime> // for time()
#include <zlib.h>

// Raw and compressed dump of vmat
// Depends on globals: vmat, outgridfn

const char* magic = "Z vmat dump file";

const String comp_suffix = ".gz";
const size_t gz_buf_size = 131072; // 128 KB
const size_t REPORT_PERIOD = 100;

bool
write_vmat_raw(std::ofstream& f, size_t xdim, size_t ydim, int map_cnt);
bool
write_vmat_compressed(const String& fn, size_t xdim, size_t ydim, int map_cnt);

bool
load_vmat_raw(std::ifstream& f, size_t xdim, size_t ydim, int map_cnt);
bool
load_vmat_compressed(const String& fn, size_t xdim, size_t ydim, int map_cnt);

String 
make_vmat_output_name(const String& proj_name, const String& fn, bool compress)
{
  // Build full path
  String vmat_fp = ChangeFileExt(ChangeFileExt(outgridfn)) + "." + fn;
  if (compress)
    vmat_fp += comp_suffix;
  return vmat_fp;
}

bool
save_vmat(const String& proj_name, const String& fn, size_t xdim, size_t ydim, int map_cnt, bool compress)
{
  // String vmat_fp = /*outdir +*/ proj_name + "." + fn;
  // if (compress)
  //   vmat_fp += comp_suffix;
  String vmat_fp = make_vmat_output_name(outgridfn, fn, compress);
  Form1->Memo1->Lines->Add("Saving vmat into: " + vmat_fp);
  std::ofstream f(vmat_fp.toUtf8().constData(), std::ios::out | std::ios::binary);
  if (!f) {
    Form1->Memo1->Lines->Add(" FATAL ERROR: failed to open vmat dump file for writing" + vmat_fp);
    return false;
  }

  Form1->Memo1->Lines->Add("Dimensions: columns: "+IntToStr(xdim)+
			   ", rows: "+IntToStr(ydim) +
			   ", features: "+IntToStr(map_cnt));
  bool r = false;
  time_t ti = time(NULL);
  if (!compress) {
    r = write_vmat_raw(f, xdim, ydim, map_cnt);
    struct stat results;
    if (r && 0 == stat(vmat_fp.toUtf8().constData(), &results))
      Form1->Memo1->Lines->Add("Finished writing vmat dump file, size: " + 
			       FloatToStrF((float)results.st_size/(1024.0*1024.0), ffFixed, 7,4)+" MB");
  } else {
    f.close();
    r = write_vmat_compressed(vmat_fp, xdim, ydim, map_cnt);
    struct stat results;
    if (r && 0 == stat(fn.toUtf8().constData(), &results))
      Form1->Memo1->Lines->Add("Finished writing compressed vmat dump file "+fn+", size: " + 
			       FloatToStrF((float)results.st_size/(1024.0*1024.0), ffFixed, 7,4)+" MB");
  }

  time_t tf = time(NULL);
  Form1->Memo1->Lines->Add("Time required to dump vmat: "+IntToStr(tf-ti) + " seconds.");
  return r;
}

bool
write_vmat_raw(std::ofstream& f, size_t xdim, size_t ydim, int map_cnt)
{

  // f << magic << xdim << ydim
  f.write(magic, strlen(magic));
  f.write((const char*)&xdim, sizeof(xdim));
  f.write((const char*)&ydim, sizeof(ydim));
  f.write((const char*)&map_cnt, sizeof(map_cnt));

  BFOC_size_t empty = 0;
  for (size_t y=0; y<ydim; y++) {
    for (size_t x=0; x<xdim; x++) {
      if (!vmat[y][x]) {
	f.write((const char*)&empty, sizeof(empty));
      } else {
	// f << size << <list of index-value pairs>
	const Biodiv_Features_Occur_Container& occ = vmat[y][x];
	BFOC_size_t n = occ.size();
	f.write((const char*)&n, sizeof(n));
	for (BFOC_size_t spp=occ.first(); spp!=occ.overflow(); spp=occ.next(spp)) {
	  //f << spp << vmat[y][x][spp];
	  f.write((const char*)&spp, sizeof(spp));
	  f.write((const char*)&(vmat[y][x][spp]), sizeof(typeof(vmat[y][x][spp])));
	}
      }
    }
    if (0 == (y+1)%REPORT_PERIOD)
      Form1->Memo1->Lines->Add(" Written "+IntToStr(y+1)+" rows...");
  }

  f.close();
  return true;
}

bool
write_vmat_compressed(const String& fn, size_t xdim, size_t ydim, int map_cnt)
{
  //Form1->Memo1->Lines->Add("saving compressed, compress bound: "+IntToStr(compressBound(64000)));
  gzFile f = gzopen(fn.toUtf8().constData(), "w");
  if (!f) {
    Form1->Memo1->Lines->Add(" FATAL ERROR: Could not open compressed vmat file: "+fn);
    return false;
  }
#if ZLIB_VERNUM >= 0x1240
  int ok = gzbuffer(f, gz_buf_size);
  if (0 != ok) {
    Form1->Memo1->Lines->Add(" FATAL ERROR: Could not set gz buffer for file: "+fn);
    return false;
  }
#endif
  //ok = gzsetparams(f, lvl, strat);

  gzwrite(f, magic, strlen(magic));
  gzwrite(f, (const char*)&xdim, sizeof(xdim));
  gzwrite(f, (const char*)&ydim, sizeof(ydim));
  gzwrite(f, (const char*)&map_cnt, sizeof(map_cnt));

  BFOC_size_t empty = 0;  

  for (size_t y=0; y<ydim; y++) {
    for (size_t x=0; x<xdim; x++) {
      if (!vmat[y][x]) {
	gzwrite(f, (const char*)&empty, sizeof(empty));
      } else {
	// f << size << <list of index-value pairs>
	const Biodiv_Features_Occur_Container& occ = vmat[y][x];
	BFOC_size_t n = occ.size();
	gzwrite(f, (const char*)&n, sizeof(n));
	for (BFOC_size_t spp=occ.first(); spp!=occ.overflow(); spp=occ.next(spp)) {
	  //f << spp << vmat[y][x][spp];
	  gzwrite(f, (const char*)&spp, sizeof(spp));
	  gzwrite(f, (const char*)&(vmat[y][x][spp]), sizeof(typeof(vmat[y][x][spp])));
	}
      }
    }
    if (0 == (y+1)%REPORT_PERIOD)
      Form1->Memo1->Lines->Add(" Written "+IntToStr(y+1)+" rows...");
  }

  gzclose(f);
  return true;
}

bool
load_vmat(bool load_direct, const String& proj_name, const String& fn, size_t xdim, size_t ydim, int map_cnt, bool compress)
{
  // String vmat_fp = /*outdir + */ proj_name + "." + fn;
  // if (compress)
  //   vmat_fp += comp_suffix;
  String vmat_fp;
  if (load_direct) {
    vmat_fp = fn;
  } else {
    vmat_fp = make_vmat_output_name(outgridfn, fn, compress);
  }
  Form1->Memo1->Lines->Add("Loading vmat from: " + vmat_fp);

  struct stat results;
  if (0 == stat(vmat_fp.toUtf8().constData(), &results))
    Form1->Memo1->Lines->Add("Size of vmat file: " + FloatToStrF((float)results.st_size/(1024.0*1024.0),ffFixed, 4,4)+" MB");

  bool r = false;
  time_t ti = time(NULL);
  if (!compress) {
    std::ifstream f((vmat_fp).toUtf8().constData(), std::ios::in | std::ios::binary);
    if (!f) {
      Form1->Memo1->Lines->Add(" FATAL ERROR: failed to open vmat dump file for reading: " + vmat_fp);
      return false;
    }
    r = load_vmat_raw(f, xdim, ydim, map_cnt);
  } else {
    r = load_vmat_compressed(vmat_fp, xdim, ydim, map_cnt);
  }

  time_t tf = time(NULL);
  Form1->Memo1->Lines->Add("Time required to load vmat: "+IntToStr(tf-ti) + " seconds.");
  return r;
}

bool
load_vmat_raw(std::ifstream& f, size_t xdim, size_t ydim, int map_cnt)
{
  char* in_magic = new char[strlen(magic)];
  f.read((char*)in_magic, strlen(magic));
  // TODO: do some check here
  delete[] in_magic;
  size_t in_xdim, in_ydim;
  int in_map_cnt;
  f.read((char*)&in_xdim, sizeof(in_xdim));
  f.read((char*)&in_ydim, sizeof(in_ydim));  
  f.read((char*)&in_map_cnt, sizeof(in_map_cnt));

  if (in_xdim != xdim || in_ydim != ydim) {
    Form1->Memo1->Lines->Add(" FATAL ERROR: dimensions do not match. In current setup, columns: "+
			     IntToStr(xdim)+", rows: "+IntToStr(ydim)+", whereas in vmat file, columns: "+
			     IntToStr(in_xdim)+", rows: "+IntToStr(in_ydim)+".");
    return false;
  } else {
    Form1->Memo1->Lines->Add("Dimensions match, columns: "+IntToStr(xdim)+
			     ", rows: "+IntToStr(ydim));    
    Form1->Memo1->Lines->Add("Features in vmat file: "+IntToStr(in_map_cnt)+ 
			     ", features in current setup: "+IntToStr(map_cnt));
  }

  Form1->Memo1->Lines->Add("Size BFOC_size_t: "+IntToStr(sizeof(BFOC_size_t))+
			   ", size key_t: "+IntToStr(sizeof(*Biodiv_Features_Occur_Container::priv_key))+
			   ", size val_t: "+IntToStr(sizeof(*Biodiv_Features_Occur_Container::priv_val)));

  BFOC_size_t len = 0;
  for (size_t y=0; y<ydim; y++) {
    for (size_t x=0; x<xdim; x++) {
      f.read((char*)&len, sizeof(len));
      if (0 == len)
	continue;

      BFOC_size_t idx;
      float val;
      // resize the Biodiv_Features_Occur_Container& and load values
      Biodiv_Features_Occur_Container& occ = vmat[y][x];
      occ.reserve(len, len);
      occ.priv_size = len;
      for (BFOC_size_t table_idx=0; table_idx<len; table_idx++) {
	f.read((char*)&(occ.priv_key[table_idx]), sizeof(*occ.priv_key));
	f.read((char*)&(occ.priv_val[table_idx]), sizeof(*occ.priv_val));
      }
    }
    if (0 == (y+1)%REPORT_PERIOD)
      Form1->Memo1->Lines->Add(" Loaded "+IntToStr(y+1)+" rows...");
  }

  f.close();
  return true;
}

bool
load_vmat_compressed(const String& fn, size_t xdim, size_t ydim, int map_cnt)
{
  gzFile f = gzopen(fn.toUtf8().constData(), "r");
  if (!f) {
    Form1->Memo1->Lines->Add(" FATAL ERROR: Could not open compressed vmat file: "+fn);
    return false;
  }
#if ZLIB_VERNUM >= 0x1240
  int ok = gzbuffer(f, gz_buf_size);
  if (0 != ok) {
    Form1->Memo1->Lines->Add(" FATAL ERROR: Could not set gz buffer for file: "+fn);
    return false;
  }
#endif

  char* in_magic = new char[strlen(magic)];
  gzread(f, (char*)in_magic, strlen(magic));
  // TODO: do some check here
  delete[] in_magic;
  size_t in_xdim, in_ydim;
  int in_map_cnt;
  gzread(f, (char*)&in_xdim, sizeof(in_xdim));
  gzread(f, (char*)&in_ydim, sizeof(in_ydim));  
  gzread(f, (char*)&in_map_cnt, sizeof(in_map_cnt));

  if (in_xdim != xdim || in_ydim != ydim) {
    Form1->Memo1->Lines->Add(" FATAL ERROR: dimensions do not match. In current setup, columns: "+
			     IntToStr(xdim)+", rows: "+IntToStr(ydim)+", whereas in vmat file, columns: "+
			     IntToStr(in_xdim)+", rows: "+IntToStr(in_ydim)+".");
    return false;
  } else {
    Form1->Memo1->Lines->Add("Dimensions match, columns: "+IntToStr(xdim)+
			     ", rows: "+IntToStr(ydim));    
    Form1->Memo1->Lines->Add("Features in vmat file: "+IntToStr(in_map_cnt)+ 
			     ", features in current setup: "+IntToStr(map_cnt));
  }

  Form1->Memo1->Lines->Add("Size BFOC_size_t: "+IntToStr(sizeof(BFOC_size_t))+
			   ", size key_t: "+IntToStr(sizeof(*Biodiv_Features_Occur_Container::priv_key))+
			   ", size val_t: "+IntToStr(sizeof(*Biodiv_Features_Occur_Container::priv_val)));

  BFOC_size_t len = 0;
  for (size_t y=0; y<ydim; y++) {
    for (size_t x=0; x<xdim; x++) {
      gzread(f, (char*)&len, sizeof(len));
      if (0 == len)
	continue;

      BFOC_size_t idx;
      float val;
      // resize the Biodiv_Features_Occur_Container& and load values
      Biodiv_Features_Occur_Container& occ = vmat[y][x];
      occ.reserve(len, len);
      occ.priv_size = len;
      for (BFOC_size_t table_idx=0; table_idx<len; table_idx++) {
	gzread(f, (char*)&(occ.priv_key[table_idx]), sizeof(*occ.priv_key));
	gzread(f, (char*)&(occ.priv_val[table_idx]), sizeof(*occ.priv_val));
      }
    }
    if (0 == (y+1)%REPORT_PERIOD)
      Form1->Memo1->Lines->Add(" Loaded "+IntToStr(y+1)+" rows...");
  }

  gzclose(f);
  return true;
}
