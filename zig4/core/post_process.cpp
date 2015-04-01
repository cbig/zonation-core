#include "Unit1.h"
#include "GridMap.h"
#include "defines.h"
#include "post_process.h"
#include "bat_run.h"
#include "PLULA.h"
#include "zig4lib/io.h"
#include "LSIdent.h"
#include <cstdio>
#include <cmath>
#include <cstring>

// LSIcnt is used for the names of output files of LSI analyses
int  LSIcnt=1;
// This is the number of comparisons done (does not have any effect on output file names)
int lsc_cnt=0;

namespace {
class PerformAnalysis : public boost::static_visitor<bool>
{
public:
	bool operator()(const ZLSI& analysis) const
	{
	  bool success = false;
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add(" **********=====-----..... Landscape identification post-processing .....-----=====**********");
	  
	  LSIRfn = ChangeFileExt(originalOutFile, name_addition + QString(".nwout.%1.ras.asc").arg(LSIcnt));
	  LSIDfn = ChangeFileExt(originalOutFile, name_addition + QString(".nwout.%1.spp_data.txt").arg(LSIcnt));

	  top_percent = analysis.f1 / 100.0f;
	  min_percent = analysis.f2 / 100.0f;
	  max_dist = analysis.d;
	  min_simil = analysis.sim;
	  
	  success = LSIdent(0);		
	  ++LSIcnt;

	  return success;
	}

	bool operator()(const ZLSC& analysis) const
	{
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("Landscape comparison post-processing");

	  LSCAnalysis(analysis.f1, analysis.f2, analysis.comp_solution, analysis.output_file);
	  
	  Form1->Memo1->Lines->Add("Landscape comparison post-processing finished successfully");
	  Form1->Memo1->Lines->Add("");

	  ++lsc_cnt;
	  return true;
	}

	bool operator()(const ZLSM& analysis) const
	{
	  bool success = false;
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("Post-processing landscape analysis for masked-in areas");
	  
	  LSIRfn = ChangeFileExt(originalOutFile, name_addition + QString(".nwout.%1.ras.asc").arg(LSIcnt));
	  LSIDfn = ChangeFileExt(originalOutFile, name_addition + QString(".nwout.%1.spp_data.txt").arg(LSIcnt));
	  
	  LSImask = analysis.mask_file;
	  min_percent = analysis.f2 / 100.0f;
	  max_dist = analysis.d;
	  min_simil = analysis.sim;
	  
	  success = LSIdent(1);
	  ++LSIcnt;

	  return success;
	}

	bool operator()(const ZLSB& analysis) const
	{
	  bool success = false;
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("Post-processing landscape analysis for top fraction inside masked-in areas");

	  LSIRfn = ChangeFileExt(originalOutFile, name_addition + QString(".nwout.%1.ras.asc").arg(LSIcnt));
	  LSIDfn = ChangeFileExt(originalOutFile, name_addition + QString(".nwout.%1.spp_data.txt").arg(LSIcnt));

	  LSImask = analysis.mask_file;
	  top_percent = analysis.f1 / 100.0f;
	  min_percent = analysis.f2 / 100.0f;
	  max_dist = analysis.d;
	  min_simil = analysis.sim;

	  success = LSIdent(2);
	  ++LSIcnt;

	  return success;
	}
};
}

int post_process(const String& PPA_fname)
{
	using namespace boost;
	ZPPAFile ppaFile;
	SuppressErrors suppress;
	if(!loadFile(loadZPPAFile, ppaFile, PPA_fname, suppress)) {
	  if (PPA_fname.isEmpty())
	    Form1->Memo1->Lines->Add("No automated post-processing file specified");
	  else 
	    Form1->Memo1->Lines->Add("Automated post-processing file ('" + 
				     PPA_fname + "') not found!");
	  return 0;
	}
	int cnt(0);
	PerformAnalysis analyze;
	foreach(const shared_ptr<ZAnalysis>& analysis, ppaFile.analysisList) {
		if(apply_visitor(analyze, *analysis)) {
			++cnt;
		}
	}

	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("From the "+IntToStr((LSIcnt-1) + lsc_cnt)+" post-processing analyses requested, "+IntToStr(cnt)+
				 " were succesfully done, including "+IntToStr(lsc_cnt)+" landscape comparison(s) (LSC).");
	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("");
	return cnt;
}
