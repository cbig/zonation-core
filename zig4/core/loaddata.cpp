#include "Unit1.h"
#include "LoadData.h"
#include "GridMap.h"
#include "defines.h"
#include "matrix_utils.h"
#include "ccsp.h"
#include "randz.h"
#include "ADMUs.h"
#include "bat_run.h"
#include "PLULA.h"
#include "LSIdent.h"
#include "smooth.h"
#include "output.h"     // need generate_output_trans_layers, etc.
#include "mem_dump.h"
#include "zig4lib/io.h"
#include "zig4lib/ini.h"
#include "zig4lib/zconf.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <map>


#include <QCoreApplication>

#ifdef WIN32
#define strnlen(str, len) strlen(str)
#endif

#define STR_TRANS_HELPER(x) #x
#define INT_TO_STR(x) STR_TRANS_HELPER(x)

const size_t MB_in_bytes = 1048576;

// interactions (fedemp: nope, only used here in loaddata.cpp)
float glob_wmax=1.0f;

float initial_wmax = 1.0f;  // wmax including all features, before dropping "0 occurrence" ones

// general settings
int resample_count = 0;
bool settings_drop_0_occur_layers = false;
// new indices of bd features after dropping some of them (0 occurrence within analysis area ones)
std::vector<int> drop_0_relocation_idx;
float random_omissions_percentage_after_vmat_simple_impl = 0.0f;
int random_omissions_max_feature = 0;

// rr4
int  new_spp_file_format=0;

// condition xxx2.1
class  GridMap condition; // also used for retention v3

// retention xxx3.0
// class  GridMap retention;

// post processing
String PPA_fname;

// community settings
struct  community_settings comm_set;
float   **sim_mat, **comm_sim_mat;
float   ret_layers_rel_weight;

// area mask
class GridMap area_mask;

// Current behavior: if the "cell area correction" layer is enabled, these 3 aspects will be corrected:
//  - occurrence values (weighted by cell area)
//  - cell cost
//  - ranking map (rank levels are for area rather than cells)
// In practice these 3 bool variables could be merged into one. Left as 3 different variables for now
// so it would be possible to adjust these 3 different aspects independently.
String glob_set_occur_size_weights_layer_fn;
bool use_occur_size_weights_correct_landscape_fraction = false;
bool use_occur_size_weights_correct_cost = false;
bool use_occur_size_weights_correct_ranking = false;
class GridMap occur_size_weights_map;

const size_t MAX_LINE_LEN = 64000;

// transformed layers generated as outputs
struct transf_layers_settings transf_layers_set;

String load_precooked_vmat_fn;
bool load_vmat_directly = false;
String save_vmat_fn;

static ZErrorCallBack zErrorCallBack;
static ZErrorStderrCallBack zErrorStderrCallBack;

void zig_legacy_alloc_error(char error_text[])
{
  Form1->Memo1->Lines->Add("ERROR: memory allocation error (see above for more details, if available). Most likely you run out of available RAM memory.");
  Form1->Memo1->Lines->Add("zig internal description of error: " + String(error_text));
  Form1->Memo1->Lines->Add("There is no solution for this, giving up.");
  graceful_exit();
}

void Fix_fname(String &fname) {} // do nothing

void Fix_output_file_names(ZConf const& conf)
{
	emffn         = conf.emfOutputPath();
	memofn        = conf.runinfoOutputPath();

	features_info_fn = conf.featuresInfoOutputPath();
	curvesfn      = conf.curvesOutputPath();

	SSI_features_info_fn  = conf.ssiFeaturesInfoOutputPath();
	SSI_curvesfn  = conf.ssiCurvesOutputPath();

	grp_curves_fn = conf.groupCurvesOutputPath();
	outgridfn     = conf.rankOutputPath();
	outgridfn2    = conf.propOutputPath();
	outgridfn3    = conf.wrscrOutputPath();

	// admu
	admu_redist_outgridfn = conf.admu_redist_rankOutputPath();

	output_transf_ig_layers_dn = conf.transf_ig_layers_DirPath();
	// output of ds transformed layers, and so on and so forth...
	output_transf_ds_layers_dn = conf.transf_ds_layers_DirPath();
	output_transf_cst_layers_dn = conf.transf_cst_layers_DirPath();
	output_transf_ct_layers_dn = conf.transf_ct_layers_DirPath();
	output_transf_rt_layers_dn = conf.transf_rt_layers_DirPath();
	output_transf_mct_layers_dn = conf.transf_mct_layers_DirPath();
	// output of interactions (IA) transformed layers
	output_transf_ia_layers_dn = conf.transf_ia_layers_DirPath();

	// output of final transformed layers
	output_transf_final_layers_dn = conf.transformed_final_layers_DirPath();
}

namespace {
void printCmdLineErr(int num_args)
{
  printf("Incorrect number of parameters (" + IntToStr(num_args) + ") on command line\n");
  printf("Expected run_mode settings_file input_grid_list_file output_grid_file IGa use_smooth smooth_mult autoclose_0_or_1\n");
  printf("Type zig4 -h or --help to get a quick help, or check the manual, part 3: The Zonation software.\n");
}
}

bool Parse_cmd_line(ZCommandLine& command)
{
	QStringList args(QCoreApplication::arguments());
	if(!parseZCommandLine(command, args, zErrorStderrCallBack)) {
	  // This would say again that it cannot run, even for -h, --help, -v and other options for which it's ok to just exit
	  // printCmdLineErr(args.size());
	  return false;
	}

	if(args[1] == "-r") {
	  run_mode = 1;
	} else if(args[1].startsWith("-l")) { // args[1] = "-lfoo.asc"
	  run_mode = 2;
	  load1fn = args[1].mid(2);
	} else {
	  printCmdLineErr(args.size());
	  return false;
	}

	// Initialize VCL
	initializeVCL(command);

	IG_set.IGa = static_cast<float>(command.uncertaintyAlpha);
	use_smoothing = command.distributionSmoothingOn ? 1 : 0;
	alpha_mult = static_cast<float>(command.dispersalKernelMultiplier);
	autoclose = command.windowLeftOpen;

	setfn     = command.datFile;
	Fix_fname(setfn);
	sppfn  = command.sppFile;
	Fix_fname(sppfn);

	outgridfn  = command.outFile;
	originalOutFile = command.outFile;

	/*
 // cmd set ip op close
 char tmp[256];

 if (_argc<9)
 {
  Form1->Memo1->Lines->Add("Incorrect number of parameters (" + IntToStr(_argc) + ")on command line");
  Form1->Memo1->Lines->Add("Expected run_mode settingsfile input_grid_file output_grid_file IGa useSmooth SmoothMult autoclo
se0_or_1");
  if (_argc==1)
   Form1->Memo1->Lines->Add("In the absence of run parameters, starting GUI parameter input form.");
  //Delay(10);
  return false;
 }

 if (strstr(_argv[1], "-r"))
  run_mode = 1; // running
 else
 {
  run_mode = 2;
  strcpy(tmp,_argv[1]+2);
  load1fn = tmp;
 }

 IG_set.IGa = static_cast<float>(atof(_argv[5]));
 use_smoothing = atoi(_argv[6]);
 alpha_mult = static_cast<float>(atof(_argv[7]));
 autoclose = atoi(_argv[8]) == 0 ? false : true;

 setfn     = _argv[2];
 Fix_fname(setfn);
 sppfn  = _argv[3];
 Fix_fname(sppfn);

 outgridfn  = _argv[4];

 char *temp = sppfn.toUtf8().constData();
*/
	return true;
}

void show_run_pars_on_form1()
{ // xxxia
	if (run_mode==1)
		Form1->RadioGroup2->ItemIndex=0;
	else
		Form1->RadioGroup2->ItemIndex=1;

	//  if (removal_rule==1)
	Form1->RadioGroup3->ItemIndex=removal_rule-1;
	//  else
	//    Form1->RadioGroup3->ItemIndex=1;

	Form1->Edit9->Text  = BQP_prof_fn;
	Form1->Edit10->Text = sppfn;
	Form1->Edit11->Text = costfn;
	Form1->Edit27->Text = maskfn;
	Form1->Edit17->Text = load1fn;
	Form1->Edit1->Text  = SSI_fname;
	Form1->Edit30->Text = PLULAfn;
	Form1->Edit40->Text = ia_fname;
	Form1->Edit41->Text = tree_fname;

	if (use_BQP)
		Form1->CheckBox8->Checked = true;
	else
		Form1->CheckBox8->Checked = false;

	if (use_edge_rem)
		Form1->CheckBox9->Checked = true;
	else
		Form1->CheckBox9->Checked = false;

	if (annotate_name)
		Form1->CheckBox1->Checked = true;
	else
		Form1->CheckBox1->Checked = false;

	if (use_smoothing)
		Form1->CheckBox10->Checked = true;
	else
		Form1->CheckBox10->Checked = false;

	if (logit)
		Form1->CheckBox11->Checked = true;
	else
		Form1->CheckBox11->Checked = false;

	if (BQP_mode==1)
		Form1->RadioGroup4->ItemIndex=0;
	else
		Form1->RadioGroup4->ItemIndex=1;

	if (use_cost)
	{
		Form1->CheckBox2->Checked=true;
		//      Form1->Chart1->BottomAxis->Title->Caption = "Proportion of landscape cost lost";
		Form1->Memo1->Lines->Add("*************** USING COST LAYER **************");
		Form1->Memo1->Lines->Add("  Cost file = "+costfn);
	}
	else
	{
		Form1->CheckBox2->Checked=false;
	}

	if (use_mask)
	{
		Form1->CheckBox4->Checked=true;
		Form1->Memo1->Lines->Add("*************** USING (HIERARCHICAL) REMOVAL MASK LAYER **************");
		Form1->Memo1->Lines->Add("  Mask file: " + maskfn);
		Form1->Memo1->Lines->Add("");
	}
	else
		Form1->CheckBox4->Checked=false;

	if (use_SSI)
	{
		Form1->CheckBox5->Checked=true;
		Form1->Chart4->Visible=true;
		Form1->Memo1->Lines->Add("*************** USING SSI features/species **************");
		Form1->Memo1->Lines->Add("  SSI list file = "+SSI_fname);
	}
	else
	{
		Form1->CheckBox5->Checked=false;
		Form1->Chart4->Visible=false;
	}

	if (use_interactions)
	{
		Form1->CheckBox15->Checked=true;
		Form1->Memo1->Lines->Add("*************** USING feature/species interactions **************");
		Form1->Memo1->Lines->Add("  IA list file = "+ia_fname);
	}
	else
	{
		Form1->CheckBox15->Checked=false;
	}

	if (corr_set.use_corr) {
	  Form1->CheckBox15->Checked=true;
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("*************** Using corridor connectivity with these settings: **************");
	  Form1->Memo1->Lines->Add("  Strength: " + FloatToStr(CCSP)); /*FloatToStrF(CCSP, ffFixed, 7, 4)*/ 
	  Form1->Memo1->Lines->Add("  Width: " + IntToStr(CCSP_thickness));
	  String var_text;
	  if (1 == CCSP_variant)
	    var_text = "Fast/approximate";
	  else if (2 == CCSP_variant)
	    var_text = "Exact/potentially slow";	    
	  Form1->Memo1->Lines->Add("  Penalty calculation variant: " + var_text);
	  if (corr_set.out_boundaries_pcs.size()>0) {
	    String fract_list = "";
	    for(size_t i=0; i<corr_set.out_boundaries_pcs.size(); i++) {
	      fract_list += FloatToStr(corr_set.out_boundaries_pcs[i])+" ";
	    }
	    Form1->Memo1->Lines->Add("  Save layers at top %: " + fract_list);
	  }
	  Form1->Memo1->Lines->Add("  Starting corridor building at top "+FloatToStr(corr_set.start_pc)+"\%");
	  if (corr_set.layers_fn.isEmpty())
	    Form1->Memo1->Lines->Add("  No domain layers file given.");
	  else
	    Form1->Memo1->Lines->Add("  Domain layers file: " + corr_set.layers_fn);
	  Form1->Memo1->Lines->Add("");
	  /*
	  if (warp_factor>1) {
	    Form1->Memo1->Lines->Add("********************* WARNING ***************");
	    Form1->Memo1->Lines->Add("   WARNING: Using warp>1 may not work properly when using the corridor loss penalty.");
	  }
	  */
	}

	if (use_PLULA)
	{
		Form1->CheckBox6->Checked=true;
		Form1->Memo1->Lines->Add("*************** USING PLULA **************");
		Form1->Memo1->Lines->Add("  Planning units file = "+PLULAfn);
		if (use_tree_conn)
		{
			Form1->Memo1->Lines->Add("*************** USING Tree connectivity measure with PLULAs ********");
			Form1->Memo1->Lines->Add("   Tree specification file = "+tree_fname);
			Form1->CheckBox16->Checked=true;
			if (use_BQP)
			{
				Form1->Memo1->Lines->Add("   WARNING: BQP switched off as tree connectivity model is used = "+tree_fname);
				use_BQP=false;
			}
		}
		else
		{
			Form1->CheckBox16->Checked=false;
		}
	}
	else
	{
		Form1->CheckBox6->Checked=false;
	}

	Form1->Edit13->Text=FloatToStrF(rem_level, ffFixed, 7,3);
	Form1->Edit14->Text=FloatToStrF(alpha_mult, ffFixed, 7,3);
	Form1->Edit18->Text=FloatToStrF(BLP, ffFixed, 7,3);
	Form1->LabeledEdit1->Text=FloatToStrF(SA_z, ffFixed, 7,3);
	Form1->LabeledEdit2->Text=IntToStr(resample_count);
	Form1->Edit25->Text=IntToStr(warp_factor);
	Form1->Edit19->Text=IntToStr(add_edge_cnt);
	//  Form1->Edit30->Text=FloatToStrF(rem_rule_par, ffFixed, 7,3);

	if(IG_set.use_IG)
		Form1->CheckBox12->Checked=true;
	else
		Form1->CheckBox12->Checked=false;

	Form1->Edit15->Text=  FloatToStrF(IG_set.IGa, ffFixed, 7,3);

	Form1->RadioGroup1->ItemIndex = IG_set.IG_proportional;

	Form1->Edit16->Text = IG_set.IGwfn;

	if(IG_set.use_IGw)
		Form1->CheckBox13->Checked=true;
	else
		Form1->CheckBox13->Checked=false;

	if (IG_set.normalize_IGw)
		Form1->CheckBox14->Checked=true;
	else
		Form1->CheckBox14->Checked=false;
}

bool Read_settings_file(ZDATFile& dat)
{
	//TIniFile *pIniFile;
	//  AnsiString tsection="Settings";
	FILE  *tmpf;
	char txt[128];

	Form1->Memo1->Lines->Add("Zonation core process parameters: ");
	if (1 == run_mode)
	  Form1->Memo1->Lines->Add("   run mode: normal run (-r)");
	else if (2==run_mode)
	  Form1->Memo1->Lines->Add("   run mode: re-load (-l)");
	else
	  Form1->Memo1->Lines->Add("   run mode: "+IntToStr(run_mode));
	Form1->Memo1->Lines->Add("   settings file: " + setfn);
	Form1->Memo1->Lines->Add("   biodiversity features list file: " + sppfn);
	Form1->Memo1->Lines->Add("   output file(s): " + ChangeFileExt(outgridfn, ""));
	Form1->Memo1->Lines->Add("   info-gap alpha: " + FloatToStr(IG_set.IGa));
	Form1->Memo1->Lines->Add("   use_smoothing: " + IntToStr(use_smoothing));
	Form1->Memo1->Lines->Add("   smoothing_alpha_multiplier: " + FloatToStr(alpha_mult));
	Form1->Memo1->Lines->Add("   autoclose (deprecated): " + IntToStr(autoclose));
	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("Reading settings file: " + setfn);

	if(!loadFile(loadZDATFile, dat, setfn, zErrorCallBack)) {
	  return false;
	}

	//  normalize  = (pIniFile->ReadInteger("Settings", "normalize", 1));
	//  show_images= (pIniFile->ReadInteger("Settings", "graphics", 1));
	//  rem_by_rarity  = (pIniFile->ReadInteger("Settings", "remove by rarity", 1));
	//  rem_step_by_rarity  = (pIniFile->ReadInteger("Settings", "removal step by rarity", 1));
	//  miss_as_zero  = (pIniFile->ReadInteger("Settings", "treat missing as zero", 0));

	normalize         = 1;
	show_images       = 1;
	rem_by_rarity     = 1;
	rem_step_by_rarity= 1;
	miss_as_zero      = 0;

	// regex:
	// \t([^\s=]*)\s*=[^,]*,[^,]*,\s*([^)]*).*
	// \1 = dat.\1 ? *dat.\1 : \2
	rem_level     = dat.rem_level ? *dat.rem_level : 0.0f;
	use_cost = dat.use_cost ? *dat.use_cost :  0 ;
	use_mask = dat.use_mask ? *dat.use_mask :  0 ;
	glob_add_to_edge_borders_between_removal_mask_levels = dat.add_to_edge_borders_between_removal_mask_levels? 
	  *dat.add_to_edge_borders_between_removal_mask_levels : true;
	use_BQP = dat.use_BQP ? *dat.use_BQP :  0 ;
	BQP_mode = dat.BQP_mode ? dat.BQP_mode->value() :  2 ;
	BLP = dat.BLP ? *dat.BLP :  0.0f ;
	SA_z = dat.SA_z ? *dat.SA_z :  0.25f ;
	use_edge_rem = dat.use_edge_rem ? *dat.use_edge_rem :  1 ;
	annotate_name = dat.annotate_name ? *dat.annotate_name : 1 ;
	resample_count = dat.resample_count ? *dat.resample_count : 0 ;
	settings_drop_0_occur_layers = dat.drop_0_occur_layers ? *dat.drop_0_occur_layers : false ;

	// removal rule: command line overrides .dat
	int rr_opt = VCLCommandLine::removal_rule();
	if ( rr_opt > 0) {
	  removal_rule = rr_opt;
	  Form1->Memo1->Lines->Add("   Option overriden in command line, removal rule: " + IntToStr(removal_rule));
	} else {
	  removal_rule = dat.removal_rule ? dat.removal_rule->value() :  1 ;
	//  rem_rule_par  = (pIniFile->ReadFloat("Settings",   "removal parameter", 1.0f)); ;
	  Form1->Memo1->Lines->Add("   Option not overriden in command line, removal rule: " + IntToStr(removal_rule));
	}

	costfn = dat.costfn ? *dat.costfn :  "" ;
	maskfn = dat.maskfn ? *dat.maskfn :  "" ;
	BQP_prof_fn = dat.BQP_prof_fn ? *dat.BQP_prof_fn :  "" ;

	// warp factor: command line overrides .dat
	int warp_opt = VCLCommandLine::warp_factor();
	if ( warp_opt > 0 ) {
	  warp_factor = warp_opt;
	  Form1->Memo1->Lines->Add("   Option overriden in command line, warp factor: " + IntToStr(warp_factor));
	} else {
	  warp_factor = dat.warp_factor ? *dat.warp_factor :  1 ;
	  Form1->Memo1->Lines->Add("   Option not overriden in command line, warp factor: " + IntToStr(warp_factor));
	}

	glob_set_output_wrscr_map = dat.output_richness? *dat.output_richness : true;
	glob_set_output_prop_rank_map = dat.output_prop_rank? *dat.output_prop_rank : false;

	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("Output settings:");
	String out_wrscr = glob_set_output_wrscr_map? "yes" : "no";
	String out_prop = glob_set_output_prop_rank_map? "yes" : "no";
	Form1->Memo1->Lines->Add("   Output weighted range size corrected richness map: "+out_wrscr);
	Form1->Memo1->Lines->Add("   Output proportional loss rank map: "+out_prop);
	Form1->Memo1->Lines->Add("");

	logit = dat.logit ? *dat.logit :  0 ;
	add_edge_cnt = dat.add_edge_cnt ? *dat.add_edge_cnt :  0 ;
	use_SSI = dat.use_SSI ? *dat.use_SSI :  0 ;
	SSI_fname = dat.SSI_fname ? *dat.SSI_fname :  "" ;
	PPA_fname = dat.PPA_fname ? *dat.PPA_fname :  "" ;
	all_zero_as_miss = dat.all_zero_as_miss ? *dat.all_zero_as_miss :  0 ;
	mem_save_mode = dat.mem_save_mode ? *dat.mem_save_mode :  0 ;
	use_groups = dat.use_groups ? *dat.use_groups :  0 ;
	groups_fname = dat.groups_fname ? *dat.groups_fname :  "" ;
	use_condition = dat.use_condition ? *dat.use_condition :  0 ;
	condition_fname = dat.condition_fname ? *dat.condition_fname :  "" ;
	use_retention = dat.use_retention ? *dat.use_retention :  0 ;
	retention_fname = dat.retention_fname ? *dat.retention_fname :  "" ;
	// use_LEC = dat.use_LEC ? *dat.use_LEC :  0 ;
	// NOTE: LEC always disabled!
	use_LEC = false;
	LEC_fname = dat.LEC_fname ? *dat.LEC_fname :  "" ;
	ret_layers_rel_weight = dat.ret_layers_rel_weight ? *dat.ret_layers_rel_weight :  1.0 ;

	mask_data = dat.mask_data ? *dat.mask_data :  0;
	area_mask_fname = dat.area_mask_fname ? *dat.area_mask_fname :  "";

	random_omissions_percentage_after_vmat_simple_impl = dat.random_omissions_percentage_after_vmat_simple_impl ? 
	  *dat.random_omissions_percentage_after_vmat_simple_impl : 0.0;
	random_omissions_max_feature = dat.random_omissions_max_feature ?
	  *dat.random_omissions_max_feature : 0;

	// xxxV3 check all gui updates!!
	// GUI, file, memo
	comm_set.load_sim_matrix = dat.load_sim_matrix ? *dat.load_sim_matrix :   0;
	comm_set.sim_matrix_name = dat.sim_matrix_name ? *dat.sim_matrix_name :  "";
	comm_set.apply_to_conn = dat.apply_to_conn ? *dat.apply_to_conn :    0;
	comm_set.comm_sim_matrix_name = dat.comm_sim_matrix_name ? *dat.comm_sim_matrix_name :  "";
	comm_set.apply_to_repr = dat.apply_to_repr ? *dat.apply_to_repr :  0;
	comm_set.edge_effect_fn = dat.edge_effect_fn ? *dat.edge_effect_fn :  "";

	ADM_set.use_ADMUs = dat.use_ADMUs ? *dat.use_ADMUs :   0;
	ADM_set.ADMU_mode = dat.ADMU_mode ? dat.ADMU_mode->value() :   1;
	ADM_set.calc_from_condition = dat.calc_from_condition ? *dat.calc_from_condition :   0;
	ADM_set.ADM_layer_file = dat.ADM_layer_file ? *dat.ADM_layer_file :   "";
	ADM_set.ADM_weights_file = dat.ADM_weights_file ? *dat.ADM_weights_file :   "";
	ADM_set.ADM_weight_matrix_file = dat.ADM_weight_matrix_file ? *dat.ADM_weight_matrix_file :   "";
	ADM_set.mode2_global_weight = dat.mode2_global_weight ? *dat.mode2_global_weight :   0.0;
	ADM_set.row_count_for_per_admu_curves = dat.ADM_row_count_for_per_admu_curves ? *dat.ADM_row_count_for_per_admu_curves : 0;

	// yes, there should be an array of these...
	transf_layers_set.output_ig_layers = dat.output_ig_layers ? *dat.output_ig_layers : 0;
	transf_layers_set.output_ds_layers = dat.output_ds_layers ? *dat.output_ds_layers : 0;
	transf_layers_set.output_cst_layers = dat.output_cst_layers ? *dat.output_cst_layers : 0;
	transf_layers_set.output_ct_layers = dat.output_ct_layers ? *dat.output_ct_layers : 0;
	transf_layers_set.output_rt_layers = dat.output_rt_layers ? *dat.output_rt_layers : 0;
	transf_layers_set.output_mct_layers = dat.output_mct_layers ? *dat.output_mct_layers : 0;
	transf_layers_set.output_ia_layers = dat.output_ia_layers ? *dat.output_ia_layers : 0;
	transf_layers_set.output_final_layers = dat.output_final_layers ? *dat.output_final_layers : 0;
	transf_layers_set.some_output = transf_layers_set.output_ig_layers || transf_layers_set.output_ds_layers
	  || transf_layers_set.output_cst_layers || transf_layers_set.output_ct_layers || transf_layers_set.output_rt_layers
	  || transf_layers_set.output_mct_layers || transf_layers_set.output_ia_layers 
	  || transf_layers_set.output_final_layers;
	// xxxPLULA
	use_PLULA = dat.use_PLULA ? *dat.use_PLULA :  0;
	PLULAfn = dat.PLULAfn ? *dat.PLULAfn :  "";

	// Tree conn
	use_tree_conn = dat.use_tree_conn ? *dat.use_tree_conn :  0;
	tree_fname = dat.tree_fname ? *dat.tree_fname :  "";

	// interactions
	use_interactions = dat.use_interactions ? *dat.use_interactions :  0;
	ia_fname = dat.ia_fname ? *dat.ia_fname :  "";

	IG_set.IG_proportional = dat.IG_proportional ? *dat.IG_proportional :  0;
	IG_set.use_IGw = dat.use_IGw ? *dat.use_IGw :  0;
	IG_set.IGwfn = dat.IGwfn ? *dat.IGwfn :  "";
	IG_set.normalize_IGw = dat.normalize_IGw ? *dat.normalize_IGw :  0;

	// arbitrary kernels
	use_arb_kernels = dat.arb_kernels_use ? *dat.arb_kernels_use : 0;
	arb_kernels_default = dat.arb_kernels_default ? *dat.arb_kernels_default : 1;
	arb_kernels_constant = dat.arb_kernels_constant ? *dat.arb_kernels_constant : 1;
	arb_kernels_prefix = dat.arb_kernels_prefix ? *dat.arb_kernels_prefix : "";

	load_precooked_vmat_fn = dat.load_precooked_vmat ? *dat.load_precooked_vmat :  "" ;
	load_vmat_directly = dat.load_vmat_directly ? *dat.load_vmat_directly : false;
	save_vmat_fn = dat.save_vmat ? *dat.save_vmat :  "" ;

	glob_set_occur_size_weights_layer_fn = dat.occur_size_weights_layer_fn? *dat.occur_size_weights_layer_fn : "";
	if (!glob_set_occur_size_weights_layer_fn.isEmpty()) {
	  use_occur_size_weights_correct_cost = true;
	  use_occur_size_weights_correct_landscape_fraction = true;
	  use_occur_size_weights_correct_ranking = true;
	}

	if (!load_precooked_vmat_fn.isEmpty()) {
	  String fn;
	  if (load_vmat_directly) {
	    fn = load_precooked_vmat_fn;
	  } else {
	    fn = make_vmat_output_name(outgridfn, load_precooked_vmat_fn);
	  }
	  Form1->Memo1->Lines->Add("*************** ENABLED LOAD-VMAT: feature layers will not be loaded from individual raster files, memory dump file '"+
				   load_precooked_vmat_fn+"' will be used instead (see below) **************");
	  Form1->Memo1->Lines->Add("");
	}
	if (!save_vmat_fn.isEmpty()) {
	  Form1->Memo1->Lines->Add("*************** ENABLED SAVE-VMAT: after loading and transforming feature layers, the resulting structure will be dumped into file '"+
				   make_vmat_output_name(outgridfn, save_vmat_fn)+"' **************");
	  Form1->Memo1->Lines->Add("");
	}

	if (use_arb_kernels) {
	  Form1->Memo1->Lines->Add("********** Enabling arbitrary kernels ********");
	  Form1->Memo1->Lines->Add("******* The default kernel is: "+IntToStr(arb_kernels_default)+" *****");
	  Form1->Memo1->Lines->Add("");
	} else {
	  Form1->Memo1->Lines->Add("********** Arbitrary kernels not enabled ********");
	  Form1->Memo1->Lines->Add("");
	}

	// corridors
	corr_set.use_corr = dat.use_corridors? *dat.use_corridors : false;
	CCSP = dat.CCSP ? *dat.CCSP :  0.0f ;
	CCSP_thickness = dat.CCSP_thickness ? *dat.CCSP_thickness :  0 ;
	CCSP_formula = dat.CCSP_formula ? *dat.CCSP_formula :  1 ;
	CCSP_variant = dat.CCSP_variant ? *dat.CCSP_variant :  1 ;
	CCSP_info_period = dat.CCSP_info_period ? *dat.CCSP_info_period :  10000 ;
	String CCSP_corridor_out_boundaries_pcs = dat.corridor_output_boundaries_pcs ? *dat.corridor_output_boundaries_pcs : "";
	corr_set.use_domain_layers = dat.use_corridor_domain_layers? *dat.use_corridor_domain_layers : false;
	corr_set.layers_fn = dat.corridor_domain_layers_fn ? *dat.corridor_domain_layers_fn : "";
	corr_set.start_pc = dat.corridor_start_pc ? *dat.corridor_start_pc : 100.0f;
	output_ccsp_fn = ChangeFileExt(outgridfn, "");	

	// What effectively enables/disables corridors is "0.0f == CCSP" (as with BLP)
	if (!corr_set.use_corr)
	  CCSP = .0f;

	if (!CCSP_corridor_out_boundaries_pcs.isEmpty()) {
	  QStringList top_fractions = CCSP_corridor_out_boundaries_pcs.split(" ");
	  for (int i = 0; i < top_fractions.size(); i++) {
	    corr_set.out_boundaries_pcs.push_back(top_fractions[i].toFloat());
	  }
	  // sort in case, descending order
	  std::sort(corr_set.out_boundaries_pcs.begin(), corr_set.out_boundaries_pcs.end(), std::greater<float>());
	}

	if (2 == run_mode && use_edge_rem) {
	  Form1->Memo1->Lines->Add("    *** Note: running in re-loading mode (-l in command line), \"edge removal\" is ignored and \"add edge points\" is set to 0 (default) even though the settings file specifies different values for these options ***");
	  Form1->Memo1->Lines->Add("");
	  use_edge_rem = false;
	  add_edge_cnt = 0;
	}

	if (!glob_set_occur_size_weights_layer_fn.isEmpty()) {
	  Form1->Memo1->Lines->Add("*************** USING MAP OF CELL AREAS / OCCURRENCE WEIGHTS **************");
	  Form1->Memo1->Lines->Add("********        Occurrence weights raster map file: " + glob_set_occur_size_weights_layer_fn);
	  Form1->Memo1->Lines->Add("********        Occurrence values of biodiversity features will be adjusted by these weights.");
	  if (use_occur_size_weights_correct_landscape_fraction)
	    Form1->Memo1->Lines->Add("********        In output files, the fraction of landscape remaining/lost will be adjusted by the cell weights/sizes.");
	  if (use_occur_size_weights_correct_ranking)
	    Form1->Memo1->Lines->Add("********        The output ranking will contain fractions of area rather than fractions of cells.");
	  Form1->Memo1->Lines->Add("");
	}

	if (mask_data) {
	  Form1->Memo1->Lines->Add("*************** MASKING ANALYSIS AREAS FROM INPUT FILES **************");
	  Form1->Memo1->Lines->Add("********        Masking according to file "+area_mask_fname);
	  Form1->Memo1->Lines->Add("");
	} else {
	  Form1->Memo1->Lines->Add("*************** NOT USING ANALYSIS AREA MASK **************");
	  Form1->Memo1->Lines->Add("");	  
	}

	if (use_condition)
	{
		if (!use_groups)
		{
			Form1->Memo1->Lines->Add("WARNING: Condition cannot be used if groups have not been specified");
			Form1->Memo1->Lines->Add("WARNING: Condition has been switched off");
			Form1->Memo1->Lines->Add("");
			use_condition=0;
		}
		Form1->Memo1->Lines->Add("");
		Form1->Memo1->Lines->Add("******* Condition layers used. Will load the list of condition layers from file: "+ condition_fname);
		Form1->Memo1->Lines->Add("");
	}

	if (ADM_set.use_ADMUs)
	{
		if (ADM_set.calc_from_condition)
		{
			if (!use_condition)
			{
				Form1->Memo1->Lines->Add("ADMU feature weights cannot be computed from condition because groups and/or condition are not used.");
				Form1->Memo1->Lines->Add("WARNING ADMU weight calculation has been switched off.");
				Form1->Memo1->Lines->Add("");
				ADM_set.calc_from_condition=0;
			}
		}
		if (ADM_set.ADMU_mode==1)
			Form1->Memo1->Lines->Add("ADMU mode = 1");
		else
			Form1->Memo1->Lines->Add("ADMU mode = 2; Wg="+FloatToStr(ADM_set.mode2_global_weight));

	}

	if (use_retention)
	{
		if (!use_groups)
		{
			Form1->Memo1->Lines->Add("Retention cannot be used if groups have not been specified");
			Form1->Memo1->Lines->Add("WARNING Retention has been switched off");
			Form1->Memo1->Lines->Add("");
			use_retention=0;
		}
		Form1->Memo1->Lines->Add("Retention used; relative weight multiplier = "+ FloatToStr(ret_layers_rel_weight));
	}

	if (logit)
	{
		Form1->Memo1->Lines->Add("********** USING logit space ************");
		Form1->Memo1->Lines->Add("");
		//  use_smoothing = (pIniFile->ReadInteger("Settings", "use smoothing", 0));
		//  alpha_mult    = (pIniFile->ReadFloat("Settings", "multiply alphas by", 1.0f));
	}

	if (mem_save_mode)
	{
		Form1->Memo1->Lines->Add("********** USING memory saving mode - limited analysis options ************");
		Form1->Memo1->Lines->Add("");
	}

	if (use_BQP)
	{
		Form1->Memo1->Lines->Add("*********** USING Boundary Quality Penalty *********");
		if (BQP_mode==1)
			Form1->Memo1->Lines->Add("  BQP MODE=1: BQP aligns missing data");
		else
			Form1->Memo1->Lines->Add("  BQP MODE=2: BQP uses feature-specific neighborhoods");
		Form1->Memo1->Lines->Add("");
	}

	if (BLP>0.0f)
	{
		Form1->Memo1->Lines->Add("********* USING the Boundary Length Penalty; BLP="+FloatToStrF(BLP, ffFixed, 7, 4));
		Form1->Memo1->Lines->Add("");
		// Old error message
		// if (warp_factor>1)
		// {
		// 	Form1->Memo1->Lines->Add("********************* ERROR ***************");
		// 	Form1->Memo1->Lines->Add("   WARNING: Using warp>1 does not work properly when using the BLP.");
		// 	Form1->Memo1->Lines->Add("");
		// }
	}

	//  if (add_edge_cnt)
	//    {
	//      Form1->Memo1->Lines->Add("*********** Adding fake edge points; count = "+IntToStr(add_edge_cnt));
	//      Form1->Memo1->Lines->Add("");
	//    }

	const int max_warp= 1000000;
	if (warp_factor<1) {
	  warp_factor = 1;
	} else if (warp_factor>max_warp) {
	  Form1->Memo1->Lines->Add("WARNING: too big warp factor requested, setting it to maximum allowed: "+IntToStr(max_warp));
	  warp_factor = max_warp;
	}

	if (IG_set.IGa!=0.0f)
		IG_set.use_IG=1;
	else
		IG_set.use_IG=0;

	//  smooth   = (pIniFile->ReadInteger("Settings", "smooth", 0));

	//  sprintf(txt, "Analysis using a=%0.2f w=%0.2f normalize=%i",
	//               alpha, weight, normalize);
	//  Form1->Memo1->Lines->Add(txt);
	val_th = StrToFloat(Form1->Edit2->Text);

	show_run_pars_on_form1();

	//delete pIniFile;

	return true;
}

bool Read_input_grid_file(int num, const char *fn, float **tmpm)
{
	int tmpx, tmpy, zfx, zfy, zf;
	String ingridfn;
	float  psum, IG_sum;

	Application->ProcessMessages();

	ingridfn=fn;

	if (miss_as_zero)
		obsmap[num].miss_as_zero=true;
	else
		obsmap[num].miss_as_zero=false;

	if (IG_set.use_IG && IG_set.use_IGw)
	{
		if (!obsmap[num].load_from_file_IG(ingridfn, IG_set,
						   IG_set.IG_spw[num], wmap[num].m, logit, psum, IG_sum, mask_data, area_mask.m, tmpm, 
						   occur_size_weights_map.m))
			return false;
		spp[num].prob_sum=psum;
		spp[num].IG_fract=IG_sum;
	}
	else
	{
		if (!obsmap[num].load_from_file_IG(ingridfn, IG_set,
						   IG_set.IG_spw[num], NULL, logit, psum, IG_sum, mask_data, area_mask.m, tmpm, 
						   occur_size_weights_map.m))
			return false;
		spp[num].prob_sum=psum;
		spp[num].IG_fract=IG_sum;
	}

	//  if (!obsmap[num].load_from_file(ingridfn, Form1->Memo1))
	//    return false; // xxx Z-old


	if (!show_images)
		return true;

	//  obsmap.set_color_for_type(1,50,50,50,1);
	//  obsmap.set_color_for_type(2,150,150,150,1);
#if 0
	Form1->Image1->Picture->Bitmap->Width = obsmap[num].cols();
	Form1->Image1->Picture->Bitmap->Height= obsmap[num].rows();
	Form1->Image1->Width = obsmap[num].cols();
	Form1->Image1->Height= obsmap[num].rows();
	obsmap[num].align_with_image(Form1->Image1->Width, Form1->Image1->Height,tmpx, tmpy);
	// xxx as far as can be seen this is only important if values need to be retrieved from scaled map.
#else
	tmpx = obsmap[num].cols();
	tmpy = obsmap[num].rows();

	Form1->Image1->Picture->Bitmap->Width = obsmap[num].cols();
	Form1->Image1->Picture->Bitmap->Height= obsmap[num].rows();
	Form1->Image1->Width = obsmap[num].cols();
	Form1->Image1->Height= obsmap[num].rows();
	obsmap[num].align_with_image(Form1->Image1->Width, Form1->Image1->Height,tmpx, tmpy);
#endif
	//  obsmap[num].draw_on_bitmap(Form1->Image1->Picture->Bitmap); // xxx2 removed
	Form1->Show();
	return true;
}

// simplified version
bool
Read_input_grid_file_just_to_count(int num, const char *fn, float **tmpm, int** occur_count_matrix, size_t& non_missing)
{
	Application->ProcessMessages();

	if (miss_as_zero)
		obsmap[num].miss_as_zero=true;
	else
		obsmap[num].miss_as_zero=false;

	String ingridfn=fn;
	bool ok = obsmap[num].load_from_file_IG_just_to_count(ingridfn, mask_data, area_mask.m, tmpm, occur_count_matrix);
	non_missing = obsmap[num].non_missing_cnt;
	if (!ok)
	  return false;

	return true;  
}


bool Read_IGw_file(int num, const String& fname, float **wtmpm)
{
	int    tmpx, tmpy, zfx, zfy, zf;
	String ingridfn;
	struct IGWrap {
		struct IG_settings *ptr;
		IGWrap(): ptr(new IG_settings()) {}
		~IGWrap() { delete ptr; }
	} dummy;
	float  weight, psum, IG_sum;

	dummy.ptr->IGa=0.0f;

	ingridfn=fname;
	//Form1->Memo1->Lines->Add("Loading IG file "+ingridfn);
	if (miss_as_zero)
		wmap[num].miss_as_zero=true;
	else
		wmap[num].miss_as_zero=false;

	wmap[num].normalize=false;

	if (!wmap[num].load_from_file_IG(ingridfn, *dummy.ptr,
					 0, NULL, 0, psum, IG_sum, mask_data, area_mask.m, wtmpm)) // error not log BT
		return false;

#if 0
	Form1->Image1->Picture->Bitmap->Width = wmap[num].cols();
	Form1->Image1->Picture->Bitmap->Height= wmap[num].rows();
	Form1->Image1->Width = wmap[num].cols();
	Form1->Image1->Height= wmap[num].rows();
	wmap[num].align_with_image(Form1->Image1->Width, Form1->Image1->Height,tmpx, tmpy);
	// xxx as far as can be seen this is only important if values need to be retrieved from scaled map.
#else
	tmpx = wmap[num].cols();
	tmpy = wmap[num].rows();

	Form1->Image1->Picture->Bitmap->Width = wmap[num].cols();
	Form1->Image1->Picture->Bitmap->Height= wmap[num].rows();
	Form1->Image1->Width = wmap[num].cols();
	Form1->Image1->Height= wmap[num].rows();
	wmap[num].align_with_image(Form1->Image1->Width, Form1->Image1->Height,tmpx, tmpy);
#endif

	return true;
}

void fake_mask_map()
{
	int x,y, masked_in, masked_out;

	masked_in=masked_out=0;
	for(y=0; y<maskmap.rows(); y++)
	{
		for(x=0; x<maskmap.cols(); x++)
		{
			if (x<100)
				maskmap.m[y][x]=0;
			else if (x<200)
				maskmap.m[y][x]=1;
			else if (x<300)
				maskmap.m[y][x]=2;
			else if (x<400)
				maskmap.m[y][x]=10;
			else
				maskmap.m[y][x]=0;
		}
	}
	//	Form1->Memo1->Lines->Add("Count of cells masked in = " +IntToStr(masked_in) );
	//	Form1->Memo1->Lines->Add("Count of cells masked out = "+IntToStr(masked_out));

	Form1->Memo1->Lines->Add("Mask map faked!");
}

void free_occur_matrix(Biodiv_Features_Occur_Container**& m, long nrl, long nrh, long ncl, long nch)
{  
#define NR_END 1
#define FREE_ARG char*
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void
free_vmat()
{
  int x,y;
  if (vmat) {
    for(y=0; y<yd; y++)
      for(x=0; x<xd; x++)
	if(vmat[y][x]) {
	    // free(vmat[y][x]);  // COMPACT_VMAT
	    // vmat[y][x]=0;
	    // this frees 'rowp's - arrays of features
	    vmat[y][x].clear();
	}
#ifdef COMPACT_VMAT  // In practice, this #ifdef could be avoided with free_void, but let's keep type safety
    free_occur_matrix(vmat, 0,yd,0,xd);
#else
    free_fpmatrix(vmat, 0,yd,0,xd);
#endif
    vmat=0;
  }
}

void
free_wrscr_Rsum()
{
  if (Rsum)
    free_matrix(Rsum, 0, yd, 0, xd);
    Rsum = NULL;
}

void
free_wrscr_Rmax()
{
  if (Rmax)
    free_matrix(Rmax, 0, yd, 0, xd);
    Rmax = NULL;
}

void
free_opt_structures()
{
  if (mat_edge_list_pos) {
    free_imatrix(mat_edge_list_pos, 0, yd, 0, xd);
    mat_edge_list_pos = NULL;
  }

  if (exl) {
    delete[] exl;
    exl = NULL;
  }
  if (eyl) {
    delete[] eyl;
    eyl = NULL;
  }
  if (edge) {
    free_cmatrix(edge, 0, yd, 0, xd);
    edge = NULL;
  }

  // Status is still needed at the end when generating output images and rasters
  // if (status) 
  //   free_cmatrix(status, 0, yd, 0, xd);
  // status = NULL;

  free_wrscr_Rsum();
  free_wrscr_Rmax();
}

void free_data()
{
	if (SSIxyPos) free_imatrix(SSIxyPos, 0,yd,0,xd);
	SSIxyPos=0;
	if (SSIxyCnt) free_imatrix(SSIxyCnt, 0,yd,0,xd);
	SSIxyCnt=0;
	if (SSI_list) {
	  delete[] SSI_list;
	  SSI_list = 0;
	}

	if (nwm) free_imatrix(nwm, 0,yd,0,xd);
	nwm=0;
	if (nbm) free_imatrix(nbm, 0,yd,0,xd);
	nbm=0;
	//  if (Rmin) free_matrix(Rmin, 0,yd,0,xd);
	//  Rmin=0;
	if (Rmax) free_matrix(Rmax, 0,yd,0,xd);
	Rmax=0;
	//  if (Rave) free_matrix(Rave, 0,yd,0,xd);
	//  Rave=0;
	if (!forget_Rsum) {
	  if (Rsum) 
	    free_matrix(Rsum, 0,yd,0,xd);
	  Rsum=0;
	}
	if (curves) free_matrix(curves, 0, CURVES_LEN,0,map_cnt+6);  // xxxia
	curves=0;
	if (SSI_indiv_curves) free_matrix(SSI_indiv_curves, 0, CURVES_LEN, 0, SSI_spp_cnt);
	SSI_indiv_curves = 0;
	if (sol) free_matrix(sol, 0,yd,0,xd);
	sol=0;
	if (sol_val) free_matrix(sol_val, 0,yd,0,xd);
	sol_val=0;
	if (cm)  free_imatrix(cm, 0,yd,0,xd);
	cm=0;

	// if (display_mat) 
	//   free_matrix(display_mat, 0,yd,0,xd);
	// display_mat=0;

	free_vmat();

	if (SSI_repr_lvls)
		delete[]SSI_repr_lvls;
	SSI_repr_lvls=0;
	if (repr_lvls)
		delete[]repr_lvls;
	repr_lvls=0;
	if (repr_after_rem)
		delete[]repr_after_rem;
	repr_after_rem=0;

	if (run_mode==2)
		if (load_order)
			delete[] load_order;
	load_order=0;

	if (use_BQP)
	{
	  int loop;
		for(loop=0;loop<radii_cnt;loop++)
		{
			if (nbms[loop])
			{
				free_imatrix(nbms[loop],0,yd,0,xd);
				nbms[loop]=0;
			}
			if (nbms_orig[loop])
			{
				free_matrix(nbms_orig[loop],0,yd,0,xd);
				nbms_orig[loop]=0;
			}
		}
		for(loop=0;loop<map_cnt;loop++)
		{
			if (!spp[loop].BQP_used)
				continue;
			if(spp[loop].BQP_dnbv)
				free(spp[loop].BQP_dnbv);
			spp[loop].BQP_dnbv = 0;

			if(spp[loop].BQP_nbv)
				free(spp[loop].BQP_nbv);
			spp[loop].BQP_nbv = 0;
		}
	}

	if (use_PLULA)
		Free_PLULA_data();

	if (sim_mat)
	{
		free_matrix(sim_mat, 0, comm_set.ct_cnt, 0, comm_set.ct_cnt);
		sim_mat = 0;
	}

	if (comm_sim_mat)
	{
		free_matrix(comm_sim_mat, 0, comm_set.comm_ct_cnt, 0, comm_set.comm_ct_cnt);
		comm_sim_mat = 0;
	}
	// Form1->Memo1->Lines->Add("So far ok.");
	if (ADM_set.use_ADMUs)
	{
		free_ADMU_data();
	}

	if (glbl_groups_info.grpv)
	  delete[] glbl_groups_info.grpv;

	if (cc_accounting) {
	   Form1->Memo1->Lines->Add("NOTE: not deleting cc_accounting.");
	   //delete cc_accounting;
	}

	if (spp)
	delete [] spp;
}

Biodiv_Features_Occur_Container**& 
occur_matrix(long nrl, long nrh, long ncl, long nch)
{
#define NR_END 1

  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  Biodiv_Features_Occur_Container** &m = vmat;
  m =  (Biodiv_Features_Occur_Container **) malloc((size_t)((nrow+NR_END)*sizeof(Biodiv_Features_Occur_Container*)));
  if (!m) Form1->Memo1->Lines->Add("Allocation error 1 in occur_matrix()");
  m += NR_END;
  m -= nrl;
	
  m[nrl]=(Biodiv_Features_Occur_Container *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(Biodiv_Features_Occur_Container)));
  if (!m[nrl]) Form1->Memo1->Lines->Add("Allocation error 2 in occur_matrix()");
  //if (!m[nrl])
  //  Form1->Memo1->Lines->Add("Allocation error 1 in matrix");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

  return m;
}

bool Check_maps_et_alloc()
{
	int loop;
#if 0
	for (loop=0; loop<(map_cnt-1); loop++)
	{
		if ( (obsmap[loop].rows()!=obsmap[loop+1].rows()) ||
				(obsmap[loop].cols()!=obsmap[loop+1].cols())   )
			return false;
	}
#endif

	SSI_repr_lvls  = new double[SSI_spp_cnt]();
	repr_lvls      = new double[map_cnt]();
	repr_after_rem = new double[map_cnt]();

	/*
	xd=tstmap.cols();
	yd=tstmap.rows();
	*/

	//  exl = new int[xd*yd]; // xxx change to max el count!, not grid size!
	//  eyl = new int[xd*yd];  // moved to clac_richness_et_al_mat
	ecnt=0;

	if (use_SSI)  // SSI data only
	{
		SSIxyPos = imatrix(0,yd,0,xd);
		if (!SSIxyPos)
			return false;

		SSIxyCnt = imatrix(0,yd,0,xd);
		if (!SSIxyCnt)
			return false;
		for(int y = 0; y < yd; y++) {
		  for(int x = 0; x < xd; x++) {
		    SSIxyCnt[y][x] = 0;
		  }
		}
	}

	if (!mem_save_mode) // LSIDENT ONLY
	{
	  // In principle these are only needed for LSI, and nwm is used as aux. matrix for BQP
	  if (use_BQP || !PPA_fname.isEmpty()) {
	    nwm = imatrix(0,yd,0,xd);
	    if (!nwm)
	      return false;
	  }
	  if (!PPA_fname.isEmpty()) {
	    cm = imatrix(0,yd,0,xd);
	    if (!cm)
	      return false;
	  }
	}
	else
	{
		Form1->Button7->Enabled=false;
	}

	nbm = imatrix(0,yd,0,xd);
	if (!nbm)
		return false;
	//  Rmin = matrix(0,yd,0,xd);
	//  if (!Rmin)
	//    return false;
	Rmax = matrix(0,yd,0,xd); // Is needed?
	if (!Rmax)
		return false;
	//  Rave = matrix(0,yd,0,xd);
	//  if (!Rave)
	//    return false;

	if (!forget_Rsum && !mem_save_mode) // only for Richness map
	{
		Rsum = matrix(0,yd,0,xd);
		if (!Rsum)
			return false;
	}
	else
	{
		SpecMap->Button3->Enabled=false;
	}

	sol = matrix(0,yd,0,xd);
	if (!sol) // rankmatrix, needed
		return false;

	sol_val = matrix(0,yd,0,xd);
	if (!sol_val) // prop remaining, not really important, but messy to remove
		return false;

	//ShowMessage("IPD 1");
	// if (!mem_save_mode)
	// {
	// 	display_mat = matrix(0,yd,0,xd);
	// 	if (!display_mat) // only needed to show spp maps and UCA maps
	// 		return false;
	// }

	curves = matrix(0, CURVES_LEN, 0, map_cnt+6);
	if (!curves) // needed - and don't take that much space either
		return false;

	SSI_indiv_curves = matrix(0, CURVES_LEN, 0, SSI_spp_cnt);
	if (!SSI_indiv_curves)
		return false;

	edge = cmatrix(0,yd,0,xd);  // partially redundant with status - candidate for removal
	if (!edge)
		return false;
	status = cmatrix(0,yd,0,xd);
	if (!status) // needed
		return false;

#ifdef COMPACT_VMAT
	vmat = occur_matrix(0, yd, 0, xd);
#else
	vmat = fpmatrix(0, yd, 0, xd); // compulsory
#endif
	if (!vmat)
	  return false;
	for(int x=0; x<=xd; x++)
	  for(int y=0;y<=yd;y++) {
	    // vmat[y][x] = 0;  // COMPACT_VMAT
	    // vmat[y][x].clear(); <- would need initialization
	    // which is not done in occur_matrix() - only malloc
	    vmat[y][x].init_constructor();
	  }
	Form1->Edit33->Text=IntToStr(SSI_spp_cnt);
	Form1->Edit39->Text=IntToStr(xd);
	Form1->Edit37->Text=IntToStr(yd);

	return true;
}

// depends on status[][], bat_run.cpp:Initialize_remove() must be called before
void
get_cost_info()
{
  int x, y;
  char txt[255];
  const size_t MAX_MSGS = 10;

  cost_in_area = 0.0f;
  cost_out     = 0.0f;

  Form1->Memo1->Lines->Add("");
  Form1->Memo1->Lines->Add("Analyzing cost layer...");
  size_t bug = 0;
  for(y=0; y<yd; y++) {
    for(x=0; x<xd; x++) {
      if (status[y][x]>0) {
	cost_in_area += costmap.m[y][x]; // xxx cost and initial removal
      } else {
	if (!isnan(costmap.m[y][x]) && costmap.m[y][x] != costmap.nodata()) {
	  cost_out += costmap.m[y][x];
	  if(bug<MAX_MSGS) {
	    Form1->Memo1->Lines->Add("   Note: cost found outside of analysis area, x: " +IntToStr(x)+ ", y: "+IntToStr(y)+", status: " +IntToStr(status[y][x]));
	  } else if(MAX_MSGS==bug) {
	    Form1->Memo1->Lines->Add("   Note: suppressing additional messages about costs outside of study area.");
	  }
	  bug++;
	}
      }
    }
  }
  
  Form1->Memo1->Lines->Add("Total cost in cost layer: "+FloatToStr(cost_in_area+cost_out)+
			   ", cost inside study/analysis area = "+FloatToStr(cost_in_area)+
			   ", cost outside = "+ FloatToStr(cost_out));
  Form1->Memo1->Lines->Add("");
}

// assumes: use_cost==true -> costmap.m is initialized
// - the main loop here could be merged with the main loop in get_cost_info()
// - but it would require a inner if
// bonus: counts # of non-missing cells into ADMU_n_non_missing
void get_per_ADMU_area_and_cost_info_and_non_missing()
{
  Form1->Memo1->Lines->Add("");
  Form1->Memo1->Lines->Add("Checking cost effective cells within administrative units");

  // assume ADMU_cost_in_area and ADMU_cost_used init. to 0.
  for(size_t y=0; y<yd; y++) {
    for(size_t x=0; x<xd; x++) {
      if (status[y][x] > 0) {
	int seq = ADMUs[y][x];
	if (seq <0 || seq >= ADM_set.count)
	  continue;

	ADMU_n_non_missing[seq]++;
	ADMU_cost_in_area[seq] += costmap.m[y][x]; // includes initial removal
	
	float cell_area = 1;
	if (use_occur_size_weights_correct_landscape_fraction)
	  cell_area = get_cell_area_correction_factor(x, y);
	ADMU_area_in[seq] += cell_area;
      }
    }
  }
  for (size_t adu=0; adu < ADM_set.count; adu++) {
    const size_t TXT_LEN = 512;
    char txt[TXT_LEN];
    snprintf(txt, TXT_LEN, "Total cost in admin. unit number %lu (id: %d) = %f,  cost outside = %f  (effective/non-missing cells: %lu", 
	     adu+1, ADM_set.descriptions[adu].id_number, ADMU_cost_in_area[adu], (cost_in_area+cost_out)-ADMU_cost_in_area[adu], ADMU_n_non_missing[adu]);
    Form1->Memo1->Lines->Add(txt);

    snprintf(txt, TXT_LEN, " Area of admin. unit number %d (id: %d) = %f", 
	     adu+1, ADM_set.descriptions[adu].id_number, ADMU_area_in[adu]);
    Form1->Memo1->Lines->Add(txt);

  }
  Form1->Memo1->Lines->Add("");
}

// only # non missing (when per admu curves files are not enabled)
void get_per_ADMU_area_and_non_missing()
{
  Form1->Memo1->Lines->Add("");
  Form1->Memo1->Lines->Add("Checking area and effective cells within administrative units...");

  // assume ADMU_cost_in_area and ADMU_cost_used init. to 0.
  for(size_t y=0; y<yd; y++) {
    for(size_t x=0; x<xd; x++) {
      if (status[y][x] > 0) {
	int seq = ADMUs[y][x];
	if (seq <0 || seq >= ADM_set.count)
	  continue;

	ADMU_n_non_missing[seq]++;
	float cell_area = 1;
	if (use_occur_size_weights_correct_landscape_fraction)
	  cell_area = get_cell_area_correction_factor(x, y);
	ADMU_area_in[seq] += cell_area;
      }
    }
  }

  for (size_t adu=0; adu < ADM_set.count; adu++) {
    const size_t LINELEN = 512;
    char txt[LINELEN];
    cost_rat = cost_out/cost_in_area;
    snprintf(txt, LINELEN, "Admin. unit %lu (id: %d), area: %f, %lu effective cells", 
	     adu, ADM_set.descriptions[adu].id_number, ADMU_area_in[adu], ADMU_n_non_missing[adu]);
    Form1->Memo1->Lines->Add(txt);
  }
  Form1->Memo1->Lines->Add("");
}

void  Check_BQP_nbms()
{
	int  x, y, mat, rad, i, j, cnt;
	int  sx, sy, ex, ey;

	if (!use_BQP)
		return;

	Form1->Memo1->Lines->Add("Checking BQP neighborhoods.");
	for(mat=0;mat<radii_cnt;mat++)
	{
		rad = radii[mat];

		//      cnt3=cnt2=0;
		cnt=0;
		for(y=0; y<yd; y++)
		{
			for(x=0; x<xd; x++)
			{
				cnt+=nbms[mat][y][x];
				//              cnt2+=cnt;
			}
		}
		Form1->Memo1->Lines->Add("Initial average count of neighbors = "+FloatToStr(static_cast<float>(cnt/(xd*yd))));
		Form1->Memo1->Lines->Add("  with radius of "+IntToStr(rad));
		//      Form1->Memo1->Lines->Add("cnt  = "+IntToStr(cnt));
		//      Form1->Memo1->Lines->Add("cnt3 = "+IntToStr(cnt3));
	}
	Form1->Memo1->Lines->Add("Done..");
	Form1->Memo1->Lines->Add("");
}

void Fix_BQP_nbms()
{
	int  x, y, mat, rad, i, j, cnt, cnt2, spnum;// max_cnt;
	int  sx, sy, ex, ey;

	if (!use_BQP)
		return;

	Form1->Memo1->Lines->Add("Fixing BQP neighborhoods");
	for(mat=0;mat<radii_cnt;mat++)
	{
		rad    = radii[mat];
		spnum  = radii_sp[mat];

		cnt2=0;
		for(y=0; y<yd; y++)
		{
			for(x=0; x<xd; x++)
			{
				sx=max(0, x-rad);
				ex=min(x+rad, xd-1);
				sy=max(0, y-rad);
				ey=min(y+rad, yd-1);

				cnt=0;
				//              max_cnt=0;
				for(int j=sy; j<=ey; ++j)
					for(int i=sx; i<=ex; ++i)
					{
						if ((i==x) && (j==y))
							continue;
						//                    max_cnt++;
						if (status[j][i]>0)
						{
							if (BQP_mode==1)
								cnt++;
							else
							{
								if (vmat[j][i][spnum]>=0.0f)
									cnt++;
							}
						}
					}
				//                  else
				//                    cnt3++;

				nbms[mat][y][x] = cnt;
				if (cnt>0)
					nbms_orig[mat][y][x] = (float)cnt; // max_cnt; // /(float)cnt;
				else
					nbms_orig[mat][y][x] = 0.0f; //1.0f; // xxxBQP is this quite right?
				//      if (Rmax[y][x]!=-1)
				//        Form1->Memo1->Lines->Add("orig_rat = "+FloatToStr(nbms_orig[mat][y][x]));

				cnt2+=cnt;
			}
		}
		Form1->Memo1->Lines->Add("Initial average count of neighbors = "+FloatToStr(cnt2/(float)(xd*yd)));
		Form1->Memo1->Lines->Add("   with radius of "+IntToStr(rad));
		//      Form1->Memo1->Lines->Add("cnt  = "+IntToStr(cnt));
		//      Form1->Memo1->Lines->Add("cnt3 = "+IntToStr(cnt3));
	}
	Form1->Memo1->Lines->Add("Done.");
	Form1->Memo1->Lines->Add("");
}

void Fix_nbm()
{
	int x, y;

	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			nbm[y][x] = 8;
		}
	}
}

int Read_IG_file_names()
{
	FILE  *f;
	// the two defined and const could be global...
#define SSCANF_MAX_FILELEN_INT 2048
#define SSCANF_MAX_FILELEN INT_TO_STR(SSCANF_MAX_FILELEN_INT)
	const int MAX_FILELEN = SSCANF_MAX_FILELEN_INT;
	char  line[MAX_FILELEN], fname[MAX_FILELEN];
	int   cnt, pos;
	float weight;

	f=fopen(IG_set.IGwfn.toUtf8().constData(), "r+t");
	if (!f)
		return 0;

	Form1->Memo1->Lines->Add("Loading error rate file list from file "+IG_set.IGwfn);

	cnt=0;
	pos=0;
	while(fgets(line, MAX_FILELEN, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<2)
			continue;
		int ns = sscanf(line, "%f %" SSCANF_MAX_FILELEN "[^\r\n]s", &weight, fname);
		if (2!=ns) {
		  Form1->Memo1->Lines->Add("ERROR: could not parse line in info-gap list file: "+String(line));
		  continue;
		}


		if (strlen(fname)>1)
		{
			if (!resamp_vec[pos])
			{
				pos++;
				continue;
			}
			pos++;
			IG_set.IG_spw[cnt]=weight;
			IG_set.fnames[cnt]=fname;
			cnt++;
		}
	}

	IG_map_cnt=cnt;
	Form1->Memo1->Lines->Add("Read files: IG weight map count = "+IntToStr(IG_map_cnt));

	fclose(f);

	return IG_map_cnt;
}

void init_resamp_vec()
{
  int loop, idx, rsc;
  if ((resample_count>0) && (resample_count<map_cnt)) {
    Form1->Memo1->Lines->Add("Note: using 'resample features/species', sampling " + 
			     IntToStr(resample_count) + " features at random (out of " + IntToStr(map_cnt) + ")");
    for(loop=0; loop<map_cnt; loop++)
      resamp_vec[loop]=0;
    rsc = 0;
    do {
      idx = static_cast<int>(map_cnt*randz());
      if (resamp_vec[idx]==0) {
	resamp_vec[idx]=1;
	rsc++;
      }
    } while(rsc<resample_count);
  } else {
    for(loop=0; loop<map_cnt; loop++)
      resamp_vec[loop]=1;
  }
}

float
get_max_spp_weight()
{
  float wmax = .0f;
  for(size_t i=0; i<map_cnt; i++)
    wmax = std::max(float(fabs(spp[i].weight)), wmax);

  if (map_cnt > 0 && .0f != wmax)
    glob_wmax = wmax;
  return wmax;
}

// Rescales the weights into [0,1] by dividing by wmax
// NOTE: updates *global* glob_wmax, whic is used for "SSI[...].weight = weight/glob_wmax;"
// So after init, calls to rescale_spp_weights should be followed by calls to
// rescale_ssi_spp_weights()
bool
rescale_spp_weights(float wmax)
{
  if (map_cnt <= 0 || .0f == wmax)
    return false;

  float wfactor = wmax / initial_wmax;

  for(int s=0;s<map_cnt;s++) {
    spp[s].weight /= wfactor;
  }

  if (map_cnt > 0 && .0f != wmax)
    glob_wmax = wmax;
  return true;
}

bool
rescale_ssi_spp_weights(float wmax)
{
  if (map_cnt <= 0 || .0f == wmax)
    return false;

  float wfactor = wmax / initial_wmax;

  for (size_t i=0; i<SSI_spp_cnt; i++)
    SSI[i].weight /= wfactor;

  return true;
}

// load spp
int Read_spp_file_names()
{
	using namespace std;
	//FILE  *f;
	//char  line[512], fname[512];
	String fname;
	float weight, alpha, wmax, sr_par;
	int   pnum, prad, pos, rcnt;
	float w2, Tj, exp1, exp2;

	Form1->Memo1->Lines->Add("Loading list of biodiversity features from file: "+sppfn);
	Form1->ComboBox1->Clear();
	Form1->ComboBox2->Clear();

	spp = new struct sp[MAX_SPP_COUNT];

	ZSPPFile sppFile;
	if(!loadFile(loadZSPPFile, sppFile, sppfn, zErrorCallBack)) {
		return 0;
	}
	if(new_spp_file_format == 0) {
		if(sppFile.format == SpeciesFormat::NewFormat) {
			ShowMessage("V1 spp file format - incorrect number of columns (6 expected) around feature/species number "+IntToStr(map_cnt));
			return -1;
		}
	} else {
		if(sppFile.format == SpeciesFormat::OldFormat) {
			ShowMessage("New -v2 spp file format - incorrect number of columns (9 expected) around feature/species number "+IntToStr(map_cnt));
			return -1;
		}
	}

	//f=fopen(sppfn.toUtf8().constData(), "r+t");
	//if (!f)
	//	return 0;

	map_cnt=0;
	wmax=0.0f;
	// weights_sum = 0.0f;  // xxxMCZ  // fedemp: weights_sum unused -> removed
	//while(fgets(line, 512, f))
	foreach(const boost::shared_ptr<ZSPPEntry>& entry, sppFile.sppList)
	{
		// v1

		weight = entry->weight();
		alpha = entry->alpha();
		pnum = entry->bqp_curve_index();
		prad = entry->bqp_buffer_size();
		sr_par = entry->additive_rule_exponent();
		fname = entry->raster();

		// v2

		w2 = entry->generalized_rule_weight();
		Tj = entry->generalized_rule_target();
		exp1 = entry->generalized_rule_exp_x();
		exp2 = entry->generalized_rule_exp_y();

		/*
  if (line[0]=='#')
   continue;
  if (strlen(line)<5)
   continue;
  if (new_spp_file_format==0)
  {
   rcnt=sscanf(line, "%f %f %i %i %f %s", &weight, &alpha, &pnum, &prad, &sr_par, fname);
   if (rcnt!=6)
   {
    ShowMessage("V1 spp file format - incorrect number of columns (6 expected) around species number "+IntToStr(map_cnt));
    fclose(f);
    return -1;
   }
  }
  else
  {
   rcnt=sscanf(line, "%f %f %i %i %f %f %f %f %s", &weight, &alpha, &pnum, &prad, &w2, &Tj, &exp1, &exp2, fname);
   if (rcnt!=9)
   {
    ShowMessage("New -v2 spp file format - incorrect number of columns (9 expected) around species number "+IntToStr(map_cnt));
    fclose(f);
    return -1;
   }
  }
  */

		if (!fname.isEmpty() && new_spp_file_format==0)
		{
			spp[map_cnt].drop  = false;
			spp[map_cnt].fname  = fname;
			spp[map_cnt].weight = weight;
			spp[map_cnt].sr_par = sr_par;
			//if (weight>0.0f)
			//	weights_sum        += weight;  // xxxMCZ
			//else 
			if (weight<0.0)
				neg_weights_used=true;
			spp[map_cnt].alpha  =alpha;
			spp[map_cnt].BQP_prof_num = pnum;
			spp[map_cnt].BQP_radius   = prad;
			spp[map_cnt].T_violation_fract   = -1.0f;
			wmax=max(weight,wmax);
			map_cnt++;
		}
		else if (!fname.isEmpty() && (new_spp_file_format>0))
		{
			spp[map_cnt].drop  = false;
			spp[map_cnt].fname  = fname;
			spp[map_cnt].weight = weight;
			spp[map_cnt].sr_par = Tj;
			//if (weight>0.0f)
			//	weights_sum        += weight;  // xxxMCZ
			//else 
			if (weight<0.0)
				neg_weights_used=true;
			spp[map_cnt].alpha  =alpha;
			spp[map_cnt].BQP_prof_num = pnum;
			spp[map_cnt].BQP_radius   = prad;
			spp[map_cnt].T_violation_fract   = -1.0f;
			spp[map_cnt].w2   = w2;
			spp[map_cnt].Tj   = Tj;
			spp[map_cnt].exp1 = exp1;
			spp[map_cnt].exp2 = exp2;
			wmax=max(weight,wmax);
			map_cnt++;
		}
	}

	// This must be done at least before loading ADMU stuff
	init_resamp_vec();

	//rewind(f);
	map_cnt=0;
	wmax=0.0f;
	//weights_sum = 0.0f;
	pos =0;
	//while(fgets(line, 512, f))
	// Sample features...
	foreach(const boost::shared_ptr<ZSPPEntry>& entry, sppFile.sppList)
	{
		// v1

		weight = entry->weight();
		alpha = entry->alpha();
		pnum = entry->bqp_curve_index();
		prad = entry->bqp_buffer_size();
		sr_par = entry->additive_rule_exponent();
		fname = entry->raster();

		// v2

		w2 = entry->generalized_rule_weight();
		Tj = entry->generalized_rule_target();
		exp1 = entry->generalized_rule_exp_x();
		exp2 = entry->generalized_rule_exp_y();

		/*
  if (line[0]=='#')
   continue;
  if (strlen(line)<5)
   continue;

  if (new_spp_file_format==0)
  {
   rcnt=sscanf(line, "%f %f %i %i %f %s", &weight, &alpha, &pnum, &prad, &sr_par, fname);
   if (rcnt!=6)
   {
    ShowMessage("V1 spp file format - incorrect number of columns (6 expected) around species number "+IntToStr(map_cnt));
    fclose(f);
    return -1;
   }
  }
  else
  {
   rcnt=sscanf(line, "%f %f %i %i %f %f %f %f %s", &weight, &alpha, &pnum, &prad, &w2, &Tj, &exp1, &exp2, fname);
   if (rcnt!=9)
   {
    ShowMessage("New -v2 spp file format - incorrect number of columns (9 expected) around species number "+IntToStr(map_cnt));
    fclose(f);
    return -1;
   }
  }
  */

		if (!fname.isEmpty())
		{
			if (!resamp_vec[pos])
			{
				pos++;
				continue;
			}
			pos++;
			if (!fname.isEmpty() && (new_spp_file_format==0))
			{
				spp[map_cnt].drop  = false;
				spp[map_cnt].fname  = fname;
				spp[map_cnt].weight = weight;
				spp[map_cnt].sr_par = sr_par;
				//if (weight>0.0f)
				//	weights_sum        += weight;  // xxxMCZ
				spp[map_cnt].alpha  =alpha;
				spp[map_cnt].BQP_prof_num = pnum;
				spp[map_cnt].BQP_radius   = prad;
				spp[map_cnt].BQP_prof_down= prad;
				spp[map_cnt].T_violation_fract   = -1.0f;
				wmax=max(fabs(weight),wmax);
				map_cnt++;
				Form1->ComboBox1->Items->Add(fname);
				Form1->ComboBox2->Items->Add(fname);
			}
			else if (!fname.isEmpty() && (new_spp_file_format>0))
			{
				spp[map_cnt].drop  = false;
				spp[map_cnt].fname  = fname;
				spp[map_cnt].weight = weight;
				spp[map_cnt].sr_par = Tj;
				//if (weight>0.0f)
				//	weights_sum        += weight;  // xxxMCZ
				spp[map_cnt].alpha  =alpha;
				spp[map_cnt].BQP_prof_num = pnum;
				spp[map_cnt].BQP_radius   = prad;
				spp[map_cnt].BQP_prof_down= prad;
				spp[map_cnt].T_violation_fract   = -1.0f;
				spp[map_cnt].w2   = w2;
				spp[map_cnt].Tj   = Tj;
				spp[map_cnt].exp1 = exp1;
				spp[map_cnt].exp2 = exp2;
				wmax=max(fabs(weight),wmax);
				map_cnt++;
				Form1->ComboBox1->Items->Add(fname);
				Form1->ComboBox2->Items->Add(fname);
			}
		}
	}

	rescale_spp_weights(wmax);
	// global
	initial_wmax = wmax;

	Form1->ComboBox1->ItemIndex=0;
	//fclose(f);
	return map_cnt;
}

int Count_rows(const char *fname)
{
  int cnt;
  FILE *f;
  char line[MAX_LINE_LEN];
  
  cnt=0;
  f=fopen(fname, "r+t");
  if (!f)
    return -1;
  
  while(fgets(line, MAX_LINE_LEN, f)) {
    if (line[0]=='#')
      continue;
    if (strlen(line)>4)
      cnt++;
  }
  
  fclose(f);
  return cnt;
}

// tokenize line to tokens separated by whitespace
// dataitems will contain tokens
// return number of tokens
int tokenize(char *dataitems[], char *line)
{
	int cnt, linelen, pos;

	linelen = static_cast<int>(strlen(line));

	cnt = pos = 0;
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
			 || (line[pos]==',') || (line[pos]=='.') || (line[pos]=='-'))
				&& (pos<linelen) )  //jumps over alphabets, digits etc.
			pos++;
		line[pos]=0;

		pos++;
	}

	return cnt;
}

int Read_SSI_file_names()
{
	FILE  *f;
	char  line[512], fname[512];
	float weight, alpha, sr_par; // wmax global
	int   pnum, prad, idx, rsc, pos, loop, rc, rcnt;
	float w2, Tj, exp1, exp2;

	f=fopen(SSI_fname.toUtf8().constData(), "r+t");
	if (!f)
		return 0;

	Form1->Memo1->Lines->Add("Loading SSI files list from "+SSI_fname);

	SSI_spp_cnt = 0;
	SSI_row_cnt = 0;
	while(fgets(line, 512, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<5)
			continue;
		if (new_spp_file_format==0)
		{
			rcnt=sscanf(line, "%f %f %i %i %f %s", &weight, &alpha, &pnum, &prad, &sr_par, fname);
			if (rcnt!=6)
			{
				ShowMessage("V1 SSI spp file format used - incorrect number of columns (6 expected) around feature/species number "+IntToStr(SSI_spp_cnt));
				fclose(f);
				return -1;
			}
		}
		else
		{
			rcnt=sscanf(line, "%f %f %i %i %f %f %f %f %s", &weight, &alpha, &pnum, &prad, &w2, &Tj, &exp1, &exp2, fname);
			if (rcnt!=9)
			{
				ShowMessage("New v2 SSI spp file format used - incorrect number of columns (9 expected) around feature/species number "+IntToStr(SSI_spp_cnt));
				fclose(f);
				return -1;
			}
		}

		if (strlen(fname)>1)
		{
			SSI[SSI_spp_cnt].fname  = fname;
			SSI[SSI_spp_cnt].weight = weight/glob_wmax;
			SSI[SSI_spp_cnt].BQP_prof_num = pnum;
			SSI[SSI_spp_cnt].BQP_radius   = prad;
			//          SSI[SSI_spp_cnt].BQP_prof_down= prad;
			SSI[SSI_spp_cnt].lost_at_fraction  = -1.0f;
			SSI[SSI_spp_cnt].T_violation_fract = -1.0f;

			if (new_spp_file_format==0)
			{
				SSI[SSI_spp_cnt].sr_par = sr_par;
			}
			else
			{
				SSI[SSI_spp_cnt].w2   = w2;
				SSI[SSI_spp_cnt].Tj   = Tj;
				SSI[SSI_spp_cnt].exp1 = exp1;
				SSI[SSI_spp_cnt].exp2 = exp2;
				SSI[SSI_spp_cnt].sr_par = Tj;
			}

			rc = Count_rows(fname);
			if (rc<=0)
			{
				ShowMessage("Error with SSI file: rowcount is zero or file not found "+(String)fname);
				Form1->Memo1->Lines->Add("Error with SSI file: rowcount is zero or file not found "+(String)fname);
				continue;
			}
			else
				SSI_row_cnt += rc;

			//          sprintf(line, "%f %f %i %i %f", weight, alpha, prad, pnum, sr_par);
			//          Form1->Memo1->Lines->Add(line);
			//          Form1->Memo1->Lines->Add("SRP"+FloatToStr(SSI[SSI_spp_cnt].sr_par));
			SSI_spp_cnt++;
		}
	}
	//-----

	fclose(f);
	Form1->Memo1->Lines->Add("Total SSI files = "+IntToStr(SSI_spp_cnt));
	Form1->Memo1->Lines->Add("Total SSI observations x features (e.g. species) = "+IntToStr(SSI_row_cnt));
	return SSI_spp_cnt;
}

const static char comment_chars = '#';
std::istream&
get_line_wo_comments(std::istream& f, std::string &line)
{  
  std::istream& is = getline(f, line);
  // cut anything coming after a comment char
  std::string::size_type i = line.find_first_of(comment_chars);
  if ( i != std::string::npos ) {
    line = line.substr(0,i);
  }
  // cut obnoxious DOS trailing \r if it's there...
  // This fixes problems with DOS files when running Roboff on gnu/unixes
  if (line.size() > 0 && '\r' == line[line.size()-1])
    line.erase(line.size()-1);
  return is;
}

bool
load_arbitrary_kernel_table_csv(String& fn, Arb_Kernel& kernel)
{
  if (!arb_kernels_prefix.isEmpty()) {
    fn = arb_kernels_prefix + "/" + fn;
  }

  std::ifstream f(fn.toUtf8().constData(), std::ifstream::in);
  if (!f.is_open()) {
    Form1->Memo1->Lines->Add("ERROR: could not open arbitary kernel file: "+fn);
    return false;
  }

  bool res = true;
  std::string line;
  size_t text_lineno = 1;
  while (get_line_wo_comments(f, line)) {
    if ( comment_chars == line[0] || line.empty()) {
      text_lineno++;
      continue;
    }
    std::string line_format = "%f %f";
    float x, y;
    unsigned int sn = sscanf(line.c_str(), line_format.c_str(), &x, &y);
    if (2 != sn) {
      Form1->Memo1->Lines->Add("ERROR: wrong format in line "+IntToStr(text_lineno)+
			       " of arbitrary kernel file: "+fn+", expecting two numeric values in each line");
      res = false;
    } else {
      if (kernel.first.size() > 1 && x <= kernel.first.back()) {
	Form1->Memo1->Lines->Add("ERROR: loading arbitrary kernel file: "+fn);
	Form1->Memo1->Lines->Add("ERROR: x coordinates are not strictly sorted from lower to higher, "+
				 IntToStr(x)+" is not strictly bigger than "+kernel.first.back());
	return false;	
      }
      kernel.first.push_back(x);
      kernel.second.push_back(y);
    }
  }
  Form1->Memo1->Lines->Add("Succesfully loaded arbitary kernel file with "+IntToStr(kernel.first.size())+" effective rows: "+fn);
  return res;
}

bool
load_arbitrary_kernel(int id)
{
  bool res = false;
  Arbitrary_Kernels_Map::iterator it;
  it = arbitrary_kernels.find(id);
  if (it != arbitrary_kernels.end()) {
    res = true;
    // already loaded
  } else {
    // id is new
    Arb_Kernel* kernel = new Arb_Kernel;
    String fn = "zonation_arbitrary_kernel_"+IntToStr(id)+".txt";

    res = load_arbitrary_kernel_table_csv(fn, *kernel);
    if (!res) {
      Form1->Memo1->Lines->Add("ERROR: could not load arbitary kernel file: "+fn);
      delete kernel;
    } else {
      arbitrary_kernels.insert(std::pair<int, Arb_Kernel>(id, *kernel));
    }
  }
  return res;
}

int load_groups()
{
	FILE  *f;
	const size_t MAX_LINE_LEN = 2048;
	char  line[MAX_LINE_LEN];
	int   row, num, cond_num, ret_num, ret_mode, arb_kernel_ds, arb_kernel_matrix, arb_kernel_ia;

	f=fopen(groups_fname.toUtf8().constData(), "r+t");
	if (!f)
	  return 0;

	Form1->Memo1->Lines->Add("*******************************************************");
	Form1->Memo1->Lines->Add("Loading feature grouping information from "+groups_fname);

	row=0;
	while(fgets(line, MAX_LINE_LEN, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<1)
			continue;
		int cnt=sscanf(line, "%i %i %i %i %i %i %i", &num, &cond_num, &ret_num, &ret_mode, 
			   &arb_kernel_ds, &arb_kernel_matrix, &arb_kernel_ia); //  xxx add retention later
		if (!use_arb_kernels) {
		  arb_kernel_ds = arb_kernel_matrix = arb_kernel_ia = -1;
		} else {
		  // make sure if optional kernels not given 
		  if (cnt < 7)
		    arb_kernel_ia = -1;
		  if (cnt < 6)
		    arb_kernel_matrix = -1;
		}

		if (cnt>4) {
		  if (row<map_cnt) {
		    int reloc_row;
		    if (settings_drop_0_occur_layers) {
		      reloc_row = drop_0_relocation_idx[row];

		      if (drop_0_relocation_idx[row]<0) {
			Form1->Memo1->Lines->Add("ERROR: in groups file, line/feature "+IntToStr(row)+" would become feature "+
						 IntToStr(reloc_row)+" which does not make sense. Something went wrong.");
			break;
		      }
		    } else {
		      reloc_row = row;
		    }

		    spp[reloc_row].grp_op1  = num;
		    spp[reloc_row].grp_cond = cond_num;
		    spp[reloc_row].grp_ret  = ret_num;
		    spp[reloc_row].grp_ret_mode = ret_mode;
		    //spp[reloc_row].grp_LEC = LEC_num;  // deprecated
		    spp[reloc_row].grp_arb_kernel_ds = arb_kernel_ds;
		    spp[reloc_row].grp_arb_kernel_matrix = arb_kernel_matrix;
		    spp[reloc_row].grp_arb_kernel_ia = arb_kernel_ia;
		    
		    row++;
		    
		    if (arb_kernel_ds>=FIRST_NUM_ARB_KERNEL) {
		      if (!load_arbitrary_kernel(arb_kernel_ds)) {
			Form1->Memo1->Lines->Add("ERROR: could not successfully load arbitary kernel for distribution smoothing in groups file ("+groups_fname+
						 ")"+", row "+IntToStr(row)+", column 5, kernel number: "+IntToStr(arb_kernel_ds));
			break;
		      }
		    }
		    if (arb_kernel_matrix>=FIRST_NUM_ARB_KERNEL) {
		      if (!load_arbitrary_kernel(arb_kernel_matrix)) {
			Form1->Memo1->Lines->Add("ERROR: could not successfully load arbitary kernel for matrix transform in groups file ("+groups_fname+
						 ")"+", row "+IntToStr(row)+", column 6, kernel number: "+IntToStr(arb_kernel_matrix));
			break;
		      }
		    }
		    if (arb_kernel_ia>=FIRST_NUM_ARB_KERNEL) {
		      if (!load_arbitrary_kernel(arb_kernel_ia)) {
			Form1->Memo1->Lines->Add("ERROR: could not successfully load arbitary kernel for interaction transform in groups file ("+groups_fname+
						 ")"+", row "+IntToStr(row)+", column 7, kernel number: "+IntToStr(arb_kernel_ia));
			break;
		      }
		    }
		    
		  }
		  else {
		    Form1->Memo1->Lines->Add("Warning: Too many rows found in groups file "+groups_fname);
		    break;
		  }
		}
	}

	fclose(f);

	if (row==map_cnt)
		return 1;
	else
	{
	  Form1->Memo1->Lines->Add("Error: number of rows successfully parsed in groups file ("+IntToStr(row)+
				   ") did not match the number of rows in the features list file ("+IntToStr(map_cnt)+
				   "); in each row 5 columns are expected at least: group, condition, retention, retention_mode and kernel_group "+groups_fname);
		return 0;
	}
}


int Read_BQP_profiles()
{
	FILE  *f;
	size_t TXT_LEN = 512;
	size_t LINE_LEN = 512;
	char  line[LINE_LEN], txt[TXT_LEN];
	int   cnt, inum, tcnt, loop, pos;
	float fnum;
	char  *tokens[100];

	f=fopen(BQP_prof_fn.toUtf8().constData(), "r+t");
	if (!f)
		return 0;

	Form1->Memo1->Lines->Add("Loading BQP profiles from file "+BQP_prof_fn);

	cnt=0;
	while(fgets(line, LINE_LEN, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<5)
			continue;
		tcnt=tokenize(tokens, line);
		if (tcnt<5)
			return 0;
		sscanf(tokens[0], "%i", &inum);
		BQPs[cnt].num=inum;
		pos=0;
		for(loop=1; loop<tcnt; loop+=2)
		{
			sscanf(tokens[loop], "%f", &fnum);
			BQPs[cnt].prop[pos]=fnum;
			sscanf(tokens[loop+1], "%f", &fnum);
			BQPs[cnt].rem[pos]=fnum;
			pos++;
		}
		BQPs[cnt].pnt_cnt=pos;
		snprintf(txt, TXT_LEN, "BQProf %i starts from (%0.2f,%0.2f) ends at (%0.2f,%0.2f)",
			BQPs[cnt].num, BQPs[cnt].prop[0], BQPs[cnt].rem[0],
			BQPs[cnt].prop[BQPs[cnt].pnt_cnt-1], BQPs[cnt].rem[BQPs[cnt].pnt_cnt-1]);
		Form1->Memo1->Lines->Add(txt);
		++cnt;
	}

	prof_cnt = cnt;
	fclose(f);

	for(cnt=0; cnt<prof_cnt; cnt++)
	{
		for(loop=0; loop<=10000; loop++)
		{
			BQPs[cnt].loss10000[loop] = BQP_interpolate(loop, 10000, BQPs[cnt].num);
		}
	}

	return 1;
}

float BQP_interpolate(int remains,  int max_nbc, int prof_num)
{
	float ratio, val;
	int   loop, pnum;
	float denum, diff;

	pnum  = -1;
	ratio = remains/(float)max_nbc;
	//   get_val(prof_num, ratio);
	for(loop=0;loop<prof_cnt;loop++)
		if (BQPs[loop].num==prof_num)
		{
			pnum=loop;
			break;
		}
	if (pnum==-1)
	{
		Form1->Memo1->Lines->Add("************** ERROR ***************");
		Form1->Memo1->Lines->Add("   BQP: wrong profile number error");
		return 0.0f;
	}

	//   Form1->Memo1->Lines->Add("prof num "+IntToStr(prof_num)+" pnum= "+IntToStr(pnum));
	for(loop=0;loop<(BQPs[pnum].pnt_cnt-1);loop++)
	{
		if ((BQPs[pnum].prop[loop]>=ratio) && (BQPs[pnum].prop[loop+1]<=ratio))
		{
			diff  = BQPs[pnum].prop[loop]-ratio;
			denum = BQPs[pnum].prop[loop]-BQPs[pnum].prop[loop+1];
			if (denum==0.0f)
			{
				Form1->Memo1->Lines->Add("************** ERROR ***************");
				Form1->Memo1->Lines->Add("pnum interp error  pntc="
							 +IntToStr(BQPs[pnum].pnt_cnt)+
							 " loop= "+IntToStr(loop));
				return 0.0f;
			}

			val = BQPs[pnum].rem[loop] -
					diff*(BQPs[pnum].rem[loop]-BQPs[pnum].rem[loop+1])/denum;
			return val;
		}
	}
	Form1->Memo1->Lines->Add("BQP: pnum interpolation error.");

	return 0.0f;
}

bool Fix_nbvms()
{
	int loop, s, max_nbc, alloc_cnt, pnum;
	float v1, v2;

	alloc_cnt=0;
	for(s=0;s<map_cnt;s++)
	{
		if (!spp[s].BQP_used)
			continue;
		max_nbc = (2*spp[s].BQP_radius+1)*(2*spp[s].BQP_radius+1);
		spp[s].BQP_dnbv =
				(float *)calloc(1+max_nbc, sizeof(float));
		if (!spp[s].BQP_dnbv)
			return false;
		spp[s].BQP_nbv =
				(float *)calloc(1+max_nbc, sizeof(float));
		if (!spp[s].BQP_nbv)
			return false;
#if 0
		if (use_PLULA && use_tree_conn)
		{
			spp[s].BQP_link2 = 0;
			for(loop=0;loop<prof_cnt;loop++)
				if (BQPs[loop].num==spp[s].BQP_prof_down)
				{
					spp[s].BQP_link2 =  BQPs[loop].loss10000; // link2 down
					continue;
				}
			if (!spp[s].BQP_link2)
				ShowMessage("Error linking downriver BQP response");
		}
#endif
		for(loop=0;loop<prof_cnt;loop++)
			if (BQPs[loop].num==spp[s].BQP_prof_num)
			{
				spp[s].BQP_link =  BQPs[loop].loss10000; // link ordinary BQP or tree up
				break;
			}

		for(loop=1;loop<=max_nbc;loop++)
		{
			v1 = BQP_interpolate(loop,   max_nbc, spp[s].BQP_prof_num);
			v2 = BQP_interpolate(loop-1, max_nbc, spp[s].BQP_prof_num);
			spp[s].BQP_dnbv[loop] = v1-v2;
			spp[s].BQP_nbv[loop]  = v1;
		}
		spp[s].BQP_dnbv[0] = 0;
		spp[s].BQP_nbv[0]  = BQP_interpolate(0, max_nbc, spp[s].BQP_prof_num);
		alloc_cnt++;
	}

	Form1->Memo1->Lines->Add(IntToStr(alloc_cnt)+" BQP profiles computed.");
	return true;
}

bool Fix_tree_profiles()
{
	int loop, s;

	for(s=0;s<map_cnt;s++)
	{
		spp[s].TREE_used=true;
		spp[s].BQP_link2 = 0;
		for(loop=0;loop<prof_cnt;loop++)
			if (BQPs[loop].num==spp[s].BQP_prof_down)
			{
				spp[s].BQP_link2 =  BQPs[loop].loss10000; // link2 down
				continue;
			}
		if (!spp[s].BQP_link2)
			ShowMessage("Error linking downriver BQP response");

		if (!spp[s].BQP_link2)
		{
			ShowMessage("Error linking downriver BQP response");
			spp[s].TREE_used=false;
		}

		spp[s].BQP_link = 0;
		for(loop=0;loop<prof_cnt;loop++)
			if (BQPs[loop].num==spp[s].BQP_prof_num)
			{
				spp[s].BQP_link =  BQPs[loop].loss10000;
				break;
			}

		if (!spp[s].BQP_link)
		{
			ShowMessage("Error linking upriver BQP response");
			spp[s].TREE_used=false;
		}

		if (!spp[s].TREE_used)
			Form1->Memo1->Lines->Add("Feature "+IntToStr(s+1)+" no tree connectivity linked.");
	}

	Form1->Memo1->Lines->Add(" Tree connectivity BQP profiles linked");
	return true;
}
bool Alloc_nbms()
{
	int  loop, nbm_idx, loop2;
	bool found;

	nbm_idx=0;
	for(loop=0;loop<map_cnt;loop++)
	{
		if (spp[loop].BQP_radius<=0)
		{
			spp[loop].BQP_used=false;
			spp[loop].nbmat_num=-1;
			continue;
		}
		spp[loop].BQP_used=true;

		found=false;
		for(loop2=0; loop2<nbm_idx; loop2++)
		{
			//         Form1->Memo1->Lines->Add(IntToStr(radii[loop2])+"  "+IntToStr(spp[loop].BQP_radius));
			if ((BQP_mode==1) && (radii[loop2]==spp[loop].BQP_radius))
			{ // xxx BQP_mode fix
				spp[loop].nbmat_num=loop2;
				nbm_num[loop]=loop2;
				found=true;
				break;
			}
		}
		if (found)
			continue;
		radii[nbm_idx]      = spp[loop].BQP_radius;
		radii_sp[nbm_idx]   = loop;
		spp[loop].nbmat_num = nbm_idx;
		nbm_num[loop]       = nbm_idx;
		nbms[nbm_idx]       = imatrix(0,yd,0,xd);
		if (!nbms[nbm_idx])
			return false;
		nbms_orig[nbm_idx]  = matrix(0,yd,0,xd);
		if (!nbms_orig[nbm_idx])
			return false;
		Form1->Memo1->Lines->Add("Neighborhood matrix allocated with radius "
					 +IntToStr(radii[nbm_idx]));
		nbm_idx++;
	}

	radii_cnt=nbm_idx;
	Fix_nbvms();
	return true;
}

// a <- a*b   (only for the positive a[][]
bool
matrix_mult_positive(float** a, float** b, int xd, int yd)
{
  for (size_t j=0; j<yd; j++) {
    for (size_t i=0; i<xd; i++) {
      if (a[j][i] >= 0)
	a[j][i] = a[j][i]*b[j][i];
    }
  }

  return true;
}

bool
check_cost_map_ok()
{
  for (int j=0; j<yd; j++) {
    for (int i=0; i<xd; i++) {
      if (status[j][i] >0 && costmap.m[j][i] <= .0f)  {
	Form1->Memo1->Lines->Add("ERROR in costs map, the cost value in coordinates x: " + IntToStr(i)+
				 ", y: " + IntToStr(j) + " is: " + FloatToStr(costmap.m[j][i]) +
				 ", but values must be strictly positive!");
	return false;
      }
    }
  }
  return true;
}

bool
load_rest_of_special_files()
{
	int  loop, x, y;
	const size_t TXT_LEN = 512;
	char txt[TXT_LEN];

	if (!glob_set_occur_size_weights_layer_fn.isEmpty()) {
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("****** Using map of occurrence weights / relative cell areas: " + glob_set_occur_size_weights_layer_fn+
				   " ******");
	  occur_size_weights_map.normalize = false;
	  // without area mask, cell size factor may be needed anywhere (e.g. at distribution centers, not necessarily inside the area mask)
	  if (!occur_size_weights_map.load_from_file(glob_set_occur_size_weights_layer_fn, 0 /*mask_data*/, NULL /*area_mask.m*/)) {
	    Form1->Memo1->Lines->Add("ERROR while loading map of occurence weights:" + glob_set_occur_size_weights_layer_fn);
	    return false;
	  } else {
	    Form1->Memo1->Lines->Add("Successfully loaded map of occurrence weights: " + glob_set_occur_size_weights_layer_fn);
	  }
	}

	if (use_cost) {
	  costmap.normalize = false;
	  if (!costmap.load_from_file(costfn, mask_data, area_mask.m)) {
	    Form1->Memo1->Lines->Add("************** ERROR ***************");
	    Form1->Memo1->Lines->Add("FAILURE when loading cost map: "+costfn);
	    return false;
	  }
	  Form1->Memo1->Lines->Add("Cost map loaded.");
	  
	  if (use_occur_size_weights_correct_cost) {
	    Form1->Memo1->Lines->Add("Correcting cost layers with relative cell areas/weights from map:"+
				     glob_set_occur_size_weights_layer_fn);
	    matrix_mult_positive(costmap.m, occur_size_weights_map.m, xd, yd);

	    if (!use_occur_size_weights_correct_landscape_fraction)
	      occur_size_weights_map.free_matrix_m();
	  }
	} else if (!glob_set_occur_size_weights_layer_fn.isEmpty()) {
	  // Cost from cell (area) weights
	  use_cost = true;
	  String w_layer_name = "cell areas/weights layer";
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("****** Note: not using standard cost layer (as in 'use cost' and 'cost file') but assigning costs from "+
				   w_layer_name+": "+ glob_set_occur_size_weights_layer_fn + " ******");
	  costmap.normalize = false;
	  // without area mask, cell size factor may be needed anywhere (e.g. at distribution centers, not necessarily inside the area mask)
	  if (!costmap.load_from_file(glob_set_occur_size_weights_layer_fn, 0 /*mask_data*/, NULL /*area_mask.m*/)) {
	    Form1->Memo1->Lines->Add("************** ERROR ***************");
	    Form1->Memo1->Lines->Add("FAILURE when loading cost map: "+glob_set_occur_size_weights_layer_fn);
	    return false;
	  }
	  Form1->Memo1->Lines->Add("Cost map loaded from "+w_layer_name+".");	  
	}

	if (use_cost && !check_cost_map_ok())
	  return false;

	if (use_mask)
	{
		maskmap.normalize=false;
		if (!maskmap.load_from_file(maskfn, mask_data, area_mask.m))
		{
			Form1->Memo1->Lines->Add("************** ERROR ***************");
			Form1->Memo1->Lines->Add("  FAILURE while trying to load removal mask map.");
			return false;
		}
		Form1->Memo1->Lines->Add("Mask map loaded.");
		//		fake_mask_map();
	}

	if (use_groups)
	{
		if (!load_groups())
		{
			Form1->Memo1->Lines->Add("************** ERROR ***************");
			Form1->Memo1->Lines->Add("  FAILURE attempting load of groups file "+groups_fname);
			return false;
		}
		Form1->Memo1->Lines->Add("Groups information loaded.");
	}

	if (use_BQP || use_tree_conn)
	{
		if (!Read_BQP_profiles())
		{
			MessageBox(0, "Could not read BQP/tree_conn profiles file.", "Data input error", MB_OK);
			return false;
		}
		if (use_tree_conn)
			Fix_tree_profiles();
	}

	if (use_BQP)
	{
		Form1->Memo1->Lines->Add("");
		for(loop=0;loop<map_cnt;loop++)	{
		  snprintf(txt, TXT_LEN, "USING BQP: SP %i, BQPProf=%i, radius=%i cells.",
			   loop, spp[loop].BQP_prof_num, spp[loop].BQP_radius);
		  Form1->Memo1->Lines->Add(txt);
		}

		if (!Alloc_nbms())
			return false;
	}

	if (run_mode==2) {
	  if (!loaded1.load_from_file(load1fn, mask_data, area_mask.m))
	    return false;
	}
	return true;
}

bool  Read_SSI_file(int num)
{
	FILE  *f;
	const int MAX_LINE_LEN = 1024;
	char  line[MAX_LINE_LEN];
	float abundance;
	int   x, y, cnt, prev_lpos, loop;
	float sum, sd, val;
	double fpx, fpy;

	f=fopen(SSI[num].fname.toUtf8().constData(), "r+t");
	if (!f)
	{
		ShowMessage("Could not open "+SSI[num].fname);
		return 0;
	}

	prev_lpos = SSI_lpos;
	while(fgets(line, MAX_LINE_LEN, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<5)
			continue;
		cnt=sscanf(line, "%lf %lf %f %f", &fpx, &fpy, &abundance, &sd);
		if (cnt==3)
			sd=0.0f;
		if (cnt<3) // xxxSSI IG, coord transform
		{
			Form1->Memo1->Lines->Add("Parameter count error in file "+SSI[num].fname);
			continue;
		}

		//      x = (int)fpx;
		//      y = (int)fpy;

		x = (int)(floor((fpx-obsmap[0].dxc)/obsmap[0].dcs));
		y = (int)(floor((fpy-obsmap[0].dyc)/obsmap[0].dcs));
		y = yd-y-1;
		x = x;
		SSI_list[SSI_lpos].spnum = num;
		SSI_list[SSI_lpos].x     = x;
		SSI_list[SSI_lpos].y     = y;
		if ((x<0) || (y<0) || (x>xd) || (y>yd))
		{
			Form1->Memo1->Lines->Add("SSI location out of features/species grid range in "+SSI[num].fname);
			snprintf(line, MAX_LINE_LEN, "  %lf %lf -> x=%i y=%i GRID:(gxc=%lf gyx=%lf cs=%lf)",
				fpx, fpy, x, y, obsmap[0].dxc, obsmap[0].dyc, obsmap[0].dcs);
			Form1->Memo1->Lines->Add(line);
			//          ShowMessage("blaah");
			SSI_err_cnt++;
			continue;
		}
		if (vmat[y][x]==0)  // checks if x,y is legal postion
		{
			Form1->Memo1->Lines->Add("SSI location in missing data area in SSI file "+SSI[num].fname);
			snprintf(line, MAX_LINE_LEN, "  %lf %lf -> %i %i GRID:(gxc=%lf gyx=%lf cs=%lf)",
				fpx, fpy, x, y, obsmap[0].dxc, obsmap[0].dyc, obsmap[0].dcs);
			Form1->Memo1->Lines->Add(line);
			SSI_err_cnt++;
			continue;
		}

		val = abundance;

		if ((IG_set.use_IG) && (IG_set.IGa>0.0f))
		{
			if (IG_set.IG_proportional)
				val /= (1.0f+IG_set.IGa*sd);
			else
				val -= (IG_set.IGa*sd);
			if (val<0.0f)
				val=0.0f;
		}
		SSI_list[SSI_lpos].val   = val;
		SSI[num].prob_sum += val;
		//      if (vmat[y][x]!=0)  // checks if x,y is legal postion
		SSI_lpos++;
	}

	fclose(f);

	snprintf(line, MAX_LINE_LEN, "SSI %s count of successfully read locations = %i, parameter=%0.3f",
		SSI[num].fname.toUtf8().constData() , SSI_lpos-prev_lpos, SSI[num].sr_par);
	Form1->Memo1->Lines->Add(line);

	sum = 0.0f;
	for(loop=prev_lpos; loop<SSI_lpos; loop++)
		if (SSI_list[loop].val>0.0f)
			sum+=SSI_list[loop].val;

	if (sum>=0.0f)
	{
		for(loop=prev_lpos; loop<SSI_lpos; loop++)
			SSI_list[loop].val /= sum;
	}
	else
	{
		SSI[num].weight=0.0f;
		Form1->Memo1->Lines->Add("Abundance sum after InfoGap is zero for SSI spp number "+IntToStr(num));
		Form1->Memo1->Lines->Add("   => complete loss of information, feature/species weight set to zero.");
	}

	return true;
}

int Load_interactions(String*& transf_layers_str_append)
{
	FILE  *f;
	char  line[512];
	float alpha, ia_par1;
	int   L1, L2, cnt, iatype;
	//  int   pnum, prad, idx, rsc, pos, loop;

	f=fopen(ia_fname.toUtf8().constData(), "r+t");
	if (!f)
		return 0;

	Form1->Memo1->Lines->Add("Loading feature/species interaction specifications from file "+ia_fname);

	ia_cnt=0;
	while(fgets(line, 512, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<5)
			continue;

		cnt=sscanf(line, "%i %i %f %i %f", &L1, &L2, &alpha, &iatype, &ia_par1);
		if (cnt==5)
		{
			if ((L1<=map_cnt) && (L2<=map_cnt) && (L1>0) && (L2>0) && (alpha>0))
			{
				Do_smoothing_interact(L1, L2, alpha, iatype, ia_par1);
				ia_cnt++;
				
				if (transf_layers_set.some_output) {
				  transf_layers_str_append[L1-1] += "_" + IntToStr(L2) + "_IA";
				}
			}
		}
	}
	
	fclose(f);
	return ia_cnt;

#if 0
	float  prob_sum;
	float  max_local_val;
	float  sr_par;
	int    nbmat_num;
#endif
}

int SSI_sortfunc( const void *a, const void *b)
{
	struct SSI_xydata *s1, *s2;
	s1 = (struct SSI_xydata *)a;
	s2 = (struct SSI_xydata *)b;

	if (s1->y<s2->y)
		return -1;
	else if (s1->y>s2->y)
		return 1;
	else
	{
		if (s1->x<s2->x)
			return -1;
		else if (s1->x>s2->x)
			return 1;
		else
			return 0;
	}
}

bool SSI_load()
{
	int  loop, prev_x, prev_y, prev_pos, cnt, x, y;
	bool ok;
	const size_t TXT_LEN = 512;
	char txt[TXT_LEN];

	SSI_list = new struct SSI_xydata[SSI_row_cnt];

	SSI_lpos    = 0;
	SSI_err_cnt = 0;

	for(loop=0;loop<SSI_spp_cnt;loop++)
	{
		ok=Read_SSI_file(loop);

		if (!ok)
		{
			ShowMessage("Error reading SSI features/species "+SSI[loop].fname);
			Form1->Memo1->Lines->Add("Error reading SSI features/species "+SSI[loop].fname);
			//          return false;
		}
	}

	Form1->Memo1->Lines->Add("* Completed load of SSI files;");
	Form1->Memo1->Lines->Add("* Count of observations successfully read = "+IntToStr(SSI_lpos));
	Form1->Memo1->Lines->Add("* Count of SSI observations with errors in coordinates = "+IntToStr(SSI_err_cnt));

	qsort((void *)SSI_list, SSI_lpos, sizeof(struct SSI_xydata), SSI_sortfunc);

	prev_x = prev_y = prev_pos = -1;
	cnt = 0;
	for(loop=0; loop<SSI_lpos; loop++)
	{
		x = SSI_list[loop].x;
		y = SSI_list[loop].y;
		if ((y!=prev_y) || (x!=prev_x))
		{
			if (prev_x!=-1)
			{
				SSIxyPos[prev_y][prev_x]= prev_pos;
				SSIxyCnt[prev_y][prev_x]= cnt;
				//              sprintf(txt, "assigned location %i %i, starting from %i, %i elements",
				//                prev_y, prev_x, prev_pos, cnt);
				//              Form1->Memo1->Lines->Add(txt);
			}
			prev_y = y;
			prev_x = x;
			prev_pos = loop;
			cnt = 0;
		}
		cnt++;
	}
	if (prev_x!=-1)
	{
		SSIxyPos[prev_y][prev_x]= prev_pos;
		SSIxyCnt[prev_y][prev_x]= cnt;
	}
	//  sprintf(txt, "assigned location %i %i, starting from %i, %i elements",
	//    prev_y, prev_x, prev_pos, cnt);
	//  Form1->Memo1->Lines->Add(txt);

#if 1
	// If SSI_lpos is big printing all the points becomes annoying and even a source of errors
	const int num_head_ssi = 10;
	const int num_tail_ssi = 10;

	if (SSI_lpos >= 20) {
	  Form1->Memo1->Lines->Add("First 10 SSI observations: x y feature_number occurrence_value");
	  for(loop=0; loop<num_head_ssi; loop++) {
	    snprintf(txt, TXT_LEN, "%i %i %i %f", SSI_list[loop].x,SSI_list[loop].y,
		    SSI_list[loop].spnum, SSI_list[loop].val);
	    Form1->Memo1->Lines->Add(txt);
	  }
	  Form1->Memo1->Lines->Add("Last 10 SSI observations: x y feature_number occurrence_value");
	  for(loop=SSI_lpos-num_tail_ssi; loop<SSI_lpos; loop++) {
	    snprintf(txt, TXT_LEN, "%i %i %i %f", SSI_list[loop].x,SSI_list[loop].y,
		    SSI_list[loop].spnum, SSI_list[loop].val);
	    Form1->Memo1->Lines->Add(txt);
	  }
	} else {
	  Form1->Memo1->Lines->Add("SSI observations read: x y feature_number occurrence_value");
	  for(loop=0; loop<SSI_lpos; loop++) {
	    snprintf(txt, TXT_LEN, "%i %i %i %f", SSI_list[loop].x,SSI_list[loop].y,
		    SSI_list[loop].spnum, SSI_list[loop].val);
	    Form1->Memo1->Lines->Add(txt);
	  }
	}
	Form1->Memo1->Lines->Add("");
#endif

	return true;
}

//-----------------------------------------------------------
int Load_and_check_condition_matrix(const String &cfname)
{
	int x, y, error_cnt, de_cnt, missing_cnt;

	condition.normalize = false;
	if (!condition.load_from_file(cfname, mask_data, area_mask.m))
	{
		Form1->Memo1->Lines->Add("************** ERROR ***************");
		Form1->Memo1->Lines->Add("FAILURE when loading condition map ");
		Form1->Memo1->Lines->Add(cfname);
		return 0;
	}

	error_cnt=0;
	de_cnt=0;
	missing_cnt=0;
	for(y=0;y<yd;y++)
		for(x=0;x<xd;x++)
		{
			if (condition.m[y][x]>1.0f)
			{
				condition.m[y][x] = 1.0f;
				error_cnt++;
			}
			else if (condition.m[y][x]>=0.0f)
				de_cnt++;
			else
				missing_cnt++;
		}

	if (error_cnt==0) {
	  Form1->Memo1->Lines->Add("Condition layer loaded without errors, count of effective elements (with values in [0,1]) = "
				   +IntToStr(de_cnt)+" , missing (<0) elements count = "+IntToStr(missing_cnt));
	} else {
	  Form1->Memo1->Lines->Add("Condition layer loaded with errors (condition values > 1) count = "+IntToStr(error_cnt));
	  return 0;
	}

	return 1;
}

//-----------------------------------------------------------
int Load_and_check_retention_matrix(const String &cfname)
{
	int x, y, error_cnt, de_cnt, missing_cnt;

	condition.normalize = false; // condition map also used for retention
	if (!condition.load_from_file(cfname, mask_data, area_mask.m))
	{
		Form1->Memo1->Lines->Add("************** ERROR ***************");
		Form1->Memo1->Lines->Add("FAILURE when loading retention map ");
		Form1->Memo1->Lines->Add(cfname);
		return 0;
	}

	error_cnt=0;
	de_cnt=0;
	missing_cnt=0;
	for(y=0;y<yd;y++)
		for(x=0;x<xd;x++)
		{
#if 0
			if (condition.m[y][x]>1.0)
			{
				condition.m[y][x] = 0.0; // loss = (1-x) from retention
				error_cnt++;
			}
			else if (condition.m[y][x]>=0.0)
			{
				condition.m[y][x] = 1.0-condition.m[y][x];
				de_cnt++;
			}
			else
				missing_cnt++;
#endif
			condition.m[y][x] = fabs(1.0-condition.m[y][x]);
		}

	if (error_cnt==0) {
	  Form1->Memo1->Lines->Add("Retention layer loaded without errors, count of effective elements (with values in [0,1]) = "
				   +IntToStr(de_cnt)+" , missing (<0) elements count = "+IntToStr(missing_cnt));
	} else {
	  Form1->Memo1->Lines->Add("Retention layer loaded with errors (retention values > 1) count = "+IntToStr(error_cnt));
	  return 0;
	}

	return 1;
}

int Load_and_check_LEC_matrix(const String &cfname)
{
	int x, y, error_cnt, de_cnt, missing_cnt;

	condition.normalize = false; // condition map also used for LEC
	if (!condition.load_from_file(cfname, mask_data, area_mask.m))
	{
		Form1->Memo1->Lines->Add("************** ERROR ***************");
		Form1->Memo1->Lines->Add("FAILURE when loading LEC map ");
		Form1->Memo1->Lines->Add(cfname);
		return 0;
	}

	error_cnt=0;
	de_cnt=0;
	missing_cnt=0;
	for(y=0;y<yd;y++)
		for(x=0;x<xd;x++)
		{
			if (condition.m[y][x]>=1.0)
			{
				error_cnt++;
			}
			else if (condition.m[y][x]>=0.0)
			{
				de_cnt++;
			}
			else
				missing_cnt++;
		}

	if (error_cnt==0)
		Form1->Memo1->Lines->Add("LEC layer loaded without errors, count of informative elements= "
					 +IntToStr(de_cnt)+" , missing (<0) elements count = "+IntToStr(missing_cnt));
	else
		Form1->Memo1->Lines->Add("LEC layer loaded with error (edge>1) count = "+IntToStr(error_cnt));

	return 1;
}

//------------------------------------------------------------
float Apply_condition_matrix(int cond_num, int cond_mode)
{
	int   x, y, s;
	float lost[MAX_SPP_COUNT], loss=0.0, val, remains[MAX_SPP_COUNT];
	const size_t TXT_LEN = 512;
	char  txt[TXT_LEN];

	for(s=0;s<map_cnt;s++)
	{
		lost[s]=0.0;
		remains[s]=0.0;
	}

	for(y=0; y<yd; y++)
	{
		for(x=0; x<xd; x++)
		{
			if (!vmat[y][x])
				continue;

			//for(s=0;s<map_cnt;s++)
			const Biodiv_Features_Occur_Container& rowp = vmat[y][x];
			for(s=rowp.first(); s != rowp.overflow(); s=rowp.next(s))
			{
				if (cond_mode)
				{
					if (spp[s].grp_cond!=cond_num)
						continue;
				}
				else
					if (spp[s].grp_ret!=cond_num)
						continue;

				val = vmat[y][x][s];
				if ((val==-1) || (condition.m[y][x]<0.0))
					continue;

				if (cond_mode)
				{
					loss = (1.0-condition.m[y][x])*val;
					vmat[y][x][s] -= loss;
					lost[s]+=loss;
				}
				else // ret_mode
				{
					lost[s]    += condition.m[y][x]*val; // condition actually has 1-Ret
					remains[s] += val;
					vmat[y][x][s] = condition.m[y][x]*val; // loss or addition
				}
			}
		}
	}

	if (cond_mode)
	{
	  Form1->Memo1->Lines->Add(" --- Report on loss of biodiversity feature proportions (by condition transformation) ---");
	  Form1->Memo1->Lines->Add("Feature#  proportion_lost");
	  int effective_cond_cnt = 0;
	  for(s=0;s<map_cnt;s++) {

	    if (spp[s].grp_cond!=cond_num)
	      continue;

	    effective_cond_cnt++;
	    
	    if ((lost[s]>0.0) && (0))
	      spp[s].weight *= 1.0/lost[s];

	    snprintf(txt, TXT_LEN, "%i\t%0.5f", s+1, lost[s]);
	    Form1->Memo1->Lines->Add(txt);
	  }
	  Form1->Memo1->Lines->Add("-----------------------");
	  if (0==effective_cond_cnt)
	    Form1->Memo1->Lines->Add("WARNING: The condition transformation with condition layer number "+IntToStr(cond_num)+" has not been applied to any any biodiversity feature. This looks suspicious. Please check settings.");
	  else
	    Form1->Memo1->Lines->Add("The condition transformation with condition layer number "+IntToStr(cond_num)+" was applied to "+IntToStr(effective_cond_cnt)+" biodiversity features.");

	} else {// retention
                int effective_ret_cnt = 0;
		for(s=0;s<map_cnt;s++) {
			if (spp[s].grp_ret!=cond_num)
				continue;

			effective_ret_cnt++;
			if (remains[s]>0.0)
			{
				if (spp[s].grp_ret_mode==1) // going down
				{
					spp[s].weight *= ret_layers_rel_weight*(lost[s]/remains[s]);

					Form1->Memo1->Lines->Add("Feature#  0=loss   remains(after cond)   diff=(1-ret)*remains   Wj*rlrw*diff/remains");
					snprintf(txt, TXT_LEN, "%i\t%i\t%0.5f\t%0.5f\t%0.4f", s+1, spp[s].grp_ret_mode, remains[s], lost[s], spp[s].weight);
					Form1->Memo1->Lines->Add(txt);
				}
				else // mgmt gain
				{
					spp[s].weight *= ret_layers_rel_weight*(lost[s]/(remains[s]+lost[s]));
					Form1->Memo1->Lines->Add("Feature#  1=gain   remains(after cond)   diff=(1-ret)*remains   Wj*rlrw*diff/(remains+diff)");
					snprintf(txt, TXT_LEN, "%i\t%i\t%0.5f\t%0.5f\t%0.4f", s+1, spp[s].grp_ret_mode, remains[s], lost[s], spp[s].weight);
					Form1->Memo1->Lines->Add(txt);
				}
			}
			else
			{
				Form1->Memo1->Lines->Add("Problem: 100% Retention makes loss zero for feature number = "+IntToStr(s));
			}
		}

		Form1->Memo1->Lines->Add("-----------------------");
		if (0==effective_ret_cnt)
		  Form1->Memo1->Lines->Add("WARNING: The retention transformation with retention layer number "+IntToStr(cond_num)+" has not had any effect on any biodiversity feature. This looks suspicious. Please check settings.");
		else
		  Form1->Memo1->Lines->Add("The retention transformation with retention layer number "+IntToStr(cond_num)+" has been applied to "+IntToStr(effective_ret_cnt)+" biodiversity features.");
		
	}
	//	for(s=0;s<map_cnt;s++)
	//		lost[s]=0.0;

	return 0.0;
}

//------------------------------------------------------------
float Apply_LEC_matrix(int LEC_num)
{
	int   x, y, s;
	char  txt[255];
	float val;
	double element_sum, pre_sum;

	for(s=0;s<map_cnt;s++)
	{
		element_sum=0.0;
		// grp_LEC deprecated
		/*
		if (spp[s].grp_LEC!=LEC_num)
			continue;
		*/

		pre_sum=0.0;
		for(y=0; y<yd; y++)
		{
			for(x=0; x<xd; x++)
			{
				if (!vmat[y][x])
					continue;
				if (vmat[y][x][s]==-1)
					continue;

				pre_sum += vmat[y][x][s];
			}
		}

		for(y=0; y<yd; y++)
		{
			for(x=0; x<xd; x++)
			{
				if (!vmat[y][x])
					continue;

				val = vmat[y][x][s];
				if ((val==-1) || (condition.m[y][x]<0.0) || (condition.m[y][x]>=1.0)) // no edge
					continue;

				vmat[y][x][s] /= (1-condition.m[y][x]);
				element_sum += vmat[y][x][s];
			}
		}

		for(y=0; y<yd; y++)
		{
			for(x=0; x<xd; x++)
			{
				if (!vmat[y][x])
					continue;
				if (vmat[y][x][s]==-1)
					continue;

				vmat[y][x][s] /= element_sum;
				vmat[y][x][s] *= pre_sum;   // make the layer sum to what it was before the edge operation (does not influence)
				// outcome of IG etc.
			}
		}
	}

	return 0.0f;
}

bool Load_and_check_condition_matrixes(int cond_mode)
{
	FILE  *f;
#define SSCANF_MAX_FILELEN_INT 2048
#define SSCANF_MAX_FILELEN INT_TO_STR(SSCANF_MAX_FILELEN_INT)
	const int MAX_FILELEN = SSCANF_MAX_FILELEN_INT;
	char  line[MAX_FILELEN], fname[MAX_FILELEN];
	int   cnt, pos, cond_num, ret_num;

	if (cond_mode)
		f=fopen(condition_fname.toUtf8().constData(), "r+t");
	else
		f=fopen(retention_fname.toUtf8().constData(), "r+t");

	if (!f) {
	  String type, name;
	  if (cond_mode) {
	    type = "condition";
	    name = condition_fname;
	  } else {	      
	    type = "retention";
	    name = retention_fname;
	  }	  
	  Form1->Memo1->Lines->Add("ERROR: could not read "+type+" layers list file: "+name);
	  return false;
	}


	Form1->Memo1->Lines->Add("");
	if (cond_mode) {
	  Form1->Memo1->Lines->Add("Loading condition file list from file "+condition_fname);
	  Form1->Memo1->Lines->Add("");
	} else 	{
	  Form1->Memo1->Lines->Add("Loading retention file list from file "+retention_fname);
	  Form1->Memo1->Lines->Add("");
	  Form1->Memo1->Lines->Add("Feature#  loss/gain   remains   change   final_retention_layer_weight");
	}

	cnt=0;
	while(fgets(line, MAX_FILELEN, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<2)
			continue;
		int ns = sscanf(line, "%i %" SSCANF_MAX_FILELEN "[^\r\n]s", &cond_num, fname);
		if (2!=ns) {
		  String sm;
		  if (cond_mode)
		    sm = "condition";
		  else
		    sm = "retention";
		  Form1->Memo1->Lines->Add("ERROR, could not parse line in "+sm+" list file: "+String(line));
		}

		if (strlen(fname)>1)
		{
			if (cond_mode) {
			  Form1->Memo1->Lines->Add(" * Applying condition layer number "+IntToStr(cond_num)+":");
			  if (Load_and_check_condition_matrix(String(fname))) {
			    Apply_condition_matrix(cond_num, 1);
			    condition.free_matrix_m();
			    Form1->Memo1->Lines->Add(" * Applied condition layer number "+IntToStr(cond_num));
			    Form1->Memo1->Lines->Add("");	
			    cnt++;
			  } else {
			    Form1->Memo1->Lines->Add("ERROR loading condition layer from file: "+String(fname));
			    fclose(f);
			    return false;
			  }
			} else {
			  Form1->Memo1->Lines->Add(" * Applying retention layer number "+IntToStr(cond_num)+":");
			  if (Load_and_check_retention_matrix(String(fname))) {
			    Apply_condition_matrix(cond_num, 0);
			    condition.free_matrix_m();
			    Form1->Memo1->Lines->Add(" * Applied retention layer number "+IntToStr(cond_num));
			    Form1->Memo1->Lines->Add("");
			    cnt++;
			  } else {
			    Form1->Memo1->Lines->Add("ERROR loading retention layer from file: "+String(fname));
			    fclose(f);
			    return false;
			  }
			}
		}
	}

	if (cond_mode)
		Form1->Memo1->Lines->Add("Read condition layers list file, number of condition layers: "+IntToStr(cnt));
	else
		Form1->Memo1->Lines->Add("Read retention layers list file, number of retention layers: "+IntToStr(cnt));

	Form1->Memo1->Lines->Add("");
	
	fclose(f);

	return true;
}

bool Load_and_check_LEC_matrixes()
{
	FILE  *f;
	char  line[512], fname[512];
	int   cnt, pos, num;

	f=fopen(LEC_fname.toUtf8().constData(), "r+t");

	if(!f)
	  return false;

	Form1->Memo1->Lines->Add("Loading local edge correction files list from file "+LEC_fname);

	cnt=0;
	while(fgets(line, 512, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<2)
			continue;
		sscanf(line, "%i %s", &num, fname);

		if (strlen(fname)>1)
		{
			if (Load_and_check_LEC_matrix(String(fname)))
			{
				Apply_LEC_matrix(num);
				condition.free_matrix_m();
				Form1->Memo1->Lines->Add("Applied local edge correction layer number "+IntToStr(num));
				cnt++;
			}
		}
	}

	Form1->Memo1->Lines->Add("Read files: edge corrections layers count = "+IntToStr(cnt));
	
	fclose(f);

	return true;
}

bool Load_sim_matrix()
{
	FILE  *f;
	char  line[MAX_LINE_LEN], fname[512], *items[MAX_SPP_COUNT];
	int   cnt, row, cnt0, cols, loop;
	float ftmp;

	f=fopen(comm_set.sim_matrix_name.toUtf8().constData(), "r+t");
	if (!f)
	  return false;

	Form1->Memo1->Lines->Add("Loading community connectivty similarity matrix from file "+comm_set.sim_matrix_name);

	cnt0=Count_rows(comm_set.sim_matrix_name.toUtf8().constData());
	if (cnt0<2)
	{
		Form1->Memo1->Lines->Add("Community connectivity similarity matrix file seems to have <2 rows"+comm_set.sim_matrix_name);
		fclose(f);
		return false;
	}

	sim_mat = matrix(0, cnt0, 0, cnt0);
	if (!sim_mat)
	{
		Form1->Memo1->Lines->Add("Out of memory allocating connectivity similarity matrix");
		fclose(f);
		return false;
	}

	row=cnt=0;
	while(fgets(line, MAX_LINE_LEN, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<2)
			continue;

		//      Form1->Memo1->Lines->Add(line);
		cols = tokenize(items, line);
		if (cols!=cnt0) {
		  ShowMessage("Inequal count of columns between rows in " + comm_set.sim_matrix_name + 
			      ", there are " + cols + " columns but " + cnt0 + " rows|");
		  fclose(f);
		  return false;
		}

		for(loop=0;loop<cols;loop++)
		{
			cnt=sscanf(items[loop],"%f", &ftmp);

			if (cnt!=1)
			{
				ShowMessage("Similarity matrix load failed: could not understand number "+(String)items[loop]);
				fclose(f);
				return false;
			}
			else if (ftmp<0)
			{
				ShowMessage("Similarity matrix load failed: negative numbers not allowed "+(String)items[loop]);
				fclose(f);
				return false;
			}

			sim_mat[row][loop] = ftmp;
		}
		row++;
	}

	comm_set.ct_cnt = row;
	Form1->Memo1->Lines->Add("Read similarity matrix with row count = " + IntToStr(comm_set.ct_cnt) + ", and same number of columns.");

	fclose(f);

	return true;
}

bool Load_comm_sim_matrix()
{
	FILE  *f;
	char  line[MAX_LINE_LEN], fname[512], *items[MAX_SPP_COUNT];
	int   cnt, row, cnt0, cols, loop;
	float ftmp;

	f=fopen(comm_set.comm_sim_matrix_name.toUtf8().constData(), "r+t");
	if (!f)
		return false;

	Form1->Memo1->Lines->Add("Loading community representation similarity matrix from file "+comm_set.comm_sim_matrix_name);

	cnt0=Count_rows(comm_set.comm_sim_matrix_name.toUtf8().constData());
	if (cnt0<2)
	{
		Form1->Memo1->Lines->Add("Community representation similarity matrix file seems to have <2 rows"+comm_set.comm_sim_matrix_name);
		return false;
	}

	comm_sim_mat = matrix(0, cnt0, 0, cnt0);
	if (!comm_sim_mat)
	{
		Form1->Memo1->Lines->Add("Out of memory allocating similarity matrix");
		return false;
	}

	row=cnt=0;
	while(fgets(line, MAX_LINE_LEN, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<2)
			continue;

		//      Form1->Memo1->Lines->Add(line);
		cols = tokenize(items, line);
		if (cols!=cnt0) {
		  ShowMessage("Column count (" + IntToStr(cols) + ") and row count (" + IntToStr(cnt0) + ") do not match in "+comm_set.comm_sim_matrix_name);
		  fclose(f);
		  return false;
		}

		for(loop=0;loop<cols;loop++)
		{
			cnt=sscanf(items[loop],"%f", &ftmp);

			if (cnt!=1)
			{
				ShowMessage("Community representation similarity matrix load failed: could not understand number "+(String)items[loop]);
				fclose(f);
				return false;
			}
			else if (ftmp<0.0)
			{
				ShowMessage("Similarity representation matrix load failed: negative numbers not allowed "+(String)items[loop]);
				fclose(f);
				return false;
			}

			comm_sim_mat[row][loop] = ftmp;
		}
		row++;
	}

	comm_set.comm_ct_cnt = row;
	Form1->Memo1->Lines->Add("Read community representation similarity matrix with row count = " + IntToStr(comm_set.comm_ct_cnt) +
				 ", and same number of columns.");

	fclose(f);

	return true;
}

int Load_ADMU_descriptions()
{
	FILE  *f;
#define SSCANF_MAX_FILELEN_INT 2048
#define SSCANF_MAX_FILELEN INT_TO_STR(SSCANF_MAX_FILELEN_INT)
	const int MAX_LINELEN = SSCANF_MAX_FILELEN_INT;
	char  line[MAX_LINELEN], fname[MAX_LINELEN];
	int   cnt, row_id, fid, in_cnt, row_count;
	float glob_weight, loc_weight;
	char  unit_name[MAX_LINELEN];

	f=fopen(ADM_set.ADM_weights_file.toUtf8().constData(), "r+t");

	if (!f)
	{
		Form1->Memo1->Lines->Add("Could not open "+ADM_set.ADM_weights_file);
		return false;
	}

	Form1->Memo1->Lines->Add("Loading ADMU descriptions and weights from file "+ADM_set.ADM_weights_file);

	// #"Rowid_","VALUE_","COUNT_","UNIT"
	// 0  106415  1.0  1.0  Northland - northern

	cnt = row_count = 0;
	while(fgets(line, MAX_LINELEN, f))
	{
		row_count++;
		if (line[0]=='#')
			continue;
		if (strlen(line)<8)
			continue;
		in_cnt = sscanf(line, "%i %f %f %" SSCANF_MAX_FILELEN "[^\r\n]s", &fid, &glob_weight, &loc_weight, &unit_name[0]);
		if (in_cnt!=4)
		{
			Form1->Memo1->Lines->Add("Wrong number - 4 expected - of elements in ADMU descriptions file on row "+IntToStr(row_count));
			fclose(f);
			return 0;
		}
		cnt++;
	}

	ADM_set.count = cnt;
	ADM_set.descriptions = new struct ADM_desc[cnt];

	rewind(f);

	cnt = 0;
	while(fgets(line, MAX_LINELEN, f))
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<8)
			continue;
		in_cnt = sscanf(line, "%i %f %f %" SSCANF_MAX_FILELEN "[^\r\n]s", &fid, &glob_weight, &loc_weight, &unit_name[0]);
		if (in_cnt!=4)
		{
			Form1->Memo1->Lines->Add("Wrong number - 4 expected - of elements in ADMU descriptions file on row "+IntToStr(row_count));
			fclose(f);
			return 0;
		}
		ADM_set.descriptions[cnt].number = cnt;
		ADM_set.descriptions[cnt].global_weight = glob_weight;
		ADM_set.descriptions[cnt].local_weight = loc_weight;
		ADM_set.descriptions[cnt].id_number = fid;
		ADM_set.descriptions[cnt].name = unit_name;
		// put id in map: ADMU_id_in_raster -> ADMU_seq_in_0_...
		ADMUs_id_to_seq[fid] = cnt;
		cnt++;
	}
	fclose(f);

	Form1->Memo1->Lines->Add("ADMU descriptions list file successfully loaded; count of ADMUs = "+IntToStr(ADM_set.count));
	return 1;
}

int Load_ADMU_weights_matrix()
{
	FILE  *f;
	char  line[100000], fname[512], *items[20000];
	int   cnt, row, cnt0, cols, loop;
	float ftmp;

	f=fopen(ADM_set.ADM_weight_matrix_file.toUtf8().constData(), "r+t");
	if (!f)
	{
		Form1->Memo1->Lines->Add("Could not open ADMU weight matrix file "+ADM_set.ADM_weight_matrix_file);
		return false;
	}

	Form1->Memo1->Lines->Add("Loading ADMU weight matrix file "+ADM_set.ADM_weight_matrix_file);


	if (!alloc_ADMU_data())
	{
		Form1->Memo1->Lines->Add("ERROR: failure allocating memory for ADMU operations.");
		fclose(f);
		return 0;
	}

	row=cnt=0;
	// This is very needed because resampling is done (and map_cnt is updated to the resample %)
	// after ADMU files are loaded.
	int expected_rows = 0;
	bool resampling_on = false;
	if ((resample_count>0) && (resample_count<map_cnt)) {
	  resampling_on = true;
	  expected_rows = resample_count;
	} else {
	  expected_rows = map_cnt;
	}
	while( fgets(line, 100000, f) && (row<expected_rows) )
	{
		if (line[0]=='#')
			continue;
		if (strlen(line)<2)
			continue;

		// Beware of the resampling monster!!!
		if (resampling_on && 0 == resamp_vec[row]) {
		  ShowMessage("contt!1: " + IntToStr(resamp_vec[row]));
		  continue;
		}
		cols = tokenize(items, line);
		if (cols!=ADM_set.count)
		{
			ShowMessage("Column count "+IntToStr(cols)+" not equal to ADMU count "+IntToStr(ADM_set.count)+" in "+ADM_set.ADM_weight_matrix_file);
			fclose(f);
			return false;
		}

		for(loop=0;loop<cols;loop++)
		{
			cnt=sscanf(items[loop],"%f", &ftmp);

			if (cnt!=1)
			{
				ShowMessage("ADMU weight matrix load failed: could not understand number "+(String)items[loop]);
				fclose(f);
				return false;
			}

			ADM_weights[row][loop] = ftmp;
		}
		row++;
	}
	if (row < expected_rows) {
	  ShowMessage("ADMU weight matrix load failed: not enough rows could be read: "+ IntToStr(row) + 
		      " (expected " + IntToStr(expected_rows) + ").");
	  fclose(f);
	  return false;
	} else if (row > map_cnt) {
	  // This excess is normal if using "resample species"
	  Form1->Memo1->Lines->Add("ADMU locations spatial data layer successfully loaded, but there is an excess of lines/rows in the ADMU weight matrix");
	  fclose(f);
	  // return 0;
	  return 1;  // If returned 0, then admu structures are not allocated... will crash sooner or later!
	}

	Form1->Memo1->Lines->Add("(" + IntToStr(row) +  ") ADMU features x " +
				 "(" + IntToStr(cols) + " )ADMUs weights matrix successfully loaded, with matching feature and ADMU counts ");
	fclose(f);
	return 1;
}


int Load_ADMU_layer()
{
	int x, y, num, loop;
	long ok_cnt, fail_cnt;

	ADMU_map.normalize = false;
	if (!ADMU_map.load_from_file(ADM_set.ADM_layer_file, 0, 0))
	{
		Form1->Memo1->Lines->Add("************** ERROR ***************");
		Form1->Memo1->Lines->Add("FAILURE when loading ADMU spatial data layer.");
		return 0;
	}

	ADMUs = imatrix(0, yd, 0, xd);
	if (!ADMUs)
	{
		ADMU_map.free_matrix_m();
		Form1->Memo1->Lines->Add("OUT OF MEMORY when loading ADMU spatial data layer.");
	}

	ok_cnt = fail_cnt = 0;
	for(y=0;y<yd;y++) {
	  for(x=0;x<xd;x++) {
	    int id = (int)(ADMU_map.m[y][x]);
	    if (id<0) {
	      fail_cnt++;
	      continue;
	    }
	    num = -1;
	    std::map<int, int>::iterator it = ADMUs_id_to_seq.find(id);
	    if (ADMUs_id_to_seq.end() == it) {
	      Form1->Memo1->Lines->Add("* ERROR, could not identify ADMU id: "+IntToStr(id)
				       +", found in coordinates x: "+IntToStr(x)+
				       +", y: "+IntToStr(y)+" in the adm. units raster map - but is not listed in the ADMU descriptions file!");
	      fail_cnt++;
	    } else {
	      num = it->second;
	    }
	    ADMUs[y][x] = num;  // after this ADMUs[][] ready to be used, with sequence number used in ADMU structures rather than the original IDs
	    ok_cnt++;
	  }
	}

	ADMU_map.free_matrix_m();

	Form1->Memo1->Lines->Add("Map of administrative units loaded successfully.");
	Form1->Memo1->Lines->Add("Locations with ok ADMU number = "+IntToStr(ok_cnt));
	Form1->Memo1->Lines->Add("Locations with missing or failed ADMU number = "+IntToStr(fail_cnt));
	Form1->Memo1->Lines->Add("");

	return 1;
}


int Load_and_check_ADMU_files()
{
	/*
[Administrative units]
use ADMUs = 1
ADMU descriptions file = ADMU_weights_all_local.txt
ADMU layer file = data\biog_unit.asc
ADMU weight matrix = ADMU_weights_matrix.txt
calculate local weights from condition = 1
# row count for per ADMU output curves = 100
*/
	Form1->Memo1->Lines->Add("");
	Form1->Memo1->Lines->Add("LOADING ADMINISTRATIVE UNITS (ADMU) INFO");

	if (!Load_ADMU_descriptions())
		return 0;
	if (!Load_ADMU_weights_matrix())
		return 0;
	if (!Load_ADMU_layer())
		return 0;

	Form1->Memo1->Lines->Add("ADMU LOAD COMPLETED");
	Form1->Memo1->Lines->Add("Identified "+IntToStr(ADMUs_id_to_seq.size())+" administrative units");
	calculate_and_output_ADMU_feature_weights();
	Form1->Memo1->Lines->Add("ADMU COMBINED WEIGHTS COMPUTED.");
	Form1->Memo1->Lines->Add("ADMU row count for per ADMU curves: "+IntToStr(ADM_set.row_count_for_per_admu_curves));
	
	if (ADM_set.row_count_for_per_admu_curves > 0) {
	  // +2 is for the extra lines (first and last)
	  const size_t MB = 1048576.0;
	  float required = float(ADM_set.count * 
				 ( (ADM_set.row_count_for_per_admu_curves + 2) * 
				   (9 + 4*7 + map_cnt*7) ) ) / MB;
	  const size_t MB_avail = 1048576.0 / 100;
	  unsigned int bytes;
	  Form1->Memo1->Lines->Add("Checking available space in " + outgridfn);
	  String where = checkAvailableDisk(outgridfn, bytes);
	  float avail =  float(bytes) / MB_avail;
	  Form1->Memo1->Lines->Add("ADMU: will generate per admin. unit output files, which (for " + IntToStr(ADM_set.count) + " adm. units and " + IntToStr(map_cnt) + " features) require over " + FloatToStrF(required, ffFixed, 4,4) + " MB of available disk space");
	  if ( avail > required ) {
	    Form1->Memo1->Lines->Add("Space available in " + where + ": " + FloatToStr(avail) +
				     " MB. OK.");		   
	  } else {
	    Form1->Memo1->Lines->Add("ERROR: not enough space available in " + where + ": " + FloatToStrF(avail, ffFixed, 4,4) + " MB.");
	    return 0;
	  }
	}

	return 1;
}

int grp_sort_function( const void *a, const void *b)
{
	int n1, n2;

	n1=*(int *)a;
	n2=*(int *)b;
	if (n1<n2)
		return -1;
	else
		return 1;
}

int get_grp_max(int *grpv)
{
	int loop, num, pos, loop2, found;

	pos=0;
	for(loop=0; loop<map_cnt; loop++)
	{
		num=spp[loop].grp_op1;
		if (num<0 || num>map_cnt)
			continue;

		found=false;

		for(loop2=0; loop2<pos; loop2++)
		{
			if (grpv[loop2]==num)
			{
				found=true;
				break;
			}
		}

		if (!found)
		{
			grpv[pos]=num;
			pos++;
		}
	}

	qsort((void *)grpv, pos, sizeof(int), grp_sort_function);

	//  for(loop=0;loop<pos;loop++)
	//    Form1->Memo1->Lines->Add("grp "+IntToStr(grpv[loop]));

	return pos;
}

// Reads the list of corridors layers and the layers
int
read_corridor_domain_layers_file_n_layers(Corr_settings& corr_set)
{
#define SSCANF_MAX_FILELEN_INT 2048
#define SSCANF_MAX_FILELEN INT_TO_STR(SSCANF_MAX_FILELEN_INT)
  const int MAX_FILELEN = SSCANF_MAX_FILELEN_INT;
  FILE  *f;
  char  line[MAX_FILELEN];
  
  f=fopen(corr_set.layers_fn.toUtf8().constData(), "r+t");
  if (!f)
    return 0;

  Form1->Memo1->Lines->Add("Loading corridor domain layers file list from file: " + corr_set.layers_fn);
  
  int line_text = 0;
  int cnt = 0;
  while(fgets(line, MAX_FILELEN, f)) {
    if (line[0]=='#') {
      line_text++;
      continue;
    }
    if (strnlen(line, MAX_FILELEN) < 2) {
      Form1->Memo1->Lines->Add("Error in corridor domain layers file, layer number " +  IntToStr(cnt) +
			       " line " + IntToStr(line_text) + ": '"+line+"', line too short, please check.");
      line_text++;
      continue;
    }

    float weight;
    char fname[MAX_FILELEN];
    int ns = sscanf(line, "%f %" SSCANF_MAX_FILELEN "[^\r\n]s", &weight, fname);
    if (2!=ns) {
      Form1->Memo1->Lines->Add("ERROR: could not parse line in corridors domain list file ("+String(fname)+"): "+String(line));
    }
	
    if (strnlen(fname, MAX_FILELEN) > 1) {
      corr_set.weights.push_back(weight);
      corr_set.file_names.push_back(String(fname));
      line_text++;
      cnt++;
    } else {
      Form1->Memo1->Lines->Add("Error in corridor domain layers file, layer number " +  IntToStr(cnt) +
			       ", in line " + IntToStr(line_text) + ": '"+fname+"', layer file name too short, please check.");
    }
  }
  fclose(f);  

  Form1->Memo1->Lines->Add("Corridor domain layers count = " + IntToStr(corr_set.file_names.size()) + 
			   ". Now trying to read layers...");

  // just to get a fake IG_settings with alpha=0 => no IG
  struct IGWrap {
    struct IG_settings *ptr;
    IGWrap(): ptr(new IG_settings()) {}
    ~IGWrap() { delete ptr; }
  } dummy;
  dummy.ptr->IGa=0.0f;

  //std::vector<String>::iterator it;
  corr_set.layers.resize(corr_set.file_names.size());
  for(size_t l=0; l < corr_set.file_names.size(); l++) { //it = corr_set.file_names.begin(); it != corr_set.file_names.end(); it++) {
    float** tmpm = NULL;
    float psum, sum;
    corr_set.layers[l].normalize = false;
    bool res = corr_set.layers[l].load_from_file(corr_set.file_names[l], 0, 0);
    if (res) {
      Form1->Memo1->Lines->Add("Successfully read corridor domain layer #"+IntToStr(l+1)+": " + corr_set.file_names[l]);
    } else {
      Form1->Memo1->Lines->Add("Error reading corridor domain layer #"+IntToStr(l+1)+": " + corr_set.file_names[l]);
      return 0;
    }

    int cols = corr_set.layers[l].cols();
    int rows = corr_set.layers[l].rows();
    if (cols != xd &&  rows != yd) {
      Form1->Memo1->Lines->Add("Error, the dimensions of the corridor domain layer #"+IntToStr(l+1)+": " + corr_set.file_names[l]+
			       ", are wrong: "+IntToStr(cols)+" columns and "+IntToStr(cols)+" rows.");
      return false;
    }
  }

  return corr_set.layers.size();
}

// Drop the spp/features identified as "empty" after initial counting.
// Must be done before populating the vmat
bool
drop_0_occur_features(size_t drop_0_count)
{
  drop_0_relocation_idx.resize(map_cnt);
  struct sp* new_spp = new struct sp[map_cnt]; //MAX_SPP_COUNT];
  int new_map_cnt = 0;
  for (int i=0; i<map_cnt; i++) {
    if (!spp[i].drop) {
      // Note: the BQP_... vectors are not yet initialized, safe to move them without care
      new_spp[new_map_cnt] = spp[i];
      drop_0_relocation_idx[i] = new_map_cnt;    // spp i in spp list file -> new_map_cnt
      new_map_cnt++;
    } else {
      drop_0_relocation_idx[i] = i;  // not really relocated, just the same idx for convenience
    }
  }
  // switch from the initial spp list to the new, filtered one
  struct sp* spp_tmp = spp;
  spp = new_spp;
  delete [] spp_tmp;

  return true;
}

// before potential memory crashes
void
flush_kludge()
{
  std::cout.flush();
  std::cerr.flush();
  fflush(stdout);
  fflush(stderr);
}

bool
preload_and_count_feature_data_layers(float**& tmpm, int** occur_count_matrix)
{
  Form1->Memo1->Lines->Add("======------ Allocating memory for input layers... ------======");
#ifdef NONCOMPACT_VMAT_FLOATP
  Form1->Memo1->Lines->Add("Using compact structure (floatp version),  sizeof(Occur_Container): " + IntToStr(sizeof(Biodiv_Features_Occur_Container)));
#elif COMPACT_VMAT_LOOKUP

  // count effective occurrencies -> we then know how to allocate compact vmat
  time_t ti, tf;
  Form1->Memo1->Lines->Add("Using compact structure (lookup version),  sizeof(Occur_Container): " + 
			   IntToStr(sizeof(Biodiv_Features_Occur_Container)) +
			   ", sizeof(allocated occurrence): " + IntToStr(Biodiv_Features_Occur_Container::allocated_occurence_size())
			   );
  Form1->Memo1->Lines->Add("Counting effective occurrencies");	
  ti = time(NULL);
  Form1->Memo1->Lines->Add("Time: " + IntToStr(ti) + " = " + zCurrentTime());

  if (NULL == occur_count_matrix) {
    occur_count_matrix = imatrix(0, yd, 0, xd);
    memset(&occur_count_matrix[0][0], 0, sizeof(int)*xd*yd);
    /*
    for(int y = 0; y < yd; y++) {	    
      for(int x = 0; x < xd; x++) {
	occur_count_matrix[y][x] = 0;
      }
      }*/
  }

  size_t non_missing_accum = 0;
  size_t drop_0_count = 0;
  for(int loop = 0; loop < map_cnt; loop++) {
    size_t non_missing = 0;
    bool ok = Read_input_grid_file_just_to_count(loop, spp[loop].fname.toUtf8().constData(), tmpm, occur_count_matrix, non_missing);
    non_missing_accum += non_missing;    // accumulates counts of non-missing cells throughout layers
    if (!tmpm)   // what could make tmpm=0 ? Only at the beginning (obsmap[0])
      tmpm = obsmap[loop].m;

    if (!ok) {
      Form1->Memo1->Lines->Add("Counting - could not read in "+spp[loop].fname);
      return false;
    }
    if ( (obsmap[loop].rows()!=yd) ||
	 (obsmap[loop].cols()!=xd)     ) {
      Form1->Memo1->Lines->Add("Error: map dimensions do not match for feature/species map "+spp[loop].fname);
      return false;
    }

    if (settings_drop_0_occur_layers && 0 == non_missing) {
      Form1->Memo1->Lines->Add(" *** Note: no occurencies were found inside the analysis area. Dropping/ignoring this feature ("+spp[loop].fname+").");
      drop_0_count++;
      spp[loop].drop = true;
    }	  

    Form1->Memo1->Lines->Add("Feature file #"+IntToStr(loop+1)+ ": "+spp[loop].fname+ 
			     "; non-missing cells: "+IntToStr(non_missing));
			     //+ "; their sum: " + FloatToStr(sum)  // not counted here
  }
	
  // Drop empty spp/features
  if (settings_drop_0_occur_layers) {
    if (drop_0_count <= 0) {
      Form1->Memo1->Lines->Add("");
      Form1->Memo1->Lines->Add("All the features have occurrencies inside the analysis area. Not dropping any feature.");
    } else {
      Form1->Memo1->Lines->Add("");
      Form1->Memo1->Lines->Add("Dropping "+IntToStr(drop_0_count)+" features...");
      bool drop_ok = drop_0_occur_features(drop_0_count); // updates map_cnt!!!
      if (drop_ok) {
	Form1->Memo1->Lines->Add(" successfully dropped "+IntToStr(drop_0_count)+" features. Remaining features: "
				 +IntToStr(drop_0_count)+" features...");
      } else {
	Form1->Memo1->Lines->Add("ERROR: something went wrong while trying to drop features.");
	Form1->Memo1->Lines->Add("There is no solution for this. Giving up.");
	graceful_exit();
      }
    }
  }

  tf = time(NULL);
  // occurrencies have been counted feature by feature in the 
  // calls to Read_input_grid_file_just_to_count(...)
  size_t global_nonzero_occur = 0;    // non-empty prows
  for(int y = 0; y < yd; y++) {
    for(int x = 0; x < xd; x++) {
      if (occur_count_matrix[y][x] > 0 )
	global_nonzero_occur++;
    }
  }
  size_t mem_base = xd*yd*sizeof(Biodiv_Features_Occur_Container*) +     // 2D matrix of pointers
    global_nonzero_occur * sizeof(Biodiv_Features_Occur_Container);      // effective containers (globally non-missing cells)
  size_t mem_occurs = non_missing_accum*Biodiv_Features_Occur_Container::allocated_occurence_size();   // total non-missing occurrencies
  size_t mem_size = mem_base + mem_occurs;
  Form1->Memo1->Lines->Add("Total # of feature layers: " + IntToStr(map_cnt));
  Form1->Memo1->Lines->Add("Total # of cells in every layer: " + IntToStr(xd*yd));
  Form1->Memo1->Lines->Add("Effective # of cells (cells with any occurrence, globally): " + IntToStr(global_nonzero_occur));
  float non_missing_avg = (float)non_missing_accum/map_cnt;
  Form1->Memo1->Lines->Add("Total # of occurrencies: " + IntToStr(non_missing_accum) + 
			   ", counted across " + IntToStr(map_cnt) + " layers; average across layers: " +
			   FloatToStr(non_missing_avg));
  Form1->Memo1->Lines->Add("Percentage of globally effective occurrencies = " +
			   FloatToStr(100.0f*non_missing_avg/global_nonzero_occur) + " % (relative to total effective cells) = " +
			   FloatToStr(100.0f*non_missing_avg/(xd*yd)) + " % (relative to total cells)");
  Form1->Memo1->Lines->Add("Memory required for biodiversity features: " + 
			   FloatToStrF((float)mem_size/MB_in_bytes, ffFixed, 4, 4) + " MB, " +
			   "of which base memory: " + FloatToStrF((float)mem_base/MB_in_bytes, ffFixed, 4, 4) + " MB, " +
			   "and occurrencies require: " + FloatToStrF((float)mem_occurs/MB_in_bytes, ffFixed, 4, 4) + " MB");

  Form1->Memo1->Lines->Add("Time now: " + IntToStr(tf) + " = " + zCurrentTime());
  Form1->Memo1->Lines->Add("Finished preload/counting effective occurrencies. Elapsed: " + IntToStr(tf-ti) + " seconds.");
  Form1->Memo1->Lines->Add("======------ Allocating memory for input layers: finished successfully ------======\n");
  /*
  // this would print some columns for debugging
  size_t print_col = 240;
  for(y = 0; y < yd; y++) {
  Form1->Memo1->Lines->Add(" " + IntToStr(print_col) + "," + IntToStr(y) + ": " + IntToStr(occur_count_matrix[y][print_col]) + " "
  + "obsmap.m: " + QString().sprintf("%.6e", obsmap[loop-1].m[y][print_col]) + " = "
  + QString().sprintf("%.6e", tmpm[y][print_col]));
  }
  */
#else
#error "Cannot find which type of vmat to use"
#endif

  return true;
}

bool
load_and_transform_feature_data_layers(float**& tmpm, int**& occur_count_matrix)
{
  Form1->Memo1->Lines->Add("******=====----- Loading feature (e.g., species) data layers -----=====**********");

  // this is to collect the _IG_DS_ST_LEC_CT_RT_MT_IT... appendix for the names of trans output layers
  String* transf_layers_str_append;
  String tr_append_ig = "_IG";
  String tr_append_ds = "_DS";
  String tr_append_cst = "_CST";
  String tr_append_ct = "_CT";
  String tr_append_rt = "_RT";
  String tr_append_mct = "_MCT";
  String tr_append_ia = "_IA";
  float* transf_buff = NULL;
  if (transf_layers_set.some_output) {
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("====== Note: generating feature (e.g., species) transformed layers =====");
    transf_layers_str_append = new String[map_cnt];
    // BEWARE: IG and DS layers are generated one-by-one => do this check here in advance -and only once.
    // For other types of layers (generated as groups) this createDir is done inside generate_output_transf_layers
    if ( transf_layers_set.output_ig_layers) {
      createSubDirIfNeeded(output_transf_ig_layers_dn, "");
    }
    if ( transf_layers_set.output_ds_layers) {
      createSubDirIfNeeded(output_transf_ds_layers_dn, "");
    }
    // allocate layer buffer for saving transformed features. Needs to be a traditional contiguous buffer
    uint64_t mem_req = xd*yd*sizeof(float);
    Form1->Memo1->Lines->Add("This requires additional memory:  "+FloatToStrF(float(mem_req)/MB_in_bytes, ffFixed, 4, 4) + " MB");
    transf_buff = new float[xd*yd]; // matrix(0, yd, 0, xd) would not work, non-contiguous rows!
    if (NULL == transf_buff) {
      Form1->Memo1->Lines->Add("Memory allocation error. Giving up!");
      delete [] transf_layers_str_append;
      return false;
    }

    Form1->Memo1->Lines->Add("Memory allocated successfully... ");
    Form1->Memo1->Lines->Add("");
  }
  
  //tmpm = NULL;  <- this would force a deallo/alloc, unnecessary
  float **wtmpm = NULL;
  // First loop: load
  for(int loop=0;loop<map_cnt;loop++) {
    if (IG_set.use_IG && IG_set.use_IGw) {
      bool ok = Read_IGw_file(loop, IG_set.fnames[loop].toUtf8().constData(), wtmpm);
      if (!wtmpm)
	wtmpm = wmap[loop].m;
      
      if (!ok){
	Form1->Memo1->Lines->Add("Error reading IG weight map "+IG_set.fnames[loop]);
	return false;
      } else {
	Form1->Memo1->Lines->Add("Info-gap robustness/opportunity analysis (DISTRIBUTION DISCOUNTING): loaded uncertainty distribution map from file"+IG_set.fnames[loop]);
      }

      if ( (wmap[loop].rows()!=yd) || (wmap[loop].cols()!=xd)     ) {
	Form1->Memo1->Lines->Add("Error: map dimensions do not match for distributional uncertainty map, file: "+
				 IG_set.fnames[loop]+" (corresponding to feature file: "+spp[loop].fname+")");
	return false;
      }
      if (IG_set.normalize_IGw)
	wmap[loop].average_normalize();
    }

    bool ok = Read_input_grid_file(loop, spp[loop].fname.toUtf8().constData(), tmpm);
    if (!tmpm)
      tmpm = obsmap[loop].m;

    if (!ok) {
      Form1->Memo1->Lines->Add("Could not read in "+spp[loop].fname);
      return false;
    }
    
    if ( (obsmap[loop].rows()!=yd) || (obsmap[loop].cols()!=xd) ) {
      Form1->Memo1->Lines->Add("Error: map dimensions do not match for feature (e.g. species) map "+spp[loop].fname);
      return false;
    }

    int nan_count = 0;
    for(int y=0; y<yd; y++) {
      for(int x=0; x<xd; x++)	{
	if (obsmap[loop].m[y][x]!=-1) {
	  //  if (obsmap[loop].m[y][x] >= 0.00001f )
	  if (!vmat[y][x]) {
	    try {
#if NONCOMPACT_VMAT_FLOATP
	      //vmat[y][x]= (float *)calloc(map_cnt, sizeof(float));    // COMPACT_VMAT
	      vmat[y][x].reserve(map_cnt, map_cnt);
#elif COMPACT_VMAT_LOOKUP
	      vmat[y][x].reserve(occur_count_matrix[y][x], map_cnt);
#else
  #error "Cannot find which type of vmat to use"
#endif
	    } catch(std::bad_alloc&) {
	      Form1->Memo1->Lines->Add("****** CRITICAL ERROR: system out of memory, unable to allocate memory for feature "+
				       IntToStr(loop)+"("+spp[loop].fname+") ******");
	      Form1->Memo1->Lines->Add("There is no solution for this. More available memory is required. Giving up.");
	      graceful_exit();
	    }
	    if (!vmat[y][x]) {
	      String msg = "out of memory loading data for spp #"+IntToStr(loop);
	      Form1->Memo1->Lines->Add(msg);
	      ShowMessage(msg);
	      return false;
	    }
	    
	    //for(z=0; z<loop; z++)
	    //        vmat[y][x][z]=-1;     // COMPACT_VMAT
	    vmat[y][x].fill_empty(loop);
	  }
	  // vmat[y][x][loop]=obsmap[loop].m[y][x]; // COMPACT_VMAT
	  vmat[y][x].insert(loop, obsmap[loop].m[y][x]);   // inserts map value into pos/feature 'loop'
#if ZCORE_DEBUG
	  if ( isnan(obsmap[loop].m[y][x]) ) {
	    ShowMessage("data_input_v2(): got NaN in obsmap / vmat! " + IntToStr(loop) + ", count: " + IntToStr(++nan_count));
	    //exit(1);
	  }
#endif
	}
	else if (vmat[y][x])
	  //vmat[y][x][loop]=-1;     // COMPACT_VMAT
	  vmat[y][x].set_empty(loop);
      }
    }
    // here I could check that the sizes of the vmat[y][x] are consistent with occur_count_matrix[y][x]

    // This cannot be done until now, when vmat is allocated and filled in...
    if (IG_set.use_IG && IG_set.use_IGw) {
      // First component of the suffix for transformed layers, may more may come later on...
      if ( transf_layers_set.some_output) {
	transf_layers_str_append[loop] += tr_append_ig; // "_IG";
	// must generate IG layers one-by-one cause DS comes immediately after, so outside of 'loop'
	// it would be too late
	if (transf_layers_set.output_ig_layers) {
	  generate_output_transf_single_layer(loop, transf_buff, xd, yd, output_transf_ig_layers_dn, transf_layers_str_append[loop]);
	}			  
      }		  
    }
    
    if (IG_set.use_IG && IG_set.use_IGw) {
      if (loop>0)
	wmap[loop].free_matrix_m();
    }
    if (loop>0)
      obsmap[loop].free_matrix_m();
    Form1->Memo1->Lines->Add("* Loaded biodiversity feature file #"+IntToStr(loop+1)+", "+spp[loop].fname +
			     ", non-missing cells:" + IntToStr(obsmap[loop].non_missing_cnt) +
			     ", their sum: " + FloatToStr(obsmap[loop].sum_orig)
			     );
  }

  if (!use_occur_size_weights_correct_landscape_fraction)
    occur_size_weights_map.free_matrix_m();
  
  flush_kludge();

#ifdef COMPACT_VMAT_LOOKUP
  if (occur_count_matrix)
    free_imatrix(occur_count_matrix, 0, yd, 0, xd);
#endif
  obsmap[0].free_matrix_m(); // the [0] are the only ones with data allocated
  obsmap[0].force_free_matrix_m();

  if (IG_set.use_IG && IG_set.use_IGw) {
    wmap[0].free_matrix_m();
    wmap[0].force_free_matrix_m();
  }

  // dist. centers are required already here (for distribution smoothing)
  // if "correcting" non-uniform cell sizes
  init_output_get_distr_centers_and_sums(false);

  if (use_smoothing) {
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("******=====----- Distribution smoothing is on. Transforming input layers... -----=====******");
    time_t t1_ds = time(NULL);
    for (int loop=0;loop<map_cnt;loop++) {
      if (spp[loop].alpha>0.0f) {
	Form1->Memo1->Lines->Add("* Applying distribution smoothing on feature #"+IntToStr(loop+1)+", "+spp[loop].fname);
	Do_smoothing_sp(loop);
	if ( transf_layers_set.some_output) {
	  transf_layers_str_append[loop] += tr_append_ds;
	  if (transf_layers_set.output_ds_layers) {
	    // Form1->Memo1->Lines->Add("* Writing layer after distrib. smoothing for input file #"+IntToStr(loop+1)+", "+spp[loop].fname);
	    generate_output_transf_single_layer(loop, transf_buff, xd, yd, output_transf_ds_layers_dn, transf_layers_str_append[loop]);
	  }
	}
      }
    }
    time_t t2_ds = time(NULL);

    Form1->Memo1->Lines->Add("Freeing structures required for distribution smoothing... ");
    // smoothing matrices also used for IA and comm connectivity
    if (!use_interactions && !((comm_set.load_sim_matrix && comm_set.apply_to_conn)))
      Free_smoothing_matrixes();
    Form1->Memo1->Lines->Add(" ...Done.");

    Form1->Memo1->Lines->Add("******=====----- Finished distribution smoothing transformations in "+
			     FloatToStr(t2_ds-t1_ds)+" seconds. -----=====******");
    Form1->Memo1->Lines->Add("");
  }

  // But if the first allocation was done inside the gridmap objects, 1 buffer would remain... kill it:
  
  Form1->Edit32->Text=IntToStr(map_cnt);
  
  if (comm_set.load_sim_matrix && comm_set.apply_to_repr) {
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("******=====----- 'load similarity matrix' is on, using similarity matrix for community similarity transform (with 'apply to representation'). Transforming input layers... -----=====******");
    time_t t1_cst = time(NULL);
    if (!Load_comm_sim_matrix())
      Form1->Memo1->Lines->Add("Error loading community representation similarity matrix"+comm_set.comm_sim_matrix_name);
    else {
      Form1->Memo1->Lines->Add("Community representation similarity matrix loaded successfully from file"+comm_set.comm_sim_matrix_name);
      if (comm_set.apply_to_repr) {
	Form1->Memo1->Lines->Add("Applying community similarity matrix to expand effective representation");
	if (!Expand_community_similarity()) {
	  ShowMessage("Failure expanding community similarity.");
	} else {
	  // Expand_community_similarity => comm_set.comm_ct_cnt
	  if ( transf_layers_set.some_output) {
	    for(size_t i=0; i<comm_set.comm_ct_cnt; i++) {
	      transf_layers_str_append[i] += tr_append_cst;
	    }
	    if (transf_layers_set.output_cst_layers) {
	      Form1->Memo1->Lines->Add("\nWriting community similarity transformed (CST) layers to disk. This can take a while...");
	      generate_output_transf_layers(transf_buff, output_transf_cst_layers_dn, transf_layers_str_append, tr_append_cst);
	    }				    
	  }
	}
      }
    }
    time_t t2_cst = time(NULL);
    Form1->Memo1->Lines->Add("******=====----- Finished community similarity transformations in "+
			   FloatToStr(t2_cst-t1_cst)+" seconds. -----=====******");
    Form1->Memo1->Lines->Add("");
  }
  
  if (use_LEC) {
    if (!use_groups) {
      Form1->Memo1->Lines->Add("Local edge corrections cannot be used if groups have not been specified");
    }
    else {
      if (Load_and_check_LEC_matrixes()) {
	Form1->Memo1->Lines->Add("Applied local edge correction layers through the landscape.");
	Form1->Memo1->Lines->Add("");
	condition.free_matrix_m();
	
	/*
	// LEC is actually disabled
	if ( transf_layers_set.output_lec_layers) {
	for(size_t i=0; i<map_cnt; i++) {
	transf_layers_str_append[i] += "_LEC";
	}
	}
	*/
      }
      else
	Form1->Memo1->Lines->Add("ERRORS applying local edge correction matrixes. ");
    }
  }
  
  if (use_condition) {
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("******=====----- 'use condition' is on. Transforming input layers... -----=====******");
    time_t t1_condition = time(NULL);
    if (!use_groups) {
      Form1->Memo1->Lines->Add("ERROR: Condition cannot be used if groups have not been specified");
    } else {
      if (Load_and_check_condition_matrixes(1)) {
	Form1->Memo1->Lines->Add("Applied condition layers through the landscape.");
	Form1->Memo1->Lines->Add("Consequently, initial remaining fractions do not necessarily start from 1.");
	condition.free_matrix_m();
	
	if ( transf_layers_set.some_output) {
	  for(size_t i=0; i<map_cnt; i++) {
	    transf_layers_str_append[i] += tr_append_ct;
	  }				  
	  if ( transf_layers_set.output_ct_layers) {
	    Form1->Memo1->Lines->Add("\nWriting condition transformed (CT) layers to disk. This can take a while...");
	    generate_output_transf_layers(transf_buff, output_transf_ct_layers_dn, transf_layers_str_append, tr_append_ct);
	  }
	}
      } else {
	Form1->Memo1->Lines->Add("ERRORS while loading and applying condition matrixes. ");
	return false;
      }
    }
    time_t t2_condition = time(NULL);
    Form1->Memo1->Lines->Add("******=====----- Finished condition transformations in "+
			   FloatToStr(t2_condition-t1_condition)+" seconds. -----=====******");
    Form1->Memo1->Lines->Add("");
  }
  
  if (use_retention) {
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("******=====----- 'use retention' is on. Transforming input layers... -----=====******");
    time_t t1_retention = time(NULL);
      if (!use_groups) {
	Form1->Memo1->Lines->Add("Retention cannot be used if groups have not been specified");
      }
      else {
	if (Load_and_check_condition_matrixes(0)) { // ret mode
	  Form1->Memo1->Lines->Add("Applied Retention layers through the landscape.");
	  Form1->Memo1->Lines->Add("");
	  condition.free_matrix_m();
	  
	  if (transf_layers_set.some_output) {
	    for(size_t i=0; i<map_cnt; i++) {
	      transf_layers_str_append[i] += tr_append_rt;
	    }
	    if (transf_layers_set.output_rt_layers) {
	      Form1->Memo1->Lines->Add("\nWriting retention transformed (RT) layers to disk. This can take a while...");
	      generate_output_transf_layers(transf_buff, output_transf_rt_layers_dn, transf_layers_str_append, tr_append_rt);
	    }
	  }
	} else {
	  Form1->Memo1->Lines->Add("ERRORS while loading and applying retention matrixes. ");
	  return false;
	}
      }
      time_t t2_retention = time(NULL);
      Form1->Memo1->Lines->Add("******=====----- Finished retention transformations in "+
			       FloatToStr(t2_retention-t1_retention)+" seconds. -----=====******");
      Form1->Memo1->Lines->Add("");
  }
  
  if (comm_set.load_sim_matrix && comm_set.apply_to_conn) {
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("******=====----- 'load similarity matrix' is on, using similarity matrix for matrix connectivity transform (with 'apply to connectivity'). Transforming input layers... -----=====******");
    time_t t1_mct = time(NULL);
    Form1->Memo1->Lines->Add("");
    if (!Load_sim_matrix())
      Form1->Memo1->Lines->Add("Error loading community similarity matrix"+comm_set.sim_matrix_name);
    else {
      Form1->Memo1->Lines->Add("Community similarity matrix loaded successfully from file"+comm_set.sim_matrix_name);
      if (comm_set.apply_to_conn) {
	Form1->Memo1->Lines->Add("Applying community similarity matrix to connectivity");
	if (!Do_matrix_smoothings())
	  ShowMessage("Failure doing matrix smoothings.");
	else {
	  if ( transf_layers_set.some_output) {
	    for(size_t i=0; i<comm_set.ct_cnt; i++) {
	      transf_layers_str_append[i] += tr_append_mct;
	    }
	    if (transf_layers_set.output_mct_layers) {
	      Form1->Memo1->Lines->Add("\nWriting matrix connectivity transformed (MCT) layers to disk. This can take a while...");
	      generate_output_transf_layers(transf_buff, output_transf_mct_layers_dn, transf_layers_str_append, tr_append_mct);
	    }
	  }
	}
      }
      if ((comm_set.apply_to_repr) && (comm_set.apply_to_conn))
	Form1->Memo1->Lines->Add("WARNING: condition with complex interpretation. You are applying community similarity matrix both on representation and connectivity.");
    }
    time_t t2_mct = time(NULL);
    Form1->Memo1->Lines->Add("******=====----- Finished matrix connectivity transformations in "+
			   FloatToStr(t2_mct-t1_mct)+" seconds. -----=====******");
    Form1->Memo1->Lines->Add("");
  }
  
  if (use_interactions) {
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("******=====----- 'use interactions' is on. Transforming input layers... -----=====******");
    time_t t1_interactions = time(NULL);
    if (!Load_interactions(transf_layers_str_append))
      Form1->Memo1->Lines->Add("Error loading feature interactions from file: "+ia_fname);
    else
      Form1->Memo1->Lines->Add("Interactions between biodiversity features loaded successfully from file: "+ia_fname);
    Form1->Memo1->Lines->Add("");
    
    if ( transf_layers_set.some_output) {
      // interactions appendices "_n_IA" have been generated before
      if (transf_layers_set.output_ia_layers) {
	Form1->Memo1->Lines->Add("\nWriting interactions (IA) transformed layers to disk. This can take a while...");
	if (generate_output_transf_layers(transf_buff, output_transf_ia_layers_dn, transf_layers_str_append)) {
	  Form1->Memo1->Lines->Add("Interactions transformed layers written succesfully.\n");
	} else {
	  Form1->Memo1->Lines->Add("WARNING: failed to write interactions transformed layers.\n");
	}
      }
    }

    time_t t2_interactions = time(NULL);
    Form1->Memo1->Lines->Add("******=====----- Finished interactions transformations in "+
			   FloatToStr(t2_interactions-t1_interactions)+" seconds. -----=====******");
    Form1->Memo1->Lines->Add("");
  }
  
  flush_kludge();

  if (smoothing_is_allocated()) {
    Form1->Memo1->Lines->Add("Freeing FFTWF structures... ");
    Free_smoothing_matrixes();
    Form1->Memo1->Lines->Add(" ...Done.");
  }

  // transformed layers
  if ( transf_layers_set.output_final_layers) {
    Form1->Memo1->Lines->Add("\nWriting transformed layers (final) to disk. This can take a while...");
    if (generate_output_transf_layers(transf_buff, output_transf_final_layers_dn, transf_layers_str_append)) {
      Form1->Memo1->Lines->Add("Final transformed layers written succesfully.\n");
    } else {
      Form1->Memo1->Lines->Add("WARNING: failed to write transformed layers.\n");
    }
    delete [] transf_layers_str_append;
    
    if (transf_buff) {
      delete [] transf_buff;
      transf_buff = NULL;
    }
  }

  Form1->Memo1->Lines->Add("");
  Form1->Memo1->Lines->Add("******=====----- Finished loading and transformation of feature data layers -----=====*********");
  Form1->Memo1->Lines->Add("");

  return true;
}

void
gen_random_omissions_percentage_simple_impl(float pc)
{
  if (random_omissions_max_feature<=0) {
    Form1->Memo1->Lines->Add("Note: not doing anything because the maximum feature number for omissions is: "+
			     IntToStr(random_omissions_max_feature));
    return;
  }

  Form1->Memo1->Lines->Add("Calulating current total repr values of features so that later on the remaining occurrencies can be rescaled back... ");
  std::vector<float> feat_repr_before;
  feat_repr_before.assign(random_omissions_max_feature, 0.0f);  // random_omissions_max_feature <= map_cnt
  for(int y=0; y<yd; y++) {
    for(int x=0; x<xd; x++) {
      const Biodiv_Features_Occur_Container& rowp = vmat[y][x];
      for(int s=rowp.first(); s != rowp.overflow() && s < random_omissions_max_feature; s=rowp.next(s)) {
	feat_repr_before[s] += rowp[s];
      }
    }
  }  

  int64_t tot_cnt = 0;
  int64_t zeroed_cnt = 0;
  pc = pc/100.0f;
  float running_pc = 0.0f;
  for(int y=0; y<yd; y++) {
    for(int x=0; x<xd; x++) {
      Biodiv_Features_Occur_Container& rowp = vmat[y][x];
      for(int s=rowp.first(); s != rowp.overflow() && s < random_omissions_max_feature; s=rowp.next(s)) {
	if (pc >= running_pc) {
	  rowp[s] = 0.0f;
	  zeroed_cnt++;
	}
	tot_cnt++;
	running_pc = float(zeroed_cnt)/float(tot_cnt);
      }
    }
  }
  Form1->Memo1->Lines->Add("Finished. Randomly zeroed "+IntToStr(zeroed_cnt)+
			   " occurrencies out of "+IntToStr(tot_cnt)+ " (including the first "+
			   IntToStr(random_omissions_max_feature)+ " features) "
			   " => "+
			   FloatToStr(100.0f*float(zeroed_cnt)/float(tot_cnt))+"%");

  
  Form1->Memo1->Lines->Add("Now rescaling values so that every feature recover its initial distrib sum... ");
  Form1->Memo1->Lines->Add(" calculating current (randomly discounted) repr values ... ");
  std::vector<float> feat_repr_after;
  feat_repr_after.assign(random_omissions_max_feature, 0.0f);   // random_omissions_max_feature <= map_cnt
  for(int y=0; y<yd; y++) {
    for(int x=0; x<xd; x++) {
      const Biodiv_Features_Occur_Container& rowp = vmat[y][x];
      for(int s=rowp.first(); s != rowp.overflow() && s < random_omissions_max_feature; s=rowp.next(s)) {
	feat_repr_after[s] += rowp[s];
      }
    }
  }  
  Form1->Memo1->Lines->Add(" calculating ratios before/after... ");
  std::vector<float> feat_repr_ratio;
  feat_repr_ratio.assign(random_omissions_max_feature, 0.0f);
  for (int i=0; i<feat_repr_ratio.size(); i++) {
    if (feat_repr_after[i] > 0.0f)
      feat_repr_ratio[i] = feat_repr_before[i] / feat_repr_after[i];
    else 
      feat_repr_ratio[i] = 0.0f;
  }
  Form1->Memo1->Lines->Add(" rescaling back to before... ");
  for(int y=0; y<yd; y++) {
    for(int x=0; x<xd; x++) {
      Biodiv_Features_Occur_Container& rowp = vmat[y][x];
      for(int s=rowp.first(); s != rowp.overflow() && s < random_omissions_max_feature; s=rowp.next(s)) {
	rowp[s] *= feat_repr_ratio[s];
      }
    }
  }
  Form1->Memo1->Lines->Add("Finished.");

}
	
// load everything except command line and dat
bool data_input_v2()
{
  if (removal_rule==4)
    new_spp_file_format=1;
  else
    new_spp_file_format=0;

  // load spp
  int spf_cnt=Read_spp_file_names();
  if (spf_cnt<=0) {
    Form1->Memo1->Lines->Add("Failure reading features/species list file - exiting.");
    return false;
  }
  
  Form1->Memo1->Lines->Add("Row count in the features/species list file: "+IntToStr(map_cnt));
  if (neg_weights_used) {
    Form1->Memo1->Lines->Add("**** NEGATIVE weights encountred in features/species list file; => running in MCZ multicriterion mode. ****");
    Form1->Memo1->Lines->Add("");
  }
  if (IG_set.use_IG && IG_set.use_IGw) { // IG used
    int igf_cnt=Read_IG_file_names();
    if (spf_cnt!=igf_cnt) {
      ShowMessage("Number of rows in spp and info-gap weights list files do not match. Exiting.");
      return false;
    }
  }

  // load the first spp file to tstmap
  if (!tstmap.load_from_file(spp[0].fname, 0, 0, true)) // set global projection
    return false;

  xd=tstmap.cols();
  yd=tstmap.rows();
  tstmap.free_matrix_m();


  if (mask_data) {
    if (!area_mask.load_from_file(area_mask_fname, 0, 0)) {
      Form1->Memo1->Lines->Add("ERROR while loading analysis area mask!");
      return false;
    } else {
      Form1->Memo1->Lines->Add("Analysis area mask loaded.");
    }
  }

  if (use_SSI) {
    SSI_spp_cnt=Read_SSI_file_names();
    if (SSI_spp_cnt<=0) {
      ShowMessage("Failure reading SSI features/species list file - exiting.");
      return false;
    }
    Form1->Memo1->Lines->Add("Number of SSI features/species = "+IntToStr(SSI_spp_cnt));
  }

  if (use_interactions) {
    ia_cnt = Count_rows(ia_fname.toUtf8().constData());
    if (ia_cnt<=0)
      use_interactions=0;
    Form1->Memo1->Lines->Add("Number of interaction pairs = "+IntToStr(ia_cnt));
  }

  Check_maps_et_alloc();
  //ShowMessage("IPD 2");

  if (use_PLULA) {
    use_PLULA=Load_PLULA(PLULAfn);
    if (!use_PLULA) {
      Form1->Memo1->Lines->Add("ERROR: load of planning unit layer failed.");
      ShowMessage("Could not successfully load planning unit layer.");
      return false;
    }
    
    if (warp_factor>1) {
      //warp_factor = 1;
      Form1->Memo1->Lines->Add("WARNING: removing planning units, warp factor set to "+IntToStr(warp_factor)+" (PLU) even though planning units are used.");
    }
    
    if (BLP>0) {
      Form1->Memo1->Lines->Add("WARNING: planning units do not presently work with BLP. Switching BLP off.");
      BLP = 0.0f;
    }
    
    if (use_tree_conn) {
      if (!Load_and_analyse_tree_file()) {
	Form1->Memo1->Lines->Add("ERROR: could not successfully load PLULA hierarchy from "+tree_fname);
	use_tree_conn=0;
      }
    }
  }
  
  Form1->Memo1->Lines->Add("");
  Form1->Memo1->Lines->Add("Matrix columns/x dimension: "+IntToStr(xd)+
			   ", rows/y dimension: "+IntToStr(yd));

  if (ADM_set.use_ADMUs) {
    if (Load_and_check_ADMU_files()) {
      Form1->Memo1->Lines->Add("Administrative units successfully loaded.");
      Form1->Memo1->Lines->Add("");
    }
    else {
      Form1->Memo1->Lines->Add("ERRORS loading information for administrative units. ");
      return false;
    }
  }

  if (!load_rest_of_special_files()) {
    //ShowMessage("Error loading cost, mask, solution file or other required files. Exiting.");
    Form1->Memo1->Lines->Add("Error loading cost, mask, solution file or other required files. Exiting.");
    return false;
  }

  // Allocate warp list (after loading PLU layer (if used))
  // In get_next_to_remove() warp_factor will be changed to <= than this one
  if (use_PLULA) {
    wrx.resize(max(warp_factor, (int)glob_plu_max_cnt), 0);
    wry.resize(max(warp_factor, (int)glob_plu_max_cnt), 0);
    wnbc.resize(max(warp_factor, (int)glob_plu_max_cnt), 0);
    wval.resize(max(warp_factor, (int)glob_plu_max_cnt), 0.0f);
  } else {
    wrx.resize(warp_factor, 0);
    wry.resize(warp_factor, 0);
    wnbc.resize(warp_factor, 0);
    wval.resize(warp_factor, 0.0f);
  }

  Form1->Memo1->Lines->Add("");
  if (removal_rule==1)
    {
      Form1->Memo1->Lines->Add("****** REMOVAL RULE: core-area Zonation (CAZ) *******");
    }
  else if (removal_rule==2)
    {
      Form1->Memo1->Lines->Add("****** REMOVAL RULE: Convex additive benefit function (ABF); r^par *******");
    }
  else if (removal_rule==3)
    {
      Form1->Memo1->Lines->Add("****** REMOVAL RULE: Targeting benefit function (TBF) *******");
    }
  else if (removal_rule==4)
    {
      Form1->Memo1->Lines->Add("****** REMOVAL RULE: Generalized benefit function (GBF) *******");
    }
  else
    Form1->Memo1->Lines->Add("****** REMOVAL RULE: Random removal *******");
  Form1->Memo1->Lines->Add("");

  if (!IG_set.use_IG || !IG_set.use_IGw)
    {
      Form1->Memo1->Lines->Add("****** NOT using Info-gap distribution discounting uncertainty analysis ****************");
      Form1->Memo1->Lines->Add("");
    }

  Form1->Memo1->Lines->Add("");
  String drop_opt_name = "\"drop 0 occurrence features\"";
  if (settings_drop_0_occur_layers)
    Form1->Memo1->Lines->Add("*** Attention: enabled "+drop_opt_name+", some features may be dropped (see below).");
  else
    Form1->Memo1->Lines->Add("*** Not using "+drop_opt_name+", all features will be processed even if they do not occur anywhere in the analysis area.");

  if ( !load_precooked_vmat_fn.isEmpty() ){
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("Note: loading precooked vmat ==> no need to allocate any memory for smoothing, interactions, matrix transforms, etc.");
  } else if ( (use_smoothing) || (use_interactions)
	      || (comm_set.load_sim_matrix && comm_set.apply_to_conn) ) {
    String feats = "";
    bool first = true;
    if (use_smoothing)
      feats = feats + "distribution smoothing ";
    if (use_interactions)
      feats = feats + ", interactions ";
    if (comm_set.load_sim_matrix && comm_set.apply_to_conn)
      feats = feats + ", matrix connectivity";	    
    Form1->Memo1->Lines->Add("===============");
    Form1->Memo1->Lines->Add("Note: using the following features which require additional memory for smoothing kernels: " + feats);
    const size_t MB_in_bytes = 1048576;
    String smooth_mem_str = FloatToStrF(mem_required_smoothing(xd, yd)/float(MB_in_bytes), ffFixed, 4,4);
    Form1->Memo1->Lines->Add("Approximately "+smooth_mem_str+" MBs of additional memory are required. Trying to allocate...");

    if (comm_set.load_sim_matrix && comm_set.apply_to_conn) {
      String cst_mem = FloatToStrF(mem_required_mct_smoothing(xd, yd)/MB_in_bytes, ffFixed, 4,4);
      Form1->Memo1->Lines->Add("Note: using 'matrix connectivity transform' with "+
			       IntToStr(comm_set.ct_cnt)+
			       " features, which requires further additional memory: "+cst_mem+" MBs");
    }
	  
    if (Do_smoothings_alloc()) {
      Form1->Memo1->Lines->Add("Memory allocated successfully!");	    
      Form1->Memo1->Lines->Add("===============");
      Form1->Memo1->Lines->Add("");
    } else {
      Form1->Memo1->Lines->Add("Failed to allocate memory for smoothing kernels. Please check that there is enough available memory in your system.");
      Form1->Memo1->Lines->Add("There is no solution for this. Giving up.");
      graceful_exit();
    }
  } else {
    Form1->Memo1->Lines->Add("");
    Form1->Memo1->Lines->Add("****** NOT using distribution smoothing, interactions or matrix connectivity *******");
    Form1->Memo1->Lines->Add("");
  }
  

  if (corr_set.use_corr && !corr_set.layers_fn.isEmpty()) {
    if (!read_corridor_domain_layers_file_n_layers(corr_set)) {
      Form1->Memo1->Lines->Add("Error loading corridor domain layers from file: " + corr_set.layers_fn);
      return false;
    } else {
      Form1->Memo1->Lines->Add("Corridor domain layers successfully loaded from file: " + corr_set.layers_fn);
    }
    Form1->Memo1->Lines->Add("");	  
  }


  if (load_precooked_vmat_fn.isEmpty()) {
    // ** normal load/transformation of biodiversity feature layers
    // buffer for input layers
    float **tmpm = NULL;
    // Count of occurrencies in every cell
    int** occur_count_matrix = imatrix(0, yd, 0, xd);

    // A) Get count of occurrencies to allocate vmat structure
    // (occur_count_matrix is freed asap)
    bool preload_ok = preload_and_count_feature_data_layers(tmpm, occur_count_matrix);
    if (!preload_ok)
      return false;
    // B) Do the real load, fill in vmat from feature layers
    bool transf_ok = load_and_transform_feature_data_layers(tmpm, occur_count_matrix);
    if (!transf_ok)
      return false;

    // Make sure it's freed (but currently it is forced-freed in load_and_transform_features()
    free_matrix(tmpm, 0, yd, 0, xd);
    tmpm = NULL;

  } else {
    // ** loading from vmat dump file
    Form1->Memo1->Lines->Add("===========================================");
    Form1->Memo1->Lines->Add("* Loading vmat, fasten your belt...");
    bool r = load_vmat(load_vmat_directly, ChangeFileExt(ChangeFileExt(outgridfn)), load_precooked_vmat_fn, xd, yd, map_cnt);
    if (!r)
      return false;
    Form1->Memo1->Lines->Add(" * vmat successfully loaded!");
    Form1->Memo1->Lines->Add("===========================================");
    Form1->Memo1->Lines->Add("");
    // fill in spp[]'s prob_sum and IG_fract (like load_from_file_IG normally would do)
    init_output_get_distr_centers_and_sums(true);
  }




  if (random_omissions_percentage_after_vmat_simple_impl > 0.0f) {
    if (random_omissions_max_feature > map_cnt) {
      Form1->Memo1->Lines->Add("");
      Form1->Memo1->Lines->Add(" * ERROR in settings: the random omissions max feature number ("+
			       IntToStr(random_omissions_max_feature)+") is bigger than the total number of features ("+
			       IntToStr(map_cnt)+")");
      Form1->Memo1->Lines->Add("");
    } else {
      Form1->Memo1->Lines->Add("");
      Form1->Memo1->Lines->Add("===========================================");
      Form1->Memo1->Lines->Add(" * Generating random omissions, percentage: "+
			       FloatToStr(random_omissions_percentage_after_vmat_simple_impl));
      Form1->Memo1->Lines->Add("   (up to feature "+IntToStr(random_omissions_max_feature)+")");
      gen_random_omissions_percentage_simple_impl(random_omissions_percentage_after_vmat_simple_impl);
      Form1->Memo1->Lines->Add("===========================================");
      Form1->Memo1->Lines->Add("");
    }
  }

  if (!save_vmat_fn.isEmpty()) {
    Form1->Memo1->Lines->Add("===========================================");
    Form1->Memo1->Lines->Add("* Saving vmat, this can take long...");
    // TODO: there we go...
    // open and save
    bool r = save_vmat(ChangeFileExt(ChangeFileExt(outgridfn)), save_vmat_fn, xd, yd, map_cnt);
    if (!r)
      return false;
    Form1->Memo1->Lines->Add(" * vmat successfully saved!");
    Form1->Memo1->Lines->Add("===========================================");
    Form1->Memo1->Lines->Add("");
  }


  if (use_SSI)
    {
      SSI_load();
    }

  if (use_groups) {
    glbl_groups_info.grpv = new int[2*map_cnt];
    glbl_groups_info.max_grp = get_grp_max(glbl_groups_info.grpv);
  } else {
    glbl_groups_info.grpv = NULL;
    glbl_groups_info.max_grp = 0;	  
  }

  if (ADM_set.use_ADMUs) {
    get_ADMU_data();
    ADMUs_output_1();
    // intra-ADMU curves
    if (ADM_set.row_count_for_per_admu_curves > 0) {
      ADMUs_output_per_admu_curves_init();
    }
  }

  // intra-ADMU groups
  if (ADM_set.use_ADMUs && use_groups) {
    ADMUs_output_per_admu_grp_curves_init(glbl_groups_info);
  }

  // not needed anymore once all the tables and restrutures have been shifted
  if (settings_drop_0_occur_layers)
    drop_0_relocation_idx.resize(0);

  return true;
}
