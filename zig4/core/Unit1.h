#ifndef Unit1H
#define Unit1H

#include "VCL.h"

class TForm1 : public TForm
{
public:
	TPageControl *PageControl1;
	TTabSheet *TabSheet1;
	TTabSheet *TabSheet2;
	TTabSheet *TabSheet3;
	TMemo *Memo1;
	TImage *Image1;
	TButton *Button3;
	TButton *Button4;
	TButton *Button5;
	TButton *Button6;
	TEdit *Edit2;
	TOpenDialog *OpenDialog1;
	TTabSheet *TabSheet4;
	TLabel *Label1;
	TLabel *Label2;
	TBevel *Bevel1;
	TLabel *Label4;
	TLabel *Label22;
	TLabel *Label32;
	TLabel *Label5;
	TLabel *Label6;
	TLabel *Label20;
	TLabel *Label21;
	TComboBox *ComboBox1;
	TLabel *Label3;
	TChart *Chart1;
	TLineSeries *Series1;
	TLineSeries *Series2;
	TPanel *Panel1;
	TImage *Image3;
	TSaveDialog *SaveDialog1;
	TEdit *Edit3;
	TTabSheet *TabSheet5;
	TLabel *Label13;
	TTabSheet *TabSheet6;
	TLabel *Label15;
	TLabel *Label16;
	TCheckBox *CheckBox2;
	TButton *Button8;
	TEdit *Edit10;
	TEdit *Edit11;
	TEdit *Edit12;
	TLabel *Label18;
	TEdit *Edit13;
	TComboBox *ComboBox2;
	TUpDown *UpDown1;
	TTabSheet *TabSheet7;
	TLabel *Label37;
	TEdit *Edit25;
	TLabel *Label39;
	TCheckBox *CheckBox3;
	TCheckBox *CheckBox4;
	TEdit *Edit27;
	TCheckBox *CheckBox11;
	TGroupBox *GroupBox2;
	TCheckBox *CheckBox12;
	TLabel *Label23;
	TEdit *Edit15;
	TRadioGroup *RadioGroup1;
	TCheckBox *CheckBox13;
	TEdit *Edit16;
	TCheckBox *CheckBox14;
	TRadioGroup *RadioGroup2;
	TEdit *Edit17;
	TGroupBox *GroupBox1;
	TLabel *Label14;
	TLabel *Label19;
	TCheckBox *CheckBox8;
	TEdit *Edit9;
	TCheckBox *CheckBox9;
	TCheckBox *CheckBox10;
	TEdit *Edit14;
	TEdit *Edit18;
	TLabel *Label25;
	TLabel *Label26;
	TLabel *Label27;
	TLabel *Label28;
	TLabel *Label29;
	TLabel *Label30;
	TEdit *Edit19;
	TChart *Chart2;
	TLineSeries *LineSeries1;
	TLineSeries *LineSeries2;
	TButton *Button9;
	TRadioGroup *RadioGroup3;
	TLineSeries *Series3;

	// fedemp: The "weighted average" curve, sister of Series2 (min) and Series3 (avg) - rep. vs. prop remaining.
	TLineSeries *Series_23;	

	TLineSeries *Series4;

	// fedemp: The "weighted average" curve, sister of LineSeries2 (min) and Series4 (avg) - val vs. cost
	TLineSeries *LineSeries2_Series4;	

	TLabel *Label45;
	TGroupBox *GroupBox3;
	TLabel *Label40;
	TEdit *Edit20;
	TLabel *Label7;
	TEdit *Edit21;
	TButton *Button2;
	TLabel *Label41;
	TEdit *Edit26;
	TGroupBox *GroupBox4;
	TLabel *Label33;
	TEdit *Edit22;
	TLabel *Label34;
	TEdit *Edit24;
	TLabel *Label35;
	TEdit *Edit23;
	TButton *Button1;
	TLabel *Label44;
	TEdit *Edit29;
	TLabel *Label36;
	TCheckBox *CheckBox7;
	TLabel *Label31;
	TProgressBar *ProgressBar1;
	TButton *Button10;
	TLabeledEdit *LabeledEdit1;
	TLabel *Label47;
	TCheckBox *CheckBox1;
	TLabeledEdit *LabeledEdit2;
	TSaveDialog *SaveDialog2;
	TChart *Chart3;
	TLineSeries *LineSeries3;
	TLineSeries *LineSeries4;
	TLabel *Label48;
	TLabel *Label49;
	TLabel *Label51;
	TImage *Image4;
	TCheckBox *CheckBox5;
	TEdit *Edit1;
	TChart *Chart4;
	TLineSeries *LineSeries5;
	TLineSeries *LineSeries6;

	// fedemp: The "weighted average" curve, sister of LineSeries5 (min) and LineSeries6 (avg) - SSI rep. vs. prop remaining.
	TLineSeries *LineSeries_56;	

	TLineSeries *LineSeries7;
	TRadioGroup *RadioGroup4;
	TSavePictureDialog *SavePictureDialog1;
	TLabel *Label38;
	TCheckBox *CheckBox6;
	TEdit *Edit30;
	TLabel *Label46;
	TEdit *Edit31;
	TChart *Chart5;
	TBarSeries *Series5;
	TLabel *Label50;
	TLabel *Label52;
	TLabel *Label53;
	TLabel *Label54;
	TLabel *Label55;
	TLabel *Label56;
	TLabel *Label57;
	TEdit *Edit32;
	TEdit *Edit33;
	TEdit *Edit34;
	TEdit *Edit35;
	TEdit *Edit37;
	TEdit *Edit38;
	TEdit *Edit39;
	TButton *Button11;
	TCheckBox *CheckBox15;
	TEdit *Edit40;
	TCheckBox *CheckBox16;
	TEdit *Edit41;
	TLabel *Label42;
	TLabel *Label59;
	TLabel *Label60;
	TLabel *Label61;
	TLabel *Label17;
	TButton *Button13;
	TButton *Button14;
	TButton *Button15;
	TButton *Button16;
	TButton *Button17;
	TButton *Button18;
	TButton *Button19;
	TButton *Button20;
	TLabel *Label24;
	TLabel *Label62;
	TGroupBox *GroupBox5;
	TLabel *Label58;
	TEdit *Edit36;
	TButton *Button12;
	TGroupBox *GroupBox6;
	TLabel *Label8;
	TEdit *Edit4;
	TLabel *Label10;
	TEdit *Edit5;
	TLabel *Label9;
	TEdit *Edit6;
	TLabel *Label11;
	TEdit *Edit7;
	TLabel *Label12;
	TEdit *Edit8;
	TLabel *Label43;
	TEdit *Edit28;
	TButton *Button7;
	TRadioGroup *RadioGroup5;
	TLabel *Label63;
	TEdit *Edit42;

	void FormClose(TObject *Sender, TCloseAction &Action);

#if 0
	void Button3Click(TObject *Sender);
	void Button4Click(TObject *Sender);
	void Button5Click(TObject *Sender);
	void Image1DblClick(TObject *Sender);
	void Button6Click(TObject *Sender);
	void ComboBox1Change(TObject *Sender);
	void Chart1DblClick(TObject *Sender);
	void Button7Click(TObject *Sender);
	void Label13Click(TObject *Sender);
	void ComboBox2Change(TObject *Sender);
	void UpDown1ChangingEx(TObject *Sender,
			       bool &AllowChange, short NewValue, TUpDownDirection Direction);
	void TabSheet4Show(TObject *Sender);
	void TabSheet4Enter(TObject *Sender);
	void PageControl1Enter(TObject *Sender);
	void TabSheet1Enter(TObject *Sender);
	void TabSheet2Enter(TObject *Sender);
	void TabSheet3Enter(TObject *Sender);
	void TabSheet5Enter(TObject *Sender);
	void TabSheet6Enter(TObject *Sender);
	void PageControl1Change(TObject *Sender);
	void Button1Click(TObject *Sender);
	void Button2Click(TObject *Sender);
	void Button9Click(TObject *Sender);
	void Button10Click(TObject *Sender);
	void Label47Click(TObject *Sender);
	void Label29Click(TObject *Sender);
	void Chart2DblClick(TObject *Sender);
	void Chart3DblClick(TObject *Sender);
	void Chart5DblClick(TObject *Sender);
	void Button11Click(TObject *Sender);
	void Button12Click(TObject *Sender);
	void Chart4DblClick(TObject *Sender);
	void Button13Click(TObject *Sender);
	void Button14Click(TObject *Sender);
	void Button15Click(TObject *Sender);
	void Button16Click(TObject *Sender);
	void Button17Click(TObject *Sender);
	void Button18Click(TObject *Sender);
	void Button19Click(TObject *Sender);
	void Button20Click(TObject *Sender);
#endif
private:	// User declarations
public:		// User declarations
	TForm1(TComponent* Owner);
};

extern String setfn, sppfn, outgridfn, outgridfn2, outgridfn3, emffn, costfn, features_info_fn, curvesfn, maskfn, BQP_prof_fn, memofn, SSI_fname, SSI_features_info_fn, SSI_curvesfn, PLULAfn, tree_fname, ia_fname, condition_fname, retention_fname, groups_fname, grp_curves_fn, area_mask_fname, LEC_fname;
extern String admu_redist_outgridfn;
extern String output_transf_ig_layers_dn, output_transf_ds_layers_dn, output_transf_cst_layers_dn, 
  output_transf_ct_layers_dn, output_transf_rt_layers_dn, output_transf_mct_layers_dn, 
  output_transf_ia_layers_dn, 
  output_transf_final_layers_dn;
extern bool   autoclose, run_finished, terminate_flag, neg_weights_used;
//float alpha, weight;
extern float rem_level, val_th, BLP, SA_z;
extern int   normalize, smooth, show_images, miss_as_zero, use_cost, use_mask, use_BQP, add_edge_cnt, removal_rule, use_SSI, use_tree_conn, use_interactions, mem_save_mode, use_condition, use_retention, use_groups, mask_data, use_LEC;
extern bool glob_add_to_edge_borders_between_removal_mask_levels;

// arbitrary kernels
extern bool use_arb_kernels;
extern int arb_kernels_default;
extern float arb_kernels_constant;
extern String arb_kernels_prefix;

// corridors
extern struct Corr_settings corr_set;

// connected components split penalty
extern float CCSP_prev_removal_threshold, CCSP_prev_prev_removal_threshold, CCSP_ma_removal_rate;
extern int CCSP_verbose;
extern float CCSP;
extern FILE* wlist_file;
extern bool use_8_connectivity;
extern String output_ccsp_fn;
// TODO: temporary, let's see what this becomes in the end.
// period in removal iterations
extern int CCSP_formula;
extern size_t CCSP_info_period;
extern int CCSP_variant, CCSP_thickness;
class Grid_CCSP;
extern Grid_CCSP* cc_accounting;

extern std::vector<float> distrib_centers_x;
extern std::vector<float> distrib_centers_y;

// outputs enable/disable
extern bool glob_set_output_wrscr_map;
extern bool glob_set_output_prop_rank_map;

extern class GridMap obsmap[], wmap[], costmap, obs_map, maskmap, tstmap;
extern struct sp* spp;

// Map of occurrence weights
extern String glob_set_occur_size_weights_layer_fn;
extern class GridMap occur_size_weights_map;
extern bool use_occur_size_weights_correct_landscape_fraction;
extern bool use_occur_size_weights_correct_cost;
extern bool use_occur_size_weights_correct_ranking;

extern float get_cell_area_correction_factor(int x, int y);

extern struct Tgroups_info glbl_groups_info;

/// number of species layers
extern int map_cnt;
// nonm1: "non-missing" cells  / removable cells, initialized in Calc_richness_et_al_matrixes()
// removed: counter for the # of cells removed so far (in Remove_site())
// current: seems redundant with 'removed' - ++ immediately after Remove_site()   (and current == removed+1)
//          the difference seems to be in their range ([0,...] vs. [1,...]) and the scope where they are used
extern int    IG_map_cnt, ia_cnt, cur_map, nonm1, m1s, removed, current;

extern float  **Rmax, **Rsum, **sol, **sol_val, **curves;
extern float **SSI_indiv_curves;
extern bool forget_Rsum;
extern char   **edge, **status;
extern float  Rmax_max, Rsum_max; // Rmin_max, **Rmin, **Rave, Rave_max
extern int    *exl, *eyl, ecnt;
extern int **mat_edge_list_pos;

extern const int FIRST_NUM_ARB_KERNEL;
typedef std::pair<std::vector<float>, std::vector<float> > Arb_Kernel;
typedef std::map<int, Arb_Kernel > Arbitrary_Kernels_Map;
extern Arbitrary_Kernels_Map arbitrary_kernels;

int sortfunc( const void *a, const void *b);
int sortfunc_with_hierarchy( const void *a, const void *b);

float z_pow(float x, float ex);

// to replace the potentially very slow std powf() implementation
#define USE_POW_APPROX_IEEE_754 1
const int powf_fast_precision = 18;
void powf_fast_init(int prec=powf_fast_precision);
void powf_fast_release();
float powf_fast(const float x, const float y);
//

extern "C" void zig_legacy_alloc_error(char error_text[]);

#endif
