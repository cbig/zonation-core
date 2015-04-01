#ifndef ZIG4POD_H
#define ZIG4POD_H

#include <QString>
#include <QList>
#include <QVector>
#include <QSet>
#include <QHash>
#include <boost/shared_ptr.hpp>
#include <boost/enum.hpp>
#include <boost/optional.hpp>
#include <boost/preprocessor.hpp>
#include <boost/variant.hpp>
#include <boost/multi_array.hpp>
#include "filepath.h"
#include "zig4lib_global.h" 

struct ZUnknownFile
{
};

//--------- Raster ---------

struct ZRasterFile
{
	int xsize;
	int ysize;
	QString projection;
	double transform[6];
	ZRasterFile() :
	        xsize(0),
	        ysize(0)
	{}
};

//--------- SPP ---------

BOOST_ENUM(SpeciesFormat,
           (OldFormat)
           (NewFormat)
           )

#define Z_SPP_VARIABLES \
	((double)(weight)(weight_)) /* > 0, first column */ \
	((double)(alpha)(alpha_)) /* > 0, second column */ \
	((int)(bqp_curve_index)(column3_)) /* third column */ \
	((int)(nqp_upstream_curve_index)(column3_)) /* third column */ \
	((int)(bqp_buffer_size)(column4_)) /* fourth column */ \
	((int)(nqp_downstream_curve_index)(column4_)) /* fourth column */ \
	((double)(additive_rule_exponent)(column5_)) /* 5th | > 0 */ \
	((double)(target_rule_proportion)(column5_)) /* 5th | [0, 1] */ \
	((double)(generalized_rule_weight)(column5_)) /* 5th | > 0 */ \
	((double)(generalized_rule_target)(generalized_rule_target_)) /* 6th | [0, 1] */ \
	((double)(generalized_rule_exp_x)(generalized_rule_exp_x_)) /* 7th */ \
	((double)(generalized_rule_exp_y)(generalized_rule_exp_y_)) /* 8th */ \
	((FilePath)(raster)(raster_)) /*  */

#define Z_SPP_VARIABLE_DEFS_REPEAT(r, data, elem) BOOST_PP_SEQ_ELEM(0, elem) & BOOST_PP_SEQ_ELEM(1, elem)() { return BOOST_PP_SEQ_ELEM(2, elem); }
//#define Z_SPP_VARIABLE_INIT_REPEAT(r, data, elem) (BOOST_PP_SEQ_ELEM(1, elem)(BOOST_PP_SEQ_ELEM(2, elem)))
//#define Z_SPP_VARIABLE_INIT_LIST BOOST_PP_SEQ_FOR_EACH(Z_SPP_VARIABLE_INIT_REPEAT, ~, Z_SPP_VARIABLES)

struct ZSPPEntry
{
	double weight_;
	double alpha_;
	int column3_;
	int column4_;
	double column5_;
	double generalized_rule_target_; // optional
	double generalized_rule_exp_x_; // optional
	double generalized_rule_exp_y_; // optional
	FilePath raster_;

	// accessor functions

	BOOST_PP_SEQ_FOR_EACH(Z_SPP_VARIABLE_DEFS_REPEAT, ~, Z_SPP_VARIABLES)

	                ZSPPEntry() :
	        weight_(1.0),
	        alpha_(1.0),
	        column3_(1),
	        column4_(1),
	        column5_(1.0),
	        generalized_rule_target_(1.0),
	        generalized_rule_exp_x_(1.0),
	        generalized_rule_exp_y_(1.0)
	{}
};

struct ZSPPFile
{
	SpeciesFormat format;
	QVector<boost::shared_ptr<ZSPPEntry> > sppList;
};

//--------- BQP ---------

struct ZBQPPoint
{
	double x;
	double y;
	ZBQPPoint() :
	        x(0.0),
	        y(0.0)
	{}
};

struct ZBQPCurve
{
	QVector<ZBQPPoint> pointList;
};

struct ZBQPFile
{
	QVector<boost::shared_ptr<ZBQPCurve> > curveList;
};

//--------- SSICoordinates ---------

struct ZSSICoordinates
{
	double x;
	double y;
	double value;
	double error;
	ZSSICoordinates() : x(0.0), y(0.0), value(0.0), error(0.0) {}
};

struct ZSSICoordinatesFile
{
	QVector<ZSSICoordinates> coordinatesList;
};

//--------- SSI ---------

#define Z_SSI_VARIABLES \
	((double)(weight)(weight_)) /* > 0, first column */ \
	((double)(additive_rule_exponent)(column5_)) /* 5th | > 0 */ \
	((double)(target_rule_proportion)(column5_)) /* 5th | [0, 1] */ \
	((double)(generalized_rule_weight)(column5_)) /* 5th | > 0 */ \
	((double)(generalized_rule_target)(generalized_rule_target_)) /* 6th | [0, 1] */ \
	((double)(generalized_rule_exp_x)(generalized_rule_exp_x_)) /* 7th */ \
	((double)(generalized_rule_exp_y)(generalized_rule_exp_y_)) /* 8th */ \
	((FilePath)(ssi)(ssi_)) /*  */

#define Z_SSI_VARIABLE_DEFS_REPEAT(r, data, elem) BOOST_PP_SEQ_ELEM(0, elem) & BOOST_PP_SEQ_ELEM(1, elem)() { return BOOST_PP_SEQ_ELEM(2, elem); }

struct ZSSIEntry
{
	double weight_;
	double column5_;
	double generalized_rule_target_; // optional
	double generalized_rule_exp_x_; // optional
	double generalized_rule_exp_y_; // optional
	FilePath ssi_;

	// accessor functions

	BOOST_PP_SEQ_FOR_EACH(Z_SSI_VARIABLE_DEFS_REPEAT, ~, Z_SSI_VARIABLES)

	ZSSIEntry() :
	        weight_(1.0),
	        column5_(1.0),
	        generalized_rule_target_(1.0),
	        generalized_rule_exp_x_(1.0),
	        generalized_rule_exp_y_(1.0)
	{}
};

struct ZSSIFile
{
	SpeciesFormat format;
	QVector<boost::shared_ptr<ZSSIEntry> > ssiList;
};

//--------- PPA ---------

struct ZLSI
{
	double f1; // percentage of lanscape [0-100]
	double f2; // inclusion minimum percentage [0-100]
	double d; // maximum distance [0-]
	double sim; // maximum difference in species composition [0-]
	ZLSI() :
	        f1(0.0),
	        f2(0.0),
	        d(0.0),
	        sim(0.0)
	{}
};

struct ZLSC
{
	double f1; // fraction of present solution [0-1]
	double f2; // fraction of comparison solution [0-1]
	FilePath comp_solution;
	FilePath output_file;
	ZLSC() :
	        f1(0.0),
	        f2(0.0)
	{}
};

// LSI to mask
struct ZLSM
{
	FilePath mask_file;
	double f2; // inclusion minimum percentage [0-100]
	double d; // maximum distance [0-]
	double sim; // maximum difference in species composition [0-]
	ZLSM() :
	        f2(0.0),
	        d(0.0),
	        sim(0.0)
	{}
};

// LSI to top fraction of a mask
struct ZLSB
{
	FilePath mask_file;
	double f1; // percentage of mask lanscape [0-100]
	double f2; // inclusion minimum percentage [0-100]
	double d; // maximum distance [0-]
	double sim; // maximum difference in species composition [0-]
	ZLSB() :
	        f1(0.0),
	        f2(0.0),
	        d(0.0),
	        sim(0.0)
	{}
};

typedef boost::variant<ZLSI, ZLSC, ZLSM, ZLSB> ZAnalysis;

struct ZPPAFile
{
	QVector<boost::shared_ptr<ZAnalysis> > analysisList;
};

//--------- Groups ---------

// loaddata.cpp: Apply_condition_matrix()
BOOST_ENUM_VALUES(RetentionMode, int,
                  (Loss)(1)
                  (Gain)(2)
		  (Unspecified)(-1)
                  )

struct ZGroupsEntry
{
	int num; // output group
	int cond_num; // condition group
	int ret_num; // retention group
	RetentionMode ret_mode; // retention mode
  // int LEC_num; // local edge correction group - deprecated
  //int arb_kernel_num;
  int arb_kernel_ds;
  int arb_kernel_matrix;
  int arb_kernel_ia;

  ZGroupsEntry():
     num(-1),
     cond_num(-1),
     ret_num(-1),
     ret_mode(RetentionMode::Loss),
     //LEC_num(-1) // always in version 3  - LEC_num deprecated
     arb_kernel_ds(-1),
     arb_kernel_matrix(-1),
     arb_kernel_ia(-1)
     {}
};

struct ZGroupsFile
{
	QVector<boost::shared_ptr<ZGroupsEntry> > groupsList;
};

//--------- Tree ---------

struct ZTreeEntry
{
	int id;
	int down_id;
	ZTreeEntry() :
	        id(-1),
	        down_id(-1)
	{}
};

struct ZTreeFile
{
	QVector<boost::shared_ptr<ZTreeEntry> > treeList;
};

//--------- Interactions ---------

BOOST_ENUM_VALUES(InteractionsMode, int,
                  (Positive)(1)
                  (Negative)(2)
                  )

struct ZInteractionsEntry
{
	int l1;
	int l2;
	double alpha;
	InteractionsMode iatype;
	double gamma;
	ZInteractionsEntry() :
	        l1(1),
	        l2(1),
	        alpha(1.0),
	        iatype(InteractionsMode::Positive),
	        gamma(1.0)
	{}
};

struct ZInteractionsFile
{
	QVector<boost::shared_ptr<ZInteractionsEntry> > interactionsList;
};

//--------- Matrix ---------

struct ZMatrixFile
{
	typedef boost::multi_array<double, 2> Type;
	Type matrix;
};

//--------- ADMU ---------

struct ZADMUEntry
{
	int fid; // identification number
	double glob_weight; // global weight
	double loc_weight; // local weight
	QString unit_name; // admu name
	ZADMUEntry() :
	        fid(0),
	        glob_weight(1.0),
	        loc_weight(1.0)
	{}
};

struct ZADMUFile
{
	QVector<boost::shared_ptr<ZADMUEntry> > admuList;
};

//--------- IndexedRasters ---------
// this is used to read retention, condition and lec lists

struct ZIndexedRastersEntry
{
	int index;
	FilePath raster;
	ZIndexedRastersEntry() : index(0) {}
};

struct ZIndexedRastersFile
{
	QVector<boost::shared_ptr<ZIndexedRastersEntry> > rasterList;
};

//--------- IG ---------

struct ZIGEntry
{
	double weight;
	FilePath raster;
	ZIGEntry() : weight(1.0) {}
};

struct ZIGFile
{
	QVector<boost::shared_ptr<ZIGEntry> > igList;
};

//--------- DAT ---------

BOOST_ENUM_VALUES(BQPMode, int,
                  (NoBQP)(0)
                  (UniformBQP)(1)
                  (NonUniformBQP)(2)
                  )

BOOST_ENUM_VALUES(RemovalRule, int,
                  (CoreRemoval)(1)
                  (AdditiveRemoval)(2)
                  (TargetRemoval)(3)
                  (GeneralizedRemoval)(4)
                  (UnknownRemoval)(5)
                  )

BOOST_ENUM_VALUES(ADMUMode, int,
                  (FirstADMU)(1)
                  (SecondADMU)(2)
                  )

#define Z_DAT_VARIABLES \
	((double)(rem_level)("Settings")("initial removal percent")(0.0)) \
	((bool)(use_cost)("Settings")("use cost")(false)) \
	((bool)(use_mask)("Settings")("use mask")(false)) \
	((bool)(use_BQP)("Settings")("use boundary quality penalty")(false)) \
	((BQPMode)(BQP_mode)("Settings")("BQP mode")(BQPMode::NonUniformBQP)) \
	((double)(BLP)("Settings")("BLP")(0.0)) \
	((double)(SA_z)("Settings")("z")(0.25)) \
	((bool)(use_edge_rem)("Settings")("edge removal")(true)) \
	((bool)(add_to_edge_borders_between_removal_mask_levels)("Settings")("add to edge borders between removal mask levels")(true)) \
	((bool)(annotate_name)("Settings")("annotate name")(true)) \
	((int)(resample_count)("Settings")("resample species")(0)) \
	((bool)(drop_0_occur_layers)("Settings")("drop 0 occurrence layers")(0)) \
	((RemovalRule)(removal_rule)("Settings")("removal rule")(RemovalRule::CoreRemoval)) \
	((FilePath)(costfn)("Settings")("cost file")("")) \
	((FilePath)(maskfn)("Settings")("mask file")("")) \
	((FilePath)(BQP_prof_fn)("Settings")("BQP profiles file")("")) \
	((int)(warp_factor)("Settings")("warp factor")(1)) \
	((bool)(logit)("Settings")("logit space")(false)) \
	((int)(add_edge_cnt)("Settings")("add edge points")(0)) \
	((bool)(use_SSI)("Settings")("use SSI")(false)) \
	((FilePath)(SSI_fname)("Settings")("SSI file name")("")) \
	((FilePath)(PPA_fname)("Settings")("post-processing list file")("")) \
	((bool)(all_zero_as_miss)("Settings")("treat zero-areas as missing data")(false)) \
	((bool)(mem_save_mode)("Settings")("memory save mode")(false)) \
	((bool)(use_groups)("Settings")("use groups")(false)) \
	((FilePath)(groups_fname)("Settings")("groups file")("")) \
	((bool)(use_condition)("Settings")("use condition layer")(false)) \
	((FilePath)(condition_fname)("Settings")("condition file")("")) \
	((bool)(use_retention)("Settings")("use retention layer")(false)) \
	((FilePath)(retention_fname)("Settings")("retention file")("")) \
	((bool)(use_LEC)("Settings")("use local edge corrections")(false)) \
	((FilePath)(LEC_fname)("Settings")("local edge corrections file")("")) \
	((double)(ret_layers_rel_weight)("Settings")("retention layers relative weight")(1.0)) \
	((bool)(mask_data)("Settings")("mask missing areas")(false)) \
	((FilePath)(area_mask_fname)("Settings")("area mask file")("")) \
	((bool)(use_PLULA)("Settings")("use planning unit layer")(false)) \
	((FilePath)(PLULAfn)("Settings")("planning unit layer file")("")) \
	((bool)(use_tree_conn)("Settings")("use tree connectivity")(false)) \
	((FilePath)(tree_fname)("Settings")("tree connectivity file")("")) \
	((bool)(use_interactions)("Settings")("use interactions")(false)) \
	((FilePath)(ia_fname)("Settings")("interaction file")("")) \
	((double)(random_omissions_percentage_after_vmat_simple_impl)("Settings")("random omission percentage after vmat simple impl")(0.0)) \
	((int)(random_omissions_max_feature)("Settings")("random omission max feature")(0.0)) \
	((bool)(load_sim_matrix)("Community analysis settings")("load similarity matrix")(false)) \
	((FilePath)(sim_matrix_name)("Community analysis settings")("connectivity similarity matrix file")("")) \
	((bool)(apply_to_conn)("Community analysis settings")("apply to connectivity")(false)) \
	((FilePath)(comm_sim_matrix_name )("Community analysis settings")("community similarity matrix file")("")) \
	((bool)(apply_to_repr)("Community analysis settings")("apply to representation")(false)) \
	((FilePath)(edge_effect_fn)("Community analysis settings")("connectivity edge effect fix file")("")) \
	((bool)(use_ADMUs)("Administrative units")("use ADMUs")(false)) \
	((ADMUMode)(ADMU_mode)("Administrative units")("ADMU mode")(ADMUMode::FirstADMU)) \
	((bool)(calc_from_condition)("Administrative units")("calculate local weights from condition")(false)) \
	((FilePath)(ADM_layer_file)("Administrative units")("ADMU layer file")("")) \
	((FilePath)(ADM_weights_file)("Administrative units")("ADMU descriptions file")("")) \
	((FilePath)(ADM_weight_matrix_file)("Administrative units")("ADMU weight matrix")("")) \
	((double)(mode2_global_weight)("Administrative units")("Mode 2 global weight")(0.0)) \
	((int)(ADM_row_count_for_per_admu_curves)("Administrative units")("row count for per ADMU output curves")(0.0)) \
	((bool)(IG_proportional)("Info-gap settings")("Info-gap proportional")(false)) \
	((bool)(use_IGw)("Info-gap settings")("use info-gap weights")(false)) \
	((FilePath)(IGwfn)("Info-gap settings")("Info-gap weights file")("")) \
	((bool)(normalize_IGw)("Info-gap settings")("normalize info-gap weights")(false)) \
	((bool)(output_final_layers)("Transformed layers")("output final transformed layers")(false)) \
	((bool)(output_ig_layers)("Transformed layers")("output info-gap transformed layers")(false)) \
	((bool)(output_ds_layers)("Transformed layers")("output distribution smoothing transformed layers")(false)) \
	((bool)(output_cst_layers)("Transformed layers")("output community similarity transformed layers")(false)) \
	((bool)(output_ct_layers)("Transformed layers")("output condition transformed layers")(false)) \
	((bool)(output_rt_layers)("Transformed layers")("output retention transformed layers")(false)) \
	((bool)(output_mct_layers)("Transformed layers")("output matrix connectivity transformed layers")(false)) \
	((bool)(output_ia_layers)("Transformed layers")("output interactions transformed layers")(false)) \
	((bool)(arb_kernels_use)("Arbitrary kernels")("use arbitrary kernels")(false)) \
	((double)(arb_kernels_constant)("Arbitrary kernels")("arbitrary kernels constant")(0.0)) \
	((int)(arb_kernels_default)("Arbitrary kernels")("default kernel")(-1)) \
	((FilePath)(arb_kernels_prefix)("Arbitrary kernels")("arbitrary kernels files prefix")("")) \
	\
	((bool)(use_corridors)("Corridor loss penalty")("use corridors")(false)) \
	((double)(CCSP)("Corridor loss penalty")("strength")(0.0))\
	((int)(CCSP_thickness)("Corridor loss penalty")("minimum width")(7))\
	((int)(CCSP_formula)("Corridor loss penalty")("penalty formula")(1))\
	((int)(CCSP_variant)("Corridor loss penalty")("variant")(1))\
	((int)(CCSP_info_period)("Corridor loss penalty")("info period")(10000))\
	((QString)(corridor_output_boundaries_pcs)("Corridor loss penalty")("output corridor boundaries layer at top percentage")("")) \
	((bool)(use_corridor_domain_layers)("Corridor loss penalty")("use domain layers")(false)) \
	((FilePath)(corridor_domain_layers_fn)("Corridor loss penalty")("domain layers list file")("")) \
	((double)(corridor_start_pc)("Corridor loss penalty")("start at top percentage")(100.0)) \
	((int)(CCSP_bfs_breadth_x)("Corridor loss penalty")("redundancy radius x")(0))\
	((int)(CCSP_bfs_breadth_y)("Corridor loss penalty")("redundancy radius y")(0))\
	\
	((FilePath)(save_vmat)("vmat")("save vmat")("")) \
	((FilePath)(load_precooked_vmat)("vmat")("load vmat")("")) \
	((bool)(load_vmat_directly)("vmat")("load vmat directly and risky")(false)) \
	\
	((FilePath)(occur_size_weights_layer_fn)("Settings")("cell area layer")("")) \
	\
	((bool)(output_richness)("Outputs")("output weighted range size corrected richness")(true)) \
	((bool)(output_occur_richness_map)("Outputs")("output occurrence richness map")(false)) \
	((bool)(output_w_occur_richness_map)("Outputs")("output weighted occurrence richness map")(false)) \
	((bool)(output_occur_level_richness_rank)("Outputs")("output occurrence level richness map")(false)) \
	((bool)(output_w_occur_level_richness_rank)("Outputs")("output weighted occurrence level richness map")(false)) \
	((bool)(output_prop_rank)("Outputs")("output proportional loss ranking")(false)) \
	((bool)(output_stop_before_ranking)("Outputs")("stop before ranking")(false)) \
	/**/
#define Z_DAT_VARIABLE_DEFS_REPEAT(r, data, elem) boost::optional<BOOST_PP_SEQ_ELEM(0, elem)> BOOST_PP_SEQ_ELEM(1, elem);

struct ZDATFile
{
	BOOST_PP_SEQ_FOR_EACH(Z_DAT_VARIABLE_DEFS_REPEAT, ~, Z_DAT_VARIABLES)
};

class FormatInfo : public QPair<QString, QString>
{
public:
	typedef QPair<QString, QString> Base;
	FormatInfo(QString const& parameter, QString const& suffix) :
	        Base(parameter, suffix) {}
	FormatInfo(QString const& parameter) :
	        Base(parameter, parameter) {}
	QString parameter() const { return first; }
	QString suffix() const { return second; }
	QString extension() const { 
	  if (suffix() == "compressed-tif") 
	    return ".compressed.tif";
	  else if (suffix() == "compressed-img") 
	    return ".compressed.img";
	  else 
	    return '.' + suffix();
	}
	// conversion to QString returns the parameter (now we're just being lazy)
	operator QString() const { return parameter(); }
};

BOOST_ENUM_VALUES(ZGridFormat, FormatInfo,
                  (ASC)("asc")
                  (TIF)("tif")
		  //                  (CTIF)(FormatInfo("compressed-tif","compressed.tif")
                  (CTIF)("compressed-tif")
                  (IMG)("img")
		  //                  (CIMG)(FormatInfo("compressed-img","compressed.img"))
                  (CIMG)("compressed-img")
                  )

BOOST_ENUM_VALUES(ZImageFormat, FormatInfo,
                  (PNG)("png")
                  (BMP)("bmp")
                  (JPG)("jpg")
                  (EMF)("emf")
                  )

//--------- BAT ---------

BOOST_ENUM_VALUES(ZCommandLineMode, QString,
                  (CalculateRank)("-r")
                  (LoadRank)("-l")
                  )

struct ZCommandLine
{
  QString line;   // text line as it is typed in the shell or .bat/project file (useful for reporting)
  DirPath workingDir;
  FilePath zigexec;
  ZCommandLineMode commandLineMode;
  FilePath rankFile;
  FilePath datFile;
  FilePath sppFile;
  FilePath outFile;
  double uncertaintyAlpha;
  bool distributionSmoothingOn;
  double dispersalKernelMultiplier;
  bool windowLeftOpen;
  boost::optional<int> numberOfThreads;
  boost::optional<int> removalRule;
  boost::optional<int> warpFactor;
  boost::optional<QSet<ZGridFormat> > gridOutputFormats;
  boost::optional<QSet<ZImageFormat> > imageOutputFormats;
};

struct ZBATFile
{
	QVector<boost::shared_ptr<ZCommandLine> > commandList;
};

#endif // ZIG4POD_H
