#include "zconf.h"
#include "enum_util.h"

#define Z_DAT_CONF_DEFS_REPEAT(r, data, elem)\
BOOST_PP_SEQ_ELEM(0, elem) ZConf::BOOST_PP_SEQ_ELEM(1, elem)() const\
{ return boost::get_optional_value_or(dat.BOOST_PP_SEQ_ELEM(1, elem), BOOST_PP_SEQ_ELEM(4, elem)); }

BOOST_PP_SEQ_FOR_EACH(Z_DAT_CONF_DEFS_REPEAT, ~, Z_DAT_VARIABLES)

ZConf::ZConf(ZCommandLine const& cmdLine, ZDATFile const& datFile) :
        cmd(cmdLine),
        dat(datFile)
{
}

QString ZConf::outputPathSuffix() const
{
	QString name_addition;
	if (!annotate_name()) {
		return name_addition;
	}
	switch (removal_rule().index()) {
	case RemovalRule::CoreRemoval:
		name_addition += "CAZ_";
		break;
	case RemovalRule::AdditiveRemoval:
		name_addition += "ABF_";
		break;
	case RemovalRule::TargetRemoval:
		name_addition += "TBF_";
		break;
	case RemovalRule::GeneralizedRemoval:
		name_addition += "GBF_";
		break;
	default: // UnknownRemoval
		name_addition += "RND_";
		break;
	}
	if (use_mask()) {
		name_addition += "M";
	}
	if (use_cost()) {
		name_addition += "C";
	}
	if (use_condition()) {
		name_addition += "D"; // xxx2.1
	}
	if (use_retention()) {
		name_addition += "R"; // xxx3.0
	}
	if (use_edge_rem()) {
		name_addition += "E";
	}
	if (add_edge_cnt() > 0) {
		name_addition += "A";
	}
	if (cmd.distributionSmoothingOn) {
		name_addition += "S";
		name_addition += QString::number(static_cast<int>(100.0 * cmd.dispersalKernelMultiplier));
	}
	if (cmd.uncertaintyAlpha != 0.0) {
		name_addition += "IG";
		name_addition += QString::number(static_cast<int>(100.0 * cmd.uncertaintyAlpha));
	}
	if (use_BQP()) {
		name_addition += "BQP";
	}
	if (BLP() > 0.0) {
		name_addition += "BLP";
		name_addition += QString::number(static_cast<int>(1000.0 * BLP()));
	}
	if(!name_addition.isEmpty()) {
		name_addition.prepend('.');
	}
	return name_addition;
}

FilePath ZConf::emfOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".emf";
}

FilePath ZConf::rankOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".rank.asc";
}

FilePath ZConf::admu_redist_rankOutputPath() const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + ".ADMU.redistributed.rank.asc";
}

FilePath ZConf::propOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".prop.asc";
}

FilePath ZConf::wrscrOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".wrscr.asc";
}

FilePath ZConf::featuresInfoOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".features_info.txt";
}

FilePath ZConf::curvesOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".curves.txt";
}

FilePath ZConf::ssiFeaturesInfoOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".SSI_features_info.txt";
}

FilePath ZConf::ssiCurvesOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".SSI_curves.txt";
}

FilePath ZConf::groupCurvesOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".grp_curves.txt";
}

const QString ADMUWeights_append = ".ADMU_weights.txt";

FilePath ZConf::ADMUWeightsOutputPath() const
{
	return cmd.outFile.removeExtension() + outputPathSuffix() + ".rank" + ADMUWeights_append;
}

FilePath ZConf::ADMUCurvesOutputPath(int i) const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + QString(".rank.per_ADMU_outputs/ADMU.%1.curves.txt").arg(i);
}

FilePath ZConf::ADMUGroupCurvesOutputPath(int i) const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + QString(".rank.per_ADMU_outputs/ADMU.%1.grp_curves.txt").arg(i);
}

FilePath ZConf::runinfoOutputPath() const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + ".run_info.txt";
}

FilePath ZConf::txtOutputPath() const
{
  return cmd.outFile;
}

FilePath ZConf::rankOutputPath(ZGridFormat format) const
{
	return rankOutputPath().changeExtension(format.value().extension());
}

FilePath ZConf::admu_redist_rankOutputPath(ZGridFormat format) const
{
  return admu_redist_rankOutputPath().changeExtension(format.value().extension());
}

FilePath ZConf::propOutputPath(ZGridFormat format) const
{
	return propOutputPath().changeExtension(format.value().extension());
}

FilePath ZConf::wrscrOutputPath(ZGridFormat format) const
{
	return wrscrOutputPath().changeExtension(format.value().extension());
}

const QString output_ig_layers_dir_suffix = "_transf_info-gap_layers";
FilePath ZConf::transf_ig_layers_DirPath() const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + output_ig_layers_dir_suffix;
}

const QString output_ds_layers_dir_suffix = "_transf_distrib_smooth_layers";
const QString output_ds_layer_suffix = "_dist_smooth";
FilePath ZConf::transf_ds_layers_DirPath() const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + output_ds_layers_dir_suffix;
}

FilePath ZConf::transf_ds_layer_DirPath(QString spp_layer, ZGridFormat format) const
{
  if (spp_layer.isNull() || spp_layer.isEmpty())
    return QString();
  // take file name (without path) and remove extension
  QString layer_name = changeFileExtension(spp_layer);
  layer_name = layer_name.mid(layer_name.lastIndexOf('/')+1);

  QString res = cmd.outFile.removeExtension() + outputPathSuffix() + output_ds_layers_dir_suffix + "/" +
    changeFileExtension( layer_name + output_ds_layer_suffix, format.value().extension()); // ".asc" and so on
  return res;
}

const QString output_cst_layers_dir_suffix = "_transf_comm_simil_layers";
FilePath ZConf::transf_cst_layers_DirPath() const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + output_cst_layers_dir_suffix;
}

const QString output_ct_layers_dir_suffix = "_transf_condition_layers";
FilePath ZConf::transf_ct_layers_DirPath() const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + output_ct_layers_dir_suffix;
}

const QString output_rt_layers_dir_suffix = "_transf_retention_layers";
FilePath ZConf::transf_rt_layers_DirPath() const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + output_rt_layers_dir_suffix;
}

const QString output_mct_layers_dir_suffix = "_transf_matrix_conn_layers";
FilePath ZConf::transf_mct_layers_DirPath() const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + output_mct_layers_dir_suffix;
}

const QString output_ia_layers_dir_suffix = "_transf_interactions_layers";
FilePath ZConf::transf_ia_layers_DirPath() const
{
  return cmd.outFile.removeExtension() + outputPathSuffix() + output_ia_layers_dir_suffix;
}

FilePath ZConf::transformed_final_layers_DirPath() const
{
  return (cmd.outFile.removeExtension() + outputPathSuffix() + "_transf_final_layers");
}

QSet<ZGridFormat> ZConf::gridOutputFormats() const
{
  // Default changed from traditional (<=3.1): ASC
  return boost::get_optional_value_or(cmd.gridOutputFormats, QSet<ZGridFormat>() << ENUM_INSTANCE(ZGridFormat, CTIF));
}

QSet<ZImageFormat> ZConf::imageOutputFormats() const
{
  // Default changed from traditional (<=3.1): EMF
  return boost::get_optional_value_or(cmd.imageOutputFormats, QSet<ZImageFormat>() << ENUM_INSTANCE(ZImageFormat, PNG));
}

const ZDATFile& 
ZConf::getDATFile() const
{
  return dat;
}
