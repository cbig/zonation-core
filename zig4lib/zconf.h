#ifndef ZCONF_H
#define ZCONF_H

/**
 * class ZConf holds instance options, which are obtained from the cmd line (example, zig4 executable)
 * and the settings file (.dat, which provides various paths).
 */

#include "pod.h"

// Compare with the pretty similar macro in pod.h, right before struct ZDATFile
#define Z_DAT_CONF_DECLS_REPEAT(r, data, elem) BOOST_PP_SEQ_ELEM(0, elem) BOOST_PP_SEQ_ELEM(1, elem)() const;

class ZIG4LIBSHARED_EXPORT ZConf
{
public:
	ZConf(ZCommandLine const& cmdLine, ZDATFile const& datFile);

	BOOST_PP_SEQ_FOR_EACH(Z_DAT_CONF_DECLS_REPEAT, ~, Z_DAT_VARIABLES)

	QString outputPathSuffix() const;

	FilePath emfOutputPath() const;
	FilePath rankOutputPath() const;  // .asc
	FilePath propOutputPath() const;  // .asc
	FilePath wrscrOutputPath() const; // .asc
	FilePath admu_redist_rankOutputPath() const;  // .asc


	FilePath txtOutputPath() const;
	FilePath runinfoOutputPath() const;
	FilePath featuresInfoOutputPath() const;
	FilePath curvesOutputPath() const;
	FilePath ssiFeaturesInfoOutputPath() const;
	FilePath ssiCurvesOutputPath() const;
	FilePath groupCurvesOutputPath() const;

	FilePath ADMUWeightsOutputPath() const;
	// get path for ADMU # i
	FilePath ADMUCurvesOutputPath(int i) const;
	FilePath ADMUGroupCurvesOutputPath(int i) const;

	FilePath rankOutputPath(ZGridFormat format) const;
	FilePath propOutputPath(ZGridFormat format) const;
	FilePath wrscrOutputPath(ZGridFormat format) const;
	FilePath admu_redist_rankOutputPath(ZGridFormat format) const;

	// Info-gap transform
	FilePath transf_ig_layers_DirPath() const;
	// distribution smoothing transform
	FilePath transf_ds_layers_DirPath() const;
	FilePath transf_ds_layer_DirPath(QString spp_layer, ZGridFormat f) const;
	// Community similarity transform
	FilePath transf_cst_layers_DirPath() const;
	// Condition transform
	FilePath transf_ct_layers_DirPath() const;
	// Retention transform
	FilePath transf_rt_layers_DirPath() const;
	// matrix connectivity transform
	FilePath transf_mct_layers_DirPath() const;

	// interactions transform
	FilePath transf_ia_layers_DirPath() const;

	// final transform
	FilePath transformed_final_layers_DirPath() const;

	QSet<ZGridFormat> gridOutputFormats() const;
	QSet<ZImageFormat> imageOutputFormats() const;

	const ZDATFile& getDATFile() const;

private:
	ZCommandLine const& cmd;
	ZDATFile const& dat;
};

#endif // ZCONF_H
