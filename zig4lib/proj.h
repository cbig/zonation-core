#ifndef PROJ_H
#define PROJ_H

#include "io.h"
#include "manager.h"
#include "zig4lib_global.h"
#include "zconf.h"
#include <QMap>

struct ZSPPEntry;

// objects

class ZProject;

class ZIG4LIBSHARED_EXPORT ZInstance
{
public:
	ZInstance(const boost::shared_ptr<ZProject>& project, int row);

	// Note this rank is not the output, that one is below!
	FileRef<ZRasterFormat> rankFile_;
	FileRef<ZDATFormat> datFile_;
	FileRef<ZBQPFormat> bqpFile_;
	FileRef<ZSSIFormat> ssiFile_;
	QList<FileRef<ZSSICoordinatesFormat> > ssiLayerFiles_;
	FileRef<ZIndexedRastersFormat> condFile_;
	QList<FileRef<ZRasterFormat> > condLayerFiles_;
	FileRef<ZIndexedRastersFormat> retFile_;
	QList<FileRef<ZRasterFormat> > retLayerFiles_;
	FileRef<ZIndexedRastersFormat> lecFile_;
	QList<FileRef<ZRasterFormat> > lecLayerFiles_;

	FileRef<ZIGFormat> igFile_;
	QList<FileRef<ZRasterFormat> > igLayerFiles_;

	FileRef<ZPPAFormat> ppaFile_;
	QList<FileRef<ZRasterFormat> > ppaLayerFiles_;
	FileRef<ZGroupsFormat> groupsFile_;
	FileRef<ZRasterFormat> costLayerFile_;
	FileRef<ZRasterFormat> maskLayerFile_;
	FileRef<ZRasterFormat> areaMaskLayerFile_;
	FileRef<ZRasterFormat> plulaLayerFile_;
	FileRef<ZRasterFormat> edgeEffectLayerFile_;
	FileRef<ZRasterFormat> admuLayerFile_;
	FileRef<ZADMUFormat> admuFile_;
	FileRef<ZMatrixFormat> admuMatrixFile_;
	FileRef<ZMatrixFormat> connMatrixFile_;
	FileRef<ZMatrixFormat> commMatrixFile_;
	FileRef<ZTreeFormat> treeFile_;
	FileRef<ZInteractionsFormat> interactionsFile_;

	FileRef<ZSPPFormat> sppFile_;
	QList<FileRef<ZRasterFormat> > sppLayerFiles_;

	// corridor domain layers - list file and layers
	FileRef<ZIGFormat> corridor_domains_list_file_;
	QList<FileRef<ZRasterFormat> > corridor_domain_layers_;

	// semantic stuff

	// output files
	// these paths depend on configuration

	FileRef<ZUnknownFormat> runinfoOutputFile_;
	FileRef<ZUnknownFormat> curvesOutputFile_;

	class OutputEntry
	{
	public:
		ZGridFormat format;
		FileRef<ZRasterFormat> file;
		OutputEntry(ZGridFormat format, FileRef<ZRasterFormat> const& file);
	};

	QList<OutputEntry> rankLayerOutputFiles_;
	QList<OutputEntry> propLayerOutputFiles_;
	QList<OutputEntry> wrscrLayerOutputFiles_;
	QList<OutputEntry> admu_redist_rankLayerOutputFiles_;

	// these FilePath ... methods use conf to get the file names
	FilePath rankLayerOutputPath(int row) const;
	FilePath propLayerOutputPath(int row) const;
	FilePath wrscrLayerOutputPath(int row) const;

	FilePath admu_redist_rankLayerOutputPath(int row) const;

	QList<OutputEntry> transf_ig_layers_files_;
	QList<OutputEntry> transf_ds_layers_files_;
	QList<OutputEntry> transf_cst_layers_files_;
	QList<OutputEntry> transf_ct_layers_files_;
	QList<OutputEntry> transf_rt_layers_files_;
	QList<OutputEntry> transf_mct_layers_files_;
	QList<OutputEntry> transf_ia_layers_files_;
	QList<OutputEntry> transformed_final_layers_files_;

	void reloadOutputFiles();

	// useless stuff begins here <- most useful stuff actually!

	boost::shared_ptr<ZConf> conf_;
	boost::weak_ptr<ZProject> project_;
	int row_; // used to access zcommandline from project

	// public api

	const boost::shared_ptr<ZCommandLine>& commandLine();
	boost::shared_ptr<ZProject> project();
};

class ZIG4LIBSHARED_EXPORT ZProject
{
public:
	DirPath dir();
	FilePath batFilePath();

	FileRef<ZBATFormat> batFile_;
	QList<boost::shared_ptr<ZInstance> > instanceList_;
};

#endif // PROJ_H
