#include "proj.h"
#include "pod.h"
#include "proj_io_util.h"

ZInstance::ZInstance(const boost::shared_ptr<ZProject>& project, int row) :
        project_(project), row_(row)
{
}

const boost::shared_ptr<ZCommandLine>& ZInstance::commandLine()
{
	return project()->batFile_.constPod().commandList.at(row_);
}

boost::shared_ptr<ZProject> ZInstance::project()
{
	return project_.lock();
}

DirPath ZProject::dir()
{
	return DirPath(batFilePath());
}

FilePath ZProject::batFilePath()
{
	return batFile_.path();
}

ZInstance::OutputEntry::OutputEntry(ZGridFormat format, FileRef<ZRasterFormat> const& file) :
        format(format), file(file)
{
}

FilePath ZInstance::rankLayerOutputPath(int row) const
{
	return conf_->rankOutputPath(rankLayerOutputFiles_.at(row).format);
}

FilePath ZInstance::propLayerOutputPath(int row) const
{
	return conf_->propOutputPath(propLayerOutputFiles_.at(row).format);
}

FilePath ZInstance::wrscrLayerOutputPath(int row) const
{
	return conf_->wrscrOutputPath(wrscrLayerOutputFiles_.at(row).format);
}

FilePath ZInstance::admu_redist_rankLayerOutputPath(int row) const
{
	return conf_->admu_redist_rankOutputPath(admu_redist_rankLayerOutputFiles_.at(row).format);
}

void ZInstance::reloadOutputFiles()
{
	SuppressErrors suppress;
	FileManager& manager(datFile_.file->manager);
	DirPath root(project()->dir());
	LoadOptionalFile loadFile(root, manager, suppress);
	loadFile(runinfoOutputFile_, conf_->runinfoOutputPath());
	loadFile(curvesOutputFile_, conf_->curvesOutputPath());
	for(int i = 0; i < rankLayerOutputFiles_.size(); ++i) {
		loadFile(rankLayerOutputFiles_[i].file, conf_->rankOutputPath(rankLayerOutputFiles_[i].format));
		loadFile(propLayerOutputFiles_[i].file, conf_->propOutputPath(propLayerOutputFiles_[i].format));
		loadFile(wrscrLayerOutputFiles_[i].file, conf_->wrscrOutputPath(wrscrLayerOutputFiles_[i].format));
		loadFile(admu_redist_rankLayerOutputFiles_[i].file, conf_->admu_redist_rankOutputPath(admu_redist_rankLayerOutputFiles_[i].format));
	}
}
