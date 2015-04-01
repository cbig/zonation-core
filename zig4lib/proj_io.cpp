#include "proj_io.h"
#include "pod.h"
#include "zconf.h"
#include "proj_io_util.h"
#include <QTime>

// for ZError class, defined in error_util.h
QString ZError::level_prefix[] = { "", "", "", "", "", "", "", "", "", "", 
				   "Notice: ", "Warning: ", "Error: ", "* Critial error: " };
QString ZError::currentProject;
QString ZError::currentInstance;
QString ZError::iline;
int ZError::proj_count = 0;
int ZError::inst_count = 0;
// qDebug()
//ErrorCallback ZError::err = DefaultErrorCallback();
zeLevel ZError::reporting_lvl = zeError;
ErrorCallback* ZError::err = new DefaultErrorCallback(); //new SuppressErrors();

ErrorCallback* ErrorCallback::config_callback = NULL;

void 
ErrorCallback::operator()(const QString& str, zeLevel lvl)  
{ 
  if (!config_callback) { 
    // nothing
    //std::cerr << str.toStdString();  
    //qDebug() << str;  // same as default
  } else {    
    (*config_callback)(str, lvl);
  }
}


void
DefaultErrorCallback::operator()(const QString& str, zeLevel lvl)
{
  //std::cerr << str.toStdString() << std::endl;
  //qDebug() << str;
}

void
ZError::setErrorCallback(ErrorCallback& ecb)
{ 
  err = &ecb;
  ErrorCallback::config_callback = &ecb; 
}

void
ZError::setErrorReportingLevel(zeLevel lvl)
{ 
  reporting_lvl = lvl;
}
  
void
ZError::startProject(QString pname)
{ 
  currentProject = pname; 
  currentInstance = "";
  proj_count = 0; 
  inst_count = 0;
}
  
void
ZError::startInstance(QString iname, QString line)
{ 
  currentInstance = iname; 
  inst_count = 0; 
  iline = line; 
}

void
ZError::err_msg(QString msg, zeLevel lvl)
{ 

  // silently ignore undesired error messages
  if (!ZError::err || lvl < ZError::reporting_lvl)
    return;

  if (0 == proj_count && !currentProject.trimmed().isEmpty())
    (*err)(QString("---------- In project '%1', at %2 ----------").
	   arg(currentProject).arg(QTime::currentTime().toString("HH:mm:ss")), zeHeader);

  if (0 == inst_count && !currentInstance.trimmed().isEmpty())
    (*err)(QString("----- In instance %1, line: '%2' -----").
	   arg(currentInstance, iline), zeHeader);
  // qDebug() << "ZError::err_msg(): typeid: " << typeid(*err).name() << std::endl;
  (*err)(msg, lvl);
  proj_count++;
  inst_count++;
}

namespace {
class GetAnalysisRaster : public boost::static_visitor<FileRef<ZRasterFormat> >
{
private:
	DirPath root;
	FileManager& manager;
	ErrorCallback err;
	FileRef<ZRasterFormat> getRaster(const FilePath& path)
	{
	  return manager.getLinked<ZRasterFormat>(root.absoluteFilePath(path), err);
	}
public:
	GetAnalysisRaster(const DirPath& root, FileManager& manager, ErrorCallback& err) :
	        root(root), manager(manager), err(err) {}
	FileRef<ZRasterFormat> operator()(const ZLSI& lsi)
	{
		return FileRef<ZRasterFormat>();
	}
	FileRef<ZRasterFormat> operator()(const ZLSC& lsc)
	{
		return getRaster(lsc.comp_solution);
	}
	FileRef<ZRasterFormat> operator()(const ZLSM& lsm)
	{
		return getRaster(lsm.mask_file);
	}
	FileRef<ZRasterFormat> operator()(const ZLSB& lsb)
	{
		return getRaster(lsb.mask_file);
	}
};
}

// I'm not responsible for this...
// this horror is needed because of the use of "BOOST_FUSION_ADAPT_STRUCT" (ppa_adapt.h)
namespace {
class GetPPAIsLSI : public boost::static_visitor<bool>
{
public:
  GetPPAIsLSI()
  {}
  bool operator()(const ZLSI& lsi) {
    return true;
  }
  bool operator()(const ZLSC& lsc) {
    return false;
  }
  bool operator()(const ZLSM& lsm) {
    return false;
  }
  bool operator()(const ZLSB& lsb) {
    return false;
  }  
};
}

/** 
 The checks done here could be done in the two operator() of the
 LoadOptionalFile template but I just kept them as they were and added
 this as a wrapper on top

 file_opt_name is the name of the 'file name' option corresponding the option given in opt_name
 example: 


 * Yes, when using this template function below (in loadZInstance), all
the option names are duplicated from pod.h. The proper thing would be
to get rid of the horrible boost defines and use an array/set of options.

*/ 
template <typename Format>
bool loadOptionalFile_w_error_check(LoadOptionalFile& loadOptionalFile,
				       const boost::optional<bool> opt, 
				       FileRef<Format>& file, 
				       const boost::optional<FilePath>& path,
				       QString opt_name,
				       QString file_opt_name
				       )
{
  bool result = false;
  if ( opt && *opt) {
    if ( !path ) {
      ZError::err_msg(opt_name + " is ON, but no file name given in " + file_opt_name, zeError);
    } else {
      result = loadOptionalFile(file, path);
      if (!result) 
	ZError::err_msg(opt_name + " is ON, but the file specified in " + file_opt_name + " cannot be loaded: '" + *path, zeError);
    }
  } else {
    result = loadOptionalFile(file, path);
    if (!result) {
      if (!opt_name.isEmpty())
	ZError::err_msg(" even though " + opt_name  + " is OFF, a file is specified in " + file_opt_name + ", and it cannot be loaded: '" + *path + "'. Does this look like a mistake?", zeNotice);
      else
	ZError::err_msg(" File is specified in " + file_opt_name + ", cannot be loaded: '" + *path + "'. Does this look like a mistake?", zeNotice);
	}
  }
  return result;
}

// ad-hoc parser for (extended, up to 7 columns) groups file
bool 
parse_groups_file(DirPath& root, const boost::optional<FilePath>& path, ErrorCallback& err, int& num_rows)
{
  if (!path)
    return false;

  FILE  *f;
  const size_t MAX_LINE_LEN = 2048;
  char  line[MAX_LINE_LEN];
  int   num, cond_num, ret_num, ret_mode, arb_kernel_ds, arb_kernel_matrix, arb_kernel_ia;

  f=fopen(root.absoluteFilePath(*path).toStdString().c_str(), "r+t");
  if (!f)
    return false;

  bool ok = true;
  
  int row = 0;
  int text_row = 0;
  while(fgets(line, MAX_LINE_LEN, f)) {
    if (line[0]=='#') {
      text_row++;
      continue;
    }
    if (strlen(line)<2) {
      text_row++;
      continue;
    }

    text_row++;
    int cnt=sscanf(line, "%i %i %i %i %i %i %i", &num, &cond_num, &ret_num, &ret_mode, 
	       &arb_kernel_ds, &arb_kernel_matrix, &arb_kernel_ia); //  xxx add retention later
    if (cnt >= 5 && cnt <=7) {
      row++;
    } else {
      ZError::err_msg(QString(" in groups file, '%1', row %2, wrong number of columns, expecting between 5 and 7").arg(*path).arg(text_row), zeError);
      ok = false;
    }
    
  }
  fclose(f);

  num_rows = row;
  return ok;
}

// Like loadOptionalFile_w_error_check but specific to the groups file
// Example call: load_optional_file_groups(dat.use_groups, temp->groupsFile_, dat.groups_fname, "'use groups'", "'groups file'");
bool
load_optional_file_groups(DirPath& root, const boost::optional<bool> opt, 
			  const boost::optional<FilePath>& path,
			  QString opt_name, QString file_opt_name,
			  ErrorCallback& err, int& num_rows)
{
  bool result = false;
  if ( opt && *opt) {
    if ( !path ) {
      ZError::err_msg(opt_name + " is ON, but no file name given in " + file_opt_name, zeError);
    } else {
      result = parse_groups_file(root, path, err, num_rows);
      if (!result) 
	ZError::err_msg(opt_name + " is ON, but the file specified in " + file_opt_name + " cannot be loaded: '" + *path, zeError);
    }
  } else {
    if (!path)
      return false;
    result = parse_groups_file(root, path, err, num_rows);
    if (!result) {
      if (!opt_name.isEmpty())
	ZError::err_msg(" even though " + opt_name  + " is OFF, a file is specified in " + file_opt_name + ", and it cannot be loaded: '" + *path + "'. Does this look like a mistake?", zeNotice);
      else
	ZError::err_msg(" File is specified in " + file_opt_name + ", cannot be opened: '" + *path + "'. Does this look like a mistake?", zeNotice);
	}
  }

  return result;
}

/**
 This is for the following (not so common) cases:

 if opt_name is an empty string, there is no 'use...' options associated to
 file_opt_name (such as "connectivity edge effect fix file"

*/
template <typename Format>
bool loadOptionalFile_w_error_check(LoadOptionalFile& loadOptionalFile,
				       FileRef<Format>& file, 
				       const boost::optional<FilePath>& path,
				       QString file_opt_name
				       )
{
  bool result = false;
  result = loadOptionalFile(file, path);
  if (!result) 
    ZError::err_msg(" The file specified in " + file_opt_name + " cannot be opened: '" + *path, zeError);
  return result;
}

boost::shared_ptr<ZInstance> loadZInstance(const boost::shared_ptr<ZProject>& project, int row, const ZCommandLine& command, FileManager& manager, ErrorCallback& err)
{
	using namespace boost;
	shared_ptr<ZInstance> instance;
	shared_ptr<ZInstance> temp(make_shared<ZInstance>(project, row));

	DirPath root(project->dir());

	ZError::startInstance(QString("#%1").arg(row+1), command.line.trimmed());

	//LineModify batPrepend(err, QObject::tr("%1: ").arg(project->batFilePath()));

	bool result = false;

	// load precalculated rank
	if(command.commandLineMode.index() == ZCommandLineMode::LoadRank) {
	  temp->rankFile_ = manager.getLinked<ZRasterFormat>(root.absoluteFilePath(command.rankFile), err);
	  if (!temp->rankFile_)
	    ZError::err_msg(" Unable to load rank file: '" + root.absoluteFilePath(command.rankFile) + "'. It does not exist or is not readable.", zeError);
	}
	// load dat
	temp->datFile_ = manager.getLinked<ZDATFormat>(root.absoluteFilePath(command.datFile), err);
	if(!temp->datFile_) {
	  ZError::err_msg(" Unable to load settings file: '" + root.absoluteFilePath(command.datFile) + "'. It does not exist or is not readable.", zeError);
	  return instance;
	}
	// load spp
	temp->sppFile_ = manager.getLinked<ZSPPFormat>(root.absoluteFilePath(command.sppFile), err);
	if(!temp->sppFile_) {
	  ZError::err_msg(" Unable to load species/features list file: '" + root.absoluteFilePath(command.sppFile) + "'. It does not exist or is not readable.", zeError);
	  return instance;
	}

	temp->conf_ = make_shared<ZConf>(command, *temp->datFile_);

	// TODO: err(QObject::tr("This is a typical in-sequence error report, loading spp file: %1").arg(temp->sppFile_.path()));

	// load spp rasters
	//LineModify sppPrepend(batPrepend, QObject::tr("%1:").arg(temp->sppFile_.path()));
	foreach(const shared_ptr<ZSPPEntry>& layer, temp->sppFile_.constPod().sppList) {
	  temp->sppLayerFiles_ << manager.getLinked<ZRasterFormat>(root.absoluteFilePath(layer->raster()), err);
	  if (!temp->sppLayerFiles_.last())
	    ZError::err_msg(" Unable to load species/features file: '" + root.absoluteFilePath(layer->raster())+ "'. It does not exist or is not readable. This instance will not run correctly!", zeError);
	}

	//LineModify datPrepend(batPrepend, QObject::tr("%1: ").arg(temp->datFile_.path()));
	//LoadOptionalFile loadOptionalFile(root, manager, datPrepend);
	LoadOptionalFile loadOptionalFile(root, manager, err);

	const ZDATFile& dat(temp->datFile_.constPod());

	// bqp
	//loadOptionalFile(temp->bqpFile_, dat.BQP_prof_fn);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_BQP, 
				       temp->bqpFile_, dat.BQP_prof_fn,
				       "'use bqp'", "'bqp profiles file'");


	// ssi
	if(dat.SSI_fname && !dat.SSI_fname->isEmpty()) {
		if(temp->ssiFile_ = manager.getLinked<ZSSIFormat>(root.absoluteFilePath(*dat.SSI_fname), err)) {
		  //LineModify ssiPrepend(datPrepend, QObject::tr("%1: ").arg(temp->ssiFile_.path()));
			foreach(const shared_ptr<ZSSIEntry>& layer, temp->ssiFile_.constPod().ssiList) {
				temp->ssiLayerFiles_ << manager.getLinked<ZSSICoordinatesFormat>(root.absoluteFilePath(layer->ssi_), err);
				if (!temp->ssiLayerFiles_.last())
				  ZError::err_msg(" Unable to load SSI layer file: '" + root.absoluteFilePath(layer->ssi_)+ "'. It does not exist or is not readable.", zeWarning);
			}
		} else {
		  ZError::err_msg(" Unable to load SSI file: '" + root.absoluteFilePath(*dat.SSI_fname)+ "'. It does not exist or is not readable.", zeWarning);
		}
	}

	// condition layers list file
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_condition, 
				       temp->condFile_, dat.condition_fname,
				       "'use condition layer'", "'condition file'");

	// condition
	if(dat.condition_fname && !dat.condition_fname->isEmpty()) {
		if(temp->condFile_ = manager.getLinked<ZIndexedRastersFormat>(root.absoluteFilePath(*dat.condition_fname), err)) {
		  //LineModify condPrepend(datPrepend, QObject::tr("%1: ").arg(temp->condFile_.path()));
			foreach(const shared_ptr<ZIndexedRastersEntry>& layer, temp->condFile_.constPod().rasterList) {
				temp->condLayerFiles_ << manager.getLinked<ZRasterFormat>(root.absoluteFilePath(layer->raster), err);
				if (!temp->condLayerFiles_.last())
				  ZError::err_msg(" Unable to load condition layer file: '" + root.absoluteFilePath(layer->raster)+ "'. It does not exist or is not readable.", zeWarning);
			}
		} else {
		  ZError::err_msg(" Unable to load condition file: '" + root.absoluteFilePath(*dat.condition_fname)+ "'. It does not exist or is not readable.", zeWarning);
		}
	}

	// retention layers list file
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_retention, 
				       temp->retFile_, dat.retention_fname,
				       "'use retention layer'", "'retention file'");

	// retention
	if(dat.retention_fname && !dat.retention_fname->isEmpty()) {
		if(temp->retFile_ = manager.getLinked<ZIndexedRastersFormat>(root.absoluteFilePath(*dat.retention_fname), err)) {
		  //LineModify retPrepend(datPrepend, QObject::tr("%1: ").arg(temp->retFile_.path()));
			foreach(const shared_ptr<ZIndexedRastersEntry>& layer, temp->retFile_.constPod().rasterList) {
				temp->retLayerFiles_ << manager.getLinked<ZRasterFormat>(root.absoluteFilePath(layer->raster), err);
				if (!temp->retLayerFiles_.last())
				  ZError::err_msg(" Unable to load retention layer file: '" + root.absoluteFilePath(layer->raster) + "'. It does not exist or is not readable.", zeWarning);
			}
		} else {
		  ZError::err_msg(" Unable to load retention file: '" + root.absoluteFilePath(*dat.retention_fname)+ "'. It does not exist or is not readable.", zeWarning);
		}
	}

	// local edge correction
	if(dat.LEC_fname && !dat.LEC_fname->isEmpty()) {
		if(temp->lecFile_ = manager.getLinked<ZIndexedRastersFormat>(root.absoluteFilePath(*dat.LEC_fname), err)) {
		  //LineModify lecPrepend(err, QObject::tr("%1: ").arg(temp->lecFile_.path()));
			foreach(const shared_ptr<ZIndexedRastersEntry>& layer, temp->lecFile_.constPod().rasterList) {
				temp->lecLayerFiles_ << manager.getLinked<ZRasterFormat>(root.absoluteFilePath(layer->raster), err);
				if (!temp->condLayerFiles_.last())
				  ZError::err_msg(" Unable to load local edge correction (LEC) layer file: '" + root.absoluteFilePath(layer->raster)+ "'. It does not exist or is not readable.", zeWarning);

			}
		} else {
		  ZError::err_msg(" Unable to load local edge correction (LEC) file: '" + root.absoluteFilePath(*dat.LEC_fname)+ "'. It does not exist or is not readable.", zeWarning);
		}
	}

	// info-gap
	if(dat.use_IGw && dat.IGwfn && !dat.IGwfn->isEmpty()) {
		if(temp->igFile_ = manager.getLinked<ZIGFormat>(root.absoluteFilePath(*dat.IGwfn), err)) {
		  //LineModify igPrepend(datPrepend, QObject::tr("%1: ").arg(temp->igFile_.path()));
			foreach(const shared_ptr<ZIGEntry>& layer, temp->igFile_->igList) {
				temp->igLayerFiles_ << manager.getLinked<ZRasterFormat>(root.absoluteFilePath(layer->raster), err);
				if (!temp->igLayerFiles_.last())
				  ZError::err_msg(" Unable to load info-gap layer file: '" + root.absoluteFilePath(layer->raster)+ "'. It does not exist or is not readable.", zeWarning);				
			}
		} else {
		  if (*dat.use_IGw)
		    ZError::err_msg(" 'use info-gap weights ' is ON, but it was not possible to load info-gap weights file: '" + root.absoluteFilePath(*dat.IGwfn) + "'. It does not exist or is not readable.", zeWarning);
		  else {
		    if (!dat.IGwfn->isEmpty())
		      ZError::err_msg(" Even though 'use info-gap weights ' is OFF, 'Info-gap weights file' is set to: '" + root.absoluteFilePath(*dat.IGwfn) + "'. Does this look like a mistake?", zeNotice);
		  }
		}
	}

	// post processing
	if(dat.PPA_fname && !dat.PPA_fname->isEmpty()) {
	  if(temp->ppaFile_ = manager.getLinked<ZPPAFormat>(root.absoluteFilePath(*dat.PPA_fname), err)) {
	    //LineModify ppaPrepend(datPrepend, QObject::tr("%1: ").arg(temp->ppaFile_.path()));
	    GetAnalysisRaster getRaster(root, manager, err);
	    // this horror is needed because of the use of "BOOST_FUSION_ADAPT_STRUCT" (ppa_adapt.h)
	    GetPPAIsLSI get_is_lsi;
	    int line_idx = 1;
	    foreach(const shared_ptr<ZAnalysis>& analysis, temp->ppaFile_->analysisList) {
	      temp->ppaLayerFiles_ << boost::apply_visitor(getRaster, *analysis);
	      if (!boost::apply_visitor(get_is_lsi, *analysis) && !temp->ppaLayerFiles_.last())
		ZError::err_msg(QString("Unable to load post processing (PPA) layer in analysis/line number %1. It does not exist or is not readable.").arg(line_idx), zeWarning);
	      line_idx++;
	    }
	  } else {
	    ZError::err_msg(" Unable to load post processing (PPA) file: '" + root.absoluteFilePath(*dat.PPA_fname)+ "'. It does not exist or is not readable.", zeNotice);
	  }
	}

	// groups
	//loadOptionalFile(temp->groupsFile_, dat.groups_fname);
	/*
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_groups, 
				       temp->groupsFile_, dat.groups_fname,
				       "'use groups'", "'groups file'");
	*/
	int num_group_rows = 0;
	if (load_optional_file_groups(root, dat.use_groups, dat.groups_fname, "'use groups'", "'groups file'", err, num_group_rows)) {
	  // Still, this needs to be done to link the FileRef object so the tree shows an "ok"
	  //loadOptionalFile(temp->groupsFile_, dat.groups_fname);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_groups, 
				       temp->groupsFile_, dat.groups_fname,
				       "'use groups'", "'groups file'");

	  if (num_group_rows != temp->sppFile_.constPod().sppList.size())
	    ZError::err_msg(QString("The number of rows in the groups file (%1) does not match the number of features listed in the features list file (%2). Please, make sure that both files are consistent: %3, %4").arg(num_group_rows).arg(temp->sppFile_.constPod().sppList.size()).arg(*dat.groups_fname).arg(command.sppFile), zeError );
	}

	// cost layer
	//loadOptionalFile(temp->costLayerFile_, dat.costfn);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_cost, 
				       temp->costLayerFile_, dat.costfn,
				       "'use cost'", "'cost file'");

	// mask
	//loadOptionalFile(temp->maskLayerFile_, dat.maskfn);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_mask, 
				       temp->maskLayerFile_, dat.maskfn,
				       "'use mask'", "'mask file'");

	// area mask
	//loadOptionalFile(temp->areaMaskLayerFile_, dat.area_mask_fname);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.mask_data, 
				       temp->areaMaskLayerFile_, dat.area_mask_fname,
				       "'mask missing areas'", "'area mask file'");

	// plula
	//loadOptionalFile(temp->plulaLayerFile_, dat.PLULAfn);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_PLULA, 
				       temp->plulaLayerFile_, dat.PLULAfn,
				       "'use planning unit layer'", "'planning unit layer file'");

	// edge effect fix
	//loadOptionalFile(temp->edgeEffectLayerFile_, dat.edge_effect_fn);
	loadOptionalFile_w_error_check(loadOptionalFile, 
				       temp->edgeEffectLayerFile_, dat.edge_effect_fn,
				       "'connectivity edge effect fix file'");

	// admu layer
	//loadOptionalFile(temp->admuLayerFile_, dat.ADM_layer_file);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_ADMUs, 
				       temp->admuLayerFile_, dat.ADM_layer_file,
				       "'use ADMUs'", "'ADMU layer file'");

	// admu
	//loadOptionalFile(temp->admuFile_, dat.ADM_weights_file);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_ADMUs, 
				       temp->admuFile_, dat.ADM_weights_file,
				       "'use ADMUs'", "'ADMU descriptions file'");


	// admu matrix
	//loadOptionalFile(temp->admuMatrixFile_, dat.ADM_weight_matrix_file);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_ADMUs, 
				       temp->admuMatrixFile_, dat.ADM_weight_matrix_file,
				       "'use ADMUs'", "'ADMU weight matrix'");

	
	/*
	[Community analysis settings]
	  load similarity matrix = 1
	  connectivity similarity matrix file = matrix_file_name.txt
	  apply to connectivity = 1
	*/
	// connectivity similarirty matrix
	//loadOptionalFile(temp->connMatrixFile_, dat.sim_matrix_name);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.apply_to_conn, 
				       temp->connMatrixFile_, dat.sim_matrix_name,
				       "'apply to connectivity'", "'connectivity similarity matrix file'");


	/*
	  Community analysis settings]
	  load similarity matrix = 1
	  community similarity matrix file = comm_sim_matrix.txt
	  apply to representation = 1 
	*/
	// community similarity matrix
	//loadOptionalFile(temp->commMatrixFile_, dat.comm_sim_matrix_name);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.apply_to_repr, 
				       temp->commMatrixFile_, dat.comm_sim_matrix_name,
				       "'apply to representation'", "'community similarity matrix file'");

	// tree
	//loadOptionalFile(temp->treeFile_, dat.tree_fname);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_tree_conn, 
				       temp->treeFile_, dat.tree_fname,
				       "'use tree connectivity'", "'tree connectivity file'");
	

	// interactions
	loadOptionalFile(temp->interactionsFile_, dat.ia_fname);
	loadOptionalFile_w_error_check(loadOptionalFile, dat.use_interactions, 
				       temp->interactionsFile_, dat.ia_fname,
				       "'use interactions'", "'interaction file'");

	// Corridor domain layers
	if (dat.use_corridors && *dat.use_corridors){
	  if ((dat.use_corridor_domain_layers && *dat.use_corridor_domain_layers) && 
	      (!dat.corridor_domain_layers_fn || (dat.corridor_domain_layers_fn && dat.corridor_domain_layers_fn->isEmpty()))) {
	    ZError::err_msg(" 'use domain layers' is ON, but 'domain layers list file' is empty'. Does this look like a mistake?", zeNotice);	  
	  } else if ((!dat.use_corridor_domain_layers || !(*dat.use_corridor_domain_layers)) &&
		     (dat.corridor_domain_layers_fn && !dat.corridor_domain_layers_fn->isEmpty())) {
	    QString val = "";
	    if (dat.corridor_domain_layers_fn)
	      val = root.absoluteFilePath(*dat.corridor_domain_layers_fn);
	    else
	      val = "";
	    ZError::err_msg(" Even though 'use domain layers' is OFF, 'domain layers list file' is set to: '" + val + "'. Does this look like a mistake?", zeNotice);
	  } else if((dat.use_corridor_domain_layers && *dat.use_corridor_domain_layers) && 
		    (dat.corridor_domain_layers_fn && !dat.corridor_domain_layers_fn->isEmpty())) {
	    if(temp->corridor_domains_list_file_ = manager.getLinked<ZIGFormat>(root.absoluteFilePath(*dat.corridor_domain_layers_fn), err)) {
	      //LineModify igPrepend(datPrepend, QObject::tr("%1: ").arg(temp->corridor_domains_list_file_.path()));
	      foreach(const shared_ptr<ZIGEntry>& layer, temp->corridor_domains_list_file_->igList) {
		temp->corridor_domain_layers_ << manager.getLinked<ZRasterFormat>(root.absoluteFilePath(layer->raster), err);
		if (!temp->corridor_domain_layers_.last())
		  ZError::err_msg(" Unable to load corridor domain layer: '" + root.absoluteFilePath(layer->raster)+ "'. It does not exist or is not readable.", zeWarning);				
	      }
	    } else {
	      if (*dat.use_corridor_domain_layers) {
		ZError::err_msg(" 'use domain layers' is ON, but it was not possible to load the corridor domain layers list file: '" + root.absoluteFilePath(*dat.corridor_domain_layers_fn) + "'. It does not exist or it is not correct.", zeWarning);
	      } else {
		if (dat.corridor_domain_layers_fn || !dat.corridor_domain_layers_fn->isEmpty())
		  ZError::err_msg(" Even though 'use domain layers' is OFF, 'domain layers list file' is set to: '" + root.absoluteFilePath(*dat.corridor_domain_layers_fn) + "'. Does this look like a mistake?", zeNotice);
	      }
	    }
	  }
	}

	// output files

	// suppress errors
	//SuppressErrors suppress;
	//LoadOptionalFile loadOutputFile(root, manager, suppress);
	LoadOptionalFile loadOutputFile(root, manager, err);

	loadOutputFile(temp->runinfoOutputFile_, temp->conf_->runinfoOutputPath());
	loadOutputFile(temp->curvesOutputFile_, temp->conf_->curvesOutputPath());
	foreach(ZGridFormat const& format, temp->conf_->gridOutputFormats()) {
		// can do anything to this crazy path manipulation?
		temp->rankLayerOutputFiles_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath(temp->conf_->rankOutputPath(format)), err));
		temp->propLayerOutputFiles_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath(temp->conf_->propOutputPath(format)), err));
		temp->wrscrLayerOutputFiles_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath(temp->conf_->wrscrOutputPath(format)), err));

		temp->admu_redist_rankLayerOutputFiles_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath(temp->conf_->admu_redist_rankOutputPath(format)), err));

		/*
		  // This block guessed the names of the DS transf. files from the species layers.
		  // Let's rather find whatever is available in the DS_transf subdir, see below...
		// Intermediate outputs, Distribution Smoothing layers
		int i = 0;
		foreach(const shared_ptr<ZSPPEntry>& layer, temp->sppFile_.constPod().sppList) {
		  FilePath spp_file = root.absoluteFilePath(layer->raster());
		  temp->interm_output_ds_files_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath(temp->conf_->transf_ds_layer_OutputPath(spp_file, format)), err));
		}
		*/
	}

	// Info-gap transform (IG) outputs:
	FilePath ig_subdir = root.absoluteFilePath(temp->conf_->transf_ig_layers_DirPath());
	QDir ig_d(ig_subdir);
	if (ig_d.exists()) {
	  ig_d.setFilter(QDir::Files);   // avoid '.', '..' and friends or crash!
	  ig_d.setSorting(QDir::Name);
	  QFileInfoList dl = ig_d.entryInfoList();
	  ZGridFormat const& format = temp->conf_->gridOutputFormats().values().first();
	  for (size_t i = 0; i < dl.size(); ++i) {
	    QFileInfo file_info = dl.at(i);
	    temp->transf_ig_layers_files_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath( ig_subdir + "/" + file_info.fileName() ), err));
	  }
        }

	// DS outputs:
	FilePath ds_subdir = root.absoluteFilePath(temp->conf_->transf_ds_layers_DirPath());
	QDir ds_d(ds_subdir);
	// load all the files in this subdir.
	if (ds_d.exists()) {
	  ds_d.setFilter(QDir::Files);   // avoid '.', '..' and friends or crash!
	  ds_d.setSorting(QDir::Name);
	  QFileInfoList dl = ds_d.entryInfoList();
	  // For every file, would need to probe all the formats, but don't
	  // do it cause files will be read correctly even if the format 
	  // specified is wrong (as long as they come in a supported format)
	  ZGridFormat const& format = temp->conf_->gridOutputFormats().values().first();
	  for (size_t i = 0; i < dl.size(); ++i) {
	    QFileInfo file_info = dl.at(i);
	    //foreach(ZGridFormat const& format, temp->conf_->gridOutputFormats()) {
	    temp->transf_ds_layers_files_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath( ds_subdir + "/" + file_info.fileName() ), err));
	    //}
	  }
        }

	// Community similarity transform (CST) outputs:
	FilePath cst_subdir = root.absoluteFilePath(temp->conf_->transf_cst_layers_DirPath());
	QDir cst_d(cst_subdir);
	if (cst_d.exists()) {
	  cst_d.setFilter(QDir::Files);   // avoid '.', '..' and friends or crash!
	  cst_d.setSorting(QDir::Name);
	  QFileInfoList dl = cst_d.entryInfoList();
	  ZGridFormat const& format = temp->conf_->gridOutputFormats().values().first();
	  for (size_t i = 0; i < dl.size(); ++i) {
	    QFileInfo file_info = dl.at(i);
	    temp->transf_cst_layers_files_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath( cst_subdir + "/" + file_info.fileName() ), err));
	  }
        }

	// Condition transform (CT) outputs:
	FilePath ct_subdir = root.absoluteFilePath(temp->conf_->transf_ct_layers_DirPath());
	QDir ct_d(ct_subdir);
	if (ct_d.exists()) {
	  ct_d.setFilter(QDir::Files);   // avoid '.', '..' and friends or crash!
	  ct_d.setSorting(QDir::Name);
	  QFileInfoList dl = ct_d.entryInfoList();
	  ZGridFormat const& format = temp->conf_->gridOutputFormats().values().first();
	  for (size_t i = 0; i < dl.size(); ++i) {
	    QFileInfo file_info = dl.at(i);
	    temp->transf_ct_layers_files_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath( ct_subdir + "/" + file_info.fileName() ), err));
	  }
        }

	// Retention transform (RT) outputs:
	FilePath rt_subdir = root.absoluteFilePath(temp->conf_->transf_rt_layers_DirPath());
	QDir rt_d(rt_subdir);
	if (rt_d.exists()) {
	  rt_d.setFilter(QDir::Files);   // avoid '.', '..' and friends or crash!
	  rt_d.setSorting(QDir::Name);
	  QFileInfoList dl = rt_d.entryInfoList();
	  ZGridFormat const& format = temp->conf_->gridOutputFormats().values().first();
	  for (size_t i = 0; i < dl.size(); ++i) {
	    QFileInfo file_info = dl.at(i);
	    temp->transf_rt_layers_files_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath( rt_subdir + "/" + file_info.fileName() ), err));
	  }
        }

	// Matrix connectivity (MCT) outputs:
	FilePath mct_subdir = root.absoluteFilePath(temp->conf_->transf_mct_layers_DirPath());
	QDir mct_d(mct_subdir);
	if (mct_d.exists()) {
	  mct_d.setFilter(QDir::Files);   // avoid '.', '..' and friends or crash!
	  mct_d.setSorting(QDir::Name);
	  QFileInfoList dl = mct_d.entryInfoList();
	  ZGridFormat const& format = temp->conf_->gridOutputFormats().values().first();
	  for (size_t i = 0; i < dl.size(); ++i) {
	    QFileInfo file_info = dl.at(i);
	    temp->transf_mct_layers_files_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath( mct_subdir + "/" + file_info.fileName() ), err));
	  }
        }

	// Interactions (IA) outputs:
	FilePath ia_subdir = root.absoluteFilePath(temp->conf_->transf_ia_layers_DirPath());
	QDir ia_d(ia_subdir);
	if (ia_d.exists()) {
	  ia_d.setFilter(QDir::Files);   // avoid '.', '..' and friends or crash!
	  ia_d.setSorting(QDir::Name);
	  QFileInfoList dl = ia_d.entryInfoList();
	  ZGridFormat const& format = temp->conf_->gridOutputFormats().values().first();
	  for (size_t i = 0; i < dl.size(); ++i) {
	    QFileInfo file_info = dl.at(i);
	    temp->transf_ia_layers_files_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath( ia_subdir + "/" + file_info.fileName() ), err));
	  }
        }
	
	// final transforms:
	FilePath trans_subdir = root.absoluteFilePath(temp->conf_->transformed_final_layers_DirPath());
	QDir d(trans_subdir);
	if (d.exists()) {
	  d.setFilter(QDir::Files);   // avoid '.', '..' and friends or crash!
	  d.setSorting(QDir::Name);
	  QFileInfoList dl = d.entryInfoList();
	  ZGridFormat const& format = temp->conf_->gridOutputFormats().values().first();
	  for (size_t i = 0; i < dl.size(); ++i) {
	    QFileInfo file_info = dl.at(i);
	    //foreach(ZGridFormat const& format, temp->conf_->gridOutputFormats()) {
	    temp->transformed_final_layers_files_ << ZInstance::OutputEntry(format, manager.getLinked<ZRasterFormat>(root.absoluteFilePath( trans_subdir + "/" + file_info.fileName() ), err));
	    //}
	  }
        }

	// all's well
	instance = temp;
	return instance;
}

boost::shared_ptr<ZProject> loadZProject(const FilePath& batFilePath, FileManager& manager,  ErrorCallback& err)
{
  using namespace boost;
  shared_ptr<ZProject> project;
  // manager.setOpeningProject(project);

  qDebug() <<"loadZProject() starts...";

  // normally, err should be consoleView->errorStream(), passed from mainwindow
  ZError::setErrorCallback(err);
  ZError::startProject(QString("%1").arg(batFilePath));

  FileRef<ZBATFormat> bat(manager.getLinked<ZBATFormat>(batFilePath, err));
  if(!bat) {
    return project;
  }

  //LineModify prepend(err, QObject::tr("%1: ").arg(bat.path()));
  shared_ptr<ZProject> temp(boost::make_shared<ZProject>());
  temp->batFile_ = bat;
  DirPath root(temp->dir());

  //err("load bat rows");

  int row = 0;
  foreach(const shared_ptr<ZCommandLine>& command, bat.constPod().commandList) {
    //err("load bat row");
    shared_ptr<ZInstance> instance(loadZInstance(temp, row, *command, manager, err));
    if(!instance) {
      return project;
    }

    temp->instanceList_ << instance;
    ++row;
  }

  // everything OK
  project = temp;
  return project;
}
