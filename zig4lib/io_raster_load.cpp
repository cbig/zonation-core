#include "io.h"
#include "pod.h"
#include <gdal/gdal_priv.h>
#include <QFile>
#include <boost/shared_ptr.hpp>

namespace {
ErrorCallback rasterErr;
void CPL_STDCALL cplErrorHandler(CPLErr eErrClass, int err_no, const char *err) {
	rasterErr(QObject::tr(err));
}
}

// error reporting is not thread safe!
bool loadZRasterFile(ZRasterFile& rasterFile, QIODevice& device, ErrorCallback& err)
{
	using namespace boost;
	// device needs to be a file
	QFile *file;
	file = qobject_cast<QFile *>(&device);
	if(!file) {
	  ZError::err_msg(QObject::tr("trying to read raster from a non-file"), zeError);
	  return false;
	}
	// install temporary error handler
	rasterErr = err;
	struct PushErrorHandler {
		PushErrorHandler() {
			CPLPushErrorHandler(::cplErrorHandler);
		}
		~PushErrorHandler() {
			CPLPopErrorHandler();
		}
	} pushErrorHandler;
	// load raster
	QString filename(QDir::toNativeSeparators(file->fileName()));
	shared_ptr<GDALDataset> dataset(static_cast<GDALDataset *>(GDALOpen(filename.toUtf8().constData(), GA_ReadOnly)), GDALClose);
	if(!dataset) {
	  ZError::err_msg(QObject::tr("unable to load raster file %1").arg(filename), zeError);
	  return false;
	}
	if(dataset->GetRasterCount() < 1) {
	  ZError::err_msg(QObject::tr("raster dataset %1 does not contain any raster layers").arg(filename), zeError);
	  return false;
	}
	rasterFile.xsize = dataset->GetRasterXSize();
	rasterFile.ysize = dataset->GetRasterYSize();
	rasterFile.projection = QString::fromUtf8(dataset->GetProjectionRef());
	dataset->GetGeoTransform(rasterFile.transform);
	return true;
}
