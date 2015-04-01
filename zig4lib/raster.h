#ifndef RASTER_H
#define RASTER_H

#include "zig4lib_global.h"
#include <gdal/gdal_priv.h>

template <typename T>
struct GDALTypeTraits {
	static const GDALDataType GDALType = GDT_Unknown;
	static const int size = 0;
};

// not really cross-platform specializations...

template <>
struct GDALTypeTraits<float>
{
	static const GDALDataType GDALType = GDT_Float32;
	static const size_t size = 4;
};

template <>
struct GDALTypeTraits<int>
{
	static const GDALDataType GDALType = GDT_Int32;
	static const size_t size = 4;
};

template <>
struct GDALTypeTraits<unsigned char>
{
	static const GDALDataType GDALType = GDT_Byte;
	static const size_t size = 1;
};

struct GlobalProjection {
	bool initialized;
	double globalTransformation[6];
	std::string globalProjection;
	GlobalProjection() :
	        initialized(false)
	{
		for(double *i = globalTransformation; i != globalTransformation + 6; ++i) {
			*i = 0.0;
		}
	}
};

class RasterBase
{
public:
	ZIG4LIBSHARED_EXPORT static GlobalProjection globalProjection;
};

template <typename T>
class Raster : public RasterBase
{
private:
	GDALDataset *dataset;
	GDALRasterBand *band;
	bool rasterLoaded_;
	bool statsNotLoaded;
	double transform_[6];
	bool transformNotLoaded;
	int xsize_, ysize_;
	T min_, max_, mean_, stddev_, nodatavalue_;
	T *data;
	CPLErr loadStats()
	{
		CPLErr err;
		double min__, max__, mean__, stddev__;
		err = band->GetStatistics(0 /* false */, 1 /* true */, &min__, &max__, &mean__, &stddev__);
		min_ = min__; max_ = max__; mean_ = mean__; stddev_ = stddev__;
		return err;
	}

public:
	explicit Raster(const char *filename, int bandindex = 1, bool setGlobalProjection = false) :
	        dataset(static_cast<GDALDataset *>(GDALOpen(filename, GA_ReadOnly))),
	        band(0),
	        rasterLoaded_(false),
	        statsNotLoaded(true),
	        transformNotLoaded(true),
	        data(0)
	{
		if(dataset) {
			band = dataset->GetRasterBand(bandindex);
			xsize_ = band->GetXSize();
			ysize_ = band->GetYSize();
			int success;
			nodatavalue_ = static_cast<T>(band->GetNoDataValue(&success));

			if(setGlobalProjection) {
				dataset->GetGeoTransform(globalProjection.globalTransformation);
				globalProjection.globalProjection = dataset->GetProjectionRef();
				globalProjection.initialized = true;
				printf("setting global projection/transform: \"%s\" / %f %f %f %f %f %f\n", globalProjection.globalProjection.c_str(),
				       globalProjection.globalTransformation[0], globalProjection.globalTransformation[1], globalProjection.globalTransformation[2],
				       globalProjection.globalTransformation[3], globalProjection.globalTransformation[4], globalProjection.globalTransformation[5]);
			}
			rasterLoaded_ = true;
		}
	}

	~Raster()
	{
		if(data)
			delete [] data;
		GDALClose(dataset);
	}

	bool rasterLoaded() const { return rasterLoaded_; }

	int xsize() const { return xsize_; }
	int ysize() const { return ysize_; }
	T nodatavalue() const { return nodatavalue_; }
	//transform doc in http://www.gdal.org/classGDALDataset.html#f9593cc241e7d140f5f3c4798a43a668
	double transform(int index) {
		if(globalProjection.initialized) {
			return globalProjection.globalTransformation[index];
		}
		if(transformNotLoaded) {
			dataset->GetGeoTransform(transform_);
			transformNotLoaded = false;
		}
		return transform_[index];
	}

	//stats
	T min()
	{
		if(statsNotLoaded)
			loadStats();
		return min_;
	}

	T max()
	{
		if(statsNotLoaded)
			loadStats();
		return max_;
	}

	T mean()
	{
		if(statsNotLoaded)
			loadStats();
		return mean_;
	}

	T stddev()
	{
		if(statsNotLoaded)
			loadStats();
		return stddev_;
	}

	CPLErr read(T *data, int linespace = 0)
	{
		return band->RasterIO(GF_Read, 0, 0, xsize_, ysize_, data, xsize_, ysize_, GDALTypeTraits<T>::GDALType, 0, linespace == 0 ? 0 : GDALTypeTraits<T>::size * linespace);
	}

	T value(int x, int y)
	{
		if(!data) {
			data = new T[xsize_ * ysize_];
			band->RasterIO(GF_Read, 0, 0, xsize_, ysize_, data, xsize_, ysize_, GDALTypeTraits<T>::GDALType, 0, 0);
		}
		return data[x + y * xsize_];
	}
};

struct NoProjection {};
static const NoProjection noProjection = NoProjection();

template <typename T>
struct Dataset {
	GDALDataset *dataset;
	// new
	Dataset(const char *driver, const char *filename, int xsize, int ysize, int bands, NoProjection) :
	        dataset(GetGDALDriverManager()->GetDriverByName(driver)->Create(filename, xsize, ysize, bands, GDALTypeTraits<T>::GDALType, 0))
	{
	}
	Dataset(const char *driver, const char *filename, int xsize, int ysize, int bands) :
	        dataset(GetGDALDriverManager()->GetDriverByName(driver)->Create(filename, xsize, ysize, bands, GDALTypeTraits<T>::GDALType, 0))
	{
		if(RasterBase::globalProjection.initialized) {
			dataset->SetGeoTransform(RasterBase::globalProjection.globalTransformation);
			dataset->SetProjection(RasterBase::globalProjection.globalProjection.c_str());
		}
	}
	// copy
	template <typename U>
	Dataset(const char *driver, const char *filename, Dataset<U>& src, char **papszOptions = 0) :
	        dataset(GetGDALDriverManager()->GetDriverByName(driver)->CreateCopy(filename, src.dataset, 0/*FALSE*/, papszOptions, 0, 0))
	{
	}
        ~Dataset() { if (dataset) GDALClose(dataset); }

	void writeBand(int band, T *data, int linespace = 0)
	{
		if(dataset) {
			int xsize = dataset->GetRasterXSize(), ysize = dataset->GetRasterYSize();
			dataset->GetRasterBand(band)->RasterIO(GF_Write, 0, 0, xsize, ysize, data, xsize, ysize, GDALTypeTraits<T>::GDALType, 0, linespace == 0 ? 0 : GDALTypeTraits<T>::size * linespace);
		}
	}

	GDALDataset *operator->() { return dataset; }
};

#endif // RASTER_H
