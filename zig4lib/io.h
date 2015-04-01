#ifndef ZIG4IO_H
#define ZIG4IO_H

#include <QMetaType>
#include <QStringList>
#include <QFile>
#include "format.h"
#include "zig4lib_global.h"
#include "error_util.h"
#include "filepath.h"
#include "enum_util.h"

class QIODevice;
class FilePath;

struct ZUnknownFile;
struct ZRasterFile;
struct ZSPPFile;
struct ZBQPFile;
struct ZSSICoordinatesFile;
struct ZSSIFile;
struct ZPPAFile;
struct ZGroupsFile;
struct ZTreeFile;
struct ZInteractionsFile;
struct ZMatrixFile;
struct ZADMUFile;
struct ZIndexedRastersFile;
struct ZIGFile;
struct ZDATFile;
struct ZBATFile;
struct ZCommandLine;

// Error callback are based on polymorphic/virtual () operators.
// They MUST be passed as pointers or by reference!

bool ZIG4LIBSHARED_EXPORT loadZUnknownFile(ZUnknownFile& unknownFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZRasterFile(ZRasterFile& rasterFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZSPPFile(ZSPPFile& sppFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZSPPFile(const ZSPPFile& sppFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZBQPFile(ZBQPFile& bqpFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZBQPFile(const ZBQPFile& bqpFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZSSICoordinatesFile(ZSSICoordinatesFile& ssiCoordinatesFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZSSICoordinatesFile(const ZSSICoordinatesFile& ssiCoordinatesFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZSSIFile(ZSSIFile& ssiFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZSSIFile(const ZSSIFile& ssiFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZPPAFile(ZPPAFile& ppaFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZPPAFile(const ZPPAFile& ppaFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZGroupsFile(ZGroupsFile& groupsFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZGroupsFile(const ZGroupsFile& groupsFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZTreeFile(ZTreeFile& treeFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZTreeFile(const ZTreeFile& treeFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZInteractionsFile(ZInteractionsFile& interactionsFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZInteractionsFile(const ZInteractionsFile& interactionsFile, QIODevice& device, ErrorCallback&);

bool ZIG4LIBSHARED_EXPORT loadZMatrixFile(ZMatrixFile& matrixFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZMatrixFile(const ZMatrixFile& matrixFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZADMUFile(ZADMUFile& admuFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZADMUFile(const ZADMUFile& admuFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZIndexedRastersFile(ZIndexedRastersFile& indexedRastersFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZIndexedRastersFile(const ZIndexedRastersFile& indexedRastersFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZIGFile(ZIGFile& igFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZIGFile(const ZIGFile& igFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZDATFile(ZDATFile& datFile, QIODevice& device, ErrorCallback& err);
bool ZIG4LIBSHARED_EXPORT saveZDATFile(const ZDATFile& datFile, QIODevice& device, ErrorCallback& err);

bool ZIG4LIBSHARED_EXPORT loadZBATFile(ZBATFile& batFile, QIODevice& device, ErrorCallback& err);
// Used for the "arg. forwarding mode":
//bool ZIG4LIBSHARED_EXPORT loadZBATFile(ZBATFile& batFile, FilePath const& file, ErrorCallback err = DefaultErrorCallback());
bool ZIG4LIBSHARED_EXPORT saveZBATFile(const ZBATFile& batFile, QIODevice& device, ErrorCallback& err);

// input: zig4 -r set.dat ...
bool ZIG4LIBSHARED_EXPORT parseZCommandLine(ZCommandLine& command, const QStringList& arguments, ErrorCallback& err);

// templates for convenience

template <typename Function, typename Object>
bool loadFile(Function loadFunction, Object& object, const FilePath& path, ErrorCallback& err)
{
  QFile file(path);
  if(!file.open(QIODevice::ReadOnly)) {
    err(QString("unable to open file %1").arg(path));
    return false;
  }
  return loadFunction(object, file, err);
}

/*
// NOTE: duplicated from above - with the prototype required by zig4load
typedef boost::function<void (const QString& error)> TF_ErrorCallback;
class TF_DefaultErrorCallback
{
public:
	void operator()(const QString& str)
	{ qDebug() << str; };
};
template <typename Function, typename Object>
bool loadFile(Function loadFunction, Object& object, const FilePath& path, TF_ErrorCallback err = TF_DefaultErrorCallback())
{
  QFile file(path);
  if(!file.open(QIODevice::ReadOnly)) {
    err(QString("unable to open file %1").arg(path));
    return false;
  }
  return loadFunction(object, file, err);
}
*/

/*
// Used for the "arg. forwarding mode":
// specialization to bat file (does not have a qiodevice interface)
template <>
inline bool loadFile<bool (*)(ZBATFile&, FilePath const&, ErrorCallback), ZBATFile>
(bool (*)(ZBATFile&, FilePath const&, ErrorCallback), ZBATFile& object, const FilePath& path, ErrorCallback err)
{
    return loadZBATFile(object, path, err);
}
*/

template <typename Function, typename Object>
bool saveFile(Function saveFunction, const Object& object, const FilePath& path, ErrorCallback& err)
{
  QFile file(path);
  if(!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
    err(QString("unable to open file %1").arg(path), zeUndef);
    return false;
  }
  return saveFunction(object, file, err);
}

// "not to be confused with" ZCommandLine...
struct CommandLine
{
	QString program;
	QStringList arguments;
};

Q_DECLARE_METATYPE(CommandLine)

CommandLine ZIG4LIBSHARED_EXPORT getZCommandLineArguments(const ZCommandLine& command);

// specializations

/* fedemp: this funny pattern will call -through template instantiation- the loadZSPPFile, loadZADMUFile, loadZGroupsFile, 
   etc. functions and their 'save' counterpart. These are defined in the io_spp_load.cpp, io_admu_load.cpp, 
   io_gropus_load.cpp, etc. and the correspoding _save.cpp files. There you go!
*/
#define Z_IO_FORMATTRAITS(pod, podname) \
template <> \
struct FormatTraits<Z##pod##Format> \
{ \
	typedef Z##pod##File Pod; \
	static const QString& name() \
	{ \
		static const QString name_( podname ); \
		return name_; \
	} \
	static bool load(Pod& file, QIODevice& device, ErrorCallback& err) \
	{ \
		return loadZ##pod##File(file, device, err); \
	} \
	static bool save(const Pod& file, QIODevice& device, ErrorCallback& err) \
	{ \
		return saveZ##pod##File(file, device, err); \
	} \
}

#define Z_IO_DEFAULTFORMATTRAITS(pod, podname) \
struct Z##pod##Format {}; \
Z_IO_FORMATTRAITS(pod, podname);

Z_IO_DEFAULTFORMATTRAITS(BQP, "zbqp/text");
Z_IO_DEFAULTFORMATTRAITS(PPA, "zppa/text");
Z_IO_DEFAULTFORMATTRAITS(SPP, "zspp/text");
Z_IO_DEFAULTFORMATTRAITS(DAT, "zdat/text");
Z_IO_DEFAULTFORMATTRAITS(BAT, "zbat/text");
Z_IO_DEFAULTFORMATTRAITS(ADMU, "zadmu/text");
Z_IO_DEFAULTFORMATTRAITS(IndexedRasters, "zindexedrasters/text");
Z_IO_DEFAULTFORMATTRAITS(IG, "zig/text");
Z_IO_DEFAULTFORMATTRAITS(SSI, "zssi/text");
Z_IO_DEFAULTFORMATTRAITS(SSICoordinates, "zssicoordinates/text");
Z_IO_DEFAULTFORMATTRAITS(Groups, "zgroups/text");
Z_IO_DEFAULTFORMATTRAITS(Tree, "ztree/text");
Z_IO_DEFAULTFORMATTRAITS(Interactions, "zinteractions/text");
Z_IO_DEFAULTFORMATTRAITS(Matrix, "zmatrix/text");

// specializations

struct ZRasterFormat {};
template <>
struct FormatTraits<ZRasterFormat>
{
	typedef ZRasterFile Pod;
	static const QString& name()
	{
		static const QString name_("raster/bin");
		return name_;
	}
	static bool load(Pod& file, QIODevice& device, ErrorCallback& err)
	{
		return loadZRasterFile(file, device, err);
	}
	static bool save(const Pod& file, QIODevice& device, ErrorCallback& err)
	{
		return false;
		//return saveZRasterFile(file, device, err);
	}
};

// The traditional ZBatFile struct is defined in pod.h

/*
// Used for the "arg. forwarding mode":
struct ZBATFormat {};

template <>
struct FormatTraits<ZBATFormat>
{
	typedef ZBATFile Pod;
	static const QString& name()
	{
		static const QString name_("zbat/text");
		return name_;
	}
	static bool load(Pod& file, QIODevice& device, ErrorCallback err = DefaultErrorCallback())
	{
		QFile *f(qobject_cast<QFile *>(&device));
		Q_ASSERT(f);
		return loadZBATFile(file, f->fileName(), err);
	}
	static bool save(const Pod& file, QIODevice& device, ErrorCallback err = DefaultErrorCallback())
	{
		return saveZBATFile(file, device, err);
	}
};
*/

struct ZUnknownFormat {};
template <>
struct FormatTraits<ZUnknownFormat>
{
	typedef ZUnknownFile Pod;
	static const QString& name()
	{
		static const QString name_("unknown/bin");
		return name_;
	}
	static bool load(Pod& file, QIODevice& device, ErrorCallback& err)
	{
		return true;
	}
	static bool save(const Pod& file, QIODevice& device, ErrorCallback& err)
	{
		return false;
	}
};

#endif // ZIG4IO_H
