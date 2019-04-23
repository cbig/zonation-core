#ifndef MANAGER_H
#define MANAGER_H

#include <typeinfo>

#include <QObject>
#include <QTextStream>
#include <QFile>
// fedemp 20120329: anything related to fs watch commented out
// #include <QFileSystemWatcher>

#include <boost/shared_ptr.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/make_shared.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/unordered_set.hpp>

#include "filepath.h"
#include "format.h"
#include "string_util.h"
#include "zig4lib_global.h"
#include "error_util.h"

struct PathTag {};
struct TypeTag {};
struct PtrTag {};
struct RefTag {};
struct FileTag {};
struct ManagerTag {};

struct PodTag {};

class FileManager;
class ZProject;

class File : public QObject
{
	Q_OBJECT

public:
	enum Status {
		UNLINKED,
		SYNCED,
		UNSYNCED,
		PARSE_ERROR, // also unsynced
		FILE_ERROR // also unsynced
	};

	FilePath path;
	Status status_;

	int refCount;
	FileManager& manager;

	File(const FilePath& path, FileManager& manager) :
	        path(path),
	        status_(UNLINKED),
	        refCount(0),
	        manager(manager)
	{
	}

	virtual ~File();

	virtual bool load(ErrorCallback& err) = 0;
	virtual bool save(ErrorCallback& err) = 0;
	virtual const std::type_info& type() const = 0;
	virtual const QString& formatName() const = 0;

	Status status() const { return status_; }

	void slotFileChanged();

signals:
	void fileChanged();
};

void ZIG4LIBSHARED_EXPORT intrusive_ptr_add_ref(File *p);
void ZIG4LIBSHARED_EXPORT intrusive_ptr_release(File *p);

template <typename Format>
class FileImpl : public File
{
public:
	typedef typename FormatTraits<Format>::Pod Pod;

	FileImpl(const FilePath& path, FileManager& manager);
	virtual ~FileImpl();
	virtual bool load(ErrorCallback& err);
	virtual bool save(ErrorCallback& err);
	virtual const std::type_info& type() const;
	virtual const QString& formatName() const;

	Pod pod;
};

typedef boost::multi_index_container<
File *,
boost::multi_index::indexed_by<
boost::multi_index::hashed_unique<
boost::multi_index::tag<PathTag>,
boost::multi_index::member<File, FilePath, &File::path>
>,
boost::multi_index::hashed_unique<
boost::multi_index::tag<FileTag>,
boost::multi_index::identity<File *>
>
>
> FileMap;


template <typename Format>
class FileRef
{
public:
	typedef typename FormatTraits<Format>::Pod Pod;

	// construct null
	FileRef() {}

	FileRef(File *file) :
	        file(static_cast<FileImpl<Format> *>(file))
	{}

	FileRef(const FilePath& path, FileManager& manager) :
	        file(new FileImpl<Format>(path, manager))
	{}

	Pod& pod() const
	{
	  if(file->status_ == File::SYNCED) {
	    file->status_ = File::UNSYNCED;
	  }
	  return file->pod;
	}

	const Pod& constPod() const
	{
	  return file->pod;
	}

	const Pod& operator*() const
	{
	  return constPod();
	}

	const Pod *operator->() const
	{
	  return &constPod();
	}

	Pod& operator*()
	{
	  return pod();
	}

	Pod *operator->()
	{
	  return &pod();
	}

	FilePath path() const
	{
	  return file->path;
	}

	operator bool() const
	{
	  return file != nullptr;
	}

	int refCount() const
	{
	  if(file)
	    return file->refCount;
	  return 0;
	}

	boost::intrusive_ptr<FileImpl<Format> > file;
};

class FileChangedCallback
{
 public:
  /* Interface for the method that reacts to a change in any file (belonging to project ?) */
  virtual void doReact() const = 0;
};

class ZIG4LIBSHARED_EXPORT FileManager : public QObject
{
  Q_OBJECT

public:

  FileManager( /*const FileChangedCallback& fccb */);

  // return null, if type error
  template <typename Format>
    FileRef<Format> get(const FilePath& path, ErrorCallback& err)
    {
      using namespace boost;
      // get absolute path
      FilePath absolutePath(path.absoluteFilePath());
      /*
      FileMap::iterator i(map.find(absolutePath));
      FileRef<Format> ret;
      if(i == map.end()) {
	ret = FileRef<Format>(absolutePath, *this);
	map.insert(ret.file.get());
	emit fileAdded(ret.file.get());
      } else if((*i)->type() == typeid(Format)) {
	ret = FileRef<Format>(*i);
      } else {
	ZError::err_msg(tr("Error loading file %1 as %2. File has already been loaded as %3.")
			.arg(path).arg(FormatTraits<Format>::name())
			.arg((*map.find(path.absoluteFilePath()))->formatName()),
			zeError);
      }
      */
      // just load it!
      FileRef<Format> ret = FileRef<Format>(absolutePath, *this);
      return ret;
    }

  // return null, if type error or if unlinked file could not be loaded
  template <typename Format>
    FileRef<Format> getLinked(const FilePath& path, ErrorCallback& err)
    {
      FileRef<Format> ref(get<Format>(path, err));
      // if(ref && ref.file->status() == File::UNLINKED) {
      if(ref) {
	if(!ref.file->load(err)) {
	  return FileRef<Format>();
	}
      }
      return ref;
    }


  // void setOpeningProject(boost::shared_ptr<ZProject> p)
  // { project_being_opened = p; };

  // void watch(File *file);
  // void remove(File *file);

  // FileMap map;
  // QFileSystemWatcher watcher;

 public slots:
  // void slotFileChanged(const QString& path) const;

 signals:
  void fileAdded(File *file);
  // void fileRemoved(File *file);

 private:
  // This would be very needed to know to which project a file belongs
  // then, when files are changed, we know which project must be updated/reloaded
  // boost::shared_ptr<ZProject> project_being_opened;

  // const FileChangedCallback& fccb;
};

// implementations

template <typename Format>
FileImpl<Format>::FileImpl(const FilePath& path, FileManager& manager) :
File(path, manager)
{
}

template <typename Format>
FileImpl<Format>::~FileImpl() {}

template <typename Format>
bool FileImpl<Format>::load(ErrorCallback& err)
{
  QFile file(path);
  if(!file.open(QIODevice::ReadOnly)) {
    //err(QObject::tr("*** Could not open file %1 for reading").arg(path));
    status_ = FILE_ERROR;
    return false;
  }

  // watch this file
  // manager.watch(this);

  //LineModify moderr(err, path + ": ");
  if(!FormatTraits<Format>::load(pod, file, err)) {
    ZError::err_msg(QObject::tr("Could not read file %1").arg(path), zeError);
    status_ = PARSE_ERROR;
    return false;
  }
  status_ = SYNCED;
  return true;
}

template <typename Format>
bool FileImpl<Format>::save(ErrorCallback& err)
{
  QFile file(path);
  if(!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
    ZError::err_msg(QObject::tr("Could not open file %1 for writing").arg(path), zeError);
    status_ = FILE_ERROR;
    return false;
  }
  LineModify moderr(err, path + ": ");
  if(!FormatTraits<Format>::save(pod, file, moderr)) {
    ZError::err_msg(QObject::tr("Could not write file %1").arg(path), zeError);
    status_ = PARSE_ERROR;
    return false;
  }
  status_ = SYNCED;
  return true;
}

template <typename Format>
const std::type_info& FileImpl<Format>::type() const
{
  return typeid(Format);
}

template <typename Format>
const QString& FileImpl<Format>::formatName() const
{
  return FormatTraits<Format>::name();
}


#endif // MANAGER_H
