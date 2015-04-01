#ifndef PROJ_IO_UTIL_H
#define PROJ_IO_UTIL_H

#include "manager.h"
#include "proj.h"

class LoadOptionalFile
{
 private:
  DirPath root;
  FileManager& manager;
  // MUST be a reference - virtual operator()
  ErrorCallback& err;

 public:
 LoadOptionalFile(DirPath& root, FileManager& manager, ErrorCallback& err): 
  root(root), manager(manager), err(err) 
  {};

  template <typename Format>
    bool operator()(FileRef<Format>& file, const FilePath& path)
    {
      if(!path.isEmpty()) {
	file = manager.getLinked<Format>(root.absoluteFilePath(path), err);
	return file;
      }
      return true;
    }
  
  template <typename Format>
    bool operator()(FileRef<Format>& file, const boost::optional<FilePath>& path)
    {
      if (path) 
	return operator()(file, *path);
      else 
	return true;
    }  
};

#endif // PROJ_IO_UTIL_H
