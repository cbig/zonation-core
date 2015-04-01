#ifndef PROJ_IO_H
#define PROJ_IO_H

#include "proj.h"
#include "zig4lib_global.h"
#include "error_util.h"

boost::shared_ptr<ZInstance> ZIG4LIBSHARED_EXPORT loadZInstance(const boost::shared_ptr<ZProject>& project, int row, const ZCommandLine& command, FileManager& manager, ErrorCallback& e);
boost::shared_ptr<ZProject> ZIG4LIBSHARED_EXPORT loadZProject(const FilePath& batFilePath, FileManager& manager, ErrorCallback& e);

#endif // PROJ_IO_H
