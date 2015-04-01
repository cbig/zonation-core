#ifndef Z_MEM_DUMP_H
#define Z_MEM_DUMP_H

#include "typedefs.h"

String make_vmat_output_name(const String& proj_name, const String& fn, bool compress=true);

// load_direct: if true, the file name in fn is a full path; otherwise it is used as a suffix to the project name (proj_name)
bool load_vmat(bool load_direct, const String& proj_name, const String& fn, size_t xdim, size_t ydim, int map_cnt, bool compress=true);

bool save_vmat(const String& proj_name, const String& fn, size_t xdim, size_t ydim, int map_cnt, bool compress=true);

#endif // Z_MEM_DUMP_H
