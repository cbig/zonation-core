#ifndef ZIG4LIB_GLOBAL_H
#define ZIG4LIB_GLOBAL_H

#include <QtCore/qglobal.h>
//#include <zig4lib/version.h>

// static compilation on Win using MinGW
#if defined(__MINGW32__) || defined(__MINGW64__)
#  define ZIG4LIBSHARED_EXPORT 
#else
#  if defined(ZIG4LIB_LIBRARY)
#    define ZIG4LIBSHARED_EXPORT Q_DECL_EXPORT
#  else
#    define ZIG4LIBSHARED_EXPORT Q_DECL_IMPORT
#endif
#endif

#include <cmath>
#ifndef _MSC_VER
using std::isnan;
#   define z_isnan isnan
#else
#pragma message {"Using _isnan in MSVC"}
#define isnan _isnan
#endif

#endif // ZIG4LIB_GLOBAL_H
