# Depends on: boost, QtCore, QtNetwork

file(GLOB ZIG4LIB_SRCS *.cpp)

include(CheckIncludeFiles)
check_include_files(malloc.h HAVE_MALLOC_H)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/zig_config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/zig4lib/zig_config.h")
include_directories(${PROJECT_BINARY_DIR})

add_definitions(-DZIG4LIB_LIBRARY)

# http://www.qtcentre.org/wiki/index.php?title=Compiling_Qt4_apps_with_CMake
# find qt, and include qt macros
find_package(Qt4 REQUIRED QtCore QtNetwork)
include(${QT_USE_FILE})

file(GLOB ZIG4LIB_MOC_HDRS core.h io_arg.h manager.h observer.h)
qt4_wrap_cpp(ZIG4LIB_MOC_SRCS ${ZIG4LIB_MOC_HDRS})

add_library(zig4lib zig_config.h ${ZIG4LIB_SRCS} ${ZIG4LIB_MOC_SRCS})

find_package(Qt4 REQUIRED QtCore QtGui QtNetwork)
include(${QT_USE_FILE})

if (NOT MINGW)
   find_package(GDAL REQUIRED)
      if(GDAL_FOUND)
   include_directories(${GDAL_INCLUDE_DIR})
   endif()
else()
   include_directories(/mingw/include/gdal)
endif()

# Needed for the enum/ dir included in the Z sources
include_directories(${CMAKE_SOURCE_DIR})
