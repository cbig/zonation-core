include(../common.pri)

TARGET = zig4lib

TEMPLATE = lib

CONFIG += qt dll thread

QT -= gui
QT += core network

DESTDIR = $$BIN_PATH

DEFINES += ZIG4LIB_LIBRARY #_HAS_ITERATOR_DEBUGGING=0

INCLUDEPATH += \
	$$BOOST_INCLUDEPATH \
	$$GDAL_INCLUDEPATH
LIBS += \
	$$BOOST_LIBS \
	$$GDAL_LIBS

QMAKE_RPATHDIR += \
	$$BOOST_LIBPATH \
	$$GDAL_LIBPATH \

SOURCES += \
	manager.cpp \
	ini_load.cpp \
	ini_save.cpp \
	io.cpp \
	spirit_util.cpp \
	proj.cpp \
	proj_io.cpp \
	core.cpp \
	observer.cpp \
	raster.cpp \
	io_spp_load.cpp \
	io_ppa_load.cpp \
	io_ppa_save.cpp \
	io_spp_save.cpp \
	io_raster_load.cpp \
	io_bqp_load.cpp \
	io_bqp_save.cpp \
	io_ssicoord_load.cpp \
	io_ssicoord_save.cpp \
	io_ssi_load.cpp \
	io_ssi_save.cpp \
	io_groups_load.cpp \
	io_groups_save.cpp \
	io_tree_load.cpp \
	io_tree_save.cpp \
	io_interactions_save.cpp \
	io_interactions_load.cpp \
	io_matrix_load.cpp \
	io_matrix_save.cpp \
	io_admu_save.cpp \
	io_admu_load.cpp \
	io_indexedrasters_load.cpp \
	io_indexedrasters_save.cpp \
	io_ig_load.cpp \
	io_ig_save.cpp \
	io_arg.cpp \
    zconf.cpp

HEADERS += \
	manager.h \
	filepath.h \
	ini.h \
	io.h \
	pod.h \
	poqstring.h \
	qstring_spirit.h \
	spirit_util.h \
	split_winmain.h \
	exception.h \
	format.h \
	proj.h \
	string_util.h \
	enum_util.h \
	proj_io.h \
	zig4lib_global.h \
	error_util.h \
	core.h \
	observer.h \
	zinterprocess.h \
	raster.h \
	ppa_adapt.h \
	ptr_util.h \
	spirit_error.h \
	spirit_enum.h \
	env.h \
	io_arg.h \
    gisnodata.h \
    zconf.h \
    proj_io_util.h







