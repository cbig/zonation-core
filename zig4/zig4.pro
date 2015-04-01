include(../common.pri)

TARGET = zig4

TEMPLATE = app

CONFIG += qt console thread
CONFIG -= app_bundle

QT += core network
QT -= gui

INCLUDEPATH += . \
	$$BOOST_INCLUDEPATH \
	$$GDAL_INCLUDEPATH \
	$$FFTW_INCLUDEPATH \
	dummy_headers \
	core

LIBS += \
	$$ZIG4_LIBS \
	$$BOOST_LIBS \
	$$GDAL_LIBS \
	$$FFTW_LIBS

QMAKE_RPATHDIR += \
	$$ZIG4_LIBPATH \
	$$BOOST_LIBPATH \
	$$GDAL_LIBPATH

#QMAKE_RPATHDIR += ../../gdal-build/lib

DESTDIR = $$BIN_PATH

HEADERS += \
	VCL.h \
	core/Unit1.h \
	core/post_process.h \
	core/output.h \
	core/nrutils.h \
        core/randz.h \
	core/LoadData.h \
	core/GridMap.h \
	core/fftw_lnk.h \
	core/defines.h \
	core/marginal_loss.h \
	core/ADMUs.h \
	zglobal.h \
	core/LSIdent.h \
	core/bat_run.h \
	core/PLULA.h \
	core/smooth.h \
	core/typedefs.h \
        core/occur_container.h

SOURCES += \
	VCL.cpp \
	main.cpp \
	core/Unit1.cpp \
	core/smooth.cpp \
	core/post_process.cpp \
	core/PLULA.cpp \
	core/output.cpp \
	core/nrutils.c \
        core/randz.c \
	core/LSIdent.cpp \
	core/loaddata.cpp \
	core/gridmap.cpp \
	core/fftw_lnk.cpp \
	core/bat_run.cpp \
	core/marginal_loss.cpp \
	core/ADMUs.cpp \
	zglobal.cpp 
        
