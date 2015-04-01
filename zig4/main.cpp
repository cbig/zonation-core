// Zonation computational core (zig4)
// Copyright (C) 2004-2013  Atte Moilanen
//
// This file and all the files contained in this directory and its 
// subdirectories are part of Zonation computational core (zig4).
//
// zig4 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// zig4 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with zig4.  If not, see <http://www.gnu.org/licenses/>.

//#include "util/zlog.h"
#include "LoadData.h"
#include "Unit1.h"
#include "zig4lib/io.h"
#include <QCoreApplication>
#include <iostream>
#include <QDebug>
#include <gdal/gdal_priv.h>
#if defined(_MSC_VER) && defined(_CRTDBG_MAP_ALLOC)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

// These includes below are needed for the "forward args" mode (env variable ZIG_FORWARD_ARGS).
#include "zig4lib/env.h"
#include <QProcessEnvironment>
#include <QLocalSocket>
#include <QTime>

void init()
{
	// init gdal
	GDALAllRegister();
	// disable aux files
	CPLSetConfigOption("GDAL_PAM_ENABLED", "NO");
	// fix ascii type to Float32 - gdal 1.8.0 option (hacked into gdal used by Z)
	// GDAL detects data type by probing the first 100KB in the ascii file looking for decimal points.
	// This could be a source of extremely obscure bugs, since Z often uses integer nodata values.
	CPLSetConfigOption("AAIGRID_DATATYPE", "Float32");
}

// abusing qt debug system...
void messageOutput(QtMsgType type, const char *msg)
{
	//uncomment to enable debugging messages
	//printf("Debug: %s\n", msg);
}

void checkArgMode(QCoreApplication const& app)
{
	QProcessEnvironment env(QProcessEnvironment::systemEnvironment());
	if(env.contains(ZEnvironment::forwardArgsVariable)) {
	  std::cerr << "argument forwarding mode\n" << std::endl;
	  QLocalSocket socket;
	  QString serverName(env.value(ZEnvironment::argServerNameVariable));
	  if(serverName.isNull()) {
	    std::cerr << QObject::tr("environment variable %1 not set").arg(ZEnvironment::argServerNameVariable).toStdString() << endl;
	    exit(-1);
	  }
	  std::cerr << "opening socket..." << std::endl;
	  socket.connectToServer(serverName, QIODevice::WriteOnly);
	  if(!socket.waitForConnected(ZEnvironment::timeoutMsecs)) {
	    std::cerr << socket.errorString().toStdString() << std::endl;
	    exit(-1);
	  }
	  std::cerr << "opening stream..." << std::endl;
	  QDataStream stream(&socket);
	  std::cerr << "writing application data..." << std::endl;
	  // QCoreApplication::applicationFilePath() not reliable on unix!
	  stream << ArgModeData(QCoreApplication::applicationFilePath(),
				QDir::currentPath(),
				app.arguments());
	  std::cerr << "flushing data..." << std::endl;
	  while(socket.bytesToWrite() > 0) {
	    if(!socket.waitForBytesWritten(ZEnvironment::timeoutMsecs)) {
	      std::cerr << socket.errorString().toStdString() << std::endl;
	      exit(-1);
	    }
	  }
	  
		// waitForDisconnected() emits an Unknown error on Windows 7 x64 / Qt 4.7.4
		// even though error doesn't occur. We'll just ignore the return value
		// from waitForDisconnected()...

	  std::cerr << "disconnecting socket..." << QTime::currentTime().toString().toStdString() << std::endl;
	  socket.disconnectFromServer();
	  std::cerr << "waiting for disconnect..." << QTime::currentTime().toString().toStdString() << std::endl;
	  socket.waitForDisconnected(ZEnvironment::timeoutMsecs);
		/*
	if(!socket.waitForDisconnected(ZEnvironment::timeoutMsecs)) {
	    std::cerr << "blimey..." << QTime::currentTime()..toStdString() << std::endl;
	    std::cerr << socket.errorString().toStdString() << std::endl;
	    exit(-1);
	}
	*/
	  std::cerr << "quitting application..." << std::endl;
	  exit(0);
	}
}

int main(int argc, char *argv[])
{
	qInstallMsgHandler(messageOutput);
	QCoreApplication app(argc, argv);

#if defined(_MSC_VER) && defined(_DEBUG)
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

        // This could be an option (would be defined in zig4lib/pod.h and so on).
        // But most likely there will not be much interest in it. 
        // For now, it's just hidden here, and this arg. fw functionality is never used.
	// --> traditional simple Z-specific parsing of .bat files is always used.
        bool enable_advancedBatchMode = false;
	// Check whether it's running in "arguments forwarding mode"/"advanced batch mode"
	// If it was activated and the option set, checkArgMode would "intercept" this
	// zonation run (explanations in env.h).
	if (enable_advancedBatchMode)
	  checkArgMode(app);
	else
	  qDebug() << "main(): argument forwarding mode disabled" << endl << endl;

	BatRun run;
	try {
	  init();
	  run.run();
	} catch(std::exception& ex) {
	  std::cerr << ex.what() << std::endl;
	} catch(...) {
	  std::cerr << "uncaught exception" << std::endl;
	}
	deleteSingletons();
	return run.retval;
}
