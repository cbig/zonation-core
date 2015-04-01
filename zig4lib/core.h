#ifndef CORE_H
#define CORE_H

#include "zig4lib_global.h"
#include "zinterprocess.h"
#include "filepath.h"
#include <QThread>
#include <QLocalServer>
#include <QFile>
#include <QTextStream>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace ZInterprocess {

class Core;

class CoreLoop : public QObject
{
	Q_OBJECT

public:
	CoreLoop(const QString& serverName, Core *core);

	void notifyCurrent(QLocalSocket *socket);

	QString serverName;
	QLocalServer *server;
	Core *core;

	QSet<QLocalSocket *> conns;

	QByteArray reusableBuffer;
	QByteArray sizeBuffer;

public slots:
	void init();
	void newConnection();
	void disconnected();
	void readyRead();
	void sendBytes(const QByteArray& bytes);
	void close();
	void quit();
};

class CoreUnicast
{
public:
	CoreLoop *loop;
	QLocalSocket *socket;
	QByteArray *array;
	QDataStream stream;
	CoreUnicast(CoreLoop *loop, QLocalSocket *socket, QByteArray *arr);
	~CoreUnicast();
};

class Core : public QThread
{
	Q_OBJECT

public:
  Core(const FilePath& outFilePath, QString& resultMsg);
  ~Core();

public:

	ZIG4LIBSHARED_EXPORT void init(const InstanceInitData& instance);
	ZIG4LIBSHARED_EXPORT void setPixel(quint32 x, quint32 y, quint32 color);
	ZIG4LIBSHARED_EXPORT void notifyInit(); // send init notification
	ZIG4LIBSHARED_EXPORT void addPlotPoint(quint32 plotIndex, float x, float y);
	ZIG4LIBSHARED_EXPORT void removeSite(quint32 x, quint32 y, quint32 color, float rank);
	ZIG4LIBSHARED_EXPORT void notifyDone();
	ZIG4LIBSHARED_EXPORT void message(const QString& msg);

	//ZIG4LIBSHARED_EXPORT const quint32 *colorMapData();

	void emitSendBytes(const QByteArray& buffer);

	ZIG4LIBSHARED_EXPORT QString compileMsgList();

	boost::shared_ptr<CoreLoop> loop;
	boost::shared_ptr<boost::interprocess::file_lock> fileLock;
	QFile outFile;
	QTextStream out;

	// local data
	QMutex mutex;
	boost::shared_ptr<InstanceData> data;
	QList<QString> msgList;
	CoreStatus status;

	QByteArray reusableBuffer;

signals:
	void init();
	void sendBytes(const QByteArray& bytes);
	void closeSignal();
};

class CoreBroadcast
{
public:
	Core *core;
	QMutexLocker lock;
	QByteArray *array;
	QDataStream stream;

	CoreBroadcast(Core *core, QByteArray *buffer);
	~CoreBroadcast();
};
}

ZIG4LIBSHARED_EXPORT ZInterprocess::Core *createCore(const FilePath& outFilePath, QString& resultMsg);
ZIG4LIBSHARED_EXPORT void releaseCore(ZInterprocess::Core *core);

#endif // CORE_H
