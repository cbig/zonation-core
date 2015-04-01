#include "core.h"
#include <QDebug>
#include <QLocalSocket>
#include <QTimer>

namespace {
const int REUSABLE_BUFFER_SIZE = 8192;
}

namespace ZInterprocess {

using namespace boost;
using namespace boost::interprocess;

// CoreLoop

CoreLoop::CoreLoop(const QString& serverName, Core *core) :
        //QObject(core),
        serverName(serverName),
        server(0),
        core(core),
        reusableBuffer(REUSABLE_BUFFER_SIZE, 0),
        sizeBuffer(sizeof(MessageSize), 0)
{
	server = new QLocalServer(this);
	moveToThread(core);
}

void CoreLoop::init() {
	qDebug() << "CoreLoop::init()";
	QLocalServer::removeServer(serverName);
	connect(server, SIGNAL(newConnection()), this, SLOT(newConnection()));
	server->listen(serverName);
}

CoreUnicast::CoreUnicast(CoreLoop *loop, QLocalSocket *socket, QByteArray *arr) :
        loop(loop),
        socket(socket),
        array(arr),
        stream(array, QIODevice::WriteOnly)
{}

CoreUnicast::~CoreUnicast()
{
	MessageSize size = array->size();
	QDataStream sstream(&loop->sizeBuffer, QIODevice::WriteOnly);
	sstream << size;
	socket->write(loop->sizeBuffer);
	socket->write(*array);
}

void CoreLoop::notifyCurrent(QLocalSocket *socket)
{
	QMutexLocker lock(&core->mutex);
	QByteArray array;
	{
		CoreUnicast unicast(this, socket, &array);
		unicast.stream << static_cast<TagType>(CORE_CURRENT_MSGLIST_NOTIFY) << core->msgList;
	}
	if(core->status != CORE_NOT_INITIALIZED) {
		CoreUnicast unicast(this, socket, &array);
		unicast.stream << static_cast<TagType>(CORE_CURRENT_DATA_NOTIFY) << *core->data;
	}
	if(core->status == CORE_DONE) {
		CoreUnicast unicast(this, socket, &array);
		unicast.stream << static_cast<TagType>(CORE_DONE);
	}
}

void CoreLoop::newConnection()
{
	while(server->hasPendingConnections()) {
		QLocalSocket *socket(server->nextPendingConnection());
		QFile out;
		out.open(stdout, QIODevice::WriteOnly);
		qDebug() << socket << "connected";
		conns.insert(socket);
		connect(socket, SIGNAL(disconnected()), this, SLOT(disconnected()));
		connect(socket, SIGNAL(readyRead()), this, SLOT(readyRead()));
		notifyCurrent(socket);
	}
}

void CoreLoop::disconnected()
{
	QLocalSocket *socket(static_cast<QLocalSocket *>(sender()));
	QFile out;
	out.open(stdout, QIODevice::WriteOnly);
	qDebug() << socket << "disconnected";
	conns.remove(socket);
	socket->deleteLater();
}

void CoreLoop::readyRead()
{
  // do something about queries!
}

void CoreLoop::sendBytes(const QByteArray& bytes)
{
	MessageSize size(bytes.size());
	QDataStream stream(&sizeBuffer, QIODevice::WriteOnly);
	stream << size;
	foreach(QLocalSocket *socket, conns) {
		socket->write(sizeBuffer);
		socket->write(bytes);
	}
}

void CoreLoop::close()
{
	qDebug() << "CoreLoop::close()";
	server->close();

	QByteArray arr;
	QDataStream stream(&arr, QIODevice::WriteOnly);
	stream << static_cast<TagType>(CORE_CLOSING_NOTIFY);
	sendBytes(arr);
	QFile out;
	out.open(stdout, QIODevice::WriteOnly);

	foreach(QLocalSocket *socket, conns) {
		qDebug() << "waiting for" << socket << "to disconnect...";
		if(!socket->waitForDisconnected(-1)) {
			qDebug() << "error disconnecting" << socket;
		}
		socket->deleteLater();
	}
	QTimer::singleShot(0, this, SLOT(quit()));
}

void CoreLoop::quit()
{
	core->quit();
}

// Core

Core::Core(const FilePath& outFilePath, QString& resultMsg):
        outFile(outFilePath),
        status(CORE_NOT_INITIALIZED),
        reusableBuffer(REUSABLE_BUFFER_SIZE, 0)
{
  bool ok = outFile.open(QFile::WriteOnly | QFile::Truncate | QFile::Text);
  // Before giving up, check if the directory path is missing and can be created
  if(!ok) {
    // dirname in the POSIX sense
    QString dirname = outFile.fileName().left(outFile.fileName().lastIndexOf('/'));
    if (!dirname.isEmpty()) {
      QDir d(dirname);
      if (!d.exists()) {
	// no way to log properly from here
	resultMsg += "Output path '" + dirname + "' does not exist, trying to create it... ";
	if (QDir("./").mkpath(dirname) )
	  resultMsg += "succeeded.";
	else
	  resultMsg += "failed!";
	// output here, cause it's going to exit soon...
	std::cout << resultMsg.toStdString() << std::endl;
      }
    }
  }
  // do a second try if needed
  if(ok || outFile.open(QFile::WriteOnly | QFile::Truncate | QFile::Text)) {
    out.setDevice(&outFile);
    fileLock = make_shared<file_lock>(QDir::toNativeSeparators(outFilePath).toStdString().c_str());
    if(!fileLock->try_lock()) {
      throw Exception("can't lock file " + outFilePath);
    }
    loop = make_shared<CoreLoop>(getServerName(outFilePath), this);
    connect(this, SIGNAL(sendBytes(const QByteArray&)), loop.get(), SLOT(sendBytes(const QByteArray&)), Qt::QueuedConnection);
    connect(this, SIGNAL(closeSignal()), loop.get(), SLOT(close()), Qt::QueuedConnection);
    connect(this, SIGNAL(init()), loop.get(), SLOT(init()), Qt::QueuedConnection);
    start();
    emit init();
  } else {
    throw Exception("can't open file " + outFilePath);
  }
}

Core::~Core()
{
	qDebug() << "Core::~Core()";
	emit closeSignal();
	//quit();
	wait();
	fileLock->unlock();
	qDebug() << "exiting Core::~Core()";
}

  /*
const quint32 *Core::colorMapData()
{
  //return data->colorMap.constData();
  return data->colorMap.data(); // vector C++11, this is for the change from QVector<> to std::vector<> in zinterprocess.h:InstanceData
}
  */
void Core::init(const InstanceInitData& instanceInit)
{
	qDebug() << "Core::init()";
	data = make_shared<InstanceData>(instanceInit);
}

void Core::setPixel(quint32 x, quint32 y, quint32 color)
{
	data->setPixel(x, y, color);
	//qDebug() << "Core::setPixel()";
}

void Core::notifyInit()
{
	qDebug() << "Core::notifyInit()";
	QByteArray array; // big temporary buffer
	CoreBroadcast sender(this, &array);
	status = CORE_INITIALIZED;
	sender.stream << static_cast<TagType>(CORE_INITIALIZED_NOTIFY) << CoreInitializedData(data->init, data->colorMap);
	//writeContainer(sender.stream, data->colorMap);
}

void Core::addPlotPoint(quint32 plotIndex, float x, float y)
{
	CoreBroadcast sender(this, &reusableBuffer);
	quint32 offset(data->addPlotPoint(plotIndex, x, y));
	sender.stream << static_cast<TagType>(CORE_PLOT_POINT_ADDED_NOTIFY) << PlotPointAddedData(plotIndex, offset, PlotPoint(x, y));
	//qDebug() << "Core::addPlotPoint()";
}

void Core::removeSite(quint32 x, quint32 y, quint32 color, float rank)
{
	CoreBroadcast sender(this, &reusableBuffer);
	quint64 offset(data->removeSite(x, y, color, rank));
	sender.stream << static_cast<TagType>(CORE_SITE_REMOVED_NOTIFY) << SiteRemovedData(data->removedCount, offset, color, rank);
	//qDebug() << "Core::removeSite()" << data->removedCount;
}

void Core::notifyDone()
{
	qDebug() << "Core::notifyDone()";
	CoreBroadcast sender(this, &reusableBuffer);
	status = CORE_DONE;
	sender.stream << static_cast<TagType>(CORE_DONE_NOTIFY);
}

void Core::message(const QString& msg)
{
	CoreBroadcast sender(this, &reusableBuffer);
	////fprintf(stdout, msg.toLatin1().data());
	//fprintf(stdout, msg.toStdString().c_str());
	//fprintf(stdout, "\n");

	out << msg << '\n'; //endl;  // Don't use endl, something stops working in QTextStreams 
	                             // redirected to stdout on windows and you'll only get the first line
	// A flush would also break everything -- out.flush();
	// could try outFile.flush()

	// This accumulates messages for 'compileMsgList' which flushes them at the end
	msgList << msg;
	MessageData data(msgList.size(), msg);
	sender.stream << static_cast<TagType>(CORE_MSG_NOTIFY) << data;
}


QString Core::compileMsgList()
{
	QMutexLocker lock(&mutex);
	QString ret;
	QTextStream s(&ret);
	foreach(const QString& msg, msgList) {
		s << msg << endl;
	}
	return ret;
}

void Core::emitSendBytes(const QByteArray& buffer)
{
	emit sendBytes(buffer);
}

CoreBroadcast::CoreBroadcast(Core *core, QByteArray *buffer) :
        core(core),
        lock(&core->mutex),
        array(buffer),
        stream(array, QIODevice::WriteOnly)
{}

CoreBroadcast::~CoreBroadcast()
{
	core->emitSendBytes(*array);
}
}

ZInterprocess::Core *createCore(const FilePath& outFilePath, QString& resultMsg)
{
  return new ZInterprocess::Core(outFilePath, resultMsg);
}

void releaseCore(ZInterprocess::Core *core)
{
  if (NULL != core)
    delete core;
}
