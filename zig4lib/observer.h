#ifndef OBSERVER_H
#define OBSERVER_H

#include "zig4lib_global.h"
#include "zinterprocess.h"
#include "pod.h"
#include "io.h"
#include <QThread>
#include <QLocalSocket>
#include <QProcess>
#include <QTimer>
#include <QQueue>
#include <QWaitCondition>
#include <QMutex>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/hashed_index.hpp>

namespace ZInterprocess {

typedef int ZInstanceId;

// minimal class describing a zonation process
class ZIG4LIBSHARED_EXPORT ZInstanceProcessInfo
{
public:
	ZInstanceId id;
	CommandLine cmdLine;
	DirPath workingDir;
	QString serverName;
	ZInstanceProcessInfo();
	// to start a process
	ZInstanceProcessInfo(ZInstanceId id,
	                     const ZCommandLine& instance);
	// to start a process
	ZInstanceProcessInfo(ZInstanceId id,
	                     const ZCommandLine& instance,
	                     const DirPath& workingDir);
};
}

Q_DECLARE_METATYPE(ZInterprocess::ZInstanceId)
Q_DECLARE_METATYPE(ZInterprocess::ZInstanceProcessInfo)

namespace ZInterprocess {

using namespace boost;
using namespace boost::multi_index;

struct SocketTag {};
struct ProcessTag {};
struct TimerTag {};
struct IdTag {};

class SocketEntry
{
public:
	enum SocketStatus
	{
		Unconnected,
		Connected
	};

	enum ConnectionState
	{
		Waiting,
		Reading
	};

	shared_ptr<QLocalSocket> socket;
	QLocalSocket* socketID;
	shared_ptr<QProcess> process;
	QProcess* processID;
	shared_ptr<QTimer> timer;
	QTimer* timerID;
	ZInstanceProcessInfo info;
	bool tryObserve; // if true, run process only after trying to observe first
	mutable SocketStatus status;
	mutable bool connectedEmitted;
	mutable bool disconnectedEmitted;
	mutable bool failedEmitted;
	mutable ConnectionState connectionState;
	mutable MessageSize messageSize;
	mutable bool processStartedHere; // true if QProcess was started by loop
	mutable bool tryingObserve; // true if process will be started only if observe fails
	SocketEntry(const shared_ptr<QLocalSocket>& socket,
	            const shared_ptr<QProcess>& process,
	            const shared_ptr<QTimer>& timer,
	            const ZInstanceProcessInfo& info,
	            bool tryObserve = false) :
	        socket(socket),
		socketID(socket.get()),
	        process(process),
		processID(process.get()),
	        timer(timer),
		timerID(timer.get()),
	        info(info),
	        tryObserve(tryObserve),
	        status(Unconnected),
	        connectedEmitted(false),
	        disconnectedEmitted(false),
	        failedEmitted(false),
	        connectionState(Waiting),
	        messageSize(0),
	        processStartedHere(false),
	        tryingObserve(false)
	{}
	~SocketEntry()
	{
		qDebug() << "deleting entry";
		socketID = NULL;
		processID = NULL;
		timerID = NULL;
	}

	ZInterprocess::ZInstanceId id() const
	{
		return info.id;
	}
};

typedef multi_index_container<
	SocketEntry,
	indexed_by<
		hashed_unique<
			tag<SocketTag>,
                        member<SocketEntry, QLocalSocket*, &SocketEntry::socketID>
		>,
		hashed_unique<
			tag<ProcessTag>,
                        member<SocketEntry, QProcess*, &SocketEntry::processID>
		>,
		hashed_unique<
			tag<TimerTag>,
                        member<SocketEntry, QTimer*, &SocketEntry::timerID>
		>,
		hashed_unique<
			tag<IdTag>,
			const_mem_fun<SocketEntry, ZInterprocess::ZInstanceId, &SocketEntry::id>
		>
	>
> SocketMap;

class Observer;

class ObserverLoop : public QObject
{
	Q_OBJECT

public:
	ObserverLoop(Observer *observer, unsigned maxProcessCount, const QString& outFile = QString());

	SocketMap::iterator init(const ZInstanceProcessInfo& info, bool createProcess = false, bool tryObserve = false);
	void parseData(SocketMap::iterator it);

	SocketMap map;
	// Note this processQueue != god.h:Process::Queue processQueue !!!
	QQueue<shared_ptr<QProcess> > processQueue;
	unsigned processCount;
	unsigned maxProcessCount;

	QByteArray reusabledBuffer; // reused buffer for tiny messages

	QWaitCondition pauseCond;
	QMutex pauseMutex;
	bool paused;

	QString outFile;

	Observer *observer;

	bool queuePaused;

public slots:
	void observe(const ZInterprocess::ZInstanceProcessInfo& info);
	void runAndObserve(const ZInterprocess::ZInstanceProcessInfo& info);
	void queueAndRun(const ZInterprocess::ZInstanceProcessInfo& info, bool hold);
	void queueAndObserveOrRun(const ZInterprocess::ZInstanceProcessInfo& info, bool hold);

private slots:
	void readyRead();
	void close();
	void pause();
	void pauseQueue();
	void resumeQueue();
	void tryStart();
	void tryListen();
	void socketConnected();
	void socketDisconnected();
	void socketError(QLocalSocket::LocalSocketError socketError);
	void processError(QProcess::ProcessError error);
	void processFinished(int exitCode, QProcess::ExitStatus exitStatus);
	void processStarted();
	void quit();
	void remove(ZInterprocess::ZInstanceId id);
	void setMaxProcesses(int n);

signals:
	// either connected and disconnected will be called
	// or only errored if there was an error launching process/listening
	void connected(ZInterprocess::ZInstanceId id);
	void disconnected(ZInterprocess::ZInstanceId id);
	void failed(ZInterprocess::ZInstanceId id);
	// fired when processcount decreases to zero
	void queueFinished();
	// fired when received a bytearray with the given tag
	void received(ZInterprocess::ZInstanceId id, ZInterprocess::TagType tag, const QByteArray& array);

};

class Observer : public QThread
{
	Q_OBJECT

public:
	Observer(unsigned maxProcessCount, const QString& outFile = QString());
	~Observer();

	// API

	/// only observe instance
	ZIG4LIBSHARED_EXPORT void observe(const ZInstanceProcessInfo& instance);

	/// fire instance up immediately, regardless of the
	/// process count and observe
	ZIG4LIBSHARED_EXPORT void runAndObserve(const ZInstanceProcessInfo& instance);

	/// queue instance for run and observe
	ZIG4LIBSHARED_EXPORT void queueAndRun(const ZInstanceProcessInfo& instance, bool hold = false);

	/// queue instance and observe or run if instance is not present
	ZIG4LIBSHARED_EXPORT void queueAndObserveOrRun(const ZInstanceProcessInfo& instance, bool hold = false);

	ZIG4LIBSHARED_EXPORT void pauseQueue();
	ZIG4LIBSHARED_EXPORT void resumeQueue();

	ZIG4LIBSHARED_EXPORT void pauseLoop();
	ZIG4LIBSHARED_EXPORT void resumeLoop();

	/// instances can be fetched from this map after receiving
	/// a signal with the corresponding id
	//InstanceMap map;
	shared_ptr<ObserverLoop> loop;

	ZIG4LIBSHARED_EXPORT bool queuePaused();

	ZIG4LIBSHARED_EXPORT void remove(ZInterprocess::ZInstanceId id);
	ZIG4LIBSHARED_EXPORT void setMaxProcesses(int n);

signals:
	// internal
	void observeSignal(const ZInterprocess::ZInstanceProcessInfo& info);
	void runAndObserveSignal(const ZInterprocess::ZInstanceProcessInfo& info);
	void queueAndRunSignal(const ZInterprocess::ZInstanceProcessInfo& info, bool hold);
	void queueAndObserveOrRunSignal(const ZInterprocess::ZInstanceProcessInfo& info, bool hold);
	void closeSignal();
	void pauseSignal();
	void tryStartSignal();
	void pauseQueueSignal();
	void resumeQueueSignal();
	void removeSignal(ZInterprocess::ZInstanceId id);
	void setMaxProcessesSignal(int n);
};
}

ZIG4LIBSHARED_EXPORT ZInterprocess::Observer *createObserver(unsigned maxProcessCount, const QString& outFile = QString());
ZIG4LIBSHARED_EXPORT void releaseObserver(ZInterprocess::Observer *obs);

#endif // OBSERVER_H
