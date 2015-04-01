#include "observer.h"
#include <QDebug>
#include <QMetaType>
#include <QTimer>
#include <QByteArray>
#include <boost/make_shared.hpp>

namespace {
const int CONNECT_RETRY_INTERVAL_MSECS = 10;
const unsigned int REUSABLE_BUFFER_SIZE = 4096;

struct RegisterMetaTypes
{
	// gets called the first time Observer is created
	RegisterMetaTypes()
	{
		qRegisterMetaType<ZInterprocess::ZInstanceProcessInfo>("ZInterprocess::ZInstanceProcessInfo");
		qRegisterMetaType<CommandLine>("CommandLine");
		qRegisterMetaType<ZInterprocess::ZInstanceId>("ZInterprocess::ZInstanceId");
		qRegisterMetaType<ZInterprocess::TagType>("ZInterprocess::TagType");
	}
};

struct UpdateStatus
{
	UpdateStatus(ZInterprocess::SocketEntry::SocketStatus status) : status(status) {}
	void operator()(ZInterprocess::SocketEntry& entry) { entry.status = status; }
private:
	ZInterprocess::SocketEntry::SocketStatus status;
};

struct DeleteQObject
{
	void operator()(QObject *obj)
	{
		qDebug() << "deleting object";
		//obj->deleteLater();
		delete obj;
		qDebug() << "deleted object";
	}
};
}

namespace ZInterprocess {

using namespace boost;

ZInstanceProcessInfo::ZInstanceProcessInfo() : id(0) {}

// This constructor is used in zig4gui/processviewmodel.h (when putting projects into the queue)
ZInstanceProcessInfo::ZInstanceProcessInfo(ZInstanceId id,
                                           const ZCommandLine& instance) :
        id(id),
        cmdLine(getZCommandLineArguments(instance)),
        workingDir(instance.workingDir),
        serverName(getServerName(workingDir.absoluteFilePath(instance.outFile)))
{
  qDebug() << "ZInstanceProcessInfo::ZInstanceProcessInfo(2 args): id: " << id << ", name: " << serverName << ", working dir: " << workingDir;
}

// This constructor is used in zig4run/main
ZInstanceProcessInfo::ZInstanceProcessInfo(ZInstanceId id,
                                           const ZCommandLine& instance,
                                           const DirPath& workingDir) :
        id(id),
        cmdLine(getZCommandLineArguments(instance)),
        workingDir(workingDir),
        serverName(getServerName(workingDir.absoluteFilePath(instance.outFile)))
{
  qDebug() << "ZInstanceProcessInfo::ZInstanceProcessInfo(3 args): id: " << id << ", name:" << serverName << ", working dir: " << workingDir;
}


// ObserverLoop

ObserverLoop::ObserverLoop(Observer *observer, unsigned maxProcessCount, const QString& outFile) :
        processCount(0),
        maxProcessCount(maxProcessCount),
        reusabledBuffer(REUSABLE_BUFFER_SIZE, 0),
        paused(false),
        outFile(outFile),
        observer(observer),
        queuePaused(false)
{}

// This creates a new SocketEntry and inserts it in the map.
// info provides for example the workingDir
SocketMap::iterator ObserverLoop::init(const ZInstanceProcessInfo& info, bool createProcess, bool tryObserve)
{
	//shared_ptr<QLocalSocket> socket(make_shared<QLocalSocket>());
	//shared_ptr<QTimer> timer(make_shared<QTimer>());
	shared_ptr<QLocalSocket> socket(new QLocalSocket, DeleteQObject());
	shared_ptr<QTimer> timer(new QTimer, DeleteQObject());
	shared_ptr<QProcess> process;
	timer->setSingleShot(true);
	connect(socket.get(), SIGNAL(connected()), this, SLOT(socketConnected()));
	connect(socket.get(), SIGNAL(disconnected()), this, SLOT(socketDisconnected()));
	connect(socket.get(), SIGNAL(error(QLocalSocket::LocalSocketError)), this, SLOT(socketError(QLocalSocket::LocalSocketError)));
	connect(socket.get(), SIGNAL(readyRead()), this, SLOT(readyRead()));
	connect(timer.get(), SIGNAL(timeout()), this, SLOT(tryListen()));
	if(createProcess) {
		//process = make_shared<QProcess>();
		process = boost::shared_ptr<QProcess>(new QProcess, DeleteQObject());
		if(!outFile.isEmpty()) {
			qDebug() << "ObserverLoop::init(): redirecting output to" << outFile;
			process->setProcessChannelMode(QProcess::MergedChannels);
			process->setStandardOutputFile(outFile, QIODevice::Append);
		}
		qDebug() << "ObserverLoop::init(): setting process working dir: " << info.workingDir;
		// NOTE: info.workingDir here is the QProcess working dir
		//         --> path to the zig binary (myself?)
		// (as opposed to the workingDir that will be used later/elsewhere
		// when QueueAndRun() -> tryStart() is called from processviewmodel, 
                //  which sets the working dir to the path of the project (.bat) 
		process->setWorkingDirectory(info.workingDir);
		connect(process.get(), SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
		connect(process.get(), SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(processFinished(int,QProcess::ExitStatus)));
		connect(process.get(), SIGNAL(started()), this, SLOT(processStarted()));
	}

	qDebug() << process.get();
	SocketMap::iterator it(map.insert(SocketEntry(socket, process, timer, info, tryObserve)).first);
	qDebug() << "ObserverLoop::init(): map size:" << map.size();
	return it;
}

void ObserverLoop::tryStart()
{
	qDebug() << "ObserverLoop::tryStart()";
	while(!queuePaused && processCount < maxProcessCount && !processQueue.isEmpty()) {
		shared_ptr<QProcess> process(processQueue.dequeue());
		CommandLine cmdLine;

			typedef SocketMap::index<ProcessTag>::type Index;
			Index& index(map.get<ProcessTag>());
			Index::const_iterator it(index.find(process.get()));
		{
			cmdLine = it->info.cmdLine;
			if(it->tryObserve) {
			  qDebug() << "ObserverLoop::tryStart(): try observe, to server: " <<
			    it->info.serverName << endl; 
				it->tryingObserve = true;
				it->socket->connectToServer(it->info.serverName);
				// fedemp 2011/12/28: increased the timeout from 0
				if(it->socket->waitForConnected(100)) {
					qDebug() << "ObserverLoop::tryStart(): already running process detected";
					continue; // next candidate
				}
			}
			qDebug() << "ObserverLoop::tryStart(): starting process";
			it->processStartedHere = true;
			++processCount;
		}
		qDebug() << "ObserverLoop::tryStart(): right before process->start(...), with program: " <<
		  cmdLine.program << " and " << cmdLine.arguments.size() << " arguments, " << 
		  "workingDir: " << it->info.workingDir << ", serverName: " << it->info.serverName << endl;
		process->start(cmdLine.program, cmdLine.arguments);
	}
	qDebug() << "ObserverLoop::tryStart(): done";
}

// public slots

void ObserverLoop::observe(const ZInstanceProcessInfo& info)
{
	SocketMap::iterator it(init(info));
	it->timer->start(0);
}

void ObserverLoop::runAndObserve(const ZInstanceProcessInfo& info)
{
	shared_ptr<QProcess> process;
	{
		SocketMap::iterator it(init(info, true));
		it->processStartedHere = true;
		++processCount;
		process = it->process;
	}
	process->start(info.cmdLine.program, info.cmdLine.arguments);
}

// queued starters

void ObserverLoop::queueAndRun(const ZInstanceProcessInfo& info, bool hold)
{
	qDebug() << "ObserverLoop::queueAndRun()";
	SocketMap::iterator it(init(info, true));
	processQueue.enqueue(it->process);
	if(!hold) {
		tryStart();
	}

}

void ObserverLoop::queueAndObserveOrRun(const ZInterprocess::ZInstanceProcessInfo& info, bool hold)
{
	SocketMap::iterator it(init(info, true, true));
	processQueue.enqueue(it->process);
	if(!hold) {
		tryStart();
	}

}

// private slots

void ObserverLoop::readyRead()
{
	QLocalSocket *socket(static_cast<QLocalSocket *>(sender()));
	SocketMap::iterator it(map.find(socket));
	while(quint64 bytes = socket->bytesAvailable()) {
		switch(it->connectionState) {
		case SocketEntry::Waiting:
			if(bytes >= sizeof(MessageSize)) {
				QDataStream stream(socket);
				MessageSize size;
				stream >> size;
				it->messageSize = size;
				it->connectionState = SocketEntry::Reading;
				break;
			}
		case SocketEntry::Reading:
			if(bytes >= it->messageSize) {
				QDataStream stream(socket);
				TagType tag;
				stream >> tag;
				if(it->messageSize > REUSABLE_BUFFER_SIZE) {
					QByteArray array(socket->read(it->messageSize - sizeof(TagType)));
					emit received(it->info.id, tag, array);
				} else {
					socket->read(reusabledBuffer.data(), it->messageSize - sizeof(TagType));
					emit received(it->info.id, tag, reusabledBuffer);
				}
				it->messageSize = 0;
				it->connectionState = SocketEntry::Waiting;
				if(tag == CORE_CLOSING_NOTIFY) {
					socket->disconnectFromServer();
				}
				break;
			}
			return; // both ifs fail (very structured indeed...)
		}
	}
}

void ObserverLoop::close()
{
	qDebug() << "ObserverLoop::close(): disconnecting";
	for(SocketMap::const_iterator i = map.begin(); i != map.end(); ++i) {
		qDebug() << "ObserverLoop::close(): disconnect socket";
		i->socket->disconnect();
		// remaining processes will not be terminated!
		if(i->process) {
			qDebug() << "ObserverLoop::close(): disconnect process signals";
			i->process->disconnect();
		}
	}
	qDebug() << "ObserverLoop::close(): disconnect loop signals";
	disconnect();
	processQueue.empty();
	map.clear();
	qDebug() << "ObserverLoop::close(): signal quit";
	QTimer::singleShot(0, this, SLOT(quit())); // calls deletes first
}

void ObserverLoop::quit()
{
	qDebug() << "ObserverLoop::quit(): observer->quit()";
	observer->quit();
}

void ObserverLoop::pauseQueue()
{
	queuePaused = true;
}

void ObserverLoop::resumeQueue()
{
	queuePaused = false;
	tryStart();
}

void ObserverLoop::pause()
{
	QMutexLocker lock(&pauseMutex);
	while(paused) {
		qDebug() << "ObserverLoop::pause(): paused";
		pauseCond.wait(&pauseMutex);
	}
	qDebug() << "ObserverLoop::pause(): resumed";
}

void ObserverLoop::remove(ZInterprocess::ZInstanceId id)
{
	qDebug() << "ObserverLoop::remove():" << id;
	typedef SocketMap::index<IdTag>::type Index;
	Index& index(map.get<IdTag>());
	Index::iterator it(index.find(id));
	if(it == index.end()) {
		return;
	}
	//check if process is still in the queue
	QQueue<shared_ptr<QProcess> >::iterator ip = qFind(processQueue.begin(), processQueue.end(), it->process);
	bool finished = false;
	if(ip != processQueue.end()) {
		processQueue.erase(ip);
	}

	else if(it->processStartedHere) {
		if(--processCount == 0 && processQueue.empty()) {
			finished = true;
		}
	}

	//it->process->kill();

	// when we forcefully delete process and socket
	// signals will not be emitted.
	// We could emit the appropriate signals but
	// suits our modle better (processqueue in processviewmodel
	// has already deleted the corresponding instance)
	/*
	if(!it->disconnectedEmitted) {
		if(it->connectedEmitted) {
			emit disconnected(id);
		} else {
			emit failed(id);
		}
	}
	*/
	// forcefully terminate process and socket
	// this might need a more graceful approach if there are
	// any problems
	index.erase(it);
	if(finished) {
		emit queueFinished();
	}
	tryStart();
	qDebug() << "ObserverLoop::remove(ZInterprocess::ZInstanceId id): done";
}

void ObserverLoop::setMaxProcesses(int n)
{
	qDebug() << "ObserverLoop::setMaxProcesses():" << n;
	maxProcessCount = n;
	tryStart();
}

void ObserverLoop::tryListen()
{
        qDebug() << "ObserverLoop::tryListen(), map size: " << map.size();
	typedef SocketMap::index<TimerTag>::type Index;
	Index& index(map.get<TimerTag>());
	QTimer *timer(static_cast<QTimer *>(sender()));
	Index::iterator it(index.find(timer));
	if(it == index.end()) {
	  qDebug() << "ObserverLoop::tryListen(): process was already destroyed";
	} else {
	  qDebug() << "ObserverLoop::tryListen(): right before ->connectToServer, with serverName:" << 
	    it->info.serverName << endl; 
	  it->socket->connectToServer(it->info.serverName);
	}
	qDebug() << "ObserverLoop::tryListen(): done";
}

// socket slots

void ObserverLoop::socketConnected()
{
	qDebug() << "ObserverLoop::socketConnected()";
	QLocalSocket *socket(static_cast<QLocalSocket *>(sender()));
	SocketMap::iterator it(map.find(socket));
	it->status = SocketEntry::Connected;
	if(!it->processStartedHere) {
		++processCount;
		qDebug() << "ObserverLoop::socketConnected(): emit connected()";
		it->connectedEmitted = true;
		emit connected(it->info.id);
	}
	qDebug() << "ObserverLoop::socketConnected(): done";
}

// DIFFERENT BEHAVIOUR ON WIN (buggy) AND NIX:
// it seems that on windows disconnected might be
// called even though socket has not been connected

void ObserverLoop::socketDisconnected()
{
	qDebug() << "ObserverLoop::socketDisconnected()";
	QLocalSocket *socket(static_cast<QLocalSocket *>(sender()));
	SocketMap::iterator it(map.find(socket));
	switch(it->status) {
	case SocketEntry::Unconnected:
		//it->timer->start(CONNECT_RETRY_INTERVAL_MSECS);
		break;
	case SocketEntry::Connected:
		if(!it->processStartedHere) {
			qDebug() << "ObserverLoop::socketDisconnected(): emit disconnected()";
			it->disconnectedEmitted = true;
			emit disconnected(it->info.id);
			map.erase(it);
			if(--processCount == 0 && processQueue.empty()) {
				qDebug() << "ObserverLoop::socketDisconnected(): emit queueFinished()";
				emit queueFinished();
			}
			tryStart();
		}
		break;
	}
	qDebug() << "ObserverLoop::socketDisconnected(): done";
}

void ObserverLoop::socketError(QLocalSocket::LocalSocketError socketError)
{
	qDebug() << "ObserverLoop::socketError(): error" << socketError;
	QLocalSocket *socket(static_cast<QLocalSocket *>(sender()));
	SocketMap::iterator it(map.find(socket));
	qDebug() << "ObserverLoop::socketError(): status: " << it->status;
	switch(it->status) {
	case SocketEntry::Unconnected:
		it->timer->start(CONNECT_RETRY_INTERVAL_MSECS);
		if(!it->processStartedHere && !it->tryingObserve) {
			qDebug() << "ObserverLoop::socketError(): emit failed()";
			it->failedEmitted = true;
			emit failed(it->info.id);
		}
		break;
	case SocketEntry::Connected:
		/*
   if(!it->processStartedHere) {
    qDebug() << "disconnecting (disconnected)";
    it->discounted = true;
    emit disconnected(it->info.id);
    map.erase(it);
    if(--processCount == 0 && processQueue.empty())
     emit queueFinished();
    tryStart();
   }
   */
		break;
	}
	qDebug() << "ObserverLoop::socketError(): done";
}

// process slots

void ObserverLoop::processStarted()
{
        qDebug() << "ObserverLoop::processStarted(), map size: " << map.size();
	typedef SocketMap::index<ProcessTag>::type Index;
	Index& index(map.get<ProcessTag>());
	QProcess *process(qobject_cast<QProcess *>(sender()));
	Q_ASSERT(process);
	Index::iterator it(index.find(process));
	Q_ASSERT(it != index.end());
	if(it->processStartedHere) { // always true
		qDebug() << "ObserverLoop::processStarted(): emit connected()";
		it->connectedEmitted = true;
		emit connected(it->info.id);
	}
	it->timer->start(0);
	qDebug() << "ObserverLoop::processStarted(): done";
}

void ObserverLoop::processFinished(int exitCode, QProcess::ExitStatus exitStatus)
{
        qDebug() << "processFinished(): exitCode: " << exitCode << 
	  ", exitStatus: " << exitStatus;
	typedef SocketMap::index<ProcessTag>::type Index;
	Index& index(map.get<ProcessTag>());
	QProcess *process(static_cast<QProcess *>(sender()));
	Index::iterator it(index.find(process));
	if(it == index.end()) {
		return;
	}

	if(it->processStartedHere) { // always true
		qDebug() << "ObserverLoop::processFinished(): emit disconnected()";
		it->disconnectedEmitted = true;
		emit disconnected(it->info.id);
		index.erase(it);
		if(--processCount == 0 && processQueue.empty()) {
			qDebug() << "ObserverLoop::processFinished(): emit queueFinished()";
			emit queueFinished();
		}
		tryStart();
	}
	qDebug() << "processFinished(): done";
}

void ObserverLoop::processError(QProcess::ProcessError error)
{
	qDebug() << "ObserverLoop::processError():" << error;
	typedef SocketMap::index<ProcessTag>::type Index;
	Index& index(map.get<ProcessTag>());
	QProcess *process(static_cast<QProcess *>(sender()));
	Index::iterator it(index.find(process));
	if(it == index.end()) {
		return;
	}

	if(it->processStartedHere) { // always true
		qDebug() << "ObserverLoop::processError(): emit failed()";
		it->failedEmitted = true;
		emit failed(it->info.id);
	}

	if(error == QProcess::FailedToStart) {
		// processStarted was never called and
		// thus processFinished will not be called
		index.erase(it);
		if(--processCount == 0 && processQueue.empty()) {
			qDebug() << "ObserverLoop::processError(): emit queueFinished()";
			emit queueFinished();
		}
	}
	qDebug() << "ObserverLoop::processError(): done";
}

// Observer

/*
 InstanceEntry::InstanceEntry() : id(0) {}
 InstanceEntry::InstanceEntry(ZInstanceId id, const shared_ptr<ZInstance>& instance) :
   id(id), instance(instance)
 {}
 */

Observer::Observer(unsigned maxProcessCount, const QString& outFile)
{
	static RegisterMetaTypes reg;
	loop = make_shared<ObserverLoop>(this, maxProcessCount, outFile);
	loop->moveToThread(this);
	connect(this, SIGNAL(observeSignal(const ZInterprocess::ZInstanceProcessInfo&)),
	        loop.get(), SLOT(observe(const ZInterprocess::ZInstanceProcessInfo&)), Qt::QueuedConnection);
	connect(this, SIGNAL(runAndObserveSignal(const ZInterprocess::ZInstanceProcessInfo&)),
	        loop.get(), SLOT(runAndObserve(const ZInterprocess::ZInstanceProcessInfo&)), Qt::QueuedConnection);
	connect(this, SIGNAL(queueAndRunSignal(const ZInterprocess::ZInstanceProcessInfo&, bool)),
	        loop.get(), SLOT(queueAndRun(const ZInterprocess::ZInstanceProcessInfo&, bool)), Qt::QueuedConnection);
	connect(this, SIGNAL(queueAndObserveOrRunSignal(const ZInterprocess::ZInstanceProcessInfo&, bool)),
	        loop.get(), SLOT(queueAndObserveOrRun(const ZInterprocess::ZInstanceProcessInfo&, bool)), Qt::QueuedConnection);
	connect(this, SIGNAL(closeSignal()), loop.get(), SLOT(close()), Qt::QueuedConnection);
	connect(this, SIGNAL(pauseSignal()), loop.get(), SLOT(pause()), Qt::QueuedConnection);
	connect(this, SIGNAL(tryStartSignal()), loop.get(), SLOT(tryStart()), Qt::QueuedConnection);
	connect(this, SIGNAL(pauseQueueSignal()), loop.get(), SLOT(pauseQueue()), Qt::QueuedConnection);
	connect(this, SIGNAL(resumeQueueSignal()), loop.get(), SLOT(resumeQueue()), Qt::QueuedConnection);
	connect(this, SIGNAL(removeSignal(ZInterprocess::ZInstanceId)),
	        loop.get(), SLOT(remove(ZInterprocess::ZInstanceId)), Qt::QueuedConnection);
	connect(this, SIGNAL(setMaxProcessesSignal(int)), loop.get(), SLOT(setMaxProcesses(int)), Qt::QueuedConnection);
	start();
}

Observer::~Observer()
{
	qDebug() << "~Observer(): emit closeSignal()";
	emit closeSignal();
	qDebug() << "~Observer(): wait()";
	wait();
	qDebug() << "~Observer(): done";
}

/// only observe instance
void Observer::observe(const ZInstanceProcessInfo& instance)
{
	emit observeSignal(instance);
}

/// fire instance up immediately, regardless of the
/// process count and observe
void Observer::runAndObserve(const ZInstanceProcessInfo& instance)
{
	emit runAndObserveSignal(instance);
}

/// queue instance for run and observe
void Observer::queueAndRun(const ZInstanceProcessInfo& instance, bool hold)
{
	qDebug() << "Observer::queueAndRun()";
	emit queueAndRunSignal(instance, hold);
}

/// queue instance and observe or run if instance is not present
void Observer::queueAndObserveOrRun(const ZInstanceProcessInfo& instance, bool hold)
{
	emit queueAndObserveOrRunSignal(instance, hold);
}

void Observer::pauseLoop()
{
	QMutexLocker lock(&loop->pauseMutex);
	loop->paused = true;
	emit pauseSignal();
}

void Observer::resumeLoop()
{
	QMutexLocker lock(&loop->pauseMutex);
	loop->paused = false;
	loop->pauseCond.wakeAll();
}

void Observer::pauseQueue()
{
	emit pauseQueueSignal();
}

void Observer::resumeQueue()
{
	emit resumeQueueSignal();
}

bool Observer::queuePaused()
{
	// this is not necessarily atomic...
	return loop->queuePaused;
}

void Observer::remove(ZInterprocess::ZInstanceId id)
{
	emit removeSignal(id);
}

void Observer::setMaxProcesses(int n)
{
	emit setMaxProcessesSignal(n);
}

}

ZInterprocess::Observer *createObserver(unsigned maxProcessCount, const QString& outFile)
{
	return new ZInterprocess::Observer(maxProcessCount, outFile);
}

void releaseObserver(ZInterprocess::Observer *obs)
{
	delete obs;
}
