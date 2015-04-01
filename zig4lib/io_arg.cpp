#include "io_arg.h"

ArgServerLoop::ArgServerLoop(QString const& serverName, ErrorCallback& err) :
        success(true),
        serverName(serverName),
        err(err)
{
	moveToThread(this);
	start();
}

void ArgServerLoop::run()
{
	QLocalServer server;
	connect(&server, SIGNAL(newConnection()), this, SLOT(on_server_newConnection()));
	server.listen(serverName);
	exec();
}

ArgServerLoop::~ArgServerLoop()
{
	close();
}

void ArgServerLoop::close()
{
	exit();
	if(!wait(1000)) {
		err(tr("Error exiting thread, terminating"));
		setTerminationEnabled(true);
		terminate();
		if(!wait(1000)) {
			err(tr("Error terminating thread, destroying"));
		}
	}
}

void ArgServerLoop::on_server_newConnection()
{
	qDebug() << "ArgServerLoop::on_server_newConnection()";
	QLocalServer *server = qobject_cast<QLocalServer *>(sender());
	Q_ASSERT(server);
	QLocalSocket *socket(server->nextPendingConnection());
	Q_ASSERT(socket);
	// insert empty bytearray
	buffer[socket];
	connect(socket, SIGNAL(disconnected()), this, SLOT(on_socket_disconnected()));
	connect(socket, SIGNAL(error(QLocalSocket::LocalSocketError)), this, SLOT(on_socket_error(QLocalSocket::LocalSocketError)));
	connect(socket, SIGNAL(readyRead()), this, SLOT(on_socket_readyRead()));
}

void ArgServerLoop::on_socket_disconnected()
{
	qDebug() << "ArgServerLoop::on_socket_disconnected()";
	QLocalSocket *socket = qobject_cast<QLocalSocket *>(sender());
	Q_ASSERT(socket);
	QMap<QLocalSocket *, QByteArray>::iterator i = buffer.find(socket);
	Q_ASSERT(i != buffer.end());
	// parse buffered data from application
	//i.value().append(socket->readAll()); // final read
	{
		QDataStream stream(&i.value(), QIODevice::ReadOnly);
		ArgModeData arg;
		stream >> arg;
		if(stream.status() == QDataStream::Ok) {
			args << arg;
		} else {
			success = false;
			err(tr("Error parsing instance arguments"));
		}
	}
	buffer.erase(i);
	socket->deleteLater();

	// QLocalserver emits an error from completeAsyncRead() after disconnecting even
	// though no apparent error is actually happening (all data is being read correctly)
	// As a workaround we disable error signaling
	// This might have something to do with (https://bugreports.qt.nokia.com/browse/QTBUG-13646)

	disconnect(socket, SIGNAL(error(QLocalSocket::LocalSocketError)), this, SLOT(on_socket_error(QLocalSocket::LocalSocketError)));
}

void ArgServerLoop::on_socket_error(QLocalSocket::LocalSocketError socketError)
{
	qDebug() << "ArgServerLoop::on_socket_error():" << socketError;
	success = false;
	err(qobject_cast<QLocalSocket *>(sender())->errorString());
}

void ArgServerLoop::on_socket_readyRead()
{
	qDebug() << "ArgServerLoop::on_socket_readyRead()";
	QLocalSocket *socket = qobject_cast<QLocalSocket *>(sender());
	Q_ASSERT(socket);
	buffer[socket].append(socket->readAll());
}

/*
void ArgServerLoop::on_socket_readChannelFinished()
{
    qDebug() << "ArgServerLoop::on_socket_readChannelFinished()";
    QLocalSocket *socket = qobject_cast<QLocalSocket *>(sender());
    Q_ASSERT(socket);
    buffer[socket].append(socket->readAll());
}
*/
