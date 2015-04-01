#ifndef IO_ARG_H
#define IO_ARG_H

#include "env.h"
#include "error_util.h"
#include <QList>
#include <QMap>
#include <QByteArray>
#include <QLocalServer>
#include <QLocalSocket>
#include <QThread>

class ArgServerLoop : public QThread
{
	Q_OBJECT

public:
	ArgServerLoop(QString const& serverName, ErrorCallback& err);
	~ArgServerLoop();

	void close();

	QList<ArgModeData> args;
	bool success;

protected:
	void run();

private slots:
	void on_server_newConnection();
	void on_socket_disconnected();
	void on_socket_error(QLocalSocket::LocalSocketError socketError);
	//void on_socket_readChannelFinished();
	void on_socket_readyRead();

private:
	QString serverName;
	ErrorCallback& err;
	QMap<QLocalSocket *, QByteArray> buffer;
};

#endif // IO_ARG_H
