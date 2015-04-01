#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <exception>
#include <string>
#include <QString>
#include <QByteArray>

struct Exception : public std::exception
{
	// in utf-8
	QByteArray err;
	Exception(const char *err) : err(err) {}
	Exception(const std::string& err) : err(err.c_str()) {}
	Exception(const QString& err) : err(err.toUtf8()) {}
	virtual ~Exception() throw() {}
	virtual const char *what() const throw()
	{
		return err.constData();
	}
	virtual const QString qwhat() const throw()
	{
		return QString::fromUtf8(err.constData());
	}
};

#endif // EXCEPTION_H
