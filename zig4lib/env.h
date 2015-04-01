#ifndef ENV_H
#define ENV_H

#include <QStringList>
#include <QString>
#include <QDataStream>

/**
 * Defines stuff needed for the "argument forwarding mode" - a special mode
 * whereby Z does not actually run, but connects to the server specified in
 * the 'ZIG_ARG_SERVER_NAME' env. variable, and sends there its command line/
 * argv.
 * 
 * The idea was that, to support complicated .bat "project" files,
 * these are send to cmd.exe rather than interpreted by Z GUI(1). cmd.exe
 * is run with the env. variable ZIG_FORWARD_ARGS defined.  Then, Z in
 * this special mode (enabled when the 'ZIG_FORWARD_ARGS'
 * env. variable is defined) takes its own argv (which has been
 * properly defined by the native .bat interpreter and sends it to the
 * arg. server.
 *
 * The ArgServerLoop class in io_arg.h/cpp implements most of this.
 * More stuff in io.cpp.
 *
 * (1) traditional parsing happens in zig4/core/loaddata.cpp:Parse_cmd_line
 *               -> zig4lib/io.cpp:parseZCommandLine
 *
 * (fedemp: this is my fuzzy interpretation...)
 */

namespace ZEnvironment {
static int const timeoutMsecs = 5000;
static int const totalTimeoutMsecs = 20000;
static QString const forwardArgsVariable = QString("ZIG_FORWARD_ARGS");
static QString const argServerNameVariable = QString("ZIG_ARG_SERVER_NAME");
}

class ArgModeData
{
public:
	QString appFilePath;
	QString workingDir;
	QStringList arguments;
	ArgModeData() {}
	ArgModeData(QString const& appFilePath, QString const& workingDir, QStringList const& arguments) :
	        appFilePath(appFilePath), workingDir(workingDir), arguments(arguments)
	{
	}
};

inline QDataStream& operator<<(QDataStream& stream, const ArgModeData& data)
{
	return stream << data.appFilePath << data.workingDir << data.arguments;
}

inline QDataStream& operator>>(QDataStream& stream, ArgModeData& data)
{
	return stream >> data.appFilePath >> data.workingDir >> data.arguments;
}

#endif // ENV_H
