#ifndef ERROR_UTIL_H
#define ERROR_UTIL_H

#include <QString>
#include <boost/function.hpp>
#include <QDebug>

// zeLevel is needed by ErrorCallback and
// LineModify. ErrorCallback is needed by ZError (static methods...)
// Note that in preferences there are 3 levels: notice, warning, error.
//                 notice implies ok
enum zeLevel { zeUndef=0, zeHeader=1, zeNotice=10, zeOk=11, zeWarning=12, zeError=13, zeCritical=14};

// template functions don't accept default parameter values!
// typedef boost::function<void (const QString& error, ZError::Level)> ErrorCallback;

class ErrorCallback
{
 public:
  virtual void operator()(const QString& str, zeLevel lvl = zeWarning);

  virtual ~ErrorCallback()
    {};

  // private:
  // horror thing...
  static ErrorCallback* config_callback;  
};

class DefaultErrorCallback: public ErrorCallback
{
public:
  virtual void operator()(const QString& str, zeLevel lvl = zeWarning);
};

class LineModify: public ErrorCallback
{
public:
 LineModify(ErrorCallback e, const QString& prepend = "", const QString& append = "") :
  err(e),
    prepend(prepend),
    append(append)
    {
    }
  
  virtual void operator()(const QString& str, zeLevel lvl = zeWarning)
  {
    err(prepend + str + append, zeUndef);
  }
  
  ErrorCallback err;
  QString prepend;
  QString append;
};

class SuppressErrors: public ErrorCallback
{
public:
  virtual void operator()(const QString&, zeLevel lvl = zeWarning)
  {
  }
};

class ZError
{
 public:
  static void setErrorCallback(ErrorCallback& ecb);

  static void setErrorReportingLevel(zeLevel lvl);
  
  static void startProject(QString pname);
  
  static void startInstance(QString iname, QString line);
  
  static void err_msg(QString msg, zeLevel lvl = zeUndef);
  
  static QString errPrefix(zeLevel lev)
  { return ""; }//zelevel_prefix[lev]; };
  
 private:
  static QString currentProject;
  static QString currentInstance;
  static QString iline;
  static int proj_count;   // errors reported for current project
  static int inst_count;   // errors reported for current instance
  static QString level_prefix[];

  static zeLevel reporting_lvl;
  static ErrorCallback* err;
};

#endif // ERROR_UTIL_H
