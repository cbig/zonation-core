#ifndef FILEPATH_H
#define FILEPATH_H

#include <QString>
#include <QFileInfo>
#include <QDir>
#include <QHash>
#include <QMetaType>
#include "exception.h"

// extension given wit a dot, eg ".jpg"
inline QString changeFileExtension(QString const& path, QString const& extension = "")
{
	int dir = path.lastIndexOf('/');
	int delim = path.lastIndexOf('.');
	if(delim > dir) { // file has a delimiter
		return path.left(delim) + extension;
	}
	return path + extension;
}

class FilePath : public QString
{
public:
	FilePath(const char *path)
	{
		QFileInfo info(path);
		absolute = info.isAbsolute();
		QString::operator=(info.filePath());
	}
	FilePath(const QString& path = QString())
	{
		QFileInfo info(path);
		absolute = info.isAbsolute();
		QString::operator=(info.filePath());
	}
	FilePath(const QFileInfo& info)
	{
		absolute = info.isAbsolute();
		QString::operator=(info.filePath());
	}

	FilePath absoluteFilePath() const
	{
		return QFileInfo(*this).absoluteFilePath();
	}

	FilePath canonicalFilePath() const
	{
		return QFileInfo(*this).canonicalFilePath();
	}

	bool isAbsolute() const
	{
		return absolute;
	}

	QString fileName() const
	{
		return QFileInfo(*this).fileName();
	}

	// fedemp: more duplication of QFileInfo...
	QString filePath() const
	{
		return QFileInfo(*this).filePath();
	}

	QString nativePathName() const
	{
		return QDir::toNativeSeparators(*this);
	}

	bool exists() const
	{
		return QFileInfo(*this).exists();
	}

	FilePath changeExtension(QString const& extension = "") const
	{
		return changeFileExtension(*this, extension);
	}

	FilePath removeExtension() const
	{
		return changeExtension();
	}

private:
	bool absolute;
};

Q_DECLARE_METATYPE(FilePath)


// absolute
class DirPath : public QString
{
public:
	DirPath(const char *path) :
	        QString(QDir(path).absolutePath())
	{}
	DirPath(const QString& path = QString()) :
	        QString(QDir(path).absolutePath())
	{}

	// parent directory
	DirPath(const FilePath& path) :
	        QString(QFileInfo(path).dir().absolutePath())
	{}

	FilePath absoluteFilePath(const FilePath& path) const
	{
		return FilePath(QDir(*this).absoluteFilePath(path));
	}

	// canonicalize

	FilePath canonicalFilePath(const FilePath& path) const
	{
		FilePath filepath(absoluteFilePath(path));
		QFileInfo info(QFileInfo(filepath).canonicalFilePath());
		if(info.isFile()) {
			return FilePath(info);
		}
		throw Exception("cannot canonicalize file: " + path);
	}

	DirPath canonicalPath() const
	{
		QString path(QDir(*this).canonicalPath());
		if(path.isEmpty()) {
			throw Exception("cannot canonicalize directory: " + *this);
		}
		return DirPath(path);
	}

	QString nativePathName() const
	{
		return QDir::toNativeSeparators(*this);
	}
};

Q_DECLARE_METATYPE(DirPath)

#endif // FILEPATH_H
