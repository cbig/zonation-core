#ifndef ZINTERPROCESS_H
#define ZINTERPROCESS_H

#include "filepath.h"
#include "gisnodata.h"
#include <QThread>
#include <QLocalServer>
#include <QVector>
#include <QList>
#include <QSharedPointer>
#include <QFileInfo>
#include <QMutex>
#include <QMutexLocker>
#include <QBuffer>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/shared_ptr.hpp>
#include <QProcess>
#include <QMutex>

namespace ZInterprocess {

// data

struct InstanceInitData {
  quint64 siteCount;
  quint32 mapWidth;
  quint32 mapHeight;
  quint32 plotCount;
  quint32 plotSize;
  InstanceInitData(quint64 siteCount = 0,
		   quint32 mapWidth = 0,
		   quint32 mapHeight = 0,
		   quint32 plotCount = 0,
		   quint32 plotSize = 0) :
  siteCount(siteCount),
    mapWidth(mapWidth),
    mapHeight(mapHeight),
    plotCount(plotCount),
    plotSize(plotSize)
  { };
};

inline QDataStream& operator<<(QDataStream& stream, const InstanceInitData& initData)
{
	return stream
	                << initData.siteCount
	                << initData.mapWidth
	                << initData.mapHeight
	                << initData.plotCount
	                << initData.plotSize
	                   ;
}

inline QDataStream& operator>>(QDataStream& stream, InstanceInitData& initData)
{
	return stream
	                >> initData.siteCount
	                >> initData.mapWidth
	                >> initData.mapHeight
	                >> initData.plotCount
	                >> initData.plotSize
	                   ;
}

struct PlotPoint
{
	float x;
	float y;
	PlotPoint(float x = 0.0f, float y = 0.0f) : x(x), y(y) {}
};

inline QDataStream& operator<<(QDataStream& stream, const PlotPoint& point)
{
	return stream << point.x << point.y;
}

inline QDataStream& operator>>(QDataStream& stream, PlotPoint& point)
{
	return stream >> point.x >> point.y;
}

typedef QList<PlotPoint> Plot;

//typedef QVector<quint32> ColorMap;
//typedef QList<quint32> ColorMap;
typedef std::vector<quint32> ColorMap;

//typedef QVector<float> RankMap;
//typedef QList<float> RankMap;
typedef std::vector<float> RankMap;

typedef QVector<Plot> Plots;
//typedef std::vector<Plot> Plots;

struct InstanceData
{
	InstanceInitData init;
	quint64 removedCount;
	ColorMap colorMap;
	RankMap rankMap;
	Plots plots;

	InstanceData(InstanceInitData initData = InstanceInitData()) :
	        init(initData),
	        removedCount(0),
		//colorMap(init.mapHeight * init.mapWidth, nodata<quint32>()),
		  colorMap(1),
	        //rankMap(init.mapHeight * init.mapWidth, nodata<float>()),
		  rankMap(1),
	        plots(init.plotCount)
       { };

	void setPixel(quint32 x, quint32 y, quint32 color)
	{
	  //colorMap[x + y * init.mapWidth] = color;
	}

	// returns offset
	quint64 removeSite(quint32 x, quint32 y, quint32 color, float rank)
	{
		quint64 pos(x + y * init.mapWidth);
		//colorMap[pos] = color;
		//rankMap[pos] = rank;
		++removedCount;
		return pos;
	}

	quint32 addPlotPoint(quint32 plotIndex, float x, float y)
	{
		Plot& plot(plots[plotIndex]);
		plot.push_back(PlotPoint(x, y));
		return plot.size();
	}
};

inline QDataStream& operator<<(QDataStream& stream, const std::vector<quint32> v)
{
  return stream;  // TODO
}

inline QDataStream& operator>>(QDataStream& stream, const std::vector<quint32> v)
{
  return stream;  // TODO
}

inline QDataStream& operator<<(QDataStream& stream, const std::vector<float> v)
{
  return stream;  // TODO
}

inline QDataStream& operator>>(QDataStream& stream, const std::vector<float> v)
{
  return stream;  // TODO
}

inline QDataStream& operator<<(QDataStream& stream, const InstanceData& data)
{
  return stream << data.init << data.removedCount /*<< data.colorMap << data.rankMap*/ << data.plots;
}

inline QDataStream& operator>>(QDataStream& stream, InstanceData& data)
{
  return stream >> data.init >> data.removedCount /*>> data.colorMap >> data.rankMap*/ >> data.plots;
}

typedef QString Msg;
typedef QList<QString> MsgList;


// tags

typedef quint64 MessageSize;

typedef qint32 TagType;

enum Tag {
	CORE_CURRENT_DATA_NOTIFY = 0,
	CORE_CURRENT_MSGLIST_NOTIFY = 1,
	CORE_MSG_NOTIFY = 2,
	CORE_INITIALIZED_NOTIFY = 3,
	CORE_SITE_REMOVED_NOTIFY = 4,
	CORE_PLOT_POINT_ADDED_NOTIFY = 5,
	CORE_DONE_NOTIFY = 6,
	CORE_CLOSING_NOTIFY = 7,
};

typedef qint32 CoreStatusType;

enum CoreStatus {
	CORE_NOT_INITIALIZED,
	CORE_INITIALIZED,
	CORE_DONE
};

// serialization

// server name from existing filename
inline QString getServerName(const FilePath& out)
{
  qDebug() << "getServerName(): " << out << endl;
	{ // make sure the file exist
		QFile file(out.absoluteFilePath());
		file.open(QIODevice::WriteOnly);
	}
  qDebug() << "getServerName(), out.canonicalFilePath(): " << out.canonicalFilePath() << endl;
  // QFileInfo::canonicalFilePath with multiple threads may fail (I've seen it returning empty string)
  // It looks like it is not thread safe in some systems. Solution: mutex.
  QMutex mutex;
  mutex.lock();
	QString ret(out.canonicalFilePath());
  mutex.unlock();
  // still, try to do something if canonicalFilePath fails
  if (ret.isEmpty()) {
    ret = out;
    qDebug() << "getServerName(), canonicalFilePath returned empty string, using: '" << ret << "' instead." << endl;  
  }
	ret.replace("_", "__");
	ret.replace("/", "_s");
	ret.replace("\\", "_b");
	ret.replace(":", "_c");
	ret.replace(".", "_p");
	qDebug() << "getServerName(), done. Result: " << ret.toStdString().c_str() << endl;
	return ret;
}


/// more types...

struct MessageData
{
	quint32 msgCount;
	QString msg;
	MessageData() {}
	MessageData(quint32 msgCount, const QString& msg) : msgCount(msgCount), msg(msg) {}
};

inline QDataStream& operator<<(QDataStream& stream, const MessageData& data)
{
	return stream << data.msgCount << data.msg;
}

inline QDataStream& operator>>(QDataStream& stream, MessageData& data)
{
	return stream >> data.msgCount >> data.msg;
}

struct CoreInitializedData
{
	InstanceInitData init;
	ColorMap colorMap;
	CoreInitializedData() {}
	CoreInitializedData(const InstanceInitData& init, const ColorMap& colorMap) :
	        init(init), colorMap(colorMap) {}
};

inline QDataStream& operator<<(QDataStream& stream, const CoreInitializedData& data)
{
  return stream << data.init /*<< data.colorMap*/;
}

inline QDataStream& operator>>(QDataStream& stream, CoreInitializedData& data)
{
  return stream >> data.init /*>> data.colorMap*/;
}


struct SiteRemovedData
{
	quint64 removedCount;
	quint64 mapOffset;
	quint32 pixelColor;
	float rank;
	SiteRemovedData() {}
	SiteRemovedData(quint64 removedCount,
	                quint64 mapOffset,
	                quint32 pixelColor,
	                float rank) :
	        removedCount(removedCount),
	        mapOffset(mapOffset),
	        pixelColor(pixelColor),
	        rank(rank)
	{}
};

inline QDataStream& operator<<(QDataStream& stream, const SiteRemovedData& data)
{
	return stream << data.removedCount << data.mapOffset << data.pixelColor << data.rank;
}

inline QDataStream& operator>>(QDataStream& stream, SiteRemovedData& data)
{
	return stream >> data.removedCount >> data.mapOffset >> data.pixelColor >> data.rank;
}

struct PlotPointAddedData
{
	quint32 plotIndex;
	quint32 plotPointCount;
	PlotPoint point;
	PlotPointAddedData() :
	        plotIndex(-1),
	        plotPointCount(0),
	        point()
	{}
	PlotPointAddedData(quint32 plotIndex,
	                   quint32 plotPointCount,
	                   PlotPoint point) :
	        plotIndex(plotIndex),
	        plotPointCount(plotPointCount),
	        point(point)
	{}
};

inline QDataStream& operator<<(QDataStream& stream, const PlotPointAddedData& data)
{
	return stream << data.plotIndex << data.plotPointCount << data.point;
}

inline QDataStream& operator>>(QDataStream& stream, PlotPointAddedData& data)
{
	return stream >> data.plotIndex >> data.plotPointCount >> data.point;
}

// parsing

struct DefaultVisitor
{
	void operator()(const ZInterprocess::MessageData& data) {}
	void operator()(const ZInterprocess::SiteRemovedData& data) {}
	void operator()(const ZInterprocess::MsgList& data) {}
	void operator()(const ZInterprocess::InstanceData& data) {}
	void operator()(const ZInterprocess::CoreInitializedData& data) {}
	void operator()(const ZInterprocess::PlotPointAddedData& data) {}
};

template <class Visitor>
inline void parseData(TagType tag, const QByteArray& array, Visitor& visitor)
{
	QDataStream stream(array);
	switch(tag) {
	case CORE_CURRENT_DATA_NOTIFY: {
		InstanceData data;
		stream >> data;
		visitor(data);
	} break;
	case CORE_CURRENT_MSGLIST_NOTIFY: {
		MsgList data;
		stream >> data;
		visitor(data);
	} break;
	case CORE_MSG_NOTIFY: {
		MessageData data;
		stream >> data;
		visitor(data);
	} break;
	case CORE_INITIALIZED_NOTIFY: {
		CoreInitializedData data;
		stream >> data;
		visitor(data);
	} break;
	case CORE_SITE_REMOVED_NOTIFY: {
		SiteRemovedData data;
		stream >> data;
		visitor(data);
	} break;
	case CORE_PLOT_POINT_ADDED_NOTIFY: {
		PlotPointAddedData data;
		stream >> data;
		visitor(data);
	} break;
	case CORE_DONE_NOTIFY:
		break;
	case CORE_CLOSING_NOTIFY:
		break;
	}
}
}

#endif // ZINTERPROCESS_H
