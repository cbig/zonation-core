#ifndef INI_H
#define INI_H

#include <QString>
#include <QIODevice>
#include <QTextStream>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include "string_util.h"
#include "zig4lib_global.h"
#include "error_util.h"

struct IniMapBySectionTag {};
struct IniMapByKeyTag {};

class IniEntry
{
public:
	QString section;
	QString key;
	QString value;
	bool quoted;
	IniEntry(QString section, QString key, QString value, bool quoted = true) :
	        section(section), key(key), value(value), quoted(quoted)
	{}
};

typedef boost::multi_index_container<
	IniEntry,
	boost::multi_index::indexed_by<
		boost::multi_index::sequenced<>,
		boost::multi_index::hashed_non_unique<
			boost::multi_index::tag<IniMapBySectionTag>,
			BOOST_MULTI_INDEX_MEMBER(IniEntry, QString, section)
		>,
		boost::multi_index::hashed_non_unique<
			boost::multi_index::tag<IniMapByKeyTag>,
			boost::multi_index::composite_key<
				IniEntry,
				BOOST_MULTI_INDEX_MEMBER(IniEntry, QString, section),
				BOOST_MULTI_INDEX_MEMBER(IniEntry, QString, key)
			>
		>
	>
> IniMap;

typedef boost::multi_index::index<IniMap, IniMapBySectionTag>::type IniMapBySection;
typedef boost::multi_index::index<IniMap, IniMapByKeyTag>::type IniMapByKey;

bool ZIG4LIBSHARED_EXPORT loadIniFile(IniMap& map, QIODevice& device, ErrorCallback&);
bool ZIG4LIBSHARED_EXPORT saveIniFile(const IniMap& map, QIODevice& device, ErrorCallback& err);

#endif // INI_H
