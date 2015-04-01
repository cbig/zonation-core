#include "ini.h"
#include "qstring_spirit.h"
#include "spirit_util.h"

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <QHash>

#include <QString>
#include <QVector>
#include <QFile>

namespace {

void StoreIniKeyValue(IniMap& map, const QString& section, const QString& key, const QString& value, bool quoted)
{
	map.push_back(IniEntry(section, key, value, quoted));
}

void onError(ErrorCallback& err, const boost::spirit::qi::info& what, int row)
{
	std::stringstream ss;
	ss << what;
	ZError::err_msg(QObject::tr("Settings file: expecting %1 on row %2").arg(QString::fromUtf8(ss.str().c_str())).arg(row), zeNotice);
}

QString trimmed(QString const& s)
{
	return s.trimmed();
}
}

bool loadIniFile(IniMap& map, QIODevice& device, ErrorCallback& err)
{
	device.setTextModeEnabled(true);
	QTextStream in(&device);
	in.setCodec("UTF-8");
	in.setAutoDetectUnicode(true);

	typedef QVector<Char> Vector;
	typedef Vector::ConstIterator Iterator;

	QString section;
	QString key;
	bool quoted;
	int row = 0;

	namespace spirit = boost::spirit;
	namespace phoenix = boost::phoenix;
	namespace qi = boost::spirit::qi;
	namespace unicode = boost::spirit::unicode;
	namespace ascii = boost::spirit::ascii;

	QuotedStringParser<Iterator> sectionRule('[', ']');
	QuotedStringParser<Iterator> quotedStringRule('"');
	typedef CommentSkipper<Iterator> SkipType;
	SkipType skipRule(";#");

	qi::rule<Iterator> keyEndRule = spirit::lit('=');
	keyEndRule.name("=");
	UnquotedStringParser<Iterator> unquotedKeyRule("=;#", true, false);
	UnquotedStringParser<Iterator> unquotedValueRule(";#", true, true);

	qi::rule<Iterator, QString()> keyRule = quotedStringRule [qi::_val = qi::_1] |
	                unquotedKeyRule [qi::_val = phoenix::bind(&::trimmed, qi::_1)];
	qi::rule<Iterator, QString()> valueRule = quotedStringRule [qi::_val = qi::_1, phoenix::ref(quoted) = phoenix::val(true)] |
	                unquotedValueRule [qi::_val = phoenix::bind(&::trimmed, qi::_1), phoenix::ref(quoted) = phoenix::val(false)];
	qi::rule<Iterator, void(), SkipType> keyValueRule = keyRule [phoenix::ref(key) = qi::_1] > keyEndRule > valueRule [phoenix::bind(&StoreIniKeyValue, phoenix::ref(map), phoenix::ref(section), phoenix::ref(key), qi::_1, phoenix::ref(quoted))];

	qi::rule<Iterator> eoiRule = qi::eoi;
	eoiRule.name("end of line/input");
	qi::rule<Iterator, void(), SkipType> lineRule = -(sectionRule [phoenix::ref(section) = qi::_1] | keyValueRule) > eoiRule;

	qi::on_error<qi::accept>(lineRule, phoenix::bind(&onError, phoenix::ref(err), qi::_4, phoenix::ref(row)));

	while(!in.atEnd()) {
		Vector line(in.readLine().toUcs4());
		++row;

		qi::phrase_parse(line.constBegin(), line.constEnd(), lineRule, skipRule);
	}

	return true;
}
