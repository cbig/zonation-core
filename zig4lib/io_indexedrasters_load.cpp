#include "io.h"
#include "pod.h"
#include "spirit_util.h"
#include "spirit_error.h"
#include "ptr_util.h"

#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(
                ZIndexedRastersEntry,
                (int, index)
                (FilePath, raster)
                )

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/home/phoenix/container.hpp>
#include <boost/spirit/home/support/attributes.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>

bool loadZIndexedRastersFile(ZIndexedRastersFile& indexedRastersFile, QIODevice& device, ErrorCallback& err)
{
	device.setTextModeEnabled(true);
	QTextStream in(&device);
	in.setCodec("UTF-8");
	in.setAutoDetectUnicode(true);

	typedef QVector<Char> Vector;
	typedef Vector::ConstIterator Iterator;

	namespace spirit = boost::spirit;
	namespace ph = boost::phoenix;
	namespace qi = boost::spirit::qi;
	namespace unicode = boost::spirit::unicode;
	using namespace boost;

	int row = 0;

	typedef CommentSkipper<Iterator> SkipType;
	SkipType skipRule("#");

	// rules
	IntParser<Iterator> int_;

	StringParser<Iterator> rasterRule('"', "#", true);

	qi::rule<Iterator, ZIndexedRastersEntry(), SkipType> indexedRastersRule =
	                int_ > rasterRule;

	qi::rule<Iterator> indexedRastersEndRule = qi::eoi;
	indexedRastersEndRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule =
	                (indexedRastersRule [ph::push_back(&ph::ref(indexedRastersFile)->*&ZIndexedRastersFile::rasterList, ph::bind(&copy_and_make_shared<ZIndexedRastersEntry>, qi::_1))] >
	                 indexedRastersEndRule) | qi::eoi;
	lineRule.name("end of line/input or index");
	qi::rule<Iterator, void(), SkipType> lineRule2 = qi::eps > lineRule;

	qi::on_error<qi::fail>(lineRule2, OnError<Iterator>(err, row));

	while(!in.atEnd()) {
		QString str(in.readLine());
		Vector line(str.toUcs4());
		++row;

		if(!qi::phrase_parse(line.constBegin(), line.constEnd(), lineRule2, skipRule)) {
			return false;
		}
	}
	return true;
}
