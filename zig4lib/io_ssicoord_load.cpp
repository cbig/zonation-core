#include "io.h"
#include "pod.h"
#include "spirit_util.h"
#include "spirit_error.h"

#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(
                ZSSICoordinates,
                (double, x)
                (double, y)
                (double, value)
                (double, error)
                )

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

bool loadZSSICoordinatesFile(ZSSICoordinatesFile& ssiCoordinatesFile, QIODevice& device, ErrorCallback& err)
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
	DoubleParser<Iterator> double_;

	qi::rule<Iterator, ZSSICoordinates(), SkipType> coordinatesRule =
	                double_ > double_ > double_ > double_;

	qi::rule<Iterator> coordinatesEndRule = qi::eoi;
	coordinatesEndRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule =
	                (coordinatesRule [ph::push_back(&ph::ref(ssiCoordinatesFile)->*&ZSSICoordinatesFile::coordinatesList, qi::_1)] >
	                 coordinatesEndRule) | qi::eoi;
	lineRule.name("end of line/input or decimal number");
	qi::rule<Iterator, void(), SkipType> lineRule2 = qi::eps > lineRule;

	qi::on_error<qi::fail>(lineRule2, OnError<Iterator>(err, row));

	while(!in.atEnd()) {
		Vector line(in.readLine().toUcs4());
		++row;
		if(!qi::phrase_parse(line.constBegin(), line.constEnd(), lineRule2, skipRule)) {
			return false;
		}
	}
	return true;
}
