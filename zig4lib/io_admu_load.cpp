#include "io.h"
#include "pod.h"
#include "spirit_util.h"
#include "spirit_error.h"
#include "ptr_util.h"

#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(
                ZADMUEntry,
                (int, fid)
                (double, glob_weight)
                (double, loc_weight)
                (QString, unit_name)
                )

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/home/phoenix/container.hpp>
#include <boost/spirit/home/support/attributes.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>

bool loadZADMUFile(ZADMUFile& admuFile, QIODevice& device, ErrorCallback& err)
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
	DoubleParser<Iterator> double_;

	StringParser<Iterator> nameRule('"', "#", true);

	qi::rule<Iterator, ZADMUEntry(), SkipType> admuRule =
	                int_ > double_ > double_ > nameRule;

	qi::rule<Iterator> admuEndRule = qi::eoi;
	admuEndRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule =
	                (admuRule [ph::push_back(&ph::ref(admuFile)->*&ZADMUFile::admuList, ph::bind(&copy_and_make_shared<ZADMUEntry>, qi::_1))] >
	                 admuEndRule) | qi::eoi;
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
