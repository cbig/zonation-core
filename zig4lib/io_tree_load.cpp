#include "io.h"
#include "pod.h"
#include "spirit_util.h"
#include "spirit_error.h"
#include "ptr_util.h"

#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(
                ZTreeEntry,
                (int, id)
                (int, down_id)
                )

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/home/support/attributes.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/make_shared.hpp>

bool loadZTreeFile(ZTreeFile& treeFile, QIODevice& device, ErrorCallback& err)
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

	qi::rule<Iterator, ZTreeEntry(), SkipType> treeRule =
	                int_ > int_;

	qi::rule<Iterator> treeEndRule = qi::eoi;
	treeEndRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule =
	                (treeRule [ph::push_back(&ph::ref(treeFile)->*&ZTreeFile::treeList, ph::bind(&copy_and_make_shared<ZTreeEntry>, qi::_1))] >
	                 treeEndRule) | qi::eoi;
	lineRule.name("end of line/input or index");
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

