#include "io.h"
#include "pod.h"
#include "spirit_util.h"
#include "spirit_error.h"
#include "spirit_enum.h"
#include "ptr_util.h"

#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(
                ZInteractionsEntry,
                (int, l1)
                (int, l2)
                (double, alpha)
                (InteractionsMode, iatype)
                (double, gamma)
                )

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/home/phoenix/container.hpp>
#include <boost/spirit/home/support/attributes.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/make_shared.hpp>

bool loadZInteractionsFile(ZInteractionsFile& interactionsFile, QIODevice& device, ErrorCallback& err)
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
	EnumParser<Iterator, InteractionsMode> iatypeRule;

	IntParser<Iterator> int_;
	DoubleParser<Iterator> double_;

	qi::rule<Iterator, ZInteractionsEntry(), SkipType> interactionsRule =
	                int_ > int_ > double_ > iatypeRule > double_;

	qi::rule<Iterator> interactionsEndRule = qi::eoi;
	interactionsEndRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule =
	                (interactionsRule [ph::push_back(&ph::ref(interactionsFile)->*&ZInteractionsFile::interactionsList, ph::bind(&copy_and_make_shared<ZInteractionsEntry>, qi::_1))] >
	                 interactionsEndRule) | qi::eoi;
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
