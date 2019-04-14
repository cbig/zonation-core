#include "io.h"
#include "pod.h"
#include "spirit_util.h"
#include "spirit_error.h"
#include "spirit_enum.h"
#include "ptr_util.h"

#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(
                ZGroupsEntry,
                (int, num)
                (int, cond_num)
                (int, ret_num)
                (RetentionMode, ret_mode)
                (int, arb_kernel_ds)
                (int, arb_kernel_matrix)
                (int, arb_kernel_ia)
                )

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>

bool loadZGroupsFile(ZGroupsFile& groupsFile, QIODevice& device, ErrorCallback& err)
{
  // assumes an external parser is run before this
  return true;

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
	EnumParser<Iterator, RetentionMode> retRule;

	IntParser<Iterator> int_;

	qi::rule<Iterator, ZGroupsEntry(), SkipType> groupsRule =
	                int_ > int_ > int_ > retRule > int_;

	qi::rule<Iterator> groupsEndRule = qi::eoi;
	groupsEndRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule =
	                (groupsRule [ph::push_back(&ph::ref(groupsFile)->*&ZGroupsFile::groupsList, ph::bind(&copy_and_make_shared<ZGroupsEntry>, qi::_1))] >
	                 groupsEndRule) | qi::eoi;
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
