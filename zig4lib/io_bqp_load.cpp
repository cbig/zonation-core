#include "io.h"
#include "pod.h"
#include "spirit_util.h"
#include "ptr_util.h"
#include "spirit_error.h"

#include <string>

#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(
                ZBQPPoint,
                (double, x)
                (double, y)
                )

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

// we should probably put this in a better place - a header, maybe?
namespace boost { namespace spirit { namespace traits {
// shared pointer
template <typename T>
struct clear_value<shared_ptr<T> >
{
	static void call(shared_ptr<T>& ptr)
	{
		qDebug() << "boost::spirit::traits::call(shared_ptr<T>&)";
		// intialize default
		ptr = make_shared<T>();
		qDebug() << "boost::spirit::traits::call(shared_ptr<T>&): done";
	}
};
}}}

bool loadZBQPFile(ZBQPFile& bqpFile, QIODevice& device, ErrorCallback& err)
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
	int index = 0;

	typedef CommentSkipper<Iterator> SkipType;
	SkipType skipRule("#");

	// rules
	qi::rule<Iterator> indexRule;

	auto setRuleName  = [&indexRule, &index]() {
	  ++ph::ref(index);
	  indexRule.name(std::to_string(index));
	};

	indexRule =
	  // set rule name (will be propagated on error)
	  qi::eps[ setRuleName ]
	  >>
	  // pass on match
	  qi::int_ [
		    qi::_pass = (qi::_1 == ph::ref(index))
		    ];

	DoubleParser<Iterator> double_;

	qi::rule<Iterator, ZBQPPoint(), SkipType> pointRule =
	                double_ > double_;

	// we basically have two ways of lazily pointing to a member variable:
	//
	// // declaration & definition
	// struct A { int b } a;
	// // getting the value of a member variable b of reference of a
	// ph::bind(&A::b, ph::ref(a)); // first one
	// &ph::ref(a)->*&A::b; // second one
	//
	// ... second one is probably better... ->*& ftw

	qi::rule<Iterator, ZBQPCurve(), SkipType> curveRule =
	                indexRule >
	                *pointRule
	                [ph::push_back(&qi::_val->*&ZBQPCurve::pointList, qi::_1)];

	qi::rule<Iterator> curveEndRule = qi::eoi;
	curveEndRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule =
	                (curveRule [ph::push_back(&ph::ref(bqpFile)->*&ZBQPFile::curveList, ph::bind(&copy_and_make_shared<ZBQPCurve>, qi::_1))] >
	                 curveEndRule) | qi::eoi;
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
