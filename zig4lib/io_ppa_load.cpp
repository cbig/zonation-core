#include "io.h"
#include "ppa_adapt.h"
#include "spirit_util.h"
#include "spirit_error.h"

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/make_shared.hpp>

namespace {
template <typename Analysis>
void addEntry(ZPPAFile& ppaFile, Analysis& analysis)
{
	ppaFile.analysisList << boost::make_shared<ZAnalysis>(analysis);
}
}

bool loadZPPAFile(ZPPAFile& ppaFile, QIODevice& device, ErrorCallback& err)
{
	device.setTextModeEnabled(true);
	QTextStream in(&device);
	in.setCodec("UTF-8");
	in.setAutoDetectUnicode(true);

	typedef QVector<Char> Vector;
	typedef Vector::ConstIterator Iterator;

	namespace spirit = boost::spirit;
	namespace phoenix = boost::phoenix;
	namespace qi = boost::spirit::qi;
	namespace unicode = boost::spirit::unicode;

	int row = 0;

	typedef CommentSkipper<Iterator> SkipType;
	SkipType skipRule(";#");
	StringParser<Iterator> stringRule('"', "#");

	qi::rule<Iterator, FilePath()> rasterRule = stringRule;
	rasterRule.name("file path");
	DoubleParser<Iterator> double_;

	// entry rules

	qi::rule<Iterator, ZLSI(), SkipType> lsiRule =
	                qi::lit("LSI") > double_ > double_ > double_ > double_;

	qi::rule<Iterator, ZLSC(), SkipType> lscRule =
	                qi::lit("LSC") > double_ > double_ > rasterRule > rasterRule;

	qi::rule<Iterator, ZLSM(), SkipType> lsmRule =
	                qi::lit("LSM") > rasterRule > double_ > double_ > double_;

	qi::rule<Iterator, ZLSB(), SkipType> lsbRule =
	                qi::lit("LSB") > rasterRule > double_ > double_ > double_ > double_;

	// line rule

	qi::rule<Iterator, void()> eoiRule = qi::eoi;
	eoiRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule =
	                -(
	                        lsiRule [phoenix::bind(&addEntry<ZLSI>, phoenix::ref(ppaFile), qi::_1)] |
	                        lscRule [phoenix::bind(&addEntry<ZLSC>, phoenix::ref(ppaFile), qi::_1)] |
	                        lsmRule [phoenix::bind(&addEntry<ZLSM>, phoenix::ref(ppaFile), qi::_1)] |
	                        lsbRule [phoenix::bind(&addEntry<ZLSB>, phoenix::ref(ppaFile), qi::_1)]
	                        ) > eoiRule;

	qi::on_error<qi::accept>(lineRule, OnError<Iterator>(err, row));

	while(!in.atEnd()) {
		Vector line(in.readLine().toUcs4());
		++row;

		qi::phrase_parse(line.constBegin(), line.constEnd(), lineRule, skipRule);
	}

	return true;
}
