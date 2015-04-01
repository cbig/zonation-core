#include "io.h"
#include "pod.h"
#include "spirit_util.h"
#include "spirit_error.h"

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/make_shared.hpp>

#define Z_SPP_FUSION_REPEAT(r, data, elem) (BOOST_PP_SEQ_ELEM(0, elem), BOOST_PP_SEQ_ELEM(1, elem))

namespace {
void addEntry(ZSPPFile& sppFile, ZSPPEntry& entry)
{
	sppFile.sppList << boost::make_shared<ZSPPEntry>(entry);
}
}

bool loadZSPPFile(ZSPPFile& sppFile, QIODevice& device, ErrorCallback& err)
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
	SkipType skipRule("#");

	StringParser<Iterator> stringRule('"', "#", true);

	// entry
	IntParser<Iterator> int_;
	DoubleParser<Iterator> double_;

	qi::rule<Iterator, ZSPPEntry(), SkipType> commonSppRule =  double_ [phoenix::bind(&ZSPPEntry::weight_, qi::_val) = qi::_1] >
	                double_ [phoenix::bind(&ZSPPEntry::alpha_, qi::_val) = qi::_1] >
	                int_    [phoenix::bind(&ZSPPEntry::column3_, qi::_val) = qi::_1] >
	                int_    [phoenix::bind(&ZSPPEntry::column4_, qi::_val) = qi::_1] >
	                double_ [phoenix::bind(&ZSPPEntry::column5_, qi::_val) = qi::_1]
	                ;
	qi::rule<Iterator, FilePath()> sppRasterRule = stringRule [qi::_val = qi::_1];
	sppRasterRule.name("file path");

	// used in determining format (notice the expectation operators)

	qi::rule<Iterator, ZSPPEntry(), SkipType> firstOldSppRule = commonSppRule [qi::_val = qi::_1] >>
	                                                                                                 sppRasterRule [phoenix::bind(&ZSPPEntry::raster_, qi::_val) = qi::_1];
	qi::rule<Iterator, ZSPPEntry(), SkipType> firstNewSppRule = commonSppRule [qi::_val = qi::_1] >>
	                                                                                                 double_ [phoenix::bind(&ZSPPEntry::generalized_rule_target_, qi::_val) = qi::_1] >>
	                                                                                                                                                                                     double_ [phoenix::bind(&ZSPPEntry::generalized_rule_exp_x_, qi::_val) = qi::_1] >
	                double_ [phoenix::bind(&ZSPPEntry::generalized_rule_exp_y_, qi::_val) = qi::_1] >
	                sppRasterRule [phoenix::bind(&ZSPPEntry::raster_, qi::_val) = qi::_1];

	// used after format has been determined

	qi::rule<Iterator, ZSPPEntry(), SkipType> oldSppRule = commonSppRule [qi::_val = qi::_1] >
	                sppRasterRule [phoenix::bind(&ZSPPEntry::raster_, qi::_val) = qi::_1];
	qi::rule<Iterator, ZSPPEntry(), SkipType> newSppRule = commonSppRule [qi::_val = qi::_1] >
	                double_ [phoenix::bind(&ZSPPEntry::generalized_rule_target_, qi::_val) = qi::_1] >
	                double_ [phoenix::bind(&ZSPPEntry::generalized_rule_exp_x_, qi::_val) = qi::_1] >
	                double_ [phoenix::bind(&ZSPPEntry::generalized_rule_exp_y_, qi::_val) = qi::_1] >
	                sppRasterRule [phoenix::bind(&ZSPPEntry::raster_, qi::_val) = qi::_1];
	qi::rule<Iterator, ZSPPEntry(), SkipType> *entryRulePtr; // we set this this to point to the correct rule after format has been detected
	qi::rule<Iterator, ZSPPEntry(), SkipType> sppFirstEntryRule = firstNewSppRule [phoenix::ref(entryRulePtr) = phoenix::val(&newSppRule),
	                phoenix::bind(&ZSPPFile::format, phoenix::ref(sppFile)) = phoenix::val(SpeciesFormat::NewFormat),
	                qi::_val = qi::_1] |
	                firstOldSppRule [phoenix::ref(entryRulePtr) = phoenix::val(&oldSppRule),
	                phoenix::bind(&ZSPPFile::format, phoenix::ref(sppFile)) = phoenix::val(SpeciesFormat::OldFormat),
	                qi::_val = qi::_1]
	                ;
	entryRulePtr = &sppFirstEntryRule;

	qi::rule<Iterator, void()> eoiRule = qi::eoi;
	eoiRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule = -qi::lazy(*phoenix::ref(entryRulePtr)) [phoenix::bind(&addEntry, phoenix::ref(sppFile), qi::_1)] > eoiRule;

	qi::on_error<qi::accept>(lineRule, OnError<Iterator>(err, row));

	while(!in.atEnd()) {
		Vector line(in.readLine().toUcs4());
		++row;

		qi::phrase_parse(line.constBegin(), line.constEnd(), lineRule, skipRule);
	}

	return true;
}
