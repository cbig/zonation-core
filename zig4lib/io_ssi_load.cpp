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

//#define Z_SSI_FUSION_REPEAT(r, data, elem) (BOOST_PP_SEQ_ELEM(0, elem), BOOST_PP_SEQ_ELEM(1, elem))

namespace {
void addEntry(ZSSIFile& ssiFile, ZSSIEntry& entry)
{
	ssiFile.ssiList << boost::make_shared<ZSSIEntry>(entry);
}
}

bool loadZSSIFile(ZSSIFile& ssiFile, QIODevice& device, ErrorCallback& err)
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

	qi::rule<Iterator, ZSSIEntry(), SkipType> commonSSIRule =  double_ [phoenix::bind(&ZSSIEntry::weight_, qi::_val) = qi::_1] >
	                double_ >
	                int_    >
	                int_    >
	                double_ [phoenix::bind(&ZSSIEntry::column5_, qi::_val) = qi::_1]
	                ;
	qi::rule<Iterator, FilePath()> ssiRule = stringRule [qi::_val = qi::_1];
	ssiRule.name("file path");

	// used in determining format (notice the expectation operators)

	qi::rule<Iterator, ZSSIEntry(), SkipType> firstOldSSIRule = commonSSIRule [qi::_val = qi::_1] >>
	                                                                                                 ssiRule [phoenix::bind(&ZSSIEntry::ssi_, qi::_val) = qi::_1];
	qi::rule<Iterator, ZSSIEntry(), SkipType> firstNewSSIRule = commonSSIRule [qi::_val = qi::_1] >>
	                                                                                                 double_ [phoenix::bind(&ZSSIEntry::generalized_rule_target_, qi::_val) = qi::_1] >>
	                                                                                                                                                                                     double_ [phoenix::bind(&ZSSIEntry::generalized_rule_exp_x_, qi::_val) = qi::_1] >
	                double_ [phoenix::bind(&ZSSIEntry::generalized_rule_exp_y_, qi::_val) = qi::_1] >
	                ssiRule [phoenix::bind(&ZSSIEntry::ssi_, qi::_val) = qi::_1];

	// used after format has been determined

	qi::rule<Iterator, ZSSIEntry(), SkipType> oldSSIRule = commonSSIRule [qi::_val = qi::_1] >
	                ssiRule [phoenix::bind(&ZSSIEntry::ssi_, qi::_val) = qi::_1];
	qi::rule<Iterator, ZSSIEntry(), SkipType> newSSIRule = commonSSIRule [qi::_val = qi::_1] >
	                double_ [phoenix::bind(&ZSSIEntry::generalized_rule_target_, qi::_val) = qi::_1] >
	                double_ [phoenix::bind(&ZSSIEntry::generalized_rule_exp_x_, qi::_val) = qi::_1] >
	                double_ [phoenix::bind(&ZSSIEntry::generalized_rule_exp_y_, qi::_val) = qi::_1] >
	                ssiRule [phoenix::bind(&ZSSIEntry::ssi_, qi::_val) = qi::_1];
	qi::rule<Iterator, ZSSIEntry(), SkipType> *entryRulePtr; // we set this this to point to the correct rule after format has been detected
	qi::rule<Iterator, ZSSIEntry(), SkipType> ssiFirstEntryRule = firstNewSSIRule [phoenix::ref(entryRulePtr) = phoenix::val(&newSSIRule),
	                phoenix::bind(&ZSSIFile::format, phoenix::ref(ssiFile)) = phoenix::val(SpeciesFormat::NewFormat),
	                qi::_val = qi::_1] |
	                firstOldSSIRule [phoenix::ref(entryRulePtr) = phoenix::val(&oldSSIRule),
	                phoenix::bind(&ZSSIFile::format, phoenix::ref(ssiFile)) = phoenix::val(SpeciesFormat::OldFormat),
	                qi::_val = qi::_1]
	                ;
	entryRulePtr = &ssiFirstEntryRule;

	qi::rule<Iterator, void()> eoiRule = qi::eoi;
	eoiRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule = -qi::lazy(*phoenix::ref(entryRulePtr)) [phoenix::bind(&addEntry, phoenix::ref(ssiFile), qi::_1)] > eoiRule;

	qi::on_error<qi::accept>(lineRule, OnError<Iterator>(err, row));

	while(!in.atEnd()) {
		Vector line(in.readLine().toUcs4());
		++row;

		qi::phrase_parse(line.constBegin(), line.constEnd(), lineRule, skipRule);
	}

	return true;
}
