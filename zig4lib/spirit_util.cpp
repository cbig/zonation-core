#include "spirit_util.h"

template <typename Iterator>
EscapedParser<Iterator>::EscapedParser(const char *escapedChars,
                                       bool acceptWhiteSpace) :
        EscapedParser::base_type(escapedRule)
{
	namespace spirit = boost::spirit;
	namespace qi = boost::spirit::qi;
	namespace unicode = boost::spirit::unicode;

	charRule = unicode::char_(escapedChars);
	escapeRule = spirit::lit('\\') >> charRule;
	if(acceptWhiteSpace) {
		escapedRule = escapeRule | (unicode::char_ - charRule);
	} else {
		escapedRule = escapeRule | (unicode::char_ - charRule - unicode::white_space);
	}
}

template <typename Iterator>
QuotedStringParser<Iterator>::QuotedStringParser(char begin, char end) :
        QuotedStringParser::base_type(stringRule)
{
	init(begin, end);
}

template <typename Iterator>
QuotedStringParser<Iterator>::QuotedStringParser(char begin) :
        QuotedStringParser::base_type(stringRule)
{
	init(begin, begin);
}

template <typename Iterator>
void QuotedStringParser<Iterator>::init(char begin, char end)
{
	namespace spirit = boost::spirit;
	namespace qi = boost::spirit::qi;
	namespace unicode = boost::spirit::unicode;

	endStringRule = spirit::lit(end);
	endStringRule.name(std::string(1, end));
	escapedRule = spirit::lit('\\') >> unicode::char_(end);
	stringCharRule = escapedRule | (unicode::char_ - spirit::lit(end));
	stringRule = spirit::lit(begin) >> *stringCharRule > endStringRule;
}

template <typename Iterator>
UnquotedStringParser<Iterator>::UnquotedStringParser(const char *escapedChars,
                                                     bool acceptWhitespace,
                                                     bool acceptZeroLength) :
        UnquotedStringParser::base_type(unquotedStringRule),
        unquotedCharRule(escapedChars, acceptWhitespace)
{
	if(acceptZeroLength) {
		unquotedStringRule = *unquotedCharRule;
	} else {
		unquotedStringRule = +unquotedCharRule;
	}
}

template <typename Iterator>
StringParser<Iterator>::StringParser(char quoteChar, const char *escapedChars, bool acceptUnquotedWs) :
        StringParser::base_type(stringRule),
        unquotedStringRule(escapedChars, acceptUnquotedWs),
        quotedStringRule(quoteChar)
{
	stringRule = quotedStringRule | unquotedStringRule;
}

template <typename Iterator>
CommentSkipper<Iterator>::CommentSkipper(const char *commentChars) :
        CommentSkipper::base_type(start)
{
	namespace unicode = boost::spirit::unicode;
	start = (unicode::char_(commentChars) >> *unicode::char_) | unicode::white_space;
}

template <typename Iterator>
DoubleParser<Iterator>::DoubleParser(const char *name) :
        DoubleParser::base_type(start)
{
	namespace qi = boost::spirit::qi;
	start = qi::double_;
	start.name(name);
}

template <typename Iterator>
IntParser<Iterator>::IntParser(const char *name) :
        IntParser::base_type(start)
{
	namespace qi = boost::spirit::qi;
	start = qi::int_;
	start.name(name);
}

// explicit instantiations

template class EscapedParser<DefaultIterator>;
template class QuotedStringParser<DefaultIterator>;
template class UnquotedStringParser<DefaultIterator>;
template class StringParser<DefaultIterator>;
template class CommentSkipper<DefaultIterator>;
template class DoubleParser<DefaultIterator>;
template class IntParser<DefaultIterator>;

/*
template boost::spirit::qi::rule<DefaultIterator> commentSkipper<DefaultIterator>(const char *commentChars);
template boost::spirit::qi::rule<DefaultIterator, double()> doubleParser<DefaultIterator>(const char *name);
template boost::spirit::qi::rule<DefaultIterator, int()> intParser<DefaultIterator>(const char *name);
*/
