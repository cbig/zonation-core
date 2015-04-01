#ifndef SPIRIT_UTIL_H
#define SPIRIT_UTIL_H

#include "qstring_spirit.h"

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>

typedef QVector<Char>::ConstIterator DefaultIterator;

template <typename Iterator>
class EscapedParser : public boost::spirit::qi::grammar<Iterator, Char()>
{
	public:
	EscapedParser(const char *escapedChars = "\"", bool acceptWhiteSpace = true);

	private:
	boost::spirit::qi::rule<Iterator, Char()> charRule;
	boost::spirit::qi::rule<Iterator, Char()> escapeRule;
	boost::spirit::qi::rule<Iterator, Char()> escapedRule;
};

template <typename Iterator>
class QuotedStringParser : public boost::spirit::qi::grammar<Iterator, QString()>
{
	public:
	QuotedStringParser(char begin, char end);
	QuotedStringParser(char begin = '"');

	private:
	void init(char begin, char end);

	boost::spirit::qi::rule<Iterator> endStringRule;
	EscapedParser<Iterator> stringEscapeRule;
	boost::spirit::qi::rule<Iterator, Char()> escapedRule;
	boost::spirit::qi::rule<Iterator, Char()> stringCharRule;
	boost::spirit::qi::rule<Iterator, QString()> stringRule;
};

template <typename Iterator>
class UnquotedStringParser : public boost::spirit::qi::grammar<Iterator, QString()>
{
	public:
	UnquotedStringParser(const char *escapedChars = ";#", bool acceptWhitespace = false, bool acceptZeroLength = false);

	private:
	EscapedParser<Iterator> unquotedCharRule;
	boost::spirit::qi::rule<Iterator, QString()> unquotedStringRule;
};

template <typename Iterator>
class StringParser : public boost::spirit::qi::grammar<Iterator, QString()>
{
	public:
	StringParser(char quoteChar = '"', const char *escapedChars = "#", bool acceptUnquotedWs = false);

	private:
	UnquotedStringParser<Iterator> unquotedStringRule;
	QuotedStringParser<Iterator> quotedStringRule;
	boost::spirit::qi::rule<Iterator, QString()> stringRule;
};

// grammar implementation

template <typename Iterator>
class CommentSkipper : public boost::spirit::qi::grammar<Iterator>
{
public:
	CommentSkipper(const char *commentChars = "#");
private:
	boost::spirit::qi::rule<Iterator> start;
};

template <typename Iterator>
class DoubleParser : public boost::spirit::qi::grammar<Iterator, double()>
{
	public:
	DoubleParser(const char *name = "decimal number");
	private:
	boost::spirit::qi::rule<Iterator, double()> start;
};

template <typename Iterator>
class IntParser : public boost::spirit::qi::grammar<Iterator, int()>
{
	public:
	IntParser(const char *name = "integer number");
	private:
	boost::spirit::qi::rule<Iterator, int()> start;
};

// template function implementation

/*
template <typename Iterator>
boost::spirit::qi::rule<Iterator> commentSkipper(const char *commentChars = "#");

template <typename Iterator>
boost::spirit::qi::rule<Iterator, double()> doubleParser(const char *name = "decimal number");
template <typename Iterator>
boost::spirit::qi::rule<Iterator, int()> intParser(const char *name = "integer number");
*/

#endif // SPIRIT_UTIL_H
