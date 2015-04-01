#ifndef SPIRIT_ENUM_H
#define SPIRIT_ENUM_H

#include <boost/enum.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <QTextStream>
#include <QString>

template <typename Iterator, typename Enum>
class EnumParser : public boost::spirit::qi::grammar<Iterator, Enum()>
{
	public:
	EnumParser() : EnumParser::base_type(enumRule)
	{
		namespace qi = boost::spirit::qi;
		namespace ph = boost::phoenix;
		namespace ln = ph::local_names;
		valueRule = qi::auto_;
		// build rule semantics
		enumRule = valueRule [
		                ph::let(ln::_a = ph::bind(&Enum::get_by_value, qi::_1))[
		                ph::if_(!!ln::_a)[
		                // NOTHING fucking works XDDD
		                //qi::_val = ph::bind(&boost::optional<Enum>::get_value_or, ln::_a, Enum())
		                //qi::_val = ph::bind(&boost::optional<Enum>::get, ln::_a)
		                //qi::_val = ph::construct<Enum>(*ln::_a)
		                //qi::_val = ph::construct<Enum>()
		                //ph::static_cast_<const Enum&>(*ln::_a)
		                qi::_val = ph::bind(&getCopy, ln::_a)
		                //qi::_pass = true
		                ].else_[
		                qi::_pass = false
		                ]
		                ]
		                ];
		// set name
		QString str;
		{
			QTextStream s(&str);
			typename Enum::const_iterator i = Enum::begin();
			if(i != Enum::end()) {
				s << i->value();
			}
			for(++i; i != Enum::end(); ++i) {
				s << QObject::tr(" or ") << i->value();
			}
		}
		enumRule.name(str.toUtf8().constData());
	}

	private:
	// must NOT be called with empty v
	static Enum getCopy(const boost::optional<Enum>& v)
	{
		return *v;
	}

	typedef typename Enum::value_type ValueType;
	boost::spirit::qi::rule<Iterator, ValueType()> valueRule;
	boost::spirit::qi::rule<Iterator, Enum()> enumRule;
};

#endif // SPIRIT_ENUM_H
