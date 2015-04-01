#ifndef QSTRING_SPIRIT_H
#define QSTRING_SPIRIT_H

#include <boost/spirit/include/support_string_traits.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <QString>
#include <QVector>

typedef uint Char;

class CharString : public QVector<Char>
{
public:
	explicit CharString(const QString& str) : QVector<Char>(str.toUcs4()) {}
};

class CharStringIterator : public boost::iterator_facade<
                CharStringIterator,
                const Char,
                boost::random_access_traversal_tag
                >
{
public:
	explicit CharStringIterator(const QString& str, bool end = false) : cstr(CharString(str)), it(end ? cstr.end() : cstr.begin()) {}
	explicit CharStringIterator(const CharString& str, bool end = false) : cstr(str), it(end ? cstr.end() : cstr.begin()) {}
	void setEnd() { it = cstr.end(); }
private:
	friend class boost::iterator_core_access;

	CharString cstr;
	CharString::const_iterator it;
	static const Char null = 0;

	// works as a c_str too
	const Char& dereference() const {
		if(it == cstr.end())
			return null;
		return *it;
	}
	bool equal(const CharStringIterator& other) const { return it == other.it; }
	void increment() { ++it; }
	void decrement() { --it; }
	void advance(difference_type n) { it += n; }
	difference_type distance_to(const CharStringIterator& other) const { return other.it - it; }
};


namespace boost { namespace spirit { namespace traits
{

// CONTAINER

// Make Qi recognize QString as a container
template <> struct is_container<QString> : mpl::true_ {};

// Expose the container's (QString's) value_type
template <> struct container_value<QString> : mpl::identity<Char> {};

// Define how to insert a new element at the end of the container (QString)
template <> struct push_back_container<QString, Char>
{
	static bool call(QString& s, const Char& v)
	{
		s.append(QString::fromUcs4(&v, 1));
		return true;
	}
};

// STRING

// Make Qi recognize QString as a string
template <> struct is_string<QString> : mpl::true_ {};

template <> struct char_type_of<QString> : mpl::identity<Char> {};

inline CharStringIterator get_c_string(const QString& str)
{ return CharStringIterator(str); }

inline CharStringIterator get_begin(const QString& str)
{ return CharStringIterator(str); }

inline CharStringIterator get_end(const QString& str)
{ return CharStringIterator(str, true); }

}}}

#endif // QSTRING_SPIRIT_H
