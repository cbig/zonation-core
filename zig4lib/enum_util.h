#ifndef ENUM_H
#define ENUM_H

#include <boost/enum.hpp>
#include <QHash>
#include <QDebug>
#include <QTextStream>
#include <QDataStream>

#define ENUM_INSTANCE(Enum, index) Enum(Enum::index)
#define ENUM_VALUE(Enum, index) ENUM_INSTANCE(Enum, index).value()

template <typename Enum>
inline typename Enum::value_type valueFromIndex(typename Enum::domain& index)
{
	return Enum(index).value();
}

template <typename Derived, typename ValueType>
inline uint qHash(const boost::detail::enum_base<Derived, ValueType>& v)
{
	return qHash(v.index());
}

template <typename Derived, typename ValueType>
inline QDebug operator<<(QDebug dbg, const boost::detail::enum_base<Derived, ValueType>& v)
{
	return dbg.nospace() << v.str();
}

template <typename Derived, typename ValueType>
inline QTextStream& operator<<(QTextStream& out, const boost::detail::enum_base<Derived, ValueType>& v)
{
	return out << v.str();
}

template <typename Derived, typename ValueType>
inline QDataStream& operator<<(QDataStream& out, const boost::detail::enum_base<Derived, ValueType>& v)
{
	qDebug() << "saving enum:" << v.str() << v.index();
	return out << static_cast<quint64>(v.index());
}

template <typename Derived, typename ValueType>
inline QDataStream& operator>>(QDataStream& in, boost::detail::enum_base<Derived, ValueType>& v)
{
	quint64 i;
	in >> i;
	qDebug() << "loading enum:" << i;
	typename Derived::optional opt(Derived::get_by_index(i));
	if(opt) {
		v = *opt;
	} else {
		in.setStatus(QDataStream::ReadCorruptData);
	}
	return in;
}

#endif // ENUM_H
