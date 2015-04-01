#ifndef PTR_UTIL_H
#define PTR_UTIL_H

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <QHash>

template <typename T>
inline uint qHash(const boost::shared_ptr<T>& ptr)
{
	return qHash(ptr.get());
}

template <typename T>
inline boost::shared_ptr<T> copy_and_make_shared(T const& t)
{
	return boost::shared_ptr<T>(new T(t));
}

template <typename T>
inline boost::shared_ptr<T> deep_copy(boost::shared_ptr<T> const& t)
{
	return boost::make_shared<T>(*t);
}

#endif // PTR_UTIL_H
