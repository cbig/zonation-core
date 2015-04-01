#ifndef GISNODATA_H
#define GISNODATA_H

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

template <typename T>
inline typename boost::enable_if<boost::is_floating_point<T>, T>::type nodata()
{
	return std::numeric_limits<T>::quiet_NaN();
}

template <typename T>
inline typename boost::enable_if<boost::is_unsigned<T>, T>::type nodata()
{
	return 0;
}

template <typename T>
inline typename boost::enable_if<boost::is_signed<T>, T>::type nodata()
{
	return -1;
}

template <typename T>
inline bool isNodata(T v, typename boost::enable_if<boost::is_floating_point<T>, T>::type *dummy = 0)
{
  //return (boost::math::isnan<T>)(v);
  return z_isnan(v);
}

template <typename T>
inline bool isNodata(T v, typename boost::enable_if<boost::is_integral<T>, T>::type *dummy = 0)
{
	return v == nodata<T>();
}

#endif // GISNODATA_H
