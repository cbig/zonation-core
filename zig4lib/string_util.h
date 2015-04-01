#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <QString>
#include <QHash>

inline std::size_t hash_value(const QString& str)
{
	return qHash(str);
}

#endif // STRING_UTIL_H
