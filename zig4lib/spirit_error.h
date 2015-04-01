#ifndef SPIRIT_ERROR_H
#define SPIRIT_ERROR_H

#include "error_util.h"
#include <boost/fusion/container/vector.hpp>

template <typename Iterator>
struct OnError
{
	typedef boost::fusion::vector<
	Iterator&       // first
	, Iterator const& // last
	, Iterator const& // err
	, boost::spirit::info const&>
	Parameters;

	ErrorCallback& callback;
	int& row;

	OnError(ErrorCallback& callback, int& row) :
	        callback(callback), row(row)
	{}

	template <typename Context, typename Action>
	void operator()(Parameters params, Context context, Action action) const
	{
		//Iterator& first(boost::fusion::at_c<0>(params));
		const Iterator& last(boost::fusion::at_c<1>(params));
		const Iterator& err(boost::fusion::at_c<2>(params));
		const boost::spirit::info& info(boost::fusion::at_c<3>(params));

		int size = last - err;
		static const int MAX_SIZE = 20;

		std::stringstream ss;
		ss << info;
		/* callback(QObject::tr("expecting %1 on row %2: %3")
		         .arg(QString::fromUtf8(ss.str().c_str())).arg(row)
		         .arg(QString::fromUcs4(err, qMin(size, MAX_SIZE)) + ((size > MAX_SIZE) ? "..." : ""))
		         );*/
		ZError::err_msg(QObject::tr("expecting %1 on row %2: %3")
		         .arg(QString::fromUtf8(ss.str().c_str())).arg(row)
		         .arg(QString::fromUcs4(err, qMin(size, MAX_SIZE)) + ((size > MAX_SIZE) ? "..." : ""))
				, zeNotice);
	}
};
#endif // SPIRIT_ERROR_H
