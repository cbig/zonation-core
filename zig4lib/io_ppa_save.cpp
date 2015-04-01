#include "io.h"
#include "ppa_adapt.h"
#include <boost/fusion/include/sequence.hpp>

namespace {

// TODO: fix this! next(i) is not going to move the iterator
/*
struct WriteAnalysis : boost::static_visitor<>
{
	WriteAnalysis(QTextStream& s) : s(s) {}
	QTextStream& s;

	template <typename Analysis>
	void operator()(const Analysis& analysis)
	{
		using namespace boost::fusion;
		typedef typename result_of::begin<const Analysis>::type Iterator;
		Iterator i = begin(analysis);
		if(i != end(analysis))
			write(deref(i));
		for(next(i); i != end(analysis); next(i)) {
			s << " ";
			write(deref(i));
		}
		s << endl;
	}

	void write(double value)
	{
		s << value;
	}

	void write(const FilePath& path)
	{
		QString str(path);
		str.replace('"', "\\\"");
		s << '"' << str << '"';
	}
};
*/
}

bool saveZPPAFile(const ZPPAFile& ppaFile, QIODevice& device, ErrorCallback& err)
{
  /*
	using namespace boost;
	device.setTextModeEnabled(true);
	QTextStream out(&device);
	out.setCodec("UTF-8");
	WriteAnalysis write(out);
	foreach(const shared_ptr<ZAnalysis>& analysis, ppaFile.analysisList)
	{
		apply_visitor(write, *analysis);
	}
	return true;
  */
  return false;
}

