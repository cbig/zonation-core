#include "io.h"
#include "pod.h"
#include "spirit_util.h"
#include "spirit_error.h"

#define BOOST_SPIRIT_UNICODE // We'll use unicode (UTF8) all throughout

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/home/support/attributes.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>

bool loadZMatrixFile(ZMatrixFile& matrixFile, QIODevice& device, ErrorCallback& err)
{
	device.setTextModeEnabled(true);
	QTextStream in(&device);
	in.setCodec("UTF-8");
	in.setAutoDetectUnicode(true);

	typedef QVector<Char> Vector;
	typedef Vector::ConstIterator Iterator;

	namespace spirit = boost::spirit;
	namespace ph = boost::phoenix;
	namespace qi = boost::spirit::qi;
	namespace unicode = boost::spirit::unicode;
	using namespace boost;

	int row = 0;
	int cols = 0; // will be set after the first line
	typedef QVector<QVector<double> > TempMatrix;
	TempMatrix tempMatrix;

	typedef CommentSkipper<Iterator> SkipType;
	SkipType skipRule("#");

	// rules
	DoubleParser<Iterator> double_;

	qi::rule<Iterator, QVector<double>(), SkipType> *rowRulePtr; // we set this this to point to the correct rule after format has been detected

	qi::rule<Iterator, QVector<double>(), SkipType> strictRowRule = qi::repeat(ph::ref(cols))[double_];
	qi::rule<Iterator, QVector<double>(), SkipType> freeRowRule = +double_;
	qi::rule<Iterator, QVector<double>(), SkipType> firstRowRule = freeRowRule
	                [
	                qi::_val = qi::_1,
	                ph::ref(cols) = ph::size(qi::_1),
	                ph::ref(rowRulePtr) = ph::val(&strictRowRule)
	                ];
	rowRulePtr = &firstRowRule;

	qi::rule<Iterator> rowEndRule = qi::eoi;
	rowEndRule.name("end of line/input");

	qi::rule<Iterator, void(), SkipType> lineRule =
	                (qi::lazy(*ph::ref(rowRulePtr)) [ph::push_back(ph::ref(tempMatrix), qi::_1)] >
	                 rowEndRule) | qi::eoi;
	lineRule.name("end of line/input or decimal number row");
	qi::rule<Iterator, void(), SkipType> lineRule2 = qi::eps > lineRule;

	qi::on_error<qi::fail>(lineRule2, OnError<Iterator>(err, row));

	while(!in.atEnd()) {
		Vector line(in.readLine().toUcs4());
		++row;
		if(!qi::phrase_parse(line.constBegin(), line.constEnd(), lineRule2, skipRule)) {
			return false;
		}
	}

	// convert temporary matrix to multi_array

	typedef ZMatrixFile::Type Matrix;
	typedef Matrix::extent_gen Extents;

	Matrix& matrix(matrixFile.matrix);
	matrix.resize(Extents()[tempMatrix.size()][cols]);

	for(size_t i = 0; i != matrix.shape()[0]; ++i) {
		for(size_t j = 0; j != matrix.shape()[1]; ++j) {
			matrix[i][j] = tempMatrix[i][j];
		}
	}

	return true;
}
