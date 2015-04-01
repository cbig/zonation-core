#ifndef POQSTRING_H
#define POQSTRING_H

#include <boost/program_options.hpp>
#include <QString>

//namespace boost { namespace program_options {
		inline void validate(boost::any& v,
					  const std::vector<std::string>& values,
					  QString*, int)
		{
			using namespace boost::program_options;
			// Make sure no previous assignment to 'v' was made.
			validators::check_first_occurrence(v);
			// Extract the first string from 'values'. If there is more than
			// one string, it's an error, and exception will be thrown.
			const std::string& s = validators::get_single_string(values);
			v = boost::any(QString::fromUtf8(s.c_str()));
		}

				/*
		inline void validate(boost::any& v,
					  const std::vector<std::string>& values,
					  FilePath*, int)
		{
			using namespace boost::program_options;
			// Make sure no previous assignment to 'v' was made.
			validators::check_first_occurrence(v);
			// Extract the first string from 'values'. If there is more than
			// one string, it's an error, and exception will be thrown.
			const std::string& s = validators::get_single_string(values);
			v = boost::any(FilePath(QString::fromUtf8(s.c_str())));
		}

		template <typename Derived, typename ValueType>
		inline void validate(boost::any& v,
					  const std::vector<std::string>& values,
					  boost::detail::enum_base<Derived, ValueType>*, int i)
		{
			typedef boost::detail::enum_base<Derived, ValueType> Base;
			ValueType *ptr;
			boost::any temp;
			validate(temp, values, ptr, i);
			v = temp
		}
		*/
//} }

#endif // POQSTRING_H
