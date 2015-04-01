#include "config.h"
#include "io.h"
#include "pod.h"
#include "ini.h"
#include "enum_util.h"
#include "poqstring.h"
#include "zinterprocess.h"
#include "io_arg.h"
#include "env.h"
#include "split_winmain.h"
#include <QFile>
#include <QList>
#include <QSet>
#include <QSettings>
#include <QVariant>
#include <QTemporaryFile>
#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>

namespace {
bool parseBool(const QString& str, bool& ok) {
	ok = true;
	if(str == "1" || str.toLower() == "true")
		return true;
	if(str == "0" || str.toLower() == "false")
		return false;
	ok = false;
	return false;
}

QString errorRow(int row) {
	return row == -1 ? QString() : QString("line %1: ").arg(row);
}

BOOST_ENUM_VALUES(ZBatFileParameter, const char *,
                  (PrintHelp)("help,h")
                  (PrintVersion)("version,v")
                  (UseThreads)("use-threads")
		  (RemovalRule)("removal-rule")
		  (WarpFactor)("warp-factor")
                  (ImageFormats)("image-output-formats")
                  (GridFormats)("grid-output-formats")
                  );

const char * batFileParam(ZBatFileParameter::domain index)
{
	return ZBatFileParameter(index).value();
}
}

void
print_help(QString command)
{
  QString vers =  QString("Zonation ") + zonation_VERSION_MAJOR + "." +
    zonation_VERSION_MINOR + "." zonation_VERSION_PATCH + ", build: " + __DATE__ + " " __TIME__;

  std::cout << "zig4 (" << command.toStdString() << ")\n"
"Usage: zig4 run_mode settingsfile spp_file output_files_prefix IGa useSmooth SmoothMult autoclose0_or_1 OPTIONS\n"
"Run Zonation core, " << vers.toStdString() <<
"\n"
"Parameters:\n"
" run_mode -r,-l               -r (run mode, calculate new solution),\n"
"                              -l FILE (reload mode, load solution in rank FILE).\n"
" settingsfile                 settings (.dat) file as specified in the manual.\n"
" spp_file                     biodiversity features list file.\n"
" output_files_prefix          name (prefix) of output files\n"
" IGa                          value of uncertainty parameter (alpha).\n"
" useSmooth                    use distribution smoothing (1) or not (0).\n"
" SmoothMult                   multiplying factor for species-specific dispersal\n"
"                              kernel widths.\n"
" autoclose0_or_1              close window at the end (1) or not (0).\n"
"\n"
"Supported options:\n"
" -h,--help                    print help and exit.\n"
" -v,--version                 print version information and exit.\n"
"\n"
" --grid-output-formats LIST   LIST is a list of formats separated by spaces.\n"
"                              Supported formats: asc tif img compressed-tif.\n"
" --image-output-formats LIST  LIST is a list of formats separated by spaces.\n"
"                              Supported formats: png bmp jpg emf.\n"
" --removal-rule N             set removal rule to N (overriding settings\n"
"                              file), 0 - CAZ, 1 - ABF, etc.\n"
" --warp-factor N              set warp factor to N (integer value).\n"
" --use-threads N              use N hardware threads (best option\n"
"                              automatically chosen if N not given).\n"
"\n";
}


void
print_version()
{
  // TODO: should there be a global version string for core, gui and everybody?
  std::cout << 
    "    zig4 version: " << zonation_VERSION_MAJOR << "." <<
    zonation_VERSION_MINOR << "." << zonation_VERSION_PATCH <<  
    ", build: " << __DATE__ << " " << __TIME__ << "\n"
"\n"
"    zig4 - Zonation computational core\n"
"    Copyright (C) 2011-2014 Conservation Biology Informatics Group\n"
"    Copyright (C) 2004-2011 Atte Moilanen\n"
"\n"
"    This program is free software: you can redistribute it and/or modify\n"
"    it under the terms of the GNU General Public License as published by\n"
"    the Free Software Foundation, either version 3 of the License, or\n"
"    (at your option) any later version.\n"
"\n"
"    This program is distributed in the hope that it will be useful,\n"
"    but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"    GNU General Public License for more details.\n"
"\n"
"    You should have received a copy of the GNU General Public License\n"
"    along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
"\n";
}

// returns false if any of the parameters fail to parse
bool parseZCommandLine(ZCommandLine& command, const QStringList& arguments, ErrorCallback& err)
{
  if(arguments.size() >= 2) {
    if (0 == arguments[1].compare("-h") || 0 == arguments[1].compare("--help")) {
      print_help(arguments[0]);
      return false;
    }
    
    if (0 == arguments[1].compare("-v") || 0 == arguments[1].compare("--version")) {
      print_version();
      return false;
    }
  }

  ZError::setErrorCallback(err);

  if(arguments.size() < 9) {
    // old, uninformative err msg.
    // ZError::err_msg(QObject::tr("not enough arguments\n"), zeError);
    ZError::err_msg(QString("Incorrect number of parameters (%1) on command line").arg(arguments.size()), zeError);
    ZError::err_msg("Expected run_mode settings_file input_grid_file output_grid_file IGa use_smooth smooth_mult autoclose_0_or_1_ignored", zeError);
    ZError::err_msg("Type zig4 -h or --help to get a quick help", zeError);
    return false;
  }
  QString binary_name = arguments[0];
#ifndef Q_OS_WIN
  // zig4.exe ---> zig4 in non-windoze systems
  // so traditional .bat files with a 'call zig4.exe ...' will work!
  QString exe_ext = ".exe";
  if (binary_name.endsWith(exe_ext))
    binary_name.chop(exe_ext.length());
  // enforce a 4 (zig4 instead of zig3, zig2 whatever else)
  binary_name.replace(binary_name.length()-1, "4");
#else
  // zig3.exe / zig3.exe -> zig4.exe
  binary_name.replace(binary_name.length()-1-4, "4.exe");
#endif
  command.zigexec = binary_name;
  QString rankModeArg = arguments[1];
  command.datFile = arguments[2];
  command.sppFile = arguments[3];
  command.outFile = arguments[4];
  bool ok = true;
  bool parsingOk = true;
  double tempDouble;
  int tempBool;
  if(rankModeArg.startsWith(ENUM_VALUE(ZCommandLineMode, CalculateRank))) {
    command.commandLineMode = ENUM_INSTANCE(ZCommandLineMode, CalculateRank);
  } else if(rankModeArg.startsWith(ENUM_VALUE(ZCommandLineMode, LoadRank))) {
    command.commandLineMode = ENUM_INSTANCE(ZCommandLineMode, LoadRank);
    command.rankFile = rankModeArg.mid(ENUM_VALUE(ZCommandLineMode, LoadRank).size());
  } else {
    ZError::err_msg(QObject::tr("could not parse rank calculation mode, defaulted to ") + ENUM_VALUE(ZCommandLineMode, CalculateRank), zeError);
    command.commandLineMode = ENUM_INSTANCE(ZCommandLineMode, CalculateRank);
    parsingOk = false;
  }
  tempDouble = arguments[5].toDouble(&ok);
  if(ok) {
    command.uncertaintyAlpha = tempDouble;
  } else {
    ZError::err_msg(QObject::tr("could not parse uncertainty alpha"), zeError);
    parsingOk = false;
  }
  tempBool = parseBool(arguments[6], ok);
  if(ok) {
    command.distributionSmoothingOn = tempBool;
  } else {
    ZError::err_msg(QObject::tr("could not parse distribution smoothing parameter"), zeError);
    parsingOk = false;
  }
  tempDouble = arguments[7].toDouble(&ok);
  if(ok) {
    command.dispersalKernelMultiplier = tempDouble;
  } else {
    ZError::err_msg(QObject::tr("could not parse dispersal kernel multiplier"), zeError);
    parsingOk = false;
  }
  tempBool = parseBool(arguments[8], ok);
  if(ok) {
    command.windowLeftOpen = tempBool;
  } else {
    ZError::err_msg(QObject::tr("could not parse window open parameter"), zeError);
    parsingOk = false;
  }

  // additional arguments

  QList<QByteArray> arrs;
  std::vector<char *> argv;

  arrs << arguments[0].toUtf8();
  for(QStringList::const_iterator i = arguments.begin() + 9; i != arguments.end(); ++i) {
    arrs << i->toUtf8();
  }
  for(QList<QByteArray>::iterator i = arrs.begin(); i != arrs.end(); ++i) {
    argv.push_back(i->data());
  }

  int argc = arrs.size();

  if(argc > 1) {

    using namespace boost::program_options;

    options_description desc("Allowed options");

    // Note ->implicit_value() is required for the following behavior:
    // --use-thread N  (with N optional, if missing: use as many as possible)
    desc.add_options()
      (batFileParam(ZBatFileParameter::PrintVersion), "Print version")
      (batFileParam(ZBatFileParameter::PrintHelp), "Print help")
      (batFileParam(ZBatFileParameter::UseThreads), value<unsigned>()->implicit_value(0),
       "Use threads. If number of threads is not specified, try to use an optimal amount")
      (batFileParam(ZBatFileParameter::WarpFactor), value<unsigned>(), "Set warp factor (overrides option in .dat file)")
      (batFileParam(ZBatFileParameter::RemovalRule), value<unsigned>(), "Set removal rule (overrides option in .dat file)")
      (batFileParam(ZBatFileParameter::ImageFormats), value<std::vector<QString> >()->multitoken()->zero_tokens(),
       "Image output formats (png, bmp, jpg or emf, default emf jpg)")
      (batFileParam(ZBatFileParameter::GridFormats), value<std::vector<QString> >()->multitoken()->zero_tokens(),
       "Grid output formats (img, compressed-img, tif, compressed-tif or asc, default asc)")
      ;

    variables_map vm;
    try {
      store(parse_command_line(argc, &argv[0], desc), vm);
    } catch(const boost::program_options::error& poerr) {
      ZError::err_msg(QObject::tr("can not parse optional arguments: %1").arg(poerr.what()), zeError);
      return false;
    }
    notify(vm);

    if(vm.count(batFileParam(ZBatFileParameter::UseThreads))) {
      command.numberOfThreads = vm[batFileParam(ZBatFileParameter::UseThreads)].as<unsigned>();
    }

    if(vm.count(batFileParam(ZBatFileParameter::RemovalRule))) {
      command.removalRule = vm[batFileParam(ZBatFileParameter::RemovalRule)].as<unsigned>();
    }

    if(vm.count(batFileParam(ZBatFileParameter::WarpFactor))) {
      command.warpFactor = vm[batFileParam(ZBatFileParameter::WarpFactor)].as<unsigned>();
    }

    if(vm.count(batFileParam(ZBatFileParameter::GridFormats))) {
      const std::vector<QString>& params(vm[batFileParam(ZBatFileParameter::GridFormats)].as<std::vector<QString> >());
      command.gridOutputFormats = QSet<ZGridFormat>();
      for(std::vector<QString>::const_iterator i = params.begin(); i != params.end(); ++i) {
	ZGridFormat::optional opt(ZGridFormat::get_by_value(*i));
	if(opt) {
	  command.gridOutputFormats->insert(*opt);
	} else {
	  ZError::err_msg(QObject::tr("unrecognized grid format %1").arg(*i), zeError);
	  parsingOk = false;
	}
      }
    }
    if(vm.count(batFileParam(ZBatFileParameter::ImageFormats))) {
      const std::vector<QString>& params(vm[batFileParam(ZBatFileParameter::ImageFormats)].as<std::vector<QString> >());
      command.imageOutputFormats = QSet<ZImageFormat>();
      for(std::vector<QString>::const_iterator i = params.begin(); i != params.end(); ++i) {
	ZImageFormat::optional opt(ZImageFormat::get_by_value(*i));
	if(opt) {
	  command.imageOutputFormats->insert(*opt);
	} else {
	  ZError::err_msg(QObject::tr("unrecognized image format %1").arg(*i), zeError);
	  parsingOk = false;
	}
      }
    }
  }

  return parsingOk;
}

// traditional loading of .bat files. Only supports simple 'call zig4...' lines
// The real stuff/parsing is done in 'parseZCommandLine'
// As a result batFile is filled in. It is just a QVector of commands
bool loadZBATFile(ZBATFile& batFile, QIODevice& device, ErrorCallback& err)
{
  using namespace boost;

  // assumed that command is in its own line and that there are no trailing arguments that
  // are not part of z command

  // utf-8 assumed (as if "chcp 65001" was set when running the bat file from command line)
  device.setTextModeEnabled(true);
  int row = 0;
  while(!device.atEnd()) {
    QByteArray l = device.readLine().constData();
    ++row;
    if (l.trimmed().isEmpty()) {
      continue;
    }
    std::string line(l);
    std::vector<std::string> words(util::split_winmain(line));
    // skip REM
    if(!words.empty() && QString::fromUtf8(words[0].c_str()).toLower() == "rem")
      continue;
    // search for "-r" beginning from second word
    typedef std::vector<std::string>::const_iterator Iterator;
    Iterator i(std::find(words.begin() + 1, words.end(), "-r"));
    // or otherwise try looking for "-l"
    if (words.end() == i) {
      std::string load_prefix("-l");
      for (i = words.begin(); i != words.end(); i++) {
	if (!i->compare(0, load_prefix.size(), load_prefix))  // found!
	  break;
      }
    }

    QStringList arguments;
    for(--i; i != words.end(); ++i) {
      arguments << QString::fromUtf8(i->c_str());
    }
    
    shared_ptr<ZCommandLine> command(make_shared<ZCommandLine>());
    // LineModify mod(err, QString("row %1: ").arg(row));
    if(parseZCommandLine(*command, arguments, err)) {
      qDebug() << "loadZBATFile(): Parsing command line, zig exec: " << command->zigexec << ", workingDir is: " << command->workingDir;
      command->line = QString::fromStdString(line);
      batFile.commandList << command;
    } else {
      ZError::err_msg(QObject::tr("Could not parse command in project (.bat) file: '%1'").
		      arg(QString::fromStdString(line).trimmed()), zeError);
    }
  }
  return true;
}

// Used for the traditional files and the "arg. forwarding mode", no difference.
bool saveZBATFile(const ZBATFile& batFile, QIODevice& device, ErrorCallback& err)
{
	using namespace boost;
	// assumed that command is in its own line and that there are no trailing arguments that
	// are not part of z command
	device.setTextModeEnabled(true);
	// utf-8 assumed (as if "chcp 65001" was set when running the bat file from command line)
	QTextStream out(&device);
	out.setCodec("UTF-8");
	foreach(const shared_ptr<ZCommandLine>& command, batFile.commandList) {
		CommandLine cmdLine(getZCommandLineArguments(*command));
		QStringList::const_iterator i(cmdLine.arguments.begin());
		out << cmdLine.program << " ";
		out << *i << " \""; // -r
		out << *++i << "\" \""; // datFile
		out << *++i << "\" \""; // sppFile
		out << *++i << "\""; // outFile
		for(++i; i != cmdLine.arguments.end(); ++i) {
			out << " " << *i;
		}
		out << endl;
	}
	return true;
}

// Used for the "arg. forwarding mode":
// It creates the arg server, calls cmd.exe with some env. variables, etc.
/*
bool loadZBATFile(ZBATFile& batFile, FilePath const& file, ErrorCallback err)
{
	using namespace boost;

	// create unique server name

	QTemporaryFile temp;
	if(!temp.open()) {
		err(QObject::tr("could not create temporary file"));
		return false;
	}
	QString argServerName(ZInterprocess::getServerName(temp.fileName()));

	// initialize socket listener

	ArgServerLoop loop(argServerName, err);

	// execute batch file

	QProcess process;
	process.setProcessChannelMode(QProcess::MergedChannels);
	process.setWorkingDirectory(DirPath(file));
	QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
	env.insert(ZEnvironment::forwardArgsVariable, "1"); // value can be anything
	env.insert(ZEnvironment::argServerNameVariable, argServerName);
	process.setProcessEnvironment(env);
#if defined(Q_OS_WIN) // on windows use cmd.exe
	process.start("cmd.exe", QStringList() << "/c" << file.absoluteFilePath().nativePathName());
#else // otherwise use bash
	// This is severely broken. Don't expect it to work
	process.start("bash", QStringList() << file.absoluteFilePath().nativePathName());
#endif

	if(!process.waitForFinished(ZEnvironment::totalTimeoutMsecs)) {
		err(QObject::tr("batch file execution timeout"));
		return false;
	}

	loop.close();
	// check process output
	//qDebug() << QTime::currentTime().toString() << QString::fromUtf8(process.readAll().constData());

	if(!loop.success) {
		err(QObject::tr("argument interception failed"));
		return false;
	}

	// parse intercepted arguments

	foreach(const ArgModeData& arg, loop.args) {

		shared_ptr<ZCommandLine> command(make_shared<ZCommandLine>());
		if(parseZCommandLine(*command, arg.arguments, err)) {
			// replace appfile
			command->zigexec = arg.appFilePath;
			command->workingDir = arg.workingDir;
			batFile.commandList << command;
		}
	}

	qDebug() << "quitting";

	return true;
}
*/

// for saveZBATFile, also used in observer and processviewmodel
CommandLine getZCommandLineArguments(const ZCommandLine& command)
{
	CommandLine cmdLine;
	cmdLine.program = command.zigexec.nativePathName();
	switch(command.commandLineMode.index()) {
	case ZCommandLineMode::CalculateRank:
		cmdLine.arguments << ENUM_VALUE(ZCommandLineMode, CalculateRank);
		break;
	case ZCommandLineMode::LoadRank:
		cmdLine.arguments << (ENUM_VALUE(ZCommandLineMode, LoadRank) + command.rankFile);
		break;
	}
	cmdLine.arguments
	                << command.datFile.nativePathName()
	                << command.sppFile.nativePathName()
	                << command.outFile.nativePathName()
	                << QString::number(command.uncertaintyAlpha)
	                << (command.distributionSmoothingOn ? "1" : "0")
	                << QString::number(command.dispersalKernelMultiplier)
	                << (command.windowLeftOpen ? "1" : "0");
	if(command.numberOfThreads) {
		cmdLine.arguments << QString("--").append(batFileParam(ZBatFileParameter::UseThreads));
		cmdLine.arguments << QString::number(*command.numberOfThreads);
	}
	if(command.removalRule) {
	  cmdLine.arguments << QString("--").append(batFileParam(ZBatFileParameter::RemovalRule));
	  cmdLine.arguments << QString::number(*command.removalRule);
	}
	if(command.warpFactor) {
	  cmdLine.arguments << QString("--").append(batFileParam(ZBatFileParameter::WarpFactor));
	  cmdLine.arguments << QString::number(*command.warpFactor);
	}
	if(command.gridOutputFormats) {
		cmdLine.arguments << QString("--").append(batFileParam(ZBatFileParameter::GridFormats));
		foreach(const ZGridFormat& format, *command.gridOutputFormats) {
			cmdLine.arguments << format.value();
		}
	}
	if(command.imageOutputFormats) {
		cmdLine.arguments << QString("--").append(batFileParam(ZBatFileParameter::ImageFormats));
		foreach(const ZImageFormat& format, *command.imageOutputFormats) {
			cmdLine.arguments << format.value();
		}
	}
	return cmdLine;
}

namespace {

// convert string to anything

template <typename T>
bool convertValueFromString(QString& str, T& t) {
	QTextStream s(&str, QIODevice::ReadOnly);
	T temp;
	s >> temp;
	if(s.status() == QTextStream::ReadCorruptData) {
		return false;
	}
	t = temp;
	return true;
}

bool convertFromString(QString& str, int& t) {
	return convertValueFromString<int>(str, t);
}

bool convertFromString(QString& str, double& t) {
	return convertValueFromString<double>(str, t);
}

bool convertFromString(QString& str, FilePath& t) {
	t = FilePath(str);
	return true;
}

bool convertFromString(QString& str, QString& t) {
	t = str;
	return true;
}


bool convertFromString(QString& str, bool& t) {
	if(str == "1" || str.toLower() == "true") {
		t = true;
		return true;
	}
	if(str == "0" || str.toLower() == "false") {
		t = false;
		return true;
	}
	return false;
}

template <typename Derived, typename T>
bool convertFromString(QString& str, boost::detail::enum_base<Derived, T>& t) {
	typedef boost::detail::enum_base<Derived, T> Base;
	QTextStream s(&str, QIODevice::ReadOnly);
	T value;
	if(convertValueFromString<T>(str, value)) {
		boost::optional<Derived> base(Base::get_by_value(value));
		if(base) {
			t = *base;
			return true;
		}
	}
	return false;
}
}

#define Z_DAT_LOAD_INIH_KEY(elem) BOOST_PP_SEQ_ELEM(2, elem) + QString(".") + BOOST_PP_SEQ_ELEM(3, elem)

// This macro will read in the options (when used in loadZDatFile)
#define Z_DAT_LOAD_INIH_REPEAT(r, data, elem) \
	if((i = index.find(boost::make_tuple(BOOST_PP_SEQ_ELEM(2, elem), BOOST_PP_SEQ_ELEM(3, elem)))) != index.end()) { \
	BOOST_PP_SEQ_ELEM(0, elem) value; \
	QString valueString(i->value); \
	if(convertFromString(valueString, value)) { \
	datFile.BOOST_PP_SEQ_ELEM(1, elem) = value; \
	} else { \
	  ZError::err_msg(QObject::tr("unable to parse option %1").arg(Z_DAT_LOAD_INIH_KEY(elem)), zeWarning); \
	} \
	index.erase(i); \
	}

bool loadZDATFile(ZDATFile& datFile, QIODevice& device, ErrorCallback& err)
{
	IniMap map;

	if(!loadIniFile(map, device, err)) {
	  ZError::err_msg(QObject::tr("could not parse dat file"), zeError);
	  return false;
	}

	IniMapByKey& index(map.get<IniMapByKeyTag>());
	IniMapByKey::iterator i;

	// convert

	BOOST_PP_SEQ_FOR_EACH(Z_DAT_LOAD_INIH_REPEAT, ~, Z_DAT_VARIABLES)

	  // unrecognized
	  
	  for(IniMap::const_iterator i = map.begin(); i != map.end(); ++i) {
	    ZError::err_msg(QObject::tr("ignoring unrecognized option in settings file, section [%1]: %2").arg(i->section).arg(i->key), zeNotice);
	  }
	return true;
}

namespace {
template <typename T>
struct ZDatVariableIniTraits
{
	static const bool quoted = false;
};

template <>
struct ZDatVariableIniTraits<FilePath>
{
	static const bool quoted = true;
};

QTextStream& operator<<(QTextStream& out, bool v)
{
	return out << (v ? "true" : "false");
}
}

#define Z_DAT_SAVE_INIH_REPEAT(r, data, elem) \
	if(datFile.BOOST_PP_SEQ_ELEM(1, elem)) { \
	QString str; \
	QTextStream stream(&str, QIODevice::WriteOnly); \
	stream << *datFile.BOOST_PP_SEQ_ELEM(1, elem); \
	map.push_back(IniEntry(BOOST_PP_SEQ_ELEM(2, elem), BOOST_PP_SEQ_ELEM(3, elem), str, ZDatVariableIniTraits<BOOST_PP_SEQ_ELEM(0, elem)>::quoted)); \
	}

bool saveZDATFile(const ZDATFile& datFile, QIODevice& device, ErrorCallback& err)
{
	//static const QSettings::Format ZIniFormat = QSettings::registerFormat("ini", readZIniFile, writeZIniFile);

	IniMap map;

	BOOST_PP_SEQ_FOR_EACH(Z_DAT_SAVE_INIH_REPEAT, ~, Z_DAT_VARIABLES)

	                return saveIniFile(map, device, err);
}

bool loadZUnknownFile(ZUnknownFile& unknownFile, QIODevice& device, ErrorCallback& err)
{
	return true;
}
