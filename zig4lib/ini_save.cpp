#include <QFile>
#include <QTextStream>
#include "ini.h"

namespace {
  /*
void writeEntry(const IniMap::const_iterator& i, QTextStream& out, QString& previousSection, bool lineEndBeforeSection) {
	if(i->section != previousSection) {
		previousSection = i->section;
		if(lineEndBeforeSection)
			out << endl;
		out << "[" << i->section << "]" << endl;
	}
	QString key(i->key);
	QString value(i->value);
	key.replace('=', "\\=");
	if(i->quoted) {
		value.replace('"', "\\\"");
		value = '"' + value + '"';
	} else {
		value.replace(';', "\\;");
	}
	out << key << " = " << value << endl;
}
  */
}

bool saveIniFile(const IniMap& map, QIODevice& device, ErrorCallback& err)
{
  /*
	device.setTextModeEnabled(true);
	QTextStream out(&device);
	out.setCodec("UTF-8"); // utf-8 without BOM
	QString previousSection;
	IniMap::const_iterator i = map.begin();
	if(i != map.end()) {
		writeEntry(i, out, previousSection, false);
		for(++i; i != map.end(); ++i) {
			writeEntry(i, out, previousSection, true);
		}
	}
	return true;
  */
  return false;
}
