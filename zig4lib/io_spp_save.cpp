#include "io.h"
#include "pod.h"

bool saveZSPPFile(const ZSPPFile& sppFile, QIODevice& device, ErrorCallback& err)
{
  /*
	using namespace boost;
	device.setTextModeEnabled(true);
	QTextStream out(&device);
	out.setCodec("UTF-8");
	foreach(const shared_ptr<ZSPPEntry>& entry, sppFile.sppList) {
		out << entry->weight_ << " " << entry->alpha_ << " " << entry->column3_ << " " << entry->column4_ << " " << entry->column5_  << " ";
		if(sppFile.format == SpeciesFormat::NewFormat) {
			out << entry->generalized_rule_target_ << " " << entry->generalized_rule_exp_x_ << " " << entry->generalized_rule_exp_y_ << " ";
		}
		QString raster(entry->raster_);
		raster.replace('"', "\\\"");
		out << '"' << raster << '"' << endl;
	}
	return true;
  */
  return false;
}
