#include "Unit1.h"
#include "VCL.h"
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <memory>
#include <set>
#include <boost/thread/thread.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/scoped_array.hpp>
#include "zig4lib/raster.h"
#include "zig4lib/zinterprocess.h"
#include "zig4lib/core.h"
#include <QDebug>
#include <QCoreApplication>
#include <QDateTime>
#include <QElapsedTimer>
#include <QHostInfo>

//Number of threads, overridden while reading command line parameters
// fedemp: now accessed through VCLCommanLine
unsigned int nthreads = 1;
//this is stupid but it can't be helped
const int PLOTSIZE = 10000; //the same number that is reserved for the global variable curves 
//const int PLOTS = 7; //to draw 7 curves in 4 windows
const int PLOTS = 10; //to draw 10 curves in 4 windows (3+3+3+1)

namespace {

const int UPDATE_FREQ(300);
const int BLOCK_SIZE(1024);
const size_t SIZEOF_DOUBLE(sizeof(double));

}

void doNothing() {}

class TForm1Delegate
{
public:
	TPageControl PageControl1_obj;
	TTabSheet TabSheet1_obj;
	TTabSheet TabSheet2_obj;
	TTabSheet TabSheet3_obj;
	TMemo Memo1_obj;
	TImage Image1_obj;
	TButton Button3_obj;
	TButton Button4_obj;
	TButton Button5_obj;
	TButton Button6_obj;
	TEdit Edit2_obj;
	TOpenDialog OpenDialog1_obj;
	TTabSheet TabSheet4_obj;
	TLabel Label1_obj;
	TLabel Label2_obj;
	TBevel Bevel1_obj;
	TLabel Label4_obj;
	TLabel Label22_obj;
	TLabel Label32_obj;
	TLabel Label5_obj;
	TLabel Label6_obj;
	TLabel Label20_obj;
	TLabel Label21_obj;
	TComboBox ComboBox1_obj;
	TLabel Label3_obj;
	TChart Chart1_obj;
	TLineSeries Series1_obj;
	TLineSeries Series2_obj;
	TLineSeries Series_23_obj;
	TPanel Panel1_obj;
	TImage Image3_obj;
	TSaveDialog SaveDialog1_obj;
	TEdit Edit3_obj;
	TTabSheet TabSheet5_obj;
	TLabel Label13_obj;
	TTabSheet TabSheet6_obj;
	TLabel Label15_obj;
	TLabel Label16_obj;
	TCheckBox CheckBox2_obj;
	TButton Button8_obj;
	TEdit Edit10_obj;
	TEdit Edit11_obj;
	TEdit Edit12_obj;
	TLabel Label18_obj;
	TEdit Edit13_obj;
	TComboBox ComboBox2_obj;
	TUpDown UpDown1_obj;
	TTabSheet TabSheet7_obj;
	TLabel Label37_obj;
	TEdit Edit25_obj;
	TLabel Label39_obj;
	TCheckBox CheckBox3_obj;
	TCheckBox CheckBox4_obj;
	TEdit Edit27_obj;
	TCheckBox CheckBox11_obj;
	TGroupBox GroupBox2_obj;
	TCheckBox CheckBox12_obj;
	TLabel Label23_obj;
	TEdit Edit15_obj;
	TRadioGroup RadioGroup1_obj;
	TCheckBox CheckBox13_obj;
	TEdit Edit16_obj;
	TCheckBox CheckBox14_obj;
	TRadioGroup RadioGroup2_obj;
	TEdit Edit17_obj;
	TGroupBox GroupBox1_obj;
	TLabel Label14_obj;
	TLabel Label19_obj;
	TCheckBox CheckBox8_obj;
	TEdit Edit9_obj;
	TCheckBox CheckBox9_obj;
	TCheckBox CheckBox10_obj;
	TEdit Edit14_obj;
	TEdit Edit18_obj;
	TLabel Label25_obj;
	TLabel Label26_obj;
	TLabel Label27_obj;
	TLabel Label28_obj;
	TLabel Label29_obj;
	TLabel Label30_obj;
	TEdit Edit19_obj;
	TChart Chart2_obj;
	TLineSeries LineSeries1_obj;
	TLineSeries LineSeries2_obj;
	TButton Button9_obj;
	TRadioGroup RadioGroup3_obj;
	TLineSeries Series3_obj;
	TLineSeries Series4_obj;

	TLineSeries LineSeries2_Series4_obj;
  
	TLabel Label45_obj;
	TGroupBox GroupBox3_obj;
	TLabel Label40_obj;
	TEdit Edit20_obj;
	TLabel Label7_obj;
	TEdit Edit21_obj;
	TButton Button2_obj;
	TLabel Label41_obj;
	TEdit Edit26_obj;
	TGroupBox GroupBox4_obj;
	TLabel Label33_obj;
	TEdit Edit22_obj;
	TLabel Label34_obj;
	TEdit Edit24_obj;
	TLabel Label35_obj;
	TEdit Edit23_obj;
	TButton Button1_obj;
	TLabel Label44_obj;
	TEdit Edit29_obj;
	TLabel Label36_obj;
	TCheckBox CheckBox7_obj;
	TLabel Label31_obj;
	TProgressBar ProgressBar1_obj;
	TButton Button10_obj;
	TLabeledEdit LabeledEdit1_obj;
	TLabel Label47_obj;
	TCheckBox CheckBox1_obj;
	TLabeledEdit LabeledEdit2_obj;
	TSaveDialog SaveDialog2_obj;
	TChart Chart3_obj;
	TLineSeries LineSeries3_obj;
	TLineSeries LineSeries4_obj;
	TLabel Label48_obj;
	TLabel Label49_obj;
	TLabel Label51_obj;
	TImage Image4_obj;
	TCheckBox CheckBox5_obj;
	TEdit Edit1_obj;
	TChart Chart4_obj;
	TLineSeries LineSeries5_obj;
	TLineSeries LineSeries6_obj;

	TLineSeries LineSeries_56_obj;

	TLineSeries LineSeries7_obj;
	TRadioGroup RadioGroup4_obj;
	TSavePictureDialog SavePictureDialog1_obj;
	TLabel Label38_obj;
	TCheckBox CheckBox6_obj;
	TEdit Edit30_obj;
	TLabel Label46_obj;
	TEdit Edit31_obj;
	TChart Chart5_obj;
	TBarSeries Series5_obj;
	TLabel Label50_obj;
	TLabel Label52_obj;
	TLabel Label53_obj;
	TLabel Label54_obj;
	TLabel Label55_obj;
	TLabel Label56_obj;
	TLabel Label57_obj;
	TEdit Edit32_obj;
	TEdit Edit33_obj;
	TEdit Edit34_obj;
	TEdit Edit35_obj;
	TEdit Edit37_obj;
	TEdit Edit38_obj;
	TEdit Edit39_obj;
	TButton Button11_obj;
	TCheckBox CheckBox15_obj;
	TEdit Edit40_obj;
	TCheckBox CheckBox16_obj;
	TEdit Edit41_obj;
	TLabel Label42_obj;
	TLabel Label59_obj;
	TLabel Label60_obj;
	TLabel Label61_obj;
	TLabel Label17_obj;
	TButton Button13_obj;
	TButton Button14_obj;
	TButton Button15_obj;
	TButton Button16_obj;
	TButton Button17_obj;
	TButton Button18_obj;
	TButton Button19_obj;
	TButton Button20_obj;
	TLabel Label24_obj;
	TLabel Label62_obj;
	TGroupBox GroupBox5_obj;
	TLabel Label58_obj;
	TEdit Edit36_obj;
	TButton Button12_obj;
	TGroupBox GroupBox6_obj;
	TLabel Label8_obj;
	TEdit Edit4_obj;
	TLabel Label10_obj;
	TEdit Edit5_obj;
	TLabel Label9_obj;
	TEdit Edit6_obj;
	TLabel Label11_obj;
	TEdit Edit7_obj;
	TLabel Label12_obj;
	TEdit Edit8_obj;
	TLabel Label43_obj;
	TEdit Edit28_obj;
	TButton Button7_obj;
	TRadioGroup RadioGroup5_obj;
	TLabel Label63_obj;
	TEdit Edit42_obj;

	void setup(TForm1 *Form1)
	{
		Form1->PageControl1 = &PageControl1_obj;
		Form1->TabSheet1 = &TabSheet1_obj;
		Form1->TabSheet2 = &TabSheet2_obj;
		Form1->TabSheet3 = &TabSheet3_obj;
		Form1->Memo1 = &Memo1_obj;
		Form1->Image1 = &Image1_obj;
		Form1->Button3 = &Button3_obj;
		Form1->Button4 = &Button4_obj;
		Form1->Button5 = &Button5_obj;
		Form1->Button6 = &Button6_obj;
		Form1->Edit2 = &Edit2_obj;
		Form1->OpenDialog1 = &OpenDialog1_obj;
		Form1->TabSheet4 = &TabSheet4_obj;
		Form1->Label1 = &Label1_obj;
		Form1->Label2 = &Label2_obj;
		Form1->Bevel1 = &Bevel1_obj;
		Form1->Label4 = &Label4_obj;
		Form1->Label22 = &Label22_obj;
		Form1->Label32 = &Label32_obj;
		Form1->Label5 = &Label5_obj;
		Form1->Label6 = &Label6_obj;
		Form1->Label20 = &Label20_obj;
		Form1->Label21 = &Label21_obj;
		Form1->ComboBox1 = &ComboBox1_obj;
		Form1->Label3 = &Label3_obj;
		Form1->Chart1 = &Chart1_obj;
		Form1->Series1 = &Series1_obj;
		Form1->Series2 = &Series2_obj;
		Form1->Panel1 = &Panel1_obj;
		Form1->Image3 = &Image3_obj;
		Form1->SaveDialog1 = &SaveDialog1_obj;
		Form1->Edit3 = &Edit3_obj;
		Form1->TabSheet5 = &TabSheet5_obj;
		Form1->Label13 = &Label13_obj;
		Form1->TabSheet6 = &TabSheet6_obj;
		Form1->Label15 = &Label15_obj;
		Form1->Label16 = &Label16_obj;
		Form1->CheckBox2 = &CheckBox2_obj;
		Form1->Button8 = &Button8_obj;
		Form1->Edit10 = &Edit10_obj;
		Form1->Edit11 = &Edit11_obj;
		Form1->Edit12 = &Edit12_obj;
		Form1->Label18 = &Label18_obj;
		Form1->Edit13 = &Edit13_obj;
		Form1->ComboBox2 = &ComboBox2_obj;
		Form1->UpDown1 = &UpDown1_obj;
		Form1->TabSheet7 = &TabSheet7_obj;
		Form1->Label37 = &Label37_obj;
		Form1->Edit25 = &Edit25_obj;
		Form1->Label39 = &Label39_obj;
		Form1->CheckBox3 = &CheckBox3_obj;
		Form1->CheckBox4 = &CheckBox4_obj;
		Form1->Edit27 = &Edit27_obj;
		Form1->CheckBox11 = &CheckBox11_obj;
		Form1->GroupBox2 = &GroupBox2_obj;
		Form1->CheckBox12 = &CheckBox12_obj;
		Form1->Label23 = &Label23_obj;
		Form1->Edit15 = &Edit15_obj;
		Form1->RadioGroup1 = &RadioGroup1_obj;
		Form1->CheckBox13 = &CheckBox13_obj;
		Form1->Edit16 = &Edit16_obj;
		Form1->CheckBox14 = &CheckBox14_obj;
		Form1->RadioGroup2 = &RadioGroup2_obj;
		Form1->Edit17 = &Edit17_obj;
		Form1->GroupBox1 = &GroupBox1_obj;
		Form1->Label14 = &Label14_obj;
		Form1->Label19 = &Label19_obj;
		Form1->CheckBox8 = &CheckBox8_obj;
		Form1->Edit9 = &Edit9_obj;
		Form1->CheckBox9 = &CheckBox9_obj;
		Form1->CheckBox10 = &CheckBox10_obj;
		Form1->Edit14 = &Edit14_obj;
		Form1->Edit18 = &Edit18_obj;
		Form1->Label25 = &Label25_obj;
		Form1->Label26 = &Label26_obj;
		Form1->Label27 = &Label27_obj;
		Form1->Label28 = &Label28_obj;
		Form1->Label29 = &Label29_obj;
		Form1->Label30 = &Label30_obj;
		Form1->Edit19 = &Edit19_obj;
		Form1->Chart2 = &Chart2_obj;
		Form1->LineSeries1 = &LineSeries1_obj;
		Form1->LineSeries2 = &LineSeries2_obj;
		Form1->Button9 = &Button9_obj;
		Form1->RadioGroup3 = &RadioGroup3_obj;
		Form1->Series3 = &Series3_obj;

		Form1->Series_23 = &Series_23_obj;

		Form1->Series4 = &Series4_obj;

		Form1->LineSeries2_Series4 = &LineSeries2_Series4_obj;

		Form1->Label45 = &Label45_obj;
		Form1->GroupBox3 = &GroupBox3_obj;
		Form1->Label40 = &Label40_obj;
		Form1->Edit20 = &Edit20_obj;
		Form1->Label7 = &Label7_obj;
		Form1->Edit21 = &Edit21_obj;
		Form1->Button2 = &Button2_obj;
		Form1->Label41 = &Label41_obj;
		Form1->Edit26 = &Edit26_obj;
		Form1->GroupBox4 = &GroupBox4_obj;
		Form1->Label33 = &Label33_obj;
		Form1->Edit22 = &Edit22_obj;
		Form1->Label34 = &Label34_obj;
		Form1->Edit24 = &Edit24_obj;
		Form1->Label35 = &Label35_obj;
		Form1->Edit23 = &Edit23_obj;
		Form1->Button1 = &Button1_obj;
		Form1->Label44 = &Label44_obj;
		Form1->Edit29 = &Edit29_obj;
		Form1->Label36 = &Label36_obj;
		Form1->CheckBox7 = &CheckBox7_obj;
		Form1->Label31 = &Label31_obj;
		Form1->ProgressBar1 = &ProgressBar1_obj;
		Form1->Button10 = &Button10_obj;
		Form1->LabeledEdit1 = &LabeledEdit1_obj;
		Form1->Label47 = &Label47_obj;
		Form1->CheckBox1 = &CheckBox1_obj;
		Form1->LabeledEdit2 = &LabeledEdit2_obj;
		Form1->SaveDialog2 = &SaveDialog2_obj;
		Form1->Chart3 = &Chart3_obj;
		Form1->LineSeries3 = &LineSeries3_obj;
		Form1->LineSeries4 = &LineSeries4_obj;
		Form1->Label48 = &Label48_obj;
		Form1->Label49 = &Label49_obj;
		Form1->Label51 = &Label51_obj;
		Form1->Image4 = &Image4_obj;
		Form1->CheckBox5 = &CheckBox5_obj;
		Form1->Edit1 = &Edit1_obj;
		Form1->Chart4 = &Chart4_obj;
		Form1->LineSeries5 = &LineSeries5_obj;
		Form1->LineSeries6 = &LineSeries6_obj;

		Form1->LineSeries_56 = &LineSeries_56_obj;

		Form1->LineSeries7 = &LineSeries7_obj;
		Form1->RadioGroup4 = &RadioGroup4_obj;
		Form1->SavePictureDialog1 = &SavePictureDialog1_obj;
		Form1->Label38 = &Label38_obj;
		Form1->CheckBox6 = &CheckBox6_obj;
		Form1->Edit30 = &Edit30_obj;
		Form1->Label46 = &Label46_obj;
		Form1->Edit31 = &Edit31_obj;
		Form1->Chart5 = &Chart5_obj;
		Form1->Series5 = &Series5_obj;
		Form1->Label50 = &Label50_obj;
		Form1->Label52 = &Label52_obj;
		Form1->Label53 = &Label53_obj;
		Form1->Label54 = &Label54_obj;
		Form1->Label55 = &Label55_obj;
		Form1->Label56 = &Label56_obj;
		Form1->Label57 = &Label57_obj;
		Form1->Edit32 = &Edit32_obj;
		Form1->Edit33 = &Edit33_obj;
		Form1->Edit34 = &Edit34_obj;
		Form1->Edit35 = &Edit35_obj;
		Form1->Edit37 = &Edit37_obj;
		Form1->Edit38 = &Edit38_obj;
		Form1->Edit39 = &Edit39_obj;
		Form1->Button11 = &Button11_obj;
		Form1->CheckBox15 = &CheckBox15_obj;
		Form1->Edit40 = &Edit40_obj;
		Form1->CheckBox16 = &CheckBox16_obj;
		Form1->Edit41 = &Edit41_obj;
		Form1->Label42 = &Label42_obj;
		Form1->Label59 = &Label59_obj;
		Form1->Label60 = &Label60_obj;
		Form1->Label61 = &Label61_obj;
		Form1->Label17 = &Label17_obj;
		Form1->Button13 = &Button13_obj;
		Form1->Button14 = &Button14_obj;
		Form1->Button15 = &Button15_obj;
		Form1->Button16 = &Button16_obj;
		Form1->Button17 = &Button17_obj;
		Form1->Button18 = &Button18_obj;
		Form1->Button19 = &Button19_obj;
		Form1->Button20 = &Button20_obj;
		Form1->Label24 = &Label24_obj;
		Form1->Label62 = &Label62_obj;
		Form1->GroupBox5 = &GroupBox5_obj;
		Form1->Label58 = &Label58_obj;
		Form1->Edit36 = &Edit36_obj;
		Form1->Button12 = &Button12_obj;
		Form1->GroupBox6 = &GroupBox6_obj;
		Form1->Label8 = &Label8_obj;
		Form1->Edit4 = &Edit4_obj;
		Form1->Label10 = &Label10_obj;
		Form1->Edit5 = &Edit5_obj;
		Form1->Label9 = &Label9_obj;
		Form1->Edit6 = &Edit6_obj;
		Form1->Label11 = &Label11_obj;
		Form1->Edit7 = &Edit7_obj;
		Form1->Label12 = &Label12_obj;
		Form1->Edit8 = &Edit8_obj;
		Form1->Label43 = &Label43_obj;
		Form1->Edit28 = &Edit28_obj;
		Form1->Button7 = &Button7_obj;
		Form1->RadioGroup5 = &RadioGroup5_obj;
		Form1->Label63 = &Label63_obj;
		Form1->Edit42 = &Edit42_obj;
	}

	void setupLineSeries(TForm1 *Form1)
	{
		Form1->Series2->init(0); //proportion of landcape lost,red
		Form1->Series3->init(1); //proportion of landcape lost, blue
		Form1->Series_23->init(2); // weighted: black

		Form1->LineSeries2->init(3); //cost needed to achieve given conservation value, red
		Form1->Series4->init(4);     //cost needed to achieve given conservation value, blue
		Form1->LineSeries2_Series4->init(5);  // weighted; black

		Form1->LineSeries3->init(6); //avg SA-extinction risk, blue

		Form1->LineSeries5->init(7); //SSI distr remaining, red
		Form1->LineSeries6->init(8); //SSI distr remaining, blue
		Form1->LineSeries_56->init(9); // weighted; black
	}
};

class VCLInterface {
private:
	TApplication Application_obj;
	TScreen Screen_obj;
	TForm1 Form1_obj;
	TAnnotForm AnnotForm_obj;
	TSpecMap SpecMap_obj;
	TAckForm AckForm_obj;
	TReferenceForm ReferenceForm_obj;
	TDisclaimerForm DisclaimerForm_obj;
	TVersionForm VersionForm_obj;
	
	TForm1Delegate Form1Delegate;

	VCLInterface(const VCLInterface&);
	VCLInterface& operator =(const VCLInterface& other);

public:
	VCLInterface():
		Form1_obj(0),
		Application(&Application_obj),
		Screen(&Screen_obj),
		Form1(&Form1_obj),
		AnnotForm(&AnnotForm_obj),
		SpecMap(&SpecMap_obj),
		AckForm(&AckForm_obj),
		ReferenceForm(&ReferenceForm_obj),
		DisclaimerForm(&DisclaimerForm_obj),
		VersionForm(&VersionForm_obj)
	{
		::Form1 = this->Form1;
		Form1Delegate.setup(Form1);
		Form1Delegate.setupLineSeries(Form1);
	}

	~VCLInterface();
	TApplication *Application;
	TScreen *Screen;
	TForm1 *Form1;
	TAnnotForm *AnnotForm;
	TSpecMap *SpecMap;
	TAckForm *AckForm;
	TReferenceForm *ReferenceForm;
	TDisclaimerForm *DisclaimerForm;
	TVersionForm *VersionForm;
	static VCLInterface& getInstance();

        static unsigned int hardware_concurrency();

	// private stuff that needs to be accessed within vcl.cpp

	std::set<ZGridFormat> gridOutputFormats;
	std::set<ZImageFormat> imageOutputFormats;

	boost::shared_ptr<ZInterprocess::Core> interprocessCore;
};

// Static members
unsigned int VCLCommandLine::wfactor= 0;
unsigned int VCLCommandLine::rrule= 0;

// Default initialization
VCLCommandLine::VCLCommandLine()
{
}

unsigned int
VCLCommandLine::hardware_or_opt_concurrency()
{
	return nthreads;
}

unsigned int
VCLCommandLine::warp_factor()
{
	return wfactor;
}

unsigned int
VCLCommandLine::removal_rule()
{
	return rrule;
}

void VCLCommandLine::readOptions(const ZCommandLine& command)
{
	VCLInterface &vcl(VCLInterface::getInstance());
	if(command.numberOfThreads) {
	  // 0: implicit value -> try to use the number of hardware threads
	  if (0 == *command.numberOfThreads)
	    nthreads = boost::thread::hardware_concurrency();
	  else if (*command.numberOfThreads <= boost::thread::hardware_concurrency())
	    nthreads = *command.numberOfThreads;
	  else
	    nthreads = boost::thread::hardware_concurrency();
	} else {
	  nthreads = 1;
	}

	if(command.removalRule) {
	  VCLCommandLine::rrule = *command.removalRule;
	}

	if(command.warpFactor) {
	  VCLCommandLine::wfactor = *command.warpFactor;
	}

	if(command.gridOutputFormats) {
		foreach(const ZGridFormat& format, *command.gridOutputFormats) {
			vcl.gridOutputFormats.insert(format);
		}
	} else {
	  // Default changed from traditional (<=3.1): ASC
	  vcl.gridOutputFormats.insert(ZGridFormat(ZGridFormat::CTIF));
	}
	if(command.imageOutputFormats) {
		foreach(const ZImageFormat& format, *command.imageOutputFormats) {
			vcl.imageOutputFormats.insert(format);
		}
	} else {
	  // Suppress EMF (very large files for large grids, more than 2 orders of magnitude bigger than jpg)
	  //vcl.imageOutputFormats.insert(ZImageFormat(ZImageFormat::EMF));
	  vcl.imageOutputFormats.insert(ZImageFormat(ZImageFormat::JPG));
	}

	/*
  using namespace boost::program_options;
  //ImageObj &image(Form1->Image1->Picture->Bitmap->Canvas->Pixels.Image);
  // Declare the supported options.
  options_description desc("Allowed options");
  desc.add_options()
    //("interprocess-id", value<String>(), "Interprocess identifier. If not specified, IPC is not used.")
    (batFileParam(ZBatFileParameter::UseThreads), value<unsigned int>(&nthreads)->implicit_value(boost::thread::hardware_concurrency()),
     "Use threads. If number of threads is not specified, try to use an optimal amount")
    (batFileParam(ZBatFileParameter::ImageFormats), value<std::vector<String> >()->multitoken()->zero_tokens(),
     "Image output formats (png, bmp, jpg or emf, default emf jpg)")
    (batFileParam(ZBatFileParameter::GridFormats), value<std::vector<String> >()->multitoken()->zero_tokens(),
     "Grid output formats (img, compressed-img, tif, compressed-tif or asc, default asc)")
    ;

  variables_map vm;
  try {
   store(parse_command_line(argc, argv, desc), vm);
   notify(vm);
  } catch(const error& err) {
   std::cout << "error reading extra parameters: " << err.what() << std::endl;
  }

  VCLInterface &vcl(VCLInterface::getInstance());
  if(vm.count(batFileParam(ZBatFileParameter::GridFormats))) {
   const std::vector<String>& params(vm[batFileParam(ZBatFileParameter::GridFormats)].as<std::vector<String> >());
   for(std::vector<String>::const_iterator i = params.begin(); i != params.end(); ++i) {
    ZGridFormat::optional opt(ZGridFormat::get_by_value(*i));
    if(opt) {
     vcl.gridOutputFormats.insert(*opt);
    } else {
     Form1->Memo1->Lines->Add("Warning: ignoring unrecognized grid format \"" + *i + "\"");
    }
   }
  } else {
   vcl.gridOutputFormats.insert(ZGridFormat(ZGridFormat::ASC));
  }
  if(vm.count(batFileParam(ZBatFileParameter::ImageFormats))) {
   const std::vector<String>& params(vm[batFileParam(ZBatFileParameter::ImageFormats)].as<std::vector<String> >());
   for(std::vector<String>::const_iterator i = params.begin(); i != params.end(); ++i) {
    ZImageFormat::optional opt(ZImageFormat::get_by_value(*i));
    if(opt) {
     vcl.imageOutputFormats.insert(*opt);
    } else {
     Form1->Memo1->Lines->Add("Warning: ignoring unrecognized image format \"" + *i + "\"");
    }
   }
  } else {
   vcl.imageOutputFormats.insert(ZImageFormat(ZImageFormat::EMF));
   vcl.imageOutputFormats.insert(ZImageFormat(ZImageFormat::JPG));
  }
  // interprocess thing is initialized always
  if(vm.count("interprocess-id")) {
   String base(vm["interprocess-id"].as<String>());
   QString name(QString::fromUtf8(base));

   vcl.interprocessCore = std::auto_ptr<ZInterprocess::Core>(new ZInterprocess::Core(name));
   if(vcl.interprocessCore->locked()) {
    qDebug() << "Server" << name << "already exists!";
    vcl.interprocessLock.reset();
    exit(ZINTERPROCESS_ERROR);
   }
  }
  */
}

// interprocess stuff

String compileMsgList()
{
  return VCLInterface::getInstance().interprocessCore->compileMsgList();
}

void ZErrorCallBack::operator()(const QString& str, zeLevel lvl)
{
  memoMessage(str + "\n");
}

void
ZErrorStderrCallBack::operator()(const QString& str, zeLevel lvl)
{
  //qerr << str << '\n'; //endl;
  fprintf(stderr, str.toStdString().c_str());
  fprintf(stderr, "\n");
  //qerr.flush();
}

void screenMessage(const String& msg)
{
  const boost::shared_ptr<ZInterprocess::Core>& core(VCLInterface::getInstance().interprocessCore);

  //  qout << msg << '\n'; //endl;
  fprintf(stdout, msg.toStdString().c_str());
  fprintf(stdout, "\n");
  if(core)
    core->message("showmessage: " + msg);
  //qout.flush();
}

void memoMessage(const String& msg)
{
  const boost::shared_ptr<ZInterprocess::Core>& core(VCLInterface::getInstance().interprocessCore);
  //// this line doesn't have any effect?
  ////qout << msg << '\n'; //endl;
  //qout << msg << endl;
  fprintf(stdout, msg.toStdString().c_str());
  fprintf(stdout, "\n");
  if(core)
    core->message(msg);
  //qout.flush();
}

void interprocessInitMap(int siteCount, int width, int height)
{
  VCLInterface::getInstance().interprocessCore->init(ZInterprocess::InstanceInitData(siteCount, width, height, PLOTS, PLOTSIZE));
}

void interprocessNotifyHardWorkBegin()
{
  VCLInterface::getInstance().interprocessCore->notifyInit();
}

void interprocessNotifyHardWorkEnd()
{
  VCLInterface::getInstance().interprocessCore->notifyDone();
}

namespace {
// interprocess stuff

void interprocessInit(const ZInterprocess::InstanceInitData& command)
{
  VCLInterface::getInstance().interprocessCore->init(command);
}

void interprocessSetPixel(quint32 x, quint32 y, quint32 color)
{
  VCLInterface::getInstance().interprocessCore->setPixel(x, y, color);
}

void interprocessNotifyInit()
{
  VCLInterface::getInstance().interprocessCore->notifyInit();
}

void interprocessAddPlotPoint(quint32 plotIndex, float x, float y)
{
  VCLInterface::getInstance().interprocessCore->addPlotPoint(plotIndex, x, y);
}

void interprocessRemoveSite(quint32 x, quint32 y, quint32 color, float rank)
{
  VCLInterface::getInstance().interprocessCore->removeSite(x, y, color, rank);
}

void interprocessNotifyDone()
{
  VCLInterface::getInstance().interprocessCore->notifyInit();
}

/*
const quint32 *colorMapData()
{
  return VCLInterface::getInstance().interprocessCore->colorMapData();
}
*/
}

// general stuff
/*
char *String:
{
 return const_cast<char *>(String:);
}

char *String::operator()()
{
 return const_cast<char *>(String:);
}

const char *String: const
{
 return String:;
}

const char *String::operator()() const
{
 return String:;
}
*/
// message stuff

void ShowMessage(const String& msg)
{
	screenMessage(msg);
}

int MessageBox(int i, const String& msg, const String& title, int b)
{
	screenMessage(msg);
	switch(i) {
	case MB_YESNO:
		return ID_YES;
	case MB_OK:
	default:
		return ID_OK;
	}
}

int MessageBox(const String& msg, const String& title, int b)
{
	screenMessage(msg);
	return ID_OK;
}

void TLines::Add(const String& msg)
{
	memoMessage(msg);
}

// convert stuff

String IntToStr(long long int i)
{
	return QString::number(i);
}

String IntToStrW(long long int i, int width)
{
  return QString("%1").arg(i, width);
}

String FloatToStr(float f)
{
  return QString::number(f);
}

int leftNumbers(float f)
{
  static const float logmult = 1.0f / std::log(10.0f);
  return (int)std::floor((std::log(f)*logmult));
}

String FloatToStrF(float f, TFloatFormat format, int precision, int digits)
{
  return FloatToStrF((double)f, format, precision, digits);
}

String FloatToStrF(double f, TFloatFormat format, int precision, int digits)
{
  //return QString::number(f, 'g', precision);
  return QString("%1").arg(f, digits, 'g', precision);
}

// extension given with dot, eg. ".jpg"
String ChangeFileExt(const String& path, const String& extension)
{
	return changeFileExtension(path, extension);
}

// Creates a directory = 'first parameter' + 'append' (if not yet there) under the path given in the first parameter
// example: createSubDirIfNeeded("/home/foo/z-test/bar_project", "_ADMU_outputs")  -> creates "/home/foo/z-test/bar_project_ADMU_outputs"
String createSubDirIfNeeded(String const& path, String const& append)
{
  // dirname and basename (in the POSIX sense)
  String dirname = path.left(path.lastIndexOf('/'));
  String basename = path.mid(path.lastIndexOf('/')+1);
  String full = dirname + "/" + basename + append;
  QDir d(full);

  if (!d.exists())
    QDir(dirname).mkdir(basename + append);
  return full;
}

// Qt: nothing to check available disk space
// boost: boost::filesystem::space(), -- but Path::initial_dir deprecated!
#include <boost/filesystem.hpp>
#ifdef WIN32
#include <direct.h>
#include <wchar.h>    // include this *and* add /Zc:wchar to cxx flags for MSVC
// otherwise: obscure linker errors related to boost methods.
#define getcwd _getcwd
#else
#include <unistd.h>
#endif

String checkAvailableDisk(String const& path, unsigned int& avail)
{
  // fedemp: yep, its the CWD instead of the initial dir, but...
  // anyway zig is often using relative paths all the time (the
  // outgridfn, memofn, etc are obtained from ZConf, which just uses
  // the command line parameters (that can be abs or relative, up to
  // the user).
  const size_t cwdMaxLength = 2048;
  char cwd[cwdMaxLength];
  if (!getcwd(cwd, cwdMaxLength)) {
    avail = 0;
    return "";
  }

  int posSlash = path.lastIndexOf('/');
  String full;
  if (-1 == posSlash )
    full = String(cwd);
  else
    full = String(cwd) + "/" + path.left(posSlash);
  boost::filesystem::space_info s = boost::filesystem::space(full.toStdString());
  avail = s.available;
  return full;
}

void printf(QString str)
{
  printf(str.toStdString().c_str());
}

// general stuff

void Randomize()
{
}

long random(long range)
{
	return ((long)std::rand()) % range;
}

// BatRun::run definition

extern bool bat_run();

void BatRun::run(void)
{
	retval = bat_run() ? 0 : -1;
	TCloseAction action;
	Form1->FormClose(0, action);
}

// VCL layer stuff

void Length::operator=(int val)
{
	*len = val;
}
bool Length::operator!=(int val)
{
	if(len) {
		return *len != val;
	}
	return true;
}

TLineSeries::TLineSeries() :
	plotIndex(-1),
	LinePen(&LinePen_obj)
{
}

void TLineSeries::init(int plotIndex)
{
	if(plotIndex >= 0) {
		this->plotIndex = plotIndex;
	}
}

void TLineSeries::AddXY(float x, float y, const char *, const TColor &) {
	if(plotIndex >= 0)
		interprocessAddPlotPoint(plotIndex, x, y);
	else
		throw Exception("TLineSeries not initialized");
}

ImageObj::ImageObj(): 
	x(0), y(0), col(0),
	width(0), height(0),
	state(SETPIXEL)
{
}

void ImageObj::setPixel(int x, int y, TColor col)
{
	switch(state) {
	case SETPIXEL:
		interprocessSetPixel(x, y, col);
		break;
	case REMOVESITE:
		//silently ignore - this happens on post process for example
		interprocessSetPixel(x, y, col);
		//throw Exception("ImageObj: setPixel called on REMOVESITE state");
		break;
	}
}

void ImageObj::removeSite(int x, int y, TColor color, float rank)
{
	switch(state) {
	case SETPIXEL:
		// initialization has been moved to bat_run()
		//interprocessNotifyInit();
		state = REMOVESITE;
	case REMOVESITE:
		interprocessRemoveSite(x, y, color, rank);
		break;
	}
}

void TJPEGImage::Assign(Graphics::TBitmap *bitmap)
{
	this->bitmap = bitmap;
}

// Does not work anymore, as colorMapData() is now empty
void Graphics::TBitmap::SaveToFiles(const String& base)
{
  return;
  /*
	VCLInterface &vcl(VCLInterface::getInstance());
	if(vcl.imageOutputFormats.empty())
		return;

	int xsize = Width, ysize = Height, size = xsize * ysize;
	String tempFilename(base + ".temp.bmp");
	String driver("BMP");

	{
		Dataset<unsigned char> tempDataset(driver.toUtf8().constData(), tempFilename.toUtf8().constData(), xsize, ysize, 3, noProjection); // temporary bmp file
		{
			const quint32 *bitmap = colorMapData();
			boost::scoped_array<unsigned char> plane(new unsigned char[xsize * ysize]);
			// red
			for(int i = 0; i != size; ++i) {
				quint32 pixel(bitmap[i]);
				plane[i] = (pixel >> 16) & 0xff;
			}
			tempDataset.writeBand(1, plane.get());
			// green
			for(int i = 0; i != size; ++i) {
				quint32 pixel(bitmap[i]);
				plane[i] = (pixel >> 8) & 0xff;
			}
			tempDataset.writeBand(2, plane.get());
			// blue
			for(int i = 0; i != size; ++i) {
				quint32 pixel(bitmap[i]);
				plane[i] = pixel & 0xff;
			}
			tempDataset.writeBand(3, plane.get());
		}

		for(std::set<ZImageFormat>::const_iterator i = vcl.imageOutputFormats.begin(); i != vcl.imageOutputFormats.end(); ++i) {
			String filename;
			switch(i->index()) {
			case ZImageFormat::EMF:
				driver = "BMP";
				break;
			case ZImageFormat::JPG:
				driver = "JPEG";
				break;
			case ZImageFormat::PNG:
				driver = "PNG";
				break;
			case ZImageFormat::BMP:
				driver = "BMP";
				break;
			}
			filename = base + i->value().extension();
			memoMessage("saving " + filename);
			Dataset<unsigned char> dataset(driver.toUtf8().constData(), filename.toUtf8().constData(), tempDataset);
		}
	}
	remove(tempFilename.toUtf8().constData()); // remove temporary bmp file
  */
}

String
get_driver_name_from_index(int idx)
{
  String driver;
  switch(idx) {
  case ZImageFormat::EMF:
    driver = "BMP";
    break;
  case ZImageFormat::JPG:
    driver = "JPEG";
    break;
  case ZImageFormat::PNG:
    driver = "PNG";
    break;
  case ZImageFormat::BMP:
    driver = "BMP";
    break;
  }
  return driver;
}

void 
save_image_to_file(String emffn, int xsize, int ysize, int** buffer)
{
  String driver = "MEM"; //get_driver_name_from_index(format_idx->index());
  Dataset<unsigned char> first_dataset(driver.toUtf8().constData(), "", xsize, ysize, 3, noProjection);

  boost::scoped_array<unsigned char> plane(new unsigned char[xsize * ysize]);
  // red
  size_t i = 0;
  for(size_t y=0; y<ysize; y++) {
    for(size_t x=0; x<xsize; x++) {
      plane[i++] = (buffer[y][x] >> 16) & 0xff;
    }
  }
  //tempDataset.writeBand(1, plane.get());
  first_dataset.writeBand(1, plane.get());
  // green
  i = 0;
  for(size_t y=0; y<ysize; y++) {
    for(size_t x=0; x<xsize; x++) {
      plane[i++] = (buffer[y][x] >> 8) & 0xff;
    }
  }
  //tempDataset.writeBand(2, plane.get());
  first_dataset.writeBand(2, plane.get());
  // blue
  i = 0;
  for(size_t y=0; y<ysize; y++) {
    for(size_t x=0; x<xsize; x++) {
      plane[i++] = buffer[y][x] & 0xff;
    }
  }
  //tempDataset.writeBand(3, plane.get());
  first_dataset.writeBand(3, plane.get());

  // Save also in all listed formats...
  VCLInterface &vcl(VCLInterface::getInstance());
  //  for(std::set<ZImageFormat>::const_iterator i = vcl.imageOutputFormats.begin(); i != vcl.imageOutputFormats.end(); ++i) {
  for(std::set<ZImageFormat>::const_iterator format_idx = vcl.imageOutputFormats.begin();
      format_idx != vcl.imageOutputFormats.end(); ++format_idx) {
    String filename;
    driver = get_driver_name_from_index(format_idx->index());
    filename = ChangeFileExt(emffn, "") + format_idx->value().extension();
    memoMessage("Saving image: " + filename + " ("+driver+" format)");
    Dataset<unsigned char> sec_dataset(driver.toUtf8().constData(), filename.toUtf8().constData(), first_dataset);
    //Dataset<unsigned char> dataset(driver.toUtf8().constData(), filename.toUtf8().constData(), tempDataset);
  }
  //remove(tmp_filename.toUtf8().constData()); // remove temporary bmp file
}

template <typename T>
void SaveToRaster(const String& base, T *data, float nodatavalue, int xsize, int ysize, int linespace)
{
  VCLInterface &vcl(VCLInterface::getInstance());
  if(vcl.gridOutputFormats.empty())
    return;

  //String tempFilename(base + ".temp.tif");
  //String driver("GTiff");
  String driver = "MEM";
  Dataset<T> tempDataset(driver.toUtf8().constData(), "" /*tempFilename.toUtf8().constData()*/, xsize, ysize, 1);
  tempDataset->GetRasterBand(1)->SetNoDataValue(nodatavalue);
  tempDataset.writeBand(1, data, linespace);

  for(std::set<ZGridFormat>::const_iterator i = vcl.gridOutputFormats.begin(); i != vcl.gridOutputFormats.end(); ++i) {
    String filename;
    char **papszOptions = 0;
    filename = base + i->value().extension();
    switch(i->index()) {
    case ZGridFormat::TIF:
      driver = "GTiff";
      papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
      break;
    case ZGridFormat::CTIF:
      driver = "GTiff";
      //papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "DEFLATE");
      papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
      papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
      // .compressed-tif -> .compressed.tif
      filename[filename.size()-4] = '.';
      break;
    case ZGridFormat::IMG:
      driver = "HFA";
      break;
    case ZGridFormat::CIMG:
      driver = "HFA";
      papszOptions = CSLSetNameValue(papszOptions, "COMPRESSED", "YES");
      // .compressed-img -> .compressed.img
      filename[filename.size()-4] = '.';
      break;
    case ZGridFormat::ASC:
      driver = "AAIGrid";
      papszOptions = CSLSetNameValue(papszOptions, "DECIMAL_PRECISION", "7");
      break;
    }
    memoMessage("Saving raster: " + filename + " ("+driver+" format)");
    Dataset<unsigned char> dataset(driver.toUtf8().constData(), filename.toUtf8().constData(), tempDataset, papszOptions);
    // avoid memory leak
    CSLDestroy(papszOptions);
  }
  //remove(tempFilename.toUtf8().constData()); // remove temporary grid file
}

// explicit instantiations
template void SaveToRaster<int>(const String& base, int *data, float nodatavalue, int xsize, int ysize, int linespace);
template void SaveToRaster<float>(const String& base, float *data, float nodatavalue, int xsize, int ysize, int linespace);


VCLInterface::~VCLInterface()
{
}

int TApplication::MessageBox(const String& msg, const String& title, int b)
{
	return ::MessageBox(msg, title, b);
}

float StrToFloat(const String& str)
{
	return str.toFloat(); // no error checking!
}

TColor RGB(TColor r, TColor g, TColor b) {
	return (r << 16) | (g << 8) | b | 0xff000000;
}

TColor BGR2ARGB(TColor bgr) {
	return RGB(bgr & 0xff, (bgr >> 8) & 0xff, (bgr >> 16) & 0xff);
}

int StrToInt(const String& str)
{
	return str.toInt(); // no error checking!
}

char DecimalSeparator = '.';

TAnnotForm *AnnotForm;
TSpecMap *SpecMap;
TAckForm *AckForm;
TReferenceForm *ReferenceForm;
TDisclaimerForm *DisclaimerForm;
TVersionForm *VersionForm;
TScreen *Screen;
TApplication *Application;

//typedef Loki::SingletonHolder<VCLInterface, Loki::CreateUsingNew, Loki::DeletableSingleton> VCLInterfaceSingleton;
boost::shared_ptr<VCLInterface> vclInterface;

String originalOutFile;

void deleteSingletons()
{
	vclInterface.reset();
	//Loki::DeletableSingleton<VCLInterface>::GracefulDelete();
}

VCLInterface& VCLInterface::getInstance()
{
	if(!vclInterface)
		vclInterface = boost::shared_ptr<VCLInterface>(new VCLInterface);
	return *vclInterface;
	//return VCLInterfaceSingleton::Instance();
}

void initializeVCL(const ZCommandLine& command)
{
	VCLInterface &vcl(VCLInterface::getInstance());
	String msg;
	vcl.interprocessCore = boost::shared_ptr<ZInterprocess::Core>
	  (createCore(command.outFile, msg), releaseCore);
	AnnotForm = vcl.AnnotForm;
	SpecMap = vcl.SpecMap;
	AckForm = vcl.AckForm;
	ReferenceForm = vcl.ReferenceForm;
	DisclaimerForm = vcl.DisclaimerForm;
	VersionForm = vcl.VersionForm;
	Screen = vcl.Screen;
	Application = vcl.Application;
	VCLCommandLine::readOptions(command);
}

String zCurrentTime(void)
{
  return QDateTime::currentDateTime().toString("h:mm:ss AP (yyyy/M/d) ");
}

static QElapsedTimer eTimer;

void zTimerStart(void)
{
  eTimer.start();
}

qint64 zElapsedTime(void)
{
  return eTimer.elapsed();
}

String zHostname(void)
{
  return QHostInfo::localHostName();
}
