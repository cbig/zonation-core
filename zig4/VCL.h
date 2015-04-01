/*
 *	Compatibility layer for Borland environment used by Zonation
 *  Jarno Leppänen (jarno.leppanen@helsinki.fi)
 *
 */
#ifndef ZIG2_VCL_H
#define ZIG2_VCL_H

#include <QtGlobal>
#include <QString>

#include <string>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <memory>
#include <boost/shared_array.hpp>
#ifndef __GNUC__
#include <direct.h>
#endif
#include "zig4lib/pod.h"
#include "zig4lib/error_util.h"
#include "typedefs.h"

// Enable only if really needed. Can slow down some parts of the code
//#define ZCORE_DEBUG

typedef int TColor;

extern unsigned int nthreads;

using std::pow;

typedef int TCursor;
const TCursor crHourGlass = 0xfffffff5;
const TCursor crDefault = 0x0;

#ifndef M_PI
const double M_PI = 3.14159265358979323846;  /* pi */
#endif

//prototype for doing nothing
void doNothing();

class EConvertError
{
};

/*
class String : public String
{
public:
 String() {}
 String(String& s): String(s) {}
 String(char* s): String(s) {}
 String(const String& s): String(s) {}
 String(const char* s): String(s) {}
 ~String() {}
 char *c_str();
 const char *c_str() const;
 char *operator()();
 const char *operator()() const;
};
*/

class TObject
{
};

class TComponent : public TObject
{
public:
	void Show() {}
};

class TPageControl
{
public:
	int ActivePageIndex;
};

class TTabSheet
{
};

class TEdit
{
public:
	String Text;
};

class TDialog
{
public:
	bool Execute() { return true; }
};

class TOpenDialog : public TDialog
{
};

class TSavePictureDialog : public TDialog
{
public:
	String FileName;
	int FilterIndex;
};

class TSaveDialog : public TDialog
{
public:
	String FileName;
};

class TAxis
{
public:
	void SetMinMax(double a, double b) {}
};

class TChart
{
private:
	TAxis BottomAxis_obj;
public:
	TAxis *BottomAxis;
	bool Visible;
	TChart(): BottomAxis(&BottomAxis_obj), Visible(true) {}
	void SaveToMetafileEnh(const String &file) {}
};

class TLabeledEdit
{
public:
	String Text;
};

class TLinePen
{
public:
	int Width;
};

class TCloseAction
{
};

class TLabel
{
public:
	String Caption;
};

class TBevel
{
};

class TItems
{
public:
	void Add(const String& str) {}
};

class TComboBox
{
private:
	TItems Items_obj;
public:
	TItems *Items;
	int ItemIndex;
	TComboBox(): Items(&Items_obj) {}
	void Clear() {}
};

class TPanel
{
};

class TGroupBox
{
};

class TUpDown
{
};

class TButton
{
public:
	bool Enabled;
};

class TSpecMap
{
	TButton Button3_obj;
public:
	TButton *Button3;
	TSpecMap(): Button3(&Button3_obj) {}
	void Show() {}
};

enum TUpDownDirection { updNone, updUp, updDown };

class TAckForm : public TComponent
{
};

class TReferenceForm : public TComponent
{
};

class TDisclaimerForm : public TComponent
{
};

class TVersionForm : public TComponent
{
};

class TVertScrollBar
{
public:
	int Position;
};

class TPen
{
public:
	int Width;
	TColor Color;
};

enum TFontStyle {fsBold, fsItalic, fsUnderline, fsStrikeOut};

class TFontStyles
{
public:
	String Name;
	int Size;
	TFontStyles& operator<<(TFontStyle style) { return *this; }
	TFontStyles& operator=(const TFontStyles& style) { return *this; }
};

class TFont
{
public:
	TFontStyles Style;
	TColor Color;
	String Name;
	int Size;
};

class TForm : public TComponent
{
private:
	TVertScrollBar VertScrollBar_obj;
public:
	TVertScrollBar *VertScrollBar;
	void (*closeAction)();
	String Caption;
	bool Visible;
	int Height;
	TForm(TComponent *owner = 0): VertScrollBar(&VertScrollBar_obj), closeAction(&doNothing) {}
	void Close() {
		(*closeAction)();
	}
};

void ShowMessage(const QString& msg);

int MessageBox(int i, const QString& msg, const QString& title, int b);
int MessageBox(const QString& msg, const QString& title, int b);

const int ID_YES = 2;
const int ID_NO = 3;
const int ID_OK = 4;

const int MB_OK = 0;
const int MB_YESNO = 1;

// in ARGB format
TColor RGB(TColor r, TColor g, TColor b);
TColor BGR2ARGB(TColor bgr);

const TColor clBlack = BGR2ARGB(0x0);
const TColor clMaroon = BGR2ARGB(0x80);
const TColor clGreen = BGR2ARGB(0x8000);
const TColor clOlive = BGR2ARGB(0x8080);
const TColor clNavy = BGR2ARGB(0x800000);
const TColor clPurple = BGR2ARGB(0x800080);
const TColor clTeal = BGR2ARGB(0x808000);
const TColor clGray = BGR2ARGB(0x808080);
const TColor clSilver = BGR2ARGB(0xc0c0c0);
const TColor clRed = BGR2ARGB(0xff);
const TColor clLime = BGR2ARGB(0xff00);
const TColor clYellow = BGR2ARGB(0xffff);
const TColor clBlue = BGR2ARGB(0xff0000);
const TColor clFuchsia = BGR2ARGB(0xff00ff);
const TColor clAqua = BGR2ARGB(0xffff00);
const TColor clLtGray = BGR2ARGB(0xc0c0c0);
const TColor clDkGray = BGR2ARGB(0x808080);
const TColor clWhite = BGR2ARGB(0xffffff);

class TRect {
public:
	int x;
	int y;
	TRect(int x1, int y1, int x2, int y2): x(x1 + 1), y(y1 + 1) {}
};

//forward declarations
class TPixels;
namespace Zig2 {
class Instance;
}

namespace ZInterprocess {
class Core;
}

class TLineSeries {
private:
	TLinePen LinePen_obj;

	// local data
	int plotIndex;

public:
	TLinePen *LinePen;
	//plot does nothing by default
	TLineSeries();
	void Clear() {}
	void init(int plotIndex); //synchronizes plot to zig2 shared memory command as plot with index plotIndex
	void AddXY(float x, float y, const char *, const TColor &);
};

class ImageObj {
public:
	// arguments
	int x;
	int y;
	TColor col;

	int width;
	int height;

	explicit ImageObj();
	ImageObj& operator[](int y) {
		this->y = y;
		return *this;
	}
	void operator=(TColor col) {
		this->col = col;
		setPixel(x, y, col);
	}

	// sets pixel color without updating
	void setPixel(int x, int y, TColor col);
	// update progress and draw colormap and rankmap
	void removeSite(int x, int y, TColor color, float rank);
private:

	enum State {
		SETPIXEL,
		REMOVESITE
	};

	State state;

};

class TPixels {
public:
	ImageObj Image;
	explicit TPixels() {}
	ImageObj& operator[](int x) {
		Image.x = x;
		return Image;
	}
};

class TProgressBarPosition {
private:
	int localPosition; //local position, if shared position is not defined
public:
	TProgressBarPosition() {}
	void operator=(int position) {}
};

class TProgressBar {
public:
	TProgressBarPosition Position;
};

class TBarSeries : public TLineSeries
{
};

class TBrush {
public:
	TColor Color;
};

class TCanvas {
private:
	TFont Font_obj;
	TPen Pen_obj;
	TBrush Brush_obj;
public:
	TFont *Font;
	TPen *Pen;
	TBrush *Brush;
	TPixels Pixels;
	explicit TCanvas(): Font(&Font_obj), Pen(&Pen_obj), Brush(&Brush_obj) {}
	void FillRect(const TRect& rect) {
		Pixels[rect.x][rect.y] = Brush->Color;
	}
	void TextOutA(int x, int y, String str) {}
	void MoveTo(int x, int y) {}
	void LineTo(int x, int y) {}
};

class Length {
private:
	int *len;
	Length();
public:
	Length(int *length): len(length) {}
	void operator=(int val);
	bool operator!=(int val);
};

namespace Graphics {
class TBitmap {
	TCanvas Canvas_obj;
public:
	TCanvas *Canvas;
	int &Width;
	int &Height;
	explicit TBitmap():
		Canvas(&Canvas_obj),
		Width(Canvas->Pixels.Image.width),
		Height(Canvas->Pixels.Image.height)
		  {};

	void SaveToFile(const String& file) {}
	// saves bitmap in formats configured in VCLInterface.outputImageFormats. not in borland vcl
	void SaveToFiles(const String& base);
};
}

// Saves image (int RGB buffer) into a file
// Requires a buffer of xsize*ysize bytes
void 
save_image_to_file(String emffn, int xsize, int ysize, int** buffer);

// saves float grid in formats configured in VCLInterface.outputImageFormats. not in borland vcl
template <typename T>
void SaveToRaster(const String& base, T *data, float nodatavalue, int xsize, int ysize, int linespace = 0);

class TJPEGImage
{
private:
	Graphics::TBitmap *bitmap;
public:
	TJPEGImage(): bitmap(0) {}
	void Assign(Graphics::TBitmap *bitmap);
	void SaveToFile(const String& file) {}
};

class TPicture {
	Graphics::TBitmap Bitmap_obj;
public:
	Graphics::TBitmap *Bitmap;
	explicit TPicture(): Bitmap(&Bitmap_obj) {}
	void SaveToFile(const String& file) {}
	void Repaint() {}
};

class TImage {
	TPicture Picture_obj;
public:
	int Width;
	int Height;
	TPicture* Picture;
	explicit TImage(): Picture(&Picture_obj) {}
	void Repaint() {}
};

class TLines {
private:
	String *Text;
public:
	void Add(const QString& msg);
	explicit TLines(String *pText): Text(pText) {}
};

class TMemo {
public:
	String Text;
private:
	TLines Lines_obj;
public:
	TLines *Lines;
	explicit TMemo(): Text(), Lines_obj(&Text), Lines(&Lines_obj) {}
};

class TCheckBox {
public:
	bool Checked;
	explicit TCheckBox() : Checked(true) {}
};

class TRadioGroup {
public:
	int ItemIndex;
	explicit TRadioGroup() : ItemIndex(0) {}
};

class TScreen {
public:
	TCursor Cursor;
	explicit TScreen(): Cursor(crDefault) {}
};

class TApplication {
public:
	int MessageBox(const String& msg, const String& title, int b);
	void ProcessMessages() {}
	void Run() {}
};

//Borland zig2 classes

class TAnnotForm : public TComponent {
	TCheckBox CheckBox1_obj;
	TCheckBox CheckBox2_obj;
	TCheckBox CheckBox3_obj;
	TCheckBox CheckBox4_obj;
	TCheckBox CheckBox5_obj;
	TCheckBox CheckBox6_obj;
	TEdit Edit1_obj;
	TEdit Edit2_obj;
	TEdit Edit3_obj;
	TEdit Edit4_obj;
	TEdit Edit5_obj;
	TEdit Edit6_obj;
	TEdit Edit7_obj;
	TEdit Edit8_obj;
	TEdit Edit9_obj;
	TEdit Edit10_obj;
	TEdit Edit11_obj;
	TRadioGroup RadioGroup1_obj;
public:
	TCheckBox *CheckBox1;
	TCheckBox *CheckBox2;
	TCheckBox *CheckBox3;
	TCheckBox *CheckBox4;
	TCheckBox *CheckBox5;
	TCheckBox *CheckBox6;
	TEdit *Edit1;
	TEdit *Edit2;
	TEdit *Edit3;
	TEdit *Edit4;
	TEdit *Edit5;
	TEdit *Edit6;
	TEdit *Edit7;
	TEdit *Edit8;
	TEdit *Edit9;
	TEdit *Edit10;
	TEdit *Edit11;
	TRadioGroup *RadioGroup1;
	explicit TAnnotForm():
		CheckBox1(&CheckBox1_obj),
		CheckBox2(&CheckBox2_obj),
		CheckBox3(&CheckBox3_obj),
		CheckBox4(&CheckBox4_obj),
		CheckBox5(&CheckBox5_obj),
		CheckBox6(&CheckBox6_obj),
		Edit1(&Edit1_obj),
		Edit2(&Edit2_obj),
		Edit3(&Edit3_obj),
		Edit4(&Edit4_obj),
		Edit5(&Edit5_obj),
		Edit6(&Edit6_obj),
		Edit7(&Edit7_obj),
		Edit8(&Edit8_obj),
		Edit9(&Edit9_obj),
		Edit10(&Edit10_obj),
		Edit11(&Edit11_obj),
		RadioGroup1(&RadioGroup1_obj)
	{}
};

class TForm1;

extern TForm1		   *Form1; // this is actually defined in unit1
extern TAckForm        *AckForm;
extern TReferenceForm  *ReferenceForm;
extern TAnnotForm      *AnnotForm;
extern TDisclaimerForm *DisclaimerForm;
extern TVersionForm    *VersionForm;
extern TSpecMap        *SpecMap;

// typedef Zig::IniFile TIniFile;

String IntToStr(long long int i);
String IntToStrW(long long int i, int width);
String FloatToStr(float f);

enum TFloatFormat { FF_FIXED };
const TFloatFormat ffFixed = FF_FIXED;

String FloatToStrF(float f, TFloatFormat format, int precision, int digits);
String FloatToStrF(double f, TFloatFormat format, int precision, int digits);
String ChangeFileExt(const String& path, const String& extension = "");
String createSubDirIfNeeded(String const& path, String const& append);
String checkAvailableDisk(String const& path, unsigned int& avail);
void printf(QString str);

void Randomize();
long random(long range);

class BatRun {
public:
	int retval;
	BatRun(): retval(0) {}
	void run();
};

float StrToFloat(const String& str);

int StrToInt(const String& str);

inline double pow(double a, float b)
{
	return pow(a, static_cast<double>(b));
}

inline double pow(float a, double b)
{
	return pow(static_cast<double>(a), b);
}

inline double pow(int a, int b)
{
	return pow(static_cast<double>(a), static_cast<double>(b));
}

inline int fabs(int x)
{
	return x < 0 ? -x : x;
}

using std::min;
using std::max;

inline void randomize()
{
}

extern TScreen *Screen;
extern TApplication *Application;
extern char DecimalSeparator;

//extern char *_argv[];
extern String originalOutFile;
//extern int _argc;

void initializeVCL(const ZCommandLine& command);
void deleteSingletons();

String compileMsgList();
void screenMessage(const String& msg);
void memoMessage(const String& msg);
void interprocessInitMap(int siteCount, int width, int height);
void interprocessNotifyHardWorkBegin();
void interprocessNotifyHardWorkEnd();

class ZErrorCallBack: public ErrorCallback
{
 public:
  virtual void operator()(const QString& str, zeLevel lvl = zeWarning);
};

class ZErrorStderrCallBack: public ErrorCallback
{
 public:
  virtual void operator()(const QString& str, zeLevel lvl = zeWarning);
};

class VCLCommandLine 
{
 public:
  static void readOptions(const ZCommandLine& command);

  // these methods return 0 if their corresponding options are not set
  static unsigned int hardware_or_opt_concurrency();
  static unsigned int removal_rule();
  static unsigned int warp_factor();

 private:
  VCLCommandLine();
  static unsigned int rrule;
  static unsigned int wfactor;
};

String zCurrentTime(void);

// start counting time
void zTimerStart(void);
// how much time elapsed since last call to zTimerStart()
qint64 zElapsedTime(void);

String zHostname(void);

#endif //ZIG_2_VCL
