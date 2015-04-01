#ifndef _RANDZ_H_
#define _RANDZ_H_

/* Simple random number generator for Z */
/* It is used in bat_run.cpp and loaddata.cpp for... */
extern "C" {
  void init_randz();
  float randz();
}
#endif
