/*
 Fast and effective random number generator (uniform distrib.)
 Partially based on [1] and references therein. 

 [1] "Good Practice in (Pseudo) Random Number Generation for Bioinformatics Applications,"
 David Jones. Online: www.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf
*/

#ifndef _MSVC_
#include <fcntl.h>
#include <sys/time.h>

// for getpid, read
#include <unistd.h>

// fallback in case /dev/urandom does not exist or fails
unsigned int get_devrand_fallback(void)
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (tv.tv_sec ^ tv.tv_usec) ^ getpid();
}
unsigned int get_devrand(void)
{
  int fn;
  unsigned int r;
  fn = open("/dev/urandom", O_RDONLY);
  if (fn == -1)
    return get_devrand_fallback();
  if (read(fn, &r, 4) != 4)
    return get_devrand_fallback();
  close(fn);
  return r;
}
#else

#define _CRT_RAND_S
#include <stdlib.h>

#endif

/* Seed variables */
static unsigned int x = 123456789, 
  y = 987654321,
  z = 43219876,
  c = 6543217;

/* Init generator using /dev/urandom */
void init_randz()
{
#if defined _MSVC_ 
  rand_s(&x);
  do {
    rand_s(&y);
  } while (!y);
  rand_s(&z);
  rand_s(&c);
#else
  x = get_devrand();
  while (!(y = get_devrand()));
  z = get_devrand();
  c = get_devrand() % 698769068 + 1;
#endif
}

float randz()
{
  // This is a kiss generator for 32-bit integers
  unsigned long long t;
  unsigned int result;
  x = 314527869 * x + 1234567;
  y ^= y << 5; y ^= y >> 7; y ^= y << 22;
  t = 4294584393ULL * z + c; c = t >> 32; z = t;
  result =  x + y + z;

  return result / (float)4294967296.0;
}
