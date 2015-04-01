#ifndef OCCUR_CONTAINER_H
#define OCCUR_CONTAINER_H

#include <stdint.h> // for int16_t
#include <cstdio>
#include <limits>

#ifdef COMPACT_VMAT
#warning "using compact vmat..."
#endif

#if defined COMPACT_VMAT_LOOKUP

//#define COMPACT_VMAT_LOOKUP
//#define COMPACT_VMAT_FLOATP
 #warning "using compact vmat - lookup version"

typedef uint16_t BFOC_size_t;  // up to 2**16 features
// typedef uint32_t BFOC_size_t;

// Contains occurrence levels (abundance, probability, presence/absence, etc.) for a list of biodiversity features at a location (cell).
// The compact version is meant for those cases where many values are -1 (or <0) - so the list is sparse
//
// Replaces float* in 'classic' Z (like in 'float *rowp' or the last * in 'float ***vmat').
// __packed__ (gcc) attributes reduces size from 24 (aligned to 8 in 64 bits) to 18
#ifdef __GNUC__
#define GCC_PACK_STRUCT_ATTRIBUTE __attribute__ ((__packed__)) 
#endif
class GCC_PACK_STRUCT_ATTRIBUTE Biodiv_Features_Occur_Container_compact
{
 public:
 Biodiv_Features_Occur_Container_compact(): priv_key(NULL), priv_val(NULL), /*priv_capacity(0),*/  priv_size(0)
    {};

 ~Biodiv_Features_Occur_Container_compact()
   { priv_size = 0; /*priv_capacity = 0;*/ priv_key = NULL; priv_val = NULL; }; // nope: if (priv_rowp) delete priv_rowp; 

  inline BFOC_size_t size() const
  { if (priv_key) return priv_size; else return 0; };

  /*
    // capacity: abandoned
  inline size_t capacity()
  { if (priv_key) return priv_capacity; else return 0; };
  */

  inline operator bool() const
  { return NULL != priv_val; };

  // semantics similar to std::
  // But note it is a single-use method. It does not free mem if previously allocated
  // map_cnt is not used at the moment, but could be useful as max allowed index or something like that
  void reserve(BFOC_size_t size, BFOC_size_t map_cnt)
  { 
    priv_key = new BFOC_size_t[size];  //new int[size];    
    priv_val = new float[size];
    /*priv_capacity = size;*/
    priv_size = 0;
  };

  // This does the real delete/free. delete is not in the destructor to allow lines like 
  // '_Container rowp = vmat[y][x]' => will call the destructor when rowp goes out of scope.
  void clear()
  { if (priv_key) delete [] priv_key; priv_key = NULL; if (priv_val) delete [] priv_val; priv_val = NULL; priv_size = 0; /*priv_capacity = 0;*/ };

  // fake constructor to use when allocating with malloc
  void init_constructor()
  { priv_key = NULL; priv_val = NULL; priv_size = 0; };

  // fill array with "empty" values (-1) up to position loop-1
  // In this compact structure, nothing is done (would ruin the whole thing)
  inline void fill_empty(int loop_up)
  {}

  // set position as empty (-1). Similar story as fill_empty
  inline void set_empty(int pos)
  {}

  // Search in the look-up table:
  // linear search: can be better for small sizes (cache)  <-- shouldn't apply here
  // binary, interpolation, ternary, etc.

  // translate linear (sparse) index into compact
  // assumes increasing sequence of keys!
  // returns an index out of range if key not found (out of range typically: priv_size)
  inline BFOC_size_t lookup_linear(BFOC_size_t i) const
  {
    int scan;
    for (scan = 0; scan < priv_size; scan++) {
      if (i == priv_key[scan])
	break;
    }
    return i;
  };

  inline BFOC_size_t lookup_interpol(BFOC_size_t i) const
  {
    int i_a = 0;
    int i_b = priv_size-1;
    while (priv_key[i_a] <= i && i <= priv_key[i_b]) {
      int center = i_a + (i-priv_key[i_a]) * (i_b-i_a) / (priv_key[i_b]-priv_key[i_a]);

      if (priv_key[center] < i)
	i_a = center + 1;
      else if (priv_key[center] > i)
	i_b = center - 1;
      else
	return center;
    }
    if (priv_key[i_a] == i)
      return i_a;
    else
      return priv_size;
  }

  // non-deferred == version
  // seems marginally faster than the deferred version ([1-1.5]/60 less time in some tests, in some machine...)
  inline BFOC_size_t lookup(BFOC_size_t i) const
  {
    int i_a = 0;
    int i_b = priv_size-1;
    while (i_b >= i_a) {
      int center = (i_a + i_b) / 2;
      if (priv_key[center] < i)
	i_a = center + 1;
      else if (priv_key[center] > i)
	i_b = center - 1;
      else 
	return center;
    }
    return priv_size;
  }

  /*
  // deferred == version
  inline int lookup(size_t i) const
  {
    int i_a = 0;
    int i_b = priv_size-1;
    while (i_b >= i_a) {
      int center = (i_a + i_b) / 2;
      if (priv_key[center] < i)
	i_a = center + 1;
      else 
	i_b = center - 1;

      if (priv_key[center] == i)
	return center;
    }
    return priv_size;
  }
  */

  // this will be used in lines like '... = rowp[s]'
  //  and should also produce an l-value because of lines like loaddata:2552: 'vmat[y][x][s] -= loss;' 
  //      (and all the reads from vmat[x][y][s] where vmat/rowp cannot be made const without changing 
  //       too many things)  => need non-const version as public
  // This returns a value, lookup returns its index
  inline const float& operator[](BFOC_size_t i) const
  { // return priv_rowp[i]; 
    int c = lookup(i);
    if (c >= priv_size)
      return priv_empty_missing;
    else
      return priv_val[c];
  };

  inline float& operator[](BFOC_size_t i) 
  { //return priv_rowp[i]; 
    /*size_t c = lookup(i);
    priv_key[c] = i;
    return priv_val[c];*/

    int c = lookup(i);
    if (c >= priv_size)
      return priv_empty_missing;
    else
      return priv_val[c];

  };

  // note it doesn't check for overflows
  inline void insert(BFOC_size_t key, float val) {
    priv_key[priv_size] = key;
    priv_val[priv_size] = val;
    priv_size++;
  }

  // pseudo-iterator methods: 'for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))'
  // their constness is arguable, but it makes it easier to have const rowp's.
  inline BFOC_size_t first() const
  { priv_it = 0; if (priv_size > 0) return priv_key[priv_it]; else return 0;};

  inline BFOC_size_t overflow() const
  { return std::numeric_limits<size_t>::max(); };

  inline BFOC_size_t next(size_t i) const
  { 
    priv_it++; 
    if (priv_it < priv_size) return priv_key[priv_it];
    else return std::numeric_limits<size_t>::max();
  };

  // An assignment operator like this is needed for example in marginal_loss.cpp:~60
  Biodiv_Features_Occur_Container_compact& operator=(Biodiv_Features_Occur_Container_compact& rhs)
    { delete [] this->priv_key; this->priv_key = NULL; delete [] this->priv_val; this->priv_val = NULL;
      this->priv_key = rhs.priv_key; this->priv_val = rhs.priv_val; this->priv_size = rhs.priv_size; /*this->priv_capacity = rhs.priv_capacity;*/ return *this; }


  // assigns from a sequence of features with consecutive indices: 0, 1,... , map_cnt
  Biodiv_Features_Occur_Container_compact& assign_seq(float* rhs, BFOC_size_t size)
    { priv_key = this->plu_idx_vec;  priv_val = rhs;  priv_size = size;
      return *this; }

  Biodiv_Features_Occur_Container_compact& assign(float* rhs, BFOC_size_t size)
    { // shallow assign: 
      // delete [] this->priv_val; this->priv_val = NULL;
      // this->priv_val = rhs; this->priv_size = size; 
      // full/deep assign
      this->clear();
      priv_key = new BFOC_size_t[size];  //new int[size];    
      for (size_t i=0; i<size; i++) {
	priv_key[i] = i;
      }
      priv_val = new float[size];
      memcpy(priv_val, rhs, size*sizeof(float));
      priv_size = size;
      return *this; }

  static size_t allocated_occurence_size() 
  { return sizeof(typeof(*priv_key))+sizeof(typeof(*priv_val)); };

 private:
  // foribd copy constructor
  Biodiv_Features_Occur_Container_compact(const Biodiv_Features_Occur_Container_compact& rhs)
    {};

  // forbid assignment? Cannot
  /*  Biodiv_Features_Occur_Container_compact& operator=(Biodiv_Features_Occur_Container_compact& rhs)
    {};
  */

 public:
  // yes, a std::vector could be better, but how to implement assign() then?
  // TODO: changed from int to int16 for now...
  // Put it here as public to simplify read/write in the vmat dump functions
  BFOC_size_t* priv_key;
  float* priv_val;
  BFOC_size_t priv_size;

 private:

  // not used at the moment, let's save some bytes... 100s of MBs maybe...
  /* size_t priv_capacity; */

  // Note this is static: rowp iterators are not thread-safe!
  //  both seem to contribute a smallish improvement over their non-static 
  //  alternatives (1-2% less time per each).
  static BFOC_size_t priv_it;
  static float priv_empty_missing;

  static BFOC_size_t* plu_idx_vec; 

 public:
  static void init_plu(size_t size) {
    // vector of indices 0, 1, 2,... , size
    if (NULL == plu_idx_vec) {
      plu_idx_vec = new BFOC_size_t[size];
    for (size_t i=0; i<size; i++) {
      plu_idx_vec[i] = i;
    }
  }

  }
  static void free_plu() {
    delete [] plu_idx_vec;
    plu_idx_vec = NULL;    
  }
};

typedef Biodiv_Features_Occur_Container_compact Biodiv_Features_Occur_Container; 

#elif defined NONCOMPACT_VMAT_FLOATP

 #warning "using compact vmat - floatp NON-COMPACT version"

// the "naive" version, just like 'float*'
class Biodiv_Features_Occur_Container_floatp
{
 public:
 Biodiv_Features_Occur_Container_floatp(): priv_val(NULL)/*, priv_size(0)*/
    {};

 ~Biodiv_Features_Occur_Container_floatp()
   { priv_val = NULL; }; // nope: if (priv_rowp) delete [] priv_rowp; 

  inline size_t size() const
  { if (priv_val) return priv_size; else return 0; };

  /*
  // in _floatp no need for any distinction between size and capacity
  // (the whole capacity is allocated and hence is also the size of priv_val)
  inline size_t capacity()
  { if (priv_val) return priv_size; else return 0; };
  */

  inline operator bool() const
  { return NULL != priv_val; };

  /*
  inline void init_zero()
  {
    priv_val = NULL; // NULL 
    }*/

  // semantics similar to std::
  // But note it is a single-use method. It does not free mem if previously allocated
  void reserve(size_t size, size_t map_cnt)
  {
    priv_val = new float[size];
    // yes, it will be repeated quite a few times...
    priv_size = size;
  };

  void clear()
  { if (priv_val) delete [] priv_val; priv_val = NULL; /*priv_size = 0;*/ };

  inline void fill_empty(int loop_up)
  { 
    for(int i=0; i<loop_up; i++) {
      priv_val[i] = -1;
    }
  }

  inline void set_empty(int pos)
  { priv_val[pos] = -1; }

  inline const float& operator[](size_t i) const
  { return priv_val[i];  };

  inline float& operator[](size_t i) 
  { return priv_val[i]; }

  inline void insert(size_t key, float val) {
    priv_val[key] = val;
  }

  // pseudo-iterator methods: 'for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))'
  // their constness is arguable, but it makes it easier to have const rowp's.
  inline size_t first() const
  { priv_it = 0; return priv_it; };

  inline size_t overflow() const
  { return priv_size; };

  inline size_t next(size_t i) const
  { 
    return priv_it++;
  };


  // NEEDS AN ASSIGNMENT OPERATOR: marginal_loss.cpp:61
  Biodiv_Features_Occur_Container_floatp& operator=(Biodiv_Features_Occur_Container_floatp& rhs)
    { delete [] priv_val;
      this->priv_val = rhs.priv_val; /*this->priv_size = rhs.priv_size;*/
      return *this; }

  Biodiv_Features_Occur_Container_floatp& assign(float* rhs, size_t size)
    { delete [] priv_val; 
      this->priv_val = rhs; /*this->priv_size = size;*/
      return *this; 
    }


 private:
  // foribd copy constructor
  Biodiv_Features_Occur_Container_floatp(const Biodiv_Features_Occur_Container_floatp& rhs)
    {};

  // forbid assignment? Cannot
  /*  Biodiv_Features_Occur_Container_floatp& operator=(Biodiv_Features_Occur_Container_floatp& rhs)
    {};
  */

 private:
  float* priv_val;

  // in the _floatp version, priv_size is the same for every container in vmat
  static size_t priv_size;

  static size_t priv_it;
  static float priv_empty_missing;
};

typedef Biodiv_Features_Occur_Container_floatp Biodiv_Features_Occur_Container; 

#endif // COMPACT_VMAT_LOOKUP elif NONCOMPACT_VMAT_FLOATP

#endif // OCCUR_CONTAINER_H
