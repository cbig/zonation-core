#ifndef CONNECTED_COMPONENTS_H
#define CONNECTED_COMPONENTS_H

#include <vector>
#include <cstddef>

class CC_Labeling
{
 public:
 CC_Labeling(int soft_max_labels, int xdim, int ydim): labels(soft_max_labels), highest_label(0), last_lbl_equiv(0), xdim(xdim), ydim(ydim) {
    clear();
  }

  int do_cc_labeling_2pass(char**& sts, int**& lbl_layer);

 private:
  void clear() {
    std::fill(labels.begin(), labels.end(), Similarity());
    highest_label = 0;
  }

  inline bool is_root_label(int id) {
    return (labels[id].sameas == id);
  }

  inline int root_of(int id) {
    while (!is_root_label(id)) {
      // link this node to its parent's parent, just to shorten
      // the tree.
      labels[id].sameas = labels[labels[id].sameas].sameas;
      
      id = labels[id].sameas;
    }
    return id;
  }

  inline bool is_equivalent(int id, int as) {
    return (root_of(id) == root_of(as));
  }

  inline bool merge(int id1, int id2) {
    if(!is_equivalent(id1, id2)) {
      labels[root_of(id1)].sameas = root_of(id2);
      return false;
    }
    return true;
  }

  inline int new_label() {
    if(highest_label+1 > labels.size())
      labels.reserve(highest_label*2);
    labels.resize(highest_label+1);
    labels[highest_label] = Similarity(highest_label);
    return highest_label++;
  }

 private:

  struct Similarity {
  Similarity() : id(0), sameas(0) {}
  Similarity(int _id, int _sameas) : id(_id), sameas(_sameas) {}
  Similarity(int _id) : id(_id), sameas(_id) {}
    int id, sameas, tag;
  };


 public:
  std::vector<Similarity> labels;
  
 private:
  size_t last_lbl_equiv; 
  int highest_label;
  size_t xdim, ydim;
};

class CC_Labeling_Chang_Chen_Lu
{
 public:
  CC_Labeling_Chang_Chen_Lu(int xdim, int ydim);

  int do_cc_labeling_Chang_Chen_Lu(char**& sts, int**& lbl_layer, bool gen_contour);

 private:
  void tracer(/*char**& contour_sts, int**& lbl_layer, */ int& cy, int& cx, int& tracingdirection);
  void contour_tracing(/*char**& sts, int**& lbl_layer, */ int cy, int cx, int label_idx, int tracingdirection, bool gen_contour);
  size_t xdim, ydim;

  char** contour_sts;
  int** contour_lbl_layer;
};

class Grid_Connected_Components
{
 public:
  Grid_Connected_Components(size_t xdim, size_t ydim);
  ~Grid_Connected_Components();

  int count_labels(char**& status, int**& lbl_layer);

  int calc_count()
  { return n_cc_lbl_calc; }

 private:
    // for typical connected components labelling
  size_t n_cc_lbl_calc;
  size_t xdim, ydim;
};

#endif //  CONNECTED_COMPONENTS_H

