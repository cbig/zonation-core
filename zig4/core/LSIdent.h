#ifndef LSIDENT_H
#define LSIDENT_H

extern float   top_percent, min_percent, max_dist, min_simil;
extern String  LSIRfn, LSIDfn, LSImask;

extern int    spot_cnt, nwc, ccnt;
extern struct spot spots[];
extern int    **cm, **nwm, nwns_cnt, s_in_nw[];
extern bool   nw_in_min[];
extern float  Ebd[];

extern class  GridMap LSI_maskmap;

int LSIdent(int LSI_mode);
void LSCAnalysis(float f1, float f2, const String& cfn, const String& comp_outfn);

#endif // LSIDENT_H
