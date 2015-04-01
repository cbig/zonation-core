#include "Unit1.h"
#include "GridMap.h"
#include "defines.h"
#include "marginal_loss.h"
#include "PLULA.h"
#include "LoadData.h"
#include "bat_run.h"
#include "grid_utils.h"
#include <cstdio>
#include <cstring>
#include <limits>

#define CAZ_ALTERNATING_EPS

// local xxxPLULA variables
float  PLULA_data_vec[MAX_SPP_COUNT], PLULA_SSI_ps[MAX_SPP_COUNT];
int    PLULA_SSI_spp[MAX_SPP_COUNT], PLULA_SSI_cnt;

float dv_initialize(int y, int x, bool single_cell_mode, Biodiv_Features_Occur_Container& rowp,
		    int &adu, float &Wjadu, int &pln, float &PLULA_cost, bool &use_PLULA_mode)
{
	int a, s;

	if (!single_cell_mode) // this is for removal, needs to evaluate single cell with BQP
		use_PLULA_mode = use_PLULA;
	else
		use_PLULA_mode = false;

	if (use_BQP || use_tree_conn)   //  Only when nbh loss, does this matrix need to be zeroed.
	{								//  Else just one row updated direct.
		for(a=0;a<ADM_set.count;++a)
		{
		        // for(s=0;s<map_cnt;++s)
                        for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s)) {
			  ADMUxspp_loss[a][s]=0.0;
			}
			loss_in_ADU[a]=false;
		}
	}

	if (ADM_set.use_ADMUs)
	{
		adu = ADMUs[y][x];
		if (adu<0 || adu>=ADM_set.count)
		  adu=0;
	}
	else
		Wjadu = 1.0;

	if (use_PLULA_mode)
	{
		pln = PLL[y][x];
		if (pln<0)
		{
			ShowMessage("GDV PLULA number error");
			return 1.0;
		}
		if (PLvec[pln].checked)  // all cells in the PLULA return same delta value
			return PLvec[pln].dv;
	}

	if (!use_PLULA_mode)  // xxx note: ADMUs work with PLUs because of rowp is fine and indep of ADM weights!
		rowp = vmat[y][x];
	else
	{
		if (use_SSI || (!PLvec[pln].allocated)) // xxx storing of SSI data with
		{                                     // plula has not been implemented
			PLULA_cost = get_PLULA_info(pln, PLULA_data_vec,
						    PLULA_SSI_spp, PLULA_SSI_ps, PLULA_SSI_cnt);
			// rowp = PLULA_data_vec;  // COMPACT_VMAT
			rowp.assign_seq(PLULA_data_vec, map_cnt);
		}
		else // spdata has been precomp + stored (but not yet SSI data)
		{
			//Form1->Memo1->Lines->Add("heer");
			// rowp = PLvec[pln].datavec;  // COMPACT_VMAT
		        rowp.assign_seq(PLvec[pln].datavec, map_cnt);
			PLULA_cost = PLvec[pln].cost;
		}
	}

	return -1;
}

void add_SSI_value(int y, int x, float &delta_value, bool use_PLULA_mode)
{
	float  SSIval, newval, oldval;
	int    s, SSIspnum, SSIpos, end;

	// xxxMCZ does not prevent wj<0 for SSI spp. SSI.wj<0 should work ok. 5/09
	if (use_PLULA_mode)
		end = SSI_spp_cnt;
	else
	{
		end    = SSIxyCnt[y][x];
		SSIpos = SSIxyPos[y][x];
	}

	for(s=0; s<end; s++)
	{

		if (use_PLULA_mode)  // note end = SSI_spp_cnt
		{
			SSIspnum = s;
			SSIval   = PLULA_SSI_ps[s];
			if (SSIval<=0.0)
				continue;
		}
		else
		{
			//          if ((SSIpos<0) || (SSIpos>SSI_lpos))
			//            ShowMessage("SSIerror");
			SSIspnum  = SSI_list[SSIpos+s].spnum;
			SSIval    = SSI_list[SSIpos+s].val;
		}

		if (removal_rule==1) // original core area zonation
		{
			//              Form1->Memo1->Lines->Add("predeltassival="+FloatToStr(delta_value));
			newval = SSIval*SSI[SSIspnum].weight/SSI_repr_lvls[SSIspnum];
			if (newval>delta_value)
				delta_value=newval;
			//              sprintf(txt, "(%i,%i) sp=%i, dv=%f %f", y,x,SSIspnum,delta_value,
			//                SSIval*SSI[SSIspnum].weight/SSI_repr_lvls[SSIspnum]);
			//              Form1->Memo1->Lines->Add(txt);
		}
		else if (removal_rule==2) // convex BF with par
		{
			//              if (SSI_repr_lvls[SSIspnum]>0.0)
			delta_value += SSI[SSIspnum].weight*        // xxx problem with pow
					( z_pow(SSI_repr_lvls[SSIspnum], SSI[SSIspnum].sr_par) -
					  z_pow(SSI_repr_lvls[SSIspnum]-SSIval, SSI[SSIspnum].sr_par)   );
			//                delta_value += SSI[SSIspnum].weight*        // xxx problem with pow
			//                  ( z_pow(SSI_repr_lvls[SSIspnum], 0.4) -
			//                    z_pow(SSI_repr_lvls[SSIspnum]-SSIval, 0.4)   );
		}
		else if (removal_rule==3) // target BF
		{
			if ((SSI_repr_lvls[SSIspnum]-SSIval)>SSI[SSIspnum].sr_par)
			{ // here derivative not ok because can be large step change in repr
				oldval = z_pow( (SSI_repr_lvls[SSIspnum]-SSI[SSIspnum].sr_par)
					      /(1.0-SSI[SSIspnum].sr_par), 0.25);
				newval = z_pow( (SSI_repr_lvls[SSIspnum]-SSIval-SSI[SSIspnum].sr_par)
					      /(1.0-SSI[SSIspnum].sr_par), 0.25);
				delta_value += (oldval-newval);
			}
			else
			{
				if (SSI_repr_lvls[SSIspnum]<SSI[SSIspnum].sr_par)
					delta_value += 0.0;
				else
					delta_value += (map_cnt+SSI_spp_cnt);
			}
		}
		else if (removal_rule==4)
		{
			delta_value += (   rr4(SSI_repr_lvls[SSIspnum],         SSI[SSIspnum].weight, SSI[SSIspnum].Tj, SSI[SSIspnum].w2, SSI[SSIspnum].exp1, SSI[SSIspnum].exp2)
					   - rr4(SSI_repr_lvls[SSIspnum]-SSIval,  SSI[SSIspnum].weight, SSI[SSIspnum].Tj, SSI[SSIspnum].w2, SSI[SSIspnum].exp1, SSI[SSIspnum].exp2) );
		}
	}  // SSI spp
}

inline float rr1(float w, float loss, float remaining)
{
	return w*(loss/remaining); // fraction of remaining
}

inline float rr2_delta(float w, float o, float n, float par)
{
#ifdef ZCORE_DEBUG
  if (isnan(f)) {
    Form1->Memo1->Lines->Add("rrw_delta(): o: %f, n: %f, par: %f",o, n, par);
  }
#endif
      
  return w*(z_pow(o,par)-z_pow(n,par));
}

inline float rr3_delta(float target, float o, float n)
{
	if ((o<target) || (o==n) || (target<=0.0))
		return 0.0;
	else if ((o>=target) && (n<target))
		return 1.0;
	else return ( z_pow((1.0-o)/(1.0-target),0.25) -z_pow((1.0-n)/(1.0-target),0.25) );
}

float add_regional_losses(int focal_adu, int x, int y)
{
	float delta_addition=0.0, old_val, loss, new_val, *arow, da2;
	int aduloop, s, adu_start, adu_end;

	if (use_BQP || use_tree_conn) // presently the only forms of conn that actively spread outside focal area
	{
		adu_start = 0;
		adu_end   = ADM_set.count;
	}
	else // only local loss from one area expected
	{
		adu_start = focal_adu;
		adu_end   = focal_adu + 1;
	}

	for(aduloop=adu_start; aduloop<adu_end; aduloop++)  // xxx change structure to (adu, sp) pairs?
	{
		if (loss_in_ADU[aduloop])
		{
			arow = &ADMUxspp_loss[aduloop][0];
			//for(s=0; s<map_cnt; s++)
			//{
			// Only need to care about the spp present in this cell
			Biodiv_Features_Occur_Container& rowp = vmat[y][x];
			for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))
			{
				loss =  arow[s];
				if (loss>0.0f && ADMUxSP_repr_orig[aduloop][s] >.0f) // avoid NaNs
				{

				  // Avoid negative remaining repr!
				  // Otherwise, for example in rr2_delta -> new_val<0
				  //           and 0<sr_par<1 => pow produces NaN
				  if (loss > ADMUxSP_repr[aduloop][s]) {
				    loss = ADMUxSP_repr[aduloop][s];
				    new_val = 0;
				  } else {
				    new_val = (ADMUxSP_repr[aduloop][s]-loss)/ADMUxSP_repr_orig[aduloop][s];
				  }
				  old_val = ADMUxSP_repr[aduloop][s]/ADMUxSP_repr_orig[aduloop][s];

					switch(removal_rule)
					{
					  int nb_count_shifted; float round_bit;

					case 1:			//
						da2 = rr1(ADM_combined_weights[s][aduloop], old_val-new_val, old_val);
#ifdef CAZ_ALTERNATING_EPS
						nb_count_shifted = get_nb8_count(status, x, y) - 3;    // nb8_count is in [0,7]
						round_bit = nb_count_shifted * da2 * std::numeric_limits<float>::epsilon();
						da2 += 4*round_bit;     // 4x because sensitivity to approx errors looks higher here (local loss)
#endif
						if (da2>delta_addition)
							delta_addition = da2;
						
						break;
					case 2:			//  addition is ok
						delta_addition += rr2_delta(ADM_combined_weights[s][aduloop], old_val, new_val, spp[s].sr_par);
						break;
					case 3:			//
						delta_addition += rr3_delta(ADM_combined_weights[s][aduloop], old_val, old_val-loss);
						break;
					case 4:			//
						delta_addition += (rr4(old_val, ADM_combined_weights[s][aduloop], spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2)
								   - rr4(new_val-loss,  ADM_combined_weights[s][aduloop], spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2) );
						break;
					};
				}
			}
		}
	}

	return delta_addition/ADM_set.count;
}

float delta_value_new(int x, int y, bool single_cell_mode)
{
        Biodiv_Features_Occur_Container rowp;
	float ret_val; // species representation values in the element
	float  delta_value, orig, nbv_val, dv_val, tmpf, SSIval, oldval, newval;
	int    pln, remaining, s;
	float  min_delta_value, PLULA_cost; // 10.12.2009; added for CAZ MCZ
	struct sp *psp;
	bool   use_PLULA_mode;
	// ADMUs
	float Wjadu=1.0;
	int adu;

	ret_val = dv_initialize(y, x, single_cell_mode, rowp, adu, Wjadu, pln, PLULA_cost, use_PLULA_mode);
	if (ret_val>=0.0)
		return ret_val; // plula value was available already; must then account for ADMUs

	if (removal_rule==1) // original core area zonation  // xxxMCZ Q resolved 10.12.2009
	{
		delta_value     = -std::numeric_limits<float>::max();
		min_delta_value = std::numeric_limits<float>::max(); // requires separate tracking for negative features
	}
	else
		delta_value=0.0;

#define ULTIMATE_RISKY_SPEEDUP_ADMU_MODE2_GLOBAL0 1
#ifdef ULTIMATE_RISKY_SPEEDUP_ADMU_MODE2_GLOBAL0
	// Beginning of block of calculations that, when in mode 2, are needed only with global weight !=0
	if (ADM_set.use_ADMUs && (ADM_set.ADMU_mode==2) && (0==ADM_set.mode2_global_weight) && !use_PLULA_mode && !use_BQP)
	{

	} else {
#endif
	//for(s=0; s<map_cnt; s++)
	for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))
	{
                float rowp_s_val = rowp[s];
		//if (!rowp[s]>=0.0)    <--- NaN trick
	        //if (rowp[s]<0.0)
		if (rowp_s_val<0.0)
	        	continue;
		if (ADM_set.use_ADMUs)
		{
			if (ADM_set.ADMU_mode==1)
				Wjadu = ADM_combined_weights[s][adu];
			else
				Wjadu = spp[s].weight;  // xxx changed for mode 2 oct 10 2010
		}
		psp = &spp[s];

		if ((removal_rule==3) && (repr_lvls[s]<psp->sr_par))
			continue; // spp under target already

		// first calculate loss, from appropriate combination of options
		if ((use_PLULA_mode) && (use_tree_conn)) // PLUs & trees
			tmpf = PLULA_tree_loss(pln, s);
		else if (!use_BQP || !psp->BQP_used) { // no BQP for spp (non-spatial)
		        // tmpf = rowp[s];
			tmpf = rowp_s_val;
		} else {
			if (use_PLULA_mode) // PLU & BQP
				tmpf = PLULA_BQP_loss(pln, s);
			else // just BQP
			{
				remaining = nbms[nbm_num[s]][y][x]; // number of nbs remaining
				orig      = nbms_orig[nbm_num[s]][y][x];
				if (orig>0)
				{
					nbv_val   = psp->BQP_link[(int)(10000*remaining/orig)]; // matches initial #nbs and init prob
					dv_val    = BQP_dv_buf(x,y,s); // dv_buf() should fill in
				}
				else
				{
					nbv_val   = 1.0;
					dv_val    = 0.0;
				}

				//tmpf      = nbv_val*rowp[s]+dv_val; // this is global loss; remaining locally + buffer loss
				tmpf      = nbv_val*rowp_s_val+dv_val; // this is global loss; remaining locally + buffer loss
				BQP_spp_loss_tmp[s]   = tmpf;       // needed in remove cell
				ADMUxspp_loss[adu][s] += nbv_val;    // = remaining locally, stored for regional losses!
				// stuff added already before in dv_buf
			}
		}

#ifdef ZCORE_DEBUG
		if (std::isnan(tmpf) )
		  Form1->Memo1->Lines->Add("delta_value_new(): --- tmpf: %f (%d) \t", tmpf,s);
#endif

		// then RR-specific aggregation of value; phase 1 - global value

		switch(removal_rule)
		{
		case 1: {

		  float this_val = tmpf*smult[s]*Wjadu;
#ifdef CAZ_ALTERNATING_EPS
		  int nb_count_shifted = get_nb8_count(status, x, y) - 3;    // nb8_count is in [0,7]
		  float round_bit = nb_count_shifted * this_val * std::numeric_limits<float>::epsilon();
		  this_val += round_bit;
#endif
		  delta_value     = max( delta_value, this_val);  // xxxMCZ; how combines with wj<0???
		  min_delta_value = min(min_delta_value, this_val);  // xxxMCZ; how combines with wj<0???
		}
		  break;
		case 2:
		  delta_value += tmpf*sp_derivative[s]*Wjadu;
		  break;

		case 3:
		  if ((repr_lvls[s]-tmpf)>psp->sr_par)
		    delta_value += tmpf*sp_derivative[s];
		  else {
		    if (repr_lvls[s]<psp->sr_par)
		      delta_value += 0.0;
		    else
		      delta_value += (map_cnt+SSI_spp_cnt);
		  }
		  break;

		case 4:
		  if (ADM_set.use_ADMUs)
		    delta_value += (   rr4(repr_lvls[s],       Wjadu, spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2)
				       - rr4(repr_lvls[s]-tmpf,  spp[s].weight, spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2) );
		  else
		    delta_value += (   rr4(repr_lvls[s],       spp[s].weight, spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2)
				       - rr4(repr_lvls[s]-tmpf,  spp[s].weight, spp[s].Tj, spp[s].w2, spp[s].exp1, spp[s].exp2) );
		  break;
		  // note: RR#5 = random; defined elsewhere does not call gdv
		};
		// note: cut loop if min_value exceeded (should be added for sure.!!!)
	} // spp

	if ((removal_rule==1) && (neg_weights_used))
	{
		if (min_delta_value<0.0)
			delta_value += min_delta_value; // MCZ 10.12.2009; mdv<0, thus addition ok.
	}										// missing 2nd variant; summing of neg effects
	
	if (use_SSI)
		add_SSI_value(y, x, delta_value, use_PLULA_mode);

#ifdef ULTIMATE_RISKY_SPEEDUP_ADMU_MODE2_GLOBAL0
	} // End of block of calculations that, when in mode 2, are needed only with global weight !=0
#endif

	if (!use_BQP && !use_tree_conn && ADM_set.use_ADMUs) // only loss in the focal adu, in cell or PLU
	{

	  // This init of all adus is very expensive!
	  //for(s=0; s<ADM_set.count; s++)
	  //		loss_in_ADU[s]=false;
	  // at the beginning all are false, change only this one here and set it false again after calling add_regional_losses(adu)...
	  loss_in_ADU[adu]=true; // set losses for one adu only;

	  // for(s=0; s<map_cnt; s++)
	  for(s = rowp.first(); s != rowp.overflow(); s = rowp.next(s))
	    ADMUxspp_loss[adu][s] = rowp[s];
	}

	if (ADM_set.use_ADMUs && (removal_rule!=5) && (ADM_set.ADMU_mode==2))
	  delta_value = ADM_set.mode2_global_weight*delta_value + (1.0-ADM_set.mode2_global_weight)*add_regional_losses(adu, x, y);
	// is low final dv cutting computations short later!!!? // xxx modified 4/2010

	if (!use_BQP && !use_tree_conn && ADM_set.use_ADMUs) {
	  loss_in_ADU[adu]=false;
	}

	if (use_cost) // ??? does cost neet to be local? No, local+global can be pre-combined in costmap?
		delta_value /= costmap.m[y][x];

	if (use_PLULA_mode)
	{
		delta_value /= PLULA_cost; // or implicitly sum of area, if cost not used
		PLvec[pln].checked = true;
		PLvec[pln].dv      = delta_value;
	}

	return delta_value;
}

