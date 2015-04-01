#ifndef MARGINAL_LOSS_H
#define MARGINAL_LOSS_H

// local xxxPLULA variables
extern float  PLULA_data_vec[], PLULA_SSI_ps[];
extern int    PLULA_SSI_spp[], PLULA_SSI_cnt;

float delta_value_new(int x, int y, bool single_cell_mode);

#endif /* MARGINAL_LOSS_H */
