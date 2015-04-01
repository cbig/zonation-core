#ifndef ADMUS_H
#define ADMUS_H

int alloc_ADMU_data();
int free_ADMU_data();
int ADMUs_output_1();
int ADMUs_output_2(float remains);
int ADMUs_output_per_admu_curves_init();
int ADMUs_output_per_admu_curves_iter(float prop_lost);
int ADMUs_output_per_admu_grp_curves_init(Tgroups_info& groups_info);
int ADMUs_output_per_admu_grp_curves_iter(float prop_lost, Tgroups_info& groups_info);

int get_ADMU_data();
int calculate_and_output_ADMU_feature_weights();

int calculate_ADMU_redistributed_rank(int num);
int output_ADMU_redistributed_rank(int num);

#endif /* ADMUS_H */
