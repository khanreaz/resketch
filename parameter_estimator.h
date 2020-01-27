#ifndef _PARAMETER_ESTIMATOR_H_
#define _PARAMETER_ESTIMATOR_H_

#ifdef __cplusplus
extern "C" void estimate_csi_parameters(double *test_data, double *return_data, int return_index);
#else
void estimate_csi_parameters(double *test_data, double *return_data, int return_index);
#endif

#endif
