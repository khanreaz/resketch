#include <lsqcpp.h>
#include <cmath>
#include <cstdlib>
#include <cstdarg>
#include <string>
#include <syslog.h>

#include "parameter_estimator.h"
#include "print.h"

#if 0
#include "csi_log.cpp"
#endif

static constexpr double SUBCARRIER_SPACING = 0.3125;

#ifdef __cplusplus
extern "C" {
double csi_fit(double observed_x, const Eigen::VectorXd& parameters)
{
    if(parameters.size() == 5) {
        auto gain_mismatch = parameters[0];
        auto timing_offset = parameters[1];
        auto phase_mismatch = parameters[2];
        auto tof_delay = parameters[3];
        auto phase_offset_error = parameters[4];
        
        return atan(gain_mismatch
            * sin(2 * M_PI * timing_offset * SUBCARRIER_SPACING * observed_x + phase_mismatch)
            / cos(2 * M_PI * timing_offset * SUBCARRIER_SPACING * observed_x))
            - 2 * M_PI * tof_delay * SUBCARRIER_SPACING * observed_x
            + phase_offset_error;
    }
    else {
        // XXX: format-security...
        std::string str = "csi_fit expects 5 parameters but number of parameters is " + std::to_string(parameters.size()) + "!\n";
        PRINT(LOG_ERR, "%s", str.c_str());
        std::exit(EXIT_FAILURE);
    }
}
}
#endif

/**
 * Convert csi phase array to vector of normalized csi phases.
 * 
 * \param observed_x Vector to contain normalized phase values.
 * \param observed_y Vector to contain subcarrier indices.
 * \param test_data Pointer to C array containing csi phase values.
 */
 #ifdef __cplusplus
extern "C" {
void csi_function(ErrorFunction::Vector& observed_x, ErrorFunction::Vector& observed_y, double *test_data)
{
    int n_observations = std::min(observed_x.size(), observed_y.size());
    int subcarrier_index = -(n_observations / 2);
    double mean = 0;
    
    for(int i = 0; i < n_observations; ++i) {
        mean += test_data[i];
    }
    mean /= n_observations;
    
    // TODO: this could be a constant.
    for(int i = 0; i < n_observations; ++i) {
        observed_y(i) = test_data[i] - mean;
        observed_x(i) = subcarrier_index;
        
        // skip 0
        if(++subcarrier_index == 0) {
            ++subcarrier_index;
        }
    }
}
}
#endif

/**
 * Estimate gain mismatch, timing offset, phase mismatch, time of flight and
 * phase offset error parameters for given csi phase values.
 * 
 * \param test_data Pointer to array containing csi phase data.
 * \param parameters Pointer to array to contain the parameters.
 */
#ifdef __cplusplus
extern "C" {
void estimate_csi_parameters(double *test_data, double *parameters)
{
    int n_observations = 56;
    ErrorFunction::Vector observed_x(n_observations);
    ErrorFunction::Vector observed_y(n_observations);
    
    csi_function(observed_x, observed_y, test_data);
    
    // Use least-squares-cpp library to fit function to observations.
    ErrorFunction fitting_error(observed_x, observed_y, &csi_fit);
    
    //~ lsq::GaussNewton<double, ErrorFunction, lsq::ArmijoBacktracking<double>> optimizer;
    //~ lsq::LevenbergMarquardt<double, ErrorFunction, lsq::ArmijoBacktracking<double>, lsq::NoCallback<double>, lsq::ForwardDifferences<double>> optimizer;
    lsq::LevenbergMarquardt<double, ErrorFunction, lsq::WolfeBacktracking<double>, lsq::NoCallback<double>, lsq::CentralDifferences<double>> optimizer;
    optimizer.setErrorFunction(fitting_error);
    optimizer.setMaxIterations(400);
    //~ optimizer.setMinGradientLength(1e-6);
    //~ optimizer.setMinStepLength(1e-6);
    optimizer.setMinError(1.49012e-8);
    //~ optimizer.setMinError(0);
    //~ optimizer.setStepSize(lsq::ArmijoBacktracking<double>(0.8, 0.1, 1e-10, 1.0, 0));
    optimizer.setVerbosity(0);

    // Initial guess for the unknown parameters of the fitting function
    Eigen::VectorXd initialGuess(5);
    //~ initialGuess << 0.512, -0.02857143, -0.006355, -0.02762, 0.1326;
    initialGuess << 0.512, -0.02812, -0.006355, -0.02762, 0.1326;

    auto result = optimizer.minimize(initialGuess);
    if(parameters != NULL) {
        parameters[0] = result.xval.transpose()(0);         // gain_mismatch
        parameters[1] = result.xval.transpose()(1);         // timing_offset
        parameters[2] = result.xval.transpose()(2);         // phase_mismatch
        parameters[3] = result.xval.transpose()(3);         // tof_delay
        parameters[4] = result.xval.transpose()(4);         // phase_offset_error
    }
}
}
#endif
