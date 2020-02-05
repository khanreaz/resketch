#include "lib/least-squares-cpp/include/lsqcpp.h"
#include <cmath>
#include <cstdlib>
#include <cstdarg>
#include <string>

#include <sstream>
#include <syslog.h>

#include "parameter_estimator.h"
#include "lib/print.h"

//#if 0
//#include "alice_bob_phase.c"
//#endif

static constexpr double SUBCARRIER_SPACING = 0.3125;

/**
 * Error function that calculates residuals.
 */
struct ErrorFunction
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

    ErrorFunction() : observed_x(0), observed_y(0), fitting_function(0) {}
    ErrorFunction(const Vector& observed_x,
                const Vector& observed_y,
                double (*fitting_function)(double, const Eigen::VectorXd&)) : observed_x(observed_x), observed_y(observed_y), fitting_function(fitting_function) {}

    void operator()(const Eigen::VectorXd& xval,
                    Eigen::VectorXd& fval,
                    Eigen::MatrixXd&) const
    {
        fval.resize(std::min(observed_x.size(), observed_y.size()));
        for(lsq::Index i = 0; i < std::min(observed_x.size(), observed_y.size()); ++i) {
            //~ fval(i) = pow(observed_y(i) - fitting_function(observed_x(i), xval), 2);
            fval(i) = observed_y(i) - fitting_function(observed_x(i), xval);
        }
    }

    /**
     * Sets data points for non-linear least squares fitting.
     *
     * @param observed_x Reference to the vector containing data from independent variable.
     * @param observed_y Reference to the vector containing data from dependent variable.
     */
    void set_observation(const Vector& observed_x, const Vector& observed_y)
    {
        this->observed_x = observed_x;
        this->observed_y = observed_y;
    }

    /**
     * Set model function and number of parametes to find.
     *
     * @param fitting_function Function pointer to the model function whose parameters are to be found.
     *        The first argument is the x variable, the second is a Vector containing the current guess
     *        for the parameters.
     */
    void set_fitting_function(double (*fitting_function)(double, const Eigen::VectorXd&))
    {
        this->fitting_function = fitting_function;
    }

private:
    Vector observed_x;
    Vector observed_y;
    double (*fitting_function)(double, const Eigen::VectorXd&);
};

double oscillation_fit(double observed_x, const Eigen::VectorXd& parameters)
{
    if(parameters.size() == 3) {
        return parameters(0) * std::sin(parameters(1) * observed_x + parameters(2));
    }
    else {
        std::cerr << "oscillation_fit expects 3 parameters but number of parameters is " << parameters.size() << "!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

/**
 * This function is for testing purposes.
 * It generates artificial observations of a sinus function with variable amplitude, frequency and phase.
 */
// void oscillation_function(ErrorFunction::Vector& observed_x,
//                           ErrorFunction::Vector& observed_y,
//                           double amplitude,
//                           double frequency,
//                           double phase)
// {
//     int n_observations = std::min(observed_x.size(), observed_y.size());
//     double x_step = 2 * M_PI / frequency / n_observations;

//     if(n_observations > 0) {
//         observed_x(0) = 0;
//         for(int i = 0; i < n_observations; ++i) {
//             observed_y(i) = amplitude * std::sin(frequency * observed_x(i) + phase);
//             if(i < n_observations - 1) {
//                 observed_x(i + 1) += observed_x(i) + x_step;
//             }
//         }
//     }
// }

// void oscillation_test()
// {
//     // Create artificial observations for an oscillation
//     double amplitude = 1.123456789;
//     double frequency = 1.89;
//     double phase = M_PI / 2.5;
//     int n_observations = 56;
//     ErrorFunction::Vector errors(n_observations);
//     ErrorFunction::Vector observed_x(n_observations);
//     ErrorFunction::Vector observed_y(n_observations);

//     oscillation_function(observed_x, observed_y, amplitude, frequency, phase);
//     std::cout << "Generated " << n_observations << " observations." << std::endl;
//     std::cout << "Amplitude: " << amplitude << ", Frequency: " << frequency << ", Phase: " << phase << std::endl;

//     // Use least-squares-cpp library to fit function to observations.
//     ErrorFunction fitting_error(observed_x, observed_y, &oscillation_fit);

//     lsq::GaussNewton<double, ErrorFunction, lsq::ArmijoBacktracking<double>> optimizer;
//     optimizer.setErrorFunction(fitting_error);
//     optimizer.setMaxIterations(200);
//     //~ optimizer.setMinGradientLength(1e-6);
//     //~ optimizer.setMinStepLength(1e-6);
//     //~ optimizer.setMinError(0);
//     optimizer.setStepSize(lsq::ArmijoBacktracking<double>(0.8, 0.1, 1e-10, 1.0, 0));
//     optimizer.setVerbosity(0);

//     // Initial guess for the unknown parameters of the fitting function
//     Eigen::VectorXd initialGuess(3);
//     initialGuess << 1.0, 1.5, 1.3;

//     auto result = optimizer.minimize(initialGuess);
//     std::cout << "Parameters: " << result.xval.transpose() << std::endl;
// }

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

void csi_function(ErrorFunction::Vector& observed_x, ErrorFunction::Vector& observed_y, double *test_data)
{
    int n_observations = std::min(observed_x.size(), observed_y.size());
    int subcarrier_index = -28;
    double mean = 0;

    for(int i = 0; i < n_observations; ++i) {
        mean += test_data[i];
    }
    mean /= n_observations;

    for(int i = 0; i < n_observations; ++i) {
        observed_y(i) = test_data[i] - mean;
        observed_x(i) = subcarrier_index;

        // skip 0
        if(++subcarrier_index == 0) {
            ++subcarrier_index;
        }
    }
}

#ifdef __cplusplus
extern "C" {
#endif
void estimate_csi_parameters(double *test_data, double *return_data, int return_index = 0)
{
    int n_observations = 56;
    ErrorFunction::Vector observed_x(n_observations);
    ErrorFunction::Vector observed_y(n_observations);

    csi_function(observed_x, observed_y, test_data);

    // Use least-squares-cpp library to fit function to observations.
    ErrorFunction fitting_error(observed_x, observed_y, &csi_fit);

    //~ lsq::GaussNewton<double, ErrorFunction, lsq::ArmijoBacktracking<double>> optimizer;
    //~ lsq::LevenbergMarquardt<double, ErrorFunction, lsq::ArmijoBacktracking<double>> optimizer;
    //~ lsq::LevenbergMarquardt<double, ErrorFunction, lsq::ArmijoBacktracking<double>, lsq::NoCallback<double>, lsq::ForwardDifferences<double>> optimizer;
    lsq::LevenbergMarquardt<double, ErrorFunction, lsq::WolfeBacktracking<double>, lsq::NoCallback<double>, lsq::CentralDifferences<double>> optimizer;
    optimizer.setErrorFunction(fitting_error);
    optimizer.setMaxIterations(100);
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

    /* to_string has fixed default precision of 6 digits, if higher precision is needed ostringstring can be used. */
    std::ostringstream gain_mismatch_sstream;
    gain_mismatch_sstream << std::fixed;
    gain_mismatch_sstream << std::setprecision(8);
    gain_mismatch_sstream << result.xval.transpose()(0);

    // XXX: format-security...
    // std::string str = "Parameters: " + gain_mismatch_sstream.str() + ","
    //                                  + std::to_string(result.xval.transpose()(1)) + ","
    //                                  + std::to_string(result.xval.transpose()(2)) + ","
    //                                  + std::to_string(result.xval.transpose()(3)) + ","
    //                                  + std::to_string(result.xval.transpose()(4)) + "\n";
    // PRINT(LOG_INFO, "%s", str.c_str());
    return_data[return_index + 0] = result.xval.transpose()(0);
    return_data[return_index + 1] = result.xval.transpose()(1);
    return_data[return_index + 2] = result.xval.transpose()(2);
    return_data[return_index + 3] = result.xval.transpose()(3);
    return_data[return_index + 4] = result.xval.transpose()(4);
}
#ifdef __cplusplus
}
#endif