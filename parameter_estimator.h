#ifndef _PARAMETER_ESTIMATOR_H_
#define _PARAMETER_ESTIMATOR_H_

#ifdef __cplusplus
extern "C" void estimate_csi_parameters(double *test_data, double *parameters);
#else
void estimate_csi_parameters(double *test_data, double *parameters);
#endif

/**
 * Error function that calculates residuals.
 */
#ifdef __cplusplus
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
#endif
#endif
