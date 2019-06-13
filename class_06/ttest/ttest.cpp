#include <iostream>
#include <iomanip>
#include "boost/math/distributions/students_t.hpp"


void two_samples_t_test_equal_sd(
        double Sm1, // Sm1 = Sample Mean 1.
        double Sd1,   // Sd1 = Sample Standard Deviation 1.
        unsigned Sn1,   // Sn1 = Sample Size 1.
        double Sm2,   // Sm2 = Sample Mean 2.
        double Sd2,   // Sd2 = Sample Standard Deviation 2.
        unsigned Sn2,   // Sn2 = Sample Size 2.
        double alpha)   // alpha = Significance Level.
{
    // A Students t test applied to two sets of data.
    // We are testing the null hypothesis that the two
    // samples have the same mean and that any difference
    // if due to chance.
    // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
    //

    // Print header:
    std::cout <<
              "_______________________________________________\n"
              "Student t test for two samples (equal variances)\n"
              "_______________________________________________\n\n";
    std::cout << std::setprecision(5);
    std::cout << std::setw(55) << std::left << "Number of Observations (Sample 1)" << "=  " << Sn1 << "\n";
    std::cout << std::setw(55) << std::left << "Sample 1 Mean" << "=  " << Sm1 << "\n";
    std::cout << std::setw(55) << std::left << "Sample 1 Standard Deviation" << "=  " << Sd1 << "\n";
    std::cout << std::setw(55) << std::left << "Number of Observations (Sample 2)" << "=  " << Sn2 << "\n";
    std::cout << std::setw(55) << std::left << "Sample 2 Mean" << "=  " << Sm2 << "\n";
    std::cout << std::setw(55) << std::left << "Sample 2 Standard Deviation" << "=  " << Sd2 << "\n";
    //
    // Now we can calculate and output some stats:
    //
    // Degrees of freedom:
    double v = Sn1 + Sn2 - 2;
    std::cout << std::setw(55) << std::left << "Degrees of Freedom" << "=  " << v << "\n";
    // Pooled variance and hence standard deviation:
    double sp = sqrt(((Sn1 - 1) * Sd1 * Sd1 + (Sn2 - 1) * Sd2 * Sd2) / v);
    std::cout << std::setw(55) << std::left << "Pooled Standard Deviation" << "=  " << sp << "\n";
    // t-statistic:
    double t_stat = (Sm1 - Sm2) / (sp * sqrt(1.0 / Sn1 + 1.0 / Sn2));
    std::cout << std::setw(55) << std::left << "T Statistic" << "=  " << t_stat << "\n";
    //
    // Define our distribution, and get the probability:
    //
    boost::math::students_t dist(v);
    double q = cdf(complement(dist, fabs(t_stat)));
    std::cout << std::setw(55) << std::left << "Probability that difference is due to chance" << "=  "
              << std::setprecision(3) << std::scientific << 2 * q << "\n\n";
    //
    // Finally print out results of alternative hypothesis:
    //
    std::cout << std::setw(55) << std::left <<
              "Results for Alternative Hypothesis and alpha" << "=  "
              << std::setprecision(4) << std::fixed << alpha << "\n\n";
    std::cout << "Alternative Hypothesis              Conclusion\n";
    std::cout << "Sample 1 Mean != Sample 2 Mean       ";
    if (q < alpha / 2)
        std::cout << "NOT REJECTED\n";
    else
        std::cout << "REJECTED\n";
    std::cout << "Sample 1 Mean <  Sample 2 Mean       ";
    if (cdf(dist, t_stat) < alpha)
        std::cout << "NOT REJECTED\n";
    else
        std::cout << "REJECTED\n";
    std::cout << "Sample 1 Mean >  Sample 2 Mean       ";
    if (cdf(complement(dist, t_stat)) < alpha)
        std::cout << "NOT REJECTED\n";
    else
        std::cout << "REJECTED\n";
    std::cout << std::endl << std::endl;
}

void two_samples_t_test_unequal_sd(
        double Sm1,   // Sm1 = Sample Mean 1.
        double Sd1,   // Sd1 = Sample Standard Deviation 1.
        unsigned Sn1,   // Sn1 = Sample Size 1.
        double Sm2,   // Sm2 = Sample Mean 2.
        double Sd2,   // Sd2 = Sample Standard Deviation 2.
        unsigned Sn2,   // Sn2 = Sample Size 2.
        double alpha)   // alpha = Significance Level.
{
    // A Students t test applied to two sets of data.
    // We are testing the null hypothesis that the two
    // samples have the same mean and
    // that any difference is due to chance.
    // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
    //
    using namespace std;

    // Print header:
    std::cout <<
              "_________________________________________________\n"
              "Student t test for two samples (unequal variances)\n"
              "_________________________________________________\n\n";
    std::cout << std::setprecision(5);
    std::cout << std::setw(55) << std::left << "Number of Observations (Sample 1)" << "=  " << Sn1 << "\n";
    std::cout << std::setw(55) << std::left << "Sample 1 Mean" << "=  " << Sm1 << "\n";
    std::cout << std::setw(55) << std::left << "Sample 1 Standard Deviation" << "=  " << Sd1 << "\n";
    std::cout << std::setw(55) << std::left << "Number of Observations (Sample 2)" << "=  " << Sn2 << "\n";
    std::cout << std::setw(55) << std::left << "Sample 2 Mean" << "=  " << Sm2 << "\n";
    std::cout << std::setw(55) << std::left << "Sample 2 Standard Deviation" << "=  " << Sd2 << "\n";
    //
    // Now we can calculate and output some stats:
    //
    // Degrees of freedom:
    double v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;
    v *= v;
    double t1 = Sd1 * Sd1 / Sn1;
    t1 *= t1;
    t1 /= (Sn1 - 1);
    double t2 = Sd2 * Sd2 / Sn2;
    t2 *= t2;
    t2 /= (Sn2 - 1);
    v /= (t1 + t2);
    std::cout << std::setw(55) << std::left << "Degrees of Freedom" << "=  " << v << "\n";
    // t-statistic:
    double t_stat = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);
    std::cout << std::setw(55) << std::left << "T Statistic" << "=  " << t_stat << "\n";
    //
    // Define our distribution, and get the probability:
    //
    boost::math::students_t dist(v);
    double q = cdf(complement(dist, fabs(t_stat)));
    std::cout << std::setw(55) << std::left << "Probability that difference is due to chance" << "=  " << std::setprecision(3) << std::scientific << 2 * q << "\n\n";


    // Calculates what T value such that the integral from [-inf,T] == 2.5%
    // This means we have 2.5% of the probability of reaching a number outside [-inf,T]
    // Because the distribution is symmetrycal, also 5% probability of a number outside [-T,T]
    double t = quantile(dist, 0.025); // -1.977e+00
    std::cout << std::setw(55) << std::left << "T value that would be enough for a 5% probability that it is due to chance (quantile) " << "=  " << std::setprecision(3) << std::scientific << t << "\n\n";
    //
    // Finally print out results of alternative hypothesis:
    //
    std::cout << std::setw(55) << std::left <<
              "Results for Alternative Hypothesis and alpha" << "=  "
              << std::setprecision(4) << std::fixed << alpha << "\n\n";
    std::cout << "Alternative Hypothesis              Conclusion\n";
    std::cout << "Sample 1 Mean != Sample 2 Mean       ";
    if (q < alpha / 2)
        std::cout << "NOT REJECTED\n";
    else
        std::cout << "REJECTED\n";
    std::cout << "Sample 1 Mean <  Sample 2 Mean       ";
    if (cdf(dist, t_stat) < alpha)
        std::cout << "NOT REJECTED\n";
    else
        std::cout << "REJECTED\n";
    std::cout << "Sample 1 Mean >  Sample 2 Mean       ";
    if (cdf(complement(dist, t_stat)) < alpha)
        std::cout << "NOT REJECTED\n";
    else
        std::cout << "REJECTED\n";
    std::cout << std::endl << std::endl;
}

int main() {
    //
    // Run tests for Car Mileage sample data
    // http://www.itl.nist.gov/div898/handbook/eda/section3/eda3531.htm
    // from the NIST website http://www.itl.nist.gov.  The data compares
    // miles per gallon of US cars with miles per gallon of Japanese cars.
    //
    two_samples_t_test_equal_sd(20.14458, 6.414700, 249, 30.48101, 6.107710, 79, 0.05);
    two_samples_t_test_unequal_sd(20.14458, 6.414700, 249, 30.48101, 6.107710, 79, 0.05);

    return 0;
}

/*
Output is:

  _______________________________________________
  Student t test for two samples (equal variances)
  _______________________________________________
  
  Number of Observations (Sample 1)                      =  249
  Sample 1 Mean                                          =  20.145
  Sample 1 Standard Deviation                            =  6.4147
  Number of Observations (Sample 2)                      =  79
  Sample 2 Mean                                          =  30.481
  Sample 2 Standard Deviation                            =  6.1077
  Degrees of Freedom                                     =  326
  Pooled Standard Deviation                              =  6.3426
  T Statistic                                            =  -12.621
  Probability that difference is due to chance           =  5.273e-030
  
  Results for Alternative Hypothesis and alpha           =  0.0500
  
  Alternative Hypothesis              Conclusion
  Sample 1 Mean != Sample 2 Mean       NOT REJECTED
  Sample 1 Mean <  Sample 2 Mean       NOT REJECTED
  Sample 1 Mean >  Sample 2 Mean       REJECTED
  
  
  _________________________________________________
  Student t test for two samples (unequal variances)
  _________________________________________________
  
  Number of Observations (Sample 1)                      =  249
  Sample 1 Mean                                          =  20.14458
  Sample 1 Standard Deviation                            =  6.41470
  Number of Observations (Sample 2)                      =  79
  Sample 2 Mean                                          =  30.48101
  Sample 2 Standard Deviation                            =  6.10771
  Degrees of Freedom                                     =  136.87499
  T Statistic                                            =  -12.94627
  Probability that difference is due to chance           =  1.571e-025
  
  Results for Alternative Hypothesis and alpha           =  0.0500
  
  Alternative Hypothesis              Conclusion
  Sample 1 Mean != Sample 2 Mean       NOT REJECTED
  Sample 1 Mean <  Sample 2 Mean       NOT REJECTED
  Sample 1 Mean >  Sample 2 Mean       REJECTED
  


*/