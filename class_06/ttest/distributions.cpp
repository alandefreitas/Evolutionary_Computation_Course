#include <iostream>
#include <iomanip>
#include "boost/math/distributions/students_t.hpp"
#include "boost/math/distributions/beta.hpp"
#include "boost/math/distributions/binomial.hpp"
#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/fisher_f.hpp"
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/normal.hpp>

int main(){
    // let's make it easier to use boost
    using namespace boost;
    // a t-distribution with 10 degrees of freedom
    {
        math::students_t dist(10);
        // Returns PDF (density) at point x of distribution my_dist.
        // Probability of getting 0.0 in a t-distribution is:
        std::cout << "pdf(dist, 0.0) = " << pdf(dist, 0.0) << std::endl;
        // Returns CDF (integral from -infinity to point x) of distribution my_dist.
        // Probability of getting a number lesser than 0.0 in a t-distribution is:
        std::cout << "cdf(dist, 0.0) = " << cdf(dist, 0.0) << std::endl;
        // Probability of getting a number lesser than 1.0 in a t-distribution is:
        std::cout << "cdf(dist, 1.0) = " << cdf(dist, 1.0) << std::endl;
        // Probability of getting a number greater than 1.0 in a t-distribution is:
        std::cout << "1.0 - cdf(dist, 0.0) = " << 1.0 - cdf(dist, 0.0) << std::endl;
        std::cout << "1.0 - cdf(dist, 0.1) = " << 1.0 - cdf(dist, 0.1) << std::endl;
        // Removing the cancellation error with the complement distribution
        std::cout << "cdf(complement(dist, 0.0)) = " << cdf(complement(dist, 0.0)) << std::endl;
        std::cout << "cdf(complement(dist, 0.1)) = " << cdf(complement(dist, 0.1)) << std::endl;
        // Returns QUANTILE the value of the random variable x such that cdf(my_dist, x) == p.
        // This is the inverse of PDF
        // Note that lower critical values are just quantiles
        // Less than 50% of the time, we have a number smaller than:
        std::cout << "quantile(dist, 0.5)) = " << quantile(dist, 0.5) << std::endl;
        // Less than 30% of the time, we have a number smaller than:
        std::cout << "quantile(dist, 0.3)) = " << quantile(dist, 0.3) << std::endl;
        // Less than 70% of the time, we have a number smaller than:
        std::cout << "quantile(dist, 0.7)) = " << quantile(dist, 0.7) << std::endl;
        // Note that upper critical values are just quantiles
        // Less than 50% of the time, we have a number greater than:
        std::cout << "quantile(complement(dist, 0.5)) = " << quantile(complement(dist, 0.5)) << std::endl;
        // Less than 30% of the time, we have a number greater than:
        std::cout << "quantile(complement(dist, 0.3)) = " << quantile(complement(dist, 0.3)) << std::endl;
        // Less than 70% of the time, we have a number greater than:
        std::cout << "quantile(complement(dist, 0.7)) = " << quantile(complement(dist, 0.7)) << std::endl;
        // Parameter finders are implemented as static member functions of the distributions
        // How many more measurements would I have to take before I would get an X% probability that the difference is real?
        auto d = math::students_t::find_degrees_of_freedom(
                1.3,        // difference from true mean to detect
                0.05,       // maximum risk of falsely rejecting the null-hypothesis.
                0.1,        // maximum risk of falsely failing to reject the null-hypothesis.
                0.13);      // sample standard deviation
        std::cout << "d = " << d << std::endl;
    }
    // beta_distribution with a = 10 and b = 20
    {
        math::beta_distribution<> dist(10, 20);
        // Returns PDF (density) at point x of distribution my_dist.
        // Probability of getting 0.0 in a t-distribution is:
        std::cout << "pdf(dist, 0.0) = " << pdf(dist, 0.0) << std::endl;
        // Returns CDF (integral from -infinity to point x) of distribution my_dist.
        // Probability of getting a number lesser than 0.0 in a t-distribution is:
        std::cout << "cdf(dist, 0.0) = " << cdf(dist, 0.0) << std::endl;
        // Probability of getting a number lesser than 1.0 in a t-distribution is:
        std::cout << "cdf(dist, 1.0) = " << cdf(dist, 1.0) << std::endl;
        // Probability of getting a number greater than 1.0 in a t-distribution is:
        std::cout << "1.0 - cdf(dist, 0.0) = " << 1.0 - cdf(dist, 0.0) << std::endl;
        std::cout << "1.0 - cdf(dist, 0.1) = " << 1.0 - cdf(dist, 0.1) << std::endl;
        // Removing the cancellation error with the complement distribution
        std::cout << "cdf(complement(dist, 0.0)) = " << cdf(complement(dist, 0.0)) << std::endl;
        std::cout << "cdf(complement(dist, 0.1)) = " << cdf(complement(dist, 0.1)) << std::endl;
        // Returns QUANTILE the value of the random variable x such that cdf(my_dist, x) == p.
        // This is the inverse of PDF
        // Note that lower critical values are just quantiles
        // Less than 50% of the time, we have a number smaller than:
        std::cout << "quantile(dist, 0.5)) = " << quantile(dist, 0.5) << std::endl;
        // Less than 30% of the time, we have a number smaller than:
        std::cout << "quantile(dist, 0.3)) = " << quantile(dist, 0.3) << std::endl;
        // Less than 70% of the time, we have a number smaller than:
        std::cout << "quantile(dist, 0.7)) = " << quantile(dist, 0.7) << std::endl;
        // Note that upper critical values are just quantiles
        // Less than 50% of the time, we have a number greater than:
        std::cout << "quantile(complement(dist, 0.5)) = " << quantile(complement(dist, 0.5)) << std::endl;
        // Less than 30% of the time, we have a number greater than:
        std::cout << "quantile(complement(dist, 0.3)) = " << quantile(complement(dist, 0.3)) << std::endl;
        // Less than 70% of the time, we have a number greater than:
        std::cout << "quantile(complement(dist, 0.7)) = " << quantile(complement(dist, 0.7)) << std::endl;
    }
    // binomial distribution, probability of success 0.3, 20 trials in total
    {
        math::binomial_distribution<long double> dist(20, 0.3);
        // Returns PDF (density) at point x of distribution my_dist.
        // Probability of getting 0.0 in a t-distribution is:
        std::cout << "pdf(dist, 0.0) = " << pdf(dist, 0.0) << std::endl;
        // Returns CDF (integral from -infinity to point x) of distribution my_dist.
        // Probability of getting a number lesser than 0.0 in a t-distribution is:
        std::cout << "cdf(dist, 0.0) = " << cdf(dist, 0.0) << std::endl;
        // Probability of getting a number lesser than 1.0 in a t-distribution is:
        std::cout << "cdf(dist, 1.0) = " << cdf(dist, 1.0) << std::endl;
        // Probability of getting a number greater than 1.0 in a t-distribution is:
        std::cout << "1.0 - cdf(dist, 0.0) = " << 1.0 - cdf(dist, 0.0) << std::endl;
        std::cout << "1.0 - cdf(dist, 0.1) = " << 1.0 - cdf(dist, 0.1) << std::endl;
        // Removing the cancellation error with the complement distribution
        std::cout << "cdf(complement(dist, 0.0)) = " << cdf(complement(dist, 0.0)) << std::endl;
        std::cout << "cdf(complement(dist, 0.1)) = " << cdf(complement(dist, 0.1)) << std::endl;
        // Returns QUANTILE the value of the random variable x such that cdf(my_dist, x) == p.
        // This is the inverse of PDF
        // Note that lower critical values are just quantiles
        // Less than 50% of the time, we have a number smaller than:
        std::cout << "quantile(dist, 0.5)) = " << quantile(dist, 0.5) << std::endl;
        // Less than 30% of the time, we have a number smaller than:
        std::cout << "quantile(dist, 0.3)) = " << quantile(dist, 0.3) << std::endl;
        // Less than 70% of the time, we have a number smaller than:
        std::cout << "quantile(dist, 0.7)) = " << quantile(dist, 0.7) << std::endl;
        // Note that upper critical values are just quantiles
        // Less than 50% of the time, we have a number greater than:
        std::cout << "quantile(complement(dist, 0.5)) = " << quantile(complement(dist, 0.5)) << std::endl;
        // Less than 30% of the time, we have a number greater than:
        std::cout << "quantile(complement(dist, 0.3)) = " << quantile(complement(dist, 0.3)) << std::endl;
        // Less than 70% of the time, we have a number greater than:
        std::cout << "quantile(complement(dist, 0.7)) = " << quantile(complement(dist, 0.7)) << std::endl;
    }


    return 0;
}
