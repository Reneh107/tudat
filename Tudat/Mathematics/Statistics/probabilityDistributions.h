/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      160429    R. Hoogendoorn    File created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_PROBABILITY_DISTRIBUTIONS_H
#define TUDAT_PROBABILITY_DISTRIBUTIONS_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Householder>
#include <Eigen/QR>
#include <Eigen/Sparse>

#include <boost/shared_ptr.hpp>
#include <tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace statistics
{

//! Base class for probability distributions.
/*!
 * Base class for general probability distributions.
 * \tparam IndependentVariableType The type of the independent (random) variable.
 */
template< typename IndependentVariableType >
class ProbabilityDistribution
{
public:

    //! Destructor
    virtual ~ProbabilityDistribution(){}

    //! Get probability density.
    /*!
     * This function computes and returns the probability density at X = x.
     * The function uses the probability density function to compute this value.
     * This function is a pure virtual function, so it must be implemented by the derived class.
     * \param IndependentVariableType x sample of random variable.
     * \return probability density.
     */
    virtual double getProbabilityDensity(IndependentVariableType x) = 0 ;

    //! Get Probability Mass?
    /*!
     * The function uses the cumulative distribution function to compute this value.
     * F(x) = P(X < x)
     * \param IndependentVariableType x sample of random variable.
     * \return probability mass.
     */
//    virtual double getProbabilityMass(IndependentVariableType x) = 0 ;

    //! Get Quantile

    //! Get inverse Cumulative
    //! F^-1(x)

protected:

private:
};

typedef boost::shared_ptr< ProbabilityDistribution<double> > ProbabilityDistributionDoublePointer;

typedef boost::shared_ptr< ProbabilityDistribution<Eigen::VectorXd> > ProbabilityDistributionXdPointer;

//! One-dimensional Gaussian distribution class.
/*!
 * One-dimensional Gaussian distribution class.
 * Source: ??
 */
class GaussianDistributiond: public ProbabilityDistribution<double>{
public:

    //! Constructor.
    GaussianDistributiond(double Mean, double StandardDeviation);

    //! Get probability density.
    double getProbabilityDensity(double x);

protected:
private:
    //! Mean value of the distribution.
    double mean_ ;

    //! Standard deviation value of the distribution.
    double standardDeviation_ ;

    //! Variance value of the distribution.
    double variance_ ;
};

//! One-dimensional Uniform distribution class.
/*!
 * One-dimensional Uniform distribution class.
 * Source:
 */
class UniformDistributiond: public ProbabilityDistribution<double> {
public:

    //! Constructor.
    UniformDistributiond(double LowerBound, double UpperBound);

    //! Compute mean and standard deviation of distribution.
    void computeMeanAndStandardDeviation();

    //! Get mean value of distribution
    double getMean(){
        computeMeanAndStandardDeviation();
        return mean_;
    }

    //! Get standard deviation of distribution
    double getStandardDeviation(){
        computeMeanAndStandardDeviation();
        return standardDeviation_;
    }

    //! Get probability density.
    double getProbabilityDensity(double x);

protected:
private:
    //! Lower bound value of the distribution.
    double lowerBound_ ;

    //! Upper bound value of the distribution.
    double upperBound_ ;

    //! Mean value of the distribution.
    double mean_ ;

    //! Standard deviation value of the distribution.
    double standardDeviation_ ;

    //! Variance value of the distribution.
    double variance_ ;

    //! Probability density value of distribution (constant).
    double probabilityDensity_ ;
};


//! Multi-dimensional Gaussian Distribution class.
/*!
 * Multi-dimensional Gaussian Distribution class.
 * Source: ??
 */
template< int Dimension >
class GaussianDistributionXd: public ProbabilityDistribution<Eigen::VectorXd> {
public:

    //! Constructor
    GaussianDistributionXd(Eigen::VectorXd Mean , Eigen::MatrixXd CovarianceMatrix ){ // constructor
        mean_ = Mean ;
        covarianceMatrix_ = CovarianceMatrix ;
        dimension_ = double(mean_.rows()) ;
        determinant_ = CovarianceMatrix.determinant() ;

        // inverse function needs to know what the dimension of the matrix is at compile time
        typedef Eigen::Matrix<double,Dimension,Dimension> MatrixNd;
        MatrixNd Covariance = covarianceMatrix_ ;
        MatrixNd Inversecovariance = Covariance.inverse() ;

        inverseCovarianceMatrix_ = Inversecovariance ;
    }

    //! Get probability density
    double getProbabilityDensity(Eigen::VectorXd x){
        using tudat::mathematical_constants::PI;
        Eigen::VectorXd u = (x - mean_) ;
        Eigen::MatrixXd location = -0.5*( u.transpose()*inverseCovarianceMatrix_*u ) ;

        double probability = std::exp( location(0,0) ) /( std::pow(2.0*PI,dimension_/2.0) * std::sqrt(determinant_) ) ;
        return probability ;
    }

private:
    //! Dimension of the random variable X
    double              dimension_           ;

    //! Mean vector of random variable X
    Eigen::VectorXd     mean_                  ;

    //! Covariance matrix of random variable X
    Eigen::MatrixXd     covarianceMatrix_    ;

    //! Determinant of covariance matrix
    double              determinant_         ;

    //! Inverse of the covariance matrix
    Eigen::MatrixXd     inverseCovarianceMatrix_ ;
};


// //! Distribution of a random variable with uniform marginals and correlation using a Gaussian copula.
//template< int Dimension >
//class UniformCorrelatedDistributionND: public Distribution {
//public:
//    using Distribution::getProbability;

//    UniformCorrelatedDistributionND(Eigen::VectorXd Mu , Eigen::MatrixXd CovarianceMatrix ){ // constructor
//        mu = Mu ;
//        covariancematrix = CovarianceMatrix ;

//        // correlationmatrix
//        correlationmatrix = Thesis::Statistics::Basics::Covariance2CorrelationMatrix(covariancematrix) ;

//        dimension = Mu.rows() ;

//        // determinant correlation matrix
//        determinantCor = correlationmatrix.determinant() ;

//        // inverse function needs to know what the dimension of the matrix is at compile time
//        typedef Eigen::Matrix<double,Dimension,Dimension> MatrixNd;
//        MatrixNd Correlation = correlationmatrix ;
//        MatrixNd Inversecorrelation = Correlation.inverse() ;
//        inversecorrelation = Inversecorrelation ;

//        // Calculate bounds
//        LeftBound.resize( Dimension ) ;
//        RightBound.resize( Dimension ) ;
//        for(int i = 0 ; i < Dimension ; i++){
//            LeftBound(i) = mu(i) - std::sqrt(3.0 * covariancematrix(i,i) ) ;
//            RightBound(i) = mu(i) + std::sqrt(3.0 * covariancematrix(i,i) ) ;
//        }
//        width = RightBound - LeftBound ;
//        volume = width.prod() ; // product of marginal PDFs

////        std::cout << "LeftBound = " << LeftBound << std::endl;
////        std::cout << "RightBound = " << RightBound << std::endl;
//    }

//    double getProbability(Eigen::VectorXd x){
//        double probability = 0.0 ;

//        int InBound = 0 ;
//        for(int i = 0 ; i < dimension ; i++){ // check in bounds
//            if( x(i) > LeftBound(i) && x(i) < RightBound(i) ){
//                InBound++ ;
//            }
//        }
//        if(InBound == dimension){
//            // Convert uniform to U = [0,1] :  F(x) = x - a / (b-a) -> U = F_X(X)
//            Eigen::VectorXd Fx = (x - LeftBound).cwiseQuotient( width ) ;

//            // Convert U[0,1] to N[0,1] using inverse CDF of standard normal distribution
//            Eigen::VectorXd y(Fx.rows()) ;
//            for(int i = 0 ; i < Fx.rows() ; i++){
//                y(i) = gsl_cdf_ugaussian_Pinv(Fx(i)) ; // Inverse standard normal CDF
//            }

//            // Calculate probability density
//            Eigen::MatrixXd location = -0.5*( y.transpose()*(inversecorrelation - Eigen::MatrixXd::Identity(Dimension,Dimension) )*y ) ;

//            probability = ( ( 1.0/( volume * std::sqrt(determinantCor) ) )*std::exp( location(0,0) ) ) ;
//        }

//        return probability ;
//    }

//private:
//    int                 dimension           ; // dimension of distribution
//    Eigen::VectorXd     mu                  ; // vector of mean values
//    Eigen::MatrixXd     covariancematrix    ; // covariance matrix of x

//    Eigen::MatrixXd     correlationmatrix   ; // correlationmatrix of x
//    double              determinantCor      ; // determinant of correlation matrix
//    Eigen::MatrixXd     inversecorrelation  ; // inverse of the correlation matrix of x

//    Eigen::VectorXd     LeftBound           ; // leftbound of uniform marginal distributions
//    Eigen::VectorXd     RightBound          ; // rightbound of uniform marginal distributions
//    Eigen::VectorXd     width               ; // width of uniform marginal distributions
//    double              volume              ; // volume of joint uniform distribution , product of marginal PDFs
//};


} // namespace statistics
} // namespace tudat

#endif // TUDAT_PROBABILITY_DISTRIBUTIONS_H
