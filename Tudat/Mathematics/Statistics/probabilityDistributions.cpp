/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      090807    J. Melman         File created.
 *      100930    D. Dirkx          Modified to comply with Tudat standards
 *      100930    J. Melman         Implemented namespace, minor comment
 *                                  changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120202    K. Kumar          Moved functions from linearAlgebra.cpp.
 *
 *    References
 *
 *    Notes
 *
 */

#include <cmath>
#include <numeric>

#include <tudat/Mathematics/Statistics/probabilityDistributions.h>

namespace tudat
{
namespace statistics
{

using tudat::mathematical_constants::PI;

//! Constructor.
GaussianDistributiond::GaussianDistributiond(double Mean, double StandardDeviation){ // constructor
    mean_ = Mean ;
    standardDeviation_ = StandardDeviation ;
    variance_ = std::pow( standardDeviation_ , 2.0 );
}

//! Get probability density of 1D Gaussian distribution
double GaussianDistributiond::getProbabilityDensity(double x){
    return (std::exp( ( - std::pow( x - mean_ , 2.0 ) ) / ( 2.0 * variance_ ) )
            /( std::sqrt( 2.0 * PI ) * standardDeviation_ ) ) ;
}

//! Constructor.
UniformDistributiond::UniformDistributiond(double LowerBound, double UpperBound){ // constructor
    lowerBound_ = LowerBound ;
    upperBound_ = UpperBound ;
    probabilityDensity_ = 1.0/( upperBound_ - lowerBound_ ) ;
}


//! Get probability density of 1D Uniform distribution
double UniformDistributiond::getProbabilityDensity(double x){
    if ( (x >= lowerBound_ && x <= upperBound_ ) ){
        return probabilityDensity_ ;
    }
    else{
        return 0.0 ;
    }
}

//! Compute mean and standard deviation 1D Uniform distribution.
void UniformDistributiond::computeMeanAndStandardDeviation(){
    mean_ = (upperBound_ + lowerBound_)/2.0 ;
    variance_ = std::pow(upperBound_ - lowerBound_,2.0) / 12.0;
    standardDeviation_ = std::sqrt( variance_ );
}

} // namespace statistics
} // namespace tudat
