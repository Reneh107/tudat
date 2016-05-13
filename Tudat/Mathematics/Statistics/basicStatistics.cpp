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
#include <algorithm>

#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

namespace tudat
{
namespace statistics
{

//! Compute average of the components of a vector.
double computeAverageOfVectorComponents( const Eigen::VectorXd& vectorOfData )
{ return vectorOfData.sum( ) / vectorOfData.rows( ); }

//! Compute standard deviation of the components of a vector.
double computeStandardDeviationOfVectorComponents( const Eigen::VectorXd& vectorOfData )
{
    // Compute average of components.
    double averageOfComponents = computeAverageOfVectorComponents( vectorOfData );

    // Declare variance of components.
    double varianceOfComponents = 0.0;

    // Compute variance of components.
    for ( int i = 0 ; i < vectorOfData.rows( ) ; i++ )
    {
        varianceOfComponents += std::pow( ( vectorOfData( i ) - averageOfComponents ), 2.0 );
    }

    varianceOfComponents /= static_cast< double >( vectorOfData.rows( ) - 1 );

    // Return square root of variance ( = standard deviation ).
    return std::sqrt( varianceOfComponents );
}

//! Compute sample mean.
double computeSampleMean( const std::vector< double >& sampleData )
{
    // Return sample mean.
    return std::accumulate( sampleData.begin( ), sampleData.end( ), 0.0 )
            / static_cast< double >( sampleData.size( ) );
}

//! Compute sample variance.
double computeSampleVariance( const std::vector< double >& sampleData )
{
    // Declare local variables.
    // Declare and compute sample mean.
    double sampleMean_ = computeSampleMean( sampleData );

    // Declare and initialize sum of residuals squared.
    double sumOfResidualsSquared_ = 0.0;

    // Compute sum of residuals of sample data squared.
    for ( unsigned int i = 0; i < sampleData.size( ); i++ )
    {
        sumOfResidualsSquared_ += std::pow( sampleData.at( i ) - sampleMean_, 2.0 );
    }

    // Return sample variance.
    return 1.0 / ( static_cast< double >( sampleData.size( ) ) - 1.0 ) * sumOfResidualsSquared_;
}

//! Compute Sample median
double computeSampleMedian( std::vector< double > sampleData )
{
    // Sort data
    std::sort( sampleData.begin() , sampleData.end());

    // Check if odd number of samples or even
    double numberOfSamples = static_cast< double >( sampleData.size() );
    int odd = tudat::basic_mathematics::computeModulo( numberOfSamples , 2.0 ) ;

    // Calculate sample median
    double sampleMedian;
    if( odd == 0 ) // even
    {
        // 0 .. 99 (100) ->
        int index = static_cast< int >( numberOfSamples / 2.0 - 0.5 ) ;
        sampleMedian = (sampleData[index] + sampleData[index+1])/2.0 ;
    }
    else // odd
    {
        // 0 .. 100 (101) -> 50
        int index = static_cast< int >( numberOfSamples / 2.0 - 0.5 ) ;
        sampleMedian = sampleData[index] ;
    }

    return sampleMedian;
}

} // namespace statistics
} // namespace tudat
