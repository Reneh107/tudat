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
 *      160520    R. Hoogendoorn    File created.
 *
 *    References
 *
 *
 *    Notes
 *
 *
 */

#ifndef TUDAT_RADIAL_BASIS_FUNCTION_INTERPOLATOR_H
#define TUDAT_RADIAL_BASIS_FUNCTION_INTERPOLATOR_H

#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/Statistics/basicStatistics.h"

namespace tudat
{
namespace interpolators
{

//! Radial basis function type
enum RadialBasisFunctionType{
    Gaussian = 0,
    Multiquadric = 1,
    InverseMultiquadric = 2,
    AssymGaussian = 3
};

//! Class that implements a radial basis function
/*!
 * ..
 *
 * Note that the types (i.e. double, float) of all independent variables must be the same.
 *
 */
class RadialBasisFunction
{
public:

    //! Default destructor
    virtual ~RadialBasisFunction(){}

    //! Evaluate function
    /*!
     * This function returns the value of the radial basis function.
     *
     * This function is pure virtual, so it must be implemented in the derived class.
     * \param
     * \param
     * \return Function value
     */
    virtual double evaluate( const Eigen::VectorXd& independentVariable, const Eigen::VectorXd& location ) = 0 ;

protected:

private:

};

typedef boost::shared_ptr< RadialBasisFunction > RadialBasisFunctionPointer;

//! Class that implements a radial basis function
/*!
 * ..
 *
 * Note that the types (i.e. double, float) of all independent variables must be the same.
 *
 */
class GaussianRadialBasisFunction: public RadialBasisFunction
{
public:

    GaussianRadialBasisFunction( double shapeParameter )
    {
        shapeParameter_ = shapeParameter;
    }

    double evaluate( const Eigen::VectorXd& independentVariable, const Eigen::VectorXd& location )
    {
        Eigen::VectorXd difference = independentVariable - location ;
        double radius = difference.norm();
        return std::exp( - shapeParameter_ * std::pow( radius , 2.0 ) );
    }

protected:

private:

    //! Shape parameter of radial basis function
    double shapeParameter_;

};

//! Class that implements a radial basis function
/*!
 * ..
 *
 * Note that the types (i.e. double, float) of all independent variables must be the same.
 *
 */
class AssymetricGaussianRadialBasisFunction: public RadialBasisFunction
{
public:

    AssymetricGaussianRadialBasisFunction( double shapeParameter , Eigen::VectorXd standardDeviationOfData )
    {
        shapeParameter_ = shapeParameter;
        standardDeviationOfData_ = standardDeviationOfData;
    }

    double evaluate( const Eigen::VectorXd& independentVariable, const Eigen::VectorXd& location )
    {
        Eigen::VectorXd difference = independentVariable - location ;
        double radius = (difference.cwiseQuotient( standardDeviationOfData_ ) ).norm();
        return std::exp( - shapeParameter_ * std::pow( radius , 2.0 ) );
    }

protected:

private:

    //! Shape parameter of radial basis function
    double shapeParameter_;

    //! Variance of to be interpolated data
    Eigen::VectorXd standardDeviationOfData_;

};

//! Class that implements a radial basis function
/*!
 * ..
 *
 * Note that the types (i.e. double, float) of all independent variables must be the same.
 *
 */
class MultiQuadricRadialBasisFunction: public RadialBasisFunction
{
public:

    MultiQuadricRadialBasisFunction( double shapeParameter )
    {
        shapeParameter_ = shapeParameter;
    }

    double evaluate( const Eigen::VectorXd& independentVariable, const Eigen::VectorXd& location )
    {
        using namespace std;
        Eigen::VectorXd difference = independentVariable - location ;
        double radius = difference.norm();
        return pow( pow( radius , 2.0 ) + pow( shapeParameter_ , 2.0 ) , 0.5 );
    }

protected:

private:

    //! Shape parameter of radial basis function
    double shapeParameter_;

};

//! Class that implements a radial basis function
/*!
 * ..
 *
 * Note that the types (i.e. double, float) of all independent variables must be the same.
 *
 */
class InverseMultiQuadricRadialBasisFunction: public RadialBasisFunction
{
public:

    InverseMultiQuadricRadialBasisFunction( double shapeParameter )
    {
        shapeParameter_ = shapeParameter;
    }

    double evaluate( const Eigen::VectorXd& independentVariable, const Eigen::VectorXd& location )
    {
        using namespace std;
        Eigen::VectorXd difference = independentVariable - location ;
        double radius = difference.norm();
        return pow( pow( radius , 2.0 ) + pow( shapeParameter_ , 2.0 ) , -0.5 );
    }

protected:

private:

    //! Shape parameter of radial basis function
    double shapeParameter_;

};

//! Class for performing radial basis function interpolation of scattered data.
/*!
 * Class for performing interpolation of scattered multi-dimensional data. Radial basis functions are
 * used to fit a function through the data and generate the interpolants.
 *
 * Note that the types (i.e. double, float) of all independent variables must be the same.
 *
 */
class RadialBasisFunctionInterpolator: public Interpolator< double , double >
{
public:

    //! Default Constructor
    /*!
     *
     */
    RadialBasisFunctionInterpolator( std::vector< Eigen::VectorXd > independentValues,
                                     std::vector< double > dependentValues,
                                     RadialBasisFunctionType radialBasisFunctionType,
                                     double shapeParameter = 1.0,
                                     double weightDecayFactor = 0.0):
                                    independentValues_( independentValues ),
                                    dependentData_( dependentValues ),
                                    radialBasisFunctionType_( radialBasisFunctionType ),
                                    shapeParameter_( shapeParameter ),
                                    weightDecayFactor_( weightDecayFactor )
    {
        numberOfDimension_ = static_cast< int >( independentValues_[0].rows() );
        numberOfDatapoints_ = static_cast< int >( independentValues_.size() );

        calculateVariance();

        // Scale shape parameter using standard deviation of independent variables
        shapeParameterScaled_ = shapeParameter_ * ( standardDeviationOfData_.maxCoeff() );

        constructRadialBasisFunction();

        generateCoefficients();
    }

    //! Constructor
    /*!
     *
     */
    RadialBasisFunctionInterpolator( RadialBasisFunctionType radialBasisFunctionType,
                                     double shapeParameter = 1.0,
                                     double weightDecayFactor = 0.0):
                                    radialBasisFunctionType_( radialBasisFunctionType ),
                                    shapeParameter_( shapeParameter ),
                                    weightDecayFactor_( weightDecayFactor )
    { }

    //! Construct the radial basis function.
    void constructRadialBasisFunction( )
    {
        // Construct radial basis function
        if( radialBasisFunctionType_ == RadialBasisFunctionType::Gaussian )
        {
            radialBasisFunction_ = boost::make_shared< GaussianRadialBasisFunction >( shapeParameterScaled_ );
        }
        else if( radialBasisFunctionType_ == RadialBasisFunctionType::Multiquadric )
        {
            radialBasisFunction_ = boost::make_shared< MultiQuadricRadialBasisFunction >( shapeParameterScaled_ );

        }
        else if( radialBasisFunctionType_ == RadialBasisFunctionType::InverseMultiquadric )
        {
            radialBasisFunction_ = boost::make_shared< InverseMultiQuadricRadialBasisFunction >(
                        shapeParameterScaled_ );
        }
        else if( radialBasisFunctionType_ == RadialBasisFunctionType::AssymGaussian )
        {
            radialBasisFunction_ = boost::make_shared< AssymetricGaussianRadialBasisFunction >(
                        shapeParameter_ , standardDeviationOfData_ );
        }
    }

    //! Default destructor
    /*!
     *  Default destructor
     */
    ~RadialBasisFunctionInterpolator( ){ }

    //! Function to perform interpolation.
    /*!
     *  This function performs t
     *  \param targetIndependentVariableValue Vector of values of independent variables at which
     *  the value of the dependent variable is to be determined.
     *  \return Interpolated value of dependent variable in all dimensions.
     */
    double interpolate( const Eigen::VectorXd& targetIndependentVariableValue )
    {
        // interpolate
        double interpolatedValue = 0.0;
        for( int i = 0 ; i < numberOfDatapoints_ ; i++ )
        {
            interpolatedValue += coefficients_( i ) * radialBasisFunction_->evaluate(
                        targetIndependentVariableValue , independentValues_[ i ] ) ;
        }
        return interpolatedValue;
    }

    //! Interpolate.
    /*!
     * Interpolation function is not implemented for an input variable with type std::vector.
     * \param targetIndependentVariableValue Target independent variable value at which point
     * the interpolation is performed.
     * \return Interpolated dependent variable value.
     */
    double interpolate( const std::vector< double > & targetIndependentVariableValue )
    {
        std::cerr << "Radial basis function interpolator works with an Eigen::VectorXd" << std::endl;
        return TUDAT_NAN;
    }

    //! Function to return the number of independent variables of the interpolation.
    /*!
     *  Function to return the number of independent variables of the interpolation, i.e. size
     *  that the vector used as input for Interpolator::interpolate should be.
     *  \return Number of independent variables of the interpolation.
     */
    int getNumberOfDimensions( )
    {
        return numberOfDimension_;
    }

    //! Get the scaled shape parameter.
    double getScaledShapeParameter( )
    {
        return shapeParameterScaled_;
    }

    //! Get the shape parameter of the radial basis functions.
    double getShapeParameter( )
    {
        return shapeParameter_;
    }

    //! Get the variance of the independent variables.
    Eigen::VectorXd getVarianceOfData( )
    {
        return varianceOfData_;
    }

    //! Get the standard deviation of the independent variables.
    Eigen::VectorXd getStandardDeviationOfData( )
    {
        return standardDeviationOfData_;
    }

    Eigen::VectorXd getCoefficients( )
    {
        return coefficients_;
    }

    void setCoefficients( Eigen::VectorXd coefficients )
    {
        coefficients_ = coefficients;
    }

    void setIndependentValues( std::vector< Eigen::VectorXd > independentValues )
    {
        independentValues_ = independentValues;
        numberOfDatapoints_ = independentValues_.size();
    }

    void computeInterpolationError( )
    {
        std::vector< double > absoluteInterpolationError( numberOfDatapoints_ );
        std::vector< double > relativeInterpolationError( numberOfDatapoints_ );
        std::vector< double > interpolationError( numberOfDatapoints_ );
        for( int i = 0 ; i < numberOfDatapoints_ ; i++ )
        {
            interpolationError[i] = ( interpolate( independentValues_[ i ] ) - dependentData_[ i ] );
            absoluteInterpolationError[i] = std::abs( interpolationError[i] );
            relativeInterpolationError[i] = absoluteInterpolationError[i] / dependentData_[ i ];
        }

        meanAbsoluteInterpolationError_ = tudat::statistics::computeSampleMean( absoluteInterpolationError );
        meanRelativeInterpolationError_ = tudat::statistics::computeSampleMean( relativeInterpolationError );
        meanInterpolationError_ = tudat::statistics::computeSampleMean( interpolationError );

        varianceAbsoluteInterpolationError_ = tudat::statistics::computeSampleVariance( absoluteInterpolationError );
        varianceRelativeInterpolationError_ = tudat::statistics::computeSampleVariance( relativeInterpolationError );
        varianceInterpolationError_ = tudat::statistics::computeSampleVariance( interpolationError );
    }

    std::vector< double > getInterpolationErrorProperties( )
    {
        std::vector< double > interpolationErrorProperties(0);
        interpolationErrorProperties.push_back( meanAbsoluteInterpolationError_ );
        interpolationErrorProperties.push_back( meanRelativeInterpolationError_ );
        interpolationErrorProperties.push_back( meanInterpolationError_ );
        interpolationErrorProperties.push_back( varianceAbsoluteInterpolationError_ );
        interpolationErrorProperties.push_back( varianceRelativeInterpolationError_ );
        interpolationErrorProperties.push_back( varianceInterpolationError_ );
        return interpolationErrorProperties;
    }

    RadialBasisFunctionPointer getRadialBasisFunction( )
    {
        return radialBasisFunction_;
    }

    //! Generate the coefficients of the radial basis functions.
    void generateCoefficients( )
    {
        if( ( weightDecayFactor_ - 0.0 ) > 1.0E-16 )
        {
            generateCoefficientsWeightDecay( );
        }
        else
        {
            generateCoefficientsExactInterpolation( );
        }
    }

private:

    //! Generate coefficients for exact interpolation (without weight decay regularization)
    void generateCoefficientsExactInterpolation( )
    {
        // Declare matrix A of : A c = y
        Eigen::MatrixXd matrixA( numberOfDatapoints_ , numberOfDatapoints_ );

        // Fill matrix
        for( int row = 0 ; row < matrixA.rows() ; row++ )
        {
            for( int col = 0 ; col < matrixA.cols() ; col++ )
            {
                matrixA( row , col ) = radialBasisFunction_->evaluate(
                            independentValues_[ col ] , independentValues_[ row ] ) ;
            }
        }

        // Generate vector y;
        Eigen::VectorXd vectorY( numberOfDatapoints_ );
        for( int i = 0 ; i < numberOfDatapoints_ ; i++ )
        {
            vectorY(i) = dependentData_[i];
        }

        // Solve coefficients
        coefficients_ = matrixA.colPivHouseholderQr( ).solve( vectorY ) ;
    }

    //! Generate coefficients for RBF interpolation with weight decay regularization.
    void generateCoefficientsWeightDecay( )
    {
        // min ( norm(Ac-y)^2 + weighDecayFactor * norm(w)^2 )
        // Declare matrix A of : A c = y
        Eigen::MatrixXd matrixA( numberOfDatapoints_ , numberOfDatapoints_ );

        // Fill matrix
        for( int row = 0 ; row < matrixA.rows() ; row++ )
        {
            for( int col = 0 ; col < matrixA.cols() ; col++ )
            {
                matrixA( row , col ) = radialBasisFunction_->evaluate(
                            independentValues_[ col ] , independentValues_[ row ] ) ;
            }
        }

        // Generate vector y;
        Eigen::VectorXd vectorY( numberOfDatapoints_ );
        for( int i = 0 ; i < numberOfDatapoints_ ; i++ )
        {
            vectorY(i) = dependentData_[i];
        }

        Eigen::VectorXd ATvectorY = matrixA.transpose() * vectorY;

        Eigen::MatrixXd identity = Eigen::MatrixXd::Identity( matrixA.rows() , matrixA.cols() );
        Eigen::MatrixXd augmentenMatrixA = matrixA.transpose() * matrixA + weightDecayFactor_ * identity ;

        // Solve coefficients
        coefficients_ = augmentenMatrixA.colPivHouseholderQr( ).solve( ATvectorY ) ;
    }

    void calculateVariance( )
    {
        // Calculate mean
        Eigen::VectorXd mean = Eigen::VectorXd::Zero( numberOfDimension_ );
        for( int i = 0 ; i < numberOfDatapoints_ ; i++ )
        {
            mean = mean + independentValues_[i];
        }
        mean = mean / ( static_cast<double>( numberOfDatapoints_ ) );

        // Calculate standard deviation
        varianceOfData_ = Eigen::VectorXd::Zero( numberOfDimension_ );
        for( int i = 0 ; i < numberOfDatapoints_ ; i++ )
        {
            varianceOfData_ = varianceOfData_ + ( (independentValues_[i] - mean).cwiseAbs2() );
        }
        varianceOfData_ = varianceOfData_ / static_cast< double >( independentValues_.size( ) - 1 );

        standardDeviationOfData_ = varianceOfData_.cwiseSqrt();
    }

    //! Vector of Eigen vectors containing independent variables.
    std::vector< Eigen::VectorXd > independentValues_;

    //! Variance vector of independent values.
    Eigen::VectorXd varianceOfData_;

    //! Standard deviation of independent values.
    Eigen::VectorXd standardDeviationOfData_;

    //! Vector of dependent data.
    std::vector< double > dependentData_;

    //! Number of dimensions
    int numberOfDimension_;

    //! Number of datapoints
    int numberOfDatapoints_;

    //! Minimum value of all the independent values.
    double minimumIndependentValue_;

    //! Maximum value of all the independent values.
    double maximumIndependentValue_;

    //! Coefficients of the radial basis functions used in the interpolation.
    Eigen::VectorXd coefficients_;

    //! Radial basis function type
    RadialBasisFunctionType radialBasisFunctionType_;

    //! Shape parameter of the radial basis functions.
    double shapeParameter_;

    //! Scaled shape parameter of the radial basis functions.
    double shapeParameterScaled_;

    //! Factor for applying weight-decay regularization (value of 0 means no weight decay)
    double weightDecayFactor_;

    //! Pointer to the radial basis function
    RadialBasisFunctionPointer radialBasisFunction_;

    double meanAbsoluteInterpolationError_;
    double meanRelativeInterpolationError_;
    double meanInterpolationError_;

    double varianceAbsoluteInterpolationError_;
    double varianceRelativeInterpolationError_;
    double varianceInterpolationError_;

};

//! A pointer type to a radial basis function interpolator
typedef boost::shared_ptr< RadialBasisFunctionInterpolator > RadialBasisFunctionInterpolatorPointer;

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_RADIAL_BASIS_FUNCTION_INTERPOLATOR_H
