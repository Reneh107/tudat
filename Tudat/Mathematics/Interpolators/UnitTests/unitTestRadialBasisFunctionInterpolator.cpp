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
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical( ) function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean( )
 *                                  and computeSampleVariance( ) functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120202    K. Kumar          Separated from unitTestBasicMathematics.cpp into new
 *                                  Interpolators sub-directory.
 *      120529    E.A.G. Heeren     Boostified unit test.
 *      120615    T. Secretin       Minor layout changes.
 *      120716    D. Dirkx          Updated with interpolator architecture.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Mathematics/Interpolators/radialBasisFunctionInterpolator.h"

#include <Eigen/Core>

namespace tudat
{
namespace unit_tests
{

//! Test function
//! z = x(0)^2 + x(0)*x(1)^2
double testFunction( Eigen::VectorXd x )
{
    return ( std::pow( x(0) , 2.0) + x(0)*std::pow( x(1) , 2.0 ) );
}

//! Test function2
double testFunction2( Eigen::VectorXd x )
{
    return ( x(0)*std::pow( std::sin( x(0) ) , 2.0 ) * std::exp( - std::pow( x(1) , 2.0 ) ) );
}

BOOST_AUTO_TEST_SUITE( test_radial_basis_function_interpolator )

// Test implementation of the gaussian radial basis function
BOOST_AUTO_TEST_CASE( test_gaussian_radial_basis_function )
{
    tudat::interpolators::GaussianRadialBasisFunction gaussianRadialBasis( 0.5 );

    Eigen::VectorXd variable(3);
    Eigen::VectorXd location(3);

    {
        variable << 1.0 , 2.0 , -1.0 ;
        location << 0.0 , 2.0 , 0.5 ;

        double computedValue = gaussianRadialBasis.evaluate( variable , location );

        // Generated using MATLAB ( Spherical Radial Basis Functions, Theory and Applications , Simon Hubbert , 2015 )
        double expectedValue = 1.969116752041940e-01 ;

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-15 );
    }

    {
        variable << 5.0 , 3.0 , -0.5 ;
        location << -1.0 , 2.0 , -0.5 ;

        double computedValue = gaussianRadialBasis.evaluate( variable , location );

        // Generated using MATLAB ( Spherical Radial Basis Functions, Theory and Applications , Simon Hubbert , 2015 )
        double expectedValue = 9.237449661970661e-09 ;

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-15 );
    }
}

// Test implementation of the multiquadric radial basis function
BOOST_AUTO_TEST_CASE( test_multiquadric_radial_basis_function )
{
    tudat::interpolators::MultiQuadricRadialBasisFunction radialBasisFuntion( 0.8 );

    Eigen::VectorXd variable(4);
    Eigen::VectorXd location(4);

    {
        variable << 1.0 , 2.0 , -1.0 , 0.5 ;
        location << 0.0 , 2.0 , 0.5 , 2.0 ;

        double computedValue = radialBasisFuntion.evaluate( variable , location );

        // Generated using MATLAB ( Spherical Radial Basis Functions, Theory and Applications , Simon Hubbert , 2015 )
        double expectedValue = 2.477902338672774e+00 ;

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-15 );
    }

    {
        variable << 5.0 , 3.0 , -0.5 , 1.0 ;
        location << -1.0 , 2.0 , -0.5 , 3.0 ;

        double computedValue = radialBasisFuntion.evaluate( variable , location );

        // Generated using MATLAB ( Spherical Radial Basis Functions, Theory and Applications , Simon Hubbert , 2015 )
        double expectedValue = 6.452906321960671e+00 ;

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-15 );
    }
}

// Test implementation of the multiquadric radial basis function
BOOST_AUTO_TEST_CASE( test_inverse_multiquadric_radial_basis_function )
{
    tudat::interpolators::InverseMultiQuadricRadialBasisFunction radialBasisFuntion( 0.8 );

    Eigen::VectorXd variable(4);
    Eigen::VectorXd location(4);

    {
        variable << 1.0 , 2.0 , -1.0 , 0.5 ;
        location << 0.0 , 2.0 , 0.5 , 2.0 ;

        double computedValue = radialBasisFuntion.evaluate( variable , location );

        // Generated using MATLAB ( Spherical Radial Basis Functions, Theory and Applications , Simon Hubbert , 2015 )
        double expectedValue = 4.035671561356308e-01 ;

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-15 );
    }

    {
        variable << 5.0 , 3.0 , -0.5 , 1.0 ;
        location << -1.0 , 2.0 , -0.5 , 3.0 ;

        double computedValue = radialBasisFuntion.evaluate( variable , location );

        // Generated using MATLAB ( Spherical Radial Basis Functions, Theory and Applications , Simon Hubbert , 2015 )
        double expectedValue = 1.549689318434359e-01 ;

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-15 );
    }
}

// Test implementation of radial basis function interpolator
//BOOST_AUTO_TEST_CASE( test_radial_basis_function_interpolator_2d )
//{
//    std::cout << "Test 2D interpolation" << std::endl;

//    using namespace tudat::interpolators;

//    // Generate random independent values
//    std::vector< Eigen::VectorXd > independentValues(0);
//    Eigen::VectorXd randomSample(2);

//    int numberOfSamples = 4E2 ;
//    for( int i = 0 ; i < numberOfSamples ; i++ )
//    {
//        randomSample = Eigen::VectorXd::Random(2) * 20 ;
//        independentValues.push_back( randomSample );
//    }

//    // Generate dependentValues
//    std::vector< double > dependentValues(0);
//    for( int i = 0 ; i < numberOfSamples ; i++ )
//    {
//        dependentValues.push_back( testFunction( independentValues[i] ) );
//    }

//    // Construct interpolator
////    double scaleParameter = 25.0;
//    double scaleParameter = 0.1;

//    RadialBasisFunctionInterpolator interpolator(
//                independentValues , dependentValues , RadialBasisFunctionType::InverseMultiquadric , scaleParameter );

//    BOOST_CHECK_CLOSE_FRACTION( scaleParameter , interpolator.getShapeParameter() , 1E-15 );
//    BOOST_CHECK_CLOSE_FRACTION( scaleParameter * 40.0 , interpolator.getScaledShapeParameter() , 2E-3 );

//    // Test if correct value at data sample
//    {
//        // Interpolate value
//        Eigen::VectorXd x(2);
//        x = independentValues[0];

//        double computedValue = interpolator.interpolate( x );
//        double expectedValue = testFunction( x );

//        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-12 );
//    }

//    // Test if correct value at data sample
//    {
//        // Interpolate value
//        Eigen::VectorXd x(2);
//        x = independentValues[30];

//        double computedValue = interpolator.interpolate( x );
//        double expectedValue = testFunction( x );

//        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-12 );
//    }

//    // Test if correct value at data sample
//    {
//        // Interpolate value
//        Eigen::VectorXd x(2);
//        x = independentValues[39];

//        double computedValue = interpolator.interpolate( x );
//        double expectedValue = testFunction( x );

//        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-12 );
//    }

//    // Test if functions correctly interpolate
//    {
//        // Interpolate value
//        Eigen::VectorXd x(2);
//        x << 2.0 , 0.0 ;

//        double computedValue = interpolator.interpolate( x );
//        double expectedValue = testFunction( x );
//        std::cout << "computed " << computedValue << std::endl;
//        std::cout << "expected " << expectedValue << std::endl;

//        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 3E-4 );
//    }

//    // Test if functions correctly interpolate
//    {
//        // Interpolate value
//        Eigen::VectorXd x(2);
//        x = independentValues[39] * 0.9 ;

//        double computedValue = interpolator.interpolate( x );
//        double expectedValue = testFunction( x );
//        std::cout << "computed " << computedValue << std::endl;
//        std::cout << "expected " << expectedValue << std::endl;

//        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 2E-3 );
//        BOOST_CHECK_CLOSE_FRACTION( computedValue , expectedValue , 1E-8 );
//    }

//}

// Test implementation of radial basis function interpolator
BOOST_AUTO_TEST_CASE( test_radial_basis_function_interpolator_2d_testfunction2 )
{
    using namespace tudat::interpolators;

    // Generate random independent values
    std::vector< Eigen::VectorXd > independentValues(0);
    Eigen::VectorXd randomSample(2);

    int numberOfSamples = 4E2 ;
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        randomSample = Eigen::VectorXd::Random(2) * 20 ;
        independentValues.push_back( randomSample );
    }

    // Generate dependentValues
    std::vector< double > dependentValues(0);
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        dependentValues.push_back( testFunction2( independentValues[i] ) );
    }

    // Construct interpolator
    double scaleParameter = 0.1;
    RadialBasisFunctionInterpolator interpolator(
                independentValues , dependentValues , RadialBasisFunctionType::InverseMultiquadric , scaleParameter );

    // Test if correct value at data sample
    {
        // Interpolate value
        Eigen::VectorXd x(2);
        x = independentValues[0];

        double computedValue = interpolator.interpolate( x );
        double expectedValue = testFunction2( x );

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-10 );
    }

    // Test if correct value at data sample
    {
        // Interpolate value
        Eigen::VectorXd x(2);
        x = independentValues[30];

        double computedValue = interpolator.interpolate( x );
        double expectedValue = testFunction2( x );

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-10 );
    }

    // Test if correct value at data sample
    {
        // Interpolate value
        Eigen::VectorXd x(2);
        x = independentValues[39];

        double computedValue = interpolator.interpolate( x );
        double expectedValue = testFunction2( x );

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 1E-10 );
    }

    // Test if functions correctly interpolate
    {
        // Interpolate value
        Eigen::VectorXd x(2);
        x << 2.0 , 0.0 ;

        double computedValue = interpolator.interpolate( x );
        double expectedValue = testFunction2( x );
        std::cout << "computed " << computedValue << std::endl;
        std::cout << "expected " << expectedValue << std::endl;

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 3E-4 );
    }

    // Test if functions correctly interpolate
    {
        // Interpolate value
        Eigen::VectorXd x(2);
        x = independentValues[39] * 0.9 ;

        double computedValue = interpolator.interpolate( x );
        double expectedValue = testFunction2( x );
//        std::cout << "computed " << computedValue << std::endl;
//        std::cout << "expected " << expectedValue << std::endl;

        BOOST_CHECK_SMALL( std::fabs( computedValue - expectedValue ) , 5E-3 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
