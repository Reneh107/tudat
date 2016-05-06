/*   Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120601    A. Ronse          File created
 *    References
 *      USNO (2011), The Astronomical Almanac Online, The U.S. Naval Observatory.
 *
 */

#include <tudat/Astrodynamics/ReferenceFrames/referenceAngles.h>
#include <tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h>
#include <tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h>

#include <cmath>

#include <tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h>

namespace tudat
{

namespace reference_frames
{

//! Compute the angle between the inertial and rotating geocentric reference frame.
double computeVernalOffset( double julianDate )
{
    using tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000;
    using tudat::basic_mathematics::computeModulo ;

    double julianDatePreviousMidnight = std::floor(julianDate-0.5)+0.5 ; // 0h UT1
    double julianDateSinceJ2000 = julianDatePreviousMidnight - JULIAN_DAY_ON_J2000 ;
    double centuries0hUT1SinceJ2000 = julianDateSinceJ2000/36525.0; // T0
    double centuriesSinceJ2000 = ( julianDate - JULIAN_DAY_ON_J2000 )/36525.0; // T

    // Compute GMST
    double greenwhichMeanSiderealTime = 24110.54841
            + 8640184.812866 * centuries0hUT1SinceJ2000
            + 1.00273790935 * (julianDate - julianDatePreviousMidnight) * 86400.0 +
            + 0.093104 * centuriesSinceJ2000 * centuriesSinceJ2000
            - centuriesSinceJ2000 * centuriesSinceJ2000 * centuriesSinceJ2000 * 6.2E-6 ;

    // Return vernal offset
    return 2.0 * tudat::mathematical_constants::PI * computeModulo<double>(greenwhichMeanSiderealTime , 86400.0) / 86400.0;
}

//! Compute aerodynamic angles
Eigen::Vector3d computeAerodynamicAngles(double longitude,
                         double latitude,
                         double time,
                         double flightPathAngle,
                         double headingAngle,
                         Eigen::Matrix3d bodyFrameToInertialFrameTransformationMatrix){

    Eigen::Quaterniond quaternion( bodyFrameToInertialFrameTransformationMatrix ) ;

    return computeAerodynamicAngles(longitude,
                                    latitude,
                                    time,
                                    flightPathAngle,
                                    headingAngle,
                                    quaternion);
}

//! Compute aerodynamic angles
Eigen::Vector3d computeAerodynamicAngles(double longitude,
                         double latitude,
                         double time,
                         double flightPathAngle,
                         double headingAngle,
                         Eigen::Quaterniond bodyFrameToInertialFrameTransformationQuaternion){
    using namespace tudat::reference_frames;

    Eigen::Matrix3d localVerticalToRotatingPlanetocentricFrameTransformationMatrix =
            getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix( longitude, latitude ) ;
    double vernalOffset = computeVernalOffset(
                tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(time) ) ;

    // Compute transformation matrix
    // rotationMatrix = C_TV * C_VR * C_RI * C_IB = Trajectory <- Body
    Eigen::Matrix3d rotationMatrix =
            getLocalVerticalFrameToTrajectoryTransformationMatrix( flightPathAngle , headingAngle ) *
             localVerticalToRotatingPlanetocentricFrameTransformationMatrix.transpose() *
            getInertialToPlanetocentricFrameTransformationMatrix( vernalOffset ) *
            bodyFrameToInertialFrameTransformationQuaternion ;

    // compute aerodynamic angles
    Eigen::Vector3d aerodynamicAngles;
    aerodynamicAngles(0) = std::atan( rotationMatrix(0,2) / rotationMatrix(0,0) );      // Angle of attack
    aerodynamicAngles(1) = std::asin( rotationMatrix(0,1) ) ;                           // Angle of sideslip
    aerodynamicAngles(2) = - std::atan( rotationMatrix(2,1) / rotationMatrix(1,1) ) ;   // bank angle
    return aerodynamicAngles;
}

} // namespace reference_frames
} // namespace tudatApp

