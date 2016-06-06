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
 *      160418    R. Hoogendoorn    File created.
 *
 *    References
 *
 *
 *    Notes
 *      The reference frame definitions/abbreviations can be found in the file
 *      referenceFrameTransformations.h.
 *
 */

#define BOOST_TEST_MAIN

#include <iostream> // cout
#include <iomanip> // precision

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/ReferenceFrames/referenceAngles.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

#include <Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

double vernalOffsetRonse(double julianDay){
    double GreenwhichSiderealTime;
    double Tu;
    double fractpart, intpart;

    fractpart = std::modf( julianDay + 0.5, &intpart );
    julianDay = julianDay - fractpart;
    Tu = (julianDay - 2451545.0)/36525;

    GreenwhichSiderealTime = 24110.54841 + Tu*(8640184.812866
                                       + Tu * (0.093104 - Tu * 6.2E-6));

    GreenwhichSiderealTime = tudat::basic_mathematics::computeModulo(
            GreenwhichSiderealTime + 86400.0*1.00273790934*fractpart, 86400.0);

    return 2.0* tudat::mathematical_constants::PI * GreenwhichSiderealTime / 86400.0;
}

namespace tudat
{
namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_vernal_offset )

using tudat::mathematical_constants::PI;

//! Test vernal offset computation
BOOST_AUTO_TEST_CASE( test_compute_vernal_offset_compare_Ronse)
{
    // Define julian date and compute vernal offset
    double julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(2013,6,11,6,25,15);
    double vernalOffset = tudat::reference_frames::computeVernalOffset(julianDate) ;

    // Comparison with matlab function: Julian date to greenwich mean sidereal time (File exchange)
    // http://www.mathworks.com/matlabcentral/fileexchange/28176-julian-date-to-greenwich-mean-sidereal-time
    BOOST_CHECK_CLOSE_FRACTION( vernalOffset*180.0/PI , 356.0725571517477 , 3E-8 );

    //! Second julian date
    julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(2015,3,11,15,25,15);
    vernalOffset = tudat::reference_frames::computeVernalOffset(julianDate) ;

    // Comparison with matlab function
    BOOST_CHECK_CLOSE_FRACTION( vernalOffset*180.0/PI , 40.28519655625729 , 3E-7 );

    //! Third julian date
    julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(2002,9,5,12,15,0);
    vernalOffset = tudat::reference_frames::computeVernalOffset(julianDate) ;

    // Comparison with matlab function
    BOOST_CHECK_CLOSE_FRACTION( vernalOffset*180.0/PI , 168.1840220843993 , 1E-7 );
}

//! test aerodynamic angles computation: aerodynamic angles zero
BOOST_AUTO_TEST_CASE( test_compute_aerodynamic_angles_1 )
{
    // Aerodynamic angles zero
    double longitude        = 0.0 ;
    double latitude         = 0.0 ;
    double flightPathAngle  = 0.0 ;
    double headingAngle     = 0.0 ; // Flying towards the north
    double time             = 0.0;

    // Compute vernal offset
    double vernalOffset = tudat::reference_frames::computeVernalOffset(
                tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(time) ) ;

    // Inertial to rotating frame transformation matrix
    Eigen::Matrix3d InertialToRotating =
            tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix( vernalOffset );

    // Define attitude of body
    // X-axis towards north (direction of flight) , Y-axis towards east ,Z-axis towards earth
    Eigen::Matrix3d bodyToInertial; // Rotation matrix
    Eigen::Vector3d rotationAxis;
    rotationAxis << 0.0 , 1.0 , 0.0 ;
    double rotationAngle = PI/2.0 ; // 90 deg around y axis
    bodyToInertial = Eigen::AngleAxisd(- rotationAngle , rotationAxis.normalized() ) ;
    bodyToInertial = InertialToRotating.transpose() * bodyToInertial ;

    Eigen::Vector3d aerodynamicAngles = tudat::reference_frames::computeAerodynamicAngles(
                longitude,latitude,time,flightPathAngle,headingAngle,bodyToInertial);

    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(0) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(1) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(2) - 0.0 ) , 1E-15 );
}

//! test aerodynamic angles computation: Only angle of attack
BOOST_AUTO_TEST_CASE( test_compute_aerodynamic_angles_2 )
{
    // Only angle of attack
    double longitude        = 0.0 ;
    double latitude         = 0.0 ;
    double flightPathAngle  = 0.0 ;
    double headingAngle     = 0.0 ; // Flying towards the north
    double time             = 0.0;

    // Compute vernal offset
    double vernalOffset = tudat::reference_frames::computeVernalOffset(
                tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(time) ) ;

    // Inertial to rotating frame transformation matrix
    Eigen::Matrix3d InertialToRotating =
            tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix( vernalOffset );

    // Define attitude of body
    // X-axis towards north (direction of flight) , Y-axis towards east ,Z-axis towards earth
    Eigen::Matrix3d bodyToInertial; // Rotation matrix
    Eigen::Vector3d rotationAxis;
    rotationAxis << 0.0 , 1.0 , 0.0 ;
    double angleOfAttack = PI/8.0 ;
    double rotationAngle = PI/2.0 - angleOfAttack ; // 90 deg around y axis - angle of attack
    bodyToInertial = Eigen::AngleAxisd(- rotationAngle , rotationAxis.normalized() ) ;
    bodyToInertial = InertialToRotating.transpose() * bodyToInertial ;

    Eigen::Vector3d aerodynamicAngles = tudat::reference_frames::computeAerodynamicAngles(
                longitude,latitude,time,flightPathAngle,headingAngle,bodyToInertial);

    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(0) - angleOfAttack ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(1) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(2) - 0.0 ) , 1E-15 );
}

//! test aerodynamic angles computation: Only bank sideslip
BOOST_AUTO_TEST_CASE( test_compute_aerodynamic_angles_3 )
{
    // Only angle of sideslip
    double longitude        = 0.0 ;
    double latitude         = 0.0 ;
    double flightPathAngle  = 0.0 ;
    double headingAngle     = 0.0 ; // Flying towards the north
    double time             = 0.0;

    // Compute vernal offset
    double vernalOffset = tudat::reference_frames::computeVernalOffset(
                tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(time) ) ;

    // Inertial to rotating frame transformation matrix
    Eigen::Matrix3d InertialToRotating =
            tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix( vernalOffset );

    // Define attitude of body
    // X-axis towards north (direction of flight) , Y-axis towards east ,Z-axis towards earth
    Eigen::Matrix3d bodyToInertial; // Rotation matrix
    Eigen::Vector3d rotationAxis;
    rotationAxis << 0.0 , 1.0 , 0.0 ;
    double rotationAngle = PI/2.0 ; // 90 deg around y axis
    bodyToInertial = Eigen::AngleAxisd(- rotationAngle , rotationAxis.normalized() ) ;

    double angleOfSideslip = -PI/6.0 ;
    rotationAxis << 0.0 , 0.0 , 1.0 ;
    Eigen::Matrix3d frameRotation;
    frameRotation = Eigen::AngleAxisd(- angleOfSideslip , rotationAxis.normalized() ) ;

    bodyToInertial = bodyToInertial * frameRotation;
    bodyToInertial = InertialToRotating.transpose() * bodyToInertial ;

    Eigen::Vector3d aerodynamicAngles = tudat::reference_frames::computeAerodynamicAngles(
                longitude,latitude,time,flightPathAngle,headingAngle,bodyToInertial);

    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(0) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(1) - angleOfSideslip ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(2) - 0.0 ) , 1E-15 );
}

//! test aerodynamic angles computation: Only bank angle
BOOST_AUTO_TEST_CASE( test_compute_aerodynamic_angles_4 )
{
    // Only bank angle
    double longitude        = 0.0 ;
    double latitude         = 0.0 ;
    double flightPathAngle  = 0.0 ;
    double headingAngle     = 0.0 ; // Flying towards the north
    double time             = 0.0;

    // Compute vernal offset
    double vernalOffset = tudat::reference_frames::computeVernalOffset(
                tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(time) ) ;

    // Inertial to rotating frame transformation matrix
    Eigen::Matrix3d InertialToRotating =
            tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix( vernalOffset );

    // Define attitude of body
    // X-axis towards north (direction of flight) , Y-axis towards east ,Z-axis towards earth
    Eigen::Matrix3d bodyToInertial; // Rotation matrix
    Eigen::Vector3d rotationAxis;
    rotationAxis << 0.0 , 1.0 , 0.0 ;
    double rotationAngle = PI/2.0 ; // 90 deg around y axis
    bodyToInertial = Eigen::AngleAxisd(- rotationAngle , rotationAxis.normalized() ) ;

    double bankAngle = PI/6.0 ;
    rotationAxis << 1.0 , 0.0 , 0.0 ;
    Eigen::Matrix3d frameRotation;
    frameRotation = Eigen::AngleAxisd(- bankAngle , rotationAxis.normalized() ) ;

    bodyToInertial = bodyToInertial * frameRotation;
    bodyToInertial = InertialToRotating.transpose() * bodyToInertial ;

    Eigen::Vector3d aerodynamicAngles = tudat::reference_frames::computeAerodynamicAngles(
                longitude,latitude,time,flightPathAngle,headingAngle,bodyToInertial);

    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(0) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(1) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(2) - bankAngle ) , 1E-15 );
}

////! test aerodynamic angles computation: Angle of attack and sideslip
//BOOST_AUTO_TEST_CASE( test_compute_aerodynamic_angles_5 )
//{
//    // Only angle of attack and sideslip
//    double longitude        = 0.0 ;
//    double latitude         = 0.0 ;
//    double flightPathAngle  = 0.0 ;
//    double headingAngle     = 0.0 ; // Flying towards the north
//    double time             = 0.0;

//    // Compute vernal offset
//    double vernalOffset = tudat::reference_frames::computeVernalOffset(
//                tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(time) ) ;

//    // Inertial to rotating frame transformation matrix
//    Eigen::Matrix3d InertialToRotating =
//            tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix( vernalOffset );

//    // Define attitude of body
//    // X-axis towards north (direction of flight) , Y-axis towards east ,Z-axis towards earth
//    Eigen::Matrix3d bodyToInertial; // Rotation matrix
//    Eigen::Vector3d rotationAxis;
//    rotationAxis << 0.0 , 1.0 , 0.0 ;
//    double angleOfAttack = PI/8.0 ;
//    double rotationAngle = PI/2.0 - angleOfAttack ; // 90 deg around y axis - angle of attack
//    bodyToInertial = Eigen::AngleAxisd(- rotationAngle , rotationAxis.normalized() ) ;

//    double angleOfSideslip = -PI/6.0 ;
//    rotationAxis << 0.0 , 0.0 , 1.0 ;
//    Eigen::Matrix3d frameRotation;
//    frameRotation = Eigen::AngleAxisd(- angleOfSideslip , rotationAxis.normalized() ) ;

//    bodyToInertial = bodyToInertial * frameRotation;
//    bodyToInertial = InertialToRotating.transpose() * bodyToInertial ;

//    Eigen::Vector3d aerodynamicAngles = tudat::reference_frames::computeAerodynamicAngles(
//                longitude,latitude,time,flightPathAngle,headingAngle,bodyToInertial);

//    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(0) - angleOfAttack ) , 1E-15 );
//    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(1) - angleOfSideslip ) , 1E-15 );
//    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(2) - 0.0 ) , 1E-15 );
//}

//! test aerodynamic angles computation: Body frame same orientation as rotating planetocentric frame
BOOST_AUTO_TEST_CASE( test_compute_aerodynamic_angles_6 )
{
    // Body frame same orientation as rotating planetocentric frame
    double longitude        = 0.0 ;
    double latitude         = 0.0 ;
    double flightPathAngle  = 0.0 ;
    double headingAngle     = 0.0 ; // Flying towards the north
    double time             = 0.0 ;

    // Compute vernal offset
    double vernalOffset = tudat::reference_frames::computeVernalOffset(
                tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(time) ) ;

    // Inertial to rotating frame transformation matrix
    Eigen::Matrix3d InertialToRotating =
            tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix( vernalOffset );

    // Define attitude of body
    // No rotation = orientation B and I frame identical
    Eigen::Matrix3d noRotation;
    noRotation = Eigen::Matrix3d::Identity(3,3);
    noRotation = InertialToRotating.transpose() * noRotation ;

    Eigen::Vector3d aerodynamicAngles = tudat::reference_frames::computeAerodynamicAngles(
                    longitude,latitude,time,flightPathAngle,headingAngle,noRotation);

    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(0) - PI/2.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(1) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(2) - 0.0 ) , 1E-15 );
}

////! test aerodynamic angles computation
//BOOST_AUTO_TEST_CASE( test_compute_aerodynamic_angles_7 )
//{
//    // Body frame same orientation as rotating planetocentric frame
//    double longitude        = 0.0 ;
//    double latitude         = 0.0 ;
//    double flightPathAngle  = 0.0 ;
//    double headingAngle     = 0.0 ; // Flying towards the north
//    double time             = 0.0;

//    // Compute vernal offset
//    double vernalOffset = tudat::reference_frames::computeVernalOffset(
//                tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(time) ) ;

//    // Inertial to rotating frame transformation matrix
//    Eigen::Matrix3d InertialToRotating =
//            tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix( vernalOffset );

//    // Define attitude of body
//    double rotationAngle = 0.0 ;
//    rotationAxis(0) = 0.0 ;
//    rotationAxis(1) = 0.0 ;
//    rotationAxis(2) = 1.0 ;
//    //    rotationAxis = rotationAxis / ( std::sin(rotationAngle/2.0) ) ;
//    Eigen::Quaterniond BodyFrameToInertialFrameTransformationMatrix;
//    BodyFrameToInertialFrameTransformationMatrix = Eigen::AngleAxisd(- rotationAngle , rotationAxis.normalized() ) ;

//    Eigen::Matrix<double,4,1> quaternion;
//    quaternion << 0.0 , 0.0 , 0.0 , 1.0 ;
//    Eigen::Matrix3d BodyFrameToInertialFrameTransformationMatrix;
//    BodyFrameToInertialFrameTransformationMatrix = tudat::reference_frames::convertQuaternionToRotationMatrix(quaternion);
//    Eigen::Matrix3d C_RI =
//    tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix(
//                tudat::reference_frames::computeVernalOffset(
//                    tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000));
//    //    Eigen::Matrix3d C_VR =

//    Eigen::Vector3d aerodynamicAngles = tudat::reference_frames::computeAerodynamicAngles(
//                    longitude,latitude,time,flightPathAngle,headingAngle,BodyFrameToInertialFrameTransformationMatrix);

//    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(0) - 0.0 ) , 1E-15 );
//    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(1) - 0.0 ) , 1E-15 );
//    BOOST_CHECK_SMALL( std::fabs( aerodynamicAngles(2) - 0.0 ) , 1E-15 );
//}

BOOST_AUTO_TEST_SUITE_END( )


} // namespace unit_tests
} // namespace tudat
