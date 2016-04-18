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

    double julianDatePreviousMidnight = std::floor(julianDate-0.5)+0.5 ;
    double julianDateSinceJ2000 = julianDatePreviousMidnight - JULIAN_DAY_ON_J2000 ;
    double centuriesSince2000 = julianDateSinceJ2000/36525.0;

    // Code: A. Ronse
//    double greenwhichMeanSiderealTime = 24110.54841 +
//            8640184.812866 * centuriesSince2000 +
//            1.00273790935 * (julianDate - julianDatePreviousMidnight) * 86400.0 +
//            0.093104 * centuriesSince2000 * centuriesSince2000 +
//            - centuriesSince2000 * centuriesSince2000 * centuriesSince2000 * 6.2E-6 ; // effect -1.50568e-008

    // Code: report A. Ronse
//    double greenwhichMeanSiderealTime = 24110.54841 +
//            236.55536 * julianDateSinceJ2000 +
//            1.00273790935 * (julianDate - julianDatePreviousMidnight) * 86400.0 +
//            0.093104 * centuriesSince2000 * centuriesSince2000;

    // Code http://aa.usno.navy.mil/faq/docs/GAST.php
    double greenwhichMeanSiderealTime = (6.697374558 +
            0.06570982441908 * julianDateSinceJ2000 +
            1.00273790935 * (julianDate - julianDatePreviousMidnight) * 24.0 +
            0.000026 * centuriesSince2000 * centuriesSince2000) * 3600;

    return 2.0 * tudat::mathematical_constants::PI * computeModulo<double>(greenwhichMeanSiderealTime , 86400.0) / 86400.0;
}

} // namespace reference_frames
} // namespace tudatApp

