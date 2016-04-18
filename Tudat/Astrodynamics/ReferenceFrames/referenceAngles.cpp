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

namespace tudat
{

namespace reference_frames
{

//! Compute the angle between the inertial and rotating geocentric reference frame.
double computeVernalOffset( double julianDay )
{
    double greenwhichSiderealTime;
    double Tu;
    double fractpart, intpart; // fractpart: part of day , intpart: day

    fractpart = std::modf( julianDay + 0.5, &intpart );
    julianDay = julianDay - fractpart;
    Tu = (julianDay - 2451545.0)/36525.0;

    greenwhichSiderealTime = 24110.54841 + Tu*(8640184.812866
                                       + Tu * (0.093104 - Tu * 6.2E-6));

    greenwhichSiderealTime = tudat::basic_mathematics::computeModulo(
            greenwhichSiderealTime + 86400.0*1.00273790934*fractpart, 86400.0);

    return 2.0*tudat::mathematics::PI * greenwhichSiderealTime / 86400.0;
}

} // namespace reference_frames
} // namespace tudatApp

