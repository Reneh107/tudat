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

#include "Tudat/Basics/testMacros.h"

#include "tudat/Astrodynamics/ReferenceFrames/referenceAngles.h"
#include "tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include <tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_vernal_offset )

BOOST_AUTO_TEST_CASE( test_compute_vernal_offset_compare_Ronse)
{
    // Do something
    // Vernal offset 2pi after 1 day
    double julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(2013,06,11,6,25,15);
    double vernalOffset = tudat::reference_frames::computeVernalOffset(julianDate) ;

    double julianDay = julianDate;
    double thetaG_JulianDay;
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

    thetaG_JulianDay = 2.0* tudat::mathematical_constants::PI * GreenwhichSiderealTime / 86400.0;

    BOOST_CHECK_SMALL( std::fabs( vernalOffset - thetaG_JulianDay ) , 1E-10 );

    // Difference in impact location due to error
    std::cout << "difference between ronse / rene = " << (vernalOffset - thetaG_JulianDay)*6378.136E3 << " meters " << std::endl;
    // Difference between Ronse code vs program Ronse = -18.00945692896082 meters
}

BOOST_AUTO_TEST_SUITE_END( )


} // namespace unit_tests
} // namespace tudat
