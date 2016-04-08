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
 *
 *    References
 *
 *    Notes
 *
 */
#include <iostream>
#include "Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h"


//! Tudat library namespace.
namespace tudat
{
namespace aerodynamics
{

void NRLMSISE00Atmosphere::computeProperties(double altitude, double longitude,
                                             double latitude, double time) {
    // Compute the hash key
    size_t hashKey = hashFunc(altitude, longitude, latitude, time);

    // If hash key is same do nothing
    if (hashKey == hashKey_) {
      return;
    }
    hashKey_ = hashKey;

    NRLMSISE00Input inputData = nrlmsise00InputFunction_(
        altitude, longitude, latitude, time);

    std::copy(inputData.apVector.begin(),
              inputData.apVector.end(), aph_.a);
    std::copy(inputData.switches.begin(),
              inputData.switches.end(), flags_.switches);

    input_.g_lat  = latitude;
    input_.g_long = longitude;
    input_.alt    = altitude*1E-3;
    input_.year   = inputData.year;
    input_.doy    = inputData.dayOfTheYear;
    input_.sec    = inputData.secondOfTheDay;
    input_.lst    = inputData.localSolarTime;
    input_.f107   = inputData.f107;
    input_.f107A  = inputData.f107a;
    input_.ap     = inputData.apDaily;
    input_.ap_a   = &aph_;

    gtd7(&input_, &flags_, &output_);

    density_ = output_.d[5]*1000.0; // GM/CM3 to kg/M3
    temperature_ = output_.t[1];
    pressure_ = TUDAT_NAN;
}

//! Overloaded ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                                 NRLMSISE00Input& nrlmsiseInput ){
    stream << "This is a NRLMSISE Input data object." << std::endl;
    stream << "The input data is stored as: " << std::endl;

    stream << "Year              = " << nrlmsiseInput.year << std::endl;
    stream << "Day of the year   = " << nrlmsiseInput.dayOfTheYear << std::endl;
    stream << "Second of the day = " << nrlmsiseInput.secondOfTheDay << std::endl;
    stream << "Local solar time  = " << nrlmsiseInput.localSolarTime << std::endl;
    stream << "f107              = " << nrlmsiseInput.f107 << std::endl;
    stream << "f107a             = " << nrlmsiseInput.f107a << std::endl;
    stream << "apDaily           = " << nrlmsiseInput.apDaily << std::endl;

    for(unsigned int i = 0 ; i < nrlmsiseInput.apVector.size() ; i++){
        stream << "apVector[ " << i << " ]     = " << nrlmsiseInput.apVector[i] << std::endl;
    }

    for(unsigned int i = 0 ; i < nrlmsiseInput.switches.size() ; i++){
        stream << "switches[ " << i << " ]     = " << nrlmsiseInput.switches[i] << std::endl;
    }

    // Return stream.
    return stream;
}

std::pair< std::vector< double >, std::vector< double >>
    NRLMSISE00Atmosphere::getFullOutput(
        double altitude, double longitude,
        double latitude, double time) {
    // Compute the properties
    computeProperties(altitude, longitude, latitude, time);
    std::pair< std::vector< double >, std::vector< double >> output;
    // Copy array members of struct to vectors on the pair.
    output.first  = std::vector< double >(output_.d,
        output_.d + sizeof output_.d / sizeof output_.d[0]);
    output.second = std::vector< double >(output_.t,
        output_.t + sizeof output_.t / sizeof output_.t[0]);
    return output;
}

}  // namespace aerodynamics
}  // namespace tudat
