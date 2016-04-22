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
 *      120701    A. Ronse          Creation of code.
 *      160324    R. Hoogendoorn    Update for use in current version of Tudat
 *
 *    References
 *      Data source for validation figures:
 *                        http://celestrak.com/SpaceData/sw19571001.txt
 *
 */

#define BOOST_TEST_MAIN

#include <istream>
#include <string>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/InputOutput/parseSolarActivityData.h"
#include "Tudat/InputOutput/extractSolarActivityData.h"
#include "Tudat/InputOutput/solarActivityData.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include <tudat/Basics/testMacros.h>

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_Solar_Activity_Parser_Extractor )

//! Test parsing and extraction process of solar dummy activity file
BOOST_AUTO_TEST_CASE( test_parsing_and_extraction )
{
tudat::input_output::parsed_data_vector_utilities::ParsedDataVectorPtr parsedDataVectorPtr;

tudat::input_output::solar_activity::ParseSolarActivityData solarActivityParser;
tudat::input_output::solar_activity::ExtractSolarActivityData solarActivityExtractor;

// Parse file
// save path of cpp file
std::string cppPath( __FILE__ );

// Strip filename from temporary string and return root-path string.
std::string folder = cppPath.substr( 0, cppPath.find_last_of("/\\")+1);
std::string filePath = folder + "testSolarActivity.txt" ;

// Open dataFile
std::ifstream dataFile;
dataFile.open(filePath.c_str(), std::ifstream::in);
parsedDataVectorPtr = solarActivityParser.parse( dataFile );
dataFile.close();

// Extract data to object of solarActivityData class
std::vector< boost::shared_ptr<tudat::input_output::solar_activity::SolarActivityData> >
        solarActivityData( 6 );
solarActivityData[0] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 0 ) );
solarActivityData[1] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 3 ) );
solarActivityData[2] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 8 ) );
solarActivityData[3] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 9 ) );
solarActivityData[4] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 12 ) );
solarActivityData[5] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 15 ) );

// Check single parameters
BOOST_CHECK_EQUAL(  solarActivityData[0]->dataType, 1 );
BOOST_CHECK_EQUAL(  solarActivityData[1]->planetaryDailyCharacterFigure, 0.0 );
BOOST_CHECK_EQUAL(  solarActivityData[1]->solarRadioFlux107Adjusted, 103.0 );
BOOST_CHECK_EQUAL(  solarActivityData[2]->month, 5 );
BOOST_CHECK_EQUAL(  solarActivityData[2]->internationalSunspotNumber, 0 );
BOOST_CHECK_EQUAL(  solarActivityData[3]->planetaryRangeIndexSum, 176 );
BOOST_CHECK_EQUAL(  solarActivityData[3]->planetaryDailyCharacterFigureConverted, 0 );
BOOST_CHECK_EQUAL(  solarActivityData[4]->bartelsSolarRotationNumber, 2441 );
BOOST_CHECK_EQUAL(  solarActivityData[5]->last81DaySolarRadioFlux107Observed, 187.8 );

// Check planetaryEquivalentAmplitudeVectors
Eigen::VectorXd correctPlanetaryEquivalentAmplitudeVector1( 8 );
correctPlanetaryEquivalentAmplitudeVector1 << 32, 27, 15, 7, 22, 9, 32, 22;
Eigen::VectorXd correctPlanetaryEquivalentAmplitudeVector2 = Eigen::VectorXd::Zero( 8 );

TUDAT_CHECK_MATRIX_CLOSE(  solarActivityData[0]->planetaryEquivalentAmplitudeVector,
                    correctPlanetaryEquivalentAmplitudeVector1, 0 );
TUDAT_CHECK_MATRIX_CLOSE(  solarActivityData[4]->planetaryEquivalentAmplitudeVector,
                    correctPlanetaryEquivalentAmplitudeVector2 , 0);
TUDAT_CHECK_MATRIX_CLOSE(  solarActivityData[5]->planetaryEquivalentAmplitudeVector,
                    correctPlanetaryEquivalentAmplitudeVector2, 0 );
}

BOOST_AUTO_TEST_CASE( test_function_readSolarActivityData ){
    using tudat::input_output::solar_activity::SolarActivityDataMap ;
    using tudat::input_output::solar_activity::SolarActivityData ;

    // Parse file
    // save path of cpp file
    std::string cppPath( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string folder = cppPath.substr( 0, cppPath.find_last_of("/\\")+1);
    std::string filePath = folder + "testSolarActivity.txt" ;

    SolarActivityDataMap SolarActivity = tudat::input_output::solar_activity::readSolarActivityData(filePath) ;

    double JulianDate ;
    JulianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                1957,
                10,
                3,
                0, 0, 0.0) ;

    SolarActivityDataMap::iterator it;
    it = SolarActivity.find(JulianDate) ;

    BOOST_CHECK_EQUAL(  it->second->day , 3 );
    BOOST_CHECK_EQUAL(  it->second->centered81DaySolarRadioFlux107Observed , 268.1 ) ;
    BOOST_CHECK_EQUAL(  it->second->planetaryEquivalentAmplitudeAverage , 19 ) ;
    BOOST_CHECK_EQUAL(  it->second->planetaryRangeIndexSum , 250 ) ;

    JulianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                2023,
                1,
                1,
                0, 0, 0.0) ;

    it = SolarActivity.find(JulianDate) ;

    BOOST_CHECK_EQUAL(  it->second->day , 1 );
    BOOST_CHECK_EQUAL(  it->second->bartelsSolarRotationNumber , 2583 ) ;
    BOOST_CHECK_EQUAL(  it->second->centered81DaySolarRadioFlux107Observed , 0.0 ) ;
    BOOST_CHECK_EQUAL(  it->second->last81DaySolarRadioFlux107Adjusted , 0.0 ) ;
    BOOST_CHECK_EQUAL(  it->second->planetaryEquivalentAmplitudeAverage , 0.0 ) ;
    BOOST_CHECK_EQUAL(  it->second->planetaryRangeIndexSum , 0.0 ) ;

    std::string filePath2 = folder + "sw19571001.txt" ;

    SolarActivityDataMap SolarActivity2 = tudat::input_output::solar_activity::readSolarActivityData(filePath2) ;

    JulianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                1993,
                12,
                11,
                0, 0, 0.0) ;

    it = SolarActivity2.find(JulianDate) ;

    BOOST_CHECK_EQUAL(  it->second->day , 11 );
    BOOST_CHECK_EQUAL(  it->second->bartelsSolarRotationNumber , 2190 ) ;
    BOOST_CHECK_EQUAL(  it->second->centered81DaySolarRadioFlux107Observed , 104.0 ) ;
    BOOST_CHECK_EQUAL(  it->second->last81DaySolarRadioFlux107Adjusted , 97.6 ) ;
    BOOST_CHECK_EQUAL(  it->second->planetaryEquivalentAmplitudeAverage , 7 ) ;
    BOOST_CHECK_EQUAL(  it->second->planetaryRangeIndexSum , 143 ) ;
    BOOST_CHECK_EQUAL(  it->second->planetaryDailyCharacterFigure , 0.4 ) ;
    BOOST_CHECK_EQUAL(  it->second->planetaryDailyCharacterFigureConverted , 2 ) ;
    BOOST_CHECK_EQUAL(  it->second->internationalSunspotNumber , 31 ) ;
}

BOOST_AUTO_TEST_CASE( test_function_readSolarActivityData_incorrect_filepath ){
    using tudat::input_output::solar_activity::SolarActivityDataMap ;
    using tudat::input_output::solar_activity::SolarActivityData ;

    // Parse file
    // save path of cpp file
    std::string cppPath( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string folder = cppPath.substr( 0, cppPath.find_last_of("/\\")+1);
    std::string filePath = folder + "incorrectFileName.txt" ;

    SolarActivityDataMap SolarActivity = tudat::input_output::solar_activity::readSolarActivityData(filePath) ;
    int isDataMapEmpty; // 1 is empty
    if( SolarActivity.empty() ){
        isDataMapEmpty = 1 ;
    }
    else{
        isDataMapEmpty = 0 ;
    }
    BOOST_CHECK_EQUAL( isDataMapEmpty , 1 );
}

BOOST_AUTO_TEST_SUITE_END( )

}   // unit_tests
}   // tudat
