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
 *      USNO (2011), Explanatory Supplement to the Astronomical Almanac, The U.S. Naval
 *          Observatory, p. 50.
 *
 */

#ifndef TUDAT_REFERENCE_ANGLES_H
#define TUDAT_REFERENCE_ANGLES_H

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace tudat
{

namespace reference_frames
{

//! Compute the angle between the inertial and rotating geocentric reference frame.
/*!
 * This function computes the angle from the inertial (J2000) frame to the rotating geocentric reference
 * frame using the Greenwich Mean Sidereal Time(GMST). Note that short-term motion
 * of the equinox due to nutation is not considered.
 *
 * References:
 *      GMST: Montenbruck, O. & Gill, E. Satellite Orbits Models, Methods and Applications, Springer, 2001 (page 167)
 *      Angle computation: https://celestrak.com/columns/v02n02/
 *
 * \param julianDate julian date at time of interest
 * \return vernal offset in radians
 */
double computeVernalOffset( double julianDate );


//! Compute aerodynamic angles
/*!
 * This function computes the aerodynamic angles.
 *
 * Note:
 *      This function does not include the difference between airspeed and groundspeed based aerodynamic angles
 *
 * Reference:
 *      Mulder, J.; van Staveren, W.; van der Vaart, J.; de Weerdt, E.; de Visser, C.; in ’t Veld, A. & Mooij, E. Flight Dynamics - Lecture Notes 2013
 *
 * \param longitude
 * \param latitude
 * \param time in seconds since the J2000 epoch
 * \param flightPathAngle
 * \param headingAngle
 * \param Body to inertial frame transformation matrix
 * \return aerodynamic angles = angle of attack, angle of sideslip , bank angle in radians
 */
Eigen::Vector3d computeAerodynamicAngles(double longitude,
                         double latitude,
                         double time,
                         double flightPathAngle,
                         double headingAngle,
                         Eigen::Matrix3d bodyFrameToInertialFrameTransformationMatrix);


//! Compute aerodynamic angles
/*!
 * This function computes the aerodynamic angles.
 *
 * Note:
 *      This function does not include the difference between airspeed and groundspeed based aerodynamic angles
 *
 * Reference:
 *      Mulder, J.; van Staveren, W.; van der Vaart, J.; de Weerdt, E.; de Visser, C.; in ’t Veld, A. & Mooij, E. Flight Dynamics - Lecture Notes 2013
 *
 * \param longitude
 * \param latitude
 * \param time in seconds since the J2000 epoch
 * \param flightPathAngle
 * \param headingAngle
 * \param Body to inertial frame transformation quaternion
 * \return aerodynamic angles = angle of attack, angle of sideslip , bank angle in radians
 */
Eigen::Vector3d computeAerodynamicAngles(double longitude,
                         double latitude,
                         double time,
                         double flightPathAngle,
                         double headingAngle,
                         Eigen::Quaterniond bodyFrameToInertialFrameTransformationQuaternion);

} // namespace reference_frame_transformations
} // namespace tudatApp

#endif // TUDAT_REFERENCE_ANGLES_H
