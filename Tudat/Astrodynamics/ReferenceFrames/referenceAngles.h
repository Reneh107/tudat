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

#ifndef TUDAT_COMPUTE_VERNAL_OFFSET_H
#define TUDAT_COMPUTE_VERNAL_OFFSET_H

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace tudat
{

namespace reference_frames
{

//! Compute the angle between the inertial and rotating geocentric reference frame.
/*!
 * This function computes the angle between the inertial and rotating geocentric reference
 * frames, based on the empirically derived relations of USNO (2011). Note that short-term motion
 * of the equinox due to nutation is not considered.
 *
 * Reference:
 *      USNO (2011), Explanatory Supplement to the Astronomical Almanac, The U.S. Naval
 *          Observatory, p. 50.
 *
 */
double computeVernalOffset( double julianDay );



} // namespace reference_frame_transformations
} // namespace tudatApp

#endif // TUDAT_COMPUTE_VERNAL_OFFSET_H
