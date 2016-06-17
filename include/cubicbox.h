/*
 * libgmxcpp
 * Copyright (C) 2015 James W. Barnett <jbarnet4@tulane.edu>
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * The full license is located in a text file titled "LICENSE" in the root
 * directory of the source.
 *
 */

/** @file
 * @author James W. Barnett jbarnet4@tulane.edu
 * @date December 5, 2014
 * @brief Header for cubicbox class
 */

#ifndef CUBICBOX_H
#define CUBICBOX_H

#include "gmxcpp/coordinates.h"
#include "gmxcpp/xdrfile.h"
#include <vector>
using namespace std;

/** @brief Box dimensions.
 * @details This is just a two dimensional array initialized to three
 * items in each dimension. To access the elements of the array use operator().
 * For example, to if the box is cubic and your have a cubicbox object named
 * mybox, to get the X dimension do mybox(0). If you it is truly a cubicbox
 * (not cubic) you can access elements with mybox(i,j).*/
class cubicbox {

private:

    array <float,3> box;

public:
    /** Constructor, makes the 2d vector 3x3 */
    cubicbox();

    /** Constructor where user provides dimensions, cubic */
    cubicbox(float x, float y, float z);

    float& operator[](int i);

    const float& operator[](int i) const;

};

#endif
