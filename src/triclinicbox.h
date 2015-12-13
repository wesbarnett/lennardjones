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
 * @brief Header for triclinicbox class
 */

#ifndef TRICLINICBOX_H
#define TRICLINICBOX_H

#include "coordinates.h"
#include "xdrfile/xdrfile.h"
#include <vector>
using namespace std;

/** @brief Box dimensions.
 * @details This is just a two dimensional vector initialized to three
 * items in each dimension. */
class triclinicbox : public vector < vector <double> > {
public:
/** Constructor, makes the 2d vector 3x3 */
triclinicbox();

/** Constructor where user provides dimensions */
triclinicbox(double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3);

/** Constructor where user provides dimensions, cubic */
triclinicbox(double x, double y, double z);

void operator=(double rhs);
};

#endif
