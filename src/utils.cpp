
/*
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

#include "utils.h"

coordinates pbc(coordinates a, cubicbox box)
{

    a[Z] -= box[Z] * nearbyint(a[Z] / box[Z]);
    a[Y] -= box[Y] * nearbyint(a[Y] / box[Y]);;
    a[X] -= box[X] * nearbyint(a[X] / box[X]);;
    return a;
}

double distance(coordinates a, coordinates b, cubicbox box)
{
    return sqrt(distance2(a, b, box));
}

double distance2(coordinates a, coordinates b, cubicbox box)
{
    coordinates c = a - b;
    c = pbc(c, box);
    return dot(c, c);
}

double dot(coordinates a, coordinates b)
{
    return a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z];
}

double magnitude(coordinates x)
{
    return sqrt(dot(x, x));
}

double volume(cubicbox box)
{
    return box[X] * box[Y] * box[Z];
}
