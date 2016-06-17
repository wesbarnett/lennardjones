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

/**
 * @file
 * @author James W. Barnett jbarnet4@tulane.edu
 * @date December 5, 2014
 * @brief Header for coordinates class
 * @see coordinates.h
 */

#include "gmxcpp/coordinates.h"

coordinates::coordinates(){ }

coordinates::coordinates(double x, double y, double z)
{
    this->r[X] = x;
    this->r[Y] = y;
    this->r[Z] = z;
}

double& coordinates::operator[](int i)
{
    return r[i];
}

const double& coordinates::operator[](int i) const
{
    return r[i];
}

void coordinates::set(double x, double y, double z)
{
    this->r[X] = x;
    this->r[Y] = y;
    this->r[Z] = z;
    return;
}

coordinates coordinates::operator-(coordinates rhs)
{
    return (coordinates (r[X] - rhs[X], r[Y] - rhs[Y], r[Z] - rhs[Z]));
}

void coordinates::operator-=(coordinates rhs)
{
    r[X] -= rhs[X];
    r[Y] -= rhs[Y];
    r[Z] -= rhs[Z];
    return;
}

coordinates coordinates::operator+(coordinates rhs)
{
    return (coordinates (r[X] + rhs[X], r[Y] + rhs[Y], r[Z] + rhs[Z]));
}

void coordinates::operator+=(coordinates rhs)
{
    r[X] += rhs[X];
    r[Y] += rhs[Y];
    r[Z] += rhs[Z];
    return;
}

coordinates coordinates::operator/(double rhs)
{
    return (coordinates (r[X] / rhs, r[Y] / rhs, r[Z] / rhs));
}

void coordinates::operator/=(double rhs)
{
    r[X] /= rhs;
    r[Y] /= rhs;
    r[Z] /= rhs;
    return;
}

coordinates operator*(coordinates lhs, double rhs)
{
    return (coordinates (lhs[X] * rhs, lhs[Y] * rhs, lhs[Z] * rhs));
}

coordinates operator*(double lhs, coordinates rhs)
{
    return (coordinates (rhs[X] * lhs, rhs[Y] * lhs, rhs[Z] * lhs));
}

void coordinates::operator*=(double rhs)
{
    r[X] *= rhs;
    r[Y] *= rhs;
    r[Z] *= rhs;
    return;
}

void coordinates::operator=(double rhs)
{
    r[X] = rhs;
    r[Y] = rhs;
    r[Z] = rhs;
    return;
}

