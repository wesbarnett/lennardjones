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
 * @brief Header for triclinicbox class
 * @see triclinicbox.h
 */

#include "triclinicbox.h"

triclinicbox::triclinicbox()
{
    this->resize(3);
    for (int i = 0; i < 3; i++)
    {
        this->at(i).resize(3);
    }
}

triclinicbox::triclinicbox(double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
{
    this->resize(3);
    for (int i = 0; i < 3; i++)
    {
        this->at(i).resize(3);
    }
    this->at(X).at(X) = x1;
    this->at(X).at(Y) = x2;
    this->at(X).at(Z) = x3;
    this->at(Y).at(X) = y1;
    this->at(Y).at(Y) = y2;
    this->at(Y).at(Z) = y3;
    this->at(Z).at(X) = z1;
    this->at(Z).at(Y) = z2;
    this->at(Z).at(Z) = z3;
}

triclinicbox::triclinicbox(double x, double y, double z)
{
    this->resize(3);
    for (int i = 0; i < 3; i++)
    {
        this->at(i).resize(3);
    }
    this->at(X).at(X) = x;
    this->at(X).at(Y) = 0.0;
    this->at(X).at(Z) = 0.0;
    this->at(Y).at(X) = 0.0;
    this->at(Y).at(Y) = y;
    this->at(Y).at(Z) = 0.0;
    this->at(Z).at(X) = 0.0;
    this->at(Z).at(Y) = 0.0;
    this->at(Z).at(Z) = z;
}

void triclinicbox::operator=(double rhs)
{
    this->at(X).at(X) = rhs;
    this->at(X).at(Y) = 0.0;
    this->at(X).at(Z) = 0.0;
    this->at(Y).at(X) = 0.0;
    this->at(Y).at(Y) = rhs;
    this->at(Y).at(Z) = 0.0;
    this->at(Z).at(X) = 0.0;
    this->at(Z).at(Y) = 0.0;
    this->at(Z).at(Z) = rhs;
    return;
}

