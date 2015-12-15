# Copyright (C) 2015 James W. Barnett <jbarnet4@tulane.edu>
#  
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# The full license is located in a text file titled "LICENSE" in the root
# directory of the source.

cmake_minimum_required(VERSION 2.8)
project(ljmd)
set (ljmd_VERSION_MAJOR 1)
set (ljmd_VERSION_MINOR 0)
set(CMAKE_CXX_FLAGS "-std=gnu++11 -Wall -fopenmp")
include_directories(include)
install(DIRECTORY include/ DESTINATION include/ljmd FILES_MATCHING PATTERN "*.h")
add_subdirectory(src)
