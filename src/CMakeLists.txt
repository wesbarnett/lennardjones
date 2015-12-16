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

add_executable(ljmd main.cpp NeighborList.cpp PdbFile.cpp
Rdf.cpp ThermodynamicVariable.cpp Thermostat.cpp Velocity.cpp)

find_library(
    LIBGMXCPP gmxcpp
    HINTS ENV LD_LIBRARY_PATH
)

if(${LIBGMXCPP} STREQUAL "LIBGMXCPP-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the gmxcpp library.")
endif()

find_path(
    GMXCPP_INCLUDES
    NAMES "gmxcpp/"
    HINTS ENV CPLUS_INCLUDE_PATH
)

message(STATUS "Found gmxcpp library at: ${LIBGMXCPP}")
message(STATUS "Found gmxcpp headrs at: ${GMXCPP_INCLUDES}")

find_package(Boost)

if(!${Boost_FOUND})
    message(FATAL_ERROR "Boost libraries not found.")
endif()

target_link_libraries(ljmd ${LIBBBOOST} ${LIBGMXCPP})
target_include_directories(ljmd PUBLIC ${BOOST_INCLUDES} ${GMXCPP_INCLUDES})

install (TARGETS ljmd DESTINATION bin)
