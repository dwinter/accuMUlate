# Copyright 2012 Tobias Marschall
# 
# This file is part of CLEVER.
# 
# CLEVER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CLEVER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CLEVER.  If not, see <http://www.gnu.org/licenses/>.


# - Try to find Bamtools (https://github.com/pezmaster31/bamtools)
# Once done this will define
#  Bamtools_FOUND
#  Bamtools_INCLUDE_DIR
#  Bamtools_LIBRARIES
#  Bamtools_DEFINITIONS

set(Bamtools_PREFIX "/usr/" "/usr/local/" "~/lib/" "~/local/lib/" "~/local/" CACHE PATH "Directory Bamtools resides in")
find_path(Bamtools_INCLUDE_DIR api/api_global.h HINTS ${Bamtools_PREFIX}/include PATH_SUFFIXES bamtools)
find_path(Bamtools_LINK_LIBRARY_DIR libbamtools.a HINTS ${Bamtools_PREFIX}/lib/ PATH_SUFFIXES bamtools)

find_library(Bamtools_LIBRARY NAMES libbamtools.a HINTS ${Bamtools_PREFIX}/lib/ PATH_SUFFIXES bamtools)
find_library(Bamtools_LIBRARY_UTILS NAMES libbamtools-utils.a HINTS ${Bamtools_PREFIX}/lib/ PATH_SUFFIXES bamtools)

set(Bamtools_LIBRARIES ${Bamtools_LIBRARIES} ${Bamtools_LIBRARY} ${Bamtools_LIBRARY_UTILS} z)


#message("Debug 
#${Bamtools_INCLUDE_DIR}
#${Bamtools_LINK_LIBRARY_DIR}
#${Bamtools_LIBRARIES}
#")


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set Bamtools_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Bamtools "Bamtools library (https://github.com/pezmaster31/bamtools) not found. If it is in a non-standard place, you have to set the variable Bamtools_PREFIX, for example by adding -DBamtools_PREFIX=<path> to your cmake call" Bamtools_LIBRARIES Bamtools_INCLUDE_DIR)

IF (BAMTOOLS_FOUND)
	set(Bamtools_FOUND TRUE)
ENDIF()

mark_as_advanced(Bamtools_INCLUDE_DIR Bamtools_LIBRARIES)
