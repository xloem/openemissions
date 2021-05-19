if(NOT PKG_CONFIG_FOUND)
    INCLUDE(FindPkgConfig)
endif()
PKG_CHECK_MODULES(PC_OPENEMISSIONS openemissions)

FIND_PATH(
    OPENEMISSIONS_INCLUDE_DIRS
    NAMES openemissions/api.h
    HINTS $ENV{OPENEMISSIONS_DIR}/include
        ${PC_OPENEMISSIONS_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    OPENEMISSIONS_LIBRARIES
    NAMES gnuradio-openemissions
    HINTS $ENV{OPENEMISSIONS_DIR}/lib
        ${PC_OPENEMISSIONS_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

include("${CMAKE_CURRENT_LIST_DIR}/openemissionsTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OPENEMISSIONS DEFAULT_MSG OPENEMISSIONS_LIBRARIES OPENEMISSIONS_INCLUDE_DIRS)
MARK_AS_ADVANCED(OPENEMISSIONS_LIBRARIES OPENEMISSIONS_INCLUDE_DIRS)
