# - Try to find CIMG
# Once done, this will define
#
# CIMG_FOUND - system has CIMG
# CIMG_INCLUDE_DIRS - the CIMG include directories
IF( CIMG_FOUND )
# in cache already
SET( CIMG_FIND_QUIETLY TRUE )
ENDIF()
include(LibFindMacros)
# Use pkg-config to get hints about paths
libfind_pkg_check_modules(CIMG_PKGCONF CIMG)
# Include dir
find_path(CImg_INCLUDE_DIR
NAMES CImg.h
PATHS ${CIMG_PKGCONF_INCLUDE_DIRS}
)
set(CIMG_PROCESS_INCLUDES CImg_INCLUDE_DIR)
IF(UNIX)
set(CIMG_PROCESS_LIBS pthread)
ELSE(UNIX)
find_library(gdi32_LIBRARY
NAMES gdi32
PATHS ${CIMG_PKGCONF_LIBRARY_DIRS}
)
find_library(user32_LIBRARY
NAMES user32
PATHS ${CIMG_PKGCONF_LIBRARY_DIRS}
)
find_library(shell32_LIBRARY
NAMES shell32
PATHS ${CIMG_PKGCONF_LIBRARY_DIRS}
)
set(CIMG_PROCESS_LIBS gdi32_LIBRARY user32_LIBRARY shell32_LIBRARY)
ENDIF(UNIX)
libfind_process(CIMG)
