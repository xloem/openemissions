if(DEFINED INCLUDED_FREESDR_UTIL_CMAKE)
  return()
endif()
set(INCLUDED_FREESDR_UTIL_CMAKE TRUE)

function(FREESDR_MODULE_UTIL)
  find_package(Pothos CONFIG)

  # todo: generate docs using https://gitlab.kitware.com/cmake/community/wikis/FAQ#how-can-i-generate-a-source-file-during-the-build
  #  will need to make an executable that enumerates blocks and outputs headerfile

  # TODO !!: only include pothos if detected

  include(CMakeParseArguments)
  CMAKE_PARSE_ARGUMENTS(FREESDR_MODULE_UTIL "NAME;CATEGORY" "SOURCES;LIBRARIES" ${ARGN})

  POTHOS_MODULE_util(
    TARGET ${FREESDR_MODULE_UTIL_NAME}
    DESTINATION ${CATEGORY}
    SOURCES ${FREESDR_MODULE_UTIL_SOURCES}
    LIBRARIES freesdr ${FREESDR_MODULE_UTIL_LIBRARIES}
  )
endfunction(FREESDR_MODULE_UTIL)
