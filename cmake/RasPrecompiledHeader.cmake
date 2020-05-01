# Copyright (c) 2017
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

#===============================================================================
# Utility for setting up precompiled header files.
#
# Control variables:
# - PCH_DISABLE: when true, target_implement_pch() does nothing.
#
# Note:
# Precompiled header file support is primarily intened for Windows/MSVC. It is
# not beneificial for gcc under Linux, though good practise on inlcude patterns
# is encouraged for all platforms. GCC setup is also somewhat flaky missing some
# compiler flags (such as -fopenmp)
#
# Recommended usage:
# 1. Create a header file to use to generate precompiled information.
# 2. Create a corresponding cpp file which only includes the header above.
# 3. Include this script in your CMakeLists.txt file.
# 3. Define your target normally: add_target() add_library()
# 4. After defining the target, invoke:
#   target_implement_pch(<target> <pch.h> FORCE)
#
# Force is not required, but is recommended. It ensures the header is included
# automatically for all source files in the target. This means it does not
# need to be explicitly included and that disabling the precompiled header
# reverts to compiling without one and without additional include statements
# for every source file.
#
# In this way, disabling precompiled header generation can be used to verify
# that the public API has all the required forward declarations and include
# statements.
#
# Using FORCE do not add the corresponding <pch.cpp> to the target source files.
# The FORCE option will do so.
#===============================================================================

# Based on https://github.com/larsch/cmake-precompiled-header/blob/master/PrecompiledHeader.cmake
macro(_combine_arguments _variable)
  set(_result "")
  foreach(_element ${${_variable}})
    set(_result "${_result} \"${_element}\"")
  endforeach()
  string(STRIP "${_result}" _result)
  set(${_variable} "${_result}")
endmacro(_combine_arguments)

function(_pch_create_compiler_FLAGS FILENAME TARGET LANG)
  set(_global_compile_flags "${CMAKE_${LANG}_FLAGS}")
  if(CMAKE_BUILD_TYPE)
    string(TOUPPER "CMAKE_${LANG}_FLAGS_${CMAKE_BUILD_TYPE}" _LANG_flags_var_name)
    set(_global_compile_flags "${_global_compile_flags} ${${_LANG_flags_var_name}}")
  endif(CMAKE_BUILD_TYPE)

  # Here we resolve flags controlled by CMake variables. These appear to be applied late in configuration and aren't
  # readily accessible. That makes this section both volatile and fragile.
  set(_behaviouralFlags)
  if(CMAKE_CXX_STANDARD AND LANG STREQUAL CXX)
    list(APPEND _behaviouralFlags "-std=gnu++${CMAKE_CXX_STANDARD}")
  endif(CMAKE_CXX_STANDARD AND LANG STREQUAL CXX)
  if(CMAKE_POSITION_INDEPENDENT_CODE)
    list(APPEND _behaviouralFlags "-fPIC")
  endif(CMAKE_POSITION_INDEPENDENT_CODE)

  if(_behaviouralFlags)
    string(REPLACE ";" "\n" _behaviouralFlags "${_behaviouralFlags}")
  endif(_behaviouralFlags)

  # Still to add:
  # - Directory properties:
  #   - COMPILE_DEFINITIONS
  #   - COMPILE_OPTIONS?
  # - Target:
  #   - LINK_FLAGS[_<CONFIG>]?

  set(_include_directories "$<TARGET_PROPERTY:${TARGET},INCLUDE_DIRECTORIES>")
  set(_compile_definitions "$<TARGET_PROPERTY:${TARGET},COMPILE_DEFINITIONS>")
  set(_compile_flags "$<TARGET_PROPERTY:${TARGET},COMPILE_FLAGS>")
  set(_compile_options "$<TARGET_PROPERTY:${TARGET},COMPILE_OPTIONS>")
  set(_include_directories "$<$<BOOL:${_include_directories}>:-I$<JOIN:${_include_directories},\n-I>\n>")
  set(_compile_definitions "$<$<BOOL:${_compile_definitions}>:-D$<JOIN:${_compile_definitions},\n-D>\n>")
  set(_compile_flags "$<$<BOOL:${_compile_flags}>:$<JOIN:${_compile_flags},\n>\n>")
  set(_compile_options "$<$<BOOL:${_compile_options}>:$<JOIN:${_compile_options},\n>\n>")

  file(GENERATE OUTPUT "${FILENAME}" CONTENT "${_global_compile_flags}\n${_compile_definitions}${_include_directories}${_compile_flags}${_compile_options}${_behaviouralFlags}\n")
endfunction(_pch_create_compiler_FLAGS)

#===============================================================================
# Helper function which applies the use of a precompiled header file. This is
# invoked for a target after defining that target.
#
# Usage:
#   ras_target_implement_pch(<target-name> <header-file> [FORCE] [VM <size>])
# Must appear after the add_library() or add_executable() statement and the
# header file must have a mathing source file already added to the target.
#
# <target-name> : The name of the target.
# <header-file> : The (pre-existing) precompiled header include file.
# FORCE : When present, forces the <header-file> to be included in into the
#   source files rathern than having an explicit include statement.
# VM <size> : VM size for Visual Studio (/Zm<size> option).
# SOURCES : Specify the sources for which the precompiled header file will be
#   included. This is faster for gcc/clang. Ignored for Visual Studio/MSVC.
#
# TODO(KS): Disable all build warnings for the pch compilation.
#===============================================================================
function(ras_target_implement_pch TARGET HEADER)
  CMAKE_PARSE_ARGUMENTS(TIP "FORCE" "VM" "C;SOURCES" "${ARGN}")
  get_filename_component(pchFileName ${HEADER} NAME_WE)
  get_filename_component(pchFilePath ${HEADER} ABSOLUTE)
  get_filename_component(pchFileDir ${HEADER} DIRECTORY)

  if(PCH_DISABLE)
    return()
  endif(PCH_DISABLE)

  if(MSVC)
    get_filename_component(pchHeader ${HEADER} NAME)
    set(PCH_OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/${pchHeader}.pch")
    set(PCH_GENERATE_FLAGS "/Yc${pchHeader} /Fp\"${PCH_OUTPUT_PATH}\"")
    set(PCH_USE_FLAGS "/Yu${HEADER} /Fp\"${PCH_OUTPUT_PATH}\"")

    if(TIP_VM)
      set(PCH_GENERATE_FLAGS "${PCH_GENERATE_FLAGS} /Zm${TIP_VM}")
      set(PCH_USE_FLAGS "${PCH_USE_FLAGS} /Zm${TIP_VM}")
    endif(TIP_VM)

    if(TIP_FORCE)
      set(PCH_USE_FLAGS "${PCH_USE_FLAGS} /FI${HEADER}")
    endif(TIP_FORCE)

    if(NOT TIP_SOURCES)
      get_target_property(TIP_SOURCES ${TARGET} SOURCES)
    endif(NOT TIP_SOURCES)
    foreach(_source ${TIP_SOURCES})
      set(PCH_COMPILE_FLAGS "")
      if(_source MATCHES \\.\(cc|cxx|cpp\)$)
        get_filename_component(_sourceWe ${_source} NAME_WE)
        # Test for PCH source file.
        if(_sourceWe STREQUAL ${pchFileName})
          # Found PCH source file. Set flags to compile to generate PCH file
          set(PCH_COMPILE_FLAGS "${PCH_GENERATE_FLAGS}")
          set(_sourceFound TRUE)
        else(_sourceWe STREQUAL ${pchFileName})
          # General source. Mark file to use PCH file.
          set(PCH_COMPILE_FLAGS "${PCH_USE_FLAGS}")
        endif(_sourceWe STREQUAL ${pchFileName})
      set_source_files_properties(${_source} PROPERTIES COMPILE_FLAGS "${PCH_COMPILE_FLAGS}")
      endif(_source MATCHES \\.\(cc|cxx|cpp\)$)
    endforeach(_source)
    if(NOT _sourceFound)
      set(pchSource)
      set(pchSourceFound FALSE)
      foreach(ext cc cxx cpp)
        if(pchFileDir)
          set(pchSource "${pchFileDir}/${pchFileName}.${ext}")
        else(pchFileDir)
          set(pchSource "${pchFileName}.${ext}")
        endif(pchFileDir)
        if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/${pchSource}")
          set(pchSourceFound TRUE)
          break()
        endif(EXISTS "${CMAKE_CURRENT_LIST_DIR}/${pchSource}")
      endforeach(ext)


      if(pchSourceFound)
        # Source file found, but not part of the target. Add it to the target.
        message(STATUS "Adding precompiled header source for ${TARGET}: ${pchSource}")
        target_sources(${TARGET} PUBLIC "${pchSource}" PUBLIC "${HEADER}")
        set_source_files_properties("${pchSource}" PROPERTIES COMPILE_FLAGS "${PCH_GENERATE_FLAGS}")
      else(pchSourceFound)
        message(FATAL_ERROR "Matching Precompiled Header source file not found for header: ${HEADER}.")
      endif(pchSourceFound)
    endif(NOT _sourceFound)
  endif(MSVC)

  if(CMAKE_COMPILER_IS_GNUCXX)
    # Based on https://github.com/larsch/cmake-precompiled-header/blob/master/PrecompiledHeader.cmake
    # Resolve just the header file name without any prefixed directories.
    get_filename_component(_name ${HEADER} NAME)
    # Resolve the original header as a full path.
    set(_pch_header "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}")
    # set(_pch_header "${CMAKE_CURRENT_LIST_DIR}/${HEADER}")
    set(_pch_binary_dir "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}_pch")
    # set(_pch_marshalled "${_pch_binary_dir}/${HEADER}")
    # set(_pch_marshalled "${CMAKE_CURRENT_LIST_DIR}/${HEADER}")
    # Generate a marshalling path.
    set(_pch_marshalled "${_pch_binary_dir}/${_name}")
    # Artefact directory aligns with the marshalled file.
    set(_outdir "${_pch_binary_dir}/${_name}.gch")
    file(MAKE_DIRECTORY "${_outdir}")
    # Artefacts
    set(_output_cxx "${_outdir}/.c++")
    set(_output_c "${_outdir}/.c")

    if(PCH_DEBUG)
      message("------------------------------------------------------------")
      message("${TARGET}")
      message("${_pch_header} -> ${_pch_marshalled}")
      message(" -> ${_outdir}")
      message("------------------------------------------------------------")
    endif(PCH_DEBUG)

    set(_pch_cxx_flags_file "${_pch_binary_dir}/compile_cxx_flags.rsp")
    _pch_create_compiler_FLAGS("${_pch_cxx_flags_file}" ${TARGET} CXX)

    set(_pch_c_flags_file "${_pch_binary_dir}/compile_c_flags.rsp")
    _pch_create_compiler_FLAGS("${_pch_c_flags_file}" ${TARGET} C)

    set(_compiler_CXX_FLAGS "@${_pch_cxx_flags_file}")
    set(_compiler_C_FLAGS "@{_pch_c_flags_file}")

    # We marshal the pch header so that it lies alongside the precompiled artefact.
    add_custom_command(
      OUTPUT "${_pch_marshalled}"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${_pch_header}" "${_pch_marshalled}"
      MAIN_DEPENDENCY "${_pch_header}"
      COMMENT "Updating ${_name}")
    add_custom_command(
      OUTPUT "${_output_cxx}"
      COMMAND "${CMAKE_CXX_COMPILER}" ${_compiler_CXX_FLAGS} -x c++-header -o "${_output_cxx}" "${_pch_marshalled}"
      DEPENDS "${_pch_marshalled}" "${_pch_cxx_flags_file}"
      COMMENT "Precompiling ${_name} for ${TARGET} (C++)")
    add_custom_command(
      OUTPUT "${_output_c}"
      COMMAND "${CMAKE_C_COMPILER}" ${_compiler_C_FLAGS} -x c-header -o "${_output_c}" "${_pch_marshalled}"
      DEPENDS "${_pch_marshalled}" "${_pch_c_flags_file}"
      COMMENT "Precompiling ${_name} for ${TARGET} (C)")

    get_property(_sources TARGET ${TARGET} PROPERTY SOURCES)
    list(REMOVE_DUPLICATES _sources)
    foreach(_source ${_sources})
      set(_pch_compile_flags "")

      if(_source MATCHES "\\.\(cc|cxx|cpp|c\)$")
        get_source_file_property(_pch_compile_flags "${_source}" COMPILE_FLAGS)
        if(NOT _pch_compile_flags)
          set(_pch_compile_flags)
        endif()
        separate_arguments(_pch_compile_flags)
        list(APPEND _pch_compile_flags "-Winvalid-pch")
        # list(APPEND _pch_compile_flags "-H")
        if(TIP_FORCE)
          list(APPEND _pch_compile_flags -include "${_pch_marshalled}")
        else(TIP_FORCE)
          list(APPEND _pch_compile_flags "-I${_pch_binary_dir}")
        endif(TIP_FORCE)

        get_source_file_property(_object_depends "${_source}" OBJECT_DEPENDS)
        if(NOT _object_depends)
          set(_object_depends)
        endif()
        list(APPEND _object_depends "${_pch_marshalled}")
        if(_source MATCHES \\.\(cc|cxx|cpp\)$)
          list(APPEND _object_depends "${_output_cxx}")
        else()
          list(APPEND _object_depends "${_output_c}")
        endif()

        _combine_arguments(_pch_compile_flags)
        # message("${_pch_compile_flags}")
        set_source_files_properties(${_source} PROPERTIES
          COMPILE_FLAGS "${_pch_compile_flags}"
          OBJECT_DEPENDS "${_object_depends}")
      endif()
    endforeach()
    # elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  #   message("xxxxxxxxxxxxxxxxxxx")
  #   get_filename_component(_name ${HEADER} NAME)

  #   set(_source "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}")
  #   set(_output "${CMAKE_CURRENT_BINARY_DIR}/${_name}.pch")

  #   # Clang precompiled header generation requires a hpp file rather than an h file
  #   # to ensure the header is treated as a C++ header, not a C header.
  #   get_filename_component(_extension ${HEADER} EXT)
  #   if(_extension MATCHES ".[Hh]")
  #     message("mutating to C++ header")
  #     get_filename_component(_nameOnly ${HEADER} NAME_WE)
  #     # C header. Need to copy and generate a C++ header.
  #     set(_source2 ${CMAKE_CURRENT_BINARY_DIR}/${_nameOnly}.hpp)
  #     add_custom_command(
  #       OUTPUT ${_source2}
  #       COMMAND ${CMAKE_COMMAND} -E copy ${_source} ${_source2}
  #       DEPENDS ${_source})
  #     set(_source ${_source2})
  #     message("Using ${_source}")
  #     endif()

  #   if(NOT CMAKE_BUILD_TYPE)
  #     set(_flags_var_name CMAKE_BUILD_TYPE_RELEASE)
  #   else(CMAKE_BUILD_TYPE)
  #    string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" _flags_var_name)
  #   endif(NOT CMAKE_BUILD_TYPE)
  #   message("CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
  #   message("CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
  #   message("${_flags_var_name}: ${${_flags_var_name}}")
  #   message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
  #   set(_compiler_FLAGS ${${_flags_var_name}})
  #   set(_compiler_FLAGS "${CMAKE_CXX_FLAGS} ${${_flags_var_name}}")
  #   message("${_compiler_FLAGS}")

  #   get_directory_property(_directory_flags INCLUDE_DIRECTORIES)
  #   foreach(item ${_directory_flags})
  #     set(_compiler_FLAGS "${_compiler_FLAGS} -I${item}")
  #     message(STATUS "set(_compiler_FLAGS \"${_compiler_FLAGS} -I${item}\")")
  #   endforeach(item)

  #   get_directory_property(_directory_flags DEFINITIONS)
  #   foreach(item ${_directory_flags})
  #     set(_compiler_FLAGS "${_compiler_FLAGS} -I${item}")
  #     message(STATUS "set(_compiler_FLAGS \"${_compiler_FLAGS} -I${item}\")")
  #   endforeach(item)

  #   # message("${CMAKE_CXX_COMPILER} -DPCHCOMPILE ${_compiler_FLAGS} -x c++-header -o {_output} ${_source}")
  #   message(STATUS "
  #   add_custom_command(
  #     OUTPUT ${_output}
  #     COMMAND ${CMAKE_CXX_COMPILER} -cc1 -emit-pch ${_compiler_FLAGS} ${_source} -o ${_output}
  #     DEPENDS ${_source} )
  #   ")
  #   separate_arguments(_compiler_FLAGS)
  #   add_custom_command(
  #     OUTPUT ${_output}
  #     COMMAND ${CMAKE_CXX_COMPILER} -cc1 -emit-pch ${_compiler_FLAGS} ${_source} -o ${_output}
  #     DEPENDS ${_source} )
  #   add_custom_target(${TARGET}_pch DEPENDS ${_output})
  #   add_dependencies(${TARGET} ${TARGET}_pch)
  #   get_target_property(_flags ${TARGET} COMPILE_FLAGS)
  #   if(COMPILE_FLAGS)
  #     list(APPEND _flags "-include-pch ${_output}")
  #   else(COMPILE_FLAGS)
  #     set(_flags "-include-pch ${_output}")
  #   endif(COMPILE_FLAGS)
  #   message("set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS ${_flags})")
  #   set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS ${_flags})
  endif(CMAKE_COMPILER_IS_GNUCXX)
endfunction(ras_target_implement_pch)
