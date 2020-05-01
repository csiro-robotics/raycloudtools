# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

# if(UTILS_DEFINED)
#   return()
# endif(UTILS_DEFINED)
# set(UTILS_DEFINED TRUE)
include(CMakeParseArguments)

# A helper function for splitting a single list into multiple lists based on the presence of modifier directives in
# the list itself.
#
# For example, the list LIST=[a b c NUMBERS 1 2 3] can be split by calling
# ras_split_list(SPLIT_LIST "NUMBERS" ${LIST})
# This results in SPLIT_LIST=[a b c] and SPLIT_LIST_NUMBERS=[1 2 3]
function(ras_split_list BASE_VAR MODIFIERS)
  set(CURRENT_MODIFIER)

  # Clear base varaiable and modifiers
  unset(${BASE_VAR})
  foreach(MOD ${MODIFIERS})
    unset(${BASE_VAR}_${MOD})
  endforeach(MOD)

  unset(CURRENT_MOD)
  foreach(ITEM ${ARGN})
    list(FIND MODIFIERS "${ITEM}" MODIFIER_INDEX)

    if(MODIFIER_INDEX GREATER -1)
      # Item is a modifier. Change the current modifier.
      # Note we prefix CURRENT_MOD with '_' so we can just dereference an empty variable to write into BASE_VAR when
      # not modified.
      set(CURRENT_MOD _${ITEM})
    else(MODIFIER_INDEX GREATER -1)
      # Item is not a modifier. Add to the current list.
      list(APPEND ${BASE_VAR}${CURRENT_MOD} "${ITEM}")
    endif(MODIFIER_INDEX GREATER -1)
  endforeach(ITEM)

  # Propagate lists to parent
  set(${BASE_VAR} "${${BASE_VAR}}" PARENT_SCOPE)
  foreach(MOD ${MODIFIERS})
    set(${BASE_VAR}_${MOD} "${${BASE_VAR}_${MOD}}" PARENT_SCOPE)
  endforeach(MOD)
endfunction(ras_split_list)

#-------------------------------------------------------------------------------
# ras_show_properties(<GLOBAL           |
#                  DIRECTORY [dir]  |
#                  TARGET <target>  |
#                  SOURCE <source>  |s
#                  INSTALL <file>   |
#                  TEST <test>      |
#                  CACHE <entry>    |
#                  VARIABLE>)
# Print all the properties and their values of the selected type. For example, target properties are listed by
# specifying ras_show_properties(TARGET <target>), while global properties are listed using ras_show_properties(GLOBAL)
function(ras_show_properties TYPE)
  # Validate TYPE values requiring an argument.
  set(TARGETED_ITEMS CACHE INSTALL SOURCE TARGET TEST)
  list(FIND TARGETED_ITEMS ${TYPE} TYPE_INDEX)
  if(TYPE_INDEX GREATER -1)
    if("${ARGV1}" STREQUAL "")
      message(SEND_ERROR "ras_show_properties(${TYPE} <item>) missing <item> argument.")
      return()
    endif("${ARGV1}" STREQUAL "")
  endif(TYPE_INDEX GREATER -1)

  execute_process(COMMAND ${CMAKE_COMMAND} --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)
  # Convert command output into a CMake list
  string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
  string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")

  # Filter the property list.
  set(ITEM_PROPERTY_LIST)
  foreach(PROPERTY ${CMAKE_PROPERTY_LIST})
    # message("${PROPERTY}...")
    if(PROPERTY MATCHES ".*<CONFIG>.*")
      # Replace <CONFIG> with each of the build configs and add those to the list.
      foreach(CONFIG Debug;MinSizeRel;Release;RelWithDebInfo)
        string(REGEX REPLACE "(.*)<CONFIG>(.*)" "\\1${CONFIG}\\2" CONFIG_PROPERTY "${PROPERTY}")
        list(APPEND ITEM_PROPERTY_LIST "${CONFIG_PROPERTY}")
        # message("...${CONFIG_PROPERTY}")
      endforeach(CONFIG)
    else(PROPERTY MATCHES "$.*<CONFIG>^")
      list(APPEND ITEM_PROPERTY_LIST "${PROPERTY}")
      # message("...${PROPERTY}")
    endif(PROPERTY MATCHES ".*<CONFIG>.*")
  endforeach(PROPERTY)

  if("${ARGV1}" STREQUAL "")
    message("${TYPE} properties:")
  else("${ARGV1}" STREQUAL "")
    set(ITEM "${ARGV1}")
    message("${ITEM} properties:")
  endif("${ARGV1}" STREQUAL "")
  foreach(PROPERTY ${ITEM_PROPERTY_LIST})
    # message ("Checking ${prop}")
    get_property(PROPERTY_VALUE ${TYPE} ${ITEM} PROPERTY ${PROPERTY} SET)
    if(PROPERTY_VALUE)
      get_property(PROPERTY_VALUE ${TYPE} ${ITEM} PROPERTY ${PROPERTY})
      message ("  ${PROPERTY}: ${PROPERTY_VALUE}")
    endif(PROPERTY_VALUE)
  endforeach(PROPERTY)
endfunction(ras_show_properties)

#-------------------------------------------------------------------------------
# ras_show_target_properties(<target>)
# Print all the properties and their values for <target>. Only lists the properties which have been set for <target>.
function(ras_show_target_properties TARGET)
  # https://stackoverflow.com/questions/32183975/how-to-print-all-the-properties-of-a-target-in-cmake
  if(NOT TARGET ${TARGET})
    message("There is no target named '${TARGET}'")
    return()
  endif(NOT TARGET ${TARGET})

  ras_show_properties(TARGET ${TARGET})
endfunction(ras_show_target_properties)

#-------------------------------------------------------------------------------
# ras_show_variables()
# List all variables and their values.
function(ras_show_variables)
  get_cmake_property(VAR_NAMES VARIABLES)
  foreach (VAR ${VAR_NAMES})
      message("${VAR}: ${${VAR}}")
  endforeach(VAR)
endfunction(ras_show_variables)

#-------------------------------------------------------------------------------
# Converts a boolean variable test to a 1 or 0.
# Warning: this overwrites the current value of VAR
#-------------------------------------------------------------------------------
macro(ras_bool_to_int VAR)
  if(${VAR})
    set(${VAR} 1)
  else(${VAR})
    set(${VAR} 0)
  endif(${VAR})
endmacro(ras_bool_to_int VAR)

#-------------------------------------------------------------------------------
# get_relative_path(<VAR> <PATH> <BASE_PATH>)
#
# Get the relative path to PATH from BASE_DIR and store in VAR.
#
# Both PATH and BASE_DIR are converted to absolute paths before evaluation.
#-------------------------------------------------------------------------------
function(get_relative_path VAR PATH RELATIVE_TO)
  # Convert PATH and RELATIVE_TO to absolute paths first.
  get_filename_component(PATH "${PATH}" ABSOLUTE)
  get_filename_component(RELATIVE_TO "${RELATIVE_TO}" ABSOLUTE)

  set(RELATIVE_PATH)
  while(PATH)
    if(PATH STREQUAL RELATIVE_TO)
      break()
    endif(PATH STREQUAL RELATIVE_TO)
    get_filename_component(TAIL "${PATH}" NAME)
    get_filename_component(PATH "${PATH}" DIRECTORY)
    if(RELATIVE_PATH)
      set(RELATIVE_PATH "${TAIL}/${RELATIVE_PATH}")
    else(RELATIVE_PATH)
      set(RELATIVE_PATH "${TAIL}")
    endif(RELATIVE_PATH)
  endwhile(PATH)

  set(${VAR} "${RELATIVE_PATH}" PARENT_SCOPE)
endfunction(get_relative_path)

