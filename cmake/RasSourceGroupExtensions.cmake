#-------------------------------------------------------------------------------
# A helper for setting up source_group()
#
# ras_make_source_groups([RELATIVE <shared-root-path>] [GENERATED <generated-source-files>] SOURCE <source>)
#-------------------------------------------------------------------------------
function(ras_make_source_groups)
  # Split into:
  # - MSG_GENERATED: generated source files
  # - MSG_UNPARSED_ARGUMENTS: source files (possibly including generated items)
  cmake_parse_arguments(MSG "" "RELATIVE" "GENERATED;SOURCE" ${ARGN})

  # message("MSG_SOURCE: ${MSG_SOURCE}")
  # message("MSG_GENERATED: ${MSG_GENERATED}")

  # Build the source list prefixed with CMAKE_CURRENT_LIST_DIR
  set(MSG_SOURCE_ABSOLUTE)
  foreach(SOURCE ${MSG_SOURCE})
    # Filter out items in the MSG_GENERATED list
    list(FIND MSG_GENERATED "${SOURCE}" GEN_IDX)
    if(GEN_IDX LESS 0)
      # Convert to absolute path for the source_group(TREE) directive
      get_filename_component(SOURCE2 "${SOURCE}" ABSOLUTE)
      # message("${SOURCE} -> ${SOURCE2}")
      # Add to list.
      list(APPEND MSG_SOURCE_ABSOLUTE ${SOURCE2})
    endif(GEN_IDX LESS 0)
  endforeach(SOURCE)

  # If MSG_RELATIVE not specified, use CMAKE_CURRENT_LIST_DIR
  if(NOT MSG_RELATIVE)
    set(MSG_RELATIVE "${CMAKE_CURRENT_LIST_DIR}")
  endif(NOT MSG_RELATIVE)

  # message("MSG_RELATIVE: ${MSG_RELATIVE}")

  # Generated file filter
  if(MSG_GENERATED)
    source_group("generated" FILES ${MSG_GENERATED})
  endif(MSG_GENERATED)
  source_group(TREE "${MSG_RELATIVE}" PREFIX source FILES ${MSG_SOURCE_ABSOLUTE})
endfunction(ras_make_source_groups)
