foreach(dim 1 2 3)
  add_executable("solution04-${dim}d" solution04.cc)
  target_compile_definitions("solution04-${dim}d" PUBLIC "GRIDDIM=${dim}")
endforeach()

dune_symlink_to_source_files(FILES tutorial04.ini)
