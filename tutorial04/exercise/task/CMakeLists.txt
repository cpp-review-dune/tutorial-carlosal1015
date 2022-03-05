foreach(dim 1 2 3)
  add_executable("exercise04-${dim}d" exercise04.cc)
  target_compile_definitions("exercise04-${dim}d" PUBLIC "GRIDDIM=${dim}")
endforeach()

dune_symlink_to_source_files(FILES exercise04.ini)
