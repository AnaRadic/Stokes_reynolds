add_executable ("stokes_reynolds" bc_extension.hh driver.hh stokes_reynolds.cc bdry_integrate.hh)
target_link_dune_default_libraries ("stokes_reynolds")
dune_symlink_to_source_files(FILES "stokes.input")
dune_symlink_to_source_files(FILES "mreza.msh")
