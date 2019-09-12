add_executable ("Stokes_reynolds" bc_extension.hh driver.hh stokes_reynolds.cc )
target_link_dune_default_libraries ("Stokes_reynolds")
