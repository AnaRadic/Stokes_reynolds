add_executable ("Stokes_reynolds" bc_extension.hh driver.hh stokes_reynolds.cc bdry_integrate.hh )
target_link_dune_default_libraries ("Stokes_reynolds")
