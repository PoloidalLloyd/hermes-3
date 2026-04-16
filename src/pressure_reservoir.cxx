#include "../include/pressure_reservoir.hxx"
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/msg_stack.hxx>

using bout::globals::mesh;

Reservoir::Reservoir(std::string name, Options& alloptions, Solver*)
    : Component({readOnly("fieldline_geometry_cell_side_area"),
                 readOnly("fieldline_geometry_cell_volume"),
                 readOnly("species:{name}:density", Regions::Interior),
                 readOnly("species:{name}:pressure", Regions::Interior),
                 readOnly("species:{name}:temperature", Regions::Interior),
                 readOnly("species:{name}:momentum", Regions::Interior),
                 readOnly("species:{name}:AA"),
                 readWrite("species:{name}:density_source"),
                 readWrite("species:{name}:energy_source"),
                 readWrite("species:{name}:momentum_source")}),
        name(name) {
    TRACE("PressureReservoir::PressureReservoir");

    const Options& units = alloptions["units"];
    const BoutReal Nnorm = units["inv_meters_cubed"];
    const BoutReal Lnorm = units["meters"];

    // Get the options for this species
    Options& options = alloptions[name];

    density_div_sol = options["density_div_sol"]
                            .doc("Set the density of the reservoir in [m^-3]. Default 1e19")
                            .withDefault<BoutReal>(1e19)
                        / Nnorm;

    density_div_pfr = options["density_div_pfr"]
                            .doc("Set the density of the reservoir in [m^-3]. Default 1e19")
                            .withDefault<BoutReal>(1e19)
                        / Nnorm;
    density_main_sol = options["density_main_sol"]
                            .doc("Set the density of the main SOL reservoir in [m^-3]. Default 1e19")
                            .withDefault<BoutReal>(1e19)
                        / Nnorm;
    temperature_div_sol = options["temperature_div_sol"]
                            .doc("Set the temperature of the divertor SOL reservoir. Default 1 eV")
                            .withDefault<BoutReal>(1)
                        / Tnorm;
    temperature_div_pfr = options["temperature_div_pfr"]
                            .doc("Set the temperature of the divertor PFR reservoir. Default 1 eV")
                            .withDefault<BoutReal>(1)
                        / Tnorm;
    temperature_main_sol = options["temperature_main_sol"]
                            .doc("Set the temperature of the main SOL reservoir. Default 1 eV")
                            .withDefault<BoutReal>(1)
                        / Tnorm;



  velocity_factor_div_sol =
      options["velocity_factor_div_sol"]
          .doc("Exchange speed of divertor sol reservoir as fraction of thermal velocity")
          .withDefault<BoutReal>(1);
  velocity_factor_div_pfr =
      options["velocity_factor_div_pfr"]
          .doc("Exchange speed of divertor pfr reservoir as fraction of thermal velocity")
          .withDefault<BoutReal>(1);
  velocity_factor_main_sol =
      options["velocity_factor_main_sol"]
          .doc("Exchange speed of main SOL reservoir as fraction of thermal velocity")
          .withDefault<BoutReal>(1);
