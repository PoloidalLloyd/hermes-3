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

