#pragma once
#ifndef RESERVOIR_H
#define RESERVOIR_H

#include "component.hxx"

/// Multi-region 1D reservoir model
///
/// Provides cross-field transport for (typically) neutrals in 1D by coupling
/// each cell to a reservoir at fixed density and temperature. The domain is split
/// by parallel position (baffle_position) into a main-SOL and a divertor-SOL
/// reservoir, each with its own density, temperature and exchange velocity factor.
/// A divertor-PFR region is stubbed out for future use.
///
/// Bi-channel behaviour: when the reservoir acts as a source (reservoir denser than
/// the local plasma) particles arrive at the reservoir temperature with no parallel
/// momentum; when it acts as a sink they leave at the local temperature carrying
/// local parallel momentum.
struct Reservoir : public Component {

  /// # Inputs
  /// - <name>
  ///   - density_main_sol / density_div_sol / density_div_pfr   BoutReal [m^-3]
  ///   - temperature_main_sol / temperature_div_sol / temperature_div_pfr  BoutReal [eV]
  ///   - velocity_factor_main_sol / velocity_factor_div_sol / velocity_factor_div_pfr
  ///                                                            BoutReal, fraction of v_th
  ///   - xpoint_position       BoutReal [m]
  ///   - baffle_position       BoutReal [m], main-SOL / div-SOL border
  ///   - reservoir_sink_only   bool, only remove particles if true
  ///   - density_floor         BoutReal [m^-3], avoids divide-by-zero
  ///   - diagnose              bool, save sources/sinks
  Reservoir(std::string name, Options& alloptions, Solver*);

  void outputVars(Options& state) override;

private:
  /// # Inputs
  /// - fieldline_geometry_cell_side_area
  /// - fieldline_geometry_cell_volume
  /// - species:<name>: density, pressure, temperature, momentum, AA
  ///
  /// # Outputs
  /// - species:<name>: density_source, energy_source, momentum_source
  void transform_impl(GuardedOptions& state) override;

  std::string name;           ///< Short name of the species e.g. h+

  BoutReal density_div_sol, density_div_pfr, density_main_sol;
  BoutReal temperature_div_sol, temperature_div_pfr, temperature_main_sol;
  BoutReal velocity_factor_div_sol, velocity_factor_div_pfr, velocity_factor_main_sol;
  BoutReal en_src_multiplier, mv_src_multiplier;

  Field3D location_div_sol, location_div_pfr, location_main_sol;

  bool diagnose, reservoir_sink_only;
  BoutReal baffle_position, xpoint_position;  ///< Parallel positions of reservoirs and x-point
  BoutReal density_floor;                     ///< Avoids divide-by-zero in source terms

  Field3D Nrate;
  Field3D density_source_main_sol, energy_source_main_sol, momentum_source_main_sol;
  Field3D density_source_div_sol, energy_source_div_sol, momentum_source_div_sol;
  Field3D density_source_div_pfr, energy_source_div_pfr, momentum_source_div_pfr;

  Field3D lpar;    ///< Parallel connection length, 0 at midplane edge
  Field3D area;    ///< Cross-sectional area in direction of cross-field transport
  Field3D volume;
};

namespace {
RegisterComponent<Reservoir> registercomponentreservoir("reservoir");
}

#endif // RESERVOIR_H