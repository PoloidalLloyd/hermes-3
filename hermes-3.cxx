/*

    Copyright B.Dudson, J.Leddy, University of York, 2016-2020
              email: benjamin.dudson@york.ac.uk

    This file is part of Hermes-3 (Hot ion, multifluid)

    Hermes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Hermes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hermes.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "hermes-3.hxx"
#include "revision.hxx"

#include <bout/constants.hxx>
#include "include/ionisation.hxx"
#include "include/neutral_mixed.hxx"
#include "include/diamagnetic_drift.hxx"
#include "include/electron_force_balance.hxx"
#include "include/evolve_density.hxx"
#include "include/isothermal.hxx"
#include "include/sheath_closure.hxx"
#include "include/sheath_boundary.hxx"
#include "include/sheath_boundary_simple.hxx"
#include "include/vorticity.hxx"
#include "include/fixed_fraction_ions.hxx"
#include "include/evolve_pressure.hxx"
#include "include/evolve_momentum.hxx"
#include "include/quasineutral.hxx"
#include "include/sound_speed.hxx"
#include "include/zero_current.hxx"
#include "include/anomalous_diffusion.hxx"
#include "include/recycling.hxx"
#include "include/collisions.hxx"
#include "include/amjuel_hyd_ionisation.hxx"
#include "include/amjuel_hyd_recombination.hxx"
#include "include/amjuel_helium.hxx"
#include "include/adas_neon.hxx"
#include "include/noflow_boundary.hxx"
#include "include/neutral_parallel_diffusion.hxx"
#include "include/solkit_neutral_parallel_diffusion.hxx"
#include "include/hydrogen_charge_exchange.hxx"
#include "include/solkit_hydrogen_charge_exchange.hxx"
#include "include/upstream_density_feedback.hxx"
#include "include/thermal_force.hxx"
#include "include/ion_viscosity.hxx"
#include "include/relax_potential.hxx"
#include "include/fixed_density.hxx"
#include "include/fixed_velocity.hxx"
#include "include/transform.hxx"
#include "include/fixed_fraction_radiation.hxx"

#include "include/loadmetric.hxx"

int Hermes::init(bool restarting) {

  auto &options = Options::root()["hermes"];
  
  output.write("\nGit Version of Hermes: {:s}\n", hermes::version::revision);
  options["revision"] = hermes::version::revision;
  options["revision"].setConditionallyUsed();

  // Save the Hermes version in the output dump files
  dump.setAttribute("", "HERMES_REVISION", hermes::version::revision);

  // Choose normalisations
  Tnorm = options["Tnorm"].doc("Reference temperature [eV]").withDefault(100.);
  Nnorm = options["Nnorm"].doc("Reference density [m^-3]").withDefault(1e19);
  Bnorm = options["Bnorm"].doc("Reference magnetic field [T]").withDefault(1.0);

  Cs0 = sqrt(SI::qe * Tnorm / SI::Mp); // Reference sound speed [m/s]
  Omega_ci = SI::qe * Bnorm / SI::Mp;  // Ion cyclotron frequency [1/s]
  rho_s0 = Cs0 / Omega_ci;             // Length scale [m]

  SAVE_ONCE(Tnorm, Nnorm, Bnorm); // Save normalisations
  SAVE_ONCE(Cs0, Omega_ci, rho_s0);

  // Put normalisation quantities into an Options to use later
  units["inv_meters_cubed"] = Nnorm;
  units["eV"] = Tnorm;
  units["Tesla"] = Bnorm;
  units["seconds"] = 1./Omega_ci;
  units["meters"] = rho_s0;

  // Put into the options tree, so quantities can be normalised
  // when creating components
  Options::root()["units"] = units;
  Options::root()["units"].setConditionallyUsed();

  /////////////////////////////////////////////////////////
  // Load metric tensor from the mesh, passing length and B
  // field normalisations
  TRACE("Loading metric tensor");

  if (options["loadmetric"]
          .doc("Load Rxy, Bpxy etc. to create orthogonal metric?")
          .withDefault(true)) {
    LoadMetric(rho_s0, Bnorm);
  } else if (options["normalise_metric"]
                 .doc("Normalise input metric tensor? (assumes input is in SI units)")
                 .withDefault<bool>(true)) {
    Coordinates *coord = mesh->getCoordinates();
    // To use non-orthogonal metric
    // Normalise
    coord->dx /= rho_s0 * rho_s0 * Bnorm;
    coord->Bxy /= Bnorm;
    // Metric is in grid file - just need to normalise
    coord->g11 /= SQ(Bnorm * rho_s0);
    coord->g22 *= SQ(rho_s0);
    coord->g33 *= SQ(rho_s0);
    coord->g12 /= Bnorm;
    coord->g13 /= Bnorm;
    coord->g23 *= SQ(rho_s0);

    coord->J *= Bnorm / rho_s0;

    coord->g_11 *= SQ(Bnorm * rho_s0);
    coord->g_22 /= SQ(rho_s0);
    coord->g_33 /= SQ(rho_s0);
    coord->g_12 *= Bnorm;
    coord->g_13 *= Bnorm;
    coord->g_23 /= SQ(rho_s0);

    coord->geometry(); // Calculate other metrics
  }

  // Tell the components if they are restarting
  options["restarting"] = restarting;
  options["restarting"].setConditionallyUsed();

  // Put pointer to restart file in global variable
  // Note: It would probably be better to pass to components as argument,
  // but awaiting DataFile refactoring
  set_restart_datafile(&restart);

  TRACE("Creating components");
  
  // Create the components
  // Here options is passed as the scheduler configuration, so that
  // settings in [hermes] are used.
  // Options::root() is passed as the root of the component options, so that
  // individual components use their own sections, rather than subsections of [hermes].
  scheduler = ComponentScheduler::create(options, Options::root(), solver);
  
  return 0;
}

int Hermes::rhs(BoutReal time) {
  // Need to reset the state, since fields may be modified in transform steps
  state = Options();
  
  set(state["time"], time);
  state["units"] = units; 

  // Call all the components
  scheduler->transform(state);

  return 0;
}

/*!
 * Preconditioner. Solves the heat conduction
 *
 * @param[in] t  The simulation time
 * @param[in] gamma   Factor in front of the Jacobian in (I - gamma*J). Related
 * to timestep
 * @param[in] delta   Not used here
 */
int Hermes::precon(BoutReal t, BoutReal gamma, BoutReal UNUSED(delta)) {
  state["time"] = t;
  scheduler->precon(state, gamma);
  return 0;
}

// Standard main() function
BOUTMAIN(Hermes);
