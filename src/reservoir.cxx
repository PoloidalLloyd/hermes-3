
#include "../include/reservoir.hxx"
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>

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
  AUTO_TRACE();

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Lnorm = units["meters"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  // Get the options for this species
  Options& options = alloptions[name];

  density_div_sol = options["density_div_sol"]
                          .doc("Set the density of the reservoir in [m^-3]. Default 1e19")
                          .withDefault<BoutReal>(1e19)
                      / Nnorm;
  en_src_multiplier = options["en_src_multiplier"]
                          .doc("Multiply the energy source by this factor. Default 1")
                          .withDefault<BoutReal>(1);
  mv_src_multiplier = options["mv_src_multiplier"]
                          .doc("Multiply the momentum source by this factor. Default 1")
                          .withDefault<BoutReal>(1);


  density_div_pfr = options["density_div_pfr"]
                          .doc("Set the density of the reservoir in [m^-3]. Default 1e19")
                          .withDefault<BoutReal>(1e19)
                      / Nnorm;
  density_main_sol = options["density_main_sol"]
                          .doc("Set the density of the reservoir in [m^-3]. Default 1e19")
                          .withDefault<BoutReal>(1e19)
                      / Nnorm;

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

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);

  reservoir_sink_only = 
    options["reservoir_sink_only"].doc("Set reservoir to only take particles away?").withDefault<bool>(true);

  density_floor = options["density_floor"].doc("Minimum density floor [m^-3]").withDefault<BoutReal>(1e-7) / Nnorm;

  xpoint_position =
      options["xpoint_position"]
          .doc("Parallel position of X-point [m]")
          .withDefault<BoutReal>(0)
      / Lnorm;  // Normalize to match lpar units

  baffle_position =
      options["baffle_position"]
          .doc("Parallel position of border between upstream-SOL and divertor-SOL reservoirs [m]")
          .withDefault<BoutReal>(xpoint_position * Lnorm)  // Default is unnormalized xpoint
      / Lnorm;  // Normalize to match lpar units

  // Get the area of the radial reservoir_region
  area = options["area"]
      .doc("Area of radial regions. Specify as an analytical function.")
      .withDefault<Field3D>(0.0);

  area = 1; // This is constant, and wrong! 
 
  BoutReal Anorm = SQ(get<BoutReal>(units["meters"]));
  area /= Anorm;

  // Get lpar
  const int MYPE = BoutComm::rank();   // Current rank
  const int NPES = BoutComm::size();    // Number of procs
  const int NYPE = NPES / mesh->NXPE;    // Number of procs in Y
  Coordinates *coord = mesh->getCoordinates();

  lpar = 0;
  BoutReal offset = 0;   // Offset to ensure ylow domain boundary starts at 0
  auto dy = coord->dy;
  auto J = coord->J;

  volume = J * dy;  // Volume of each cell


  lpar(0,0,0) = 0.5 * dy(0,0,0);
  for (int id = 0; id <= NYPE-1; ++id) {   // Iterate through each proc
    for (int j = 0; j <= mesh->LocalNy; ++j) {
          if (j > 0) {
            lpar(0, j, 0) = lpar(0, j-1, 0) + 0.5 * dy(0, j-1, 0) + 0.5 * dy(0, j, 0);
          }
    }
    mesh->communicate(lpar);  // Communicate across guard cells so other procs keep counting
    
    if (MYPE == 0) {
      offset = (lpar(0,1,0) + lpar(0,2,0)) / 2;  // Offset lives on first proc
    }
  }

  MPI_Bcast(&offset, 1, MPI_DOUBLE, 0, BoutComm::get());  // Ensure all procs get offset
  lpar -= offset;
  lpar /= Lnorm;  // Normalize lpar to match baffle_position units

  // Capture reservoir locations and regions
  // location_div_sol = 0;

  // Region<Ind3D>::RegionIndices indices;
  // BOUT_FOR_SERIAL(i, location_div_sol.getRegion("RGN_NOBNDRY")) {
  //   if (lpar[i] > baffle_position) {
  //     // Add this cell to the iteration
  //     indices.push_back(i);
  //     location_div_sol[i] = 1;
  //   }
  // }
  // region_div_sol = Region<Ind3D>(indices);

  substitutePermissions("name", {name});
}

void Reservoir::transform_impl(GuardedOptions& state) {
  AUTO_TRACE();

  // Initialize output fields to zero (needed for outputVars even if we return early)
 
  density_source_main_sol = 0;
  density_source_div_sol = 0;
  density_source_div_pfr = 0;
  energy_source_main_sol = 0;
  energy_source_div_sol = 0;
  energy_source_div_pfr = 0;
  momentum_source_main_sol = 0;
  momentum_source_div_sol = 0;
  momentum_source_div_pfr = 0;
  location_main_sol = 0;
  location_div_sol = 0;
  location_div_pfr = 0;

  // We are operating on only one species
  GuardedOptions species = state["species"][name];

  // If fieldline geometry is available in the state, use it for area/volume (should add error handling)
  area = get<Field3D>(state["fieldline_geometry_cell_side_area"]);
  volume = get<Field3D>(state["fieldline_geometry_cell_volume"]);

  // Get conditions. Boundary conditions do not need to be set
  // because we don't use the state in the boundary cells.
  Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
  Field3D P = GET_NOBOUNDARY(Field3D, species["pressure"]);
  Field3D T = GET_NOBOUNDARY(Field3D, species["temperature"]);
  Field3D NV = GET_NOBOUNDARY(Field3D, species["momentum"]);
  const BoutReal AA = get<BoutReal>(species["AA"]);

  Field3D vth = 0.25 * sqrt( (8*T) / (PI*AA) );    // 1D maxwellian particle flux in 3D system (Stangeby)

  // Pressure and momentum mass flows proportional to density
  // When flow is reversed and these become sources, the new particles have
  // the same pressure and momentum as the local particles.

  // Apply flooring to density to prevent division by zero
  // No flooring needed for P and NV - they can be zero or negative
  Field3D Nfloor = softFloor(N, density_floor);

  BOUT_FOR(i, N.getRegion("RGN_NOBNDRY")) {
    // Main SOL reservoir
    //////////////////////////////////////////
    if (lpar[i] <= baffle_position) {
      BoutReal Nrate = (density_main_sol - Nfloor[i]) * area[i] * vth[i] * velocity_factor_main_sol;
      if (reservoir_sink_only && Nrate > 0) {
        Nrate = 0;
      };

      BoutReal Prate  = P[i]  / Nfloor[i] * Nrate;
      BoutReal NVrate = NV[i] / Nfloor[i] * Nrate;

      density_source_main_sol[i]  += Nrate / volume[i];
      energy_source_main_sol[i]   += (3. / 2) * Prate / volume[i] * en_src_multiplier;
      momentum_source_main_sol[i] += NVrate / volume[i] * mv_src_multiplier;
      location_main_sol[i] = 1;
    }

    // Divertor SOL reservoir
    //////////////////////////////////////////
    if (lpar[i] > baffle_position) {
      BoutReal Nrate = (density_div_sol - Nfloor[i]) * area[i] * vth[i] * velocity_factor_div_sol;
      if (reservoir_sink_only && Nrate > 0) {
        Nrate = 0;
      };

      BoutReal Prate  = P[i]  / Nfloor[i] * Nrate;
      BoutReal NVrate = NV[i] / Nfloor[i] * Nrate;

      density_source_div_sol[i]  += Nrate / volume[i];
      energy_source_div_sol[i]   += ((3. / 2) * Prate / volume[i]) * en_src_multiplier;
      momentum_source_div_sol[i] += (NVrate / volume[i]) * mv_src_multiplier;
      location_div_sol[i] = 1;
    }

    // Divertor PFR reservoir
    //////////////////////////////////////////
    // Note: PFR region requires additional logic to define its spatial extent
    // Currently using velocity_factor_div_pfr = 0 to disable

  }

  // Add all source terms to the species
  add(species["density_source"], density_source_main_sol + density_source_div_sol + density_source_div_pfr);
  add(species["energy_source"], energy_source_main_sol + energy_source_div_sol + energy_source_div_pfr);
  add(species["momentum_source"], momentum_source_main_sol + momentum_source_div_sol + momentum_source_div_pfr);
}

void Reservoir::outputVars(Options& state) {
  AUTO_TRACE();

  // Normalisations - get them from state with the correct keys
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  auto Lnorm = get<BoutReal>(state["rho_s0"]);  // Use rho_s0 for length
  BoutReal Anorm = SQ(Lnorm);  // Area normalization
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  // Always save reservoir location (time independent)
  
  if (diagnose) {

    // Area
    ///////////////////////////

    set_with_attrs(state[std::string("area_") + name], area,
          {{"units", "m^2"},
          {"conversion", Anorm},
          {"standard_name", "area"},
          {"long_name", name + std::string(" area")},
          {"source", "reservoir"}});



    // Reservoir rate
    set_with_attrs(
      state[{std::string("Nrate") + name + std::string("_rsv")}], Nrate,
      {{"time_dimension", "t"},
      {"units", "m^-3 s^-1"},
      {"conversion", Omega_ci},
      {"standard_name", "particle transfer rate"},
      {"long_name", name + std::string(" particle transfer rate from reservoir")},
      {"source", "reservoir"}});


    // Main sol reservoir
    ///////////////////////////
    set_with_attrs(state[{std::string("rsv_main_sol_") + name}], location_main_sol,
                 {{"standard_name", name + std::string("reservoir location")},
                  {"long_name", name + std::string("main sol reservoir location")},
                  {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("S") + name + std::string("_rsv_main_sol")}], density_source_main_sol,
          {{"time_dimension", "t"},
          {"units", "m^-3 s^-1"},
          {"conversion", Nnorm * Omega_ci},
          {"standard_name", "particle transfer"},
          {"long_name", name + std::string(" density transfer from div sol reservoir")},
          {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("E") + name + std::string("_rsv_main_sol")}], energy_source_main_sol,
          {{"time_dimension", "t"},
          {"units", "W m^-3"},
          {"conversion", Pnorm * Omega_ci},
          {"standard_name", "energy transfer"},
          {"long_name", name + std::string(" energy transfer from div sol reservoir")},
          {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("F") + name + std::string("_rsv_main_sol")}], momentum_source_main_sol,
          {{"time_dimension", "t"},
          {"units", "kg m^-2 s^-2"},
          {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
          {"standard_name", "momentum transfer"},
          {"long_name", name + std::string(" momentum transfer from div sol reservoir")},
          {"source", "reservoir"}});



    // Div sol reservoir
    ///////////////////////////
    set_with_attrs(state[{std::string("rsv_div_sol_") + name}], location_div_sol,
                 {{"standard_name", name + std::string("reservoir location")},
                  {"long_name", name + std::string("div sol reservoir location")},
                  {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("S") + name + std::string("_rsv_div_sol")}], density_source_div_sol,
          {{"time_dimension", "t"},
          {"units", "m^-3 s^-1"},
          {"conversion", Nnorm * Omega_ci},
          {"standard_name", "particle transfer"},
          {"long_name", name + std::string(" density transfer from div sol reservoir")},
          {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("E") + name + std::string("_rsv_div_sol")}], energy_source_div_sol,
          {{"time_dimension", "t"},
          {"units", "W m^-3"},
          {"conversion", Pnorm * Omega_ci},
          {"standard_name", "energy transfer"},
          {"long_name", name + std::string(" energy transfer from div sol reservoir")},
          {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("F") + name + std::string("_rsv_div_sol")}], momentum_source_div_sol,
          {{"time_dimension", "t"},
          {"units", "kg m^-2 s^-2"},
          {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
          {"standard_name", "momentum transfer"},
          {"long_name", name + std::string(" momentum transfer from div sol reservoir")},
          {"source", "reservoir"}});

      


    // Div pfr reservoir
    ///////////////////////////
    set_with_attrs(state[{std::string("rsv_div_pfr_") + name}], location_div_pfr,
                 {{"standard_name", name + std::string("reservoir location")},
                  {"long_name", name + std::string("div pfr reservoir location")},
                  {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("S") + name + std::string("_rsv_div_pfr")}], density_source_div_pfr,
          {{"time_dimension", "t"},
          {"units", "m^-3 s^-1"},
          {"conversion", Nnorm * Omega_ci},
          {"standard_name", "particle transfer"},
          {"long_name", name + std::string(" density transfer from div sol reservoir")},
          {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("E") + name + std::string("_rsv_div_pfr")}], energy_source_div_pfr,
          {{"time_dimension", "t"},
          {"units", "W m^-3"},
          {"conversion", Pnorm * Omega_ci},
          {"standard_name", "energy transfer"},
          {"long_name", name + std::string(" energy transfer from div sol reservoir")},
          {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("F") + name + std::string("_rsv_div_pfr")}], momentum_source_div_pfr,
          {{"time_dimension", "t"},
          {"units", "kg m^-2 s^-2"},
          {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
          {"standard_name", "momentum transfer"},
          {"long_name", name + std::string(" momentum transfer from div sol reservoir")},
          {"source", "reservoir"}});

    // Other stuff
    set_with_attrs(    // Doesn't need unnormalising for some reason
          state[std::string("lpar")], lpar,
          {{"units", "m"},
          {"standard_name", "Parallel distance from midplane"},
          {"long_name", "Parallel distance from midplane"},
          {"source", "reservoir"}});


  }
}
