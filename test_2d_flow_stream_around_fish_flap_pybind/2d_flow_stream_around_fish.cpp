/**
 * @file 2d_flow_stream_around_fish.cpp
 * @brief fish swimming driven by active muscles
 * @author Mai Ye, Yaru Ren and Xiangyu Hu
 */
#include <pybind11/pybind11.h>
#include "2d_flow_stream_around_fish.h"
#include "sphinxsys.h"

using namespace SPH;
namespace py = pybind11;

class SphBasicSystemSetting : public SphBasicGeometrySetting
{
  protected:
    BoundingBox system_domain_bounds;
    SPHSystem sph_system;
    IOEnvironment io_environment;

  public:
    SphBasicSystemSetting(int episode_env) :
        system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW)),
        sph_system(system_domain_bounds, particle_spacing_ref),
        io_environment(sph_system, episode_env) {}
};

class SphFishRelaxationEnvironment : public SphBasicSystemSetting
{
  protected:
    FluidBody water_block;
    SolidBody fish_body, flap_body;
    ObserverBody fish_observer;

  public:
    SphFishRelaxationEnvironment(int episode_env) :
        SphBasicSystemSetting(episode_env),
        water_block(sph_system, makeShared<WaterBlock>("WaterBody")),
        fish_body(sph_system, makeShared<FishBody>("FishBody")),
        flap_body(sph_system, makeShared<FlapBody>("FlapBody")),
        fish_observer(sph_system, "FishObserver")
    {
        water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
        water_block.generateParticles<ParticleGeneratorLattice>();

        fish_body.defineAdaptationRatios(1.15, 2.0);
        fish_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
        fish_body.defineParticlesAndMaterial<ElasticSolidParticles, FishBodyComposite>();
        fish_body.generateParticles<ParticleGeneratorLattice>();

        flap_body.defineAdaptationRatios(1.15, 2.0);
        flap_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
        flap_body.defineParticlesAndMaterial<ElasticSolidParticles, FishBodyComposite>();
        flap_body.generateParticles<ParticleGeneratorLattice>();

        fish_observer.generateParticles<ObserverParticleGenerator>(observation_locations);
    }
};

class SphFishRelaxation : public SphFishRelaxationEnvironment
{
  protected:
    InnerRelation fish_inner, flap_inner, water_block_inner;
    ContactRelation water_block_contact, fish_contact, flap_contact, fish_observer_contact_water, fish_observer_contact_fish;
    ComplexRelation water_block_complex;
    SimpleDynamics<RandomizeParticlePosition> random_fish_body_particles, random_flap_body_particles;
    BodyStatesRecordingToVtp write_fish_body, write_flap_body;
    ReloadParticleIO write_fish_particle_reload_files, write_flap_particle_reload_files;
    relax_dynamics::RelaxationStepInner relaxation_step_inner_fish, relaxation_step_inner_flap;
    int ite_p = 0;

  public:
    explicit SphFishRelaxation(int episode_env) :
        SphFishRelaxationEnvironment(episode_env),
        fish_inner(fish_body),
        flap_inner(flap_body),
        water_block_inner(water_block),
        water_block_contact(water_block, {&fish_body, &flap_body}),
        fish_contact(fish_body, {&water_block}),
        flap_contact(flap_body, {&water_block}),
        fish_observer_contact_water(fish_observer, {&water_block}),
        fish_observer_contact_fish(fish_observer, {&fish_body}),
        water_block_complex(water_block_inner, water_block_contact),
        random_fish_body_particles(fish_body),
        random_flap_body_particles(flap_body),
        write_fish_body(io_environment, fish_body),
        write_flap_body(io_environment, flap_body),
        write_fish_particle_reload_files(io_environment, {&fish_body}),
        write_flap_particle_reload_files(io_environment, {&flap_body}),
        relaxation_step_inner_fish(fish_inner),
        relaxation_step_inner_flap(flap_inner)
    {
        random_fish_body_particles.exec(0.25);
        relaxation_step_inner_fish.SurfaceBounding().exec();
        write_fish_body.writeToFile();
        while (ite_p < 1000)
        {
            relaxation_step_inner_fish.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the fish body N = " << ite_p << "\n";
                write_fish_body.writeToFile(ite_p);
            }
        }
        write_fish_particle_reload_files.writeToFile();

        ite_p = 0;
        random_flap_body_particles.exec(0.25);
        relaxation_step_inner_flap.SurfaceBounding().exec();
        write_flap_body.writeToFile();
        while (ite_p < 1000)
        {
            relaxation_step_inner_flap.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the flap body N = " << ite_p << "\n";
                write_flap_body.writeToFile(ite_p);
            }
        }
        write_flap_particle_reload_files.writeToFile();

        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
    }
};

class SphFishReloadTrainEnvironment : public SphBasicSystemSetting
{
  protected:
    FluidBody water_block;
    SolidBody fish_body, flap_body;
    ObserverBody fish_observer;

  public:
    SphFishReloadTrainEnvironment(int episode_env) :
        SphBasicSystemSetting(episode_env),
        water_block(sph_system, makeShared<WaterBlock>("WaterBody")),
        fish_body(sph_system, makeShared<FishBody>("FishBody")),
        flap_body(sph_system, makeShared<FlapBody>("FlapBody")),
        fish_observer(sph_system, "FishObserver")
    {
        water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
        water_block.generateParticles<ParticleGeneratorLattice>();

        fish_body.defineAdaptationRatios(1.15, 2.0);
        fish_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
        fish_body.defineParticlesAndMaterial<ElasticSolidParticles, FishBodyComposite>();
        fish_body.generateParticles<ParticleGeneratorReload>(io_environment, fish_body.getName());

        flap_body.defineAdaptationRatios(1.15, 2.0);
        flap_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
        flap_body.defineParticlesAndMaterial<ElasticSolidParticles, FishBodyComposite>();
        flap_body.generateParticles<ParticleGeneratorReload>(io_environment, flap_body.getName());

        fish_observer.generateParticles<ObserverParticleGenerator>(observation_locations);
    }
};

class SphFishReloadTrain : public SphFishReloadTrainEnvironment
{
  protected:
    InnerRelation fish_inner, flap_inner, water_block_inner;
    ContactRelation water_block_contact, fish_contact, flap_contact, fish_observer_contact_water, fish_observer_contact_fish;
    ComplexRelation water_block_complex;
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step;
    BodyAlignedBoxByParticle emitter;
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection;
    BodyAlignedBoxByCell emitter_buffer, disposer;
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>> emitter_buffer_inflow_condition;
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion;
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator;
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density;
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size;
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size;
    SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint;
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation;
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation;
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration;
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction;
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity;
    //----------------------------------------------------------------------
    //	Algorithms of flap moving.
    //----------------------------------------------------------------------
    SimpleDynamics<FlapInitialVelocity> flap_velocity;
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> fish_body_normal_direction, flap_body_normal_direction;
    InteractionWithUpdate<KernelCorrectionMatrixInner> fish_body_corrected_configuration, flap_body_corrected_configuration;
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_fish, viscous_force_on_flap;
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluidRiemann> fluid_force_on_fish_update, fluid_force_on_flap_update;
    solid_dynamics::AverageVelocityAndAcceleration average_fish_velocity_and_acceleration;

    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    SimpleDynamics<FishMaterialInitialization> composite_material_id_fish;
    SimpleDynamics<FlapMaterialInitialization> composite_material_id_flap;
    SimpleDynamics<ImposingActiveStrain> imposing_active_strain;
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> fish_body_computing_time_step_size;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> fish_body_stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> fish_body_stress_relaxation_second_half;
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> fish_body_update_normal, flap_body_update_normal;
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states;
    RestartIO restart_io;
    ObservedQuantityRecording<Real> fish_pressure_probe;
    ObservedQuantityRecording<Vecd> fish_velocity_probe;
    ObservedQuantityRecording<Vecd> fish_position_probe;
    //ObservedQuantityRecording<Real> elastic_energy_probe;
    //ObservedQuantityRecording<Real> elastic_energy_probe;
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int write_episode_vtp = 10;
    Real D_Time = 0.01;  /**< time stamps for output. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;

  public: 
    explicit SphFishReloadTrain(int episode_env) :
        SphFishReloadTrainEnvironment(episode_env),
        fish_inner(fish_body), flap_inner(flap_body), water_block_inner(water_block),
        water_block_contact(water_block, {&fish_body, &flap_body}),
        fish_contact(fish_body, {&water_block}),
        flap_contact(flap_body, {&water_block}),
        fish_observer_contact_water(fish_observer, {&water_block}),
        fish_observer_contact_fish(fish_observer, {&fish_body}),
        water_block_complex(water_block_inner, water_block_contact),
        initialize_a_fluid_step(water_block, makeShared<TimeDependentAcceleration>(Vec2d::Zero())),
        emitter(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize)),
        emitter_inflow_injection(emitter, 10, 0),
        emitter_buffer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize)),
        emitter_buffer_inflow_condition(emitter_buffer),
        disposer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize)),
        disposer_outflow_deletion(disposer, 0),
        free_stream_surface_indicator(water_block_inner, water_block_contact),
        update_fluid_density(water_block_inner, water_block_contact),
        get_fluid_advection_time_step_size(water_block, U_f),
        get_fluid_time_step_size(water_block),
        velocity_boundary_condition_constraint(water_block),
        pressure_relaxation(water_block_inner, water_block_contact),
        density_relaxation(water_block_inner, water_block_contact),
        viscous_acceleration(water_block_inner, water_block_contact),
        transport_velocity_correction(water_block_inner, water_block_contact),
        compute_vorticity(water_block_inner),
        flap_velocity(flap_body),
        fish_body_normal_direction(fish_body),
        flap_body_normal_direction(flap_body),
        fish_body_corrected_configuration(fish_inner),
        flap_body_corrected_configuration(flap_inner),
        viscous_force_on_fish(fish_contact),
        viscous_force_on_flap(flap_contact),
        fluid_force_on_fish_update(fish_contact, viscous_force_on_fish),
        fluid_force_on_flap_update(flap_contact, viscous_force_on_flap),
        average_fish_velocity_and_acceleration(fish_body),
        composite_material_id_fish(fish_body),
        composite_material_id_flap(flap_body),
        imposing_active_strain(fish_body),
        fish_body_computing_time_step_size(fish_body),
        //write_fish_body_elastic_energy(io_environment, fish_body),
        //write_fish_body_kinetic_energy(io_environment, fish_body),
        fish_body_stress_relaxation_first_half(fish_inner),
        fish_body_stress_relaxation_second_half(fish_inner),
        fish_body_update_normal(fish_body),
        flap_body_update_normal(flap_body),
        write_real_body_states(io_environment, sph_system.real_bodies_),
        restart_io(io_environment, sph_system.real_bodies_),
        fish_pressure_probe("Pressure", io_environment, fish_observer_contact_water),
        fish_velocity_probe("Velocity", io_environment, fish_observer_contact_fish),
        fish_position_probe("Position", io_environment, fish_observer_contact_fish)
    {
        GlobalStaticVariables::physical_time_ = 0.0;
        water_block.addBodyStateForRecording<Real>("Pressure");
        water_block.addBodyStateForRecording<int>("Indicator");
        fish_body.addBodyStateForRecording<Real>("Density");
        fish_body.addBodyStateForRecording<int>("MaterialID");
        fish_body.addBodyStateForRecording<Matd>("ActiveStrain");

        /** initialize cell linked lists for all bodies. */
        sph_system.initializeSystemCellLinkedLists();
        /** initialize configurations for all bodies. */
        sph_system.initializeSystemConfigurations();
        /** computing surface normal direction for the fish. */
        fish_body_normal_direction.exec();
        /** computing surface normal direction for the flap. */
        flap_body_normal_direction.exec();
        /** computing linear reproducing configuration for the fish. */
        fish_body_corrected_configuration.exec();
        /** computing linear reproducing configuration for the flap. */
        flap_body_corrected_configuration.exec();
        /** initialize material ids for the fish. */
        composite_material_id_fish.exec();
        /** initialize material ids for the flap. */
        composite_material_id_flap.exec(); 

        write_real_body_states.writeToFile();
        fish_pressure_probe.writeToFile();
        fish_velocity_probe.writeToFile();
        fish_position_probe.writeToFile();
    }

    virtual ~SphFishReloadTrain(){};

    //----------------------------------------------------------------------
    //	Action from python environment.
    //----------------------------------------------------------------------
    void setLambdafromPython(Real lambda_from_python)
    {
        lambda_ = lambda_from_python;
    }

    void setFreqfromPython(Real freq_from_python)
    {
        frequency_ = freq_from_python;
    }
    //----------------------------------------------------------------------
    //	observation for python environment.
    //----------------------------------------------------------------------
    Real getFishPressure(int number)
    {
        return fish_pressure_probe.getCurrentPressure(number);
    }

    Real getFishVelocityX(int number)
    {
        return fish_velocity_probe.getCurrentVelocityX(number);
    }

    Real getFishVelocityY(int number)
    {
        return fish_velocity_probe.getCurrentVelocityY(number);
    }

    Real getFishPositionX(int number)
    {
        return fish_position_probe.getCurrentPositionX(number);
    }

    Real getFishPositionY(int number)
    {
        return fish_position_probe.getCurrentPositionY(number);
    }

    void runCase(int episode, Real pause_time_from_python)
    {
        //----------------------------------------------------------------------
        //	Main loop starts here.
        //----------------------------------------------------------------------
        while (GlobalStaticVariables::physical_time_ < pause_time_from_python)
        {
            Real integration_time = 0.0;

            /** Integrate time (loop) until the next output time. */
            while (integration_time < D_Time)
            {
                initialize_a_fluid_step.exec();
                Real Dt = get_fluid_advection_time_step_size.exec();
                free_stream_surface_indicator.exec();
                update_fluid_density.exec();
                viscous_acceleration.exec();
                transport_velocity_correction.exec();

                /** FSI for viscous force. */
                viscous_force_on_fish.exec();
                viscous_force_on_flap.exec();
                /** Update normal direction on elastic body.*/
                fish_body_update_normal.exec();
                flap_body_update_normal.exec();
                size_t inner_ite_dt = 0;
                size_t inner_ite_dt_s = 0;
                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    Real dt = get_fluid_time_step_size.exec();
                    /** Fluid pressure relaxation, first half. */
                    pressure_relaxation.exec(dt);
                    /** FSI for fluid force on solid body. */
                    fluid_force_on_fish_update.exec();
                    fluid_force_on_flap_update.exec();
                    /** Fluid pressure relaxation, second half. */
                    density_relaxation.exec(dt);

                    /** Solid dynamics. */
                    inner_ite_dt_s = 0;
                    Real dt_s_sum = 0.0;
                    average_fish_velocity_and_acceleration.initialize_displacement_.exec();
                    while (dt_s_sum < dt)
                    {
                        Real dt_s = SMIN(fish_body_computing_time_step_size.exec(), dt - dt_s_sum);
                        imposing_active_strain.exec();
                        fish_body_stress_relaxation_first_half.exec(dt_s);
                        fish_body_stress_relaxation_second_half.exec(dt_s);
                        dt_s_sum += dt_s;
                        inner_ite_dt_s++;
                    }
                    average_fish_velocity_and_acceleration.update_averages_.exec(dt);

                    flap_velocity.exec(dt);

                    relaxation_time += dt;
                    integration_time += dt;
                    GlobalStaticVariables::physical_time_ += dt;
                    emitter_buffer_inflow_condition.exec(dt);
                    inner_ite_dt++;
                }

                fish_pressure_probe.writeToFile(number_of_iterations);
                fish_velocity_probe.writeToFile(number_of_iterations);
                fish_position_probe.writeToFile(number_of_iterations);

                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                            << GlobalStaticVariables::physical_time_
                            << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
                }
                number_of_iterations++;

                /** Water block configuration and periodic condition. */
                emitter_inflow_injection.exec();
                disposer_outflow_deletion.exec();

                water_block.updateCellLinkedListWithParticleSort(100);
                fish_body.updateCellLinkedList();
                flap_body.updateCellLinkedList();
                /** one need update configuration after periodic condition. */
                water_block_complex.updateConfiguration();
                /** one need update configuration after periodic condition. */
                fish_contact.updateConfiguration();
                /** one need update configuration after periodic condition. */
                flap_contact.updateConfiguration();
            }

            TickCount t2 = TickCount::now();
            /** write run-time observation into file */
            compute_vorticity.exec();
            write_real_body_states.writeToFile(number_of_iterations);
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
    }
};

PYBIND11_MODULE(test_2d_flow_stream_around_fish_flap_pybind, m)
{
	py::class_<SphFishRelaxation>(m, "from_sph_relaxation")
        .def(py::init<const int&>());
    
    py::class_<SphFishReloadTrain>(m, "from_sph_reload_and_train")
		.def(py::init<const int&>())
        .def("SetLambda", &SphFishReloadTrain::setLambdafromPython)
        .def("SetFreq", &SphFishReloadTrain::setFreqfromPython)
        .def("GetFishPressurePoint", &SphFishReloadTrain::getFishPressure)
        .def("GetFishVelocityX", &SphFishReloadTrain::getFishVelocityX)
        .def("GetFishVelocityY", &SphFishReloadTrain::getFishVelocityY)
        .def("GetFishPositionX", &SphFishReloadTrain::getFishPositionX)
        .def("GetFishPositionY", &SphFishReloadTrain::getFishPositionY)
		.def("RunCase", &SphFishReloadTrain::runCase);
}