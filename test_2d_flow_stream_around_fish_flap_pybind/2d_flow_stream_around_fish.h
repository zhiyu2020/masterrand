#include "2d_fish_and_bones.h"
#include "active_model.h"
#include "sphinxsys.h"
#define PI 3.14159265359

using namespace SPH;
class SphBasicGeometrySetting : public SphFish, public SphFlap
{
  public:
    //----------------------------------------------------------------------
    //	Basic geometry parameters and numerical setup.
    //----------------------------------------------------------------------
    Real DL = 0.64;                                /**< Channel length. */
    Real DH = 0.4;                                /**< Channel height. */
    Real particle_spacing_ref = 0.004;           /**< Initial reference particle spacing. */
    Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose inflow condition. */
    Real BW = particle_spacing_ref * 4.0;         /**< Boundary width, determined by specific layer of boundary particles. */

    Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
    Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
    Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
    Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
    Vec2d emitter_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
    Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;
    Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
    Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;

    Real cx = 0.3 * DL;           /**< Center of fish in x direction. */
    Real cy = DH / 4;             /**< Center of fish in y direction. */
    Real fish_length = 0.2;       /**< Length of fish. */
    Real fish_thickness = 0.03;   /**< The maximum fish thickness. */
    Real muscle_thickness = 0.022; /**< The maximum fish thickness. */
    Real head_length = 0.03;      /**< Length of fish bone. */
    Real bone_thickness = 0.003;  /**< Length of fish bone. */
    Real fish_shape_resolution = particle_spacing_ref * 0.5;
    Real frequency_ = 3.0;
    Real lambda_ = 2.0;

    //----------------------------------------------------------------------
    //	Material properties of the fluid.
    //----------------------------------------------------------------------
    Real rho0_f = 1000.0;                /**< Density. */
    Real U_f = 1.0;                      /**< freestream velocity. */
    Real c_f = 10.0 * U_f;               /**< Speed of sound. */
    Real mu_f = 0.001;                   /**< Dynamics viscosity. */
    //----------------------------------------------------------------------
    //	Global parameters on the solid properties
    //----------------------------------------------------------------------
    Real rho0_s = 1050.0;
    Real Youngs_modulus1 = 0.8e6;
    Real Youngs_modulus2 = 0.5e6;
    Real Youngs_modulus3 = 1.1e6;
    Real poisson = 0.49;
    //----------------------------------------------------------------------
    //	Global parameters on the flap
    //----------------------------------------------------------------------
    Real gap = 0.13;
    Real flap_cord = 0.067;
    Real flap_thickness = 0.0081;
    Real flap_frequency = 2.0;
    Real flap_yaw_amplitude = 20.0 * PI / 180.0 ;
    Real flap_sway_amplitude = 0.01;
    Real flap_shape_resolution = particle_spacing_ref * 0.5;
    Real flap_cx = cx - gap - flap_cord;
    Real flap_cy = cy + flap_sway_amplitude + 0.2 * DH;
    Real flap_rotation_center_x = flap_cx + 0.019;
    Real flap_rotation_center_y = flap_cy;

    Real a1 = 1.22 * fish_thickness / fish_length;
    Real a2 = 3.19 * fish_thickness / fish_length / fish_length;
    Real a3 = -15.73 * fish_thickness / pow(fish_length, 3);
    Real a4 = 21.87 * fish_thickness / pow(fish_length, 4);
    Real a5 = -10.55 * fish_thickness / pow(fish_length, 5);

    Real b1 = 1.22 * muscle_thickness / fish_length;
    Real b2 = 3.19 * muscle_thickness / fish_length / fish_length;
    Real b3 = -15.73 * muscle_thickness / pow(fish_length, 3);
    Real b4 = 21.87 * muscle_thickness / pow(fish_length, 4);
    Real b5 = -10.55 * muscle_thickness / pow(fish_length, 5);

    Real obx[16] = {0.05, 0.1, 0.15, 0.2,
                    0.25, 0.3, 0.35, 0.4,
                    0.5, 0.6, 0.7, 0.75, 
                    0.8, 0.85, 0.9, 0.95};
    StdVec<Vec2d> observation_locations;
    SphBasicGeometrySetting() : observation_locations()
    {
        for (size_t i = 0; i < 16; ++i)
        {
            Real oby = RevertOutline(obx[i] * fish_length, fish_thickness, fish_length) + cy;
            Real oby_neg = -RevertOutline(obx[i] * fish_length, fish_thickness, fish_length) + cy;

            observation_locations.push_back(Vec2d(cx + obx[i] * fish_length, oby));
            observation_locations.push_back(Vec2d(cx + obx[i] * fish_length, oby_neg));
        }
    }

    std::vector<Vecd> createWaterBlockShape()
    {
        // geometry
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
        water_block_shape.push_back(Vecd(-DL_sponge, DH));
        water_block_shape.push_back(Vecd(DL, DH));
        water_block_shape.push_back(Vecd(DL, 0.0));
        water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

        return water_block_shape;
    }
};

//----------------------------------------------------------------------
//	SPH bodies with cases dependent geometries.
//----------------------------------------------------------------------
class FishBody : public MultiPolygonShape, public SphBasicGeometrySetting
{
  public:
    explicit FishBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> fish_shape = CreatFishShape(cx, cy, fish_length, fish_thickness, fish_shape_resolution);
        multi_polygon_.addAPolygon(fish_shape, ShapeBooleanOps::add);
    }
};

class FlapBody : public MultiPolygonShape, public SphBasicGeometrySetting
{
  public:
    explicit FlapBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> flap_shape = CreatFlapShape(flap_cx, flap_cy, flap_cord, flap_thickness, flap_shape_resolution);
        multi_polygon_.addAPolygon(flap_shape, ShapeBooleanOps::add);
    }
};

class WaterBlock : public ComplexShape, public SphBasicGeometrySetting
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon fish(CreatFishShape(cx, cy, fish_length, fish_thickness, fish_shape_resolution));
        subtract<MultiPolygonShape>(fish);
        MultiPolygon flap(CreatFlapShape(flap_cx, flap_cy, flap_cord, flap_thickness, flap_shape_resolution));
        subtract<MultiPolygonShape>(flap);
    }
};
//----------------------------------------------------------------------
//	Define case dependent bodies material, constraint and boundary conditions.
//----------------------------------------------------------------------
struct FreeStreamVelocity
{
    Real u_ref_, t_ref_;

    template <class BoundaryConditionType>
    FreeStreamVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(1.0), t_ref_(5) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;
        target_velocity[0] = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
    Real t_ref_, u_ref_, du_ave_dt_;

  public:
    explicit TimeDependentAcceleration(Vecd gravity_vector)
        : Gravity(gravity_vector), t_ref_(5), u_ref_(1.0), du_ave_dt_(0) {}

    virtual Vecd InducedAcceleration(const Vecd &position) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);
        return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
    }
};
//----------------------------------------------------------------------
//	Case dependent composite material
//----------------------------------------------------------------------
class FishBodyComposite : public SphBasicGeometrySetting, public CompositeSolid
{
  public:
    FishBodyComposite() : CompositeSolid(rho0_s)
    {
        add<ActiveModelSolid>(rho0_s, Youngs_modulus1, poisson);
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus2, poisson);
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus3, poisson);
    };
};

class FlapBodyComposite : public SphBasicGeometrySetting, public CompositeSolid
{
  public:
    FlapBodyComposite() : CompositeSolid(rho0_s)
    {
        add<ActiveModelSolid>(rho0_s, Youngs_modulus1, poisson);
    };
};
//----------------------------------------------------------------------
//	Case dependent initialization material ids
//----------------------------------------------------------------------
class FishMaterialInitialization
    : public MaterialIdInitialization, public SphBasicGeometrySetting
{
  public:
    explicit FishMaterialInitialization(SolidBody &solid_body)
        : MaterialIdInitialization(solid_body){};

    void update(size_t index_i, Real dt = 0.0)
    {
        Real x = pos0_[index_i][0] - cx;
        Real y = pos0_[index_i][1];

        Real x1 = abs(pos0_[index_i][0] - (cx + fish_length));
        Real y1 = a1 * pow(x1, 0 + 1) + a2 * pow(x1, 1 + 1) + a3 * pow(x1, 2 + 1) + a4 * pow(x1, 3 + 1) + a5 * pow(x1, 4 + 1);
        if (x >= head_length && y > (y1 - 0.004 + cy) && y > (cy + bone_thickness / 2))
        {
            material_id_[index_i] = 0; // region for muscle
        }
        else if (x >= head_length && y < (-y1 + 0.004 + cy) && y < (cy - bone_thickness / 2))
        {
            material_id_[index_i] = 0; // region for muscle
        }
        else if ((x < head_length) || ((y < (cy + bone_thickness / 2)) && (y > (cy - bone_thickness / 2))))
        {
            material_id_[index_i] = 2;
        }
        else
        {
            material_id_[index_i] = 1;
        }
    };
};

class FlapMaterialInitialization
    : public MaterialIdInitialization, public SphBasicGeometrySetting
{
  public:
    explicit FlapMaterialInitialization(SolidBody &solid_body)
        : MaterialIdInitialization(solid_body){};

    void update(size_t index_i, Real dt = 0.0)
    {
        material_id_[index_i] = 0; // region for muscle
    };
};
//----------------------------------------------------------------------
//	imposing active strain to fish muscle
//----------------------------------------------------------------------
class ImposingActiveStrain : public solid_dynamics::ElasticDynamicsInitialCondition, public SphBasicGeometrySetting
{
  public:
    explicit ImposingActiveStrain(SolidBody &solid_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(solid_body),
          material_id_(*particles_->getVariableByName<int>("MaterialID")),
          pos0_(*particles_->getVariableByName<Vecd>("InitialPosition")),
          active_strain_(*particles_->getVariableByName<Matd>("ActiveStrain")){};
    virtual void update(size_t index_i, Real dt = 0.0)
    {
        if (material_id_[index_i] == 0)
        {
            Real x = pos0_[index_i][0] - cx;
            Real y = pos0_[index_i][1];

            Real Am = 0.12;
            Real frequency = frequency_;
            Real w = 2 * Pi * frequency;
            Real lambda = lambda_ * fish_length;
            Real wave_number = 2 * Pi / lambda;
            Real hx = -(pow(x, 2) - pow(fish_length, 2)) / pow(fish_length, 2);
            Real start_time = 0.2;
            Real current_time = GlobalStaticVariables::physical_time_;
            Real strength = 1 - exp(-current_time / start_time);

            Real phase_shift = y > (cy + bone_thickness / 2) ? 0 : Pi / 2;
            active_strain_[index_i](0, 0) =
                -Am * hx * strength * pow(sin(w * current_time / 2 + wave_number * x / 2 + phase_shift), 2);
        }
    };

  protected:
    StdLargeVec<int> &material_id_;
    StdLargeVec<Vecd> &pos0_;
    StdLargeVec<Matd> &active_strain_;
};

class FlapInitialVelocity : public solid_dynamics::ElasticDynamicsInitialCondition, public SphBasicGeometrySetting
{
  public:
    explicit FlapInitialVelocity(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body),
          material_id_(*particles_->getVariableByName<int>("MaterialID")),
          pos0_(*particles_->getVariableByName<Vecd>("InitialPosition")){};

    void update(size_t index_i, Real dt)
    {
        Real flap_angle = flap_yaw_amplitude * sin(2 * PI * flap_frequency * GlobalStaticVariables::physical_time_);
        Real flap_distance = flap_sway_amplitude * sin(2 * PI * flap_frequency * GlobalStaticVariables::physical_time_ + PI / 2) - flap_sway_amplitude;
        //pos_[index_i][1] = (pos0_[index_i][0] - flap_rotation_center_x) * sin(flap_angle) + (pos0_[index_i][1] - flap_rotation_center_y) * cos(flap_angle) 
        //                    + flap_rotation_center_y + flap_distance;
        pos_[index_i][1] = (pos0_[index_i][0] - flap_rotation_center_x) * sin(flap_angle) + (pos0_[index_i][1] - flap_rotation_center_y) * cos(flap_angle) + flap_rotation_center_y;
        pos_[index_i][0] = (pos0_[index_i][0] - flap_rotation_center_x) * cos(flap_angle) - (pos0_[index_i][1] - flap_rotation_center_y) * sin(flap_angle) 
                            + flap_rotation_center_x;
    };

  protected:
    StdLargeVec<int> &material_id_;
    StdLargeVec<Vecd> &pos0_;
};
