#include "elastic_energy.h"

#include <numeric>

namespace SPH
{
	namespace solid_dynamics
	{
	//=================================================================================================//
	ElasticEnergy::ElasticEnergy(SPHBody& sph_body)
		: LocalDynamicsReduce<Real, ReduceSum<Real>>(sph_body, Real(0)),
		ElasticSolidDataSimple(sph_body), 
		Vol_(particles_->Vol_), F_(particles_->F_),
		stress_PK1_B_(*particles_->getVariableByName<Matd>("CorrectedStressPK1")),
		active_strain_(*particles_->getVariableByName<Matd>("ActiveStrain")),
		materail_id_(*particles_->getVariableByName<int>("MaterialID")),
		active_strain_rate_(*particles_->getVariableByName<Matd>("ActiveStrainRate"))
	{
		quantity_name_ = "ElasticEnergy";
	}
	//=================================================================================================//
	Real ElasticEnergy::reduce(size_t index_i, Real dt)
	{
		Matd active_stress = Matd:: Zero();
		Matd elastic_stress = Matd::Zero();
		Matd F_e = Matd::Zero();
		Real lambda0_(0.0), G0_(0.0);

		Matd F0 = (2.0 * active_strain_[index_i] + Matd::Identity()).llt().matrixL();

		F_e = F_[index_i] * F0.inverse();
		Matd strain = 0.5 * (F_e.transpose() * F_e - Matd::Identity());

		if (materail_id_[index_i] == 0)
		{
			lambda0_ = 13154362.42;
			G0_ = 268456.3758;
		}
		else if(materail_id_[index_i] == 1)
		{
			lambda0_ = 8221476.51;
			G0_ = 167785.2349;
		}
		else
		{
			lambda0_ = 18087248.32;
			G0_ = 369127.5168;
		}

		elastic_stress = F_e * (lambda0_ * strain.trace() * Matd::Identity() + 2.0 * G0_ * strain);
		active_stress = stress_PK1_B_[index_i] - elastic_stress;

		return  0.5 * CalculateBiDotProduct(active_stress, active_strain_rate_[index_i]) * Vol_[index_i];
	}
	//=================================================================================================//
	SolidKinecticEnergy::SolidKinecticEnergy(SPHBody& sph_body)
		: LocalDynamicsReduce<Real, ReduceSum<Real>>(sph_body, Real(0)),
		ElasticSolidDataSimple(sph_body), vel_(particles_->vel_), mass_(particles_->mass_)
	{
		quantity_name_ = "SolidKinecticEnergy";
	}
	//=================================================================================================//
	Real SolidKinecticEnergy::reduce(size_t index_i, Real dt)
	{

		return 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
	}

	}
}
