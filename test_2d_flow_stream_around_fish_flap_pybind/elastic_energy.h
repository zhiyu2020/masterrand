/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	general_solid_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid dynamics. 
 * @details We consider here a weakly compressible solids.  
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef ELASTIC_ENERGY_H
#define ELASTIC_ENERGY_H

#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "all_body_relations.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "elastic_solid.h"

namespace SPH
{
	namespace solid_dynamics
	{
		//----------------------------------------------------------------------
		//		for general solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<ElasticSolidParticles> ElasticSolidDataSimple;

		/**
        * @class ElasticEnergy
        */
		class ElasticEnergy
			: public LocalDynamicsReduce<Real, ReduceSum<Real>>,
			public ElasticSolidDataSimple
		{
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Matd>& F_, &stress_PK1_B_, &active_strain_, &active_strain_rate_;
			StdLargeVec<int>& materail_id_;

		public:
			ElasticEnergy(SPHBody& sph_body);
			virtual ~ElasticEnergy() {};

			Real reduce(size_t index_i, Real dt = 0.0);
		};

		/**
	    * @class KineticEnergy
	    */
		class SolidKinecticEnergy
			: public LocalDynamicsReduce<Real, ReduceSum<Real>>,
			public ElasticSolidDataSimple
		{
		protected:
			StdLargeVec<Vecd>&vel_;
			StdLargeVec<Real>& mass_;

		public:
			SolidKinecticEnergy(SPHBody& sph_body);
			virtual ~SolidKinecticEnergy() {};

			Real reduce(size_t index_i, Real dt = 0.0);
		};
	}
}
#endif //ELASTIC_ENERGY_H