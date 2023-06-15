/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseFluidUpdates.cpp
 */

#include "CompositionalMultiphaseFluid.hpp"

namespace geos
{
namespace constitutive
{
CompositionalMultiphaseFluid::CompositionalMultiphaseFluidUpdates
  ::CompositionalMultiphaseFluidUpdates(
  bool const useMass,
  arrayView1d< real64 const > const & componentMolarWeight,
  arrayView1d< real64 const > const & componentCriticalPressure,
  arrayView1d< real64 const > const & componentCriticalTemperature,
  arrayView1d< real64 const > const & componentAcentricFactor,
  arrayView1d< real64 const > const & componentVolumeShift,
  arrayView2d< real64 const > const & componentBinaryCoeff,
  PhaseProp::ViewType phaseFraction,
  PhaseProp::ViewType phaseMolarDensity,
  PhaseProp::ViewType phaseMassDensity,
  PhaseProp::ViewType phaseViscosity,
  PhaseProp::ViewType phaseEnthalpy,
  PhaseProp::ViewType phaseInternalEnergy,
  PhaseComp::ViewType phaseComponentFraction,
  FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( componentMolarWeight,
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseMolarDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseEnthalpy ),
                                   std::move( phaseInternalEnergy ),
                                   std::move( phaseComponentFraction ),
                                   std::move( totalDensity ) ),
  m_componentCriticalPressure( componentCriticalPressure ),
  m_componentCriticalTemperature( componentCriticalTemperature ),
  m_componentAcentricFactor( componentAcentricFactor ),
  m_componentVolumeShift( componentVolumeShift ),
  m_componentBinaryCoeff( componentBinaryCoeff )
{}

CompositionalMultiphaseFluid::CompositionalMultiphaseFluidUpdates
CompositionalMultiphaseFluid::CompositionalMultiphaseFluidUpdates
  ::createKernelWrapper(
  bool const useMass,
  arrayView1d< real64 const > const & componentMolarWeight,
  arrayView1d< real64 const > const & componentCriticalPressure,
  arrayView1d< real64 const > const & componentCriticalTemperature,
  arrayView1d< real64 const > const & componentAcentricFactor,
  arrayView1d< real64 const > const & componentVolumeShift,
  arrayView2d< real64 const > const & componentBinaryCoeff,
  PhaseProp::ViewType phaseFraction,
  PhaseProp::ViewType phaseMolarDensity,
  PhaseProp::ViewType phaseMassDensity,
  PhaseProp::ViewType phaseViscosity,
  PhaseProp::ViewType phaseEnthalpy,
  PhaseProp::ViewType phaseInternalEnergy,
  PhaseComp::ViewType phaseComponentFraction,
  FluidProp::ViewType totalDensity )
{
  return CompositionalMultiphaseFluidUpdates(
    useMass,
    componentMolarWeight,
    componentCriticalPressure,
    componentCriticalTemperature,
    componentAcentricFactor,
    componentVolumeShift,
    componentBinaryCoeff,
    phaseFraction,
    phaseMolarDensity,
    phaseMassDensity,
    phaseViscosity,
    phaseEnthalpy,
    phaseInternalEnergy,
    phaseComponentFraction,
    totalDensity
    );
}

GEOS_HOST_DEVICE
void CompositionalMultiphaseFluid::CompositionalMultiphaseFluidUpdates
  ::compute( real64 const pressure,
             real64 const temperature,
             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMolarDensity,
             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
             arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseComponentFraction,
             real64 & totalDensity ) const
{
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );
  GEOS_UNUSED_VAR( composition );
  GEOS_UNUSED_VAR( phaseEnthalpy );
  GEOS_UNUSED_VAR( phaseInternalEnergy );

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  for( integer ip = 0; ip < numPhase; ++ip )
  {
    real64 const phaseMolecularWeight = 80.0 + ip*40.0;

    phaseFraction[ip] = 0.35;
    phaseMassDensity[ip] = 500.0 + 300.0*ip;
    phaseMolarDensity[ip] = phaseMassDensity[ip] / phaseMolecularWeight;
    phaseViscosity[ip] = 0.001;
    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseComponentFraction[ip][jc] = 1.0 / numComp;
    }
  }

  computeTotalDensity< maxNumComp, maxNumPhase >( phaseFraction,
                                                  phaseMassDensity,
                                                  totalDensity );
}

GEOS_HOST_DEVICE
void CompositionalMultiphaseFluid::CompositionalMultiphaseFluidUpdates
  ::compute( real64 const pressure,
             real64 const temperature,
             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
             PhaseProp::SliceType const phaseFraction,
             PhaseProp::SliceType const phaseMolarDensity,
             PhaseProp::SliceType const phaseMassDensity,
             PhaseProp::SliceType const phaseViscosity,
             PhaseProp::SliceType const phaseEnthalpy,
             PhaseProp::SliceType const phaseInternalEnergy,
             PhaseComp::SliceType const phaseComponentFraction,
             FluidProp::SliceType const totalDensity ) const
{
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );
  GEOS_UNUSED_VAR( composition );
  GEOS_UNUSED_VAR( phaseEnthalpy );
  GEOS_UNUSED_VAR( phaseInternalEnergy );

  using Deriv = multifluid::DerivativeOffset;

  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  for( integer ip = 0; ip < numPhase; ++ip )
  {
    real64 const phaseMolecularWeight = 80.0 + ip*40.0;

    phaseFraction.value[ip] = (1-ip)*0.35 + ip*0.65;
    phaseFraction.derivs[ip][Deriv::dP] = 0.0;
    phaseFraction.derivs[ip][Deriv::dT] = 0.0;

    phaseMolarDensity.value[ip] = 500.0 + 300.0*ip;
    phaseMolarDensity.derivs[ip][Deriv::dP] = 0.0;
    phaseMolarDensity.derivs[ip][Deriv::dT] = 0.0;

    phaseMassDensity.value[ip] = phaseMolarDensity.value[ip] / phaseMolecularWeight;
    phaseMassDensity.derivs[ip][Deriv::dP] = 0.0;
    phaseMassDensity.derivs[ip][Deriv::dT] = 0.0;

    phaseViscosity.value[ip] = 0.001;
    phaseViscosity.derivs[ip][Deriv::dP] = 0.0;
    phaseViscosity.derivs[ip][Deriv::dT] = 0.0;

    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseFraction.derivs[ip][Deriv::dC+jc] = 0.0;
      phaseMolarDensity.derivs[ip][Deriv::dC+jc] = 0.0;
      phaseMassDensity.derivs[ip][Deriv::dC+jc] = 0.0;
      phaseViscosity.derivs[ip][Deriv::dC+jc] = 0.0;

      phaseComponentFraction.value[ip][jc] = 1.0 / numComp;
      phaseComponentFraction.derivs[ip][jc][Deriv::dP] = 0.0;
      phaseComponentFraction.derivs[ip][jc][Deriv::dT] = 0.0;

      for( integer ic = 0; ic < numComp; ++ic )
      {
        phaseComponentFraction.derivs[ip][ic][Deriv::dC+jc] = 0.0;
      }
    }
  }

  computeTotalDensity( phaseFraction,
                       phaseMolarDensity,
                       totalDensity );
}

GEOS_HOST_DEVICE
void CompositionalMultiphaseFluid::CompositionalMultiphaseFluidUpdates
  ::update( localIndex const k,
            localIndex const q,
            real64 const pressure,
            real64 const temperature,
            arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  compute( pressure,
           temperature,
           composition,
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseEnthalpy( k, q ),
           m_phaseInternalEnergy( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

} /* namespace constitutive */
} /* namespace geos */
