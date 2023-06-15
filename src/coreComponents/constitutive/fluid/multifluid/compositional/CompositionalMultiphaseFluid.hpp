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
 * @file CompositionalMultiphaseFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUID_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{
namespace constitutive
{

class CompositionalMultiphaseFluid : public MultiFluidBase
{
public:
  using exec_policy = serialPolicy;

  class CompositionalMultiphaseFluidUpdates : public MultiFluidBase::KernelWrapper
  {
public:
    GEOS_HOST_DEVICE
    void compute( real64 const pressure,
                  real64 const temperature,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                  arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseComponentFraction,
                  real64 & totalDensity ) const override;
    GEOS_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          PhaseProp::SliceType const phaseFraction,
                          PhaseProp::SliceType const phaseMolarDensity,
                          PhaseProp::SliceType const phaseMassDensity,
                          PhaseProp::SliceType const phaseViscosity,
                          PhaseProp::SliceType const phaseEnthalpy,
                          PhaseProp::SliceType const phaseInternalEnergy,
                          PhaseComp::SliceType const phaseComponentFraction,
                          FluidProp::SliceType const totalDensity ) const override;
    GEOS_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

    static CompositionalMultiphaseFluidUpdates createKernelWrapper(
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
      FluidProp::ViewType totalDensity );

protected:
    CompositionalMultiphaseFluidUpdates(
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
      FluidProp::ViewType totalDensity );

    arrayView1d< real64 const > m_componentCriticalPressure;
    arrayView1d< real64 const > m_componentCriticalTemperature;
    arrayView1d< real64 const > m_componentAcentricFactor;
    arrayView1d< real64 const > m_componentVolumeShift;
    arrayView2d< real64 const > m_componentBinaryCoeff;
  };

  using KernelWrapper = CompositionalMultiphaseFluidUpdates;

public:
  CompositionalMultiphaseFluid( string const & name, Group * const parent );
  ~CompositionalMultiphaseFluid() override = default;

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName() { return "CompositionalMultiphaseFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual integer getWaterPhaseIndex() const override final;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * equationsOfStateString() { return "equationsOfState"; }
    static constexpr char const * componentCriticalPressureString() { return "componentCriticalPressure"; }
    static constexpr char const * componentCriticalTemperatureString() { return "componentCriticalTemperature"; }
    static constexpr char const * componentAcentricFactorString() { return "componentAcentricFactor"; }
    static constexpr char const * componentVolumeShiftString() { return "componentVolumeShift"; }
    static constexpr char const * componentBinaryCoeffString() { return "componentBinaryCoeff"; }
  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

protected:

  virtual void postProcessInput() override;

  virtual void initializePostSubGroups() override;

  void createFluid();

  class IFluidProps
  {};

  template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR, typename VISC_TYPE_LIQUID, typename VISC_TYPE_VAPOUR >
  class FluidProps;

private:
  // names of equations of state to use for each phase
  string_array m_equationsOfState;

  // standard EOS component input
  array1d< real64 > m_componentCriticalPressure;
  array1d< real64 > m_componentCriticalTemperature;
  array1d< real64 > m_componentAcentricFactor;
  array1d< real64 > m_componentVolumeShift;
  array2d< real64 > m_componentBinaryCoeff;

  /// Fluid
  std::unique_ptr< IFluidProps > m_fluid{};
};

} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_
