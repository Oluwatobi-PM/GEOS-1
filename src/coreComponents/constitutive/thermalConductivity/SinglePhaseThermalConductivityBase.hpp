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
 * @file SinglePhaseThermalConductivityBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYBASE_HPP
#define GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief The abstract base class to perform the thermal conductivity
 */
class SinglePhaseThermalConductivityBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_effectiveConductivity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_effectiveConductivity.size( 1 ); }

protected:

  /**
   * @brief Constructor for the class performing the thermal conductivity updates
   * @param effectiveConductivity the array of cell-wise effective conductivities in the subregion
   * @param dEffectiveConductivity_dPhaseVolFrac the array of cell-wise derivatives of effective conductivities wrt phase vol fractions in
   * the subregion
   */
  SinglePhaseThermalConductivityBaseUpdate( arrayView3d< real64 > const & effectiveConductivity )
    : m_effectiveConductivity( effectiveConductivity )
  {}

  /// View on the cell-wise effective conductivities
  arrayView3d< real64 > m_effectiveConductivity;

private:

  /**
   * @brief Pointwise update function called from the solver
   * @param[in] k index of the cell in the subRegion
   * @param[in] q constitutive index (equal to one in this class)
   * @param[in] laggedPorosity lagged porosity in the cell (for fractures, this will be unused)
   */
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & laggedPorosity ) const = 0;
};

/**
 * @brief The abstract base class for thermal conductivity
 */
class SinglePhaseThermalConductivityBase : public ConstitutiveBase
{
public:

  /**
   * @brief Constructor for the abstract base class
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  SinglePhaseThermalConductivityBase( string const & name, dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Initialize the thermal conductivity state (needed when thermal conductivity depends on porosity and phase volume fraction)
   * @param[in] initialPorosity the initial porosity field after reservoir initialization
   *
   * Note: this is needed because for now, the porosity and phase volume fraction are treated **explictly**
   */
  virtual void initializeRockFluidState( arrayView2d< real64 const > const & initialPorosity ) const
  { GEOS_UNUSED_VAR( initialPorosity ); }

  /**
   * @brief Save the thermal conductivity state (needed when thermal conductivity depends on porosity and phase volume fraction)
   * @param[in] convergedPorosity the converged porosity field after reservoir initialization
   * @param[in] convergedPhaseVolumeFraction the converged phase volume fraction field
   *
   * Note: this is needed because for now, the porosity and phase volume fraction are treated **explictly**
   */
  virtual void saveConvergedRockFluidState( arrayView2d< real64 const > const & convergedPorosity ) const
  { GEOS_UNUSED_VAR( convergedPorosity ); }

  /**
   * @brief Update the thermal conductivity state
   * @param[in] porosity the  porosity field after reservoir initialization
   *
   * Note: this is needed because of the fracture subregions which do not exist at initialization
   */
  virtual void update( arrayView2d< real64 const > const & porosity ) const
  { GEOS_UNUSED_VAR( porosity ); }

  /**
   * @brief Getter for the effective conductivities in the subRegion
   * @return an arrayView of effective conductivities
   */
  arrayView3d< real64 const > effectiveConductivity() const { return m_effectiveConductivity; }

private:

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void resizeFields( localIndex const size, localIndex const numPts );

protected:

  virtual void postProcessInput() override;

  /// cell-wise effective conductivities in the subregion
  array3d< real64 > m_effectiveConductivity;

};

} // namespace constitutive

} // namespace geos


#endif //GEOS_CONSTITUTIVE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYBASE_HPP
