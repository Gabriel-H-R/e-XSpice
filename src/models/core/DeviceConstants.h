/**
 * @author Gabriel Henrique Silva
 * @file DeviceConstants.h
 * @brief Shared SPICE physical constants (single definition).
 *
 * Values follow Xyce N_DEV_Const.h / SPICE3: electron charge, Boltzmann
 * constant, silicon bandgap vs temperature, reference temperature and the
 * exponent clamp. Centralized here so models (DiodeModel, ...) and devices
 * (Diode, ...) can never disagree on them.
 */

#ifndef DEVICE_CONSTANTS_H
#define DEVICE_CONSTANTS_H

namespace e_XSpice {
namespace Const {

static constexpr double Q        = 1.6021918e-19;   // electron charge (C)
static constexpr double Boltz    = 1.3806226e-23;   // Boltzmann (J/K)
static constexpr double KoverQ   = Boltz / Q;       // 8.6170344e-05 V/K
static constexpr double RefTemp  = 300.15;          // 27 degC in Kelvin
static constexpr double Eg0      = 1.16;            // Si bandgap @ 0 K (eV)
static constexpr double AlphaEg  = 7.02e-4;         // eV/K
static constexpr double BetaEg   = 1108.0;          // K
static constexpr double Eg300    = 1.1150877;       // Si bandgap @ 300.15 K
static constexpr double E        = 2.71828182845904523536028747135;
static constexpr double Root2    = 1.41421356237309504880168872421;
static constexpr double MaxExpArg = 100.0;          // exp() argument clamp
static constexpr double CtoK     = 273.15;          // Celsius -> Kelvin

} // namespace Const
} // namespace e_XSpice

#endif // DEVICE_CONSTANTS_H
