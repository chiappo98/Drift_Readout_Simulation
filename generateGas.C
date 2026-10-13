// ============================================================
// Ex1_GenerateGas.C  [Student Version]
//
// PURPOSE:
//   Generate Magboltz gas transport tables for two Ar-CO2
//   mixtures at atmospheric pressure and room temperature:
//     - Ar 90% + CO2 10%   [Worked Example]
//     - Ar 70% + CO2 30%   [YOUR TASK - 5 TODOs]
//
//   No magnetic field is included; this is the baseline case.
//   The resulting .gas files are prerequisites for Ex3 and Ex6.
//
// OUTPUT FILES:
//   ar_co2_90_10_Production.gas   (worked example)
//   ar_co2_70_30_Production.gas   (your output)
//
// DIFFERENT GAS FILE VERSIONS:
//   For Demo, nE = 10, ncoll = 5.
//   For Production, nE = 20, ncoll = 10.
//
// TIMING NOTE:
//   With ncoll=10, each mixture takes ~90 minutes on a desktop.
//   Pre-generated .gas files will be provided for the school;
//   you do not need to run this during the session.
//   The exercise is to understand and write the code.
//
// USAGE:
//   ./Ex1_GenerateGas
// ============================================================

#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/MediumMagboltz.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {

  // --- Conditions -------------------------------------------
  const double pressure    = AtmosphericPressure; // 760 Torr
  const double temperature = 298.15;              // 20 deg C in Kelvin

  // --- Electric field grid ----------------------------------
  // For Demo, nE = 10, for Production, nE = 20.
  // Logarithmically spaced: 100 V/cm to 100 kV/cm.
  const size_t nE      = 20;
  const double eMin    = 10.;    // V/cm
  const double eMax    = 4e4; // V/cm
  constexpr bool useLog = true;

  // --- Magboltz collision count -----------------------------
  // ncoll * 1e7 collisions per field point.
  // For Demo, ncoll = 5, for Production, ncoll = 10.
  const int ncoll = 10;

  // ==========================================================
  // Mixture 1: Ar 90% + CO2 10%  [Worked Example]
  // ==========================================================
  std::cout << "\n=== Generating table for Ar/CO2 90:10 ===\n";
  MediumMagboltz gas1("Ar", 90., "CO2", 10.);
  gas1.SetTemperature(temperature);
  gas1.SetPressure(pressure);
  gas1.SetFieldGrid(eMin, eMax, nE, useLog);
  gas1.GenerateGasTable(ncoll);
  gas1.WriteGasFile("ar_co2_90_10.gas");
  std::cout << "Written: ar_co2_90_10.gas\n";

  // ==========================================================
  // Mixture 2: Ar 70% + CO2 30%  [YOUR TASK]
  // Study the worked example above, then complete the five
  // steps below for the second mixture.
  // ==========================================================
  std::cout << "\n=== Generating table for Ar/CO2 70:30 ===\n";

  MediumMagboltz gas2("Ar", 70., "CO2", 30.);
  gas2.SetTemperature(temperature);
  gas2.SetPressure(pressure);
  gas2.SetFieldGrid(eMin, eMax, nE, useLog);
  gas2.GenerateGasTable(ncoll);
  gas2.WriteGasFile("ar_co2_70_30.gas");

  std::cout << "Written: ar_co2_70_30.gas\n";

  std::cout << "\n=== Generating table for Ar/CO2 85:15 ===\n";

  MediumMagboltz gas3("Ar", 85., "CO2", 15.);
  gas3.SetTemperature(temperature);
  gas3.SetPressure(pressure);
  gas3.SetFieldGrid(eMin, eMax, nE, useLog);
  gas3.GenerateGasTable(ncoll);
  gas3.WriteGasFile("ar_co2_85_15.gas");

  std::cout << "Written: ar_co2_85_15.gas\n";

  return 0;
}
