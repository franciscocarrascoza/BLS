#pragma once

#include <cstddef>

namespace bls {

std::size_t currentRSSBytes();

// Get total available system RAM in bytes
std::size_t availableSystemRAMBytes();

// Estimate memory required for grid allocation (occupancy + visited arrays)
std::size_t estimateGridMemoryBytes(std::size_t nx, std::size_t ny, std::size_t nz);

// Calculate maximum box dimension for a given grid spacing and RAM limit
double maxBoxDimensionForRAM(double gridSpacing, std::size_t maxRAMBytes);

}  // namespace bls

