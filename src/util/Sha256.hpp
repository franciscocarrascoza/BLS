// Created: 2026-09-26T16:40+02:00 | by: CC audit | purpose: SHA-256 (FIPS 180-4) for the per-frame grid fingerprint
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace bls {

// Lower-case hex SHA-256 of `len` bytes. Written from FIPS 180-4; checked against the
// standard test vectors in tests/audit/test_audit_harness.cpp.
std::string sha256Hex(const std::uint8_t* data, std::size_t len);

// SHA-256 of the occupancy BITS: voxel i (linear index, x-major, z fastest) is bit
// (i % 8) of byte (i / 8), least significant bit first; any non-zero byte counts as 1.
// Two methods label the same grid only if this, the dims, the per-axis spacing, the
// origin and the periodic flags agree (brief §5 "grid fingerprint").
std::string occupancySha256(const std::vector<std::uint8_t>& occupancy);

}  // namespace bls
