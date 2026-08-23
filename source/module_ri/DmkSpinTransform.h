#ifndef DMK_SPIN_TRANSFORM_H
#define DMK_SPIN_TRANSFORM_H

#include <complex>
#include <vector>

class Parallel_2D;

namespace DmkSpinTransform
{

/**
 * Convert spin-major {D_up, D_down} streams to one packed [D0 | Dz]
 * stream per physical k point.
 */
template <typename Tdata>
std::vector<std::vector<Tdata>> to_pauli_nspin2(const std::vector<std::vector<Tdata>>& physical);

/**
 * Convert packed [D0 | Dz] streams back to spin-major {D_up, D_down}.
 */
template <typename Tdata>
std::vector<std::vector<Tdata>> to_physical_nspin2(
    const std::vector<const std::vector<Tdata>*>& pauli);

/**
 * Convert physical spinor blocks {D_uu, D_ud, D_du, D_dd} to packed
 * [D0 | Dx | Dy | Dz] streams without changing the local vector length.
 */
std::vector<std::vector<std::complex<double>>> to_pauli_nspin4(
    const std::vector<std::vector<std::complex<double>>>& physical,
    const Parallel_2D& pv,
    bool column_major);

/**
 * Convert packed [D0 | Dx | Dy | Dz] streams back to physical spinor
 * matrices in the local layout described by pv.
 */
std::vector<std::vector<std::complex<double>>> to_physical_nspin4(
    const std::vector<const std::vector<std::complex<double>>*>& pauli,
    const Parallel_2D& pv,
    bool column_major);

} // namespace DmkSpinTransform

#endif
