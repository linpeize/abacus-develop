#include "DmkSpinTransform.h"

#include "module_base/parallel_2d.h"
#include "module_base/tool_quit.h"

#include <array>
#include <cstddef>
#include <utility>

namespace
{

using Complex = std::complex<double>;
using SpinorOffsets = std::array<std::size_t, 4>;

void require_transform_shape(const bool condition, const char* message)
{
    if (!condition)
    {
        ModuleBase::WARNING_QUIT("DmkSpinTransform", message);
    }
}

std::size_t flat_offset(const int local_row,
                        const int local_col,
                        const Parallel_2D& pv,
                        const bool column_major)
{
    const std::size_t row = static_cast<std::size_t>(local_row);
    const std::size_t col = static_cast<std::size_t>(local_col);
    const std::size_t row_size = static_cast<std::size_t>(pv.get_row_size());
    const std::size_t col_size = static_cast<std::size_t>(pv.get_col_size());
    return column_major ? col * row_size + row : row * col_size + col;
}

std::vector<SpinorOffsets> build_spinor_offsets(const Parallel_2D& pv, const bool column_major)
{
    const int global_rows = pv.get_global_row_size();
    const int global_cols = pv.get_global_col_size();
    require_transform_shape(global_rows >= 0 && global_cols >= 0 && global_rows % 2 == 0
                                && global_cols % 2 == 0,
                            "nspin=4 DMK requires even global matrix dimensions");
    require_transform_shape(pv.get_row_size() >= 0 && pv.get_col_size() >= 0 && pv.get_local_size() >= 0,
                            "nspin=4 DMK has an invalid local matrix shape");

    std::vector<std::pair<int, int>> local_row_pairs;
    for (int mu = 0; mu < global_rows / 2; ++mu)
    {
        const int up = pv.global2local_row(2 * mu);
        const int down = pv.global2local_row(2 * mu + 1);
        const bool up_is_local = up >= 0;
        const bool down_is_local = down >= 0;
        require_transform_shape(up_is_local == down_is_local,
                                "nspin=4 row spin pair crosses MPI ranks");
        if (up_is_local)
        {
            local_row_pairs.push_back(std::make_pair(up, down));
        }
    }

    std::vector<std::pair<int, int>> local_col_pairs;
    for (int nu = 0; nu < global_cols / 2; ++nu)
    {
        const int up = pv.global2local_col(2 * nu);
        const int down = pv.global2local_col(2 * nu + 1);
        const bool up_is_local = up >= 0;
        const bool down_is_local = down >= 0;
        require_transform_shape(up_is_local == down_is_local,
                                "nspin=4 column spin pair crosses MPI ranks");
        if (up_is_local)
        {
            local_col_pairs.push_back(std::make_pair(up, down));
        }
    }

    const std::size_t local_size = static_cast<std::size_t>(pv.get_local_size());
    std::vector<SpinorOffsets> offsets;
    offsets.reserve(local_size / 4);

    // Keep the packed spatial-matrix order consistent with the KS solver.
    if (column_major)
    {
        for (const std::pair<int, int>& col : local_col_pairs)
        {
            for (const std::pair<int, int>& row : local_row_pairs)
            {
                offsets.push_back({flat_offset(row.first, col.first, pv, true),
                                   flat_offset(row.first, col.second, pv, true),
                                   flat_offset(row.second, col.first, pv, true),
                                   flat_offset(row.second, col.second, pv, true)});
            }
        }
    }
    else
    {
        for (const std::pair<int, int>& row : local_row_pairs)
        {
            for (const std::pair<int, int>& col : local_col_pairs)
            {
                offsets.push_back({flat_offset(row.first, col.first, pv, false),
                                   flat_offset(row.first, col.second, pv, false),
                                   flat_offset(row.second, col.first, pv, false),
                                   flat_offset(row.second, col.second, pv, false)});
            }
        }
    }

    require_transform_shape(offsets.size() * std::size_t{4} == local_size,
                            "nspin=4 spin quartets do not cover the local matrix");

    std::vector<bool> covered(local_size, false);
    for (const SpinorOffsets& quartet : offsets)
    {
        for (const std::size_t offset : quartet)
        {
            require_transform_shape(offset < local_size && !covered[offset],
                                    "nspin=4 spin quartet offsets are invalid or duplicated");
            covered[offset] = true;
        }
    }
    return offsets;
}

template <typename Tdata>
std::size_t validate_nspin2_physical(const std::vector<std::vector<Tdata>>& physical)
{
    require_transform_shape(!physical.empty() && physical.size() % 2 == 0,
                            "nspin=2 physical DMK requires complete up/down stream pairs");
    const std::size_t local_size = physical.front().size();
    for (const std::vector<Tdata>& stream : physical)
    {
        require_transform_shape(stream.size() == local_size,
                                "nspin=2 physical DMK streams must have equal lengths");
    }
    return local_size;
}

template <typename Tdata>
std::size_t validate_nspin2_pauli(const std::vector<const std::vector<Tdata>*>& pauli)
{
    require_transform_shape(!pauli.empty(), "nspin=2 Pauli DMK requires at least one stream");
    for (const std::vector<Tdata>* stream : pauli)
    {
        require_transform_shape(stream != nullptr, "nspin=2 Pauli DMK contains a null stream");
    }
    const std::size_t packed_size = pauli.front()->size();
    require_transform_shape(packed_size % 2 == 0,
                            "nspin=2 Pauli DMK stream length must be divisible by two");
    for (const std::vector<Tdata>* stream : pauli)
    {
        require_transform_shape(stream->size() == packed_size,
                                "nspin=2 Pauli DMK streams must have equal lengths");
    }
    return packed_size / 2;
}

void validate_nspin4_physical(const std::vector<std::vector<Complex>>& physical,
                              const std::size_t local_size)
{
    require_transform_shape(!physical.empty(), "nspin=4 physical DMK requires at least one stream");
    for (const std::vector<Complex>& stream : physical)
    {
        require_transform_shape(stream.size() == local_size,
                                "nspin=4 physical DMK stream length does not match Parallel_2D");
    }
}

void validate_nspin4_pauli(const std::vector<const std::vector<Complex>*>& pauli,
                           const std::size_t local_size)
{
    require_transform_shape(!pauli.empty(), "nspin=4 Pauli DMK requires at least one stream");
    for (const std::vector<Complex>* stream : pauli)
    {
        require_transform_shape(stream != nullptr, "nspin=4 Pauli DMK contains a null stream");
        require_transform_shape(stream->size() == local_size,
                                "nspin=4 Pauli DMK stream length does not match Parallel_2D");
    }
}

} // namespace

namespace DmkSpinTransform
{

template <typename Tdata>
std::vector<std::vector<Tdata>> to_pauli_nspin2(const std::vector<std::vector<Tdata>>& physical)
{
    const std::size_t component_size = validate_nspin2_physical(physical);
    const std::size_t nk = physical.size() / 2;
    std::vector<std::vector<Tdata>> pauli(nk, std::vector<Tdata>(2 * component_size));

    for (std::size_t ik = 0; ik < nk; ++ik)
    {
        const std::vector<Tdata>& up = physical[ik];
        const std::vector<Tdata>& down = physical[ik + nk];
        for (std::size_t i = 0; i < component_size; ++i)
        {
            pauli[ik][i] = up[i] + down[i];
            pauli[ik][component_size + i] = up[i] - down[i];
        }
    }
    return pauli;
}

template <typename Tdata>
std::vector<std::vector<Tdata>> to_physical_nspin2(
    const std::vector<const std::vector<Tdata>*>& pauli)
{
    const std::size_t component_size = validate_nspin2_pauli(pauli);
    const std::size_t nk = pauli.size();
    std::vector<std::vector<Tdata>> physical(2 * nk, std::vector<Tdata>(component_size));

    for (std::size_t ik = 0; ik < nk; ++ik)
    {
        const std::vector<Tdata>& packed = *pauli[ik];
        for (std::size_t i = 0; i < component_size; ++i)
        {
            physical[ik][i] = Tdata(0.5) * (packed[i] + packed[component_size + i]);
            physical[ik + nk][i] = Tdata(0.5) * (packed[i] - packed[component_size + i]);
        }
    }
    return physical;
}

std::vector<std::vector<Complex>> to_pauli_nspin4(const std::vector<std::vector<Complex>>& physical,
                                                   const Parallel_2D& pv,
                                                   const bool column_major)
{
    const std::vector<SpinorOffsets> offsets = build_spinor_offsets(pv, column_major);
    const std::size_t local_size = static_cast<std::size_t>(pv.get_local_size());
    validate_nspin4_physical(physical, local_size);

    const std::size_t component_size = offsets.size();
    const Complex imaginary(0.0, 1.0);
    std::vector<std::vector<Complex>> pauli(physical.size(), std::vector<Complex>(local_size));
    for (std::size_t ik = 0; ik < physical.size(); ++ik)
    {
        for (std::size_t i = 0; i < component_size; ++i)
        {
            const SpinorOffsets& quartet = offsets[i];
            const Complex& uu = physical[ik][quartet[0]];
            const Complex& ud = physical[ik][quartet[1]];
            const Complex& du = physical[ik][quartet[2]];
            const Complex& dd = physical[ik][quartet[3]];
            pauli[ik][i] = uu + dd;
            pauli[ik][component_size + i] = ud + du;
            pauli[ik][2 * component_size + i] = imaginary * (ud - du);
            pauli[ik][3 * component_size + i] = uu - dd;
        }
    }
    return pauli;
}

std::vector<std::vector<Complex>> to_physical_nspin4(
    const std::vector<const std::vector<Complex>*>& pauli,
    const Parallel_2D& pv,
    const bool column_major)
{
    const std::vector<SpinorOffsets> offsets = build_spinor_offsets(pv, column_major);
    const std::size_t local_size = static_cast<std::size_t>(pv.get_local_size());
    validate_nspin4_pauli(pauli, local_size);

    const std::size_t component_size = offsets.size();
    const Complex imaginary(0.0, 1.0);
    std::vector<std::vector<Complex>> physical(pauli.size(), std::vector<Complex>(local_size));
    for (std::size_t ik = 0; ik < pauli.size(); ++ik)
    {
        const std::vector<Complex>& packed = *pauli[ik];
        for (std::size_t i = 0; i < component_size; ++i)
        {
            const Complex& d0 = packed[i];
            const Complex& dx = packed[component_size + i];
            const Complex& dy = packed[2 * component_size + i];
            const Complex& dz = packed[3 * component_size + i];
            const SpinorOffsets& quartet = offsets[i];
            physical[ik][quartet[0]] = 0.5 * (d0 + dz);
            physical[ik][quartet[1]] = 0.5 * (dx - imaginary * dy);
            physical[ik][quartet[2]] = 0.5 * (dx + imaginary * dy);
            physical[ik][quartet[3]] = 0.5 * (d0 - dz);
        }
    }
    return physical;
}

template std::vector<std::vector<double>> to_pauli_nspin2(
    const std::vector<std::vector<double>>& physical);
template std::vector<std::vector<Complex>> to_pauli_nspin2(
    const std::vector<std::vector<Complex>>& physical);
template std::vector<std::vector<double>> to_physical_nspin2(
    const std::vector<const std::vector<double>*>& pauli);
template std::vector<std::vector<Complex>> to_physical_nspin2(
    const std::vector<const std::vector<Complex>*>& pauli);

} // namespace DmkSpinTransform
