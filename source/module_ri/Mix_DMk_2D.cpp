//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#include "Mix_DMk_2D.h"
#include "module_base/module_mixing/plain_mixing.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"
#include "module_parameter/parameter.h"

#include <cassert>

template <typename Tdata>
Mix_DMk_2D<Tdata>::~Mix_DMk_2D<Tdata>()
{
    if(this->flag_del_mixing)
        delete this->mixing;
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::set_nks(const int nks)
{
    this->mix_DMk.clear();
    this->mix_DMk.resize(nks);
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::set_mixing(Base_Mixing::Mixing* mixing_in)
{
    if(this->flag_del_mixing)
        delete this->mixing;
    this->mixing = mixing_in;
    this->flag_del_mixing = false;
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::set_mixing_plain(const double mixing_beta)
{
    if(this->flag_del_mixing)
        delete this->mixing;
    this->mixing = new Base_Mixing::Plain_Mixing(mixing_beta);
    this->flag_del_mixing = true;
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::mix(const std::vector<std::vector<Tdata>>& dm, const bool flag_restart)
{
    ModuleBase::TITLE("Mix_DMk_2D", "mix");
    this->validate_input(dm, flag_restart);
    if (flag_restart)
        { this->restart_all(dm); }
    else
        { this->mix_all(dm); }
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::mix(const std::vector<std::vector<Tdata>>& dm,
                            const bool flag_restart,
                            const std::size_t beta_split_index)
{
    ModuleBase::TITLE("Mix_DMk_2D", "mix_split");
    this->validate_input(dm, flag_restart);
    for (const std::vector<Tdata>& stream : dm)
    {
        if (beta_split_index > stream.size())
        {
            ModuleBase::WARNING_QUIT("Mix_DMk_2D", "beta split index exceeds a DMK stream length");
        }
    }

    if (flag_restart)
        { this->restart_all(dm); }
    else
        { this->mix_all(dm, beta_split_index); }
}

template <typename Tdata>
std::vector<const std::vector<Tdata>*> Mix_DMk_2D<Tdata>::get_DMk_out() const
{
    std::vector<const std::vector<Tdata>*> DMk_out(this->mix_DMk.size());
    for (int ik = 0; ik < this->mix_DMk.size(); ++ik)
        { DMk_out[ik] = &this->mix_DMk[ik].data_out; }
    return DMk_out;
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::validate_input(const std::vector<std::vector<Tdata>>& data_in,
                                       const bool flag_restart) const
{
    if (this->mixing == nullptr)
    {
        ModuleBase::WARNING_QUIT("Mix_DMk_2D", "mixing engine is not configured");
    }
    if (data_in.size() != this->mix_DMk.size())
    {
        ModuleBase::WARNING_QUIT("Mix_DMk_2D", "DMK stream count does not match set_nks");
    }

    if (!flag_restart)
    {
        for (std::size_t ik = 0; ik < data_in.size(); ++ik)
        {
            if (data_in[ik].size() != this->mix_DMk[ik].data_out.size()
                || data_in[ik].size() != this->mix_DMk[ik].mixing_data.length)
            {
                ModuleBase::WARNING_QUIT("Mix_DMk_2D", "DMK stream shape changed without resetting history");
            }
        }
    }
    else
    {
        for (std::size_t ik = 0; ik < data_in.size(); ++ik)
        {
            // Mixing_Data::resize cannot safely shrink an allocated history to
            // zero length. A rank with no local elements is supported from the
            // initial seed; changing rank ownership requires rebuilding set_nks.
            if (data_in[ik].empty() && this->mix_DMk[ik].mixing_data.data != nullptr)
            {
                ModuleBase::WARNING_QUIT("Mix_DMk_2D",
                                         "an allocated DMK history cannot be reset to an empty stream");
            }
        }
    }
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::restart_all(const std::vector<std::vector<Tdata>>& data_in)
{
    ModuleBase::TITLE("Mix_DMk_2D", "restart_all");
    assert(this->mix_DMk.size() == data_in.size());
    assert(this->mixing != nullptr);
    for (int ik = 0; ik < data_in.size(); ++ik)
    {
        this->mix_DMk[ik].data_out = data_in[ik];
        this->mixing->init_mixing_data(this->mix_DMk[ik].mixing_data, data_in[ik].size(), sizeof(Tdata));
    }
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::mix_all(const std::vector<std::vector<Tdata>>& data_in)
{
    ModuleBase::TITLE("Mix_DMk_2D", "mix_all");
    assert(this->mix_DMk.size() == data_in.size());
    assert(this->mixing != nullptr);
    for (int ik = 0; ik < data_in.size(); ++ik)
    {
        // A rank may own no local matrix elements. Its empty history carries
        // no numerical data and does not need a push.
        if (data_in[ik].empty())
            { continue; }
        this->mixing->push_data(this->mix_DMk[ik].mixing_data, this->mix_DMk[ik].data_out.data(), data_in[ik].data(), nullptr, false);
        this->mixing->mix_data(this->mix_DMk[ik].mixing_data, this->mix_DMk[ik].data_out.data());
    }
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::mix_all(const std::vector<std::vector<Tdata>>& data_in,
                                const std::size_t beta_split_index)
{
    ModuleBase::TITLE("Mix_DMk_2D", "mix_all_split");
    assert(this->mix_DMk.size() == data_in.size());
    assert(this->mixing != nullptr);
    for (std::size_t ik = 0; ik < data_in.size(); ++ik)
    {
        if (data_in[ik].empty())
            { continue; }

        const std::size_t length = data_in[ik].size();
        this->mixing->push_data(
            this->mix_DMk[ik].mixing_data,
            this->mix_DMk[ik].data_out.data(),
            data_in[ik].data(),
            nullptr,
            [this, beta_split_index, length](Tdata* out, const Tdata* in, const Tdata* residual) {
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 4096 / sizeof(Tdata))
                #endif
                for (std::size_t i = 0; i < length; ++i)
                {
                    const double beta = i < beta_split_index ? this->mixing->mixing_beta
                                                             : PARAM.inp.mixing_beta_mag;
                    out[i] = in[i] + beta * residual[i];
                }
            },
            false);
        this->mixing->mix_data(this->mix_DMk[ik].mixing_data, this->mix_DMk[ik].data_out.data());
    }
}

template class Mix_DMk_2D<double>;
template class Mix_DMk_2D<std::complex<double>>;
