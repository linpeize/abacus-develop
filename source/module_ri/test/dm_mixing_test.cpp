#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/module_mixing/broyden_mixing.h"
#include "module_base/parallel_2d.h"
#include "module_ri/DmkSpinTransform.h"
#include "module_ri/Mix_DMk_2D.h"
#define private public
#include "module_parameter/parameter.h"
#undef private

#include <cstddef>

/************************************************
 *  unit test of charge_mixing.cpp & Mix_DMk_2D.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Mix_DMk_2D::mix:
 *      mix the density matrix data according to the set mixing mode.
 *
 */

class DM_Mixing_Test : public ::testing::Test
{
  public:
    DM_Mixing_Test()
    {
        mixing = new Base_Mixing::Broyden_Mixing(ndim, mixing_beta);
        mix_data_vector = std::vector<std::vector<double>>(2);
        mix_complexdata_vector = std::vector<std::vector<std::complex<double>>>(3);
        mix_data_vector[0].resize(nr * nc);
        mix_data_vector[1].resize(nr * nc);
        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j < nc; ++j)
            {
                mix_data_vector[0][i * nc + j] = i * nc + j;
                mix_data_vector[1][i * nc + j] = i * nc + j + 0.2;
            }
        }
        mix_complexdata_vector[0].resize(nr * nc);
        mix_complexdata_vector[1].resize(nr * nc);
        mix_complexdata_vector[2].resize(nr * nc);
        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j < nc; ++j)
            {
                mix_complexdata_vector[0][i * nc + j] = std::complex<double>{ double(i), double(j) };
                mix_complexdata_vector[1][i * nc + j] = std::complex<double>{ double(i), double(j) + 0.2 };
                mix_complexdata_vector[2][i * nc + j] = std::complex<double>{ double(i) + 0.8, double(j) };
            }
        }
    };
    ~DM_Mixing_Test()
    {
        delete mixing;
    };
    Base_Mixing::Mixing* mixing = nullptr;
    const int nr = 2;
    const int nc = 2;
    const int ndim = 1;
    const double mixing_beta = 0.3;

  protected:
    std::vector<std::vector<double>> mix_data_vector;
    std::vector<std::vector<std::complex<double>>> mix_complexdata_vector;
};

namespace
{

using Complex = std::complex<double>;

template <typename Tdata>
std::vector<const std::vector<Tdata>*> make_pointer_view(const std::vector<std::vector<Tdata>>& data)
{
    std::vector<const std::vector<Tdata>*> view(data.size());
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        view[i] = &data[i];
    }
    return view;
}

void expect_complex_vector_near(const std::vector<Complex>& actual,
                                const std::vector<Complex>& expected,
                                const double tolerance = 1e-14)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        EXPECT_NEAR(actual[i].real(), expected[i].real(), tolerance);
        EXPECT_NEAR(actual[i].imag(), expected[i].imag(), tolerance);
    }
}

class ScopedMixingBetaMag
{
  public:
    explicit ScopedMixingBetaMag(const double mixing_beta_mag)
        : saved_mixing_beta_mag_(PARAM.inp.mixing_beta_mag)
    {
        PARAM.input.mixing_beta_mag = mixing_beta_mag;
    }

    ~ScopedMixingBetaMag()
    {
        PARAM.input.mixing_beta_mag = this->saved_mixing_beta_mag_;
    }

  private:
    const double saved_mixing_beta_mag_;
};

class SplitSpinParallel : public Parallel_2D
{
  public:
    void split_row_pair()
    {
        set_serial(2, 2);
        // Mimic a block-cyclic layout in which only one spin partner is local.
        global2local_row_[1] = -1;
    }

    void split_col_pair()
    {
        set_serial(2, 2);
        global2local_col_[1] = -1;
    }

    void set_empty_serial()
    {
        is_serial = true;
        nrow = 0;
        ncol = 0;
        nloc = 0;
        global2local_row_.clear();
        global2local_col_.clear();
        local2global_row_.clear();
        local2global_col_.clear();
    }
};

} // namespace

TEST_F(DM_Mixing_Test, Mix_DMk_2D)
{
    //Gamma only
    Mix_DMk_2D<double> mix_dmk_gamma;
    mix_dmk_gamma.set_nks(1);
    mix_dmk_gamma.set_mixing_plain(1.0);
    std::vector<std::vector<std::vector<double>>> dm_gamma(2);
    dm_gamma[0] = std::vector<std::vector<double>>(1);
    dm_gamma[0][0] = mix_data_vector[0];
    dm_gamma[1] = std::vector<std::vector<double>>(1);
    dm_gamma[1][0] = mix_data_vector[1];
    for (int istep = 0; istep < 2; ++istep)
    {
        mix_dmk_gamma.mix(dm_gamma[istep], (istep == 0));
    }
    std::vector<const std::vector<double>*> dm_gamma_out = mix_dmk_gamma.get_DMk_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            EXPECT_DOUBLE_EQ(dm_gamma_out[0][0][i * nc + j], mix_data_vector[1][i * nc + j]);
        }
    }

    // not Gamma only
    Mix_DMk_2D<std::complex<double>> mix_dmk;
    mix_dmk.set_nks(1);
    mix_dmk.set_mixing_plain(1.0);
    std::vector<std::vector<std::vector<std::complex<double>>>> dm(2);
    dm[0] = std::vector<std::vector<std::complex<double>>>(1);
    dm[0][0] = mix_complexdata_vector[0];
    dm[1] = std::vector<std::vector<std::complex<double>>>(1);
    dm[1][0] = mix_complexdata_vector[1];
    for (int istep = 0; istep < 2; ++istep)
    {
        mix_dmk.mix(dm[istep], (istep == 0));
    }
    std::vector<const std::vector<std::complex<double>>*> dm_out = mix_dmk.get_DMk_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            EXPECT_DOUBLE_EQ(dm_out[0][0][i * nc + j].real(), mix_complexdata_vector[1][i * nc + j].real());
            EXPECT_DOUBLE_EQ(dm_out[0][0][i * nc + j].imag(), mix_complexdata_vector[1][i * nc + j].imag());
        }
    }

    // Shared Broyden mix
    Mix_DMk_2D<std::complex<double>> mix_dmk_broyden;
    mix_dmk_broyden.set_nks(1);
    mix_dmk_broyden.set_mixing(mixing);
    mixing->coef = { 1.1, -0.1 };
    std::vector<std::vector<std::vector<std::complex<double>>>> dm_broyden(3);
    for (int istep = 0; istep < 3; ++istep)
    {
        dm_broyden[istep] = std::vector<std::vector<std::complex<double>>>(1);
        dm_broyden[istep][0] = mix_complexdata_vector[istep];
        mix_dmk_broyden.mix(dm_broyden[istep], (istep == 0));
    }
    std::vector<const std::vector<std::complex<double>>*> dm_broyden_out = mix_dmk_broyden.get_DMk_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            std::complex<double> first_step_result
                = (1 - mixing_beta) * mix_complexdata_vector[0][i * nc + j]
                  + mixing_beta * mix_complexdata_vector[1][i * nc + j];
            std::complex<double> second_step_result
                = (1 - mixing_beta) * first_step_result + mixing_beta * mix_complexdata_vector[2][i * nc + j];
            std::complex<double> ref = second_step_result * mixing->coef[1] + first_step_result * mixing->coef[0];
            EXPECT_DOUBLE_EQ(dm_broyden_out[0][0][i * nc + j].real(), ref.real());
            EXPECT_DOUBLE_EQ(dm_broyden_out[0][0][i * nc + j].imag(), ref.imag());
        }
    }
}

TEST(DmkSpinTransformTest, Nspin2PairsSpinMajorStreamsAndRoundTrips)
{
    // Physical streams are [up(k0), up(k1), down(k0), down(k1)].
    const std::vector<std::vector<double>> physical = {
        {1.0, 2.0}, {10.0, 20.0}, {3.0, 4.0}, {30.0, 40.0}};

    const std::vector<std::vector<double>> pauli = DmkSpinTransform::to_pauli_nspin2(physical);

    ASSERT_EQ(pauli.size(), 2);
    EXPECT_EQ(pauli[0], (std::vector<double>{4.0, 6.0, -2.0, -2.0}));
    EXPECT_EQ(pauli[1], (std::vector<double>{40.0, 60.0, -20.0, -20.0}));
    EXPECT_EQ(DmkSpinTransform::to_physical_nspin2(make_pointer_view(pauli)), physical);
}

TEST(DmkSpinTransformTest, Nspin2RejectsInvalidShapesAndNullPointers)
{
    const std::vector<std::vector<double>> empty_physical;
    EXPECT_EXIT(DmkSpinTransform::to_pauli_nspin2(empty_physical), ::testing::ExitedWithCode(1), "");

    const std::vector<std::vector<double>> odd_outer = {{1.0}};
    EXPECT_EXIT(DmkSpinTransform::to_pauli_nspin2(odd_outer), ::testing::ExitedWithCode(1), "");

    const std::vector<std::vector<double>> unequal_lengths = {{1.0}, {2.0, 3.0}};
    EXPECT_EXIT(DmkSpinTransform::to_pauli_nspin2(unequal_lengths), ::testing::ExitedWithCode(1), "");

    const std::vector<double> odd_packed = {1.0, 2.0, 3.0};
    const std::vector<const std::vector<double>*> odd_view = {&odd_packed};
    EXPECT_EXIT(DmkSpinTransform::to_physical_nspin2(odd_view), ::testing::ExitedWithCode(1), "");

    const std::vector<const std::vector<double>*> null_view = {nullptr};
    EXPECT_EXIT(DmkSpinTransform::to_physical_nspin2(null_view), ::testing::ExitedWithCode(1), "");

    const std::vector<const std::vector<double>*> empty_view;
    EXPECT_EXIT(DmkSpinTransform::to_physical_nspin2(empty_view), ::testing::ExitedWithCode(1), "");

    const std::vector<double> packed_two = {1.0, 2.0};
    const std::vector<double> packed_four = {1.0, 2.0, 3.0, 4.0};
    const std::vector<const std::vector<double>*> unequal_view = {&packed_two, &packed_four};
    EXPECT_EXIT(DmkSpinTransform::to_physical_nspin2(unequal_view), ::testing::ExitedWithCode(1), "");
}

TEST(MixDmk2DTest, PlainMixWithoutSplitUsesOneBetaForWholeStream)
{
    const ScopedMixingBetaMag magnetic_beta(0.75);
    Mix_DMk_2D<double> mixer;
    mixer.set_nks(1);
    mixer.set_mixing_plain(0.25);

    mixer.mix({{0.0, 0.0, 0.0, 0.0}}, true);
    mixer.mix({{4.0, 8.0, 4.0, 8.0}}, false);

    const std::vector<const std::vector<double>*> output = mixer.get_DMk_out();
    ASSERT_EQ(output.size(), 1);
    EXPECT_EQ(*output[0], (std::vector<double>{1.0, 2.0, 1.0, 2.0}));
}

TEST(MixDmk2DTest, SplitBetaUsesEngineBetaForPrefixAndInputMagneticBetaForSuffix)
{
    const ScopedMixingBetaMag magnetic_beta(0.75);
    Mix_DMk_2D<double> mixer;
    mixer.set_nks(1);
    mixer.set_mixing_plain(0.25);

    mixer.mix({{0.0, 0.0, 0.0, 0.0}}, true, 2);
    mixer.mix({{4.0, 8.0, 4.0, 8.0}}, false, 2);

    const std::vector<const std::vector<double>*> output = mixer.get_DMk_out();
    ASSERT_EQ(output.size(), 1);
    EXPECT_EQ(*output[0], (std::vector<double>{1.0, 2.0, 3.0, 6.0}));
}

TEST_F(DM_Mixing_Test, BorrowedChargeEngineUsesInputMagneticBeta)
{
    const ScopedMixingBetaMag magnetic_beta(0.8);

    Mix_DMk_2D<double> mixer;
    mixer.set_nks(1);
    mixer.set_mixing(mixing);
    mixer.mix({{0.0, 0.0, 0.0, 0.0}}, true, 2);
    mixer.mix({{1.0, 1.0, 1.0, 1.0}}, false, 2);

    EXPECT_EQ(*mixer.get_DMk_out()[0], (std::vector<double>{mixing_beta, mixing_beta, 0.8, 0.8}));
}

TEST(MixDmk2DTest, EmptyLocalStreamsAreValidFromInitialSeed)
{
    Mix_DMk_2D<double> standard;
    standard.set_nks(1);
    standard.set_mixing_plain(0.5);
    standard.mix({{}}, true);
    standard.mix({{}}, false);
    ASSERT_EQ(standard.get_DMk_out().size(), 1);
    EXPECT_TRUE(standard.get_DMk_out()[0]->empty());

    Mix_DMk_2D<double> split;
    split.set_nks(1);
    split.set_mixing_plain(0.5);
    split.mix({{}}, true, 0);
    split.mix({{}}, false, 0);
    ASSERT_EQ(split.get_DMk_out().size(), 1);
    EXPECT_TRUE(split.get_DMk_out()[0]->empty());
}

TEST(MixDmk2DTest, RejectsResetFromAllocatedHistoryToEmptyStream)
{
    Mix_DMk_2D<double> mixer;
    mixer.set_nks(1);
    mixer.set_mixing_plain(0.5);
    mixer.mix({{1.0}}, true);
    EXPECT_EXIT(mixer.mix({{}}, true), ::testing::ExitedWithCode(1), "");
}

TEST(MixDmk2DTest, RestartReseedsHistoryAndAlwaysChecksSplitBoundary)
{
    const ScopedMixingBetaMag magnetic_beta(0.25);
    Mix_DMk_2D<double> mixer;
    mixer.set_nks(1);
    mixer.set_mixing_plain(0.5);

    const std::vector<std::vector<double>> seed = {{0.0, 0.0}};
    EXPECT_EXIT(mixer.mix(seed, true, 3), ::testing::ExitedWithCode(1), "");

    mixer.mix(seed, true, 1);
    mixer.mix({{2.0, 4.0}}, false, 1);
    EXPECT_EQ(*mixer.get_DMk_out()[0], (std::vector<double>{1.0, 1.0}));

    // A restart is a reset/seed of this SCF's DM history, not a new engine.
    mixer.mix({{10.0, 20.0}}, true, 1);
    EXPECT_EQ(*mixer.get_DMk_out()[0], (std::vector<double>{10.0, 20.0}));
    mixer.mix({{14.0, 28.0}}, false, 1);
    EXPECT_EQ(*mixer.get_DMk_out()[0], (std::vector<double>{12.0, 22.0}));
}

TEST(MixDmk2DTest, RejectsUnsafeEngineCountAndHistoryShapes)
{
    const std::vector<std::vector<double>> one_stream = {{1.0, 2.0}};
    Mix_DMk_2D<double> no_engine;
    no_engine.set_nks(1);
    EXPECT_EXIT(no_engine.mix(one_stream, true), ::testing::ExitedWithCode(1), "");

    Mix_DMk_2D<double> wrong_count;
    wrong_count.set_nks(2);
    wrong_count.set_mixing_plain(1.0);
    EXPECT_EXIT(wrong_count.mix(one_stream, true), ::testing::ExitedWithCode(1), "");

    Mix_DMk_2D<double> changed_shape;
    changed_shape.set_nks(1);
    changed_shape.set_mixing_plain(1.0);
    changed_shape.mix(one_stream, true);
    const std::vector<std::vector<double>> longer_stream = {{1.0, 2.0, 3.0}};
    EXPECT_EXIT(changed_shape.mix(longer_stream, false), ::testing::ExitedWithCode(1), "");
}

TEST(DmkMixingPipelineTest, Nspin2MixesPauliAndReturnsPhysicalStreams)
{
    const ScopedMixingBetaMag magnetic_beta(0.5);
    const std::vector<std::vector<double>> initial_physical = {{1.0}, {1.0}};
    const std::vector<std::vector<double>> next_physical = {{5.0}, {3.0}};

    Mix_DMk_2D<double> mixer;
    mixer.set_nks(1);
    mixer.set_mixing_plain(0.25);
    const std::vector<std::vector<double>> initial_pauli
        = DmkSpinTransform::to_pauli_nspin2(initial_physical);
    const std::vector<std::vector<double>> next_pauli
        = DmkSpinTransform::to_pauli_nspin2(next_physical);
    mixer.mix(initial_pauli, true, 1);
    mixer.mix(next_pauli, false, 1);

    // get_DMk_out remains a Pauli view until the explicit inverse transform.
    ASSERT_EQ(*mixer.get_DMk_out()[0], (std::vector<double>{3.5, 1.0}));
    const std::vector<std::vector<double>> mixed_physical
        = DmkSpinTransform::to_physical_nspin2(mixer.get_DMk_out());
    EXPECT_EQ(mixed_physical, (std::vector<std::vector<double>>{{2.25}, {1.25}}));
}

TEST(DmkSpinTransformTest, Nspin4SupportsRowAndColumnMajorLayouts)
{
    Parallel_2D pv;
    pv.set_serial(2, 2); // nb=1 is valid because both spin partners are local.

    const Complex uu(1.0, 2.0);
    const Complex ud(3.0, 4.0);
    const Complex du(5.0, 6.0);
    const Complex dd(7.0, 8.0);
    const std::vector<Complex> expected_pauli = {
        Complex(8.0, 10.0), Complex(8.0, 10.0), Complex(2.0, -2.0), Complex(-6.0, -6.0)};

    const std::vector<std::vector<Complex>> row_major = {{uu, ud, du, dd}};
    const std::vector<std::vector<Complex>> row_pauli
        = DmkSpinTransform::to_pauli_nspin4(row_major, pv, false);
    expect_complex_vector_near(row_pauli[0], expected_pauli);
    const std::vector<std::vector<Complex>> row_roundtrip
        = DmkSpinTransform::to_physical_nspin4(make_pointer_view(row_pauli), pv, false);
    expect_complex_vector_near(row_roundtrip[0], row_major[0]);

    // Column-major physical storage swaps the ud and du flat positions.
    const std::vector<std::vector<Complex>> column_major = {{uu, du, ud, dd}};
    const std::vector<std::vector<Complex>> column_pauli
        = DmkSpinTransform::to_pauli_nspin4(column_major, pv, true);
    expect_complex_vector_near(column_pauli[0], expected_pauli);
    const std::vector<std::vector<Complex>> column_roundtrip
        = DmkSpinTransform::to_physical_nspin4(make_pointer_view(column_pauli), pv, true);
    expect_complex_vector_near(column_roundtrip[0], column_major[0]);
}

TEST(DmkSpinTransformTest, Nspin4RejectsInvalidAndNonColocatedLayouts)
{
    Parallel_2D odd_shape;
    odd_shape.set_serial(3, 2);
    const std::vector<std::vector<Complex>> odd_physical(1, std::vector<Complex>(6));
    EXPECT_EXIT(DmkSpinTransform::to_pauli_nspin4(odd_physical, odd_shape, false),
                ::testing::ExitedWithCode(1),
                "");

    Parallel_2D even_shape;
    even_shape.set_serial(2, 2);
    const std::vector<std::vector<Complex>> wrong_length(1, std::vector<Complex>(3));
    EXPECT_EXIT(DmkSpinTransform::to_pauli_nspin4(wrong_length, even_shape, false),
                ::testing::ExitedWithCode(1),
                "");

    const std::vector<const std::vector<Complex>*> null_view = {nullptr};
    EXPECT_EXIT(DmkSpinTransform::to_physical_nspin4(null_view, even_shape, false),
                ::testing::ExitedWithCode(1),
                "");

    const std::vector<std::vector<Complex>> empty_physical;
    EXPECT_EXIT(DmkSpinTransform::to_pauli_nspin4(empty_physical, even_shape, false),
                ::testing::ExitedWithCode(1),
                "");

    const std::vector<const std::vector<Complex>*> empty_view;
    EXPECT_EXIT(DmkSpinTransform::to_physical_nspin4(empty_view, even_shape, false),
                ::testing::ExitedWithCode(1),
                "");

    const std::vector<Complex> wrong_inverse_length(3);
    const std::vector<const std::vector<Complex>*> wrong_inverse_view = {&wrong_inverse_length};
    EXPECT_EXIT(DmkSpinTransform::to_physical_nspin4(wrong_inverse_view, even_shape, false),
                ::testing::ExitedWithCode(1),
                "");

    SplitSpinParallel split;
    split.split_row_pair();
    const std::vector<std::vector<Complex>> local_physical(1, std::vector<Complex>(4));
    EXPECT_EXIT(DmkSpinTransform::to_pauli_nspin4(local_physical, split, false),
                ::testing::ExitedWithCode(1),
                "");

    SplitSpinParallel col_split;
    col_split.split_col_pair();
    EXPECT_EXIT(DmkSpinTransform::to_pauli_nspin4(local_physical, col_split, false),
                ::testing::ExitedWithCode(1),
                "");
}

TEST(DmkSpinTransformTest, ZeroLocalSizeRoundTripsAreValid)
{
    const std::vector<std::vector<double>> nspin2_physical(2);
    const std::vector<std::vector<double>> nspin2_pauli
        = DmkSpinTransform::to_pauli_nspin2(nspin2_physical);
    ASSERT_EQ(nspin2_pauli.size(), 1);
    EXPECT_TRUE(nspin2_pauli[0].empty());
    EXPECT_EQ(DmkSpinTransform::to_physical_nspin2(make_pointer_view(nspin2_pauli)), nspin2_physical);

    SplitSpinParallel empty_pv;
    empty_pv.set_empty_serial();
    const std::vector<std::vector<Complex>> nspin4_physical(1);
    const std::vector<std::vector<Complex>> nspin4_pauli
        = DmkSpinTransform::to_pauli_nspin4(nspin4_physical, empty_pv, false);
    ASSERT_EQ(nspin4_pauli.size(), 1);
    EXPECT_TRUE(nspin4_pauli[0].empty());
    EXPECT_EQ(DmkSpinTransform::to_physical_nspin4(make_pointer_view(nspin4_pauli), empty_pv, false),
              nspin4_physical);
}

TEST(DmkMixingPipelineTest, Nspin4UsesChargeAndMagneticBetasBeforeInverse)
{
    const ScopedMixingBetaMag magnetic_beta(0.5);
    Parallel_2D pv;
    pv.set_serial(2, 2);
    const std::vector<std::vector<Complex>> initial_physical(1, std::vector<Complex>(4));
    const std::vector<std::vector<Complex>> next_physical = {
        {Complex(4.0, 0.0), Complex(2.0, 0.0), Complex(6.0, 0.0), Complex(0.0, 0.0)}};

    Mix_DMk_2D<Complex> mixer;
    mixer.set_nks(1);
    mixer.set_mixing_plain(0.25);
    const std::vector<std::vector<Complex>> initial_pauli
        = DmkSpinTransform::to_pauli_nspin4(initial_physical, pv, false);
    const std::vector<std::vector<Complex>> next_pauli
        = DmkSpinTransform::to_pauli_nspin4(next_physical, pv, false);
    mixer.mix(initial_pauli, true, 1);
    mixer.mix(next_pauli, false, 1);

    expect_complex_vector_near(*mixer.get_DMk_out()[0],
                               {Complex(1.0, 0.0),
                                Complex(4.0, 0.0),
                                Complex(0.0, -2.0),
                                Complex(2.0, 0.0)});
    const std::vector<std::vector<Complex>> mixed_physical
        = DmkSpinTransform::to_physical_nspin4(mixer.get_DMk_out(), pv, false);
    expect_complex_vector_near(mixed_physical[0],
                               {Complex(1.5, 0.0),
                                Complex(1.0, 0.0),
                                Complex(3.0, 0.0),
                                Complex(-0.5, 0.0)});
}
