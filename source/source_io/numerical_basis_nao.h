#ifndef NUMERICAL_BASIS_NAO_H
#define NUMERICAL_BASIS_NAO_H

#include "source_lcao/module_hcontainer/hcontainer.h"

#include <array>
#include <vector>

	namespace ModuleESolver{ template<typename TK, typename TR> class ESolver_KS_LCAO; }
	class Parallel_Orbitals;
	namespace hamilt{ template<typename TK, typename TR> class HamiltLCAO; }
	class K_Vectors;
	namespace psi{ template<typename TK, typename Device> class Psi; }
	class UnitCell;

template<typename TK, typename TR>
class Numerical_Basis_Nao
{
  public:
	Numerical_Basis_Nao(
		const UnitCell &ucell_in,
		const K_Vectors &kv_in,
		const Parallel_Orbitals &pv_in,
		const psi::Psi<TK, base_device::DEVICE_CPU> &psi_in,
		const hamilt::HamiltLCAO<TK,TR> &hamilt_in,
		const ModuleESolver::ESolver_KS_LCAO<TK,TR> &es_in);

	void output_overlap() const;


  private:

	class Matrix_TK
	{
	  public:
		Matrix_TK(){}
		Matrix_TK(const std::size_t nr_in, const std::size_t nc_in)
			:nr(nr_in), nc(nc_in) { v.resize(nr*nc); }
		Matrix_TK& operator=(Matrix_TK &&m)
			{ nr=m.nr; nc=m.nc; v=std::move(m.v); }
		std::size_t nr=0;
		std::size_t nc=0;
		std::vector<TK> v;
		TK* data() { return v.data(); }
		const TK* data() const { return v.data(); }
		TK operator()(const std::size_t ir, const std::size_t ic) const { return v[ir*nc+ic]; }
	};

	const UnitCell &ucell;
	const K_Vectors &kv;
	const Parallel_Orbitals &pv;
	const psi::Psi<TK, base_device::DEVICE_CPU> &psi;
	const ModuleESolver::ESolver_KS_LCAO<TK,TR> &es;		// for calculate T(R)
	const hamilt::HamiltLCAO<TK,TR> &hamilt;				// for get H(R) and S(R)

	std::array<int,9> desc_nb_nb;
	std::array<int,9> desc_nb_nw;
	std::array<int,9> desc_nw_nw;

	std::array<int,9> init_desc(const int &N, const int &M) const;
	Matrix_TK init_matrix(const int desc[9]) const;

	static void pgemm(
		const char &transA, const char &transB,
		const int m, const int n, const int k,
		const TK*const A, const int*const descA,
		const TK*const B, const int*const descB,
		TK*const C, const int*const descC);

	std::vector<Matrix_TK> cal_Sk(const hamilt::HContainer<TR> &SR) const;
	std::vector<Matrix_TK> cal_Qk(const std::vector<Matrix_TK> &Sk) const;
	std::vector<Matrix_TK> cal_Vk(const std::vector<Matrix_TK> &Qk) const;
	hamilt::HContainer<TR> cal_Tr() const;

	std::vector<std::array<int,4>> get_mu_index() const;
	void output(
		const std::string &file_name,
		const std::vector<Matrix_TK> &Q,
		const std::vector<Matrix_TK> &S,
		const std::vector<Matrix_TK> &V) const;
	void output_info(std::ofstream &ofs, const std::size_t nwfc) const;
};

#endif