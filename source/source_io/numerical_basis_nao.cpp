#include "numerical_basis_nao.h"
#include "source_base/element_basis_index.h"
#include "source_esolver/esolver_ks_lcao.h"
#include "source_cell/unitcell.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_operator_lcao/ekinetic_new.h"
#include "source_cell/klist.h"
#include "source_psi/psi.h"
#include "source_base/tool_title.h"

#include <cassert>

template<typename TK, typename TR>
Numerical_Basis_Nao<TK,TR>::Numerical_Basis_Nao(
		const UnitCell &ucell_in,
		const K_Vectors &kv_in,
		const Parallel_Orbitals &pv_in,
		const psi::Psi<TK, base_device::DEVICE_CPU> &psi_in,
		const hamilt::HamiltLCAO<TK,TR> &hamilt_in,
		const ModuleESolver::ESolver_KS_LCAO<TK,TR> &es_in)
	:ucell(ucell_in),
	 kv(kv_in),
	 pv(pv_in),
	 psi(psi_in),
	 hamilt(hamilt_in),
	 es(es_in)
{
	this->desc_nb_nb = init_desc(PARAM.inp.nbands, PARAM.inp.nbands);
	this->desc_nb_nw = init_desc(PARAM.inp.nbands, PARAM.globalv.nlocal);
	this->desc_nw_nw = init_desc(PARAM.globalv.nlocal, PARAM.globalv.nlocal);
}

template<typename TK, typename TR>
void Numerical_Basis_Nao<TK,TR>::output_overlap() const
{
	ModuleBase::TITLE("Numerical_Basis_Nao", "output_overlap");
	assert(this->pv.desc[1] == this->pv.desc_wfc[1]);
	
	const std::vector<Matrix_TK> Sk = cal_Sk(*(this->hamilt.getSR()));
	const std::vector<Matrix_TK> Qk = cal_Qk(Sk);
	const std::vector<Matrix_TK> Vk = cal_Vk(Qk);
	output("orb_matrix.0.dat", Qk, Sk, Vk);

	const hamilt::HContainer<TR> Tr = cal_Tr();
	const std::vector<Matrix_TK> Sk_T = cal_Sk(Tr);
	const std::vector<Matrix_TK> Qk_T = cal_Qk(Sk_T);
	const std::vector<Matrix_TK> Vk_T = cal_Vk(Qk_T);
	output("orb_matrix.1.dat", Qk_T, Sk_T, Vk_T);

	const std::vector<Matrix_TK> Sk_H = cal_Sk(*(this->hamilt.getHR()));
	const std::vector<Matrix_TK> Qk_H = cal_Qk(Sk_H);
	const std::vector<Matrix_TK> Vk_H = cal_Vk(Qk_H);
	output("orb_matrix.2.dat", Qk_H, Sk_H, Vk_H);
}


template<typename TK, typename TR>
std::array<int,9>
Numerical_Basis_Nao<TK,TR>::init_desc(const int &N, const int &M) const
{
	int nprow, npcol, myprow, mypcol;
	Cblacs_gridinfo(this->pv.desc_wfc[1], &nprow, &npcol, &myprow, &mypcol);
	const int mlocal = numroc_(&M, &this->pv.desc_wfc[4], &myprow, &this->pv.desc_wfc[6], &nprow);

	std::array<int,9> desc;
	int info;
	descinit_( desc.data(), &M, &N, &this->pv.desc_wfc[4], &this->pv.desc_wfc[5], &this->pv.desc_wfc[6], &this->pv.desc_wfc[7], &this->pv.desc_wfc[1], &mlocal, &info );
	assert(info==0);
	return desc;
};

template<typename TK, typename TR>
auto Numerical_Basis_Nao<TK,TR>::init_matrix(const int desc[9]) const
	-> Numerical_Basis_Nao::Matrix_TK
{
	int nprow, npcol, myprow, mypcol;
	Cblacs_gridinfo(desc[1], &nprow, &npcol, &myprow, &mypcol);
	const int mlocal = numroc_(&desc[2], &desc[4], &myprow, &desc[6], &nprow);
	const int nlocal = numroc_(&desc[3], &desc[5], &mypcol, &desc[7], &npcol);
	return Matrix_TK(mlocal,nlocal);
};


template<typename TK, typename TR>
void Numerical_Basis_Nao<TK,TR>::pgemm(
	const char &transA, const char &transB,
	const int m, const int n, const int k,
	const TK*const A, const int*const descA,
	const TK*const B, const int*const descB,
	TK*const C, const int*const descC)
{
	ScalapackConnector::gemm(
		transB, transA,
		n, m, k,
		TK(1.0),
		B, 1, 1, descB,
		A, 1, 1, descA,
		TK(0.0),
		C, 1, 1, descC);
};


template<typename TK, typename TR>
auto Numerical_Basis_Nao<TK,TR>::cal_Sk(const hamilt::HContainer<TR> &SR) const
	-> std::vector<Numerical_Basis_Nao<TK,TR>::Matrix_TK>
{
	ModuleBase::TITLE("Numerical_Basis_Nao", "cal_Sk");
	std::vector<Matrix_TK> Sk(this->kv.get_nkstot());
	const int major = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver)
		? this->pv.get_row_size() : this->pv.get_col_size();
	const int hk_type = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver)
		? 1 : 0;
	for (int ik=0; ik<this->kv.get_nkstot(); ++ik)
	{
		Sk[ik] = init_matrix(this->desc_nw_nw.data());
		hamilt::folding_HR(SR, Sk[ik].data(), this->kv.kvec_d[ik], major, hk_type);
	}
	return Sk;
}

template<typename TK, typename TR>
auto Numerical_Basis_Nao<TK,TR>::cal_Qk(const std::vector<Matrix_TK> &Sk) const
	-> std::vector<Numerical_Basis_Nao<TK,TR>::Matrix_TK>
{
	ModuleBase::TITLE("Numerical_Basis_Nao", "cal_Qk");
	std::vector<Matrix_TK> Qk(this->kv.get_nkstot());
	for (int ik=0; ik<this->kv.get_nkstot(); ++ik)
	{
		// C_S(ib,iwt2) = C(ib,iwt1) * S(iwt1,iwt2)
		Qk[ik] = init_matrix(this->desc_nb_nw.data());
		pgemm(
			'N', 'N',
			PARAM.inp.nbands, PARAM.globalv.nlocal, PARAM.globalv.nlocal,
			this->psi.get_pointer(ik), this->pv.desc_wfc,
			Sk[ik].data(), this->pv.desc,
			Qk[ik].data(), this->desc_nb_nw.data());
	}
	return Qk;
}

template<typename TK, typename TR>
auto Numerical_Basis_Nao<TK,TR>::cal_Vk(const std::vector<Matrix_TK> &Qk) const
	-> std::vector<Numerical_Basis_Nao<TK,TR>::Matrix_TK>
{
	ModuleBase::TITLE("Numerical_Basis_Nao", "cal_Vk");
	std::vector<Matrix_TK> Vk(this->kv.get_nkstot());
	for (int ik=0; ik<this->kv.get_nkstot(); ++ik)
	{
		// C_S_C(ib1,ib2) = C_S(ib1,iwt) * C(ib2,iwt)
		Vk[ik] = init_matrix(this->desc_nb_nb.data());
		pgemm(
			'N', 'T',
			PARAM.inp.nbands, PARAM.inp.nbands, PARAM.globalv.nlocal,
			Qk[ik].data(), this->desc_nb_nw.data(),
			this->psi.get_pointer(ik), this->pv.desc_wfc,
			Vk[ik].data(), this->desc_nb_nb.data());
	}
	return Vk;
}

template<typename TK, typename TR>
hamilt::HContainer<TR> Numerical_Basis_Nao<TK,TR>::cal_Tr() const
{
	ModuleBase::TITLE("Numerical_Basis_Nao", "cal_Tr");
	hamilt::HContainer<TR> Tr(this->ucell, &this->pv);
	hamilt::HS_Matrix_K<TK> Tk_tmp(&this->pv);
	hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>> ekinetic(
		&Tk_tmp,
		this->kv.kvec_d,
		&Tr,
		&this->ucell,
		this->es.get_orb().cutoffs(),
		&this->es.get_gd(),
		this->es.get_two_center_bundle().kinetic_orb.get());
	ekinetic.contributeHR();
	return Tr;	
}

template<typename TK, typename TR>
std::vector<std::array<int,4>> Numerical_Basis_Nao<TK,TR>::get_mu_index() const
{
	std::vector<std::array<int,4>> mu_index;
	for (int it=0; it<this->ucell.ntype; ++it) {
		for (int ia=0; ia<this->ucell.atoms[it].na; ++ia) {
			for (int il=0; il<this->ucell.atoms[it].nwl+1; ++il) {
					for (int im=0; im<2*il+1; ++im) {
						mu_index.push_back({it, ia, il, im});
	}}}}
	return mu_index;
}


template<typename TK, typename TR>
void Numerical_Basis_Nao<TK,TR>::output(
	const std::string &file_name,
	const std::vector<Matrix_TK> &Q,
	const std::vector<Matrix_TK> &S,
	const std::vector<Matrix_TK> &V) const
{
	ModuleBase::TITLE("Numerical_Basis_Nao", "output");

	assert(GlobalV::NPROC==1);				// parallel output unsupported yet

	std::ofstream ofs(file_name);

	const std::vector<std::array<int,4>> mu_index = get_mu_index();

	const ModuleBase::Element_Basis_Index::Range iw_range = ModuleBase::Element_Basis_Index::construct_range( this->es.get_orb() );
	const ModuleBase::Element_Basis_Index::IndexLNM iw_index = ModuleBase::Element_Basis_Index::construct_index( iw_range );

	output_info(ofs, mu_index.size());
	ofs << std::scientific << std::setprecision(15);

	ofs<<"<OVERLAP_Q>"<<std::endl;
	for (int ik=0; ik<this->kv.get_nkstot(); ++ik) {
		for (int ib=0; ib<PARAM.inp.nbands; ++ib) {
			for (const std::array<int,4> &mu : mu_index) {
				for (int in=0; in<this->ucell.atoms[mu[0]].l_nchi[mu[2]]; ++in)
				{
					const int iwt = ucell.itiaiw2iwt(mu[0], mu[1], int(iw_index[mu[0]][mu[2]][in][mu[3]]));
					ofs << std::real(Q[ik](ib,iwt)) <<"\t"<< std::imag(Q[ik](ib,iwt)) <<std::endl;
				}}}}
	ofs<<"</OVERLAP_Q>"<<std::endl<<std::endl;

	ofs<<"<OVERLAP_Sq>"<<std::endl;
	for (int ik=0; ik<this->kv.get_nkstot(); ++ik) {
		for (const std::array<int,4> &mu1 : mu_index) {
			for (const std::array<int,4> &mu2 : mu_index) {
				for (int in1=0; in1<this->ucell.atoms[mu1[0]].l_nchi[mu1[2]]; ++in1) {
					for (int in2=0; in2<this->ucell.atoms[mu2[0]].l_nchi[mu2[2]]; ++in2)
					{
						const int iwt1 = ucell.itiaiw2iwt(mu1[0], mu1[1], int(iw_index[mu1[0]][mu1[2]][in1][mu1[3]]));
						const int iwt2 = ucell.itiaiw2iwt(mu2[0], mu2[1], int(iw_index[mu2[0]][mu2[2]][in2][mu2[3]]));
						ofs << std::real(S[ik](iwt1,iwt2)) <<"\t"<< std::imag(S[ik](iwt1,iwt2)) <<std::endl;
					}}}}}
	ofs<<"</OVERLAP_Sq>"<<std::endl<<std::endl;
					
	ofs<<"<OVERLAP_V>"<<std::endl;
	for (int ik=0; ik<this->kv.get_nkstot(); ++ik) {
		for (int ib=0; ib<PARAM.inp.nbands; ++ib)
		{
			ofs << V[ik](ib,ib) <<std::endl;
		}}
	//for (int ik=0; ik<this->kv.get_nkstot(); ++ik) {
	//	for (int ib1=0; ib1<PARAM.inp.nbands; ++ib1)
	//	{
	//		for (int ib2=0; ib2<PARAM.inp.nbands; ++ib2)
	//		{
	//			ofs << V[ik](ib1,ib2) <<"\t";
	//		}
	//		ofs<<std::endl;
	//	}}
	ofs<<"</OVERLAP_V>"<<std::endl<<std::endl;
}


template<typename TK, typename TR>
void Numerical_Basis_Nao<TK,TR>::output_info(std::ofstream &ofs, const std::size_t nwfc) const
{
	ModuleBase::TITLE("Numerical_Basis_Nao", "output_info");
	if (GlobalV::MY_RANK == 0)
	{
		const std::streamsize precision_old = ofs.precision(10);
		ofs << this->ucell.lat0 << std::endl;

		ofs << this->ucell.latvec.e11 << " " << this->ucell.latvec.e12 << " " << this->ucell.latvec.e13 << std::endl;
		ofs << this->ucell.latvec.e21 << " " << this->ucell.latvec.e22 << " " << this->ucell.latvec.e23 << std::endl;
		ofs << this->ucell.latvec.e31 << " " << this->ucell.latvec.e32 << " " << this->ucell.latvec.e33 << std::endl;

		ofs << this->ucell.ntype << " ntype" << std::endl;
		for (int it=0; it<this->ucell.ntype; ++it)
		{
			ofs << this->ucell.atoms[it].label << " label" << std::endl;
			ofs << this->ucell.atoms[it].na << " na" << std::endl;
			for (int ia=0; ia<this->ucell.atoms[it].na; ++ia)
			{
				ofs << this->ucell.atoms[it].tau[ia].x << " " << this->ucell.atoms[it].tau[ia].y << " " << this->ucell.atoms[it].tau[ia].z
					<< std::endl;
			}
		}
		
		ofs << PARAM.inp.ecutwfc << " ecutwfc" << std::endl;
		ofs << PARAM.inp.ecutwfc << " ecutwfc_jlq" << std::endl;							// meaningless
		ofs << this->es.get_orb().get_rcutmax_Phi() << " rcut_Jlq" << std::endl;			// error when Rcut of different elements and different L are different
		ofs << 0 << " smooth" << std::endl;													// meaningless
		ofs << 0 << " sigma" << std::endl;													// meaningless
		ofs << 0 << " tolerence" << std::endl;												// meaningless

		ofs << this->es.get_orb().get_lmax() << " lmax" << std::endl;

		// NOTICE: ofs_warning << "\n The precison may affect the optimize result.";
		ofs << this->kv.get_nkstot() << " nks" << std::endl;
		ofs << PARAM.inp.nbands << " nbands" << std::endl;
		ofs << nwfc << " nwfc" << std::endl;
		ofs << this->es.get_orb().get_nchimax() << " ne " << std::endl;						// error when nchi of different elements and different L are different

		ofs << "<WEIGHT_OF_KPOINTS>" << std::endl;
		for (int ik=0; ik<this->kv.get_nkstot(); ++ik)
			{ ofs << this->kv.kvec_c[ik].x << " " << this->kv.kvec_c[ik].y << " " << this->kv.kvec_c[ik].z << " " << this->kv.wk[ik] * 0.5 << std::endl; }
        ofs << "</WEIGHT_OF_KPOINTS>" << std::endl<<std::endl;

		ofs.precision(precision_old);
	}
}




template class Numerical_Basis_Nao<double, double>;
template class Numerical_Basis_Nao<std::complex<double>, double>;
template class Numerical_Basis_Nao<std::complex<double>, std::complex<double>>;