#ifndef XC_FUNCTIONAL_LIBXC_H
#define XC_FUNCTIONAL_LIBXC_H

#ifdef USE_LIBXC

#include "module_base/matrix.h"
#include "module_base/vector3.h"

#include <xc.h>

#include <tuple>
#include <vector>

	class Charge;

namespace XC_Functional_Libxc
{
//-------------------
//  xc_functional_libxc.cpp
//-------------------	

	// sets functional type, which allows combination of LIBXC keyword connected by "+"
	//		for example, "XC_LDA_X+XC_LDA_C_PZ"
	std::pair<int,std::vector<int>> set_xc_type_libxc(const std::string xc_func_in);

	// converts func_id into corresponding xc_func_type vector
	std::vector<xc_func_type> init_func(const std::vector<int> &func_id, const int xc_polarized);

	void finish_func(std::vector<xc_func_type> &funcs);

//-------------------
//  xc_functional_libxc_vxc.cpp
//-------------------	

	std::tuple<double,double,ModuleBase::matrix> v_xc_libxc(
        const std::vector<int> &func_id,
		const int &nrxx, // number of real-space grid
		const double &omega, // volume of cell
		const double tpiba,
		const Charge* const chr); // charge density

	// for mGGA functional
	std::tuple<double,double,ModuleBase::matrix,ModuleBase::matrix> v_xc_meta(
        const std::vector<int> &func_id,
		const int &nrxx, // number of real-space grid
		const double &omega, // volume of cell
		const double tpiba,
		const Charge* const chr);

//-------------------
//  xc_functional_libxc_wrapper_xc.cpp
//-------------------

	void xc_spin_libxc(
        const std::vector<int> &func_id,
		const double &rhoup, const double &rhodw,
		double &exc, double &vxcup, double &vxcdw);


//-------------------
//  xc_functional_libxc_wrapper_gcxc.cpp
//-------------------

	// the entire GGA functional, for nspin=1 case
	void gcxc_libxc(
        const std::vector<int> &func_id,
		const double &rho, const double &grho,
		double &sxc, double &v1xc, double &v2xc);

	// the entire GGA functional, for nspin=2 case
	void gcxc_spin_libxc(
        const std::vector<int> &func_id,
		double rhoup, double rhodw, 
		ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
		double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud);


//-------------------
//  xc_functional_libxc_wrapper_tauxc.cpp
//-------------------

	// wrapper for the mGGA functionals
	void tau_xc(
        const std::vector<int> &func_id,
		const double &rho, const double &grho, const double &atau, double &sxc,
		double &v1xc, double &v2xc, double &v3xc);

	void tau_xc_spin(
        const std::vector<int> &func_id,
		double rhoup, double rhodw, 
		ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
		double tauup, double taudw,
		double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud,
		double &v3xcup, double &v3xcdw);
}

#endif

#endif