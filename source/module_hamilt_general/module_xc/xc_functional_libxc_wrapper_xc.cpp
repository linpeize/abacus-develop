#ifdef USE_LIBXC

#include "xc_functional_libxc.h"

void XC_Functional_Libxc::xc_spin_libxc(
        const std::vector<int> &func_id,
        const double &rhoup, const double &rhodw,
		double &exc, double &vxcup, double &vxcdw)
{
    double e, vup, vdw;
    double *rho_ud, *vxc_ud;
    exc = vxcup = vxcdw = 0.0;

    rho_ud = new double[2];
    vxc_ud = new double[2];
    rho_ud[0] = rhoup;
    rho_ud[1] = rhodw;

    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(func_id, XC_POLARIZED);

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_LDA)
        {
            // call Libxc function: xc_lda_exc_vxc
            xc_lda_exc_vxc( &func, 1, rho_ud, &e, vxc_ud);
        }
        exc += e;
        vxcup += vxc_ud[0];
        vxcdw += vxc_ud[1];
    }    

    XC_Functional_Libxc::finish_func(funcs);
    delete[] rho_ud;
    delete[] vxc_ud;
}

#endif