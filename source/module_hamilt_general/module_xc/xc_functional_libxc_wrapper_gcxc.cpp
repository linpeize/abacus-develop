#ifdef USE_LIBXC

#include "xc_functional_libxc.h"

#include <xc.h>

void XC_Functional_Libxc::gcxc_libxc(
        const std::vector<int> &func_id,
    const double &rho, const double &grho, double &sxc,
          double &v1xc, double &v2xc)
{
    const double small = 1.e-6;
    const double smallg = 1.e-10;
    double s,v1,v2;
    sxc = v1xc = v2xc = 0.0;

    if (rho <= small || grho < smallg)
    {
        return;
    }

    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(func_id, XC_UNPOLARIZED);
    

    for(xc_func_type &func : funcs)
    {
        xc_gga_exc_vxc(&func, 1, &rho, &grho, &s, &v1, &v2);
        
        sxc += s * rho;
        v2xc += v2 * 2.0;
        v1xc += v1;
    }
    XC_Functional_Libxc::finish_func(funcs);
} // end subroutine gcxc_libxc

void XC_Functional_Libxc::gcxc_spin_libxc(
        const std::vector<int> &func_id,
        double rhoup, double rhodw, 
        ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
        double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud)
{
	std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(func_id, XC_POLARIZED);
    double *rho, *grho, *v1xc, *v2xc, *sgn, s;
    sxc = v1xcup = v1xcdw = 0.0;
    v2xcup = v2xcdw = v2xcud = 0.0;
    rho = new double[2];
    grho= new double[3];
    v1xc= new double[2];
    v2xc= new double[3];
    sgn = new double[2];
    
    rho[0] = rhoup;
    rho[1] = rhodw;
    grho[0] = gdr1.norm2();
    grho[1] = gdr1 * gdr2;
    grho[2] = gdr2.norm2();

    const double rho_threshold = 1E-6;
    const double grho_threshold = 1E-10;

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_GGA || func.info->family == XC_FAMILY_HYB_GGA)
        {
            sgn[0] = sgn[1] = 1.0;
            // call Libxc function: xc_gga_exc_vxc
            xc_gga_exc_vxc( &func, 1, rho, grho, &s, v1xc, v2xc);
            if(func.info->kind==XC_CORRELATION)
            {
                if ( rho[0]<rho_threshold || sqrt(std::abs(grho[0]))<grho_threshold )
                    sgn[0] = 0.0;
                if ( rho[1]<rho_threshold || sqrt(std::abs(grho[2]))<grho_threshold )
                    sgn[1] = 0.0;
            }
            sxc += s * (rho[0] * sgn[0] + rho[1] * sgn[1]);
            v1xcup += v1xc[0] * sgn[0];
            v1xcdw += v1xc[1] * sgn[1];
            v2xcup += 2.0 * v2xc[0] * sgn[0];
            v2xcud += v2xc[1] * sgn[0] * sgn[1];
            v2xcdw += 2.0 * v2xc[2] * sgn[1];
        }
    }
    delete[] grho;
    delete[] rho;
    delete[] v1xc;
    delete[] v2xc;
    delete[] sgn;
    XC_Functional_Libxc::finish_func(funcs);
}

#endif