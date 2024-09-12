#ifdef USE_LIBXC

#include "write_libxc_r.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_general/module_xc/xc_functional_libxc.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

#include <xc.h>

#include <vector>

void ModuleIO::write_libxc_r(
        const int order,
		const std::vector<int> &func_id,
        const int &nrxx, // number of real-space grid
        const double &omega, // volume of cell
        const double tpiba,
        const Charge* const chr)
{
    ModuleBase::TITLE("ModuleIO","write_libxc_r");
    ModuleBase::timer::tick("ModuleIO","write_libxc_r");

    const int nspin =
        (GlobalV::NSPIN == 1 || ( GlobalV::NSPIN ==4 && !GlobalV::DOMAG && !GlobalV::DOMAG_Z))
        ? 1 : 2;

    //----------------------------------------------------------
    // xc_func_type is defined in Libxc package
    // to understand the usage of xc_func_type,
    // use can check on website, for example:
    // https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/
    //----------------------------------------------------------

    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func( func_id, (1==nspin) ? XC_UNPOLARIZED:XC_POLARIZED );

    const bool is_gga = [&funcs]()
    {
        for( xc_func_type &func : funcs )
        {
            switch( func.info->family )
            {
                case XC_FAMILY_GGA:
                case XC_FAMILY_HYB_GGA:
                    return true;
            }
        }
        return false;
    }();

    // converting rho
    std::vector<double> rho;
    std::vector<double> amag;
    if(1==nspin || 2==GlobalV::NSPIN)
    {
        rho = XC_Functional_Libxc::convert_rho(nspin, nrxx, chr);
    }
    else
    {
        std::tuple<std::vector<double>,std::vector<double>> rho_amag = XC_Functional_Libxc::convert_rho_amag_nspin4(nspin, nrxx, chr);
        rho = std::get<0>(std::move(rho_amag));
        amag = std::get<1>(std::move(rho_amag));
    }

    std::vector<std::vector<ModuleBase::Vector3<double>>> gdr;
    std::vector<double> sigma;
    if(is_gga)
    {
        gdr = XC_Functional_Libxc::cal_gdr(nspin, rho, tpiba, chr);
        sigma = XC_Functional_Libxc::convert_sigma(gdr);
    }

    for( xc_func_type &func : funcs )
    {
        // jiyy add for threshold
        constexpr double rho_threshold = 1E-6;
        constexpr double grho_threshold = 1E-10;

        xc_func_set_dens_threshold(&func, rho_threshold);

        // sgn for threshold mask
        const std::vector<double> sgn = XC_Functional_Libxc::cal_sgn(rho_threshold, grho_threshold, func, nspin, nrxx, rho, sigma);

		std::vector<double> exc;
		std::vector<double> vrho;
		std::vector<double> vsigma;
		std::vector<double> v2rho2;
		std::vector<double> v2rhosigma;
		std::vector<double> v2sigma2;
		std::vector<double> v3rho3;
		std::vector<double> v3rho2sigma;
		std::vector<double> v3rhosigma2;
		std::vector<double> v3sigma3;
		std::vector<double> v4rho4;
		std::vector<double> v4rho3sigma;
		std::vector<double> v4rho2sigma2;
		std::vector<double> v4rhosigma3;
		std::vector<double> v4sigma4;
		// attention: order 4321 don't break
		switch( order )
		{
			case 4:		v4rho4.resize( nrxx * ((1==nspin)?1:5) );
			case 3:		v3rho3.resize( nrxx * ((1==nspin)?1:4) );
			case 2:		v2rho2.resize( nrxx * ((1==nspin)?1:3) );
			case 1:		vrho  .resize( nrxx * nspin            );
			case 0:		exc   .resize( nrxx                    );
						break;
			default:	throw std::domain_error("order ="+std::to_string(order)
							+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
						break;
		}
		if(is_gga)
		{
			switch( order )
			{
				case 4:		v4rho3sigma .resize( nrxx * ((1==nspin)?1:12) );
							v4rho2sigma2.resize( nrxx * ((1==nspin)?1:15) );
							v4rhosigma3 .resize( nrxx * ((1==nspin)?1:20) );
							v4sigma4    .resize( nrxx * ((1==nspin)?1:15) );
				case 3:		v3rho2sigma .resize( nrxx * ((1==nspin)?1:9) );
							v3rhosigma2 .resize( nrxx * ((1==nspin)?1:12) );
							v3sigma3    .resize( nrxx * ((1==nspin)?1:10) );
				case 2:		v2rhosigma  .resize( nrxx * ((1==nspin)?1:6) );
							v2sigma2    .resize( nrxx * ((1==nspin)?1:6) );
				case 1:		vsigma      .resize( nrxx * ((1==nspin)?1:3) );
				case 0:		break;
				default:	throw std::domain_error("order ="+std::to_string(order)
								+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
							break;
			}
			break;
		}

		// call libxc function
		// attention: order 432 don't break
        switch( func.info->family )
        {
            case XC_FAMILY_LDA:
			{
				switch( order )
				{
					case 4:		xc_lda_lxc    ( &func, nrxx, rho.data(), v4rho4.data() );
					case 3:		xc_lda_kxc    ( &func, nrxx, rho.data(), v3rho3.data() );
					case 2:		xc_lda_fxc    ( &func, nrxx, rho.data(), v2rho2.data() );
					case 1:		xc_lda_exc_vxc( &func, nrxx, rho.data(), exc.data(), vrho.data() );
								break;
					case 0:		xc_lda_exc    ( &func, nrxx, rho.data(), exc.data() );
								break;
					default:	throw std::domain_error("order ="+std::to_string(order)
									+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
								break;
				}
                break;
			}
            case XC_FAMILY_GGA:
            case XC_FAMILY_HYB_GGA:
			{
				switch( order )
				{
					case 4:		xc_gga_lxc    ( &func, nrxx, rho.data(), sigma.data(), v4rho4.data(), v4rho3sigma.data(), v4rho2sigma2.data(), v4rhosigma3.data(), v4sigma4.data() );
					case 3:		xc_gga_kxc    ( &func, nrxx, rho.data(), sigma.data(), v3rho3.data(), v3rho2sigma.data(), v3rhosigma2.data(), v3sigma3.data() );
					case 2:		xc_gga_fxc    ( &func, nrxx, rho.data(), sigma.data(), v2rho2.data(), v2rhosigma.data(), v2sigma2.data() );
					case 1:		xc_gga_exc_vxc( &func, nrxx, rho.data(), sigma.data(), exc.data(), vrho.data(), vsigma.data() );
								break;
					case 0:		xc_gga_exc    ( &func, nrxx, rho.data(), sigma.data(), exc.data() );
								break;
					default:	throw std::domain_error("order ="+std::to_string(order)
									+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
								break;
				}
				break;
			}
            default:
			{
                throw std::domain_error("func.info->family ="+std::to_string(func.info->family)
                    +" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
                break;
			}
        }
    } // end for( xc_func_type &func : funcs )

    XC_Functional_Libxc::finish_func(funcs);

    ModuleBase::timer::tick("ModuleIO","write_libxc_r");
}

#endif