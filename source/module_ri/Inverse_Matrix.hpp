//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef INVERSE_MATRIX_HPP
#define INVERSE_MATRIX_HPP

#include "Inverse_Matrix.h"
#include "module_base/lapack_connector.h"
#include <RI/global/Blas_Interface-Contiguous.h>
#include <RI/global/Blas_Interface-Tensor.h>
#include <RI/global/Lapack_Interface-Tensor.h>

#include <cassert>

template<typename Tdata>
void Inverse_Matrix<Tdata>::cal_inverse( const Method &method )
{
	switch(method)
	{
		case Method::potrf:		using_potrf();			break;
		case Method::heev:		using_heev(1E-6);		break;
	}
}

template<typename Tdata>
void Inverse_Matrix<Tdata>::using_potrf()
{
	int info;
	LapackConnector::potrf('U', this->A.shape[0], this->A.ptr(), this->A.shape[0], info);
	if(info)
		throw std::range_error("info="+std::to_string(info)+"\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));

	LapackConnector::potri('U', this->A.shape[0], this->A.ptr(), this->A.shape[0], info);
	if(info)
		throw std::range_error("info="+std::to_string(info)+"\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));

	copy_down_triangle();
}

template<typename Tdata>
void Inverse_Matrix<Tdata>::using_heev( const RI::Global_Func::To_Real_t<Tdata> &threshold_condition_number )
{
	using Tdata_real = RI::Global_Func::To_Real_t<Tdata>;
	std::vector<Tdata_real> eigen_value(this->A.shape[1]);

	const int info = RI::Lapack_Interface::heev('V', 'U', this->A, eigen_value);
	if(info)
		throw std::range_error("info="+std::to_string(info)+"\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));

	Tdata_real eigen_value_max = 0;
	for( const Tdata_real &ie : eigen_value )
		eigen_value_max = std::max( ie, eigen_value_max );
	const RI::Global_Func::To_Real_t<Tdata> threshold = eigen_value_max * threshold_condition_number;

	RI::Tensor<Tdata> eA( this->A.shape );
	int ie=0;
	for( int i=0; i!=this->A.shape[1]; ++i )
		if( eigen_value[i] > threshold )
		{
			RI::Blas_Interface::axpy( this->A.shape[1], Tdata(std::sqrt(1.0/eigen_value[i])), this->A.ptr()+i*this->A.shape[1], eA.ptr()+ie*eA.shape[1] );
			++ie;
		}
	this->A = RI::Blas_Interface::gemm('T', 'N', Tdata(1), eA, eA);
}

template<typename Tdata>
void Inverse_Matrix<Tdata>::input( const RI::Tensor<Tdata> &m )
{
	assert(m.shape.size()==2);
	assert(m.shape[0]==m.shape[1]);
	this->A = m.copy();
}


template<typename Tdata>
void Inverse_Matrix<Tdata>::input(const std::vector<std::vector<RI::Tensor<Tdata>>> &ms)
{
	const size_t N0 = ms.size();
	assert(N0>0);
	const size_t N1 = ms[0].size();
	assert(N1>0);
	for(size_t Im0=0; Im0<N0; ++Im0)
		assert(ms[Im0].size()==N1);

	for(size_t Im0=0; Im0<N0; ++Im0)
		for(size_t Im1=0; Im1<N1; ++Im1)
			assert(ms[Im0][Im1].shape.size()==2);

	std::vector<size_t> n0(N0);
	for(size_t Im0=0; Im0<N0; ++Im0)
		n0[Im0] = ms[Im0][0].shape[0];
	std::vector<size_t> n1(N1);
	for(size_t Im1=0; Im1<N1; ++Im1)
		n1[Im1] = ms[0][Im1].shape[1];

	for(size_t Im0=0; Im0<N0; ++Im0)
		for(size_t Im1=0; Im1<N1; ++Im1)
			assert((ms[Im0][Im1].shape[0]==n0[Im0]) && (ms[Im0][Im1].shape[1]==n1[Im1]));

	const size_t n_all = std::accumulate(n0.begin(), n0.end(), 0);
	assert(n_all == std::accumulate(n1.begin(), n1.end(), 0));
	this->A = RI::Tensor<Tdata>({n_all, n_all});

	std::vector<size_t> n0_partial(N0+1);
	std::partial_sum(n0.begin(), n0.end(), n0_partial.begin()+1);
	std::vector<size_t> n1_partial(N1+1);
	std::partial_sum(n1.begin(), n1.end(), n1_partial.begin()+1);

	for(size_t Im0=0; Im0<N0; ++Im0)
		for(size_t Im1=0; Im1<N1; ++Im1)
		{
			const RI::Tensor<Tdata> &m_tmp = ms.at(Im0).at(Im1);
			for(size_t im0=0; im0<m_tmp.shape[0]; ++im0)
				for(size_t im1=0; im1<m_tmp.shape[1]; ++im1)
					this->A(im0+n0_partial[Im0], im1+n1_partial[Im1]) = m_tmp(im0,im1);
		}
}


template<typename Tdata>
RI::Tensor<Tdata> Inverse_Matrix<Tdata>::output() const
{
	return this->A.copy();
}


template<typename Tdata>
std::vector<std::vector<RI::Tensor<Tdata>>>
Inverse_Matrix<Tdata>::output(const std::vector<size_t> &n0, const std::vector<size_t> &n1) const
{
	assert( std::accumulate(n0.begin(), n0.end(), 0) == this->A.shape[0] );
	assert( std::accumulate(n1.begin(), n1.end(), 0) == this->A.shape[1] );

	const size_t N0 = n0.size();
	const size_t N1 = n1.size();

	std::vector<size_t> n0_partial(N0+1);
	std::partial_sum(n0.begin(), n0.end(), n0_partial.begin()+1);
	std::vector<size_t> n1_partial(N1+1);
	std::partial_sum(n1.begin(), n1.end(), n1_partial.begin()+1);

	std::vector<std::vector<RI::Tensor<Tdata>>> ms(N0, std::vector<RI::Tensor<Tdata>>(N1));
	for(size_t Im0=0; Im0<N0; ++Im0)
		for(size_t Im1=0; Im1<N1; ++Im1)
		{
			RI::Tensor<Tdata> &m_tmp = ms[Im0][Im1] = RI::Tensor<Tdata>({n0[Im0], n1[Im1]});
			for(size_t im0=0; im0<n0[Im0]; ++im0)
				for(size_t im1=0; im1<n1[Im1]; ++im1)
					m_tmp(im0,im1) = this->A(im0+n0_partial[Im0], im1+n1_partial[Im1]);
		}
	return ms;
}


template<typename Tdata>
void Inverse_Matrix<Tdata>::copy_down_triangle()
{
	for( size_t i0=0; i0<this->A.shape[0]; ++i0 )
		for( size_t i1=0; i1<i0; ++i1 )
			this->A(i0,i1) = this->A(i1,i0);
}

#endif