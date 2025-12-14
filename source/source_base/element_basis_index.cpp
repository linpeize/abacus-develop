//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#include "element_basis_index.h"

namespace ModuleBase
{

Element_Basis_Index::Index_T
Element_Basis_Index::construct_index( const std::vector<NM> &range )
{
	Index_T index;
	std::size_t count=0;
	index.resize( range.size() );
	for( std::size_t L=0; L!=range.size(); ++L )
	{
		index[L].resize( range[L].N );
		for( std::size_t N=0; N!=range[L].N; ++N )
		{
			index[L][N].resize( range[L].M );
			for( std::size_t M=0; M!=range[L].M; ++M )
			{
				index[L][N][M] = count;
				++count;
			}
		}
		index[L].N = range[L].N;
		index[L].M = range[L].M;
	}
	index.count_size = count;
	return index;
}

Element_Basis_Index::IndexLNM
Element_Basis_Index::construct_index( const Range &range )
{
	IndexLNM index(range.size());
	for( std::size_t T=0; T!=range.size(); ++T )
		{ index[T] = construct_index(range[T]); }
	return index;
}

}