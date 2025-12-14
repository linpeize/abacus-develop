//==========================================================
// AUTHOR : Peize Lin
// DATE : 2025-12-14
//==========================================================

#pragma once

#include "element_basis_index.h"

namespace ModuleBase
{

template<typename Tkey>
std::map<Tkey, Element_Basis_Index::Index_T>
Element_Basis_Index::construct_index( const std::map<Tkey, std::vector<NM>> &range )
{
	std::map<Tkey, Element_Basis_Index::Index_T> index(range.size());
	for( const auto &range_T : range )
		{ index[range_T.first] = construct_index(range_T.second); }
	return index;
}

}