//-----------------------------------------------------------------------------
// boost/variant/detail/bool_trait_def.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2003
// Eric Friedman
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

// Needed until move-related traits incorporated into type_traits library.
// no include guards, the header is intended for multiple inclusion!

// should be the last #include
#include "boost/type_traits/detail/bool_trait_def.hpp"

#define BOOST_VARIANT_TT_AUX_BOOL_TRAIT_DEF1(trait,T,C) \
template< typename T > struct trait \
    BOOST_TT_AUX_BOOL_C_BASE(C) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1,trait,(T)) \
}; \
/**/

#define BOOST_VARIANT_TT_AUX_TRAIT_SUFFIX(arity, name) \
BOOST_TT_AUX_TEMPLATE_ARITY_SPEC(arity, name) \
/**/
