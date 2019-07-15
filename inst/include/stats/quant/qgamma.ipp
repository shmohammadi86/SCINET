/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
  ##
  ##   This file is part of the StatsLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * quantile function of the gamma distribution
 */

//
// single input

template<typename T>
statslib_constexpr
T
qgamma_int(const T p, const T shape_par, const T scale_par)
{
    return ( scale_par*gcem::incomplete_gamma_inv(shape_par,p) );
}

template<typename T>
statslib_constexpr
T
qgamma_check(const T p, const T shape_par, const T scale_par)
{
    return ( STLIM<T>::epsilon() > p ? T(0) :
             //
             qgamma_int(p,shape_par,scale_par) );
}

template<typename Ta, typename Tb>
statslib_constexpr
Ta
qgamma(const Ta p, const Tb shape_par, const Tb scale_par)
{
    return qgamma_check<Ta>(p,shape_par,scale_par);
}

//
// matrix/vector input

template<typename Ta, typename Tb, typename Tc>
statslib_inline
void
qgamma_int(const Ta* __stats_pointer_settings__ vals_in, const Tb shape_par, const Tb scale_par, 
                 Tc* __stats_pointer_settings__ vals_out, const uint_t num_elem)
{
#ifdef STATS_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint_t j=0U; j < num_elem; j++)
    {
        vals_out[j] = qgamma(vals_in[j],shape_par,scale_par);
    }
}

#ifdef STATS_USE_ARMA
template<typename Ta, typename Tb, typename Tc>
statslib_inline
ArmaMat<Tc>
qgamma(const ArmaMat<Ta>& X, const Tb shape_par, const Tb scale_par)
{
    ArmaMat<Tc> mat_out(X.n_rows,X.n_cols);

    qgamma_int<Ta,Tb,Tc>(X.memptr(),shape_par,scale_par,mat_out.memptr(),mat_out.n_elem);

    return mat_out;
}
#endif

#ifdef STATS_USE_BLAZE
template<typename Ta, typename Tb, typename Tc, bool To>
statslib_inline
BlazeMat<Tc,To>
qgamma(const BlazeMat<Ta,To>& X, const Tb shape_par, const Tb scale_par)
{
    BlazeMat<Tc,To> mat_out(X.rows(),X.columns());

    qgamma_int<Ta,Tb,Tc>(X.data(),shape_par,scale_par,mat_out.data(),X.rows()*X.spacing());

    return mat_out;
}
#endif

#ifdef STATS_USE_EIGEN
template<typename Ta, typename Tb, typename Tc, int iTr, int iTc>
statslib_inline
EigMat<Tc,iTr,iTc>
qgamma(const EigMat<Ta,iTr,iTc>& X, const Tb shape_par, const Tb scale_par)
{
    EigMat<Tc,iTr,iTc> mat_out(X.rows(),X.cols());

    qgamma_int<Ta,Tb,Tc>(X.data(),shape_par,scale_par,mat_out.data(),mat_out.size());

    return mat_out;
}
#endif
