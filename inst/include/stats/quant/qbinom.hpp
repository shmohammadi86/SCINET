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
 * quantile function of the Binomial distribution
 */

#ifndef _statslib_qbinom_HPP
#define _statslib_qbinom_HPP

//
// single input

template<typename Ta, typename Tb = uint_t>
statslib_constexpr
Tb qbinom(const Ta p, const uint_t n_trials_par, const Ta prob_par);

//
// matrix/vector input

template<typename Ta, typename Tb, typename Tc = Tb>
statslib_inline
void qbinom_int(const Ta* __stats_pointer_settings__ vals_in, const uint_t n_trials_par, const Tb prob_par, 
                      Tc* __stats_pointer_settings__ vals_out, const uint_t num_elem);

#ifdef STATS_USE_ARMA
template<typename Ta, typename Tb, typename Tc = Tb>
statslib_inline
ArmaMat<Tc> qbinom(const ArmaMat<Ta>& X, const uint_t n_trials_par, const Tb prob_par);
#endif

#ifdef STATS_USE_BLAZE
template<typename Ta, typename Tb, typename Tc = Tb, bool To = blaze::columnMajor>
statslib_inline
BlazeMat<Tc,To> qbinom(const BlazeMat<Ta,To>& X, const uint_t n_trials_par, const Tb prob_par);
#endif

#ifdef STATS_USE_EIGEN
template<typename Ta, typename Tb, typename Tc = Tb, int iTr = Eigen::Dynamic, int iTc = Eigen::Dynamic>
statslib_inline
EigMat<Tc,iTr,iTc> qbinom(const EigMat<Ta,iTr,iTc>& X, const uint_t n_trials_par, const Tb prob_par);
#endif

//
// include implementation files

#include "qbinom.ipp"

#endif
