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
 * Sample from a normal distribution
 */

#ifndef _statslib_rnorm_HPP
#define _statslib_rnorm_HPP

template<typename T>
statslib_inline
return_t<T> rnorm(const T mu_par, const T sigma_par, rand_engine_t& engine);

template<typename T>
statslib_inline
return_t<T> rnorm(const T mu_par, const T sigma_par, uint_t seed_val = std::random_device{}());

template<typename T = double>
statslib_inline
T rnorm();

template<typename T>
statslib_inline
void rnorm_int(const T mu_par, const T sigma_par, T* vals_out, const uint_t num_elem);

template<typename mT, typename eT = double>
statslib_inline
mT rnorm(const uint_t n, const uint_t k, const eT mu_par = eT(0), const eT sigma_par = eT(1));

#include "rnorm.ipp"

#endif
