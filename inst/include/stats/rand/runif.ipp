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
 * Sample from a uniform distribution
 */

template<typename T>
statslib_inline
T
runif_int(const T a_par, const T b_par, rand_engine_t& engine)
{
    // convert from [a,b) to (a,b)
    T a_par_adj = std::nextafter(a_par, b_par);
    std::uniform_real_distribution<T> unif_dist(a_par_adj, b_par);

    return unif_dist(engine);
}

template<typename T>
statslib_inline
return_t<T>
runif(const T a_par, const T b_par, rand_engine_t& engine)
{
    return runif_int<return_t<T>>(a_par,b_par,engine);
}

template<typename T>
statslib_inline
return_t<T>
runif(const T a_par, const T b_par, uint_t seed_val)
{
    rand_engine_t engine(seed_val);
    return runif_int<return_t<T>>(a_par,b_par,engine);
}

template<typename T>
statslib_inline
T
runif()
{
    return runif<T>(T(0),T(1));
}

//

template<typename T>
statslib_inline
void
runif_int(const T a_par, const T b_par, T* vals_out, const uint_t num_elem, const uint_t seed_val)
{
    rand_engine_t engine(seed_val);

    for (uint_t j=0U; j < num_elem; j++)
    {
        vals_out[j] = runif(a_par,b_par,engine);
    }
}

#ifdef STATS_WITH_MATRIX_LIB
template<typename mT, typename eT>
statslib_inline
mT
runif(const uint_t n, const uint_t k, const eT a_par, const eT b_par, const uint_t seed_val)
{
    mT mat_out(n,k);

    runif_int(a_par,b_par,mat_ops::get_mem_ptr(mat_out),n*mat_ops::spacing(mat_out), seed_val);

    return mat_out;
}
#endif
