#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4LGTe_mod) {


    class_<rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> >("model_LGTe")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_LGTe_namespace::model_LGTe, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4SGTe_mod) {


    class_<rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> >("model_SGTe")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_SGTe_namespace::model_SGTe, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4lgt_mod) {


    class_<rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> >("model_lgt")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_lgt_namespace::model_lgt, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4msgt_mod) {


    class_<rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> >("model_msgt")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_msgt_namespace::model_msgt, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sgt_mod) {


    class_<rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> >("model_sgt")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sgt_namespace::model_sgt, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4trend_mod) {


    class_<rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> >("model_trend")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_trend_namespace::model_trend, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
