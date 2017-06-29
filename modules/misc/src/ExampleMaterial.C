/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ExampleMaterial.h"

template <>
InputParameters
validParams<ExampleMaterial>()
{
  InputParameters params = validParams<Material>();
  //params.addParam<Real>("initial_diffusivity", 1.0, "The Initial Diffusivity");
  params.addCoupledVar("temperature", 0, "tempature to be coupled");

  params.addParam<Real>("k0_param", 1.0, "k0 parameters");
  params.addParam<Real>("Q_param", 1.0, "Q parameters");
  params.addParam<Real>("n_param", 1.0, "n parameters");

  return params;
}

ExampleMaterial::ExampleMaterial(const InputParameters & parameters)
  : Material(parameters),

     _k0(getParam<Real>("k0_param")),
     _Q(getParam<Real>("Q_param")),
     _n(getParam<Real>("n_param")),
     _kB(getMaterialPropertyOld<Real>("kB_param")),

    // Get a parameter value for the diffusivity
//    _initial_diffusivity(getParam<Real>("initial_diffusivity")),
    _temp(coupledValue("temperature")),
    // Declare that this material is going to have a Real
    // valued property named "diffusivity" that Kernels can use.
    _S_prop(declareProperty<Real>("S_prop"))
    // Retrieve/use an old value of diffusivity.
    // Note: this is _expensive_ as we have to store values for each
    // qp throughout the mesh. Only do this if you REALLY need it!
   // _diffusivity_old(getMaterialPropertyOld<Real>("diffusivity"))
{
}

void
ExampleMaterial::initQpStatefulProperties()
{
  // init the diffusivity property (this will become
  // _diffusivity_old in the first call of computeProperties)
  _S_prop[_qp] = 0.0;
}

void
ExampleMaterial::computeQpProperties()
{

  Real _k_prop = _k0 * std::exp( - _Q/_kB[_qp] * _temp[_qp]  );

  Real _funckt = _k_prop*_t;

  _S_prop[_qp] = 1 - std::exp( - std::pow(_funckt,_n) );

}
