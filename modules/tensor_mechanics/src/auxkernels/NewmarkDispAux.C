/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "NewmarkDispAux.h"

template<>
InputParameters validParams<NewmarkDispAux>()
{
  InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("acceleration","acceleration variable");
    params.addRequiredCoupledVar("velocity","velocity variable");
    params.addRequiredParam<Real>("gamma","gamma parameter");
    params.addRequiredParam<Real>("beta","beta parameter");

  return params;
}

NewmarkDispAux::NewmarkDispAux(const InputParameters & parameters) :
  AuxKernel(parameters),
   _accel_old(coupledValueOld("acceleration")),
   _accel(coupledValue("acceleration")),
   _vel_old(coupledValueOld("velocity")),
   _vel(coupledValue("velocity")),
   _beta(getParam<Real>("beta")),
   _gamma(getParam<Real>("gamma"))
{
}

Real
NewmarkDispAux::computeValue()
{
  Real disp_old = _u_old[_qp];
  if (!isNodal())
    mooseError("must run on a nodal variable");
  // Calculates Velocity using Newmark time integration scheme
  //return vel_old + (_dt*(1-_gamma))*_accel_old[_qp] + _gamma*_dt*_accel[_qp];
   
  return disp_old + _dt*_vel_old[_qp] + 0.5*_dt*_dt*( (1-2*_beta)*_accel_old[_qp] + 2*_beta*_accel[_qp]);
}
