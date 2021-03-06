/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "SecondDerivativeImplicitEuler.h"
#include "SubProblem.h"

template <>
InputParameters
validParams<SecondDerivativeImplicitEuler>()
{
  InputParameters params = validParams<TimeKernel>();
  return params;
}

SecondDerivativeImplicitEuler::SecondDerivativeImplicitEuler(const InputParameters & parameters)
  : TimeKernel(parameters), _u_old(valueOld()), _u_older(valueOlder())
{
}

Real
SecondDerivativeImplicitEuler::computeQpResidual()
{
  return _test[_i][_qp] * ((_u[_qp] - 2 * _u_old[_qp] + _u_older[_qp]) / (_dt * _dt));
}

Real
SecondDerivativeImplicitEuler::computeQpJacobian()
{
  return _test[_i][_qp] * (_phi[_j][_qp] / (_dt * _dt));
}
