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
#include "CoefMassReaction.h"

template <>
InputParameters
validParams<CoefMassReaction>()
{
  InputParameters params = validParams<Reaction>();
  params.addParam<Real>("coefficient", 1.0, "Coefficient of the term");
  return params;
}

CoefMassReaction::CoefMassReaction(const InputParameters & parameters)
  : Reaction(parameters), _coef(getParam<Real>("coefficient"))
{
}

Real
CoefMassReaction::computeQpResidual()
{
  return _coef * Reaction::computeQpResidual();
}

Real
CoefMassReaction::computeQpJacobian()
{
  return _coef * Reaction::computeQpJacobian();
}
