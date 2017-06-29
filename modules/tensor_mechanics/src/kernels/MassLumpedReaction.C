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

#include "MassLumpedReaction.h"
#include "Assembly.h"
#include "MooseVariable.h"


// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<MassLumpedReaction>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

MassLumpedReaction::MassLumpedReaction(const InputParameters & parameters) :
    Kernel(parameters),
    _density(getMaterialProperty<Real>("density")),
    _u_nodal(_var.nodalValue())
{
}

Real
MassLumpedReaction::computeQpResidual()
{
  return _density[_qp] * _test[_i][_qp] * _u_nodal[_i];
}

Real
MassLumpedReaction::computeQpJacobian()
{
  return _density[_qp] * _test[_i][_qp];
}

void
MassLumpedReaction::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());

  for (_i = 0; _i < _test.size(); _i++)
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
}

