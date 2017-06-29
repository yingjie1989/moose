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

#ifndef MASSLUMPEDREACTION_H
#define MASSLUMPEDREACTION_H

#include "Kernel.h"

// Forward Declaration
class MassLumpedReaction;

template<>
InputParameters validParams<MassLumpedReaction>();

class MassLumpedReaction : public Kernel
{
public:
  MassLumpedReaction(const InputParameters & parameters);

  virtual void computeJacobian() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  const MaterialProperty<Real> & _density;
  const VariableValue & _u_nodal;
};

#endif // MASSLUMPEDREACTION_H
