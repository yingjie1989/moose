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
#ifndef COEFMASSREACTION_H
#define COEFMASSREACTION_H

#include "Reaction.h"

// Forward Declarations
class CoefMassReaction;

template <>
InputParameters validParams<CoefMassReaction>();

class CoefMassReaction : public Reaction
{
public:
  CoefMassReaction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  Real _coef;
};

#endif // COEFMASSREACTION_H
