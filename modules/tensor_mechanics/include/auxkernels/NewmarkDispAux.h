/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NEWMARKDISPAUX_H
#define NEWMARKDISPAUX_H

#include "AuxKernel.h"

class NewmarkDispAux;

template<>
InputParameters validParams<NewmarkDispAux>();

class NewmarkDispAux : public AuxKernel
{
public:

  /**
   * Calcualtes velocity using Newmark time integration scheme
   */
  NewmarkDispAux(const InputParameters & parameters);

  virtual ~NewmarkDispAux() {}

protected:
  virtual Real computeValue();

  const VariableValue & _accel_old;
  const VariableValue & _accel;
  const VariableValue & _vel_old;
  const VariableValue & _vel;
  Real _beta;
  Real _gamma;
};

#endif //NEWMARKDISPAUX_H
