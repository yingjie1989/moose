/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CONTACTTRACTIONAUX_H
#define CONTACTTRACTIONAUX_H

#include "AuxKernel.h"

class NodalArea;
class PenetrationLocator;

class ContactTractionAux : public AuxKernel
{
public:
  ContactTractionAux(const InputParameters & parameters);

  virtual ~ContactTractionAux();

protected:
  virtual Real computeValue();

  const VariableValue & _nodal_area;
  int _component;
  const PenetrationLocator & _penetration_locator;
};

template <>
InputParameters validParams<ContactTractionAux>();

#endif
