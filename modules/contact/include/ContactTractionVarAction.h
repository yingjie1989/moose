/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CONTACTTRACTIONVARACTION_H
#define CONTACTTRACTIONVARACTION_H

#include "Action.h"
#include "MooseTypes.h"

class ContactTractionVarAction : public Action
{
public:
  ContactTractionVarAction(const InputParameters & params);

  virtual void act();
};

template <>
InputParameters validParams<ContactTractionVarAction>();

#endif
