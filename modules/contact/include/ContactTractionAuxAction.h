/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CONTACTTRACTIONAUXACTION_H
#define CONTACTTRACTIONAUXACTION_H

#include "Action.h"
#include "MooseTypes.h"
#include "MooseEnum.h"

class ContactTractionAuxAction : public Action
{
public:
  ContactTractionAuxAction(const InputParameters & params);

  virtual void act();

private:
  const BoundaryName _master;
  const BoundaryName _slave;
  const MooseEnum _order;
};

template <>
InputParameters validParams<ContactTractionAuxAction>();

#endif
