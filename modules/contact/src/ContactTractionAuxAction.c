/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ContactTractionAuxAction.h"

#include "Factory.h"
#include "FEProblem.h"
#include "Parser.h"
#include "MooseApp.h"
#include "Conversion.h"

static unsigned int counter = 0;

template <>
InputParameters
validParams<ContactTractionAuxAction>()
{
  MooseEnum orders("FIRST SECOND THIRD FOURTH", "FIRST");

  InputParameters params = validParams<Action>();
  params.addRequiredParam<BoundaryName>("master", "The master surface");
  params.addRequiredParam<BoundaryName>("slave", "The slave surface");
  params.addParam<MooseEnum>("order", orders, "The finite element order: " + orders.getRawNames());
  return params;
}

ContactTractionAuxAction::ContactTractionAuxAction(const InputParameters & params)
  : Action(params),
    _master(getParam<BoundaryName>("master")),
    _slave(getParam<BoundaryName>("slave")),
    _order(getParam<MooseEnum>("order"))
{
}

void
ContactTractionAuxAction::act()
{
  if (!_problem->getDisplacedProblem())
  {
    mooseError("Contact requires updated coordinates.  Use the 'displacements = ...' line in the "
               "Mesh block.");
  }

  {
    InputParameters params = _factory.getValidParams("ContactTractionAux");

    // Extract global params
    if (isParamValid("parser_syntax"))
      _app.parser().extractParams(getParam<std::string>("parser_syntax"), params);

    params.set<std::vector<BoundaryName>>("boundary") = {_slave};
    params.set<BoundaryName>("paired_boundary") = _master;
    params.set<AuxVariableName>("variable") = "contact_traction";
    params.addRequiredCoupledVar("nodal_area", "The nodal area");
    params.set<std::vector<VariableName>>("nodal_area") = {"nodal_area_" + _name};
    params.set<MooseEnum>("order") = _order;

    params.set<bool>("use_displaced_mesh") = true;

    std::stringstream name;
    name << _name;
    name << "_contact_traction_";
    name << counter++;

    MultiMooseEnum execute_options(SetupInterface::getExecuteOptions());
    execute_options = "nonlinear timestep_end timestep_begin";
    params.set<MultiMooseEnum>("execute_on") = execute_options;
    _problem->addAuxKernel("ContactTractionAux", name.str(), params);
  }
}
