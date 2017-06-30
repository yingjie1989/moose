/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ContactTractionAux.h"

#include "NodalArea.h"
#include "PenetrationLocator.h"

#include "libmesh/string_to_enum.h"

template <>
InputParameters
validParams<ContactTractionAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("nodal_area", "The nodal area");
  params.addParam<int>("component",-1,"The direction of the traction");
  params.addRequiredParam<BoundaryName>("paired_boundary", "The boundary to be penetrated");
  params.set<MultiMooseEnum>("execute_on") = "nonlinear";
  return params;
}

ContactTractionAux::ContactTractionAux(const InputParameters & params)
  : AuxKernel(params),
    _nodal_area(coupledValue("nodal_area")),
    _component(getParam<int>("component")),
    _penetration_locator(
        getPenetrationLocator(getParam<BoundaryName>("paired_boundary"),
                              getParam<std::vector<BoundaryName>>("boundary")[0],
                              Utility::string_to_enum<Order>(getParam<MooseEnum>("order"))))
{
}

ContactTractionAux::~ContactTractionAux() {}

Real
ContactTractionAux::computeValue()
{
  Real value(0);
  const Real area = _nodal_area[_qp];
  const PenetrationInfo * pinfo(NULL);

  const auto it = _penetration_locator._penetration_info.find(_current_node->id());
  if (it != _penetration_locator._penetration_info.end())
    pinfo = it->second;

  if (pinfo && area != 0){
    Real contact_force_normal = pinfo->_contact_force * pinfo->_normal;

    RealVectorValue traction = pinfo->_contact_force - (pinfo->_normal * contact_force_normal);

    if (_component > -1)
        value = traction(_component) / area;
   else
    //RealVectorValue
    value = traction.norm() / area;
    //traction.norm2();

  }

  return value;
}
