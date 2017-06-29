/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergenceExpPFFracTensors.h"

template<>
InputParameters validParams<StressDivergenceExpPFFracTensors>()
{
  InputParameters params = validParams<StressDivergenceTensors>();
  params.addClassDescription("Stress divergence kernel for phase-field fracture: Additionally computes off diagonal damage dependent Jacobian components");
  params.addCoupledVar("c", 0.0, "Phase field damage variable: Used to indicate calculation of Off Diagonal Jacobian term");
  return params;
}


StressDivergenceExpPFFracTensors::StressDivergenceExpPFFracTensors(const InputParameters & parameters) :
    DerivativeMaterialInterface<StressDivergenceTensors>(parameters),
    _c_coupled(isCoupled("c")),
    _c_var(_c_coupled ? coupled("c") : 0),
    _d_stress_dc(getMaterialPropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name()))
{
}

Real
StressDivergenceExpPFFracTensors::computeQpJacobian()
{
	return 0.0;
}


Real
StressDivergenceExpPFFracTensors::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.0;
}
