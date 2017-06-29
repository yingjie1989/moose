/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ContactAugLagMulProblem.h"

// MOOSE includes
#include "AuxiliarySystem.h"
#include "DisplacedProblem.h"
#include "MooseApp.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "NearestNodeLocator.h"
#include "NonlinearSystem.h"
#include "PenetrationLocator.h"


#include "SystemBase.h"
#include "Assembly.h"
#include "Executioner.h"
#include "AddVariableAction.h"

// libMesh includes
#include "libmesh/string_to_enum.h"
#include "libmesh/sparse_matrix.h"


#include <limits>

template <>
InputParameters
validParams<ContactAugLagMulProblem>()
{
  InputParameters params = validParams<ReferenceResidualProblem>();
  params.addRequiredParam<std::vector<int>>(
      "master", "IDs of the master surfaces for which the slip should be calculated");
  params.addRequiredParam<std::vector<int>>(
      "slave", "IDs of the slave surfaces for which the slip should be calculated");
  params.addParam<Real>(
          "penalty",
          1e8,
          "The penalty to apply.  This can vary depending on the stiffness of your materials");
  params.addParam<bool>(
              "normalize_penalty",
              false,
              "Whether to normalize the penalty parameter with the nodal area for penalty contact.");
  params.addRequiredParam<NonlinearVariableName>("disp_x",
                                                 "Variable containing the x displacement");
  params.addRequiredParam<NonlinearVariableName>("disp_y",
                                                 "Variable containing the y displacement");
  params.addParam<NonlinearVariableName>("disp_z", "Variable containing the z displacement");
  params.addRequiredParam<AuxVariableName>("residual_x",
                                           "Auxiliary variable containing the saved x residual");
  params.addRequiredParam<AuxVariableName>("residual_y",
                                           "Auxiliary variable containing the saved y residual");
  params.addParam<AuxVariableName>("residual_z",
                                   "Auxiliary variable containing the saved z residual");

  params.addParam<int>("minimum_update_iterations", 1, "Minimum number of update Lagrangian Multiplier iterations per step");
  params.addParam<int>(
      "maximum_update_iterations", 100, "Maximum number of update Lagrangian Multiplier iterations per step");

  params.addParam<Real>("target_contact_residual",
                        "Augmented Lagrangian Multiplier contact residual convergence criterion");

  params.addParam<Real>("target_relative_contact_residual",
                                              "Augmented Lagrangian Multiplier contact residual convergence criterion");
  params.addParam<Real>("contact_lagmul_tolerance_factor",
                                              "Augmented Lagrangian Multiplier tolerance factor");

  params.addParam<Real>("target_relative_contact_residual",
                        "Frictional contact relative residual convergence criterion");
  params.addParam<std::vector<std::string>>(
      "contact_reference_residual_variables",
      "Set of variables that provide reference residuals for relative contact convergence check");


  return params;
}


ContactAugLagMulProblem::ContactAugLagMulProblem(const InputParameters & params)
  : ReferenceResidualProblem(params),
    _penalty(getParam<Real>("penalty")),
    _normalize_penalty(getParam<bool>("normalize_penalty")),
    //_nodal_area_var(getVar("nodal_area", 0)),
    //_nodal_area_var(0),
    //_aux_system(getAuxiliarySystem()),
    //_aux_solution(_aux_system.currentSolution()),
    _refResidContact(0.0),
    _do_lagmul_update(false),
    _num_lagmul_iterations(0),
    _target_contact_residual(0.0),
    _target_relative_contact_residual(0.0),
    _contact_lagmul_tol_factor(1.0),
    _num_nl_its_since_contact_update(0),
    _num_contact_nodes(0)
{

  _console << "Initialize the LM Problem\n";

//  AuxiliarySystem & aux_sys = getAuxiliarySystem();


  std::vector<int> master = params.get<std::vector<int>>("master");
  std::vector<int> slave = params.get<std::vector<int>>("slave");

  unsigned int dim = getNonlinearSystemBase().subproblem().mesh().dimension();

  _disp_x = params.get<NonlinearVariableName>("disp_x");
  _residual_x = params.get<AuxVariableName>("residual_x");

  _disp_y = params.get<NonlinearVariableName>("disp_y");
  _residual_y = params.get<AuxVariableName>("residual_y");

  if (dim == 3)
  {
    if (!params.isParamValid("disp_z"))
      mooseError("Missing disp_z in FrictionalContactProblem");
    if (!params.isParamValid("residual_z"))
      mooseError("Missing residual_z in FrictionalContactProblem");
    _disp_z = params.get<NonlinearVariableName>("disp_z");
    _residual_z = params.get<AuxVariableName>("residual_z");
  }

  unsigned int num_interactions = master.size();
  if (num_interactions != slave.size())
    mooseError(
        "Sizes of master surface and slave surface lists must match in FrictionalContactProblem");

  for (unsigned int i = 0; i < master.size(); ++i)
  {
    std::pair<int, int> ms_pair(master[i], slave[i]);
    InteractionParams ip;

    _interaction_params[ms_pair] = ip;
  }

  _min_lagmul_iters = params.get<int>("minimum_update_iterations");
  _max_lagmul_iters = params.get<int>("maximum_update_iterations");
  //_lagmul_updates_per_iter = params.get<int>("lagmul_updates_per_iteration");

  bool have_target = false;
  bool have_target_relative = false;
  if (params.isParamValid("target_contact_residual"))
  {
    _target_contact_residual = params.get<Real>("target_contact_residual");
    have_target = true;
  }
  if (params.isParamValid("target_relative_contact_residual"))
  {
    _target_relative_contact_residual = params.get<Real>("target_relative_contact_residual");
    have_target_relative = true;
  }
  if (!(have_target || have_target_relative))
    mooseError("Must specify either target_contact_residual or target_relative_contact_residual");

  //_contact_lagmul_tol_factor = params.get<Real>("contact_lagmul_tolerance_factor");

  if (params.isParamValid("contact_reference_residual_variables"))
    _contactRefResidVarNames =
        params.get<std::vector<std::string>>("contact_reference_residual_variables");
}

void
ContactAugLagMulProblem::initialSetup()
{
  _console << "the LM Problem initial\n";

  ReferenceResidualProblem::initialSetup();
}

void
ContactAugLagMulProblem::timestepSetup()
{

  _console << "the LM Problem timestepSetup\n";

  _do_lagmul_update = false;
  _num_lagmul_iterations = 0;
  _num_nl_its_since_contact_update = 0;
  _refResidContact = 0.0;
  ReferenceResidualProblem::timestepSetup();
}

void
ContactAugLagMulProblem::updateContactReferenceResidual()
{

  _console << "the LM Problem update\n";

  if (_contactRefResidVarIndices.size() > 0)
  {
    _refResidContact = 0.0;
    for (unsigned int i = 0; i < _contactRefResidVarIndices.size(); ++i)
      _refResidContact += _refResid[i] * _refResid[i];

    _refResidContact = std::sqrt(_refResidContact);
  }
  else if (_refResid.size() > 0)
  {
    _refResidContact = 0.0;
    for (unsigned int i = 0; i < _refResid.size(); ++i)
      _refResidContact += _refResid[i] * _refResid[i];

    _refResidContact = std::sqrt(_refResidContact);
  }
  _console << "Contact reference convergence residual: " << _refResidContact << std::endl;
}



MooseNonlinearConvergenceReason
ContactAugLagMulProblem::checkNonlinearConvergence(std::string & msg,
                                                    const PetscInt it,
                                                    const Real xnorm,
                                                    const Real snorm,
                                                    const Real fnorm,
                                                    const Real rtol,
                                                    const Real stol,
                                                    const Real abstol,
                                                    const PetscInt nfuncs,
                                                    const PetscInt /*max_funcs*/,
                                                    const Real ref_resid,
                                                    const Real /*div_threshold*/)
{


  Real my_max_funcs = std::numeric_limits<int>::max();
  Real my_div_threshold = std::numeric_limits<Real>::max();

  MooseNonlinearConvergenceReason reason =
      ReferenceResidualProblem::checkNonlinearConvergence(msg,
                                                          it,
                                                          xnorm,
                                                          snorm,
                                                          fnorm,
                                                          rtol,
                                                          stol,
                                                          abstol,
                                                          nfuncs,
                                                          my_max_funcs,
                                                          ref_resid,
                                                          my_div_threshold);

  _refResidContact = ref_resid; // use initial residual if no reference variables are specified
  updateContactReferenceResidual();

  int min_nl_its_since_contact_update = 1;
  ++_num_nl_its_since_contact_update;

  if ((reason > 0) ||                         // converged
      (reason == MOOSE_NONLINEAR_ITERATING && // iterating and converged within factor
       (fnorm < abstol * _contact_lagmul_tol_factor ||
        checkConvergenceIndividVars(
            fnorm, abstol * _contact_lagmul_tol_factor, rtol * _contact_lagmul_tol_factor, ref_resid))))
  {
    _console << "Augmentd Lagrangian Multiplier iteration " << _num_lagmul_iterations << " ";
    if (_num_lagmul_iterations < _max_lagmul_iters)
    { // do a slip update if there is another iteration
        _do_lagmul_update = true;

        NonlinearSystemBase & nonlinear_sys = getNonlinearSystemBase();
        nonlinear_sys.update();
        const NumericVector<Number> *& ghosted_solution = nonlinear_sys.currentSolution();

        //calculateSlip(*ghosted_solution, nullptr); // Just to calculate slip residual
        //Real _pene_residual = calculatePenetration(*ghosted_solution);

        bool _augLM_repeat_step;

        if (_displaced_problem != NULL)
              _augLM_repeat_step = nonlinear_sys.updateLagMul(true);
        else
              _augLM_repeat_step = nonlinear_sys.updateLagMul(false);

        if (_augLM_repeat_step)
        { // force it to keep iterating
          reason = MOOSE_NONLINEAR_ITERATING;
          _console << "Contact penetration > tolerance" << std::endl;
          _num_lagmul_iterations ++;
        }
        else
        {
          //_do_slip_update = false; //maybe we want to do this
          _console << "Contact penetration pene_resid < tolerance" << std::endl;
        }
    }
    else
    { // maxed out
      _console << "Max lagmul iterations" << std::endl;
      reason = MOOSE_DIVERGED_FUNCTION_COUNT;
    }
  }

  return reason;
}
