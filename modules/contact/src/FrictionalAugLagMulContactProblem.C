/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "FrictionalAugLagMulContactProblem.h"

// MOOSE includes
#include "AuxiliarySystem.h"
#include "DisplacedProblem.h"
#include "MooseApp.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "NearestNodeLocator.h"
#include "NonlinearSystem.h"
#include "PenetrationLocator.h"

#include <limits>

template <>
InputParameters
validParams<FrictionalAugLagMulContactProblem>()
{
  InputParameters params = validParams<ReferenceResidualProblem>();
  params.addRequiredParam<std::vector<int>>(
      "master", "IDs of the master surfaces for which the slip should be calculated");
  params.addRequiredParam<std::vector<int>>(
      "slave", "IDs of the slave surfaces for which the slip should be calculated");
  params.addRequiredParam<std::vector<Real>>(
      "friction_coefficient", "Coefficient of friction for sliding contact for each interaction");
  params.addRequiredParam<std::vector<Real>>(
      "slip_factor", "Fraction of calculated slip to be applied for each interaction");
  params.addRequiredParam<std::vector<Real>>("slip_too_far_factor",
                                             "Fraction of calculated slip to be applied for each "
                                             "interaction in the slipped-too-far state");
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
  params.addRequiredParam<AuxVariableName>(
      "diag_stiff_x", "Auxiliary variable containing the saved x diagonal stiffness");
  params.addRequiredParam<AuxVariableName>(
      "diag_stiff_y", "Auxiliary variable containing the saved y diagonal stiffness");
  params.addParam<AuxVariableName>("diag_stiff_z",
                                   "Auxiliary variable containing the saved z diagonal stiffness");
  params.addRequiredParam<AuxVariableName>("inc_slip_x",
                                           "Auxiliary variable to store the x incremental slip");
  params.addRequiredParam<AuxVariableName>("inc_slip_y",
                                           "Auxiliary variable to store the y incremental slip");
  params.addParam<AuxVariableName>("inc_slip_z",
                                   "Auxiliary variable to store the z incremental slip");
  params.addParam<int>("minimum_slip_iterations", 1, "Minimum number of slip iterations per step");
  params.addParam<int>(
      "maximum_slip_iterations", 100, "Maximum number of slip iterations per step");
  params.addParam<int>(
      "slip_updates_per_iteration", 1, "Number of slip updates per contact iteration");
  params.addParam<Real>("contact_slip_tol",
                            "Slipping tolerance for stickking part");
  params.addParam<Real>("target_contact_residual",
                        "Frictional contact residual convergence criterion");
  params.addParam<Real>("target_relative_contact_residual",
                        "Frictional contact relative residual convergence criterion");
  params.addParam<Real>("contact_slip_tolerance_factor",
                        10.0,
                        "Multiplier on convergence criteria to determine when to start slipping");
  params.addParam<std::vector<std::string>>(
      "contact_reference_residual_variables",
      "Set of variables that provide reference residuals for relative contact convergence check");
  return params;
}

FrictionalAugLagMulContactProblem::SlipData::SlipData(const Node * node, unsigned int dof, Real slip)
  : _node(node), _dof(dof), _slip(slip)
{
}

FrictionalAugLagMulContactProblem::SlipData::SlipData(const SlipData & sd)
  : _node(sd._node), _dof(sd._dof), _slip(sd._slip)
{
}

FrictionalAugLagMulContactProblem::SlipData::~SlipData() {}

FrictionalAugLagMulContactProblem::FrictionalAugLagMulContactProblem(const InputParameters & params)
  : ReferenceResidualProblem(params),
    _refResidContact(0.0),
    _slip_residual(0.0),
    _pene_residual(0.0),
    _do_slip_update(false),
    _do_lagmul_update(false),
    _num_iterations(0),
    _stick_slip_tol(0.0),
    _target_contact_residual(0.0),
    _target_relative_contact_residual(0.0),
    _num_nl_its_since_contact_update(0),
    _num_contact_nodes(0),
    _num_slipping(0),
    _num_slipped_too_far(0),
    _inc_slip_norm(0.0),
    _it_slip_norm(0.0)
{
  std::vector<int> master = params.get<std::vector<int>>("master");
  std::vector<int> slave = params.get<std::vector<int>>("slave");
  std::vector<Real> friction_coefficient = params.get<std::vector<Real>>("friction_coefficient");
  std::vector<Real> slip_factor = params.get<std::vector<Real>>("slip_factor");
  std::vector<Real> slip_too_far_factor = params.get<std::vector<Real>>("slip_too_far_factor");

  unsigned int dim = getNonlinearSystemBase().subproblem().mesh().dimension();

  _disp_x = params.get<NonlinearVariableName>("disp_x");
  _residual_x = params.get<AuxVariableName>("residual_x");
  _diag_stiff_x = params.get<AuxVariableName>("diag_stiff_x");
  _inc_slip_x = params.get<AuxVariableName>("inc_slip_x");

  _disp_y = params.get<NonlinearVariableName>("disp_y");
  _residual_y = params.get<AuxVariableName>("residual_y");
  _diag_stiff_y = params.get<AuxVariableName>("diag_stiff_y");
  _inc_slip_y = params.get<AuxVariableName>("inc_slip_y");

  if (dim == 3)
  {
    if (!params.isParamValid("disp_z"))
      mooseError("Missing disp_z in FrictionalAugLagMulContactProblem");
    if (!params.isParamValid("residual_z"))
      mooseError("Missing residual_z in FrictionalAugLagMulContactProblem");
    if (!params.isParamValid("diag_stiff_z"))
      mooseError("Missing diag_stiff_z in FrictionalAugLagMulContactProblem");
    if (!params.isParamValid("inc_slip_z"))
      mooseError("Missing inc_slip_z in FrictionalAugLagMulContactProblem");
    _disp_z = params.get<NonlinearVariableName>("disp_z");
    _residual_z = params.get<AuxVariableName>("residual_z");
    _diag_stiff_z = params.get<AuxVariableName>("diag_stiff_z");
    _inc_slip_z = params.get<AuxVariableName>("inc_slip_z");
  }

  unsigned int num_interactions = master.size();
  if (num_interactions != slave.size())
    mooseError(
        "Sizes of master surface and slave surface lists must match in FrictionalAugLagMulContactProblem");
  if (num_interactions != friction_coefficient.size())
    mooseError(
        "Must have friction coefficient defined for every interaction in FrictionalAugLagMulContactProblem");
  if (num_interactions != slip_factor.size())
    mooseError("Must have slip factor defined for every interaction in FrictionalAugLagMulContactProblem");
  if (num_interactions != slip_too_far_factor.size())
    mooseError(
        "Must have slip too far factor defined for every interaction in FrictionalAugLagMulContactProblem");

  for (unsigned int i = 0; i < master.size(); ++i)
  {
    std::pair<int, int> ms_pair(master[i], slave[i]);
    InteractionParams ip;
    ip._friction_coefficient = friction_coefficient[i];
    ip._slip_factor = slip_factor[i];
    ip._slip_too_far_factor = slip_too_far_factor[i];

    _interaction_params[ms_pair] = ip;
  }

  _min_iters = params.get<int>("minimum_slip_iterations");
  _max_iters = params.get<int>("maximum_slip_iterations");
  _slip_updates_per_iter = params.get<int>("slip_updates_per_iteration");

  bool have_target = false;
  bool have_target_relative = false;

  if (params.isParamValid("stick_slip_tol"))
  {
    _stick_slip_tol = params.get<Real>("stick_slip_tol");
  }
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

  _contact_slip_tol_factor = params.get<Real>("contact_slip_tolerance_factor");

  if (params.isParamValid("contact_reference_residual_variables"))
    _contactRefResidVarNames =
        params.get<std::vector<std::string>>("contact_reference_residual_variables");
}

void
FrictionalAugLagMulContactProblem::initialSetup()
{
  ReferenceResidualProblem::initialSetup();

  _contactRefResidVarIndices.clear();
  for (unsigned int i = 0; i < _contactRefResidVarNames.size(); ++i)
  {
    bool foundMatch = false;
    for (unsigned int j = 0; j < _refResidVarNames.size(); ++j)
    {
      if (_contactRefResidVarNames[i] == _refResidVarNames[j])
      {
        _contactRefResidVarIndices.push_back(j);
        foundMatch = true;
        break;
      }
    }
    if (!foundMatch)
      mooseError("Could not find variable '",
                 _contactRefResidVarNames[i],
                 "' in reference_residual_variables");
  }

  //  if (_contactRefResidVarIndices.size()>0)
  //  {
  //    Moose::out << "Contact reference convergence variables:" << std::endl;
  //    for (unsigned int i=0; i<_contactRefResidVarIndices.size(); ++i)
  //    {
  //      _console << _contactRefResidVarNames[i] << std::endl;;
  //    }
  //  }
}

void
FrictionalAugLagMulContactProblem::timestepSetup()
{
  _do_slip_update = false;
  _do_lagmul_update = false;
  _num_iterations = 0;
  _num_nl_its_since_contact_update = 0;
  _refResidContact = 0.0;
  ReferenceResidualProblem::timestepSetup();
}

void
FrictionalAugLagMulContactProblem::updateContactReferenceResidual()
{
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
  ;

}


bool
FrictionalAugLagMulContactProblem::calculateSlip(const NumericVector<Number> & ghosted_solution,
                                        std::vector<SlipData> * iterative_slip)
{
  NonlinearSystemBase & nonlinear_sys = getNonlinearSystemBase();
  unsigned int dim = nonlinear_sys.subproblem().mesh().dimension();

  MooseVariable * residual_x_var = &getVariable(0, _residual_x);
  MooseVariable * residual_y_var = &getVariable(0, _residual_y);
  MooseVariable * residual_z_var = nullptr;
  MooseVariable * diag_stiff_x_var = &getVariable(0, _diag_stiff_x);
  MooseVariable * diag_stiff_y_var = &getVariable(0, _diag_stiff_y);
  MooseVariable * diag_stiff_z_var = nullptr;
  MooseVariable * inc_slip_x_var = &getVariable(0, _inc_slip_x);
  MooseVariable * inc_slip_y_var = &getVariable(0, _inc_slip_y);
  MooseVariable * inc_slip_z_var = nullptr;
  if (dim == 3)
  {
    residual_z_var = &getVariable(0, _residual_z);
    diag_stiff_z_var = &getVariable(0, _diag_stiff_z);
    inc_slip_z_var = &getVariable(0, _inc_slip_z);
  }

  bool updatedSolution = false;
  _slip_residual = 0.0;
  _it_slip_norm = 0.0;
  _inc_slip_norm = 0.0;

  if (iterative_slip)
    iterative_slip->clear();

  if (getDisplacedProblem() && _interaction_params.size() > 0)
  {
    computeResidual(ghosted_solution, getNonlinearSystemBase().RHS());

    _num_contact_nodes = 0;
    _num_slipping = 0;
    _num_slipped_too_far = 0;

    GeometricSearchData & displaced_geom_search_data = getDisplacedProblem()->geomSearchData();
    std::map<std::pair<unsigned int, unsigned int>, PenetrationLocator *> * penetration_locators =
        &displaced_geom_search_data._penetration_locators;

    AuxiliarySystem & aux_sys = getAuxiliarySystem();
    const NumericVector<Number> & aux_solution = *aux_sys.currentSolution();

    for (PLIterator plit = penetration_locators->begin(); plit != penetration_locators->end();
         ++plit)
    {
      PenetrationLocator & pen_loc = *plit->second;

      bool frictional_contact_this_interaction = false;

      std::map<std::pair<int, int>, InteractionParams>::iterator ipit;
      std::pair<int, int> ms_pair(pen_loc._master_boundary, pen_loc._slave_boundary);
      ipit = _interaction_params.find(ms_pair);
      if (ipit != _interaction_params.end())
        frictional_contact_this_interaction = true;

      if (frictional_contact_this_interaction)
      {

        InteractionParams & interaction_params = ipit->second;
        Real slip_factor = interaction_params._slip_factor;
        Real slip_too_far_factor = interaction_params._slip_too_far_factor;
        Real friction_coefficient = interaction_params._friction_coefficient;

        std::vector<dof_id_type> & slave_nodes = pen_loc._nearest_node._slave_nodes;

        for (unsigned int i = 0; i < slave_nodes.size(); i++)
        {
          dof_id_type slave_node_num = slave_nodes[i];

          if (pen_loc._penetration_info[slave_node_num])
          {
            PenetrationInfo & info = *pen_loc._penetration_info[slave_node_num];
            const Node * node = info._node;

            if (node->processor_id() == processor_id())
            {

              if (info.isCaptured())
              {
                _num_contact_nodes++;

                VectorValue<dof_id_type> residual_dofs(
                    node->dof_number(aux_sys.number(), residual_x_var->number(), 0),
                    node->dof_number(aux_sys.number(), residual_y_var->number(), 0),
                    (residual_z_var
                         ? node->dof_number(aux_sys.number(), residual_z_var->number(), 0)
                         : 0));

                VectorValue<dof_id_type> diag_stiff_dofs(
                    node->dof_number(aux_sys.number(), diag_stiff_x_var->number(), 0),
                    node->dof_number(aux_sys.number(), diag_stiff_y_var->number(), 0),
                    (diag_stiff_z_var
                         ? node->dof_number(aux_sys.number(), diag_stiff_z_var->number(), 0)
                         : 0));

                VectorValue<dof_id_type> inc_slip_dofs(
                    node->dof_number(aux_sys.number(), inc_slip_x_var->number(), 0),
                    node->dof_number(aux_sys.number(), inc_slip_y_var->number(), 0),
                    (inc_slip_z_var
                         ? node->dof_number(aux_sys.number(), inc_slip_z_var->number(), 0)
                         : 0));

                RealVectorValue res_vec;
                RealVectorValue stiff_vec;
                RealVectorValue slip_inc_vec;

                for (unsigned int i = 0; i < dim; ++i)
                {
                  res_vec(i) = aux_solution(residual_dofs(i));
                  stiff_vec(i) = aux_solution(diag_stiff_dofs(i));
                  slip_inc_vec(i) = aux_solution(inc_slip_dofs(i));
                }

                RealVectorValue slip_iterative(0.0, 0.0, 0.0);
                Real interaction_slip_residual = 0.0;
                //  _console << "inc  slip: " << slip_inc_vec << std::endl;
                //  _console << "info slip: " << info._incremental_slip << std::endl;
                //  ContactState state = calculateInteractionSlip(slip_iterative,
                //  interaction_slip_residual, info._normal, res_vec, info._incremental_slip,
                //  stiff_vec, friction_coefficient, slip_factor, slip_too_far_factor, dim);
                ContactState state = calculateInteractionSlip(slip_iterative,
                                                              interaction_slip_residual,
                                                              info._normal,
                                                              info._lagrange_multiplier,
                                                              res_vec,
                                                              slip_inc_vec,
                                                              stiff_vec,
                                                              friction_coefficient,
                                                              slip_factor,
                                                              slip_too_far_factor,
                                                              dim);
                // _console << "iter slip: " << slip_iterative << std::endl;
                _slip_residual += interaction_slip_residual * interaction_slip_residual;

                if (state == SLIPPING || state == SLIPPED_TOO_FAR)
                {
                  _num_slipping++;
                  if (state == SLIPPED_TOO_FAR)
                    _num_slipped_too_far++;
                  for (unsigned int i = 0; i < dim; ++i)
                  {
                    SlipData sd(node, i, slip_iterative(i));
                    if (iterative_slip)
                      iterative_slip->push_back(sd);
                    _it_slip_norm += slip_iterative(i) * slip_iterative(i);
                    _inc_slip_norm += (slip_inc_vec(i) + slip_iterative(i)) *
                                      (slip_inc_vec(i) + slip_iterative(i));
                  }
                }
              }
            }
          }
        }
      }
    }

    _communicator.sum(_num_contact_nodes);
    _communicator.sum(_num_slipping);
    _communicator.sum(_num_slipped_too_far);
    _communicator.sum(_slip_residual);
    _slip_residual = std::sqrt(_slip_residual);
    _communicator.sum(_it_slip_norm);
    _it_slip_norm = std::sqrt(_it_slip_norm);
    _communicator.sum(_inc_slip_norm);
    _inc_slip_norm = std::sqrt(_inc_slip_norm);
    if (_num_slipping > 0)
      updatedSolution = true;
  }

  return updatedSolution;
}


bool
FrictionalAugLagMulContactProblem::calculatePenetration(const NumericVector<Number> & ghosted_solution)
{
  //NonlinearSystemBase & nonlinear_sys = getNonlinearSystemBase();
  //unsigned int dim = nonlinear_sys.subproblem().mesh().dimension();

  _pene_residual = 0.0;

  if (getDisplacedProblem() && _interaction_params.size() > 0)
  {
    computeResidual(ghosted_solution, getNonlinearSystemBase().RHS());
    _num_contact_nodes = 0;

    GeometricSearchData & displaced_geom_search_data = getDisplacedProblem()->geomSearchData();
    std::map<std::pair<unsigned int, unsigned int>, PenetrationLocator *> * penetration_locators =
        &displaced_geom_search_data._penetration_locators;

  //  AuxiliarySystem & aux_sys = getAuxiliarySystem();
  //  const NumericVector<Number> & aux_solution = *aux_sys.currentSolution();

    for (PLIterator plit = penetration_locators->begin(); plit != penetration_locators->end();
         ++plit)
    {
      PenetrationLocator & pen_loc = *plit->second;

      bool frictional_contact_this_interaction = false;

      std::map<std::pair<int, int>, InteractionParams>::iterator ipit;
      std::pair<int, int> ms_pair(pen_loc._master_boundary, pen_loc._slave_boundary);
      ipit = _interaction_params.find(ms_pair);
      if (ipit != _interaction_params.end())
        frictional_contact_this_interaction = true;

      if (frictional_contact_this_interaction)
      {

      //  InteractionParams & interaction_params = ipit->second;

        std::vector<dof_id_type> & slave_nodes = pen_loc._nearest_node._slave_nodes;

        for (unsigned int i = 0; i < slave_nodes.size(); i++)
        {
          dof_id_type slave_node_num = slave_nodes[i];

          if (pen_loc._penetration_info[slave_node_num])
          {
            PenetrationInfo & info = *pen_loc._penetration_info[slave_node_num];
            const Node * node = info._node;

            if (node->processor_id() == processor_id())
            {

              if (info.isCaptured())
              {
                _num_contact_nodes++;
                Real interaction_pene_residual = info._normal * (info._closest_point - _mesh.nodeRef(slave_node_num));
                //  _console << "inc  slip: " << slip_inc_vec << std::endl;
                //  _console << "info slip: " << info._incremental_slip << std::endl;
                //  ContactState state = calculateInteractionSlip(slip_iterative,
                //  interaction_slip_residual, info._normal, res_vec, info._incremental_slip,
                //  stiff_vec, friction_coefficient, slip_factor, slip_too_far_factor, dim);

                // _console << "iter slip: " << slip_iterative << std::endl;
                _pene_residual += interaction_pene_residual * interaction_pene_residual;

              }
            }
          }
        }
      }
    }

    _communicator.sum(_num_contact_nodes);
    _communicator.sum(_pene_residual);
    _pene_residual = std::sqrt(_pene_residual);
    if (_pene_residual > _target_contact_residual)
      return true;
  }

  return false;
}


FrictionalAugLagMulContactProblem::ContactState
FrictionalAugLagMulContactProblem::calculateInteractionSlip(RealVectorValue & slip,
                                                   Real & slip_residual,
                                                   const RealVectorValue & normal,
                                                   const Real lagrange_multiplier,
                                                   const RealVectorValue & residual,
                                                   const RealVectorValue & incremental_slip,
                                                   const RealVectorValue & stiffness,
                                                   const Real friction_coefficient,
                                                   const Real slip_factor,
                                                   const Real slip_too_far_factor,
                                                   const int dim)
{




  ContactState state = STICKING;

  RealVectorValue normal_residual = normal * (normal * residual);
  //Real normal_force = normal_residual.norm();

  // _console << "normal=" << info._normal << std::endl;
  // _console << "normal_force=" << normal_force << std::endl;
  // _console << "residual=" << residual << std::endl;
  // _console << "stiffness=" << stiff_vec << std::endl;

  RealVectorValue tangential_force =
      normal_residual - residual; // swap sign to make the code more manageable
  Real tangential_force_magnitude = tangential_force.norm();

  Real normal_force = lagrange_multiplier + normal_residual.norm();

  Real capacity = normal_force * friction_coefficient;
  if (capacity < 0.0)
    capacity = 0.0;

  Real slip_inc = incremental_slip.norm();

  slip(0) = 0.0;
  slip(1) = 0.0;
  slip(2) = 0.0;

  if (slip_inc > _stick_slip_tol)
  {
    state = SLIPPING;
    Real slip_dot_tang_force = incremental_slip * tangential_force / slip_inc;
    if (slip_dot_tang_force < capacity)
    {
      state = SLIPPED_TOO_FAR;
      // _console << "STF slip_dot_force: " << slip_dot_tang_force << " capacity: " << capacity <<
      // std::endl;
    }
  }

  Real excess_force = tangential_force_magnitude - capacity;
  if (excess_force < 0.0)
    excess_force = 0;

  if (state == SLIPPED_TOO_FAR)
  {
    RealVectorValue slip_inc_direction = incremental_slip / slip_inc;
    Real tangential_force_in_slip_dir = slip_inc_direction * tangential_force;
    slip_residual = capacity - tangential_force_in_slip_dir;

    RealVectorValue force_from_unit_slip(0.0, 0.0, 0.0);
    for (int i = 0; i < dim; ++i)
      force_from_unit_slip(i) = stiffness(i) * slip_inc_direction(i);

    Real stiffness_slipdir =
        force_from_unit_slip *
        slip_inc_direction; // k=f resolved to slip dir because f=kd and d is a unit vector

    Real slip_distance =
        slip_too_far_factor * (capacity - tangential_force_in_slip_dir) / stiffness_slipdir;
    //    _console << "STF dist: " << slip_distance << " inc: " << slip_inc << std::endl;
    //    _console << "STF  cap: " << capacity << " tfs: " << tangential_force_in_slip_dir <<
    //    std::endl;
    if (slip_distance < slip_inc)
    {
      //      _console << "STF" << std::endl;
      slip = slip_inc_direction * -slip_distance;
    }
    else
    {
      //      _console << "STF max" << std::endl;
      slip = -incremental_slip;
    }
  }
  else if (excess_force > 0)
  {
    state = SLIPPING;
    Real tangential_force_magnitude = tangential_force.norm();
    slip_residual = excess_force;

    RealVectorValue tangential_direction = tangential_force / tangential_force_magnitude;
    RealVectorValue excess_force_vector = tangential_direction * excess_force;

    for (int i = 0; i < dim; ++i)
      slip(i) = slip_factor * excess_force_vector(i) / stiffness(i);

    // zero out component of slip in normal direction
    RealVectorValue slip_normal_dir = normal * (normal * slip);
    slip = slip - slip_normal_dir;
  }
  return state;
}



MooseNonlinearConvergenceReason
FrictionalAugLagMulContactProblem::checkNonlinearConvergence(std::string & msg,
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

  _console << "check whether to run the Augmented LM loop\n";


  ++_num_nl_its_since_contact_update;

  if ((reason > 0) ||                         // converged
      (reason == MOOSE_NONLINEAR_ITERATING && // iterating and converged within factor
       (fnorm < abstol * _contact_slip_tol_factor ||
        checkConvergenceIndividVars(
            fnorm, abstol * _contact_slip_tol_factor, rtol * _contact_slip_tol_factor, ref_resid))))
  {

    _console << "Augmentd Lagrangian Multiplier iteration " << _num_iterations << " ";
    if (_num_iterations < _max_iters)
    { // do a slip update if there is another iteration
        NonlinearSystemBase & nonlinear_sys = getNonlinearSystemBase();
        nonlinear_sys.update();
        const NumericVector<Number> *& ghosted_solution = nonlinear_sys.currentSolution();

      //  bool _augLM_repeat_step = false;
        _do_lagmul_update = false;
        _do_slip_update = false;


            if (_displaced_problem != NULL)
                  _augLM_repeat_step = nonlinear_sys.updateLagMul(true);
            else
                  _augLM_repeat_step = nonlinear_sys.updateLagMul(false);


        //   calculatePenetartion(*ghost_solution);
        //   calculateSlip(*ghosted_solution, nullptr);


      if (_augLM_repeat_step)
      {
        reason = MOOSE_NONLINEAR_ITERATING;

        _console << "iteration " << _num_iterations << " ";

       // do a slip update if there is another iteration
        _do_lagmul_update = true;
        _do_slip_update = true;
        _num_iterations++;

        if (_displaced_problem != NULL)
              _augLM_repeat_step = nonlinear_sys.updateLagMul(true);
        else
              _augLM_repeat_step = nonlinear_sys.updateLagMul(false);

        _num_iterations ++;

      // Just to calculate slip residual

       // force it to keep iterating
          _console << "Force slip update slip_resid > target: " << _slip_residual << std::endl;
        }
        else
        {          //_do_slip_update = false; //maybe we want to do this
          _console << "Not forcing slip update slip_resid <= target: " << _slip_residual
                   << std::endl;
        }
      }
      else
      {
        _console << "Max slip iterations" << std::endl;
        reason = MOOSE_DIVERGED_FUNCTION_COUNT;
      }
    }

  return reason;
}
