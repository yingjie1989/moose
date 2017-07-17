/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FRICTIONALAUGLAGMULCONTACTPROBLEM_H
#define FRICTIONALAUGLAGMULCONTACTPROBLEM_H

#include "ReferenceResidualProblem.h"

class FrictionalAugLagMulContactProblem;

template <>
InputParameters validParams<FrictionalAugLagMulContactProblem>();

/**
 * FEProblemBase derived class for frictional contact-specific callbacks
 */
class FrictionalAugLagMulContactProblem : public ReferenceResidualProblem
{
public:
  FrictionalAugLagMulContactProblem(const InputParameters & params);
  virtual ~FrictionalAugLagMulContactProblem() {}

  struct InteractionParams;
  struct SlipData;
  enum ContactState
  {
    STICKING,
    SLIPPING,
    SLIPPED_TOO_FAR
  };

  virtual void initialSetup();
  virtual void timestepSetup();

  void updateContactReferenceResidual();

  virtual MooseNonlinearConvergenceReason checkNonlinearConvergence(std::string & msg,
                                                                    const PetscInt it,
                                                                    const Real xnorm,
                                                                    const Real snorm,
                                                                    const Real fnorm,
                                                                    const Real rtol,
                                                                    const Real stol,
                                                                    const Real abstol,
                                                                    const PetscInt nfuncs,
                                                                    const PetscInt max_funcs,
                                                                    const Real ref_resid,
                                                                    const Real div_threshold);


protected:
  std::map<std::pair<int, int>, InteractionParams> _interaction_params;
  NonlinearVariableName _disp_x;
  NonlinearVariableName _disp_y;
  NonlinearVariableName _disp_z;
  AuxVariableName _residual_x;
  AuxVariableName _residual_y;
  AuxVariableName _residual_z;
  AuxVariableName _diag_stiff_x;
  AuxVariableName _diag_stiff_y;
  AuxVariableName _diag_stiff_z;
  AuxVariableName _inc_slip_x;
  AuxVariableName _inc_slip_y;
  AuxVariableName _inc_slip_z;

  std::vector<std::string> _contactRefResidVarNames;
  std::vector<unsigned int> _contactRefResidVarIndices;
  Real _refResidContact;

  Real _slip_residual;
  Real _pene_residual;
  bool _do_slip_update;
  bool _do_lagmul_update;
  int _num_iterations;
  int _min_iters;
  int _max_iters;
  int _slip_updates_per_iter;
  Real _stick_slip_tol;
  Real _target_contact_residual;
  Real _target_relative_contact_residual;
  Real _contact_slip_tol_factor;
  int _num_nl_its_since_contact_update;
  int _num_contact_nodes;
  int _num_slipping;
  int _num_slipped_too_far;
  Real _inc_slip_norm;
  Real _it_slip_norm;

  /// Convenient typedef for frequently used iterator
  typedef std::map<std::pair<unsigned int, unsigned int>, PenetrationLocator *>::iterator
      PLIterator;
};

struct FrictionalAugLagMulContactProblem::InteractionParams
{
  Real _friction_coefficient;
  Real _slip_factor;
  Real _slip_too_far_factor;
};

struct FrictionalAugLagMulContactProblem::SlipData
{
  SlipData(const Node * node, unsigned int dof, Real slip);
  SlipData(const SlipData & sd);
  ~SlipData();

  const Node * _node;
  unsigned int _dof;
  Real _slip;
};

#endif /* FrictionalAugLagMulContactProblem_H */
