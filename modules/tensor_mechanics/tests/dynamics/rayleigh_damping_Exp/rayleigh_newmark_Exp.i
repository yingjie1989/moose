# test for rayleigh damping implemented using newmark time integration

# the test is for an 1d bar element of  unit length fixed on one end
# with a ramped pressure boundary condition applied to the other end.
# zeta and eta correspond to the stiffness and mass proportional rayleigh damping
# beta and gamma are newmark time integration parameters
# the equation of motion in terms of matrices is:
#
# m*accel + eta*m*vel + zeta*k*vel + k*disp = p*area
#
# here m is the mass matrix, k is the stiffness matrix, p is the applied pressure
#
# this equation is equivalent to:
#
# density*accel + eta*density*vel + zeta*d/dt(div stress) + div stress = p
#
# the first two terms on the left are evaluated using the inertial force kernel
# the next two terms on the left involving zeta are evaluated using the
# dynamicstressdivergencetensors kernel
# the residual due to pressure is evaluated using pressure boundary condition
#
# the system will come to steady state slowly after the pressure becomes constant.
# the store_stress_old flag in the computestressbase material model needs to be
# turned on to store stress old. in this example, this flag is turned on using
# the child class computelinearelasticstress.

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmin = 0.0
  xmax = 0.1
  ymin = 0.0
  ymax = 1.0
  zmin = 0.0
  zmax = 0.1
  displacements = 'disp_x disp_y disp_z'
[]


[Variables]
  [./accel_x]
  [../]
  [./accel_y]
  [../]
  [./accel_z]
  [../]
[]

[AuxVariables]
  [./vel_x]
  [../]
  [./disp_x]
  [../]
  [./vel_y]
  [../]
  [./disp_y]
  [../]
  [./vel_z]
  [../]
  [./disp_z]
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
   [./solid_x]
     type = StressDivergenceExpTensors
     variable = accel_x
     displacements = 'disp_x disp_y disp_z'
     component = 0
 [../]
  [./solid_y]
    type = StressDivergenceExpTensors
     variable = accel_y
     displacements = 'disp_x disp_y disp_z'
    component = 1
  [../]
  [./solid_z]
     type = StressDivergenceExpTensors
     variable = accel_z
     displacements = 'disp_x disp_y disp_z'
    component = 2
  [../]
  [./inertia_x]
    type = MassLumpedReaction
#  type = CoefMassReaction
    variable = accel_x
#    coefficient = 7750
  [../]
  [./inertia_y]
    type = MassLumpedReaction
#    type = CoefMassReaction
    variable = accel_y
#    coefficient = 7750
  [../]
  [./inertia_z]
    type = MassLumpedReaction
#    type = CoefMassReaction
    variable = accel_z
#    coefficient = 7750
  [../]
[]

[AuxKernels]
   [./disp_x]
    type = NewmarkDispAux
    variable = disp_x
    acceleration = accel_x
    velocity = vel_x
    beta = 0.5
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./disp_y]
    type = NewmarkDispAux
    variable = disp_y
    acceleration = accel_y
    velocity = vel_y
    beta = 0.5
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./disp_z]
    type = NewmarkDispAux
    variable = disp_z
    acceleration = accel_z
    velocity = vel_z
    beta = 0.5
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./vel_z]
    type = NewmarkVelAux
    variable = vel_z
    acceleration = accel_z
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./stress_yy]
   type = RankTwoAux
   rank_two_tensor = stress
   variable = stress_yy
   index_i = 1
    index_j = 1
    execute_on = timestep_end
 [../]
 [./strain_yy]
   type = RankTwoAux
   rank_two_tensor = total_strain
    variable = strain_yy
    index_i = 1
   index_j = 1
   execute_on = timestep_end
 [../]
[]


[BCs]
  [./top_y]
    type = DirichletBC
    variable = accel_y
    boundary = top
    value=0.0
  [../]
 [./top_x]
   type = DirichletBC
   variable = accel_x
   boundary = top
   value=0.0
 [../]
 [./top_z]
    type = DirichletBC
   variable = accel_z
       boundary = top
    value=0.0
  [../]
  [./bottom_x]
    type = DirichletBC
    variable = accel_x
    boundary = bottom
    value=0.0
  [../]
  [./bottom_z]
    type = DirichletBC
    variable = accel_z
    boundary = bottom
    value=0.0
  [../]
  [./bottom_y]
    type = NeumannBC
    variable = accel_y
    boundary = bottom
    value = 1e8
  [../]
#  [./Pressure]
#    [./Side1]
#      boundary = bottom
#      function = pressure
#      disp_x = disp_x
#      disp_y = disp_y
#      disp_z = disp_z
#      factor = 1
#    [../]
#  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    fill_method = symmetric_isotropic
    C_ijkl = '210e9 0'
  [../]

  [./strain]
    type = ComputeSmallStrain
    block = 0
#    ifOld = true
    displacements = 'disp_x disp_y disp_z'
  [../]

  [./stress]
    type = ComputeLinearElasticStress
    store_stress_old = True
    block = 0
  [../]

  [./density]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'density'
    prop_values = '7750'
  [../]
[]

[Executioner]
  type = Transient

 petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
 petsc_options_value = 'lu     mumps'

  start_time = 0
  end_time = 2
  dt = 1e-4
  dtmin = 1e-5
  nl_abs_tol = 1e-9
#  num_steps = 10
[]


[Functions]
  [./pressure]
    type = PiecewiseLinear
    x = '0.0 0.1 0.2 1.0 2.0 5.0'
    y = '0.0 0.1 0.2 1.0 1.0 1.0'
    scale_factor = 1e9
  [../]
[]

#[Postprocessors]
# [./_dt]
#   type = TimestepSize
#  [../]
#  [./disp]
#    type = NodalMaxValue
#    variable = disp_y
#    boundary = bottom
#  [../]
#  [./vel]
#    type = NodalMaxValue
#    variable = vel_y
#    boundary = bottom
#  [../]
#  [./accel]
#    type = NodalMaxValue
#    variable = accel_y
#    boundary = bottom
#  [../]
#  [./stress_yy]
#    type = ElementAverageValue
#    variable = stress_yy
#  [../]
#  [./strain_yy]
#    type = ElementAverageValue
#    variable = strain_yy
#  [../]
#[]

[Outputs]
  exodus = true
  interval = 100
  print_perf_log = true
[]
