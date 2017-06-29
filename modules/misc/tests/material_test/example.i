[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmax = 1.0
  ymax = 1.0
[]

[Variables]
[]

[AuxVariables]
  [./temp]
     order = FIRST
     family = LAGRANGE
     initial_condition = 1
  [../]
[]

[Kernels]
[]

#[AuxKernel]
#  [./TempAux]
#    type = FunctionAux
#    function = uniform    
#  [../]
#[]



[Materials]
  [./example_material]
    type = ExampleMaterial
    block = 0
    temperature = temp
    k0_param = 2.04e11
    Q_param = 0.9
    n_param = 1
    kB_param = 'kB'
  [../]
  [./constant]
    type = GenericConstantMaterial
    prop_names = 'kB'
    prop_values = '1.38e-23' #add the value of Boltzman constant
    block = 0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'None'
#  num_steps = 10
  dt = 1
  start_time = 0.0
  end_time = 10
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
