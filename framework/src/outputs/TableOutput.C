/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "TableOutput.h"

// MOOSE includes
#include "Conversion.h"
#include "FEProblem.h"
#include "Executioner.h"
#include "MooseApp.h"
#include "MooseVariableScalar.h"
#include "PetscSupport.h"
#include "Postprocessor.h"
#include "SystemBase.h"

// libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/string_to_enum.h"

template <>
InputParameters
validParams<TableOutput>()
{
  // Fit mode selection Enum
  MooseEnum pps_fit_mode(FormattedTable::getWidthModes());

  // Base class parameters
  InputParameters params = validParams<AdvancedOutput>();
  params += AdvancedOutput::enableOutputTypes("postprocessor scalar vector_postprocessor");

  // Option for writing vector_postprocessor time file
  params.addParam<bool>("time_data",
                        false,
                        "When true and VecptorPostprocessor data exists, write "
                        "a csv file containing the timestep and time "
                        "information.");

  // Add option for appending file on restart
  params.addParam<bool>("append_restart", false, "Append existing file on restart");

  params.addParam<bool>(
      "time_column",
      true,
      "Whether or not the 'time' column should be written for Postprocessor CSV files");

  return params;
}

TableOutput::TableOutput(const InputParameters & parameters)
  : AdvancedOutput(parameters),
    _tables_restartable(getParam<bool>("append_restart")),
    _postprocessor_table(_tables_restartable
                             ? declareRestartableData<FormattedTable>("postprocessor_table")
                             : declareRecoverableData<FormattedTable>("postprocessor_table")),
    _vector_postprocessor_time_tables(
        _tables_restartable ? declareRestartableData<std::map<std::string, FormattedTable>>(
                                  "vector_postprocessor_time_table")
                            : declareRecoverableData<std::map<std::string, FormattedTable>>(
                                  "vector_postprocessor_time_table")),
    _scalar_table(_tables_restartable ? declareRestartableData<FormattedTable>("scalar_table")
                                      : declareRecoverableData<FormattedTable>("scalar_table")),
    _all_data_table(_tables_restartable ? declareRestartableData<FormattedTable>("all_data_table")
                                        : declareRecoverableData<FormattedTable>("all_data_table")),
    _time_data(getParam<bool>("time_data")),
    _time_column(getParam<bool>("time_column"))

{
}

void
TableOutput::outputPostprocessors()
{
  // List of names of the postprocessors to output
  const std::set<std::string> & out = getPostprocessorOutput();

  // Loop through the postprocessor names and extract the values from the PostprocessorData storage
  for (const auto & out_name : out)
  {
    PostprocessorValue value = _problem_ptr->getPostprocessorValue(out_name);

    _postprocessor_table.outputTimeColumn(_time_column);
    _postprocessor_table.addData(out_name, value, time());

    _all_data_table.outputTimeColumn(_time_column);
    _all_data_table.addData(out_name, value, time());
  }
}

void
TableOutput::outputVectorPostprocessors()
{
  // List of names of the postprocessors to output
  const std::set<std::string> & out = getVectorPostprocessorOutput();

  // Loop through the postprocessor names and extract the values from the VectorPostprocessorData
  // storage
  for (const auto & vpp_name : out)
  {
    if (_problem_ptr->vectorPostprocessorHasVectors(vpp_name))
    {
      const auto & vectors = _problem_ptr->getVectorPostprocessorVectors(vpp_name);

      auto table_it = _vector_postprocessor_tables.lower_bound(vpp_name);
      if (table_it == _vector_postprocessor_tables.end() || table_it->first != vpp_name)
        table_it = _vector_postprocessor_tables.emplace_hint(table_it, vpp_name, FormattedTable());

      FormattedTable & table = table_it->second;

      table.clear();
      table.outputTimeColumn(false);

      for (const auto & vec_it : vectors)
      {
        const auto & vector = *vec_it.second.current;

        for (auto i = beginIndex(vector); i < vector.size(); ++i)
          table.addData(vec_it.first, vector[i], i);
      }

      if (_time_data)
      {
        FormattedTable & t_table = _vector_postprocessor_time_tables[vpp_name];
        t_table.addData("timestep", _t_step, _time);
      }
    }
  }
}

void
TableOutput::outputScalarVariables()
{
  // List of scalar variables
  const std::set<std::string> & out = getScalarOutput();

  // Loop through each variable
  for (const auto & out_name : out)
  {
    // Get reference to the variable (0 is for TID)
    MooseVariableScalar & scalar_var = _problem_ptr->getScalarVariable(0, out_name);

    // Make sure the value of the variable is in sync with the solution vector
    scalar_var.reinit();

    // Next we need to make sure all processors agree on the value of
    // the variable - not all processors may be able to see all
    // scalars!

    // Make a copy rather than taking a reference to the MooseArray,
    // because if a processor can't see that scalar variable's values
    // then we'll need to do our own communication of them.
    VariableValue value = scalar_var.sln();

    // libMesh *does* currently guarantee that all processors can
    // calculate all scalar DoF indices, so this is a const reference
    const std::vector<dof_id_type> & dof_indices = scalar_var.dofIndices();
    const unsigned int n = dof_indices.size();

    // In dbg mode, if we don't see a scalar we might not even have
    // its array allocated to full length yet.
    value.resize(n);

    // Finally, let's just let the owners broadcast their values.
    // There's probably lots of room to optimize this communication
    // via merging broadcasts and making them asynchronous, but this
    // code path shouldn't be hit often enough for that to matter.

    const DofMap & dof_map = scalar_var.sys().dofMap();
    for (unsigned int i = 0; i != n; ++i)
    {
      const processor_id_type pid = dof_map.dof_owner(dof_indices[i]);
      this->comm().broadcast(value[i], pid);
    }

    // If the variable has a single component, simply output the value with the name
    if (n == 1)
    {
      _scalar_table.addData(out_name, value[0], time());
      _all_data_table.addData(out_name, value[0], time());
    }

    // Multi-component variables are appended with the component index
    else
      for (unsigned int i = 0; i < n; ++i)
      {
        std::ostringstream os;
        os << out_name << "_" << i;
        _scalar_table.addData(os.str(), value[i], time());
        _all_data_table.addData(os.str(), value[i], time());
      }
  }
}
