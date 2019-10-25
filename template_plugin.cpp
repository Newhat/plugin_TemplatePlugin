/*
 * Copyright (c) 2011-2016:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "bridge/util.h"

// if possible, replace this with util_domain_dependent.h or
// util_algebra_dependent.h to speed up compilation time
#include "bridge/util_domain_algebra_dependent.h"

// begin test 10/24/2019 ---------------------------------------------------------------------------

//header from grid_function_util.h-----------------------------------------------------------------
#include <vector>
#include <string>
#include <cmath>  // for isinf, isnan
#include <boost/function.hpp>
#include "common/util/file_util.h"
#include "lib_algebra/cpu_algebra/sparsematrix_print.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/debug_writer.h"
#include "lib_algebra/operator/vector_writer.h"
#include "lib_algebra/common/matrixio/matrix_io_mtx.h"
#include "lib_algebra/common/connection_viewer_output.h"
#include "lib_algebra/common/csv_gnuplot_output.h"
#include "lib_grid/algorithms/debug_util.h"  // for ElementDebugInfo
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "lib_disc/io/vtkoutput.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/dof_position_util.h"



using namespace std;
using namespace ug::bridge;

namespace ug{
namespace TemplatePlugin{


//	'TemplateSampleClass' and 'TemplateSampleFunction' serve as examples. Please
//	delete them if you use this file as a template for your own plugin.
//
//	If you want to test the registry with those classes, please uncomment the
//	commented sections in 'Functionality::DomainAlgebra' and 'Functionality::Common'.

///	a sample class that is used to show how class-groups can be registered.
template <class TDomain, class TAlgebra>
class TemplateSampleClass {
	public:
		TemplateSampleClass ()		{}
		void print_hello () const	{UG_LOG("hello\n");}
};

///	a sample function that us used to show how a simple function can be registered
void TemplateSampleFunction () {
	UG_LOG("TemplateSampleFunction executed.\n");
}

// for test
// from connection_viewer_input.h
// with additional checks
template<typename vector_type>
bool ReadVector(std::string filename, vector_type &vec,int dim)
{
    Progress p;
	std::cout << " Reading std::vector from " <<  filename << "... ";
	std::fstream matfile(filename.c_str(), std::ios::in);
	if(matfile.is_open() == false) { std::cout << "failed.\n"; return false; }

	int version=-1, dimension=-1, gridsize;

	matfile >> version;
	matfile >> dimension;
	matfile >> gridsize;

	assert(version == 1);
	assert(dimension == dim);
	// todo check positions and not just size
	assert(gridsize == (int)vec.size());


	PROGRESS_START(prog, gridsize*2, "ReadVector " << dimension << "d from " << filename << " , " << gridsize << " x " << gridsize);
	for(int i=0; i<gridsize; i++)
	{
		if(i%100) { PROGRESS_UPDATE(prog, i); }
		if(matfile.eof())
		{
			std::cout << " failed.\n";
			assert(0);
			return false;
		}
		double x, y, z;
		matfile >> x >> y;
		if(dimension==3) matfile >> z;
	}

	int printStringsInWindow;
	matfile >> printStringsInWindow;

	// vec.resize(gridsize);
	bool bEOF = matfile.eof();
	while(!bEOF)
	{
		int from, to; double value;
		char c = matfile.peek();
		if(c == -1 || c == 'c' || c == 'v' || matfile.eof())
			break;

		matfile >> from >> to >> value;
		assert(from == to);
		vec[from] = value;
		if(from%100) { PROGRESS_UPDATE(prog, from); }
		bEOF = matfile.eof();
	}
	return true;
}

// load vector that has been saved in connection viewer format and write it
// into grid function
template<typename TGridFunction>
void LoadVector22(TGridFunction& u,const char* filename){
	PROFILE_FUNC();
	typename TGridFunction::algebra_type::vector_type b;
	b.resize(u.num_indices());
	ReadVector(filename,b,TGridFunction::dim);
	u.assign(b);
}
//test end


/** 
 *  \defgroup plugin_template Plugin Template
 *  \ingroup plugins_experimental
 *  This is a template for new plugins.
 *  \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	//	typedef
		//static const int dim = TDomain::dim;
		//typedef typename TAlgebra::vector_type vector_type;
		//typedef typename TAlgebra::matrix_type matrix_type;
		typedef GridFunction<TDomain, TAlgebra> function_type;

//	The code below illustrates how a template-dependend class
//	can be registered as a class-group.
	/*
	 {
	 	typedef TemplateSampleClass<TDomain, TAlgebra> T;
	 	string name = string("TemplateSampleClass").append(suffix);
	 	reg.add_class_<T>(name, grp)
	 		.add_constructor()
	 		.add_method("print_hello", &T::print_hello, "", "", "prints hello")
	 		.set_construct_as_smart_pointer(true);
	 	reg.add_class_to_group(name, "TemplateSampleClass", tag);
	 }
	 */
		//	WriteGridToVTK
			{
				reg.add_function("LoadVector22",
								 &LoadVector22<function_type>, grp,
									"", "GridFunction#Filename|save-dialog|endings=[\"vtk\"];description=\"VTK-Files\"",
									"Saves GridFunction to *.vtk file", "No help");
			}
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */

static void Common(Registry& reg, string grp)
{
//	The code below shows how a simple function can be registered
	 reg.add_function("TemplateSampleFunction", &TemplateSampleFunction, grp,
	 				 "", "", "Prints a short message");
}

}; // end Functionality

// end group plugin_template
/// \}

} // end namespace TemplatePlugin


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_TemplatePlugin(Registry* reg, string grp)
{
	grp.append("TemplatePlugin");
	typedef TemplatePlugin::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
