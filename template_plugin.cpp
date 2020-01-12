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

#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

#include <fstream>

// test
#include "common/common.h"
#include "common/math/ugmath.h"
#include "common/util/smart_pointer.h"
#include "lib_disc/dof_manager/dof_distribution.h"

//#include "dof_position_util.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/domain.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/domain_traits.h"

#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/reference_element/reference_element_util.h"
// end test
#include <limits>
#include <map>


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

// for test 10/25/2019
// from connection_viewer_input.h
// with additional checks
template<typename TGridFunction>
bool LoadVector22(TGridFunction& vec,const char* filename)
{
	//Progress p;
	std::cout << " Reading std::vector from " <<  filename << "... ";
	std::ifstream matfile;
	matfile.open(filename);
	if(matfile.is_open() == false) { std::cout << "failed.\n"; return false; }

	int version=-1, dimension=-1, gridsize;

	matfile >> version;
	matfile >> dimension;
	matfile >> gridsize;

	assert(version == 1);
	// todo check positions and not just size
	assert(gridsize == (int)vec.size());

	PROGRESS_START(prog, gridsize*2, "ReadVector "
			<< dimension << "d from " << filename << " , " << gridsize << " x " << gridsize);
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
	std::cout<<"\n Me ---------------"<<std::endl;
	return true;
}

template<typename TGridFunction, typename TAlgebra>
void Assign22(TGridFunction& vec,TAlgebra &b){
	//PROFILE_FUNC();
	int gridsize, size2;
	gridsize = (int)vec.size();
	size2 = (int)b.size();
	cout<<"u: "<<endl;
	for(int i=0; i<gridsize; i++)
	{
		//vec[i] = b[i];
		//cout<<vec[i]<<endl;
	}
	cout<<"b: "<<endl;
	for(int i=0; i<size2; i++)
	{
		cout<<b[i]<<endl;
	}
}
//test end
template<typename TAlgebra>
void disp(TAlgebra &vec){
	//this function is to disply vector
	int gridsize;
	gridsize = (int)vec.size();
	//cout<<"This vector is: "<<endl;
	for(int i=0; i<gridsize; i++)
	{
		cout<<vec[i]<<endl;
	}
}

template<typename TGridFunction>
void Restriction3dto1d(TGridFunction &b3to1,TGridFunction &b11,TGridFunction &b31to1){
	//this function is to disply vector
	int gridsize;
	gridsize = (int)b11.size();
	ug::DenseVector<ug::FixedArray1<double, 2> > temp1, temp2, temp3;
	for(int i=2; i<gridsize; i++)
	{
		temp1 = b3to1[i];
		temp2 = b11[i];
		temp3 = b31to1[i];
		b3to1[i] = temp1[0]/temp3[0]*temp2[0];	
	}
}
// test --------------10/30/2019

template <typename TDomain, typename TGridFunction, typename TAlgebra>
void ExtractPositionsV3dto1d(TGridFunction &u,TAlgebra &b, SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd,
		std::vector<MathVector<TDomain::dim> >& vPos)
{
	//	number of total dofs
	int ns, nr = dd->num_indices();
	ug::DenseVector<ug::FixedArray1<double, 2> > temp;
	//	resize positions
	vPos.resize(nr);
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	std::vector<size_t> ind;
	//	get iterators for subsets
	int MarkerSize = (int)b.size();
	int marker;
	for(marker = 3;marker<=MarkerSize;marker++)
	{	//cout<<"marker is: ";
		//cout<<marker<<endl;
		iter = dd->begin<Vertex>(marker, SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(marker, SurfaceView::ALL);
		temp[0] = 0.0;
		temp[1] = 0.0;
		ns = 0;
		//	loop vertices with certian marker
		for(;iter != iterEnd; ++iter)
		{
			//	get vertex
			Vertex* v = *iter;		
			//	load indices associated with vertex
			dd->inner_algebra_indices(v, ind);
			//	write position
			//cout<<ind.size()<<endl;
			for(size_t i = 0; i < ind.size(); ++i)
			{
				const size_t index = ind[i];
				vPos[index] = aaPos[v];
				//cout<<vPos[index]<<endl;
				//cout<<u[index]<<endl;
				//u[index] = 0;
				temp=temp+u[index];
				//cout<<temp<<endl;
				ns += 1;
			}
		}
		//cout<<"double vector"<<endl;
		b[marker-1] = temp[0];
		//cout<<temp[0]/((double)(1))<<endl;
	}
}

template <typename TDomain, typename TGridFunction, typename TAlgebra>
void ExtractPositionsV1dto3d(TGridFunction &u,TAlgebra &b, SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd,
		std::vector<MathVector<TDomain::dim> >& vPos)
{
	//	number of total dofs
	int ns, nr = dd->num_indices();
	ug::DenseVector<ug::FixedArray1<double, 2> > temp;
	//	resize positions
	vPos.resize(nr);
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	std::vector<size_t> ind;
	//	get iterators for subsets
	int MarkerSize = (int)b.size();
	int marker;
	for(marker = 3;marker<=MarkerSize;marker++)
	{	//cout<<"marker is: ";
		//cout<<marker<<endl;
		iter = dd->begin<Vertex>(marker, SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(marker, SurfaceView::ALL);
		temp[0] = 0.0;
		temp[1] = 0.0;
		ns = 0;
		//	loop vertices with certian marker
		for(;iter != iterEnd; ++iter)
		{
			//	get vertex
			Vertex* v = *iter;		
			//	load indices associated with vertex
			dd->inner_algebra_indices(v, ind);
			//	write position
			//cout<<ind.size()<<endl;
			for(size_t i = 0; i < ind.size(); ++i)
			{
				const size_t index = ind[i];
				vPos[index] = aaPos[v];
				//cout<<vPos[index]<<endl;
				//cout<<u[index]<<endl;
				u[index] = b[marker-1];
				//temp=temp+u[index];
				//cout<<temp<<endl;
				ns += 1;
			}
		}
		//cout<<"double vector"<<endl;
		//b[marker-1] = temp[0];
		//cout<<temp[0]/((double)(1))<<endl;
	}
}


template<class TGridFunction,typename TAlgebra>
void ExtractPositions3dto1d(TGridFunction &u,TAlgebra &b)
{

	// 	get positions of vertices
	const static int dim = TGridFunction::domain_type::dim;
	std::vector<MathVector<dim> > vPos;
	ExtractPositionsV3dto1d(u,b,u.domain(),u.dof_distribution(), vPos);
	//	cout<<vPos[0]<<endl;
	//	cout<<vPos[0][0]<<endl;
}

template<class TGridFunction,typename TAlgebra>
void ExtractPositions1dto3d(TAlgebra &b,TGridFunction &u)
{

	// 	get positions of vertices
	const static int dim = TGridFunction::domain_type::dim;
	std::vector<MathVector<dim> > vPos;
	ExtractPositionsV1dto3d(u,b,u.domain(),u.dof_distribution(), vPos);
	//	cout<<vPos[0]<<endl;
	//	cout<<vPos[0][0]<<endl;
}

// test 2020 ------------------------------------------------------------------------

// 1/7/2020 -----------------------------------------------1214
template <typename TDomain>
Vertex* Vnext(Vertex* v,SmartPtr<TDomain> domain, char lr)
{
	double temp = std::numeric_limits<double>::max();
	double temp1 = 0.0;
	int Nt,Nv;
	Vertex* vt = NULL;
	Vertex* vn = NULL;
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	MultiGrid::edge_traits::secure_container edges;
	domain->grid()->associated_elements(edges, v);

	//cout<<"how many edges associated with v"<<endl;
	//cout<<edges.size()<<endl;

	for(size_t i = 0; i < edges.size(); ++i)
	{
		Edge* e = edges[i];
		//cout<<"-------------------"<<endl;
		for(size_t j=0;j<2;j++)
		{
			if(e->vertex(j)!=v)
			{
				vt = e->vertex(j);
				Nt = Dsh->get_subset_index(vt);
				Nv = Dsh->get_subset_index(v) ;
				if(lr=='l')
				{
					if( Nt<Nv && Nt>1)
					{
						//cout<<i<<" "<<aaPos[vt]<<endl;
						//cout<<"marker = "<<Dsh->get_subset_index(vt)<<endl;
						temp1 = sqrt(   pow(aaPos[v][0]-aaPos[vt][0],2)
								+pow(aaPos[v][1]-aaPos[vt][1],2)
								+pow(aaPos[v][2]-aaPos[vt][2],2) );
						//cout<<"the dist = "<<temp1<<endl;
						if(temp1<temp)
						{
							temp = temp1;
							vn = vt;
						}
					}
				}
				else
				{
					if( Nt>Nv && Nt>1 )
					{
						//cout<<i<<" "<<aaPos[vt]<<endl;
						//cout<<"marker = "<<Dsh->get_subset_index(vt)<<endl;
						temp1 = sqrt(   pow(aaPos[v][0]-aaPos[vt][0],2)
								+pow(aaPos[v][1]-aaPos[vt][1],2)
								+pow(aaPos[v][2]-aaPos[vt][2],2) );
						//cout<<"the dist = "<<temp1<<endl;
						if(temp1<temp)
						{
							temp = temp1;
							vn = vt;
						}
					}
				}

			}
		}
	}
	//cout<<"min dist = "<<temp<<endl;
	//cout<<"the vn = "<<aaPos[vn]<<endl;
	//cout<<"the v  = "<<aaPos[v]<<endl;
	return vn;

}
// 1/7/2020 -----------------------------------------------1214
template <typename TDomain,typename TAlgebra>
ug::DenseMatrix<ug::FixedArray2<double, 2, 2> >
CoarseFromFine(Vertex* V,Vertex* Vn,SmartPtr<DoFDistribution> dd,
		SmartPtr<TDomain> domain,
		TAlgebra &N,
		map<Vertex*, size_t> Mmap,
		vector<size_t> markers)
		{
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	ug::DenseMatrix<ug::FixedArray2<double, 2, 2> > temp;
	Vertex* v11;
	Vertex* v12;
	vector<Vertex*> Vrs;
	vector<Vertex*> Vls;
	vector<Vertex*> Vnrs;
	vector<Vertex*> Vnls;
	vector<Vertex*> Vs;
	vector<Vertex*> Vns;
	map<Vertex*, double> CV,CVn;

	//cout<<"marker = "<<Dsh->get_subset_index(V)<<endl;

	std::vector<size_t> ind, ind1,ind2;
	Dsh->get_subset_index(V);
	//Dsh->get_subset_index(Vrs[0]);
	aaPos[V];
	//1, basis for V point -------------------------------
	v12 = V;
	for(size_t i=markers[Mmap[V]];i<markers[Mmap[V]+1];i++)
	{
		v11 = Vnext(v12,domain,'r');
		Vrs.push_back(v11);
		v12 = v11;
	}
	v12 = V;
	for(size_t i=markers[Mmap[V]-1];i<markers[Mmap[V]];i++)
	{
		v11 = Vnext(v12,domain,'l');
		Vls.push_back(v11);
		v12 = v11;
	}
	//2, basis for Vn point -------------------------------
	v12 = Vn;
	for(size_t i=markers[Mmap[Vn]];i<markers[Mmap[Vn]+1];i++)
	{
		v11 = Vnext(v12,domain,'r');
		Vnrs.push_back(v11);
		v12 = v11;
		//cout<<aaPos[v11]<<endl;
	}
	v12 = Vn;
	for(size_t i=markers[Mmap[Vn]-1];i<markers[Mmap[Vn]];i++)
	{
		v11 = Vnext(v12,domain,'l');
		Vnls.push_back(v11);
		v12 = v11;
		//cout<<aaPos[v11]<<endl;
	}
	// the coefficients of V's basis functions ---------------
	CV[V] = 1.0;
	Vs.push_back(V);
	//cout<<"V: "<<aaPos[V]<<endl;
	//cout<<"Vrs.size(): "<<Vrs.size()<<endl;
	//cout<<"Vrs.back(): "<<aaPos[Vrs.back()]<<endl;
	for(size_t j=0;j<Vrs.size();j++)
	{
		CV[Vrs[j]] =((double)Vrs.size()-1.0-(double)j )/(double)Vrs.size();
		Vs.push_back(Vrs[j]);
		//cout<<aaPos[Vrs[j]]<<";  "<<CV[Vrs[j]]<<endl;
	}
	//cout<<"Vls.size(): "<<Vls.size()<<endl;
	//cout<<"Vls.back(): "<<aaPos[Vls.back()]<<endl;
	for(size_t j=0;j<Vls.size();j++)
	{
		CV[Vls[j]] =((double)Vls.size()-1.0-(double)j )/(double)Vls.size();
		Vs.push_back(Vls[j]);
		//cout<<aaPos[Vls[j]]<<";  "<<CV[Vls[j]]<<endl;
	}
	// the coefficients of Vn's basis functions ---------------
	CVn[Vn] = 1.0;
	Vns.push_back(Vn);
	for(size_t j=0;j<Vnrs.size();j++)
	{
		CVn[Vnrs[j]] =((double)Vnrs.size()-1.0-(double)j )/(double)Vnrs.size();
		Vns.push_back(Vnrs[j]);
	}
	for(size_t j=0;j<Vnls.size();j++)
	{
		CVn[Vnls[j]] =((double)Vnls.size()-1.0-(double)j )/(double)Vnls.size();
		Vns.push_back(Vnls[j]);
	}
	// Matrix -----------------------------------------------------------------
	//cout<<N(ind2[0],ind2[0])<<endl;
	//cout<<"Vs.size(): "<<Vs.size()<<endl;
	//cout<<"Vns.size(): "<<Vns.size()<<endl;
	temp = 0.0;
	//N(0,0) = 0.0;
	for(size_t k=0;k<Vs.size();k++)
	{
		dd->inner_algebra_indices(Vs[k],  ind1);
		for(size_t j=0;j<Vns.size();j++)
		{
			dd->inner_algebra_indices(Vns[j], ind2);
			temp+=CV[Vs[k]]*CVn[Vns[j]]*N(ind1[0],ind2[0]);
		}
	}

	return temp;
		}
// 1/9/2020 -----------------------------------------------0836

template <typename TDomain, typename TGridFunction, typename TAlgebra>
void Matrix3dFineToCoarse(TGridFunction &u,TAlgebra &M,TAlgebra &N,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd)
{
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//  get marker accessor
	typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	std::vector<size_t> ind;
	// 12301------------
	map<Vertex*, size_t> umap, Mmap;
	//	The unkowns collected by chosen markers
	size_t MarkerSize = Dsh->num_subsets();
	cout<<"marker max size is "<<MarkerSize<<endl;
	size_t num = MarkerSize-3;
	size_t Ms = 10;
	vector<size_t> markers;
	markers.push_back(3);
	for(size_t i=1; i<=(num-num % Ms)/Ms;i++)
	{
		markers.push_back(2+Ms*i);
		//cout << i <<"  "<< 2+Ms*i<<endl;
	}
	markers.push_back(MarkerSize-1);
	// the map from Vertex* to index of coarse matrix
	size_t jc = 0;
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		iter = dd->begin<Vertex>(markers[mk], SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(markers[mk], SurfaceView::ALL);
		for(;iter != iterEnd; ++iter)
		{
			//	get vertex
			//Vertex* v = *iter;
			Mmap[*iter] = mk;
			umap[*iter] = jc;
			jc+=1;
			//cout<<umap[*iter]<<endl;
		}
	}
	cout<<"jc is: "<<jc<<endl;
	M.resize_and_clear(jc, jc);
	//cout<<markers.size()<<endl;
	//cout<<markers[0]<<endl;
	size_t marker;
	//int mk = 1;
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		marker = markers[mk];
		iter = dd->begin<Vertex>(marker, SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(marker, SurfaceView::ALL);
		//iter = dd->begin<Vertex>(SurfaceView::ALL);
		//iterEnd = dd->end<Vertex>(SurfaceView::ALL);
		for(;iter != iterEnd; ++iter)
		{

			//	get vertex---------------------------------------------12301
			Vertex* V = *iter;
			Vertex* v11;
			Vertex* v12;
			//vector<Vertex*> Vrs;
			//vector<Vertex*> Vls;
			vector<Vertex*> Vos;
			//cout<<aaPos[V]<<endl;
			aaPos[V];
			// dd->inner_algebra_indices(V, ind);
			// cout<<" :"<<u[ind[0]]<<endl;
			// cout<<"coarse index = "<<umap[V]<<endl;
			// 1/7/2020 -----------------------------------------------1214
			Vos.push_back(V);

			v12 = V;
			for(size_t i=markers[mk];i<markers[mk+1];i++)
			{
				v11 = Vnext(v12,domain,'r');
				//Vrs.push_back(v11);
				v12 = v11;
				//cout<<"marker = "<<Dsh->get_subset_index(v11)<<endl;
				//dd->inner_algebra_indices(v11, ind);
				//cout<<"index = "<<ind[0]<<endl;
				//cout<<aaPos[v11]<<endl;
			}
			if((size_t)Dsh->get_subset_index(v11)<markers.back())
			{
				Vos.push_back(v11);
				//cout<<aaPos[Vrs.back()]<<endl;
			}
			v12 = V;
			for(size_t i=markers[mk-1];i<markers[mk];i++)
			{
				v11 = Vnext(v12,domain,'l');
				//Vls.push_back(v11);
				v12 = v11;
				//cout<<"marker = "<<Dsh->get_subset_index(v11)<<endl;
				//dd->inner_algebra_indices(v11, ind);
				//cout<<"index = "<<ind[0]<<endl;
				//cout<<aaPos[v11]<<endl;
			}
			if((size_t)Dsh->get_subset_index(v11)>markers[0])
			{
				Vos.push_back(v11);
				//cout<<aaPos[Vls.back()]<<endl;
			}

			// 1/7/2020 -----------------------------------------------1214
			//  collected edges
			MultiGrid::edge_traits::secure_container edges;
			domain->grid()->associated_elements(edges, V);
			//cout<<"how many edges associated with V"<<endl;
			//cout<<edges.size()<<endl;
			for(size_t i = 0; i < edges.size(); ++i)
			{
				//vector<Vertex*> Vnrs;
				//vector<Vertex*> Vnls;
				Vertex* Vn;
				//cout<<i<<endl;
				Edge* e = edges[i];
				//cout<<"-------------------"<<endl;
				for(size_t j=0;j<2;j++)
				{
					if((e->vertex(j)!=V)
							&&(Dsh->get_subset_index(e->vertex(j))==Dsh->get_subset_index(V)))
					{
						Vn = e->vertex(j);
						Vos.push_back(Vn);
						//dd->inner_algebra_indices(e->vertex(j), ind);
						//cout<<"nby verts: "<<aaPos[e->vertex(j)]<<endl;
						//cout<<"marker = "<<markers[Mmap[Vn]]<<endl;
						//cout<<"index = "<<ind[0]<<endl;
						v12 = Vn;
						for(size_t i=markers[mk];i<markers[mk+1];i++)
						{
							v11 = Vnext(v12,domain,'r');
							//Vnrs.push_back(v11);
							v12 = v11;
							//cout<<"marker = "<<Dsh->get_subset_index(v11)<<endl;
							//dd->inner_algebra_indices(v11, ind);
							//cout<<"index = "<<ind[0]<<endl;
							//cout<<aaPos[v11]<<endl;
						}
						if((size_t)Dsh->get_subset_index(v11)<markers.back())
						{
							Vos.push_back(v11);
							//cout<<"right: "<<aaPos[v11]<<endl;
						}
						//cout<<"nearby verts: "<<aaPos[e->vertex(j)]<<endl;
						v12 = Vn;
						for(size_t i=markers[mk-1];i<markers[mk];i++)
						{
							v11 = Vnext(v12,domain,'l');
							//Vnls.push_back(v11);
							v12 = v11;
							//cout<<"marker = "<<Dsh->get_subset_index(v11)<<endl;
							//dd->inner_algebra_indices(v11, ind);
							//cout<<"index = "<<ind[0]<<endl;
							//cout<<aaPos[v11]<<endl;
						}
						if((size_t)Dsh->get_subset_index(v11)>markers[0])
						{
							Vos.push_back(v11);
							//cout<<"left: "<<aaPos[v11]<<endl;
						}
						//cout<<"center v : "<<aaPos[v]<<endl;
						//cout<<"nby verts: "<<aaPos[Vn]<<endl;
					}
				}
			}
			//--------------------------------------------------------12301
			//assemble coarse matrix from fine matrix
			ug::DenseMatrix<ug::FixedArray2<double, 2, 2> > temp;
			for(size_t j=0;j<Vos.size();j++)
			{
				temp = CoarseFromFine(V,Vos[j],dd,domain,N,Mmap,markers);
				//cout<<"temp :"<<temp<<endl;
				M(umap[V],umap[Vos[j]]) = temp(0,0);
			}
			//cout<<M(umap[V],umap[Vos[j]])<<endl;
		}

	}
}

template <typename TDomain, typename TGridFunction, typename TAlgebra>
void Vector3dFineToCoarse(TGridFunction &u,TAlgebra &vec,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd)
{
	//ug::DenseVector<ug::FixedArray1<double, 2> > temp;
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//  get marker accessor
	typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	std::vector<size_t> ind;
	// 12301------------
	map<Vertex*, size_t> umap, Mmap;
	//	The unkowns collected by chosen markers
	size_t MarkerSize = Dsh->num_subsets();
	cout<<"marker max size is "<<MarkerSize<<endl;
	size_t num = MarkerSize-3;
	size_t Ms = 10;
	vector<size_t> markers;
	markers.push_back(3);
	for(size_t i=1; i<=(num-num % Ms)/Ms;i++)
	{
		markers.push_back(2+Ms*i);
		//cout << i <<"  "<< 2+Ms*i<<endl;
	}
	markers.push_back(MarkerSize-1);
	// the map from Vertex* to index of coarse matrix
	size_t jc = 0;
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		iter = dd->begin<Vertex>(markers[mk], SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(markers[mk], SurfaceView::ALL);
		for(;iter != iterEnd; ++iter)
		{
			//	get vertex
			//Vertex* v = *iter;
			Mmap[*iter] = mk;
			umap[*iter] = jc;
			jc+=1;
			//cout<<umap[*iter]<<endl;
		}
	}
	cout<<"jc is: "<<jc<<endl;
	vec.resize(jc);
	//cout<<markers.size()<<endl;
	//cout<<markers[0]<<endl;
	size_t marker;
	//int mk = 1;
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		marker = markers[mk];
		iter = dd->begin<Vertex>(marker, SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(marker, SurfaceView::ALL);
		//iter = dd->begin<Vertex>(SurfaceView::ALL);
		//iterEnd = dd->end<Vertex>(SurfaceView::ALL);
		for(;iter != iterEnd; ++iter)
		{

			//	get vertex---------------------------------------------12301
			Vertex* V = *iter;
			Vertex* v11;
			Vertex* v12;
			vector<Vertex*> Vrs;
			vector<Vertex*> Vls;
			vector<Vertex*> Vs;
			ug::DenseVector<ug::FixedArray1<double, 2> > temp;
			std::vector<size_t> ind;
			aaPos[V];
			// cout<<"aaPos[V] :"<<aaPos[V]<<endl;
			// cout<<"marker   :"<<Dsh->get_subset_index(V)<<endl;
			// cout<<"coarse index = "<<umap[V]<<endl;
			// 1/7/2020 -----------------------------------------------1214

			v12 = V;
			for(size_t i=markers[mk];i<markers[mk+1];i++)
			{
				v11 = Vnext(v12,domain,'r');
				Vrs.push_back(v11);
				v12 = v11;
			}
			v12 = V;
			for(size_t i=markers[mk-1];i<markers[mk];i++)
			{
				v11 = Vnext(v12,domain,'l');
				Vls.push_back(v11);
				v12 = v11;
			}
			// 1/7/2020 -----------------------------------------------1214
			map<Vertex*, double> CV;
			// the coefficients of V's basis functions ---------------
			CV[V] = 1.0;
			Vs.push_back(V);
			//cout<<"V: "<<aaPos[V]<<endl;
			//cout<<"Vrs.size(): "<<Vrs.size()<<endl;
			//cout<<"Vrs.back(): "<<aaPos[Vrs.back()]<<endl;
			for(size_t j=0;j<Vrs.size();j++)
			{
				CV[Vrs[j]] =((double)Vrs.size()-1.0-(double)j )/(double)Vrs.size();
				Vs.push_back(Vrs[j]);
				//cout<<aaPos[Vrs[j]]<<";  "<<CV[Vrs[j]]<<endl;
			}
			//cout<<"Vls.size(): "<<Vls.size()<<endl;
			//cout<<"Vls.back(): "<<aaPos[Vls.back()]<<endl;
			for(size_t j=0;j<Vls.size();j++)
			{
				CV[Vls[j]] =((double)Vls.size()-1.0-(double)j )/(double)Vls.size();
				Vs.push_back(Vls[j]);
				//cout<<aaPos[Vls[j]]<<";  "<<CV[Vls[j]]<<endl;
			}
			//--------------------------------------------------------12301
			//assemble coarse vector from fine vector
			temp = 0;
			for(size_t j=0;j<Vs.size();j++)
			{
				dd->inner_algebra_indices(Vs[j],  ind);
				temp += CV[Vs[j]]*u[ind[0]];
				//cout<<"CV[Vs[j]] :"<<CV[Vs[j]]<<endl;
			}
			vec[umap[V]] = temp(0,0);
			//cout<<M(umap[V],umap[Vos[j]])<<endl;
		}

	}
}
//------------------------------------------------------------------------------
template <typename TDomain, typename TGridFunction, typename TAlgebra>
void Vector3dCompare(TGridFunction &u,TAlgebra &vec,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd)
{
	//ug::DenseVector<ug::FixedArray1<double, 2> > temp;
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//  get marker accessor
	typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	std::vector<size_t> ind;
	// 12301------------
	map<Vertex*, size_t> umap, Mmap;
	//	The unkowns collected by chosen markers
	size_t MarkerSize = Dsh->num_subsets();
	cout<<"marker max size is "<<MarkerSize<<endl;
	size_t num = MarkerSize-3;
	size_t Ms = 10;
	vector<size_t> markers;
	markers.push_back(3);
	for(size_t i=1; i<=(num-num % Ms)/Ms;i++)
	{
		markers.push_back(2+Ms*i);
		//cout << i <<"  "<< 2+Ms*i<<endl;
	}
	markers.push_back(MarkerSize-1);
	// the map from Vertex* to index of coarse matrix
	size_t jc = 0;
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		iter = dd->begin<Vertex>(markers[mk], SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(markers[mk], SurfaceView::ALL);
		for(;iter != iterEnd; ++iter)
		{
			//	get vertex
			//Vertex* v = *iter;
			Mmap[*iter] = mk;
			umap[*iter] = jc;
			jc+=1;
			//cout<<umap[*iter]<<endl;
		}
	}
	//cout<<"jc is: "<<jc<<endl;
	//vec.resize(jc);
	//cout<<markers.size()<<endl;
	//cout<<markers[0]<<endl;
	size_t marker;
	//int mk = 1;
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		marker = markers[mk];
		iter = dd->begin<Vertex>(marker, SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(marker, SurfaceView::ALL);
		//iter = dd->begin<Vertex>(SurfaceView::ALL);
		//iterEnd = dd->end<Vertex>(SurfaceView::ALL);
		cout<<marker<<"---------------------------------------"<<endl;
		for(;iter != iterEnd; ++iter)
		{
			//	get vertex---------------------------------------------12301
			Vertex* V = *iter;
			ug::DenseVector<ug::FixedArray1<double, 2> > temp;
			std::vector<size_t> ind;
			aaPos[V];
			dd->inner_algebra_indices(V, ind);
			cout<<u[ind[0]]-vec[umap[V]]<<endl;
			// 1/7/2020 -----------------------------------------------1214
		}
	}
}
//------------------------------------------------------------------------------

template<class TGridFunction,typename TAlgebra>
void CoarseM(TGridFunction &v,TAlgebra &M,TAlgebra &N ){
	//this function is to disply matrix
	//const static int dim = TGridFunction::domain_type::dim;
	/*
	ug::DenseMatrix<ug::FixedArray2<double, 2, 2> > temp;
	//-----------------------------------------------------
	//M.print();
	int gridsize = 10;
	M.resize_and_clear(gridsize, gridsize);
	v.resize(gridsize);
	for(int i = 0; i < gridsize; ++i)
	{   v[i] = (double)(i);
		M(i,i) = 1.52310;
	}
	temp= M(0,0);
	cout<<"M(0,0): "<<M(0,0)<<endl;
	cout<<"temp(0,0): "<<temp(0,0)<<endl;
	cout<<"temp: "<<temp<<endl;
	M(1,1) = temp(0,0);
	cout<<"M(1,1): "<<M(1,1)<<endl;
	 */
	Matrix3dFineToCoarse(v,M,N,v.domain(),v.dof_distribution());

	//N.print();
	//cout<<M(i,j)<<endl;
	//M.print();
	//M.num_rows();
	//v.resize(M.num_rows());
	//N.print();
}

template<class TGridFunction,typename TAlgebra>
void CoarseV(TGridFunction &r3d,TAlgebra &vec){
	Vector3dFineToCoarse(r3d,vec,r3d.domain(),r3d.dof_distribution());
	//-----------------------------------------------------

}

template<class TGridFunction,typename TAlgebra>
void ComV(TGridFunction &r3d,TAlgebra &vec){
	Vector3dCompare(r3d,vec,r3d.domain(),r3d.dof_distribution());
	//-----------------------------------------------------

}
//end test ------------------------------------------------------

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
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename TAlgebra::matrix_type matrix_type;
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
			reg.add_function("LoadVector22", &LoadVector22<function_type>, grp);
			reg.add_function("Assign22", &Assign22<function_type, vector_type>, grp);
			reg.add_function("Extract3dto1d", &ExtractPositions3dto1d<function_type, vector_type>, grp);
			reg.add_function("Prolongation1dto3d", &ExtractPositions1dto3d<function_type, vector_type>, grp);
			reg.add_function("Restriction", &Restriction3dto1d<function_type>, grp);

			reg.add_function("CoarseM", &CoarseM<function_type,matrix_type>, grp);
			reg.add_function("CoarseV", &CoarseV<function_type,vector_type>, grp);
			reg.add_function("ComV", &ComV<function_type,vector_type>, grp);
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
		typedef typename TAlgebra::vector_type vector_type;
		//typedef typename TAlgebra::matrix_type matrix_type;

		reg.add_function("disp", &disp<vector_type>, grp);
		//reg.add_function("dispM", &dispM<matrix_type>, grp);

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
