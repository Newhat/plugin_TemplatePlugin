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
#include<set>
#include <string>
#include <cmath>  // for isinf, isnan
#include <algorithm>
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
	vector<Vertex*> Vlists;
	//vector<Vertex*> Vtest[10];
	vector<size_t> Vmarkers;
	void print_hello () const	{UG_LOG("hello\n");}
	vector<Vertex*> get_Vlists() {return Vlists;}
	vector<size_t> get_Vmarkers() {return Vmarkers;}
	//vector<Vertex*> get_Vtest() {return Vtest;}
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
			<< dimension << "d from " << filename
			<< " , " << gridsize << " x " << gridsize);
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
	int marker; //all inner points, no boundary points
	for(marker = 3;marker<=MarkerSize;marker++)
	{	//cout<<"marker is: "; //----------------------test here
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

//-----------------------------------------------------------------------------
// assign subsets for Non-zero boundary set
//-----------------------------------------------------------------------------
template <typename TDomain, typename TGridFunction>
void Assign_Subsets_Non0_Boundary(TGridFunction &u,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd,int Meshcounts,double Meshbegin,double Meshsize,vector<Vertex*>* Vlis)
{
	//Meshcounts is the number of edges for 1D mesh
	//ug::DenseVector<ug::FixedArray1<double, 2> > temp;
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//  get marker accessor
	//typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	int i, j;
	double tempX;
	map<int, int> umap;
	map<Vertex*, int> Vmap;
	for(j=0;j<(Meshcounts+1);j++)
	{
		i=round( (Meshbegin+(double)(j)*Meshsize)*(double)round(2/Meshsize) );
		umap[i] = j;
	}
	int counter = 0;
	iter = dd->begin<Vertex>(1, SurfaceView::ALL);
	iterEnd = dd->end<Vertex>(1, SurfaceView::ALL);
	for(;iter != iterEnd; ++iter)
	{
			counter+=1;
			tempX = aaPos[*iter].x()*(double)round(2/Meshsize);
			Vmap[*iter] = umap[round(tempX)];
			//cout << counter<<", "<<aaPos[V] << ", "<<Vmap[V] <<endl;
	}

	//vector<Vertex*> Vlis[Meshcounts+1]; //-------------
	map<Vertex*, int>::iterator it;
	for(it = Vmap.begin();it != Vmap.end();++it)
	{
		Vlis[it->second].push_back(it->first);

	}
/*
	for(int i=0;i<Meshcounts+1;++i)
	{
		for(size_t j=0;j<Vlis[i].size();++j)
		{
			cout<<i<<", "<<aaPos[Vlis[i][j]]<<endl;
		}
	}
*/
}
//----------------------------------------------------------------------------------
template<class TGridFunction,typename TAlgebra>
//vector<vector<Vertex*>>
void Subsets_Non0_Boundary(TGridFunction &u,TAlgebra &vec,
		size_t Meshcounts,double Meshbegin,double Meshsize)
{
	//vec can be used later such as non-uniform Meshzise
	vector<Vertex*> Vlists[Meshcounts+1];
	Assign_Subsets_Non0_Boundary(u,u.domain(),u.dof_distribution(),Meshcounts,Meshbegin,Meshsize,Vlists);
	//cout<<Vlists[Meshcounts].size()<<endl;
	//-----------------------------------------------------
}

//-----------------------------------------------------------------------------
// matrix 3d to 1d, Non-zero boundary
//-----------------------------------------------------------------------------
template <typename TDomain, typename TGridFunction, typename TAlgebra>
void Matrix3dFineTo1dNon0(TGridFunction &u,TAlgebra &M,TAlgebra &N,
		SmartPtr<TDomain> domain,SmartPtr<DoFDistribution> dd,
		vector<Vertex*>* Vlists)
{
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//  get marker accessor
	typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	std::vector<size_t> ind1,ind2;
	size_t Ms=1;
	// 12301------------
	map<Vertex*, size_t> umap;// Mmap;
	//	The unkowns collected by chosen markers
	size_t MarkerSize = Dsh->num_subsets();
	//cout<<"marker max size is "<<MarkerSize<<endl;
	size_t num = MarkerSize-3;
	//size_t Ms = 10;
	vector<size_t> markers;
	markers.push_back(3);
	for(size_t i=1; i<=(num-num % Ms)/Ms;i++)
	{
		markers.push_back(2+Ms*i);
		//cout << i <<"  "<< 2+Ms*i<<endl;
	}
	markers.push_back(MarkerSize-1); //here comes the problem 2+Ms*i can be MarkerSize-1
	// the map from Vertex* to index of coarse matrix
	size_t jc = 0;
	for(size_t i = 0;i<MarkerSize-1;i++)
	{
		for(size_t j=0;j<Vlists[i].size();j++)
		{
			umap[Vlists[i][j]] = jc;
		}
		jc+=1;
	}
	//cout<<"jc is: "<<jc<<endl;
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		//cout<<markers[mk]<<endl;
		iter = dd->begin<Vertex>(markers[mk], SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(markers[mk], SurfaceView::ALL);
		for(;iter != iterEnd; ++iter)
		{
			umap[*iter] = jc;
		}
		jc+=1;
	}
	cout<<"jc is: "<<jc<<endl;
	M.resize_and_clear(jc, jc);
	for(size_t k = 0;k<MarkerSize-1;k++)
	{
		M(k,k) = 1.0;
	}

	for(size_t mk = 0;mk<MarkerSize-1;mk++)
	{
		ug::DenseMatrix<ug::FixedArray2<double, 2, 2> > temp,temp1;
		temp=0.0;
		temp1=0.0;
		for(size_t jk=0;jk<Vlists[mk].size();jk++)
		{
			//	get vertex---------------------------------------------
			Vertex* V = Vlists[mk][jk];
			vector<Vertex*> Vos;
			//Vos collects points around V and V
			//cout<<aaPos[V]<<endl;
			aaPos[V];
			//Vos.push_back(V);
			//  collected edges
			MultiGrid::edge_traits::secure_container edges;
			//nearest points around point V--------------------------------
			domain->grid()->associated_elements(edges, V);
			for(size_t i = 0; i < edges.size(); ++i)
			{
				Vertex* Vn=NULL;
				Edge* e = edges[i];
				for(size_t j=0;j<2;j++)
				{
					if((e->vertex(j)!=V)
							&&(Dsh->get_subset_index(e->vertex(j))>Dsh->get_subset_index(V)))
					{
						Vn = e->vertex(j);
						Vos.push_back(Vn);
						//cout<<"nby verts: "<<aaPos[Vn]<<endl;
					}
				}
			}

			//assemble boundary parts of coarse matrix from fine matrix
			for(size_t j=0;j<Vos.size();j++)
			{
				dd->inner_algebra_indices(V, ind1);
				dd->inner_algebra_indices(Vos[j], ind2);
				M(umap[Vos[j]],umap[V])+=N(ind2[0],ind1[0]);
			}
		}
		//cout<<"---------------------------------"<<endl;
	}
}

template<class TGridFunction,typename TAlgebra>
void Matrix3dto1dNon0(TGridFunction &v,TAlgebra &M,TAlgebra &N,TGridFunction &vec,
		size_t Meshcounts,double Meshbegin,double Meshsize)
{
	vector<Vertex*> Vlists[Meshcounts+1];
	Assign_Subsets_Non0_Boundary(v,v.domain(),v.dof_distribution(),Meshcounts,Meshbegin,Meshsize,Vlists);
	Matrix3dFineTo1dNon0(v,M,N,v.domain(),v.dof_distribution(),Vlists);
	vec.resize(2*Meshcounts);
	for(size_t i=0;i<vec.size();i++)
		vec[i] = double(i);
}
//-----------------------------------------------------------------------------
// matrix 3d to 1d, Non-zero boundary end
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//-------------------------- matrix 3d to 1d, zero boundary -------------------
//-----------------------------------------------------------------------------
template <typename TDomain, typename TGridFunction, typename TAlgebra>
void Matrix3dFineTo1d(TGridFunction &u,TAlgebra &M,TAlgebra &N,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd)
{
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//  get marker accessor
	typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	std::vector<size_t> ind1,ind2;
	size_t Ms=1;
	// 12301------------
	map<Vertex*, size_t> umap;// Mmap;
	//	The unkowns collected by chosen markers
	size_t MarkerSize = Dsh->num_subsets();
	//cout<<"marker max size is "<<MarkerSize<<endl;
	size_t num = MarkerSize-3;
	//size_t Ms = 10;
	vector<size_t> markers;
	markers.push_back(3);
	for(size_t i=1; i<=(num-num % Ms)/Ms;i++)
	{
		markers.push_back(2+Ms*i);
		//cout << i <<"  "<< 2+Ms*i<<endl;
	}
	markers.push_back(MarkerSize-1); //here comes the problem 2+Ms*i can be MarkerSize-1
	// the map from Vertex* to index of coarse matrix
	size_t jc = 2;
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		//cout<<markers[mk]<<endl;
		iter = dd->begin<Vertex>(markers[mk], SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(markers[mk], SurfaceView::ALL);
		for(;iter != iterEnd; ++iter)
		{
			//	get vertex
			//Vertex* v = *iter;
			//Mmap[*iter] = mk;
			umap[*iter] = jc;
			//cout<<umap[*iter]<<endl;
		}
		jc+=1;
	}
	cout<<"jc is: "<<jc<<endl;
	M.resize_and_clear(jc, jc);
	M(0,0)=1.0;
	M(1,1)=1.0;
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
		ug::DenseMatrix<ug::FixedArray2<double, 2, 2> > temp,temp1;
		temp=0.0;
		temp1=0.0;
		//M(mk+1,mk+1)=0;
		//M(mk+1,mk+2)=0;
		//M(mk+2,mk+1)=0;
		for(;iter != iterEnd; ++iter)
		{

			//	get vertex---------------------------------------------12301
			Vertex* V = *iter;
			//Vertex* v11 = NULL;
			//Vertex* v12 = NULL;
			//vector<Vertex*> Vrs;
			//vector<Vertex*> Vls;
			vector<Vertex*> Vos;
			//Vos collects points around V and V
			//cout<<aaPos[V]<<endl;
			aaPos[V];
			// dd->inner_algebra_indices(V, ind);
			// 1/7/2020 -----------------------------------------------1214
			Vos.push_back(V);
			// 1/7/2020 -----------------------------------------------1214
			//  collected edges
			//set<Vertex*> setVs;
			MultiGrid::edge_traits::secure_container edges;
			//nearest points around point V--------------------------------
			domain->grid()->associated_elements(edges, V);
			//cout<<"how many edges associated with V"<<endl;
			//cout<<edges.size()<<endl;
			for(size_t i = 0; i < edges.size(); ++i)
			{
				//vector<Vertex*> Vnrs;
				//vector<Vertex*> Vnls;
				Vertex* Vn=NULL;
				//cout<<i<<endl;
				Edge* e = edges[i];
				//cout<<"-------------------"<<endl;
				for(size_t j=0;j<2;j++)
				{
					if((e->vertex(j)!=V)
							&&(Dsh->get_subset_index(e->vertex(j))>=Dsh->get_subset_index(V)))
					{
						Vn = e->vertex(j);
						//setVs.insert(Vn);
						//setVs.clear();
						Vos.push_back(Vn);
						//cout<<"nby verts: "<<aaPos[Vn]<<endl;
					}
				}
			}
			//--------------------------------------------------------12301
			//assemble coarse matrix from fine matrix

			for(size_t j=0;j<Vos.size();j++)
			{
				dd->inner_algebra_indices(V, ind1);
				dd->inner_algebra_indices(Vos[j], ind2);
				if(umap[Vos[j]]==umap[V])
				{
					temp+= N(ind2[0],ind1[0]);
				}
				if(umap[Vos[j]]>umap[V])
				{
					temp1 += N(ind2[0],ind1[0]);
				}
			}

		}
		M(mk+1,mk+1)=temp(0,0);
		if(mk+2<jc)
		{
			M(mk+1,mk+2)=temp1(0,0);
			M(mk+2,mk+1)=temp1(0,0);
		}
		//cout<<mk+2<<endl;
	}
	//cout<<"mk+2 is done  --------"<<endl;
}

template<class TGridFunction,typename TAlgebra>
void Matrix3dto1d(TGridFunction &v,TAlgebra &M,TAlgebra &N)
{
	Matrix3dFineTo1d(v,M,N,v.domain(),v.dof_distribution());
	cout<<" "<<endl;
}
//-----------------------------------------------------------------------------
// matrix 3d to 1d, zero boundary end -----------------------------------------
//-----------------------------------------------------------------------------

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
	if(v!=NULL)
	{
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


	Vertex* v11=NULL;
	Vertex* v12=NULL;
	vector<Vertex*> Vrs;
	vector<Vertex*> Vls;
	vector<Vertex*> Vnrs;
	vector<Vertex*> Vnls;
	vector<Vertex*> Vs;
	vector<Vertex*> Vns;
	map<Vertex*, double> CV,CVn;

	Dsh->get_subset_index(V);

	std::vector<size_t> ind, ind1,ind2;
	//Dsh->get_subset_index(V);
	//Dsh->get_subset_index(Vrs[0]);
	aaPos[V];
	//cout<<V<<endl;

	//1, basis for V point -------------------------------
	v12 = V;

	//cout<<Mmap[V]+1<<endl;

	for(size_t i=markers[Mmap[V]];i<markers[Mmap[V]+1];i++)
	{
		v11 = Vnext(v12,domain,'r');
		if(v11!=NULL)
			Vrs.push_back(v11);
		v12 = v11;
	}

	v11 = NULL;
	v12 = V;
	for(size_t i=markers[Mmap[V]-1];i<markers[Mmap[V]];i++)
	{
		v11 = Vnext(v12,domain,'l');
		if(v11!=NULL)
			Vls.push_back(v11);
		v12 = v11;
	}
	v11 = NULL;

	//2, basis for Vn point -------------------------------
	v12 = Vn;
	for(size_t i=markers[Mmap[Vn]];i<markers[Mmap[Vn]+1];i++)
	{
		v11 = Vnext(v12,domain,'r');
		if(v11!=NULL)
			Vnrs.push_back(v11);
		v12 = v11;
		//cout<<aaPos[v11]<<endl;
	}
	v11 = NULL;
	v12 = Vn;
	for(size_t i=markers[Mmap[Vn]-1];i<markers[Mmap[Vn]];i++)
	{
		v11 = Vnext(v12,domain,'l');
		if(v11!=NULL)
			Vnls.push_back(v11);
		v12 = v11;
		//cout<<aaPos[v11]<<endl;
	}
	v11 = NULL;
	// the coefficients of V's basis functions ---------------
	CV[V] = 1.0;
	Vs.push_back(V);
	size_t Len;

	//cout<<"Vrs.size(): "<<Vrs.size()<<endl;
	//cout<<"Vrs.back(): "<<aaPos[Vrs.back()]<<endl;
	for(size_t j=0;j<Vrs.size();j++)
	{
		Len = Vrs.size();
		if(Mmap[V]==markers.size()-2){Len +=1;}
		CV[Vrs[j]] =((double)Len-1.0-(double)j )/(double)Len;
		Vs.push_back(Vrs[j]);
		//cout<<aaPos[Vrs[j]]<<";  "<<CV[Vrs[j]]<<endl;
	}
	//cout<<"Vls.size(): "<<Vls.size()<<endl;
	//cout<<"Vls.back(): "<<aaPos[Vls.back()]<<endl;
	//cout<<"V: "<<aaPos[V]<<endl;
	for(size_t j=0;j<Vls.size();j++)
	{
		Len = Vls.size();
		if(Mmap[V]==1){Len +=1;}
		CV[Vls[j]] =((double)Len-1.0-(double)j )/(double)Len;
		Vs.push_back(Vls[j]);
		//if(Mmap[V]==1)
		//cout<<aaPos[Vls[j]]<<";  "<<CV[Vls[j]]<<endl;
	}
	// the coefficients of Vn's basis functions ---------------
	CVn[Vn] = 1.0;
	Vns.push_back(Vn);
	for(size_t j=0;j<Vnrs.size();j++)
	{
		Len = Vnrs.size();
		if(Mmap[Vn]==markers.size()-2){Len +=1;}
		CVn[Vnrs[j]] =((double)Len-1.0-(double)j )/(double)Len;
		Vns.push_back(Vnrs[j]);
	}
	//if(Mmap[Vn]==1){cout<<"Vn: "<<aaPos[Vn]<<endl;}
	for(size_t j=0;j<Vnls.size();j++)
	{
		Len = Vnls.size();
		if(Mmap[Vn]==1){Len +=1;}
		CVn[Vnls[j]] =((double)Len-1.0-(double)j )/(double)Len;
		Vns.push_back(Vnls[j]);
		//if(Mmap[Vn]==1){cout<<aaPos[Vnls[j]]<<";  "<<CVn[Vnls[j]]<<endl;}
	}
	// Matrix -----------------------------------------------------------------
	//cout<<N(ind2[0],ind2[0])<<endl;
	//cout<<"Vs.size(): "<<Vs.size()<<endl;
	//cout<<"Vns.size(): "<<Vns.size()<<endl;

	temp = 0.0;
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
		SmartPtr<DoFDistribution> dd,size_t Ms)
{
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//  get marker accessor
	typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	std::vector<size_t> ind,ind1,ind2;
	// 12301------------
	map<Vertex*, size_t> umap, Mmap;
	//	The unkowns collected by chosen markers
	size_t MarkerSize = Dsh->num_subsets();
	cout<<"marker max size is "<<MarkerSize<<endl;
	size_t num = MarkerSize-3;
	//size_t Ms = 10;
	vector<size_t> markers;
	markers.push_back(3);
	for(size_t i=1; i<=(num-num % Ms)/Ms;i++)
	{
		markers.push_back(2+Ms*i);
		//cout << i <<"  "<< 2+Ms*i<<endl;
	}
	markers.push_back(MarkerSize-1); //here comes the problem 2+Ms*i can be MarkerSize-1
	// the map from Vertex* to index of coarse matrix
	size_t jc = 0;
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		//cout<<markers[mk]<<endl;
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
			Vertex* v11 = NULL;
			Vertex* v12 = NULL;
			//vector<Vertex*> Vrs;
			//vector<Vertex*> Vls;
			vector<Vertex*> Vos;
			//Vos collects points around V and V
			//cout<<aaPos[V]<<endl;
			aaPos[V];
			// dd->inner_algebra_indices(V, ind);
			// cout<<" :"<<u[ind[0]]<<endl;
			/*
			if(umap[V]==260 || umap[V]==121) //test 1/19/2020
			cout<<umap[V]<<", "<<"marker :"<<marker<<", "<<aaPos[V]<<endl;
			 */
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
			if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)<=markers.back())
			{
				Vos.push_back(v11);
				//cout<<aaPos[Vrs.back()]<<endl;
			}
			v11 = NULL;
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
			if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)>=markers[0])
			{
				Vos.push_back(v11);
				//cout<<aaPos[Vls.back()]<<endl;
			}
			v11 = NULL;
			// 1/7/2020 -----------------------------------------------1214
			//  collected edges
			set<Vertex*> setVs;
			MultiGrid::edge_traits::secure_container edges;
			//nearest points around point V--------------------------------
			domain->grid()->associated_elements(edges, V);
			//cout<<"how many edges associated with V"<<endl;
			//cout<<edges.size()<<endl;
			for(size_t i = 0; i < edges.size(); ++i)
			{
				//vector<Vertex*> Vnrs;
				//vector<Vertex*> Vnls;
				Vertex* Vn=NULL;
				//cout<<i<<endl;
				Edge* e = edges[i];
				//cout<<"-------------------"<<endl;
				for(size_t j=0;j<2;j++)
				{
					if((e->vertex(j)!=V)
							&&(Dsh->get_subset_index(e->vertex(j))==Dsh->get_subset_index(V)))
					{
						Vn = e->vertex(j);
						setVs.insert(Vn);
						//setVs.clear();
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
						if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)<=markers.back())
						{
							Vos.push_back(v11);
							//cout<<"right: "<<aaPos[v11]<<endl;
						}
						v11 = NULL;
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
						if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)>=markers[0])
						{
							Vos.push_back(v11);
							//cout<<"left: "<<aaPos[v11]<<endl;
						}
						v11 = NULL;
						//cout<<"center v : "<<aaPos[v]<<endl;
						//cout<<"nby verts: "<<aaPos[Vn]<<endl;
					}
					/*
					//-----------------------------------------------------------
					//-----------add more Vn 1/28/2020
					//-----------------------------------------------------------
					if((e->vertex(j)!=V)
							&&(Dsh->get_subset_index(e->vertex(j))!=Dsh->get_subset_index(V))) //------working
					{
						Vn=NULL;
						v11=NULL;
						int Vnr, Vnl;
						Vnr=0;
						Vnl=0;
						//cout<<", nby verts: "<<aaPos[e->vertex(j)]<<endl;
						if(Vnext(e->vertex(j),domain,'r')!=NULL)
							Vnr=Dsh->get_subset_index(Vnext(e->vertex(j),domain,'r'));
						if(Vnext(e->vertex(j),domain,'l')!=NULL)
							Vnl=Dsh->get_subset_index(Vnext(e->vertex(j),domain,'l'));

						if(Vnr==Dsh->get_subset_index(V))
						{
							Vn = Vnext(e->vertex(j),domain,'r');
							if(std::find(Vos.begin(), Vos.end(), Vn) != Vos.end()) {
								//cout<<"Found X"<<endl;
								Vn=NULL;
							} else{

								//cout<<"can't find X"<<endl;
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
								if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)<=markers.back())
								{
									Vos.push_back(v11);
									//cout<<"right: "<<aaPos[v11]<<endl;
								}
								v11 = NULL;

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
								if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)>=markers[0])
								{
									Vos.push_back(v11);
									//cout<<"left: "<<aaPos[v11]<<endl;
								}
								v11 = NULL;
								//cout<<"center v : "<<aaPos[v]<<endl;
								//cout<<"nby verts: "<<aaPos[Vn]<<endl;

							}
						}

						if(Vnl==Dsh->get_subset_index(V))
						{
							Vn = Vnext(e->vertex(j),domain,'l');
							if(std::find(Vos.begin(), Vos.end(), Vn) != Vos.end()) {
								//cout<<"Found X"<<endl;
								Vn=NULL;
							} else {
								//cout<<"can't find X"<<endl;
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
								if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)<=markers.back())
								{
									Vos.push_back(v11);
									//cout<<"right: "<<aaPos[v11]<<endl;
								}
								v11 = NULL;
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
								if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)>=markers[0])
								{
									Vos.push_back(v11);
									//cout<<"left: "<<aaPos[v11]<<endl;
								}
								v11 = NULL;
								//cout<<"center v : "<<aaPos[v]<<endl;
								//cout<<"nby verts: "<<aaPos[Vn]<<endl;
							}
						}
					}
					 */
					//-----------------------------------------------------------
					//-----------add more Vn 1/28/2020 End
					//---------------------------------------------------------
				}
			}
			/*
			//------------------------------------------------------------------
			//put second nearest points in to Vos
			//------------------------------------------------------------------
			for(set<Vertex*>::iterator it=setVs.begin(); it!=setVs.end(); ++it)
			{
				domain->grid()->associated_elements(edges, *it);
				//cout<<"how many edges associated with V"<<endl;
				//cout<<edges.size()<<endl;
				for(size_t i = 0; i < edges.size(); ++i)
				{
					//vector<Vertex*> Vnrs;
					//vector<Vertex*> Vnls;
					Vertex* Vn=NULL;
					//cout<<i<<endl;
					Edge* e = edges[i];
					//cout<<"-------------------"<<endl;
					for(size_t j=0;j<2;j++)
					{
						if((e->vertex(j)!=*it)
								&&(Dsh->get_subset_index(e->vertex(j))==Dsh->get_subset_index(V)))
						{
							Vn=e->vertex(j);
							v11=NULL;
							//cout<<", nby verts: "<<aaPos[e->vertex(j)]<<endl;
							if(std::find(Vos.begin(), Vos.end(), Vn) != Vos.end())
							{
								//cout<<"Found X"<<endl;
								Vn=NULL;
							} else{

								//cout<<"can't find X"<<endl;
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
								if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)<=markers.back())
								{
									Vos.push_back(v11);
									//cout<<"right: "<<aaPos[v11]<<endl;
								}
								v11 = NULL;

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
								if(v11!=NULL && (size_t)Dsh->get_subset_index(v11)>=markers[0])
								{
									Vos.push_back(v11);
									//cout<<"left: "<<aaPos[v11]<<endl;
								}
								v11 = NULL;
								//cout<<"center v : "<<aaPos[v]<<endl;
								//cout<<"nby verts: "<<aaPos[Vn]<<endl;
							}
						}
					}
				}
			}
			//------------------------------------------------------------------
			//put second nearest points in to Vos End
			//------------------------------------------------------------------
			*/
			//assemble coarse matrix from fine matrix
			ug::DenseMatrix<ug::FixedArray2<double, 2, 2> > temp;
			for(size_t j=0;j<Vos.size();j++)
			{
				if(umap[Vos[j]]>=umap[V])
				{	temp=0.0;
				temp = CoarseFromFine(V,Vos[j],dd,domain,N,Mmap,markers);
				M(umap[V],umap[Vos[j]]) = temp(0,0);
				M(umap[Vos[j]],umap[V]) = temp(0,0);
				//cout<<"V :"<<aaPos[V]<<", Vos["<<j<<"] :"<<aaPos[Vos[j]]<<endl;
				dd->inner_algebra_indices(V, ind1);
				dd->inner_algebra_indices(Vos[j], ind2);
				/*
					if(M(umap[V],umap[Vos[j]])-N(ind2[0],ind1[0])!=0)
					cout<<umap[V]<<","<<umap[Vos[j]]<<endl;
					//cout<<M(umap[V],umap[Vos[j]])-N(ind1[0],ind2[0])<<endl;
					if(ind1[0]-umap[V]!=10** || ind2[0]-umap[Vos[j]]!=10**)
					cout<<ind1[0]-umap[V]<<", "<<ind2[0]-umap[Vos[j]]<<endl;
				 */
				}
			}

		}

	}
}

template <typename TDomain, typename TGridFunction, typename TAlgebra>
void Vector3dFineToCoarse(TGridFunction &u,TAlgebra &vec,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd,size_t Ms)
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
	//cout<<"marker max size is "<<MarkerSize<<endl;
	size_t num = MarkerSize-3;
	//size_t Ms = 10;
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
			Vertex* v11=NULL;
			Vertex* v12=NULL;
			vector<Vertex*> Vrs;
			vector<Vertex*> Vls;
			vector<Vertex*> Vs;
			ug::DenseVector<ug::FixedArray1<double, 2> > temp;
			std::vector<size_t> ind;
			size_t Len;
			aaPos[V];
			// cout<<"aaPos[V] :"<<aaPos[V]<<endl;
			// cout<<"marker   :"<<Dsh->get_subset_index(V)<<endl;
			// cout<<"coarse index = "<<umap[V]<<endl;
			// 1/7/2020 -----------------------------------------------1214

			v12 = V;
			for(size_t i=markers[mk];i<markers[mk+1];i++)
			{
				v11 = Vnext(v12,domain,'r');
				if(v11!=NULL)
					Vrs.push_back(v11);
				v12 = v11;
			}
			v11 = NULL;
			v12 = V;
			for(size_t i=markers[mk-1];i<markers[mk];i++)
			{
				v11 = Vnext(v12,domain,'l');
				if(v11!=NULL)
					Vls.push_back(v11);
				v12 = v11;
			}
			v11 = NULL;
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
				Len = Vrs.size();
				if(Mmap[V]==markers.size()-2){Len +=1;}
				CV[Vrs[j]] =((double)Len-1.0-(double)j )/(double)Len;
				Vs.push_back(Vrs[j]);
				//cout<<aaPos[Vrs[j]]<<";  "<<CV[Vrs[j]]<<endl;
			}
			//cout<<"Vls.size(): "<<Vls.size()<<endl;
			//cout<<"Vls.back(): "<<aaPos[Vls.back()]<<endl;
			for(size_t j=0;j<Vls.size();j++)
			{
				Len = Vls.size();
				if(Mmap[V]==1){Len +=1;}
				CV[Vls[j]] =((double)Len-1.0-(double)j )/(double)Len;
				Vs.push_back(Vls[j]);
				//if(Mmap[V]==1)
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
			//cout<<umap[V]<<", :"<<vec[umap[V]]<<endl; //test
			//dd->inner_algebra_indices(V,  ind);
			/*
			if(vec[umap[V]]-u[ind[0]]!=0)
			{
				cout<<vec[umap[V]]-u[ind[0]]<<endl;
				cout<<"V :"<<aaPos[V]<<endl;
			}
			 */
		}

	}
}
//------------------------------------------------------------------------------
template <typename TDomain, typename TGridFunction, typename TAlgebra>
void Vector3dCoarseToFine(TGridFunction &u,TAlgebra &vec,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd,size_t Ms)
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
	//cout<<"marker max size is "<<MarkerSize<<endl;
	size_t num = MarkerSize-3;
	//size_t Ms = 10;
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
	// Main part below----------------------------------------------------------------------
	for(size_t mk = 1;mk<(markers.size()-1);mk++)
	{
		marker = markers[mk];
		iter = dd->begin<Vertex>(marker, SurfaceView::ALL);
		iterEnd = dd->end<Vertex>(marker, SurfaceView::ALL);
		//iter = dd->begin<Vertex>(SurfaceView::ALL);
		//iterEnd = dd->end<Vertex>(SurfaceView::ALL);
		//cout<<marker<<"---------------------------------------"<<endl;
		for(;iter != iterEnd; ++iter)
		{
			//	get vertex---------------------------------------------12301
			Vertex* V = *iter;
			Vertex* v11=NULL;
			Vertex* v12=NULL;
			vector<Vertex*> Vrs;
			vector<Vertex*> Vls;

			ug::DenseVector<ug::FixedArray1<double, 2> > temp;
			std::vector<size_t> ind;
			aaPos[V];
			dd->inner_algebra_indices(V, ind);
			//the error between A3d^-1*res and M1^-1*Restriction(res)-------------error test
			//cout<<u[ind[0]]-vec[umap[V]]<<endl;
			// 1/7/2020 -----------------------------------------------1214
			v12 = V;
			for(size_t i=markers[mk];i<markers[mk+1];i++)
			{
				v11 = Vnext(v12,domain,'r');
				if(v11!=NULL)
					Vrs.push_back(v11);
				v12 = v11;
			}
			v11 = NULL;
			if(mk==1)
			{
				v12 = V;
				for(size_t i=markers[mk-1];i<markers[mk];i++)
				{
					v11 = Vnext(v12,domain,'l');
					if(v11!=NULL)
						Vls.push_back(v11);
					v12 = v11;
				}
			}
			v11 = NULL;
			// 1/7/2020 -----------------------------------------------1214
			map<Vertex*, double> CV,CVL;
			size_t Len;
			// the coefficients of V's basis functions ---------------
			CV[V] = 1.0;
			CVL[V]= 0.0;
			//cout<<aaPos[V]<<";---"<<CV[V]<<endl;
			//cout<<"V: "<<aaPos[V]<<endl;
			//cout<<"Vrs.size(): "<<Vrs.size()<<endl;
			//cout<<"Vrs.back(): "<<aaPos[Vrs.back()]<<endl;
			for(size_t j=0;j<Vrs.size();j++)
			{
				Len = Vrs.size();
				if(mk==markers.size()-2){Len +=1;}
				CV[Vrs[j]] =((double)Len-1.0-(double)j )/(double)Len;
				CVL[Vrs[j]] =1.0-((double)Len-1.0-(double)j )/(double)Len;
				dd->inner_algebra_indices(Vrs[j], ind);
				if(mk==markers.size()-2){
					u[ind[0]] = CV[Vrs[j]]*vec[umap[V]];
				}
				else{
					u[ind[0]] = CV[Vrs[j]]*vec[umap[V]]+CVL[Vrs[j]]*vec[umap[Vrs.back()]];
				}
				//cout<<aaPos[Vrs[j]]<<";  "<<CV[Vrs[j]]<<endl;
				//cout<<aaPos[Vrs[j]]<<";  "<<CVL[Vrs[j]]<<endl;
			}
			if(mk==1) //L R sides uncorrect
			{
				dd->inner_algebra_indices(V, ind);
				u[ind[0]] = CV[V]*vec[umap[V]];
				//cout<<"V :"<<aaPos[V]<<endl;
				for(size_t j=0;j<Vls.size();j++)
				{
					Len = Vls.size()+1;
					CVL[Vls[j]] =((double)Len-1.0-(double)j )/(double)Len;
					dd->inner_algebra_indices(Vls[j], ind);
					u[ind[0]] = CVL[Vls[j]]*vec[umap[Vls[j]]];
					//cout<<aaPos[Vls[j]]<<";  "<<CVL[Vls[j]]<<endl;
				}
			}
		}
	}
}
//------------------------------------------------------------------------------

template<class TGridFunction,typename TAlgebra>
void CoarseM(TGridFunction &v,TAlgebra &M,TAlgebra &N,size_t Ms){
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
	Matrix3dFineToCoarse(v,M,N,v.domain(),v.dof_distribution(),Ms);

	//N.print();
	//cout<<M(i,j)<<endl;
	//M.print();
	//M.num_rows();
	//v.resize(M.num_rows());
	//N.print();
}

template<class TGridFunction,typename TAlgebra>
void CoarseV(TGridFunction &r3d,TAlgebra &vec,size_t Ms){
	Vector3dFineToCoarse(r3d,vec,r3d.domain(),r3d.dof_distribution(),Ms);
	//-----------------------------------------------------

}

template<class TGridFunction,typename TAlgebra>
void ComV(TGridFunction &r3d,TAlgebra &vec,size_t Ms){
	Vector3dCoarseToFine(r3d,vec,r3d.domain(),r3d.dof_distribution(),Ms);
	//-----------------------------------------------------

}
//-------------------------------------------------------------------------------------------------
// 1d-3d, then multigrid tet
//--------------------------------------------------------------------------------------------------
template <typename TDomain, typename TGridFunction, typename TAlgebra>
void Refined_Domain3dTo1d(TGridFunction &u,TAlgebra &b,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd,vector<Vertex*> Vlists)
{
	ug::DenseVector<ug::FixedArray1<double, 2> > temp;
	//	get position accessor
	//typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//	algebra indices vector
	std::vector<size_t> ind;
	//	Sz is size of each cross section, Meshcounts is how many subsets
	size_t Sz,Meshcounts;
	Meshcounts = b.size()-2;
	Sz = Vlists.size()/Meshcounts;
	for(size_t i=0;i<Meshcounts;++i)
	{
		temp[0] = 0.0;
		temp[1] = 0.0;
		for(size_t j=0;j<Sz;++j)
		{
			//cout<<Sz<<", "<<i+3<<", "<<aaPos[Vlists[i*Sz+j]]<<endl;
			//	load indices associated with vertex
			dd->inner_algebra_indices(Vlists[i*Sz+j], ind);
			temp=temp+u[ind[0]];
		}
		b[i+2] = temp[0];
	}


}
template<class TGridFunction,typename TAlgebra>
void RefinedSubsetsTet3Dto1D(TGridFunction &u,TAlgebra &vec,vector<Vertex*> Vlists)
{
	//the 3d to 1d map-----------------------------------------------------
	Refined_Domain3dTo1d(u,vec,u.domain(),u.dof_distribution(),Vlists);

}

//--------------------------------------------------------------------------------------------------
template <typename TDomain, typename TGridFunction, typename TAlgebra>
void Refined_Domain1dTo3d(TGridFunction &u,TAlgebra &b,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd,vector<Vertex*> Vlists)
{
	ug::DenseVector<ug::FixedArray1<double, 2> > temp;
	//	get position accessor
	//typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//	algebra indices vector
	std::vector<size_t> ind;
	//	Sz is size of each cross section, Meshcounts is how many subsets
	size_t Sz,Meshcounts;
	Meshcounts = b.size()-2;
	Sz = Vlists.size()/Meshcounts;
	for(size_t i=0;i<Meshcounts;++i)
	{
		temp[0] = 0.0;
		temp[1] = 0.0;
		for(size_t j=0;j<Sz;++j)
		{
			//cout<<Sz<<", "<<i+3<<", "<<aaPos[Vlists[i*Sz+j]]<<endl;
			//	load indices associated with vertex
			dd->inner_algebra_indices(Vlists[i*Sz+j], ind);
			u[ind[0]]=b[i+2];
		}
	}


}
template<class TGridFunction,typename TAlgebra>
void RefinedSubsetsTet1Dto3D(TGridFunction &u,TAlgebra &vec,vector<Vertex*> Vlists)
{
	//the 1d to 3d map-----------------------------------------------------
	Refined_Domain1dTo3d(u,vec,u.domain(),u.dof_distribution(),Vlists);

}

//-----------------------------------------------------------------------------
// assign subsets for refined domain in multigrid 3d tetrahendron mesh
//-----------------------------------------------------------------------------
template <typename TDomain, typename TGridFunction, typename TAlgebra>
map<Vertex*, int> Assign_Subsets_Refined_Domain(TGridFunction &u,TAlgebra &vec,SmartPtr<TDomain> domain,
		SmartPtr<DoFDistribution> dd,int Meshcounts,double Meshbegin,double Meshsize)
		{
	//ug::DenseVector<ug::FixedArray1<double, 2> > temp;
	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();
	//  get marker accessor
	typename TDomain::subset_handler_type* Dsh = domain->subset_handler().get();
	//	iterator
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	//	algebra indices vector
	//map<Vertex*, size_t> umap, Mmap;
	int i, j;
	double tempX;
	map<int, int> umap;
	map<Vertex*, int> Vmap;
	// inserting values by using [] operator
	//rf = 0,1,2...; //refine number
	for(j=0;j<(Meshcounts-1);j++)
	{
		i=round( (Meshbegin+(double)(j+1)*Meshsize)*(double)round(2/Meshsize) );
		umap[i] = j+3;
	}
	int counter = 0;
	iter = dd->begin<Vertex>(SurfaceView::ALL);
	iterEnd = dd->end<Vertex>(SurfaceView::ALL);
	for(;iter != iterEnd; ++iter)
	{
		if(Dsh->get_subset_index(*iter)!=1)
		{
			counter+=1;
			Vertex* V = *iter;
			tempX = aaPos[*iter].x()*(double)round(2/Meshsize);
			Vmap[V]=umap[round(tempX)];
			//cout << counter<<", "<<aaPos[V] << ", "<<Vmap[V] <<endl;
		}
	}
	return Vmap;
		}
//----------------------------------------------------------------------------------
template<class TGridFunction,typename TAlgebra>
vector<Vertex*> RefinedSubsetsTet(TGridFunction &u,TAlgebra &vec,size_t Meshcounts,double Meshbegin,double Meshsize,
		vector<Vertex*> Vlists,vector<size_t> Vmarkers)
		{
	map<Vertex*, int> Vmap;
	Vmap = Assign_Subsets_Refined_Domain(u,vec,u.domain(),u.dof_distribution(),Meshcounts,Meshbegin,Meshsize);
	//-----------------------------------------------------
	vector<Vertex*> Vlis[Meshcounts-1];
	// Vlists[0] is corresponding to marker 3 -------------
	map<Vertex*, int>::iterator it;
	for(it = Vmap.begin();it != Vmap.end();++it)
	{
		Vlis[it->second-3].push_back(it->first);

	}
	size_t counter=0;
	//Vmarkers.push_back(0);
	for(size_t i=0;i<Meshcounts-1;++i)
	{
		for(size_t j=0;j<Vlis[i].size();++j)
		{
			//cout<<i+3<<", "<<aaPos[Vlis[i][j]]<<endl;
			Vlists.push_back(Vlis[i][j]);
			//cout<<i+3<<", "<<Vlists[counter]<<endl;
			counter+=1;
		}
		//Vmarkers.push_back(counter);
	}
	//-------------------------------------------------------------------
	//Refined_Domain1dTo3d(u,vec,u.domain(),u.dof_distribution(),Vlists);
	return Vlists;
		}

template<class TGridFunction,typename TAlgebra>
vector<size_t> RefinedSubsetsTet2(TGridFunction &u,TAlgebra &vec,size_t Meshcounts,double Meshbegin,double Meshsize,
		vector<Vertex*> Vlists,vector<size_t> Vmarkers)
		{
	map<Vertex*, int> Vmap;
	Vmap = Assign_Subsets_Refined_Domain(u,vec,u.domain(),u.dof_distribution(),Meshcounts,Meshbegin,Meshsize);
	//-----------------------------------------------------
	vector<Vertex*> Vlis[Meshcounts-1];
	// Vlists[0] is corresponding to marker 3 -------------
	map<Vertex*, int>::iterator it;
	for(it = Vmap.begin();it != Vmap.end();++it)
	{
		Vlis[it->second-3].push_back(it->first);

	}
	size_t counter=0;
	Vmarkers.push_back(0);
	for(size_t i=0;i<Meshcounts-1;++i)
	{
		for(size_t j=0;j<Vlis[i].size();++j)
		{
			//cout<<i+3<<", "<<aaPos[Vlis[i][j]]<<endl;
			//Vlists.push_back(Vlis[i][j]);
			//cout<<i+3<<", "<<Vlists[counter]<<endl;
			counter+=1;
		}
		Vmarkers.push_back(counter);
	}
	//-------------------------------------------------------------------
	//Refined_Domain3dTo1d(u,vec,u.domain(),u.dof_distribution(),Vlists);
	return Vmarkers;
		}
//-----------------------------------------------------------------------------
// assign subsets end----------------------------------------------------------
//-----------------------------------------------------------------------------
template<typename TAlgebra>
void dispm(TAlgebra &M){
	M.print();
}

template<typename TAlgebra>
void CompareMatrix(TAlgebra &M,TAlgebra &N)
{
	size_t Mr=M.num_rows();
	size_t Nr=N.num_rows();
	if(Mr>=Nr)
	{
		for(size_t i=0;i<Nr;i++)
		{
			for(size_t j=0;j<Nr;j++)
			{
				if(M(i+Mr-Nr,j+Mr-Nr)-N(i,j) != 0)
					cout<<"("<<i+Mr-Nr<<","<<j+Mr-Nr<<"), M: "<<M(i+Mr-Nr,j+Mr-Nr)<<" ,N: "<<N(i,j) <<endl;
			}
		}
	}
	else
	{
		for(size_t i=0;i<Mr;i++)
		{
			for(size_t j=0;j<Mr;j++)
			{
				if(N(i-Mr+Nr,j-Mr+Nr)-M(i,j) != 0)
					cout<<N(i-Mr+Nr,j-Mr+Nr)-M(i,j) <<endl;
			}
		}
	}
}

template<typename TAlgebra>
void AssignMatrix(TAlgebra &M,TAlgebra &N)
{
	size_t Mr=M.num_rows();
	size_t Nr=N.num_rows();
	if(Mr>=Nr)
	{
		for(size_t i=0;i<Nr;i++)
		{
			for(size_t j=0;j<Nr;j++)
			{
				if(N(i,j) != 0)
					M(i+Mr-Nr,j+Mr-Nr)=N(i,j);
			}
		}
	}
	else
	{
		for(size_t i=0;i<Mr;i++)
		{
			for(size_t j=0;j<Mr;j++)
			{
				if(M(i,j) != 0)
					N(i-Mr+Nr,j-Mr+Nr)=M(i,j);
			}
		}
	}
}

template<typename TGridFunction, typename TAlgebra>
void AssignVectorNon0(TGridFunction& vec,TAlgebra &b)
{
	//PROFILE_FUNC();
	int gridsize, size2;
	gridsize = (int)vec.size();
	size2 = (int)b.size();
	if(gridsize>=size2)
	{
		//cout<<"u: "<<endl;
		for(int i=2; i<size2; i++)
		{
			vec[i+gridsize-size2] = b[i];
			//cout<<i+gridsize-size2<<", "<<vec[i+gridsize-size2]<<"; "<<b[i]<<endl; //test
		}
	}
	else
	{
		vec[0]=0;
		vec[1]=0;
		for(int i=2; i<gridsize; i++)
		{
			vec[i] = b[i-gridsize+size2];
		}
	}
}
template<typename TGridFunction>
void AssignVectorNon01(TGridFunction& vec,double value,size_t i)
{
	vec[i] = value;
}

template<typename TGridFunction, typename TAlgebra>
void Assign22(TGridFunction& vec,TAlgebra &b){
	//PROFILE_FUNC();
	int gridsize, size2;
	gridsize = (int)vec.size();
	size2 = (int)b.size();
	if(gridsize>=size2)
	{
		//cout<<"u: "<<endl;
		for(int i=0; i<size2; i++)
		{
			vec[i+gridsize-size2] = b[i];
			//cout<<vec[i+gridsize-size2]<<"; "<<b[i]<<endl; //test
		}
	}
	/*else
	{
		cout<<'\n'<<"Assign22: dimension not match"<<endl;
	}*/
}

template<typename TGridFunction, typename TAlgebra>
void Compare22(TGridFunction& vec,TAlgebra &b){
	//PROFILE_FUNC();
	int gridsize, size2;
	gridsize = (int)vec.size();
	size2 = (int)b.size();
	if(gridsize>=size2)
	{
		for(int i=0; i<size2; i++)
		{
			if(vec[i+gridsize-size2]-b[i]!=0)
			{
				cout<<i<<": "<<vec[i+gridsize-size2]<<"; "<<b[i]<<". err: "<<vec[i+gridsize-size2]-b[i]<<endl;
			}//test
		}
	}

	/*else
	{
		cout<<'\n'<<"Assign22: dimension not match"<<endl;
	}*/
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

		{
			typedef TemplateSampleClass<TDomain, TAlgebra> T;
			string name = string("TemplateSampleClass").append(suffix);
			reg.add_class_<T>(name, grp)
	 							.add_constructor()
								.add_method("print_hello",  &T::print_hello, "", "", "prints hello")
								.add_method("get_Vlists",   &T::get_Vlists, "", "", "get Vlists")
								.add_method("get_Vmarkers", &T::get_Vmarkers, "", "", "markers")
								//.add_method("get_Vtest", &T::get_Vtest, "", "", "tests")
								.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "TemplateSampleClass", tag);
		}

		//	WriteGridToVTK
		{
			reg.add_function("LoadVector22", &LoadVector22<function_type>, grp);
			reg.add_function("AssignVector", &Assign22<function_type, vector_type>, grp);
			reg.add_function("CompareVector", &Compare22<function_type, vector_type>, grp);
			reg.add_function("Extract3dto1d", &ExtractPositions3dto1d<function_type, vector_type>, grp);
			reg.add_function("Prolongation1dto3d", &ExtractPositions1dto3d<function_type, vector_type>, grp);
			reg.add_function("Restriction", &Restriction3dto1d<function_type>, grp);
			reg.add_function("Matrix3dto1d", &Matrix3dto1d<function_type, matrix_type>, grp);
			reg.add_function("Matrix3dto1dNon0", &Matrix3dto1dNon0<function_type, matrix_type>, grp);
			reg.add_function("Subsets_Non0_Boundary",&Subsets_Non0_Boundary<function_type,vector_type>, grp);
			reg.add_function("AssignVectorNon0", &AssignVectorNon0<function_type, vector_type>, grp);
			reg.add_function("AssignVectorNon0", &AssignVectorNon01<function_type>, grp);

			reg.add_function("CoarseM", &CoarseM<function_type,matrix_type>, grp);
			reg.add_function("CoarseV", &CoarseV<function_type,vector_type>, grp);
			reg.add_function("ComV", &ComV<function_type,vector_type>, grp);
			reg.add_function("RefinedSubsetsTet", &RefinedSubsetsTet<function_type,vector_type>, grp);
			reg.add_function("RefinedSubsetsTet2", &RefinedSubsetsTet2<function_type,vector_type>, grp);
			reg.add_function("RefinedSubsetsTet3Dto1D",&RefinedSubsetsTet3Dto1D<function_type,vector_type>, grp);
			reg.add_function("RefinedSubsetsTet1Dto3D",&RefinedSubsetsTet1Dto3D<function_type,vector_type>, grp);
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
		typedef typename TAlgebra::matrix_type matrix_type;

		reg.add_function("disp", &disp<vector_type>, grp);
		reg.add_function("disp", &dispm<matrix_type>, grp);
		reg.add_function("CompareMatrix", &CompareMatrix<matrix_type>, grp);
		reg.add_function("AssignMatrix",  &AssignMatrix<matrix_type>, grp);
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
