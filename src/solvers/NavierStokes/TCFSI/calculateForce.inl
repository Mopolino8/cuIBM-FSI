/***************************************************************************//**
 * \file calculateForce.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c TCFSISolver
 *        to calculate forces acting on each immersed body.
 */

#include <solvers/NavierStokes/kernels/calculateForce.h>

/**
 * \brief Calculates forces acting on each immersed body.
 *
 * For each immersed body, it spreads the Lagrangian forces on the Eulerian grid
 * and sum the entire field to get the x- and y- components of the force.
 *
 */
template <typename memoryType>
void TCFSISolver<memoryType>::calculateForce()
{
	typedef typename cusp::array1d<real, memoryType>::iterator ValueIterator;
	typedef typename cusp::array1d_view<ValueIterator>         View;

	int     nx = NavierStokesSolver<memoryType>::domInfo->nx,
	        ny = NavierStokesSolver<memoryType>::domInfo->ny,
	        numBodies = NSWithBody<memoryType>::B.numBodies,
	        totalPoints = NSWithBody<memoryType>::B.totalPoints,
	        ETRows = (nx-1)*ny + nx*(ny-1);
	real    dx, dy;
	View    f, fView;

	cusp::array1d<real, memoryType> F(ETRows), fTemp(2*totalPoints);

	dx = NavierStokesSolver<memoryType>::domInfo->dx[ NSWithBody<memoryType>::B.I[0] ],
	dy = dx;

	//f = View(NavierStokesSolver<memoryType>::lambda.begin() + nx*ny, NavierStokesSolver<memoryType>::lambda.end());
	f = View(NavierStokesSolver<memoryType>::lambdak.begin() + nx*ny, NavierStokesSolver<memoryType>::lambdak.end());
	// loop through bodies
	for(int l=0; l < numBodies; l++)
	{
		dx = NavierStokesSolver<memoryType>::domInfo->dx[ NSWithBody<memoryType>::B.I[l]];
		dy = dx;
		// x-component of the force
		fTemp = f;
		fView = View(fTemp.begin(), fTemp.begin() + NSWithBody<memoryType>::B.offsets[l]);
		cusp::blas::fill(fView, 0.0);
		if(l < numBodies-1)
		{
			fView = View(fTemp.begin() + NSWithBody<memoryType>::B.offsets[l+1], fTemp.end());
			cusp::blas::fill(fView, 0.0);
		}
		cusp::multiply(ET, fTemp, F);
		NSWithBody<memoryType>::B.forceX[l] = (dx*dy)/dx*thrust::reduce(F.begin(), F.end());

		// y-component of the force
		fTemp = f;
		fView = View(fTemp.begin(), fTemp.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l]);
		cusp::blas::fill(fView, 0.0);
		if(l < numBodies-1)
		{
			fView = View(fTemp.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l+1], fTemp.end());
			cusp::blas::fill(fView, 0.0);
		}
		cusp::multiply(ET, fTemp, F);
		NSWithBody<memoryType>::B.forceY[l] = (dx*dy)/dy*thrust::reduce(F.begin(), F.end());
	}
}

//not used
/*
template <>
void TCFSISolver<device_memory>::calculateForce2()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	real dt = (*paramDB)["simulation"]["dt"].get<real>(),
	     nu = (*paramDB)["flow"]["nu"].get<real>();
	
	real *q_r      = thrust::raw_pointer_cast(&qk[0]),
		 *qOld_r   = thrust::raw_pointer_cast(&qOld[0]),
		 *lambda_r = thrust::raw_pointer_cast(&lambdak[0]),
		 *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
		 *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	// Calculating drag
	cusp::array1d<real, device_memory>
		FxX(B.numCellsY[0]),
		FxY(B.numCellsX[0]+1),
		FxU((B.numCellsX[0]+1)*B.numCellsY[0]);
	
	real *FxX_r = thrust::raw_pointer_cast(&FxX[0]),
	     *FxY_r = thrust::raw_pointer_cast(&FxY[0]),
	     *FxU_r = thrust::raw_pointer_cast(&FxU[0]);

	const int blockSize = 256;
	dim3 dimGrid( int((B.numCellsX[0]+B.numCellsY[0]+1-0.5)/blockSize)+1, 1 );
	dim3 dimBlock(blockSize, 1);
	
	kernels::dragLeftRight <<<dimGrid, dimBlock>>> (FxX_r, q_r, lambda_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::dragBottomTop <<<dimGrid, dimBlock>>> (FxY_r, q_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridX( int( ( (B.numCellsX[0]+1)*B.numCellsY[0]-0.5 )/blockSize )+1, 1 );
	kernels::dragUnsteady <<<dimGridX, dimBlock>>> (FxU_r, q_r, qOld_r, dxD, dyD, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	B.forceX[0] = thrust::reduce(FxX.begin(), FxX.end()) + thrust::reduce(FxY.begin(), FxY.end()) + thrust::reduce(FxU.begin(), FxU.end());
	
	// Calculating lift
	cusp::array1d<real, device_memory>
		FyX(B.numCellsY[0]+1),
		FyY(B.numCellsX[0]),
		FyU((B.numCellsX[0]+1)*B.numCellsY[0]);
	
	real *FyX_r = thrust::raw_pointer_cast(&FyX[0]),
	     *FyY_r = thrust::raw_pointer_cast(&FyY[0]),
	     *FyU_r = thrust::raw_pointer_cast(&FyU[0]);

	kernels::liftLeftRight <<<dimGrid, dimBlock>>> (FyX_r, q_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::liftBottomTop <<<dimGrid, dimBlock>>> (FyY_r, q_r, lambda_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridY( int( ( B.numCellsX[0]*(B.numCellsY[0]+1)-0.5 )/blockSize )+1, 1 );
	kernels::liftUnsteady <<<dimGridY, dimBlock>>> (FyU_r, q_r, qOld_r, dxD, dyD, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	B.forceY[0] = thrust::reduce(FyX.begin(), FyX.end()) + thrust::reduce(FyY.begin(), FyY.end()) + thrust::reduce(FyU.begin(), FyU.end());
}





template <>
void TCFSISolver<host_memory>::calculateForce2()
{
}*/
