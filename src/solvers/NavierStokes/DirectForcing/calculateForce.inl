/***************************************************************************//**
 * \file calculateForce.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver
 *        to calculate forces acting on each immersed body.
 */

#include <solvers/NavierStokes/kernels/calculateForce.h>

template <>
void DirectForcingSolver<device_memory>::calculateForce()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	real dt = (*paramDB)["simulation"]["dt"].get<real>(),
	     nu = (*paramDB)["flow"]["nu"].get<real>();

	real *q_r      = thrust::raw_pointer_cast(&q[0]),
		 *qOld_r   = thrust::raw_pointer_cast(&qOld[0]),
		 *lambda_r = thrust::raw_pointer_cast(&lambda[0]),
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
void DirectForcingSolver<host_memory>::calculateForce()
{
}

template<typename memoryType>
void DirectForcingSolver<memoryType>::calcForce()
{
	real 	*q_is_r      = thrust::raw_pointer_cast( &(q_is[0]) ),
		*A_is_r      = thrust::raw_pointer_cast( &(A_is[0]) ),
		*rhs1_r      = thrust::raw_pointer_cast( &(NavierStokesSolver<memoryType>::rhs1[0]) ),
		*lambda_is_r = thrust::raw_pointer_cast( &(lambda_is[0]) ),
		*coef1_r     = thrust::raw_pointer_cast( &(coeffsD[0]) ),
		*coef2_r     = thrust::raw_pointer_cast( &(coeffs2D[0]) ),
		*dx_r        = thrust::raw_pointer_cast( &(NavierStokesSolver<memoryType>::domInfo->dx[0]) ),
		*dy_r        = thrust::raw_pointer_cast( &(NavierStokesSolver<memoryType>::domInfo->dy[0]) );

	int 	*tags_r = thrust::raw_pointer_cast( &(tagsD[0]) );

	int	numX = 0,
	   	numY = 0,
		totalPoints = NSWithBody<memoryType>::B.totalPoints,
		nx = NavierStokesSolver<memoryType>::domInfo->nx,
         	ny = NavierStokesSolver<memoryType>::domInfo->ny,
		numU = (nx-1)*ny;
	/*
	vecH	*xuGrid = NavierStokesSolver<memoryType>::domInfo->xu,
		*yuGrid = NavierStokesSolver<memoryType>::domInfo->yu,
		*xvGrid = NavierStokesSolver<memoryType>::domInfo->xv,
		*yvGrid = NavierStokesSolver<memoryType>::domInfo->yv;
	real 	yMax,
		xMax,
		yMin,
		xMin,
		radius = 0.5;
	cusp::array1d<real, memoryType>
		tempForceX,
		tempForceY;
	*/

	//set q_is to q
	q_is = NavierStokesSolver<memoryType>::qOld;
	//make A for qkp1
	NavierStokesSolver<memoryType>::generateA(NavierStokesSolver<memoryType>::intgSchm.alphaImplicit[NavierStokesSolver<memoryType>::subStep]);
	//remake rn without tagged 0 values
	NavierStokesSolver<memoryType>::assembleRHS1();	
	/*looks like this isn't needed
	//calculate vector size
	yMax = NSWithBody<memoryType>::B.y[0] + radius;
	xMax = NSWithBody<memoryType>::B.x[0];
	yMin = NSWithBody<memoryType>::B.y[0] - radius;
	xMin = NSWithBody<memoryType>::B.x[0] - 2*radius; //this doesn't work if the body is rotating...
	//find number of intersections in x
	for (int i=0; i < xuGrid.size; i ++)
	{
		if (xuGrid[i] > xMin && xuGrid[i] < xMax)
			numX ++;
	}
	//find number of intersections in y
	for (int i=0; i < yuGrid.size; i ++)
	{
		if (yvGrid[i] > yMin && yvGrid[i] < yMax)
			numY ++;
	}

	//make vectors to hold temporary forces
	tempForceX.resize(numX);
	tempForceY.resize(numY);
	*/

	//interpolate vB, q, rn, Gphi and Aqlp1 to the surface
	const int blocksize = 256;
	dim3 dimGrid( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 dimBlock(blocksize, 1);

	//q
	//qintersec = q[i]/(1-coef1[i]-coef2[i])
	kernels::qinterpX<<<dimGrid,dimBlock>>>(q_is_r, coef1_r, coef2_r, tags_r, nx, ny);
	kernels::qinterpY<<<dimGrid,dimBlock>>>(q_is_r, coef1_r, coef2_r, tags_r, nx, ny);

	//A
	//Aintersectionu[I] = 0.5*(A[i] + A[i-1])/(1-coef1[I]-coef2[I])
	cusp::multiply(NavierStokesSolver<memoryType>::A, NavierStokesSolver<memoryType>::q, A_is);
	kernels::interpX<<<dimGrid,dimBlock>>>(A_is_r, dx_r, coef1_r, coef2_r, tags_r, nx, ny);		//x
	kernels::interpY<<<dimGrid,dimBlock>>>(A_is_r, dy_r, coef1_r, coef2_r, tags_r, nx, ny);		//y
	//RHS
	//rhsintersectionx[I] = 0.5*(RHS[I] + RHS[I-1])/(1-coef1[I]-coef2[I])
	kernels::interpX<<<dimGrid,dimBlock>>>(rhs1_r, dx_r, coef1_r, coef2_r, tags_r, nx, ny);		//x
	kernels::interpY<<<dimGrid,dimBlock>>>(rhs1_r, dy_r, coef1_r, coef2_r, tags_r, nx, ny);		//y
	//lambda
	//lambdaintersectionx[I] = 0.5*(lambda[I]+lambda[I-1])/(1-coef1[I]-coef2[I])
	cusp::multiply(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::lambda, lambda_is);
	kernels::interpX<<<dimGrid,dimBlock>>>(lambda_is_r, dx_r, coef1_r, coef2_r, tags_r, nx, ny);		//x
	kernels::interpY<<<dimGrid,dimBlock>>>(lambda_is_r, dy_r, coef1_r, coef2_r, tags_r, nx, ny);		//y

	//solve for the force at the intersections
		//vB - q
		cusp::blas::axpby(uvFlex,       q_is,      tempForce, 1.0, -1.0);
		//... - RHS1
		cusp::blas::axpby(tempForce, NavierStokesSolver<memoryType>::rhs1,   tempForce, 1.0, -1.0);
		//... + lambda
		cusp::blas::axpby(tempForce, lambda_is, tempForce, 1.0,  1.0);
		//... - Aulp1
		cusp::blas::axpby(tempForce, A_is,      tempForce,1.0, -1.0);
	//interpolate the forces at the intersections back to the lagrangian body points ... going to try skipping this for now
	//reduce forces
	NSWithBody<memoryType>::B.forceX[0] = thrust::reduce(tempForce.begin(), tempForce.begin() + numU);
	NSWithBody<memoryType>::B.forceY[0] = thrust::reduce(tempForce.begin()+numU+1, tempForce.end());
	//profit

	//fix generate C so it doesn't break when the cylinder moves too far
	//interpolate vB from lagrangian points to intersections
}
