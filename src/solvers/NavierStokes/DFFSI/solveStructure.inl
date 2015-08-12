/***************************************************************************//**
 * \file solveStructure.inl
 * \author Christopher Minar
 * \brief Implementation of the methods of the class \c DFFSI
 *        to solve the structure equation on the bodies and solve for the new position
 */

#include <solvers/NavierStokes/kernels/fluidStructureInteraction.h>

template <typename memoryType>
void DFFSI<memoryType>::solveStructure()
{
	NavierStokesSolver<memoryType>::logger.startTimer("solveStructure");

	parameterDB  &db = *NavierStokesSolver<memoryType>::paramDB;

	real forcey,
	     Mred,
	     Ured,
	     Cy,
	     dt,
	     alpha_,
	     NumP,
	     tol;

	forcey = NSWithBody<memoryType>::B.forceY[0];
	Mred = 2;
	Ured = 3;
	Cy = forcey*2;
	dt  = db["simulation"]["dt"].get<real>();
	alpha_ = 1.0;
	NumP = NSWithBody<memoryType>::B.numPoints[0]; //only looks at the first immeresed body
	tol = 0.000001;

	real *vb_old        = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.vBk[0]),
	     *vb_iterator   = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.vB[0]),
	     *y_iterator    = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.y[0]),
	     *y_old         = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.yk[0]),
	     *y_updated     = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.ykp1[0]),
	     *test_r        = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.test[0]);

	bool *con_r  = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.converged[0]);

	//call fsi kernel
	dim3 dimGrid(1,1);
	dim3 dimBlock(NumP,1);
	kernels::vorticityInducedVibrationsSC<<<dimGrid,dimBlock>>>(vb_iterator, vb_old, y_old, y_updated, forcey, Mred, Ured, Cy, dt, alpha_);

	printFSI();
	//check for convergence
	kernels::checkConvergencePositionDF<<<dimGrid,dimBlock>>>(tol, con_r, y_iterator, y_updated);

	NavierStokesSolver<memoryType>::logger.stopTimer("solveStructure");
}
