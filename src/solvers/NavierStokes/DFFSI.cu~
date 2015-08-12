/***************************************************************************//**
 * \file DFFSI.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DFFSI.
 */


#include "DFFSI.h"
#include <sys/stat.h>
#include <thrust/extrema.h>
#include <cusp/io/matrix_market.h>
#include <ctime>

/**
 * \brief Constructor. Initializes the simulation parameters and the domain info.
 */
template <typename memoryType>
DFFSI<memoryType>::DFFSI(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

/**
 * \brief Initialize the vectors used in the simulation.
 */
template <typename memoryType>
void DFFSI<memoryType>::initialise()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;

	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	NavierStokesSolver<memoryType>::initialiseCommon();
	
	NSWithBody<memoryType>::initialiseBodies();
	
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP);
	
	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");

	tags.resize(numUV);
	tagsD.resize(numUV);
	tags2.resize(numUV);
	tags2D.resize(numUV);
	coeffs.resize(numUV);
	coeffsD.resize(numUV);
	coeffs2.resize(numUV);
	coeffs2D.resize(numUV);
	uv.resize(numUV);
	uvD.resize(numUV);
	q_is.resize(numUV);
	A_is.resize(numUV);
	lambda_is.resize(numUV);
	tempForce.resize(numUV);

	pressure.resize(numP);
	cusp::blas::fill(pressure, 0.0);

	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
	
	tagPoints();//tagged for update
	std::cout << "Done tagging points!" << std::endl;
	
	NavierStokesSolver<memoryType>::assembleMatrices();
}

template <typename memoryType>
void DFFSI<memoryType>::stepTime()
{
	NavierStokesSolver<memoryType>::qOld = NavierStokesSolver<memoryType>::q;
	NSWithBody<memoryType>::B.yk = NSWithBody<memoryType>::B.y;
	//NSWithBody<memoryType>::B.xk = NSWithBody<memoryType>::B.x;
	NSWithBody<memoryType>::B.vBk = NSWithBody<memoryType>::B.vB;
	//NSWithBody<memoryType>::B.uBk = NSWithBody<memoryType>::B.uB;
	int k = 1; 

	if (NavierStokesSolver<memoryType>::timeStep < 1) //no fsi
	{
		DFFSI<memoryType>::updateSolverState();
		NavierStokesSolver<memoryType>::generateRN();
		NavierStokesSolver<memoryType>::generateBC1();
		DFFSI<memoryType>::assembleRHS1();
		NavierStokesSolver<memoryType>::solveIntermediateVelocity();

		// Set up and solve the Poisson system
		NavierStokesSolver<memoryType>::generateBC2();
		NavierStokesSolver<memoryType>::assembleRHS2();
		NavierStokesSolver<memoryType>::solvePoisson();

		// Projection step
		DFFSI<memoryType>::projectionStep(); //tagged for update

		DFFSI<memoryType>::calculateForce();
		printFSI();
	}
	else //fsi
	{
		do
		{
		//updates required for a moving body
		DFFSI<memoryType>::updateSolverState(); //qt q c pc2

		//intermediate velocity
		NavierStokesSolver<memoryType>::generateRN();
		NavierStokesSolver<memoryType>::generateBC1();
		DFFSI<memoryType>::assembleRHS1();
		NavierStokesSolver<memoryType>::solveIntermediateVelocity();

		// Set up and solve the Poisson system
		NavierStokesSolver<memoryType>::generateBC2();
		NavierStokesSolver<memoryType>::assembleRHS2();
		NavierStokesSolver<memoryType>::solvePoisson();

		// Projection step
		DFFSI<memoryType>::projectionStep();

		//FSI
		DFFSI<memoryType>::calculateForce();
		NSWithBody<memoryType>::B.converged[0] = true;
		DFFSI<memoryType>::solveStructure();

		k++;
		}while(k<1);
		//}while (NSWithBody<memoryType>::B.converged[0] == false);
	}
	NavierStokesSolver<memoryType>::timeStep++;
}

/**
 * \brief Updates the matrices every time the body is moved.
 *
 * Change the location of the body points, re-tags all the points on the
 * velocity grid to locate the new forcing nodes. Reassembles the
 * matrices. Moving bodies have not been tested for DFFSI.
 */
template <typename memoryType>
void DFFSI<memoryType>::updateSolverState()
{
	//if (fabs(NSWithBody<memoryType>::B.y[0]) > 0.001)
	{
	//calculateCellIndicies
	NavierStokesSolver<memoryType>::logger.startTimer("updateSolverState");
	NSWithBody<memoryType>::B.calculateCellIndices(*NavierStokesSolver<memoryType>::domInfo);
	NSWithBody<memoryType>::B.calculateBoundingBoxes(*NavierStokesSolver<memoryType>::paramDB,  *NavierStokesSolver<memoryType>::domInfo);
	// tag points
	NavierStokesSolver<memoryType>::logger.startTimer("updateSolverState_tagPoints");
	tagPoints();
	NavierStokesSolver<memoryType>::logger.stopTimer("updateSolverState_tagPoints");
	
	//assembleMatrices
	//generateL
	generateL();
	//generateA
	generateA(NavierStokesSolver<memoryType>::intgSchm.alphaImplicit[NavierStokesSolver<memoryType>::subStep]);
	//update PC1
	NavierStokesSolver<memoryType>::PC1 ->update(NavierStokesSolver<memoryType>::A);
	//generateQT
	generateQT();
	//generateC
	NavierStokesSolver<memoryType>::logger.startTimer("updateSolverState_generateC");
	generateC(); // QT*BN*Q
	NavierStokesSolver<memoryType>::logger.stopTimer("updateSolverState_generateC");
	//Update PC2
	NavierStokesSolver<memoryType>::PC2 ->update(NavierStokesSolver<memoryType>::C);

	NavierStokesSolver<memoryType>::logger.stopTimer("updateSolverState");
	}
}

template<typename memoryType>
void DFFSI<memoryType>::printFSI()
{
	std::cout<<NavierStokesSolver<memoryType>::timeStep*0.01<<"\t"<<NSWithBody<memoryType>::B.ykp1[0]<<"\t"<<NSWithBody<memoryType>::B.forceY[0]<<"\n";
	NSWithBody<memoryType>::output<<NavierStokesSolver<memoryType>::timeStep*0.01<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.ykp1[0]<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.y[0]<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.yk[0]<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.forceY[0]<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.forceX[0]<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.vBk[0]<<"   \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.vB[0]<<"   \n";
}

/**
 * \brief Assembles the matrix rhs1 for DFFSI.
 *
 * This function first calls the function assembleRHS1 from NavierStokesSolver.
 * Then it called the function updateRHS1 to modify only the elements of the 
 * vector that correspond to the interpolation nodes.
 */
template <typename memoryType>
void DFFSI<memoryType>::assembleRHS1()
{
	NavierStokesSolver<memoryType>::assembleRHS1();
	
	NavierStokesSolver<memoryType>::logger.startTimer("updateRHS1");
	updateRHS1();
	NavierStokesSolver<memoryType>::logger.startTimer("updateRHS1");
}

/**
 * \brief Prints the min, max and sum of the divergences of the velocity field 
 *        in every cell of the domain.
 *
 * The divergence is calculated as QT*q, which is technically the sum of the 
 * mass fluxes in every cell. This QT is also differs depending on the Direct 
 * Forcing method used.
 */
template <typename memoryType>
void DFFSI<memoryType>::writeMassFluxInfo()
{
	parameterDB  &db = *NavierStokesSolver<memoryType>::paramDB;
	int     nx = NavierStokesSolver<memoryType>::domInfo->nx,
	        ny = NavierStokesSolver<memoryType>::domInfo->ny,
	        timeStep = NavierStokesSolver<memoryType>::timeStep;

	cusp::array1d<real, memoryType> fluxes(nx*ny);
	cusp::multiply(NavierStokesSolver<memoryType>::QT, NavierStokesSolver<memoryType>::q, fluxes);
	int minPosition = thrust::min_element(fluxes.begin(), fluxes.end()) - fluxes.begin(),
	    maxPosition = thrust::max_element(fluxes.begin(), fluxes.end()) - fluxes.begin();
	real minFlux = fluxes[minPosition],
	     maxFlux = fluxes[maxPosition],
	     globalSum = thrust::reduce(fluxes.begin(), fluxes.end());

	std::ofstream fluxInfoFile;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/massFlux";
	
	if(timeStep==1)
		fluxInfoFile.open(out.str().c_str());
	else
		fluxInfoFile.open(out.str().c_str(), std::ios::out | std::ios::app);
		
	fluxInfoFile << timeStep << '\t' << minFlux << '\t' << maxFlux << '\t' << globalSum << std::endl;
	fluxInfoFile.close();
}

/**
 * \brief Projects the pressure gradient on to the intermediate velocity field
 *        to obtain the divergence-free velocity field at the next time step.
 */
template <typename memoryType>
void DFFSI<memoryType>::projectionStep()
{
	NavierStokesSolver<memoryType>::projectionStep();

	NavierStokesSolver<memoryType>::logger.startTimer("projectionStep");
	cusp::blas::axpy(NavierStokesSolver<memoryType>::lambda, pressure , 1.0);
	NavierStokesSolver<memoryType>::logger.stopTimer("projectionStep");
}

/**
 * \brief Writes the velocity, pressure, force and mass flux data at every save point.
 */
template <typename memoryType>
void DFFSI<memoryType>::writeData()
{	
	NavierStokesSolver<memoryType>::logger.startTimer("output");

	parameterDB  &db = *NavierStokesSolver<memoryType>::paramDB;
	real         dt  = db["simulation"]["dt"].get<real>();
	int          timeStep = NavierStokesSolver<memoryType>::timeStep;

	NSWithBody<memoryType>::writeCommon();
	
	// Print forces calculated using the CV approach
	//NSWithBody<memoryType>::calculateForce();
	NSWithBody<memoryType>::forceFile << timeStep*dt << '\t' << NSWithBody<memoryType>::forceX << '\t' << NSWithBody<memoryType>::forceY << std::endl;

	writeMassFluxInfo();
	
	NavierStokesSolver<memoryType>::logger.stopTimer("output");
}

/**
 * \brief Generates the right-hand side matrix in the Poisson step.
 *
 * Because the fully discrete direct forcing method separates the domains 
 * inside and outside the immersed boundary, and Neumann boundary conditions
 * are enforced at the immersed boundary, a point also needs to be fixed 
 * inside the immersed boundary, so that Poisson system does not have any 
 * zero eigenvalues.
 *
 * In this function, \a phi at a point that corresponds to the center of the 
 * grid is fixed as zero. This is the cell with indices (nx/2, ny/2), and is 
 * the center in terms of the indices and not the physical location in space.
 * Of course, this is invalid if the interior of the body does not include 
 * this point, and the solution obtained will be unphysical.
 *
 * This needs to be done in a better way (i.e. by locating points that 
 * are inside the immersed boundaries, and fixing them to zero)
 */
template <typename memoryType>
void DFFSI<memoryType>::generateC()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
	    ny = NavierStokesSolver<memoryType>::domInfo->ny;
	int index = 5*(ny/2)*nx - nx - ny + 5*(nx/2) - 1 + 2;
	int row = (ny/2)*nx+nx/2;

	NavierStokesSolver<memoryType>::generateC();
	bool flag = true;
	while(flag)
	{
		if(NavierStokesSolver<memoryType>::C.row_indices[index]==NavierStokesSolver<memoryType>::C.column_indices[index] && NavierStokesSolver<memoryType>::C.column_indices[index]==row)
		{
			NavierStokesSolver<memoryType>::C.values[index] += NavierStokesSolver<memoryType>::C.values[index];
			flag = false;
		}
		index++;
	}
}

// inline files in the folder "DFFSI"
#include "DFFSI/tagPoints.inl"
#include "DFFSI/generateL.inl"
#include "DFFSI/generateA.inl"
#include "DFFSI/updateRHS1.inl"
#include "DFFSI/generateQT.inl"
#include "DFFSI/calculateForce.inl"
#include "DFFSI/solveStructure.inl"

// specialization of the class
template class DFFSI<host_memory>;
template class DFFSI<device_memory>;
