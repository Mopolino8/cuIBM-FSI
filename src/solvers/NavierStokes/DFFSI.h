/***************************************************************************//**
 * \file  DFFSI.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c DFFSI.
 */


#pragma once

#include "NSWithBody.h"


/**
 * \class DFFSI
 * \brief A fully discrete formulation of the direct forcing method.
 *
 * Based on the method first proposed by Fadlun et al (2000)
 * with modifications by Kim et al (2001).
 * It does not follow the same equations that they used, but use a 
 * fractional step method starting with the discretized equations,
 * based on the idea by Perot (1993).
 *
 * This method does not use an explicit pressure term in the step where
 * the intermediate velocity is calculated, and the pressure is directly 
 * obtained from the Poisson equation.
 */
template <typename memoryType>
class DFFSI : public NSWithBody<memoryType>
{
protected:
	// For the 1D interpolations
	cusp::array1d<int, host_memory>    tags,  tags2;
	cusp::array1d<int, device_memory>  tagsD, tags2D;
	cusp::array1d<real, device_memory>  q_is, A_is, lambda_is, tempForce;
	vecH coeffs, coeffs2;
	vecD coeffsD, coeffs2D;
	vecH uv;
	vecD uvD;
	cusp::array1d<real, memoryType> pressure;
	
	void tagPoints();
	
	//void tagPoints(real *bx, real *by);
	
	void tagPoints(real *bx, real *by, real *uB, real *vB);

	virtual void generateA(real alpha);

	virtual void generateL();
	void generateQT(int *QTRows, int *QTCols, real *QTVals){} // check this!
	virtual void generateQT();
	void updateQ();

	virtual void generateC();
	
	virtual void assembleRHS1();
	void updateRHS1();

	virtual void updateSolverState();
	virtual void projectionStep();

	void writeMassFluxInfo();
	
public:
	DFFSI(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual void initialise();
	virtual void writeData();

	void printFSI();

	virtual void calculateForce();
	void calcForce();
	void solveStructure();

	//time stepping with fsi fun!
	virtual void stepTime();
	
	/**
	 * \brief Returns the name of the solver as a string.
	 */
	virtual std::string name()
	{
		return "Direct Forcing";
	}
};
