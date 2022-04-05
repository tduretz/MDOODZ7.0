#include "cholmod.h"
#include "cs.h"
#include "mdoodz.h"
//---------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ MACROS DEFINITIONS -------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------//
#define MINV(a,b) ((a)<=(b)? (a):(b))
#define MAXV(a,b) ((a)>=(b)? (a):(b))
#define ABSV(a)   ((a)>=0? (a):-(a))
#define _TRUE_  1
#define _FALSE_ 0
#define DoodzFP double

//---------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------- STRUCTURE DEFINITIONS -----------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------//

// Stucture scale contains Stokes sparse matrix
typedef struct _OutputSparseMatrix OutputSparseMatrix ;
struct _OutputSparseMatrix  {
    double params[4];
    double *V, *eta_cell, *b, *F;
    int *Ic, *J, *eqn_u, *eqn_v, *eqn_p;
};

// Contains discrete systems of equation
typedef struct _SparseMat SparseMat;
struct _SparseMat {
	double   *A, *x, *b, *F, *d, *bbc;
	int      *Ic, *J, neq;
    int   *eqn_u, *eqn_v, *eqn_p;
    int   nnz, neq_mom, neq_cont;
};

// Nparams contains numerical parameters of the non-linear solver
typedef struct _n_params Nparams;
struct _n_params {
	int    nit, nit_max, stagnated;
    double abs_tol_u, rel_tol_u, abs_tol_p, rel_tol_p;
	double resx, resz, resp, rest;
    double resx0, resz0, resp0;
    double resx_f, resz_f, resp_f;
	double vrlx,  prlx, trlx;
    int    Picard2Newton, let_res_grow, nit_Pic_max, *LogIsNewtonStep;
    double Pic2NewtCond;
    double *rx_abs, *rz_abs, *rp_abs, *rx_rel, *rz_rel, *rp_rel;
};

// Contains information needed for the direct solver
typedef struct _DirectSolver DirectSolver;
struct _DirectSolver {
    int mtype;           /* Real unsymmetric matrix */
    double res, res0;
    int nrhs;             /* Number of right hand sides */
    int maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    int i, j;
    double ddum;                  /* Double dummy */
    int idum;                 /* Integer dummy */
    char *uplo;
    int n_th;
    int flag;
    int     Analyze;
    cholmod_factor *Lfact;
    cholmod_common c;
};

//---------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ FUNCTION PROTOTYPES ------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------//

// Miscellaneous functions
void Initialise1DArrayDouble( double*, int, double );

// Input/Output
void ScaleMe( scale* );
void    ReadInputFile( char*, int*, int*, int*, int*, params*, scale*, mat_prop*, markers*, Nparams* );
void    UpdateInputFile( char[], int );
int     ReadInt2( FILE*, char[], int );
double  ReadDou2( FILE*, char[], double );
//float   ReadFlo2( FILE*, char[], float );
double  ReadMatProps( FILE* , char[], int, double );
//double  ReadDou21( FILE*, char[], double );
char*   ReadChar( FILE*, char[], char[] );
char*   ReadPhaseDiagram( FILE*, char FieldName[] );
double* ReadBin( char[], int, int, double );
//void    WriteBin(double*, int, const char[]);

// Allocation
void GridAlloc( grid* , params*  );
void GridFree( grid* , params*  );
void PartAlloc( markers*, params*  );
void PartFree( markers*, params* );

// Initialisation
void InitialiseSolutionFields( grid* , params* );
void GridIndices( grid* );
void SetGridCoordinates( grid*, params*, int, int );

// Memory
void* DoodzMalloc(  size_t );
void* DoodzCalloc( int, size_t );
void* DoodzRealloc( void*, size_t );
void  DoodzFree( void* );
void  AllocMat( SparseMat*, int );
void  FreeMat( SparseMat* );
void  FreeSparseSystems( int, int, SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat* );

// General function prototypes
void ArrayPlusArray( double*, double*, int );
void ArrayMinusArray( double*, double*, int );
void ArrayEqualArray( double*, double*, int );
void ArrayEqualArrayI( int*, int*, int );
//void ArrayPlusScalar( double*, double, int );
void ArrayTimesScalar( double*, double, int );
void ArrayTimesScalarArray( double*, double, double*, int );
void ArrayDividedScalarArray( double*, double, double*, int );
void ArrayPlusScalarArray( double*, double, double*, int );
void Initialise2DArray( double*, int, int, double );
void Initialise2DArrayInt( int*, int, int, int );
void MinMaxArrayVal( DoodzFP*, int, double*, double* );
void MinMaxArrayTag( DoodzFP*, double, int, char*, char* );
void MinMaxArrayTagInt( int*, double, int, char*, char* );
void MinMaxArrayPart( DoodzFP*, double, int, char*, int* ) ;
double SumArray( double*, double, int, char*);
//void MinMaxArrayF( float*, double, int, char*);
void MinMaxArrayI( int*, int, int, char*);
//void SumArrayF( float*, double, int, char*);
//void SumArrayI( int*, double, int, char* );
//void Initialise1DArrayDouble( double*, int, double );
void Initialise1DArrayChar( char*, int, char );
//void Print2DArrayDoubleTag( DoodzFP*, int, int, double, char* );
//void Print2DArrayChar( char*, int, int, double);
void Initialise1DArrayInt( int*, int, int );
void IsNanArray2DFP( DoodzFP*, int );
void IsInfArray2DFP( DoodzFP*, int );
void InterpCentroidsToVerticesDouble( double*, double*, grid*, params* );
void InterpVerticesToCentroidsDouble( double*, double*, grid*, params* );
//
//
// Grid initialisation and boundary conditions
//void MInitialiseSolutionFields( grid*, params* );
//void MinitMG( grid*, paramsOutputSparseMatrix  );
void SetBCs_new( grid*, params*, scale, markers*, mat_prop* );
void SetBCs_user( grid*, params*, scale, markers*, mat_prop* );
//void MSetRes( grid*, paramsOutputSparseMatrix  );
//void gridFillLevels( grid*, paramsOutputSparseMatrix  );
//void gridSetCoords( grid*, params*OutputSparseMatrix );
//void MSetGrid( grid*, params, int, int, int);
//void eval_anal_Dani( double*, double*, double*, double*, double*, double*, double, double, int, double, double, double );
void ComputeLithostaticPressure( grid*, params*, double, scale, int );

//// Particles
void PutPartInBox( markers*, grid*, params, surface, scale );
void PartInit( markers*, params* );
void Interp_P2U( markers, DoodzFP*, grid*, double*, double*, double*, int, int, int, char* , params* );
void Interp_P2N( markers, DoodzFP*, grid*, double*, double*, double*, int, int, params* );
void Interp_P2C( markers, DoodzFP*, grid*, double*, double*, double*, int, int );
void Interp_Grid2P( markers, DoodzFP*, grid*, double*, double*, double*, int, int, char* );
void Interp_Grid2P_strain( markers, DoodzFP*, grid*, double*, double*, double*, int, int, char* );
//void FreeP2Mesh( grid* );
void Interp_Phase2VizGrid( markers, int*, grid*, char*, double*, double*, int, int, params, surface );
void ParticleInflowCheck ( markers*, grid*, params, surface, int);
void P2Mastah ( params*, markers , DoodzFP*, grid*, double*, char*, int, int, int, int, int );
//
//// Stokes
//void ResidualCalc2( grid*OutputSparseMatrix *, params, int, double*, double*, double*, double*, double*, double*, int, scale );
//void Gauss_Seidel( grid*OutputSparseMatrix *, params, scale );
//void SaveRHS0( grid* );
void EvaluateRHS( grid*, params, scale, double );
//void DefineResidualScales( Mparams*, grid* );
//void ResidualNorm( gridOutputSparseMatrix  * );
//void StrainStressCalc( grid* );
//void PressureScaling( Mparams, grid*, mat_prop, Eparams );
//
//// Interpolation fine to coarse (restriction)
//void RestrictAll( grid *, int, int );
//void Interp_F2C_U5( grid*, int, double*, double* );
//void Interp_F2C_V5( grid*, int, double*, double* );
//void Interp_F2C_Peta2( grid*, int, double*, double*, double*, double* );
//void Interp_F2C_C4( grid*, int, double*, double*, int );
//void Interp_F2C_N4( grid*, int, double*, double*, int );
//
//// Interpolation coarse to fine (prolongation)
//void ProlongAll( grid*, intOutputSparseMatrix * );
//void Interp_C2F_C( grid*, int, double*, double* );
//void Interp_C2F_U( grid*, int, double*, double* );
//void Interp_C2F_V( grid*, int, double*, double* );
//
//// Visualisation prototypes
void CreateOutputHDF5( const char[] );
void AddGroupToHDF5( const char[], const char[] );
void AddFieldToGroup( const char[], const char[], const char[], const char, int, void*, int );
//void Myfopen( char*, FILE** );
//void MViz_vtk( grid*, char*OutputSparseMatrix  );
void WriteOutputHDF5( grid*, markers*, surface*, markers*, params, Nparams, char*, mat_prop, scale );
void WriteOutputHDF5Particles( grid*, markers*, surface*, markers*, surface*, markers* , params, char*, mat_prop, scale );
void WriteResiduals( grid, params, Nparams, scale );
//
//
//// Phase diagrams
//void AllocatePhaseDiagrams( params* );
//void FreePhaseDiagrams( params* );
//
// Output
void MakeBreakpointParticles( markers*, grid*, markers*, markers*, params , surface*, surface*, scale );
void LoadBreakpointParticles( markers*, grid*, markers*, markers*, params*, surface*, surface*, scale );
void LoadIniParticles       ( char*, markers*, grid* , markers*, markers*, params*, scale );
void DeletePreviousBreakpoint( int, int );
//
// Direct solver
void EvalNumberOfEquations( grid*, SparseMat* );
void SAlloc( SparseMat*, int );
void SFree( SparseMat* );
void BuildStokesOperator( grid*, params, int, double*, double*, double*, SparseMat*, int );
void EvaluateStokesResidual( SparseMat*, Nparams*, grid*, params, scale, int );
//void StokesDirectSolve( grid*OutputSparseMatrix *, params, int, double*, double*, double*, double*, double*, double* );
//void StokesDirectSolveCoarse( grid*OutputSparseMatrix *, params, int, double*, double*, double*, double*, double*, double* );
//void DirectSolverCall( double*, int*, int*, grid*, int );
void SolveStokes( SparseMat*, DirectSolver* );
void SolveStokesDefect( SparseMat*, DirectSolver*, Nparams*, grid*, params*, markers*, markers*, surface*, mat_prop, scale );
void DirectStokes( SparseMat*, DirectSolver*, double* , double* );
void ExtractSolutions( SparseMat*, grid*, params* );
void InitialiseSolutionVector( grid*, SparseMat*, params* );

// Viscoelastoplasticity
void RotateStresses( grid, markers*, params, scale* );
void UpdateParticleStress( grid*, markers*, params*, mat_prop*, scale* );
void ShearModCompExpGrid( grid*, mat_prop, params, scale );
void CohesionFrictionDilationGrid( grid*, markers*, mat_prop, params, scale  );

// Non-Newtonian rheology
void UpdateNonLinearity( grid*, markers*, markers*, surface*, mat_prop, params*, Nparams*, scale, int, double );
double LineSearch( SparseMat*, double*, grid*, params*, Nparams*, markers*, markers*, surface*, mat_prop, scale );
void NonNewtonianViscosityGrid( grid*, mat_prop*, params*, Nparams, scale* );
void StrainRateComponents( grid*, scale, params* );
void GenerateDeformationMaps( grid*, mat_prop*, params*, Nparams, scale*);
void UpdateParticleGrainSize( grid*, scale, params, markers*, mat_prop* );
void UpdateParticleDensity( grid*, scale, params, markers*, mat_prop* );
void UpdateParticleX( grid*, scale, params, markers*, mat_prop* );
void UpdateParticlePhi( grid*, scale, params, markers*, mat_prop* );
// Advection
void DefineInitialTimestep( params*, grid*, markers, mat_prop, scale );
void EvaluateCourantCriterion( double*, double*, params*, scale, grid*, int);
void Check_dt_for_advection( double*, double*, params*, scale, grid*, int);
void RogerGunther( markers*, params, grid, int, scale );
void isout( markers*, params );
void isoutPart( markers*, params*, int );
void CountPartCell    ( markers*, grid* , params, surface, surface, int, scale );
void CountPartCell_Old( markers*, grid* , params, surface, int, scale );
void CountPartCell2   ( markers*, grid* , params, surface, surface, int, scale );

void AccumulatedStrain( grid*, scale , params, markers* );
void PureShearALE( params*,  grid*, markers*, scale );
void VelocitiesOnCenters( double*, double*, double*, double*, int, int, scale );
void VelocitiesToParticles( grid*, markers*, DoodzFP*, DoodzFP*, params, scale );
void DeformationGradient ( grid, scale, params , markers * );


// Energy
void UpdateParticleEnergy( grid*, scale, params, markers*, mat_prop* );
void EnergyDirectSolve( grid*, params, double*, double*, double*, double*, markers*, double, int, int, scale, int );
cholmod_factor* FactorEnergyCHOLMOD( cholmod_common*, cs_di*, double*, int*, int*, int, int, int );
cs_di* TransposeA( cholmod_common*, double*, int*, int*, int, int );
void SolveEnergyCHOLMOD( cholmod_common*, cs_di*, cholmod_factor*, double*, double*, int, int, int );
void ThermalSteps( grid*, params, double*, double*, double*, double*, markers*, double, scale );
void SetThermalPert( grid*, params, scale );
void UpdateMaxPT ( scale, params, markers* );

// Free surface routines
void SetTopoChainHorizontalCoords( surface*, markers*, params, grid, scale );
void AllocateMarkerChain( surface*, markers*, params );
void FreeMarkerChain( surface*, markers* );
void CellFlagging( grid*, params, surface, scale );
void ProjectTopography( surface*, markers*, params, grid, scale, double*, int );
void MarkerChainPolyFit( surface*, markers*, params, grid );
void CleanUpSurfaceParticles( markers*, grid*, surface, scale );
void RemeshMarkerChain( markers*, surface*, params, scale, grid*, int );
void SurfaceDensityCorrection( grid*, params, surface, scale);
void SurfaceVelocity( grid*, params, surface*, markers*, scale );
void UpdateDensity( grid*, markers*, mat_prop*, params*, scale* );
void DiffuseAlongTopography( grid*, params, scale, double*, double*, int, double, double );
void AddPartSed( markers *, mat_prop , markers *, surface *, params , scale , grid *);
void CorrectTopoIni( markers *, mat_prop , markers *, surface *, params , scale , grid *);

// Decoupled solver
void KillerSolver( SparseMat*,  SparseMat*,  SparseMat*,  SparseMat*, DirectSolver*, double*, double*, double*, params, grid*, scale, SparseMat*, SparseMat*, SparseMat*,  SparseMat*,  SparseMat* );
void KSPStokesDecoupled( SparseMat*,  SparseMat*,  SparseMat*,  SparseMat*, DirectSolver*, double*, double*, double*, params, grid*, scale, SparseMat*, SparseMat*, SparseMat*,  SparseMat*,  SparseMat* );
double LineSearchDecoupled( SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, double*, grid*, params*, Nparams*, markers*, markers*, surface*, mat_prop, scale );
void EvaluateStokesResidualDecoupled( SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, Nparams*, grid*, params, scale, int );
void BuildStokesOperatorDecoupled( grid*, params,int, double*, double*, double*, double*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, int );
void SolveStokesDecoupled( SparseMat*, SparseMat*, SparseMat*,  SparseMat*, SparseMat*, DirectSolver*, params, grid*, scale );
void SolveStokesDefectDecoupled( SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, DirectSolver*, Nparams*, grid*, params*, markers*, markers*, surface*, mat_prop, scale, SparseMat*, SparseMat*, SparseMat* );
void AddCoeff2( int*, double*, int, int, int*, double, int, double, double* );
void MergeParallelMatrix( SparseMat*, double**, int**, int**, grid*, int*, int*, int*, int*, int*, int, char*, int* );
void DirectStokesDecoupled( SparseMat*, SparseMat*, SparseMat*,  SparseMat*, DirectSolver*, double*, double*, double*, params, grid*, scale, SparseMat* );
void DirectStokesDecoupledComp( SparseMat*, SparseMat*, SparseMat*,  SparseMat*, DirectSolver*, double*, double*, double*, params, grid*, scale, SparseMat* );
//void ConvertTo1Based( int*, int );
void DecompressCSRtoTriplets( int, int*, int* );
//void  ArrayEqualScalarArray( DoodzFP*, DoodzFP, DoodzFP*, int );
//void ApplyBC( grid*, params );
//void GridIndices( grid* );
//void ScaleBack          (float*,  double , int );
void ScaleBackD         (double*, double , int );
//void DoubleToFloat      (double*,  float*, int );
//
// Newton
void BuildJacobianOperatorDecoupled( grid*, params, int, double*, double*, double*, double*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, int );
//
//
// Flow law Database
void ReadDataPowerLaw   ( mat_prop*, params*, int, int, scale* );
void ReadDataLinear     ( mat_prop*, params*, int, int, scale* );
void ReadDataGBS        ( mat_prop*, params*, int, int, scale* );
void ReadDataExponential( mat_prop*, params*, int, int, scale* );
void ReadDataGSE        ( mat_prop*, params*, int, int, scale* );
void ReadDataKinetics   ( mat_prop*, params*, int, int, scale* );

void AllocatePhaseDiagrams( params* );
void FreePhaseDiagrams( params* );

double DotProduct( DoodzFP*, DoodzFP*, int  );
void BackToSolutionVector( cholmod_dense*, cholmod_dense*, double*, grid*, SparseMat* );
void NormResidualCholmod( double*, double*, cholmod_dense*, cholmod_dense*, int, int, params, scale, int );
void BuildInitialSolutions( double*, double*, grid* );
void copy_cholmod_to_cs_matrix( cholmod_sparse*, cs* );
void copy_cs_to_cholmod_matrix( cholmod_sparse*, cs* );
void copy_vec_to_cholmod_dense( cholmod_dense*, DoodzFP* );
void copy_cholmod_dense_to_cholmod_dense( cholmod_dense*, cholmod_dense* );
void cholmod_dense_plus_cholmod_dense( cholmod_dense*, cholmod_dense* );

void ApplyBC( grid*, params* );
void AssignMarkerProperties (markers*, int, int, params*, grid*, int );


// GLOBAL
//void Interp_P2G( markers, DoodzFP*, grid*, double*, double*, double*, int, int, double, double, int, int, params*, char*  );

void AdvectFreeSurf_BEN( markers*, params, scale );
void BuildInitialTopography_BEN( surface*, markers*, params, grid, scale );
void SetTopoChainHorizontalCoords_BEN( surface*, markers*, params, grid, scale );
void AllocateMarkerChain_BEN( surface*, markers*, params );
void FreeMarkerChain_BEN( surface*, markers* );
void CellFlagging_BEN( grid*, params, surface, scale );
void ProjectTopography_BEN( surface*, markers*, params, grid, scale, double*, int );
//double TopoFun( double, int, surface, scale );
void MarkerChainPolyFit_BEN( surface*, markers*, params, grid );
void CleanUpSurfaceParticles_BEN( markers*, grid*, surface, scale );
void RemeshMarkerChain_BEN( markers*, surface*, params, scale, grid*, int );
void SurfaceDensityCorrection_BEN( grid*, params, surface, scale);
void SurfaceVelocity_BEN( grid*, params, surface*, markers*, scale );
void UpdateDensity_BEN( grid*, markers*, mat_prop*, params*, scale* );
//void AdvectFreeSurf( markers*, params, scale );
//void PhaseGrowth( markers*, params, grid* );
void DiffuseAlongTopography_BEN( grid*, params, scale, double*, int, double, double );
void AddPartSed_BEN( markers *, mat_prop , markers *, surface *, params , scale , grid *);
void CorrectTopoIni_BEN( markers *, mat_prop , markers *, surface *, params , scale , grid *);
void EvaluateRHS_BEN( grid*, params, scale, double );
void UpdateNonLinearity_BEN( grid*, markers*, markers*, surface*, mat_prop, params*, Nparams*, scale, int, double );
void PressureScaling_BEN(  grid *, mat_prop , params  );
void StrainStressCalc_BEN( grid* mesh );
void UpdateParticleVelocity_BEN( grid*, scale, params, markers* );

void BuildStokesOperator_BEN( grid*, params, int, double*, double*, double*, SparseMat*, int );
int EvalNumberOfEquations_BEN( grid* mesh, SparseMat *Stokes );
void ExtractSolutions_BEN( SparseMat*, grid*, params );
void EvaluateStokesResidual_BEN( SparseMat *, Nparams *, grid *, params , scale , int  );
void EvaluateCourantCriterion_BEN( double* , double* , params *, scale , grid*, int  );
void RogerGunther_BEN( markers *, params , grid , int , scale  );

void Interp_P2N_BEN ( markers , DoodzFP* , grid *, double* , double* , double* , int , int, params*  ) ;
void Interp_P2C_BEN ( markers , DoodzFP* , grid *, double* , double* , double* , int , int   ) ;
void Interp_P2U_BEN ( markers , DoodzFP* , grid *, double* , double* , double* , int , int , int , char*  ) ;
void Interp_Grid2P_BEN ( markers , DoodzFP* , grid *, double* , double* , double* , int , int , char * ) ;
void CountPartCell_BEN ( markers* , grid *, params , surface, int, scale  );
void PutPartInBox_BEN( markers *, grid *, params , surface , scale  );
void SetBCs_BEN( grid *, params *, scale , markers* , mat_prop * );

void SetParticles_BEN( markers *, scale , params , mat_prop*   );
void BuildInitialTopography_BEN( surface *, markers *, params , grid , scale  );
void SolveStokes_BEN( SparseMat*, DirectSolver* );

void V2P( double*, double*, markers*, double*,  double*, double*, double*, double*, double*, int, int, int, int, char*, char*, double, double, int, int );
double Vertices2Particle( markers*, double*, double*, double*, int, int , char *, double, double, int );
double Centers2Particle( markers*, double*, double*, double*, int, int, char*, double, double, int, int );
void RogerGuntherII( markers*, params, grid, int, scale );
void AccumulatedStrainII( grid*, scale, params, markers*, double*, double*, int, int, char * );
void AdvectFreeSurf( markers*, params, scale );

void InitialiseDirectorVector (grid*, markers*, params*, mat_prop*);
void NormalizeDirector ( grid* , DoodzFP*, DoodzFP*, DoodzFP*, DoodzFP*, params*  );
void RotateDirectorVector( grid, markers*, params, scale* );
void UpdateParticlePressure( grid*, scale, params, markers*, mat_prop* );
void DetectCompressibleCells ( grid* , params*  );
void ScaleVelocitiesRHSBack(SparseMat*, SparseMat*, double*);
void ExtractDiagonalScale(SparseMat *, SparseMat *, SparseMat *, SparseMat * );
void ScaleMatrix(SparseMat *, SparseMat *, SparseMat *, SparseMat * ) ;

void RheologicalOperators( grid*, params*, scale*, int );
void ComputeViscosityDerivatives_FD( grid*, mat_prop*, params*, Nparams, scale* );
void SetUpModel_NoMarkers ( grid*, params*, scale* );
void Diffuse_X( grid*, params*, scale* );
void FiniteStrainAspectRatio ( grid*, scale, params, markers* );
void Print2DArrayDouble( DoodzFP*, int, int, double );

void OldDeviatoricStressesPressure( grid*, markers*, scale, params*  );
void TotalStresses( grid*, markers*, scale, params* );

void Interp_P2G ( markers*, DoodzFP*, grid*, double*, double*, double*, int, int, params*, char*, int, int );
void Interp_Grid2P_centroids ( markers, DoodzFP* , grid *, double* , double* , double* , int , int , char *, params* );
void Interp_Grid2P_centroids2( markers, DoodzFP* , grid *, double* , double* , double* , int , int , char *, params* );
void ExpandCentroidArray( double*, double*, grid*, params* );
void ComputeIncrementsOnParticles( grid*, markers*, params*, mat_prop*, scale* );
void UpdateGridFields( grid*, markers*, params*, mat_prop*, scale* );

void RogerGunther( markers *, params, grid, int, scale );
void CheckSym( DoodzFP*, double, int, int, char*, int, int );
void ChemicalDirectSolve( grid*, params, markers*, mat_prop*, double, scale );
void InitialiseGrainSizeParticles( markers*, mat_prop* );

void ViscosityDerivatives( grid*, mat_prop*, params*, Nparams, scale*);
double ViscosityConcise( int , double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, mat_prop*, params*, scale*, double*, double*, double*, double*, double*, double*, double*, double*, double* , double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double, double, double, double, double, double,  double*, double*, double*, double*, double, double, double*, double*, double*, int, int );
double EvaluateDensity( int, double, double, double, params*, mat_prop* );
void ComputeMeanQuantitesForTimeSeries( grid *mesh );
void LogTimeSeries( grid*, params, scale );
