#pragma once


#include "CGAL/CGAL_Mesh.h"
#include "Utils/GMM_Macros.h"
#include "Parser/Parser.h"
#include <unordered_set>
#include "CGAL/CGAL_Borders.h"
#include "Eigen/SparseCore"
#include "Utils/PardisoLinearSolver.h"


class GIF
{

public:
	GIF();
	~GIF() {}
	bool run();

private:

	PardisoLinearSolver mPardisoSolver;

	std::string mObjpath;
	Mesh mCgalMesh;
	double mRadiusOfInitialTutte;
	std::vector<std::complex<double>> mInitialUVsOfCones;
	MeshBuffer mMeshBuffer;
	bool mUseCurvatureMetaVertices; // this option adds the vertices with the largest Gaussian curvature to the meta vertices set
	double mPercentCurvatureMeta; // mPrecenteCurvatureMeta is the percent of the boundary that goes to curvature meta vertices
	double mOuterRateTerminationCondition;
	int mSegSize;
	int mNumInternalFaces;
	bool mUseElimination;
	bool mVisBoundary;
	int mOperationMode;
	int mSizeOfSystemVar, mSizeOfScaffoldSystemVar, mNumDOF, mSizeOfMatrix, mNumScaffoldDOF;
	GMMDenseColMatrix mUVs;
	std::vector<Vertex_handle> mConesAndMetaMap; 
	std::vector<Vertex_handle> mBoundaryVerticesMap; 
	std::vector<Vertex_handle> mScaffoldConesAndMetaMap; 

	Eigen::SparseMatrix<double, Eigen::RowMajor> mKKtEigen;//KKT matrix
	GMMSparseRowMatrix mL12_elimination;
	GMMSparseRowMatrix mQ_mat_elimination;
	GMMSparseRowMatrix mSampling_matrix_elimination;
	std::vector<int> mInternalVerticesIndices;
	std::vector<int> mBoundaryVerticesIndices;
	GMMSparseComplexRowMatrix mHarmonicBasisInAllV;//Harmonic basis matrix in the size of mSizeOfSystemVar X mNumDOF - but only rows that are cones or near cones are filled (the others are just 0 rows)
	std::vector<std::complex<double>> mAllBoundaryUVs;
	std::vector<int> mIndexOfHEinSystem;//this array maps from halfedges indices to the index in the system
	std::vector<int> mIndexOfVertexByIndexOfHEinSystem;
	std::vector<int> mIndexOfHEinScaffoldSystem;
	std::vector<int> mIndicesOfMetaVerticesInUVbyGeneralIndex;
	std::vector<int> mIndicesOfAllBorderVerticesInUVbyGeneralIndex;
	std::vector<int> mGeneralIndicesUVbyIndicesOfAllBorderVertices;
	std::vector<std::vector<Halfedge_handle>> mMetaVerticesInBorderByHalfEdge;//for each border - the halfedges corresponding to meta vertices in this border
	std::vector<std::vector<Halfedge_handle> > mScaffoldMetaVerticesInBorderByHalfEdge;
	std::vector<bool> mIsMeta;
	std::vector<bool> mIsNearMeta;
	std::vector<bool> mIsMetaScaffold;
	std::unordered_set<int> mConesAndNearConesMapOfRowsInKKT;
	std::unordered_set<int> mConesVertices;
	std::unordered_set<int> mNearConesVertices;
	std::unordered_set<int> mNearBoundaryVertices;
	std::unordered_set<int> mBoundaryVertices;
	std::vector<int> mNearConesVerticesVec;
	std::vector<int> mNearBoundaryVerticesVec;
	std::vector<int> mMetaIndicesVec;
	std::vector<int> mNonBoundaryVerticesNearConesVerticesVec;
	std::vector<int> mNonMetaBoundaryVerticesNearConesVerticesVec;
	GMMDenseColMatrix mParameters_matrix;


	bool initialize(CGAL_Borders& borders);
	bool loadMesh();
	void getBordersMapAndSetMetaVertices(CGAL_Borders& borders);
	void getScaffoldBordersMapAndSetMetaVertices(CGAL_Borders& borders1, int numVertices);
	void getVerticesThatMustBeMeta(CGAL_Borders& borders);
	void UpdateIndexOfHEinSystem();
	void UpdateIndexOfHEinScaffoldSystem(Mesh& mesh1);
	bool TransferBorderFacesAndVerticesToMatlab();
	bool TransferScaffoldBorderFacesAndVerticesToMatlab(Mesh& mesh1);

	bool runAlgorithm(CGAL_Borders& borders);
	void constructKKTmatrix(CGAL_Borders& borders, int& conesConstraintsStartRow);
	void calculateHarmonicBasisUsingEliminationWithoutMetaVertices(int& conesConstraintsStartRow);
	void calculateHarmonicBasisUsingEliminationWithMetaVertices(int& conesConstraintsStartRow);
	void FillLaplacianInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues);
	void updateTermInLaplacianByHalfEdge(Halfedge_handle h, const double& term, std::vector<Eigen::Triplet<double>>& tripletListValues);
	void updateTermInLaplacianByIndices(int i, int j, const double& term, std::vector<Eigen::Triplet<double>>& tripletListValues);
	void FillMetaVerticesConstraintsInKKT(CGAL_Borders& borders, std::vector<Eigen::Triplet<double>>& tripletListValues, int& rowInKKT);
	void FillConesConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow);
	void SetElementsForPARDISO(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow);
	bool constructHarmonicBasisAndSendToMATLAB(int conesConstraintsStartRow);
	bool calculateHarmonicBasisInPARDISO(int conesConstraintsStartRow);
	bool constructScaffoldHarmonicBasisAndSendToMATLAB();
	bool calculateScaffoldHarmonicBasis();
	bool getTutteInitialValue(CGAL_Borders& borders);
	bool runNewtonWithRemeshedScaffold(int conesConstraintsStartRow, GMMDenseColMatrix& RHS);
	bool getAllUVsByRHS(GMMDenseColMatrix& RHS);
	void updatHalfedgeUVs();

	bool getResult();
	void checkForFoldovers(std::vector<Facet_handle>& flippedTriangles); //returns a vector of flipped triangles
	void fixCotFoldoversWithMeanValue(std::vector<Facet_handle> flippedTriangles);//this function builds a scaffold around the result, glues the scaffold and the result into one mesh, and maps the new mesh to its own boundary by using mean value weights
	void openReportWindow(bool isFoldoversFree);
};

