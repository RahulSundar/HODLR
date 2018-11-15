#ifndef __HODLR_Tree__
#define __HODLR_Tree__

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <vector>

#include "HODLR_Matrix.hpp"
#include "HODLR_Node.hpp"

class HODLR_Tree 
{
private:
	int N;
	int n_levels;
	double tolerance;
	std::vector<int> nodes_in_level;
	HODLR_Matrix* A;
	
    // Vector of levels(which contain nodes) thereby giving the tree:
    std::vector<std::vector<HODLR_Node*>> tree;
	void createTree();
	void createRoot();
	void createChildren(int j, int k);

	// Variables and methods needed for HODLR solver
	void factorizeLeaf(int k);
	void factorizeNonLeaf(int j, int k);
	Eigen::MatrixXd solveLeaf(int k, Eigen::MatrixXd b);
	Eigen::MatrixXd solveNonLeaf(int j, int k, Eigen::MatrixXd b);

public:
	HODLR_Tree(int n_levels, double tolerance, HODLR_Matrix* A);
	~HODLR_Tree();

	//  Methods for HODLR solver
	void assembleTree();
	void factorize();
	double determinant();
	void matmatProduct(Eigen::MatrixXd x, Eigen::MatrixXd& b);
	Eigen::MatrixXd solve(Eigen::MatrixXd b);
};

#endif /*__HODLR_Tree__*/