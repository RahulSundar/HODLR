#ifndef __HODLR_Tree__
#define __HODLR_Tree__

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <vector>

#include "HODLR_Matrix.hpp"
#include "HODLR_Node.hpp"

class HODLR_Tree {
private:
	int N;
	int nLevels;
	double tolerance;
	std::vector<int> nodesInLevel;
	HODLR_Matrix* A;
public:

	HODLR_Tree(int nLevels, double tolerance, HODLR_Matrix* A);
	~HODLR_Tree();
	virtual double get_Matrix_Element(int j, int k) {
		return j*j+k*k;
	}
	std::vector< std::vector<HODLR_Node*> > tree;
	void createTree();
	void createRoot();
	void createChildren(int j, int k);
	void assembleTree();
	void matmat_Product(Eigen::MatrixXd x, Eigen::MatrixXd& b);
	void factorize();
	void factorize_Leaf(int k);
	void factorize_Non_Leaf(int j, int k);
	Eigen::MatrixXd solve_Leaf(int k, Eigen::MatrixXd b);
	Eigen::MatrixXd solve_Non_Leaf(int j, int k, Eigen::MatrixXd b);
	Eigen::MatrixXd solve(Eigen::MatrixXd b);
};

HODLR_Tree::HODLR_Tree(int nLevels, double tolerance, HODLR_Matrix* A) {
	// // std::cout << "\nStart HODLR_Tree\n";
	this->nLevels	=	nLevels;
	this->tolerance	=	tolerance;
	this->A			=	A;
	this->N			=	A->N;
	nodesInLevel.push_back(1);
	for (int j=1; j<=nLevels; ++j) {
		nodesInLevel.push_back(2*nodesInLevel.back());
	}
	// // std::cout << "\nDone HODLR_Tree\n";
	createTree();
}

void HODLR_Tree::createRoot() {
	// // std::cout << "\nStart createRoot\n";
	HODLR_Node* root	=	new HODLR_Node(0, 0, 0, 0, N, tolerance);
	std::vector<HODLR_Node*> level;
	level.push_back(root);
	tree.push_back(level);
	// // std::cout << "\nDone createRoot\n";
}

void HODLR_Tree::createChildren(int j, int k) {
	// // std::cout << "\nStart createChildren\n";
	//	Adding left child
	HODLR_Node* left	=	new HODLR_Node(2*j, k+1, 0, tree[j][k]->cStart[0], tree[j][k]->cSize[0], tolerance);
	tree[j+1].push_back(left);

	//	Adding right child
	HODLR_Node* right	=	new HODLR_Node(2*j+1, k+1, 1, tree[j][k]->cStart[1], tree[j][k]->cSize[1], tolerance);
	tree[j+1].push_back(right);
	// // std::cout << "\nDone createChildren\n";
}

void HODLR_Tree::createTree() {
	// // std::cout << "\nStart createTree\n";
	createRoot();
	for (int j=0; j<nLevels; ++j) {
		std::vector<HODLR_Node*> level;
		tree.push_back(level);
		for (int k=0; k<nodesInLevel[j]; ++k) {
			createChildren(j,k);
		}
	}
	// // std::cout << "\nDone createTree\n";
}

void HODLR_Tree::assembleTree() {
	// // std::cout << "\nStart assembleTree\n";
	for (int j=0; j<nLevels; ++j) {
		#pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			tree[j][k]->assemble_Non_Leaf_Node(A);
		}
	}
	#pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		tree[nLevels][k]->assemble_Leaf_Node(A);
	}
	// // std::cout << "\nDone assembleTree\n";
}

void HODLR_Tree::matmat_Product(Eigen::MatrixXd x, Eigen::MatrixXd& b) {
	// // std::cout << "\nStart matmat_Product\n";
	int r	=	x.cols();
	b		=	Eigen::MatrixXd::Zero(N,r);
	for (int j=0; j<nLevels; ++j) {
		#pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			tree[j][k]->matmat_Product_Non_Leaf(x, b);
		}
	}
	#pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		tree[nLevels][k]->matmat_Product_Leaf(x, b);
	}
	// // std::cout << "\nDone matmat_Product\n";
}

//	Factorization begins from here

void HODLR_Tree::factorize_Leaf(int k) {
	// std::cout << "\nStart factorize Leaf: " << k << "\n";
	tree[nLevels][k]->Kfactor.compute(tree[nLevels][k]->K);
	int parent	=	k;
	int child	=	k;
	int size	=	tree[nLevels][k]->nSize;
	int tstart, r;
	#pragma omp parallel for
	for (int l=nLevels-1; l>=0; --l) {
		child	=	parent%2;
		parent	=	parent/2;
		tstart	=	tree[nLevels][k]->nStart-tree[l][parent]->cStart[child];
		r		=	tree[l][parent]->rank[child];
		// std::cout << size << "\t" << child << "\t" << parent << "\t" << tstart << "\t" << r << "\n";
		// std::cout << tree[l][parent]->Ufactor[child].rows() << "\t" << tree[l][parent]->Ufactor[child].cols() << "\n";
		// std::cout << l << "\n";
		tree[l][parent]->Ufactor[child].block(tstart,0,size,r)	=	solve_Leaf(k, tree[l][parent]->Ufactor[child].block(tstart,0,size,r));
	}
	// std::cout << "\nDone factorize Leaf: " << k << "\n";
}

Eigen::MatrixXd HODLR_Tree::solve_Leaf(int k, Eigen::MatrixXd b) {
	// std::cout << "\nStart solve Leaf: " << k << "\n";
	Eigen::MatrixXd x	=	tree[nLevels][k]->Kfactor.solve(b);
	// std::cout << "\nDone solve Leaf: " << k << "\n";
	return x;
}

void HODLR_Tree::factorize_Non_Leaf(int j, int k) {
	// std::cout << "\nStart factorize; Level: " << j << "Node: " << k << "\n";
	int r0	=	tree[j][k]->rank[0];
	int r1	=	tree[j][k]->rank[1];
	tree[j][k]->K	=	Eigen::MatrixXd::Identity(r0+r1, r0+r1);
	// std::cout << tree[j][k]->K.rows() << "\t" << tree[j][k]->K.cols() << "\n";
	tree[j][k]->K.block(0, r0, r0, r1)	=	tree[j][k]->Vfactor[1].transpose()*tree[j][k]->Ufactor[1];
	tree[j][k]->K.block(r0, 0, r1, r0)	=	tree[j][k]->Vfactor[0].transpose()*tree[j][k]->Ufactor[0];
	tree[j][k]->Kfactor.compute(tree[j][k]->K);
	int parent	=	k;
	int child	=	k;
	int size	=	tree[j][k]->nSize;
	int tstart, r;
	#pragma omp parallel for
	for (int l=j-1; l>=0; --l) {
		child	=	parent%2;
		parent	=	parent/2;
		tstart	=	tree[j][k]->nStart-tree[l][parent]->cStart[child];
		r		=	tree[l][parent]->rank[child];
		// std::cout << "\n" << size << "\n";
		tree[l][parent]->Ufactor[child].block(tstart,0,size,r)	=	solve_Non_Leaf(j, k, tree[l][parent]->Ufactor[child].block(tstart,0,size,r));
	}
	// std::cout << "\nDone factorize; Level: " << j << "Node: " << k << "\n";
}

Eigen::MatrixXd HODLR_Tree::solve_Non_Leaf(int j, int k, Eigen::MatrixXd b) {
	// std::cout << "\nStart Solve non leaf; Level: " << j << "Node: " << k << "\n";
	int r0	=	tree[j][k]->rank[0];
	int r1	=	tree[j][k]->rank[1];
	int n0	=	tree[j][k]->cSize[0];
	int n1	=	tree[j][k]->cSize[1];
	int r	=	b.cols();
	Eigen::MatrixXd temp(r0+r1, r);
	// std::cout << "\nHere\n";
	// std::cout << tree[j][k]->Vfactor[1].rows() << "\t" << tree[j][k]->Vfactor[1].cols() << "\n";
	// std::cout << tree[j][k]->Vfactor[0].rows() << "\t" << tree[j][k]->Vfactor[0].cols() << "\n";
	// std::cout << b.rows() << "\t" << b.cols() << "\n";
	// std::cout << n0 << "\t" << n1 << "\t" << r0 << "\t" << r1 << "\n";
	temp << tree[j][k]->Vfactor[1].transpose()*b.block(n0,0,n1,r),
			tree[j][k]->Vfactor[0].transpose()*b.block(0,0,n0,r);
	// std::cout << "\nHere\n";
	// std::cout << "\n" << temp.rows() << "\t" << temp.cols() << "\n";
	// std::cout << "\n" << tree[j][k]->K.rows() << "\t" << tree[j][k]->K.cols() << "\n";
	temp	=	tree[j][k]->Kfactor.solve(temp);
	// std::cout << "\nHere\n";
	Eigen::MatrixXd y(n0+n1, r);
	y << tree[j][k]->Ufactor[0]*temp.block(0,0,r0,r), tree[j][k]->Ufactor[1]*temp.block(r0,0,r1,r);
	// std::cout << "\nEnd Solve non leaf; Level: " << j << "Node: " << k << "\n";
	return (b-y);
}

void HODLR_Tree::factorize() {
	// std::cout << "\nStart factorize...\n";
	//	Set things for factorization
	// std::cout << "Number of levels: " << nLevels << "\n";
	for (int j=0; j<=nLevels; ++j) {
		#pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			for (int l=0; l<2; ++l) {
				tree[j][k]->Ufactor[l]	=	tree[j][k]->U[l];
				tree[j][k]->Vfactor[l]	=	tree[j][k]->V[l];
				// std::cout << "Level " << j << "; Node " << k << "; Matrix Size" << tree[j][k]->Ufactor[l].rows() << "\t" << tree[j][k]->Ufactor[l].cols() << "\n";
			}
		}
	}

	//	Factorize the leaf nodes
	#pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		factorize_Leaf(k);
	}
	//	Factorize the non-leaf nodes
	for (int j=nLevels-1; j>=0; --j) {
		#pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			factorize_Non_Leaf(j, k);
		}
	}
	// std::cout << "\nEnd factorize...\n";
}

Eigen::MatrixXd HODLR_Tree::solve(Eigen::MatrixXd b) {
	// std::cout << "\nStart solve...\n";
	int start, size;
	Eigen::MatrixXd x	=	Eigen::MatrixXd::Zero(b.rows(),b.cols());
	int r	=	b.cols();
	// #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		x.block(start, 0, size, r)	=	solve_Leaf(k, b.block(start, 0, size, r));
	}
	b=x;
	for (int j=nLevels-1; j>=0; --j) {
		// #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			x.block(start, 0, size, r)	=	solve_Non_Leaf(j, k, b.block(start, 0, size, r));
		}
		b=x;
	}
	// std::cout << "\nEnd solve...\n";
	return x;
}

//Symmetric Factorization

void HODLR_Tree::assembleSymmetricTree() {
	// // std::cout << "\nStart assembleTree\n";
	for (int j=0; j<nLevels; ++j) {
		#pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; k=k+2) {
			tree[j][k]->assemble_Symmetric_Non_Leaf_Node(A);
           /* //Try using pointers instead of storing the same array twice
			tree[j][k+1]->Q[0] = tree[j][k]->Q[1];
			tree[j][k+1]->Q[1] = tree[j][k]->Q[0];*/
		}
	}
	#pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {

		tree[nLevels][k]->assemble_Leaf_Node(A);
	}
	// // std::cout << "\nDone assembleTree\n";
}



void symmetric_factorize(){

    /*for (int j=0; j<=nLevels; ++j) {
		#pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			for (int l=0; l<2; ++l) {
				tree[j][k]->Qfactor[l]	=	tree[j][k]->Q[l];
//				tree[j][k]->Vfactor[l]	=	tree[j][k]->V[l];
				// std::cout << "Level " << j << "; Node " << k << "; Matrix Size" << tree[j][k]->Ufactor[l].rows() << "\t" << tree[j][k]->Ufactor[l].cols() << "\n";
			}
		}
	}*/
//	Factorize the leaf nodes
	#pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		factorize_Symmetric_Leaf(k);
	}

    qr(nLevels-1);

//	Factorize the non-leaf nodes
	for (int j=nLevels-1; j>=0; --j) {
		#pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			factorize_Symmetric_Non_Leaf(j, k);
		}
		qr(j);
	}
	// std::cout << "\nEnd factorize...\n";
}


void HODLR_Tree::qr_Level(int level){

    for(int k=0; k<nodesInLevel[level]; k++){
            qr(level, k);
        }
    }

}

void HODLR_Tree::qr(int j, int k){

    //tree[j][k]->Q[0];
    //tree[j][k]->Q[1];
    //tree[j][k]->R;
    int min[2];
    min[0] = std::min(tree[j][k]->Q[0].rows(), tree[j][k]->Q[0].cols());
    min[1] = std::min(tree[j][k]->Q[1].rows(), tree[j][k]->Q[1].cols());

    HouseholderQR<MatrixXd> qr(Q[0]);
    tree[j][k]->Q[0] = qr.householderQ()*(Eigen::MatrixXd::Identity((Q[0].rows(), m[0]));
    if(R!=null){
        tree[j][k]->R = qr.matrixQR().block(0,0,m[0],Q[0].cols()).triangularView<Eigen::Upper>()*tree[j][k]->R;
    } else {
        tree[j][k]->R = qr.matrixQR().block(0,0,m[0],Q[0].cols()).triangularView<Eigen::Upper>();
    }
    }

    qr(Q[1]);
    tree[j][k]->Q[1] = qr.householderQ()*(Eigen::MatrixXd::Identity((Q[1].rows(), m[1]));
    tree[j][k]->R *= qr.matrixQR().block(0,0,m[1],Q[1].cols()).triangularView<Eigen::Upper>().transpose();

}

void HODLR_Tree::factorize_Symmetric_Leaf(int k) {
	// std::cout << "\nStart factorize Leaf: " << k << "\n";

    //Cholesky of W:
	tree[nLevels][k]->W = tree[nLevels][k]->K.llt().matrixL();
	int parent	=	k;
	int child	=	k;
	int size	=	tree[nLevels][k]->nSize;
	int tstart, r;
	#pragma omp parallel for
	for (int l=nLevels-1; l>=0; --l) {
		child	=	parent%2;
		parent	=	parent/2;
		tstart	=	tree[nLevels][k]->nStart-tree[l][parent]->cStart[child];
		r		=	tree[l][parent]->sym_rank;
		// std::cout << size << "\t" << child << "\t" << parent << "\t" << tstart << "\t" << r << "\n";
		// std::cout << tree[l][parent]->Ufactor[child].rows() << "\t" << tree[l][parent]->Ufactor[child].cols() << "\n";
		// std::cout << l << "\n";
		tree[l][parent]->Q[child].block(tstart,0,size,r)	=	solve_Symmetric_Leaf(k, tree[l][parent]->Q[child].block(tstart,0,size,r));
	}
	// std::cout << "\nDone factorize Leaf: " << k << "\n";
}

Eigen::MatrixXd HODLR_Tree::solve_Symmetric_Leaf(int k, Eigen::MatrixXd b) {
	// std::cout << "\nStart solve Leaf: " << k << "\n";
	Eigen::MatrixXd x	=	tree[nLevels][k]->W.solve(b);

	// std::cout << "\nDone solve Leaf: " << k << "\n";
	return x;
}

void HODLR_Tree::factorize_Non_Leaf(int j, int k) {
	// std::cout << "\nStart factorize; Level: " << j << "Node: " << k << "\n";
	tree[j][k] X = symmetric_Cholesky(j,k,tree[j][k]->sym_rank);
	int parent	=	k;
	int child	=	k;
	int size	=	tree[j][k]->nSize;
	int tstart, r;
	#pragma omp parallel for
	for (int l=j-1; l>=0; --l) {
		child	=	parent%2;
		parent	=	parent/2;
		tstart	=	tree[j][k]->nStart-tree[l][parent]->cStart[child];
		r		=	tree[l][parent]->sym_rank;
		// std::cout << "\n" << size << "\n";
		tree[l][parent]->Qfactor[child].block(tstart,0,size,r)	=	solve_Symmetric_Non_Leaf(j, k, tree[l][parent]->Qfactor[child].block(tstart,0,size,r));
	}
	// std::cout << "\nDone factorize; Level: " << j << "Node: " << k << "\n";
}


Eigen::MatrixXd HODLR_Tree::symmetric_Cholesky(int j, int k,int r) {

    // Finding X using symmetric factorization
    // what should be the rank of R1*R2??
    tree[j][k]->K	=	Eigen::MatrixXd::Identity(2*r, 2*r);
    tree[j][k]->K.block(0, r, r, r)	=	tree[j][k]->R;
    tree[j][k]->K.block(r, 0, r, r)	=	tree[j][k]->R.transpose();
    Eigen::MatrixXd x = tree[j][k]->K.llt().matrixL();
    x -= Eigen::MatrixXd::Identity(2*r, 2*r);
	return x;
}

Eigen::MatrixXd HODLR_Tree::solve_Symmetric_Non_Leaf(int j, int k, Eigen::MatrixXd b) {
	// std::cout << "\nStart Solve non leaf; Level: " << j << "Node: " << k << "\n";
	/*int r	=	tree[j][k]->sym_rank;
	int n0	=	tree[j][k]->cSize[0];
	int n1	=	tree[j][k]->cSize[1];
	Eigen::MatrixXd temp(2*r, r);
	temp << tree[j][k]->Vfactor[1].transpose()*b.block(n0,0,n1,r),
			tree[j][k]->Vfactor[0].transpose()*b.block(0,0,n0,r);
	// std::cout << "\nHere\n";
	// std::cout << "\n" << temp.rows() << "\t" << temp.cols() << "\n";
	// std::cout << "\n" << tree[j][k]->K.rows() << "\t" << tree[j][k]->K.cols() << "\n";
	temp	=	tree[j][k]->Kfactor.solve(temp);
	// std::cout << "\nHere\n";
	Eigen::MatrixXd y(n0+n1, r);
	y << tree[j][k]->Ufactor[0]*temp.block(0,0,r0,r), tree[j][k]->Ufactor[1]*temp.block(r0,0,r1,r);
	// std::cout << "\nEnd Solve non leaf; Level: " << j << "Node: " << k << "\n";
	return (b-y);*/
}


#endif /*__HODLR_Tree__*/
