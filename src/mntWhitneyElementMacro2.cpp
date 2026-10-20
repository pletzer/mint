#include <cmath>
#include <deque>
#include <set>
#include <string>
#include <mntWhitneyElementMacro2.h>
#include <mntGlobal.h>
#include <mntLogger.h>

// cell positions (p, q) in a block, in the order used by blockCells
static const int POS_P[] = {0, 1, 1, 0};
static const int POS_Q[] = {0, 0, 1, 1};

// local node (0-3) of the cell at each position that is the centre of the block
static const int CENTRE_LOCAL_NODE[] = {2, 3, 0, 1};

static int posIndex(int p, int q) {
    return (q == 0)? p: 3 - p;
}

/**
 * Find the 4 cells of the block centred on node centre
 * @param faceNodes face-node connectivity
 * @param node2cells cells around each node
 * @param centre candidate block centre
 * @param cells cells at positions (0,0), (1,0), (1,1), (0,1) (output)
 * @return true if the cells form a 2x2 block with a shared parametric frame
 */
static bool findBlock(const std::vector<std::size_t>& faceNodes,
                      const std::vector< std::vector<vtkIdType> >& node2cells,
                      std::size_t centre, vtkIdType cells[]) {

    const std::vector<vtkIdType>& around = node2cells[centre];
    if (around.size() != 4) {
        return false;
    }

    for (int k = 0; k < 4; ++k) cells[k] = -1;
    for (vtkIdType c : around) {
        int loc = -1;
        for (int n = 0; n < 4; ++n) {
            if (faceNodes[4*c + n] == centre) loc = n;
        }
        // the local node of the centre fixes the position of the cell in the block
        int pos = -1;
        for (int k = 0; k < 4; ++k) {
            if (CENTRE_LOCAL_NODE[k] == loc) pos = k;
        }
        if (pos < 0 || cells[pos] >= 0) {
            return false;
        }
        cells[pos] = c;
    }

    // the cells share their parametric frame if the mid nodes of the block's
    // outer edges are seen consistently from both adjacent cells
    const std::size_t* n00 = &faceNodes[4*cells[0]];
    const std::size_t* n10 = &faceNodes[4*cells[1]];
    const std::size_t* n11 = &faceNodes[4*cells[2]];
    const std::size_t* n01 = &faceNodes[4*cells[3]];
    return n10[0] == n00[1] && n01[0] == n00[3] && n11[1] == n10[2] && n11[3] == n01[2];
}

int
WhitneyElementMacro2::build(Grid_t* grid) {

    std::string msg;
    const std::vector<std::size_t>& faceNodes = grid->faceNodeConnectivity;
    if (faceNodes.size() == 0) {
        msg = "the 2x2 block element needs the face-node connectivity, load the grid from UGRID";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    const std::size_t numCells = faceNodes.size() / MNT_NUM_VERTS_PER_QUAD;
    std::size_t numNodes = 0;
    for (std::size_t n : faceNodes) numNodes = std::max(numNodes, n + 1);

    std::vector< std::vector<vtkIdType> > node2cells(numNodes);
    for (std::size_t c = 0; c < numCells; ++c) {
        for (int k = 0; k < MNT_NUM_VERTS_PER_QUAD; ++k) {
            node2cells[faceNodes[4*c + k]].push_back((vtkIdType) c);
        }
    }

    this->cellBlock.assign(numCells, -1);
    this->cellPos.assign(numCells, -1);
    this->blockCells.clear();

    // Grow the tiling block by block from a seed block. The block across each
    // macro edge is centred on the node shared by the two cells across that
    // edge, other than the macro edge's mid node. A wrong seed (the centre of
    // a 2x2 group that is not a block of the tiling) produces overlapping
    // blocks, in which case the next seed is tried.
    while (true) {

        // first cell not in a block
        std::size_t firstCell = numCells;
        for (std::size_t c = 0; c < numCells; ++c) {
            if (this->cellBlock[c] < 0) {
                firstCell = c;
                break;
            }
        }
        if (firstCell == numCells) {
            break;
        }

        bool success = false;
        for (int k = 0; k < MNT_NUM_VERTS_PER_QUAD && !success; ++k) {

            std::vector<long long> cellBlockTrial = this->cellBlock;
            std::vector<int> cellPosTrial = this->cellPos;
            std::vector<vtkIdType> blockCellsTrial = this->blockCells;

            std::deque<std::size_t> queue;
            std::set<std::size_t> seen;
            std::size_t seed = faceNodes[4*firstCell + k];
            queue.push_back(seed);
            seen.insert(seed);
            success = true;

            while (!queue.empty()) {

                std::size_t centre = queue.front();
                queue.pop_front();

                vtkIdType cells[4];
                if (!findBlock(faceNodes, node2cells, centre, cells)) {
                    success = false;
                    break;
                }
                bool overlaps = false;
                for (int i = 0; i < 4; ++i) {
                    if (cellBlockTrial[cells[i]] >= 0) overlaps = true;
                }
                if (overlaps) {
                    success = false;
                    break;
                }
                long long blockId = (long long) (blockCellsTrial.size() / 4);
                for (int i = 0; i < 4; ++i) {
                    cellBlockTrial[cells[i]] = blockId;
                    cellPosTrial[cells[i]] = i;
                    blockCellsTrial.push_back(cells[i]);
                }

                // mid nodes of the macro edges: south, east, north, west
                const std::size_t mids[] = {faceNodes[4*cells[0] + 1], faceNodes[4*cells[1] + 2],
                                            faceNodes[4*cells[2] + 3], faceNodes[4*cells[0] + 3]};
                for (std::size_t mid : mids) {
                    std::vector<vtkIdType> outside;
                    for (vtkIdType c : node2cells[mid]) {
                        if (c != cells[0] && c != cells[1] && c != cells[2] && c != cells[3]) {
                            outside.push_back(c);
                        }
                    }
                    if (outside.size() != 2) {
                        // domain boundary
                        continue;
                    }
                    for (int i = 0; i < 4; ++i) {
                        std::size_t n = faceNodes[4*outside[0] + i];
                        if (n == mid) continue;
                        for (int j = 0; j < 4; ++j) {
                            if (faceNodes[4*outside[1] + j] == n && seen.count(n) == 0) {
                                seen.insert(n);
                                queue.push_back(n);
                            }
                        }
                    }
                }
            }

            if (success) {
                this->cellBlock = cellBlockTrial;
                this->cellPos = cellPosTrial;
                this->blockCells = blockCellsTrial;
            }
        }

        if (!success) {
            std::size_t numUntiled = 0;
            for (long long b : this->cellBlock) {
                if (b < 0) numUntiled++;
            }
            msg = "could not tile the grid into 2x2 blocks of cells, cell " + std::to_string(firstCell) +
                  " is not in any block (" + std::to_string(numUntiled) + " cells not tiled)";
            mntlog::error(__FILE__, __func__, __LINE__, msg);
            return (int) numUntiled;
        }
    }

    return 0;
}

void
WhitneyElementMacro2::getBlockCells(std::size_t blockId, vtkIdType cellIds[]) const {
    for (int i = 0; i < 4; ++i) {
        cellIds[i] = this->blockCells[4*blockId + i];
    }
}

int
WhitneyElementMacro2::getCellBlock(vtkIdType cellId, std::size_t* blockId, int* p, int* q) const {
    if (cellId < 0 || (std::size_t) cellId >= this->cellBlock.size() || this->cellBlock[cellId] < 0) {
        return 1;
    }
    *blockId = (std::size_t) this->cellBlock[cellId];
    *p = POS_P[this->cellPos[cellId]];
    *q = POS_Q[this->cellPos[cellId]];
    return 0;
}

void
WhitneyElementMacro2::getDofEdge(int dofIndex, int* p, int* q, int* edgeIndex) {
    if (dofIndex < 6) {
        // X-directed edges at Y = 0, 1/2, 1
        int i = dofIndex % 2;
        int j = dofIndex / 2;
        *p = i;
        *q = (j == 2)? 1: 0;
        *edgeIndex = (j == 0)? 0: 2;
    }
    else {
        // Y-directed edges at X = 0, 1/2, 1
        int i = (dofIndex - 6) / 2;
        int j = (dofIndex - 6) % 2;
        *p = (i == 2)? 1: 0;
        *q = j;
        *edgeIndex = (i == 0)? 3: 1;
    }
}

void
WhitneyElementMacro2::evalMacroBasis(double x, double y, double a[], double b[]) {

    // duals of the half-interval integrals
    const double phiX[] = {3.0 - 4.0*x, 4.0*x - 1.0};
    const double phiY[] = {3.0 - 4.0*y, 4.0*y - 1.0};

    // quadratic Lagrange polynomials at 0, 1/2, 1
    const double lagX[] = {(1.0 - x)*(1.0 - 2.0*x), 4.0*x*(1.0 - x), x*(2.0*x - 1.0)};
    const double lagY[] = {(1.0 - y)*(1.0 - 2.0*y), 4.0*y*(1.0 - y), y*(2.0*y - 1.0)};

    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 2; ++i) {
            a[2*j + i] = phiX[i]*lagY[j];
            b[2*j + i] = 0.0;
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            a[6 + 2*i + j] = 0.0;
            b[6 + 2*i + j] = lagX[i]*phiY[j];
        }
    }
}

void
WhitneyElementMacro2::evalBasis(vtkIdType cellId, double xi, double eta, double a[], double b[]) const {

    std::size_t blockId;
    int p, q;
    if (this->getCellBlock(cellId, &blockId, &p, &q) != 0) {
        for (int k = 0; k < 12; ++k) {
            a[k] = 0.0;
            b[k] = 0.0;
        }
        return;
    }

    evalMacroBasis(0.5*(p + xi), 0.5*(q + eta), a, b);

    // dX = dxi/2 and dY = deta/2
    for (int k = 0; k < 12; ++k) {
        a[k] *= 0.5;
        b[k] *= 0.5;
    }
}

void
WhitneyElementMacro2::integrateBasis(vtkIdType cellId, const Vec3& xia, const Vec3& xib, double w[]) const {

    for (int k = 0; k < 12; ++k) w[k] = 0.0;

    std::size_t blockId;
    int p, q;
    if (this->getCellBlock(cellId, &blockId, &p, &q) != 0) {
        return;
    }

    // the integrands are at most cubic along a straight segment, the 2-point
    // Gauss-Legendre rule is exact
    const double ts[] = {0.5 - 0.5/std::sqrt(3.0), 0.5 + 0.5/std::sqrt(3.0)};

    const double xa = 0.5*(p + xia[0]);
    const double ya = 0.5*(q + xia[1]);
    const double dx = 0.5*(xib[0] - xia[0]);
    const double dy = 0.5*(xib[1] - xia[1]);

    double a[12], b[12];
    for (double t : ts) {
        evalMacroBasis(xa + t*dx, ya + t*dy, a, b);
        for (int k = 0; k < 12; ++k) {
            w[k] += 0.5*(a[k]*dx + b[k]*dy);
        }
    }
}

int
WhitneyElementMacro2::getCounterclockwiseSign(int dofIndex) const {
    int p, q, edgeIndex;
    getDofEdge(dofIndex, &p, &q, &edgeIndex);
    // the north and west edges of a cell point clockwise
    return (edgeIndex < 2)? 1: -1;
}

int
WhitneyElementMacro2::getDofDataIndex(Grid_t* grid, vtkIdType cellId, int dofIndex, int placement,
                                      std::size_t* index, int* sign) const {

    std::size_t blockId;
    int p, q;
    if (this->getCellBlock(cellId, &blockId, &p, &q) != 0) {
        return 1;
    }

    // the cell and edge holding this degree of freedom
    int pe, qe, edgeIndex;
    getDofEdge(dofIndex, &pe, &qe, &edgeIndex);
    vtkIdType edgeCellId = this->blockCells[4*blockId + posIndex(pe, qe)];

    if (placement == MNT_CELL_BY_CELL_DATA) {
        *index = edgeCellId*MNT_NUM_EDGES_PER_QUAD + edgeIndex;
        *sign = 1;
        return 0;
    }
    return mnt_grid_getEdgeId(&grid, edgeCellId, edgeIndex, index, sign);
}
