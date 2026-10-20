#define _USE_MATH_DEFINES // M_PI for Visual Studio
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

// quadratic Lagrange polynomials at 0, 1/2, 1 and their derivatives
static void lagrange2(double t, double l[], double dl[]) {
    l[0] = (1.0 - t)*(1.0 - 2.0*t);  dl[0] = 4.0*t - 3.0;
    l[1] = 4.0*t*(1.0 - t);          dl[1] = 4.0 - 8.0*t;
    l[2] = t*(2.0*t - 1.0);          dl[2] = 4.0*t - 1.0;
}

int
WhitneyElementMacro2::build(Grid_t* grid, bool biquadratic) {

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

    // block geometry
    const std::size_t numBlocks = this->getNumBlocks();
    this->blockBiquadratic.assign(numBlocks, 0);
    this->blockNodes.assign(numBlocks*9*2, 0.0);
    this->cellVerts.assign(numCells*4*2, 0.0);
    const double* verts = grid->verts;
    if (!verts) {
        // no vertex coordinates, the macro coordinates of the cells are used
        return 0;
    }
    const bool periodic = grid->fixLonAcrossDateline;
    const double period = grid->periodX;
    const double poleLat = grid->degrees? 90.0: 0.5*M_PI;
    // corner offsets of the local vertices of a cell
    const int cornerA[] = {0, 1, 1, 0};
    const int cornerB[] = {0, 0, 1, 1};
    for (std::size_t blockId = 0; blockId < numBlocks; ++blockId) {

        const double lonRef = verts[3*MNT_NUM_VERTS_PER_QUAD*this->blockCells[4*blockId]];
        bool ok = true;
        bool nodeSet[9] = {false, false, false, false, false, false, false, false, false};

        for (int pos = 0; pos < 4; ++pos) {
            vtkIdType c = this->blockCells[4*blockId + pos];
            for (int k = 0; k < MNT_NUM_VERTS_PER_QUAD; ++k) {
                double lon = verts[3*(MNT_NUM_VERTS_PER_QUAD*c + k) + LON_INDEX];
                double lat = verts[3*(MNT_NUM_VERTS_PER_QUAD*c + k) + LAT_INDEX];
                if (periodic) {
                    // continuous longitudes within the block
                    while (lon - lonRef > 0.5*period) lon -= period;
                    while (lon - lonRef < -0.5*period) lon += period;
                    if (std::abs(lat) > poleLat*(1.0 - 1.e-10)) ok = false;
                }
                this->cellVerts[8*c + 2*k + 0] = lon;
                this->cellVerts[8*c + 2*k + 1] = lat;

                int a = POS_P[pos] + cornerA[k];
                int b = POS_Q[pos] + cornerB[k];
                double* node = &this->blockNodes[(9*blockId + 3*b + a)*2];
                if (nodeSet[3*b + a]) {
                    // the cells of the block must agree on their shared nodes
                    double scale = std::max(1.0, std::abs(node[0]) + std::abs(node[1]));
                    if (std::abs(node[0] - lon) + std::abs(node[1] - lat) > 1.e-10*scale) ok = false;
                }
                node[0] = lon;
                node[1] = lat;
                nodeSet[3*b + a] = true;
            }
        }
        this->blockBiquadratic[blockId] = (biquadratic && ok)? 1: 0;
    }

    return 0;
}

std::size_t
WhitneyElementMacro2::getNumBiquadraticBlocks() const {
    std::size_t n = 0;
    for (char b : this->blockBiquadratic) n += b;
    return n;
}

void
WhitneyElementMacro2::evalMacroMap(std::size_t blockId, double x, double y, double g[], double jac[][2]) const {
    double lx[3], dlx[3], ly[3], dly[3];
    lagrange2(x, lx, dlx);
    lagrange2(y, ly, dly);
    g[0] = g[1] = 0.0;
    jac[0][0] = jac[0][1] = jac[1][0] = jac[1][1] = 0.0;
    for (int b = 0; b < 3; ++b) {
        for (int a = 0; a < 3; ++a) {
            const double* node = &this->blockNodes[(9*blockId + 3*b + a)*2];
            for (int d = 0; d < 2; ++d) {
                g[d] += node[d]*lx[a]*ly[b];
                jac[d][0] += node[d]*dlx[a]*ly[b];
                jac[d][1] += node[d]*lx[a]*dly[b];
            }
        }
    }
}

int
WhitneyElementMacro2::getMacroCoords(vtkIdType cellId, double xi, double eta, double xy[], double dxy[][2]) const {

    std::size_t blockId;
    int p, q;
    if (this->getCellBlock(cellId, &blockId, &p, &q) != 0) {
        return 1;
    }

    // macro coordinates of the piecewise bilinear cell maps, also the initial guess
    xy[0] = 0.5*(p + xi);
    xy[1] = 0.5*(q + eta);
    dxy[0][0] = 0.5; dxy[0][1] = 0.0;
    dxy[1][0] = 0.0; dxy[1][1] = 0.5;
    if (!this->blockBiquadratic[blockId]) {
        return 0;
    }

    // position of the point and its derivatives with the bilinear cell map
    const double* v = &this->cellVerts[8*cellId];
    double pt[2], dpt[2][2];
    for (int d = 0; d < 2; ++d) {
        pt[d] = v[d]*(1.0 - xi)*(1.0 - eta) + v[2 + d]*xi*(1.0 - eta) + v[4 + d]*xi*eta + v[6 + d]*(1.0 - xi)*eta;
        dpt[d][0] = (v[2 + d] - v[d])*(1.0 - eta) + (v[4 + d] - v[6 + d])*eta;
        dpt[d][1] = (v[6 + d] - v[d])*(1.0 - xi) + (v[4 + d] - v[2 + d])*xi;
    }

    // invert the biquadratic map. Newton converges quadratically from the
    // bilinear cell coordinates; the last steps are limited by round-off in the
    // residual (coordinates / cell size), so stop once the step no longer decreases
    double x = xy[0], y = xy[1];
    double g[2], jac[2][2], det = 0.0;
    double step = 1.0, prevStep = 2.0;
    for (int iter = 0; iter < 20 && step > 1.e-15 && step < prevStep; ++iter) {
        this->evalMacroMap(blockId, x, y, g, jac);
        det = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
        if (det == 0.0) break;
        double r0 = pt[0] - g[0], r1 = pt[1] - g[1];
        double dx = ( jac[1][1]*r0 - jac[0][1]*r1)/det;
        double dy = (-jac[1][0]*r0 + jac[0][0]*r1)/det;
        x += dx;
        y += dy;
        prevStep = step;
        step = std::abs(dx) + std::abs(dy);
    }
    if (det == 0.0 || step > 1.e-8) {
        // keep the macro coordinates of the bilinear cell maps
        return 0;
    }

    // d(X, Y)/d(xi, eta) = (dG/d(X, Y))^-1 d(lon, lat)/d(xi, eta)
    this->evalMacroMap(blockId, x, y, g, jac);
    det = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
    const double inv[2][2] = {{ jac[1][1]/det, -jac[0][1]/det},
                              {-jac[1][0]/det,  jac[0][0]/det}};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            dxy[i][j] = inv[i][0]*dpt[0][j] + inv[i][1]*dpt[1][j];
        }
    }
    xy[0] = x;
    xy[1] = y;
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

    double xy[2], dxy[2][2];
    if (this->getMacroCoords(cellId, xi, eta, xy, dxy) != 0) {
        for (int k = 0; k < 12; ++k) {
            a[k] = 0.0;
            b[k] = 0.0;
        }
        return;
    }

    double am[12], bm[12];
    evalMacroBasis(xy[0], xy[1], am, bm);

    // a dX + b dY in terms of dxi and deta
    for (int k = 0; k < 12; ++k) {
        a[k] = am[k]*dxy[0][0] + bm[k]*dxy[1][0];
        b[k] = am[k]*dxy[0][1] + bm[k]*dxy[1][1];
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

    if (this->blockBiquadratic[blockId]) {
        // the segment is a curve in macro coordinates, 5-point Gauss-Legendre rule
        const double s = std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0))/3.0;
        const double t = std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0))/3.0;
        const double gaussX[] = {-t, -s, 0.0, s, t};
        const double gaussW[] = {(322.0 - 13.0*std::sqrt(70.0))/900.0, (322.0 + 13.0*std::sqrt(70.0))/900.0,
                                 128.0/225.0,
                                 (322.0 + 13.0*std::sqrt(70.0))/900.0, (322.0 - 13.0*std::sqrt(70.0))/900.0};
        const double dxi = xib[0] - xia[0];
        const double deta = xib[1] - xia[1];
        double a[12], b[12];
        for (int ig = 0; ig < 5; ++ig) {
            double tg = 0.5*(1.0 + gaussX[ig]);
            this->evalBasis(cellId, xia[0] + tg*dxi, xia[1] + tg*deta, a, b);
            for (int k = 0; k < 12; ++k) {
                w[k] += 0.5*gaussW[ig]*(a[k]*dxi + b[k]*deta);
            }
        }
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
