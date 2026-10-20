#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGlobal.h>
#include <mntGrid.h>
#include <mntPolylineIntegral.h>
#include <mntWhitneyElementMacro2.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <cstdio>
#include <map>
#include <array>
#include <random>
#include <vector>

/**
 * Line integrals of the gradient of a smooth potential on equiangular
 * cubed-sphere grids of increasing resolution, along random paths that cross
 * the panel edges. The edge data are the exact edge integrals, i.e. differences
 * of the potential between the edge nodes, and the exact line integral is the
 * difference of the potential between the path end points.
 *
 * For such a field the discrete line integral is the difference of the
 * element's discrete potential at the path end points, so the error is an
 * interpolation error: O(h^2) for the lowest order element, O(h^3) for the
 * order 2 element with its biquadratic block geometry. With the piecewise
 * bilinear cell geometry instead, the order 2 element is limited to O(h^2) by
 * the kink of the four cell maps inside each block.
 */

// cube faces: outward normal n and in-face axes e1, e2 with e1 x e2 = n
const int FACE_N[6][3] =  {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
const int FACE_E1[6][3] = {{0, 1, 0}, {0, 0, 1},  {0, 0, 1}, {1, 0, 0},  {1, 0, 0}, {0, 1, 0}};
const int FACE_E2[6][3] = {{0, 0, 1}, {0, 1, 0},  {1, 0, 0}, {0, 0, 1},  {0, 1, 0}, {1, 0, 0}};

void lonLatToXyz(double lon, double lat, double xyz[]) {
    double lam = lon*M_PI/180., the = lat*M_PI/180.;
    xyz[0] = cos(the)*cos(lam);
    xyz[1] = cos(the)*sin(lam);
    xyz[2] = sin(the);
}

// smooth potential on the sphere
double potential(const double p[]) {
    return p[0]*p[1] + 0.5*p[2]*p[2] - 0.3*p[0]*p[2] + 0.7*p[1] + 0.2*sin(2.*p[0] + p[2]);
}

double potentialLonLat(double lon, double lat) {
    double p[3];
    lonLatToXyz(lon, lat, p);
    return potential(p);
}

/**
 * Equiangular cubed sphere with n x n cells per panel (n even) as UGRID data
 * @param n number of cells along each panel edge
 * @param lonLat node coordinates (lon, lat, 0) in degrees (output)
 * @param xyz node Cartesian coordinates on the unit sphere (output)
 * @param f2n face-node connectivity (output)
 * @param e2n edge-node connectivity (output)
 */
void createCubedSphere(int n, std::vector<double>& lonLat, std::vector<double>& xyz,
                       std::vector<std::size_t>& f2n, std::vector<std::size_t>& e2n) {

    // nodes are points of the integer lattice [0, n]^3 on the cube surface
    std::map< std::array<int, 3>, std::size_t > nodeIds;
    auto nodeId = [&](const std::array<int, 3>& lattice) {
        auto it = nodeIds.find(lattice);
        if (it != nodeIds.end()) return it->second;
        std::size_t id = nodeIds.size();
        nodeIds[lattice] = id;
        double p[3], norm = 0.;
        for (int d = 0; d < 3; ++d) {
            p[d] = tan(0.25*M_PI*(2.0*lattice[d]/n - 1.0));
            norm += p[d]*p[d];
        }
        norm = sqrt(norm);
        for (int d = 0; d < 3; ++d) {
            p[d] /= norm;
            xyz.push_back(p[d]);
        }
        lonLat.push_back(atan2(p[1], p[0])*180./M_PI);
        lonLat.push_back(asin(p[2])*180./M_PI);
        lonLat.push_back(0.);
        return id;
    };

    std::map< std::pair<std::size_t, std::size_t>, std::size_t > edgeIds;
    auto addEdge = [&](std::size_t n0, std::size_t n1) {
        auto key = std::make_pair(std::min(n0, n1), std::max(n0, n1));
        if (edgeIds.find(key) == edgeIds.end()) {
            edgeIds[key] = edgeIds.size();
            e2n.push_back(n0);
            e2n.push_back(n1);
        }
    };

    const int h = n/2;
    const int di[] = {0, 1, 1, 0};
    const int dj[] = {0, 0, 1, 1};
    for (int face = 0; face < 6; ++face) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                std::size_t ids[4];
                for (int k = 0; k < 4; ++k) {
                    std::array<int, 3> lattice;
                    for (int d = 0; d < 3; ++d) {
                        lattice[d] = h + h*FACE_N[face][d] + (i + di[k] - h)*FACE_E1[face][d] +
                                     (j + dj[k] - h)*FACE_E2[face][d];
                    }
                    ids[k] = nodeId(lattice);
                }
                for (int k = 0; k < 4; ++k) {
                    f2n.push_back(ids[k]);
                    addEdge(ids[k], ids[(k + 1) % 4]);
                }
            }
        }
    }
}

typedef std::vector<double> Path; // lon, lat pairs in degrees

// element variants
enum Method {ORDER1 = 0, ORDER2_BILINEAR_CELLS, ORDER2};
const char* METHOD_NAMES[] = {"order 1", "order 2, bilinear cells", "order 2"};

/**
 * @return line integral errors for each path
 */
std::vector<double> lineIntegralErrors(int n, Method method, const std::vector<Path>& paths) {

    std::vector<double> lonLat, xyz;
    std::vector<std::size_t> f2n, e2n;
    createCubedSphere(n, lonLat, xyz, f2n, e2n);
    const std::size_t numEdges = e2n.size()/2;

    Grid_t* grid = NULL;
    mnt_grid_new(&grid);
    mnt_grid_setFlags(&grid, 1, 1, 1);
    int ier = mnt_grid_loadFromUgrid2DData(&grid, f2n.size()/4, numEdges, lonLat.size()/3,
                                           &lonLat[0], &f2n[0], &e2n[0]);
    assert(ier == 0);

    // exact edge integrals of the gradient, from the first to the second edge node
    std::vector<double> data(numEdges);
    for (std::size_t e = 0; e < numEdges; ++e) {
        data[e] = potential(&xyz[3*e2n[2*e + 1]]) - potential(&xyz[3*e2n[2*e]]);
    }

    PolylineIntegral_t* pli = NULL;
    mnt_polylineintegral_new(&pli);
    mnt_polylineintegral_setGrid(&pli, grid);
    ier = mnt_polylineintegral_setOrder(&pli, (method == ORDER1)? 1: 2);
    assert(ier == 0);
    if (method == ORDER2_BILINEAR_CELLS) {
        WhitneyElementMacro2* elem = new WhitneyElementMacro2();
        assert(elem->build(grid, false) == 0);
        assert(elem->getNumBiquadraticBlocks() == 0);
        delete pli->element;
        pli->element = elem;
    }
    else if (method == ORDER2) {
        // all blocks but the 4 around each pole (the pole node has a different
        // longitude in each cell) use the biquadratic map
        const WhitneyElementMacro2* elem = dynamic_cast<const WhitneyElementMacro2*>(pli->element);
        assert(elem);
        assert(elem->getNumBiquadraticBlocks() + 8 == elem->getNumBlocks());
    }
    mnt_polylineintegral_buildLocator(&pli, 100, 360., 1, 0);

    std::vector<double> errors;
    for (const Path& path : paths) {
        int npts = (int) path.size()/2;
        std::vector<double> pts(3*npts, 0.);
        for (int i = 0; i < npts; ++i) {
            pts[3*i] = path[2*i];
            pts[3*i + 1] = path[2*i + 1];
        }
        ier = mnt_polylineintegral_computeWeights(&pli, npts, &pts[0], 0);
        assert(ier == 0);
        double result;
        ier = mnt_polylineintegral_getIntegral(&pli, &data[0], MNT_UNIQUE_EDGE_DATA, &result);
        assert(ier == 0);
        double exact = potentialLonLat(path[2*(npts - 1)], path[2*(npts - 1) + 1]) -
                       potentialLonLat(path[0], path[1]);
        errors.push_back(std::abs(result - exact));
    }

    mnt_polylineintegral_del(&pli);
    mnt_grid_del(&grid);
    return errors;
}

/**
 * Random polylines away from the poles and the dateline, crossing panel edges
 */
std::vector<Path> randomPaths(int numPaths, int numSegments) {
    std::mt19937 gen(1234);
    std::uniform_real_distribution<double> lon(-150., 150.), lat(-70., 70.);
    std::vector<Path> paths;
    for (int i = 0; i < numPaths; ++i) {
        Path path;
        for (int k = 0; k <= numSegments; ++k) {
            path.push_back(lon(gen));
            path.push_back(lat(gen));
        }
        paths.push_back(path);
    }
    return paths;
}

/**
 * Random closed polylines (16-gons of radius 0.5 to 8 degrees)
 */
std::vector<Path> randomLoops(int numLoops) {
    std::mt19937 gen(5);
    std::uniform_real_distribution<double> lon(-150., 150.), lat(-60., 60.), rad(0.5, 8.);
    std::vector<Path> loops;
    for (int i = 0; i < numLoops; ++i) {
        double lo = lon(gen), la = lat(gen), r = rad(gen);
        Path path;
        for (int k = 0; k < 16; ++k) {
            path.push_back(lo + r*cos(2*M_PI*k/16.));
            path.push_back(la + r*sin(2*M_PI*k/16.));
        }
        path.push_back(path[0]);
        path.push_back(path[1]);
        loops.push_back(path);
    }
    return loops;
}

void testConvergence() {

    const std::vector<int> resolutions = {8, 16, 32, 64};
    const std::vector<Path> paths = randomPaths(30, 3);

    // RMS and max errors at each resolution, for each method
    double rms[3][4], mx[3][4];
    for (int method : {ORDER1, ORDER2_BILINEAR_CELLS, ORDER2}) {
        printf("%-24s     RMS error (order)      max error (order)\n", METHOD_NAMES[method]);
        for (std::size_t ir = 0; ir < resolutions.size(); ++ir) {
            std::vector<double> errors = lineIntegralErrors(resolutions[ir], (Method) method, paths);
            double r = 0., m = 0.;
            for (double e : errors) {
                r += e*e;
                m = std::max(m, e);
            }
            rms[method][ir] = sqrt(r/errors.size());
            mx[method][ir] = m;
            printf("  C%-4d                  %10.3e", resolutions[ir], rms[method][ir]);
            if (ir > 0) printf(" (%4.2f)", log2(rms[method][ir - 1]/rms[method][ir])); else printf("       ");
            printf("  %10.3e", mx[method][ir]);
            if (ir > 0) printf(" (%4.2f)", log2(mx[method][ir - 1]/mx[method][ir]));
            printf("\n");
        }
    }

    // average order from C16 to C64
    auto order = [](const double e[]) { return 0.5*log2(e[1]/e[3]); };
    printf("orders C16 -> C64 (RMS/max): order 1 %4.2f/%4.2f, order 2 bilinear cells %4.2f/%4.2f, order 2 %4.2f/%4.2f\n",
           order(rms[ORDER1]), order(mx[ORDER1]),
           order(rms[ORDER2_BILINEAR_CELLS]), order(mx[ORDER2_BILINEAR_CELLS]),
           order(rms[ORDER2]), order(mx[ORDER2]));

    assert(order(rms[ORDER1]) > 1.8 && order(mx[ORDER1]) > 1.7);
    // the bilinear cell geometry limits the order 2 element to O(h^2)
    assert(order(rms[ORDER2_BILINEAR_CELLS]) > 1.8 && order(rms[ORDER2_BILINEAR_CELLS]) < 2.5);
    // O(h^3) with the biquadratic block geometry
    assert(order(rms[ORDER2]) > 2.8 && order(mx[ORDER2]) > 2.7);
    assert(rms[ORDER2][3] < 0.1*rms[ORDER1][3]);
    printf("testConvergence: OK\n");
}

void testClosedLoops() {
    // the discrete field is a gradient: loop integrals vanish, up to the O(h^4)
    // mismatch between the straight cell edges and the curved block edges
    const std::vector<Path> loops = randomLoops(100);
    double mx[3] = {0., 0., 0.};
    for (int method : {ORDER1, ORDER2_BILINEAR_CELLS, ORDER2}) {
        for (double e : lineIntegralErrors(32, (Method) method, loops)) {
            mx[method] = std::max(mx[method], e);
        }
        printf("  C32 %-24s max |loop integral| %10.3e\n", METHOD_NAMES[method], mx[method]);
    }
    assert(mx[ORDER1] < 1.e-9);
    assert(mx[ORDER2_BILINEAR_CELLS] < 1.e-9);
    assert(mx[ORDER2] < 1.e-5);
    printf("testClosedLoops: OK\n");
}

int main(int argc, char** argv) {
    setvbuf(stdout, NULL, _IONBF, 0);
    testConvergence();
    testClosedLoops();
    return 0;
}
