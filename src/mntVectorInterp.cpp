#include <mntVectorInterp.h>
#include <vtkGenericCell.h>

double crossDotZHat(const Vec3& a, const Vec3& b) {
    return a[0]*b[1] - a[1]*b[0];
}

extern "C" 
int mnt_vectorinterp_new(VectorInterp_t** self) {
    *self = new VectorInterp_t();
    (*self)->grid = NULL;
    (*self)->locator = NULL;
    return 0;
}

extern "C"
int mnt_vectorinterp_del(VectorInterp_t** self) {
    delete *self;
    return 0;
}

extern "C"
int mnt_vectorinterp_set(VectorInterp_t** self, Grid_t* grid, vmtCellLocator* locator) {
    (*self)->grid = grid;
    (*self)->locator = locator;
    return 0;
}

extern "C"
int mnt_vectorinterp_find(VectorInterp_t** self, const double targerPoint[], double tol2) {
    vtkGenericCell *cell;
    double weights[8];
    (*self)->cellId = (*self)->locator->FindCell(targerPoint, tol2, cell, &(*self)->pcoords[0], weights);
    if ((*self)->cellId < 0) {
        // failed
        return 1;
    }
    return 0;
}

extern "C"
int mnt_vectorinterp_getVector(VectorInterp_t** self, const double data[], double vector[]) {

    int ier;

    if ((*self)->cellId < 0) {
        // failed, maybe the cell was not found?
        // be sure to call mnt_vectorinterp_find before
        return 1;
    }

    // get the cell vertices, assumes edges are in the positive pcoords direction
    Vec3 v0, v1, v2, v3;
    ier = mnt_grid_getPoints(&(*self)->grid, (*self)->cellId, 0, &v0[0], &v1[0]);
    ier = mnt_grid_getPoints(&(*self)->grid, (*self)->cellId, 2, &v3[0], &v2[0]);

    double xsi = (*self)->pcoords[0];
    double eta = (*self)->pcoords[1];
    double isx = 1.0 - xsi;
    double ate = 1.0 - eta;

    // compute the Jacobians attached at each vertex
    double a013 = crossDotZHat(v1 - v0, v3 - v0);
    double a120 = crossDotZHat(v2 - v1, v0 - v1);
    double a231 = crossDotZHat(v3 - v2, v1 - v2);
    double a302 = crossDotZHat(v0 - v3, v2 - v3);

    // interpolate, assume edges in positive pcoords direction
    double jac = isx*ate*a013 + xsi*ate*a120 + xsi*eta*a231 + isx*eta*a302;
    Vec3 drdXsi = ate*(v1 - v0) + eta*(v2 - v3);
    Vec3 drdEta = isx*(v3 - v0) + xsi*(v2 - v1);

    Vec3 gradXsi; 
    gradXsi[0] = + drdEta[1]/jac;
    gradXsi[1] = - drdEta[0]/jac;
    gradXsi[2] = 0.0;

    Vec3 gradEta; 
    gradEta[0] = - drdXsi[1]/jac;
    gradEta[1] = + drdXsi[0]/jac;
    gradEta[2] = 0.0;

    // interpolate
    size_t k = 4*(*self)->cellId;
    double data0 = data[k + 0];
    double data1 = data[k + 1];
    double data2 = data[k + 2];
    double data3 = data[k + 3];
    for (std::size_t i = 0; i < 3; ++i) {
        vector[i] = (data0*ate + data2*eta)*gradXsi[i] + 
                    (data3*isx + data1*xsi)*gradEta[i];
    }

    return 0;
}


