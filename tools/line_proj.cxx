#include <mntLogger.h>
#include <mntVecN.h>
#include <mntPolylineParser.h>
#include <mntGrid.h>
#include <mntPolylineIntegral.h>
#include <CmdLineArgParser.h>
#include <vtkCellData.h>
#include <string>
#include <iostream>
#include <cstdio>

/**
 * Compute line integral of a 2d edge field
 *
 */

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Project a line segment.");
    args.set("-i", std::string(""), "Source grid file in VTK format");
    args.set("-p", std::string("(0., 0.),(360., 0.)"), "Points defining the path.");
    args.set("-v", std::string("edgeData"), "Edge variable name.");
    args.set("-N", 1, "Average number of cells per bucket.");
    args.set("-P", 0.0, "Periodicity length of x (0 if not periodic)");
    args.set("-verbose", false, "Turn on verbosity");
    args.set("-counter", false, "Whether the edge values correspond to edges oriented in counterclockwise direction");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {

        std::string srcFile = args.get<std::string>("-i");

        const size_t ndims = 2;
        PolylineParser pp(ndims);
        pp.parse(args.get<std::string>("-p"));
        const std::vector< Vec3 >& points = pp.getPoints();
        size_t npts = points.size();
        // PolylineIntegrals expects a flat array
        std::vector<double> xyz(npts * 3); // 3D
        std::cout << "Path:\n";
        for (size_t i = 0; i < npts; ++i) {
            std::cout << i << ": " << points[i] << '\n';
            for (size_t j = 0; j < points[i].size(); ++j) {
                xyz[3*i + j] = points[i][j];
            }
        }

        // checks
        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-i)\n";
            return 1;
        }
        if (npts < 2) {
            std::cerr << "ERROR: must have at least two points\n";
            return 2;
        }
        // compute length of path and check it is != 0
        double pathLengthSq = 0;
        for (size_t i = 0; i < npts - 1; ++i) {
            double dx = points[i + 1][0] - points[i + 0][0];
            double dy = points[i + 1][1] - points[i + 0][1];
            double dz = points[i + 1][2] - points[i + 0][2];
            pathLengthSq += dx*dx + dy*dy + dz*dz;
        }
        if (pathLengthSq == 0) {
            std::cerr << "ERROR: path length must be > 0\n";
            return 3;
        }

        // read/build the src grid
        Grid_t* srcGrid = NULL;
        mnt_grid_new(&srcGrid);
        mnt_grid_load(&srcGrid, srcFile.c_str());

        // get the grid pointers
        vtkUnstructuredGrid* vgrid = NULL;
        mnt_grid_get(&srcGrid, &vgrid);

        // build the interpolator
        // create the flux calculator object
        PolylineIntegral_t* fluxCalc = NULL;
        mnt_polylineintegral_new(&fluxCalc);
        int counterclock = 0;
        if (args.get<bool>("-counter")) {
          std::cout << "Warning: edge values are assumed to be stored in counter clockwise direction\n";
          counterclock = 1;
        }
        double periodX = args.get<double>("-P");
        ier = mnt_polylineintegral_build(&fluxCalc, srcGrid, (int) npts, &xyz[0], counterclock, periodX);
        if (ier != 0) {
            std::cerr << "ERROR: after calling mnt_lineintegral_build ier = " << ier << '\n';
        }

        std::cout << "no of cells " << vgrid->GetNumberOfCells() << " no of points " << vgrid->GetNumberOfPoints() << '\n';

        // get the velocity field
        vtkCellData* cellData = vgrid->GetCellData();
        std::string varname = args.get<std::string>("-v");
        vtkAbstractArray* aa = cellData->GetAbstractArray(varname.c_str());
        if (!aa) {
            std::cerr << "ERROR: must specify a valid field to project (-v)\n";
            return 3;
        }
        double* srcData = (double*) aa->GetVoidPointer(0);

        double totFlux = 0.0;
        ier = mnt_polylineintegral_getIntegral(&fluxCalc, (const double*) srcData, &totFlux);
        std::cout << "total flux is: " << totFlux << '\n';

        // cleanup
        mnt_polylineintegral_del(&fluxCalc);
        mnt_grid_del(&srcGrid);


    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
        return 3;
    }

    if (args.get<bool>("-verbose")) mnt_printLogMessages();
    std::string logname = "line_proj_logs.txt";
    mnt_writeLogMessages(logname.c_str(), logname.size());

    return 0;
}
