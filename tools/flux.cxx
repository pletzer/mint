#include <GrExprParser.h>
#include <GrExprAdaptor.h>
#include <mntGrid.h>
#include <mntPolylineIntegral.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <string>
#include <iostream>
#include <cstdio>

#define NDIMS 3

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Project a line segment.");
    args.set("-i", std::string(""), "Source grid file in VTK format");
    args.set("-v", std::string("edge_integrated_velocity"), "Specify edge staggered field variable name in VTK file");
    args.set("-xline", std::string("180./4. + 0.2*180.*cos(5*pi*t/3. - 0.2)/pi"), "Parametric x coordinate (lon in [rad]) expression of 0 <= t <= 1");
    args.set("-yline", std::string("180./4. + 0.2*180.*sin(5*pi*t/3. - 0.2)/pi"), "Parametric y coordinate (lat in [rad]) expression of 0 <= t <= 1");
    args.set("-nline", 2, "Number of points defining the line (>= 2)");
    args.set("-counter", false, "Whether the edge values correspond to edges oriented in counterclockwise direction");
    args.set("-P", 0.0, "Periodicity length of x (0 if not periodic)");


    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        
        // get the variable name of the field to project
        std::string varname = args.get<std::string>("-v");

        // construct the projection line from the x and y expressions
        std::string xLineExpr = args.get<std::string>("-xline");
        std::string yLineExpr = args.get<std::string>("-yline");
        int nline = args.get<int>("-nline");
        Vec ts(nline);
        ts.space(0., 1.);

        GrExprParser xExpr(ts.size(), GrExprAdaptor(xLineExpr).getPrefixExpr());
        xExpr.defineVariable("t", &ts);
        
        GrExprParser yExpr(ts.size(), GrExprAdaptor(yLineExpr).getPrefixExpr());
        yExpr.defineVariable("t", &ts);

        // hold the line coordinates
        Vec* xs = xExpr.eval();
        Vec* ys = yExpr.eval();
        size_t npts = (*xs).size();
        std::vector<double> xyz(npts * 3); //3D

        for (size_t i = 0; i < npts; ++i) {
            double x = (*xs)[i];
            double y = (*ys)[i];
            xyz[3*i + 0] = x;
            xyz[3*i + 1] = y;
            xyz[3*i + 2] = 0.; 
            std::cout << i << " x = " << x << ", " << " y = " << y << '\n';
        }

        std::string srcFile = args.get<std::string>("-i");
        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-i)\n";
            return 1;
        }

        if (ts.size() < 2) {
            std::cerr << "ERROR: need at least two points to create a path\n";
            return 2;
        }

        std::cout << "Start point " << (*xs)[0] << " [deg east], " << (*ys)[0] << " [deg north] -> End point " 
                  << (*xs)[nline - 1] << " [deg east], " << (*ys)[nline -1] << " [deg north]\n";

        // read/build the src grid
        Grid_t* srcGrid = NULL;
        mnt_grid_new(&srcGrid);
        mnt_grid_load(&srcGrid, srcFile.c_str());

        // get the grid pointers
        vtkUnstructuredGrid* grid = NULL;
        mnt_grid_get(&srcGrid, &grid);

        std::cout << "Grid: no of cells " << grid->GetNumberOfCells() << " no of points " << grid->GetNumberOfPoints() << '\n';

        // get the velocity field
        vtkCellData* cellData = grid->GetCellData();
        vtkAbstractArray* aa = cellData->GetAbstractArray(varname.c_str());
        if (!aa) {
            std::cerr << "ERROR: must specify a a valid field to project (-v)\n";
            return 3;
        }
        double* srcData = (double*) aa->GetVoidPointer(0);

        // create the flux calculator object
        PolylineIntegral_t* fluxCalc = NULL;
        mnt_polylineintegral_new(&fluxCalc);
        int counterclock = 0;
        if (args.get<bool>("-counter")) {
          std::cout << "Warning: edge values are assumed to be stored in counterclockwise direction\n";
          counterclock = 1;
        }
        double periodX = args.get<double>("-P");
        ier = mnt_polylineintegral_build(&fluxCalc, srcGrid, (int) npts, &xyz[0], 
                                         counterclock, periodX);
        if (ier != 0) {
            std::cerr << "ERROR: after calling mnt_lineintegral_build ier = " << ier << '\n';
        }

        double flux = 0.0;
        ier = mnt_polylineintegral_getIntegral(&fluxCalc, (const double*) srcData, &flux);

        printf("Flux across the line: %5.18f\n", flux);

        // cleanup
        mnt_polylineintegral_del(&fluxCalc);
        mnt_grid_del(&srcGrid);

    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
    }

    return 0;
}
