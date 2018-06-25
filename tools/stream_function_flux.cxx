#include <GrExprParser.h>
#include <GrExprAdaptor.h>
#include <mntGrid.h>
#include <mntPolysegmentIter.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <string>
#include <iostream>
#include <cstdio>

#define NDIMS 3

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Project a line segment.");
    args.set("-i", std::string(""), "Source grid file in VTK format");
    args.set("-xline", std::string("pi/4. + 0.2*cos(5*pi/3. - 0.2"), "Parametric x coordinate (lon in [rad]) expression of 0 <= t <= 1");
    args.set("-yline", std::string("pi/4. + 0.2*sin(5*pi/3. - 0.2"), "Parametric y coordinate (lat in [rad]) expression of 0 <= t <= 1");
    args.set("-nline", 2, "Number of points defining the line (>= 2)");


    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {

        std::string xLineExpr = args.get<std::string>("-xline");
        std::string yLineExpr = args.get<std::string>("-yline");
        int nline = args.get<int>("-nline");
        Vec ts(nline);
        ts.space(0., 1.);

        GrExprParser xExpr(ts.size(), GrExprAdaptor(xLineExpr).getPrefixExpr());
        xExpr.defineVariable("t", &ts);
        Vec* xs = xExpr.eval();

        GrExprParser yExpr(ts.size(), GrExprAdaptor(yLineExpr).getPrefixExpr());
        yExpr.defineVariable("t", &ts);
        Vec* ys = yExpr.eval();

        std::string srcFile = args.get<std::string>("-i");
        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-i)\n";
            return 1;
        }

        if (ts.size() < 2) {
            std::cerr << "ERROR: need at least two points to create a path\n";
            return 2;
        }

        // read/build the src grid
        Grid_t* srcGrid = NULL;
        mnt_grid_new(&srcGrid);
        mnt_grid_load(&srcGrid, srcFile.c_str());

        // get the grid pointers
        vtkUnstructuredGrid* grid = NULL;
        mnt_grid_get(&srcGrid, &grid);

        std::cout << "Grid: no of cells " << grid->GetNumberOfCells() << " no of points " << grid->GetNumberOfPoints() << '\n';

        // build locator
        vtkCellLocator* loc = vtkCellLocator::New();
        loc->SetDataSet(grid);
        loc->BuildLocator();

        // iterate over segments
        for (size_t ip0 = 0; ip0 < ts.size() - 1; ++ip0) {

            double p0[] = {(*xs)[ip0], (*ys)[ip0], 0.};

            size_t ip1 = ip0 + 1;
            double p1[] = {(*xs)[ip1], (*ys)[ip1], 0.};

            PolysegmentIter polyseg(grid, loc, p0, p1);

            size_t numSubSegs = polyseg.getNumberOfSegments();
            polyseg.reset();

            // iterate over the sub-segments 
            for (size_t i = 0; i < numSubSegs; ++i) {

                vtkIdType cellId = polyseg.getCellId();
                double ta = polyseg.getBegLineParamCoord();
                double tb = polyseg.getEndLineParamCoord();
                const std::vector<double>& xia = polyseg.getBegCellParamCoord();
                const std::vector<double>& xib = polyseg.getEndCellParamCoord();
                double coeff = polyseg.getCoefficient();

                std::cout << "cell " << cellId << " t = " << ta  << " -> " << tb << 
                            " xi = " << xia[0] << ',' << xia[1] << " -> " << xib[0] << ',' << xib[1] << 
                            " coeff = " << coeff << '\n';
                polyseg.next();
            }
            double tTotal = polyseg.getIntegratedParamCoord();

        }

        // cleanup
        mnt_grid_del(&srcGrid);
        loc->Delete();


    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
    }

    return 0;
}