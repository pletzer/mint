#include <GrExprParser.h>
#include <GrExprAdaptor.h>
#include <mntGrid.h>
#include <mntPolysegmentIter.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
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

        // hold the line cooridnates
        Vec* xs = xExpr.eval();
        Vec* ys = yExpr.eval();

        for (size_t i = 0; i < (*xs).size(); ++i) {
            std::cout << i << " x = " << (*xs)[i] << ", " << " y = " << (*ys)[i] << '\n';
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

        // build locator
        vtkCellLocator* loc = vtkCellLocator::New();
        loc->SetDataSet(grid);
        loc->BuildLocator();

        double flux = 0.0;

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

                std::vector<double> dxi({xib[0] - xia[0], xib[1] - xia[1]});
                std::vector<double> xiMid({0.5*(xia[0] + xib[0]), 0.5*(xia[1] + xib[1])});

                // compute the weight contributions from each src edge
                double ws[] = {+ dxi[0] * (1.0 - xiMid[1]) * coeff,
                               + dxi[1] * (0.0 + xiMid[0]) * coeff,
                               - dxi[0] * (0.0 + xiMid[1]) * coeff,
                               - dxi[1] * (1.0 - xiMid[0]) * coeff};

                // project
                const size_t numEdges = 4;
                for (size_t j = 0; j < numEdges; ++j) {
                    flux += ws[j]*srcData[cellId*numEdges + j];
                }

                //std::cout << "cell " << cellId << " t = " << ta  << " -> " << tb << 
                //            " xi = " << xia[0] << ',' << xia[1] << " -> " << xib[0] << ',' << xib[1] << 
                //            " coeff = " << coeff << '\n';

                polyseg.next();
            }

            double tTotal = polyseg.getIntegratedParamCoord();
            if (std::abs(tTotal - 1.0) > 1.e-10) {
                std::cerr << "Warning: total t of segment " << ip0
                          << " != 1 (diff=" << tTotal - 1.0 << ")\n";
            }

        }

        printf("Flux across the line: %5.18f\n", flux);

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
