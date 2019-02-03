#include <mntRegridEdges3d.h>
#include <mntGrid.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <iostream>
#include <limits>
#include <cmath>

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Regrid an edge centred field.");
    args.set("-s", std::string(""), "Source grid file in VTK format");
    args.set("-v", std::string("edge_integrated_velocity"), "Specify edge staggered field variable name in source VTK file");
    args.set("-d", std::string(""), "Destination grid file in VTK format");
    args.set("-w", std::string(""), "Write interpolation weights to file");
    args.set("-o", std::string(""), "Specify output VTK file where regridded edge data is saved");
    args.set("-N", 1024, "Average number of cells per bucket");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        std::string srcFile = args.get<std::string>("-s");
        std::string dstFile = args.get<std::string>("-d");
        std::string weightsFile = args.get<std::string>("-w");
        std::string regridFile = args.get<std::string>("-o");

        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-s)\n";
            return 1;
        }
        if (dstFile.size() == 0) {
            std::cerr << "ERROR: must specify a destination grid file (-d)\n";
            return 2;
        }

        vtkUnstructuredGrid* sg = NULL;
        vtkUnstructuredGrid* dg = NULL;

        // read/build the src grid
        Grid_t* srcGrid = NULL;
        mnt_grid_new(&srcGrid);
        mnt_grid_load(&srcGrid, srcFile.c_str());

        // read/build the dst grid
        Grid_t* dstGrid = NULL;
        mnt_grid_new(&dstGrid);
        mnt_grid_load(&dstGrid, dstFile.c_str());

        // get the grid pointers
        mnt_grid_get(&srcGrid, &sg);
        mnt_grid_get(&dstGrid, &dg);

        // compute the interpolation weights
        RegridEdges3d_t* rge = NULL;
        mnt_regridedges3d_new(&rge);
        ier = mnt_regridedges3d_setSrcGrid(&rge, sg);
        if (ier != 0) return 1;
        ier = mnt_regridedges3d_setDstGrid(&rge, dg);
        if (ier != 0) return 2;
        ier = mnt_regridedges3d_build(&rge, args.get<int>("-N"));
        if (ier != 0) return 3;

        if (weightsFile.size() != 0) {
            std::cout << "INFO saving weights in file " << weightsFile << '\n';
            ier = mnt_regridedges3d_dump(&rge, weightsFile.c_str());
        }

        // regrid
        std::string varname = args.get<std::string>("-v");
        vtkCellData* cellData = sg->GetCellData();
        if (varname.size() > 0) {
            vtkAbstractArray* aa = cellData->GetAbstractArray(varname.c_str());

            if (aa) {
                // found the array
                std::cout << "Found variable '" << varname << "'\n";

                double* srcData = (double*) aa->GetVoidPointer(0);


                int numDstCells, numEdgesPerCell;
                mnt_regridedges3d_getNumDstCells(&rge, &numDstCells);
                mnt_regridedges3d_getNumEdgesPerCell(&rge, &numEdgesPerCell);
                std::vector<double> dstData(numDstCells * numEdgesPerCell);

                // regrid
                mnt_regridedges3d_applyWeights(&rge, srcData, &dstData[0]);

                // compute loop integrals for each cell
                std::vector<double> loop_integrals(numDstCells);
                double minAbsLoop = std::numeric_limits<double>::max();
                double maxAbsLoop = - std::numeric_limits<double>::max();
                double avgAbsLoop = 0.0;
                for (size_t i = 0; i < numDstCells; ++i) {
                    size_t k = i*numEdgesPerCell;
                    double loop = 0.0;
                    for (size_t j = 0; j < numEdgesPerCell; ++j) {
                        loop += dstData[k + j];
                    }
                    loop_integrals[i] = loop;
                    loop = std::abs(loop);
                    minAbsLoop = std::min(loop, minAbsLoop);
                    maxAbsLoop = std::max(loop, maxAbsLoop);
                    avgAbsLoop += loop;
                }
                avgAbsLoop /= double(numDstCells);
                std::cout << "Min/avg/max cell loop integrals: " << minAbsLoop << "/" << avgAbsLoop << "/" << maxAbsLoop << '\n';

                if (regridFile.size() > 0) {
                	// attach field to grid so we can save the data in file
                	mnt_grid_attach(&dstGrid, varname.c_str(), 1, &dstData[0]);
                	std::string loop_integral_varname = std::string("loop_integrals_of_") + varname;
                	mnt_grid_attach(&dstGrid, loop_integral_varname.c_str(), 1, &loop_integrals[0]);

                	std::cout << "Writing " << varname << " to " << regridFile << '\n';
                	mnt_grid_dump(&dstGrid, regridFile.c_str());
                }
            }
        }

        // cleanup
        mnt_grid_del(&dstGrid);
        mnt_grid_del(&srcGrid);
        mnt_regridedges3d_del(&rge);

    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
    }

    return 0;
}
