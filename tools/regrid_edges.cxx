#include <mntRegridEdges.h>
#include <mntGrid.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <iostream>
#include <limits>

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Regrid an edge centred field.");
    args.set("-s", std::string(""), "Source grid file in VTK format");
    args.set("-v", std::string("edge_integrated_velocity"), "Specify edge staggered field variable name in source VTK file");
    args.set("-d", std::string(""), "Destination grid file in VTK format");
    args.set("-o", std::string(""), "Write interpolation weights to file");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        std::string srcFile = args.get<std::string>("-s");
        std::string dstFile = args.get<std::string>("-d");
        std::string weightsFile = args.get<std::string>("-o");

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
        RegridEdges_t* rge = NULL;
        mnt_regridedges_new(&rge);
        ier = mnt_regridedges_setSrcGrid(&rge, sg);
        if (ier != 0) return 1;
        ier = mnt_regridedges_setDstGrid(&rge, dg);
        if (ier != 0) return 2;
        ier = mnt_regridedges_build(&rge);
        if (ier != 0) return 3;

        if (weightsFile.size() != 0) {
            std::cout << "INFO saving weights in file " << weightsFile << '\n';
            ier = mnt_regridedges_dump(&rge, weightsFile.c_str());
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
                mnt_regridedges_getNumDstCells(&rge, &numDstCells);
                mnt_regridedges_getNumEdgesPerCell(&rge, &numEdgesPerCell);
                std::vector<double> dstData(numDstCells * numEdgesPerCell);

                // regrid
                mnt_regridedges_applyWeights(&rge, srcData, &dstData[0]);

                // compute loop integrals for each cell
                double minAbsLoop = std::numeric_limits<double>::max();
                double maxAbsLoop = - std::numeric_limits<double>::max();
                double avgAbsLoop = 0.0;
                for (size_t i = 0; i < numDstCells; ++i) {
                    size_t k = i*numEdgesPerCell;
                    double loop = 0.0;
                    for (size_t j = 0; j < numEdgesPerCell; ++j) {
                        loop += dstData[k + j];
                    }
                    loop = std::abs(loop);
                    minAbsLoop = std::min(loop, minAbsLoop);
                    maxAbsLoop = std::max(loop, maxAbsLoop);
                    avgAbsLoop += loop;
                }
                avgAbsLoop /= double(numDstCells);
                std::cout << "Min/avg/max cell loop integrals: " << minAbsLoop << "/" << avgAbsLoop << "/" << maxAbsLoop << '\n';
            }
        }

        // cleanup
        mnt_grid_del(&dstGrid);
        mnt_grid_del(&srcGrid);
        mnt_regridedges_del(&rge);

    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
    }

    return 0;
}