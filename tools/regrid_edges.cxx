#include <mntRegridEdges.h>
#include <mntGrid.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Regrid an edge centred field.");
    args.set("-s", std::string(""), "UGRID source grid file and mesh name, specified as \"filename:meshname\"");
    args.set("-v", std::string(""), "Specify edge staggered field variable name in source UGRID file, varname[@filename:meshname]");
    args.set("-d", std::string(""), "UGRID destination grid file name");
    args.set("-w", std::string(""), "Write interpolation weights to file");
    args.set("-W", std::string(""), "Load interpolation weights from file");
    args.set("-o", std::string(""), "Specify output VTK file where regridded edge data are saved");
    args.set("-O", std::string(""), "Specify output 2D UGRID file where regridded edge data are saved");
    args.set("-S", 1, "Set to zero if you want to disable source grid regularization. This might be required for uniform lon-lat grids");
    args.set("-D", 1, "Set to zero if you want to disable destination grid regularization. This might be required for uniform lon-lat grids");
    args.set("-N", 10, "Average number of cells per bucket");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        std::string srcFile = args.get<std::string>("-s");
        std::string dstFile = args.get<std::string>("-d");
        std::string weightsFile = args.get<std::string>("-w");
        std::string loadWeightsFile = args.get<std::string>("-W");
        std::string regridFile = args.get<std::string>("-o");
        std::string dstEdgeDataFile = args.get<std::string>("-O");

        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-s)\n";
            return 1;
        }
        if (dstFile.size() == 0) {
            std::cerr << "ERROR: must specify a destination grid file (-d)\n";
            return 2;
        }

        RegridEdges_t* rg;
        mnt_regridedges_new(&rg);

        // defaults are suitable for cubed-sphere 
        int fixLonAcrossDateline = 1;
        int averageLonAtPole = 1;
        if (args.get<int>("-S") == 0) {
            fixLonAcrossDateline = 0;
            averageLonAtPole = 0;
            std::cout << "info: no regularization applied to source grid\n";
        }
        ier = mnt_regridedges_setSrcGridFlags(&rg, fixLonAcrossDateline, averageLonAtPole);

        // ...destination griod
        fixLonAcrossDateline = 1;
        averageLonAtPole = 1;
        if (args.get<int>("-D") == 0) {
            fixLonAcrossDateline = 0;
            averageLonAtPole = 0;
            std::cout << "info: no regularization applied to destination grid\n";
        }
        ier = mnt_regridedges_setDstGridFlags(&rg, fixLonAcrossDateline, averageLonAtPole);

        // read the source grid
        ier = mnt_regridedges_loadSrcGrid(&rg, srcFile.c_str(), srcFile.size());
        if (ier != 0) {
            std::cerr << "ERROR: could not read file \"" << srcFile << "\"\n";
            return 3;
        }

        // read the destination grid
        ier = mnt_regridedges_loadDstGrid(&rg, dstFile.c_str(), dstFile.size());
        if (ier != 0) {
            std::cerr << "ERROR: could not read file \"" << dstFile << "\"\n";
            return 4;
        }

        if (loadWeightsFile.size() == 0) {
            ier = mnt_regridedges_build(&rg, args.get<int>("-N"));
            if (ier != 0) return 5;
        
            if (weightsFile.size() != 0) {
                std::cout << "info: saving weights in file " << weightsFile << '\n';
                ier = mnt_regridedges_dumpWeights(&rg, weightsFile.c_str(), (int) weightsFile.size());
            }

        }
        else {
            std::cout << "info: loading weights from file " << loadWeightsFile << '\n';
            ier = mnt_regridedges_loadWeights(&rg, loadWeightsFile.c_str(), (int) loadWeightsFile.size());
            if (ier != 0) return 6;
        }


        // regrid
        size_t numSrcEdges, numDstEdges;
        mnt_regridedges_getNumSrcEdges(&rg, &numSrcEdges);
        mnt_regridedges_getNumDstEdges(&rg, &numDstEdges);
        std::vector<double> srcEdgeData(numSrcEdges);
        std::vector<double> dstEdgeData(numDstEdges);

        std::string varAtFileMesh = args.get<std::string>("-v");
        if (varAtFileMesh.size() > 0) {

            // by default the variable is stored in srcFile
            std::string fileAndMeshName = srcFile;
            std::string vname = varAtFileMesh;

            size_t columnAt = varAtFileMesh.find('@');
            if (columnAt < std::string::npos) {
                // user specified the file and mesh names
                fileAndMeshName = varAtFileMesh.substr(columnAt + 1);
                vname = varAtFileMesh.substr(0, columnAt);
            }

            std::cout << "info: loading field " << vname << " from file \"" << fileAndMeshName << "\"\n";
            ier = mnt_regridedges_loadEdgeField(&rg, 
                                                fileAndMeshName.c_str(), fileAndMeshName.size(),
                                                vname.c_str(), vname.size(),
                                                numSrcEdges, &srcEdgeData[0]);
            if (ier != 0) {
                std::cerr << "ERROR: could not load edge centred data \"" << vname << "\" from file \"" << fileAndMeshName << "\"\n";
                return 6;
            }

            std::cerr << "info: num src edge = " << numSrcEdges << " num dst edges = " << numDstEdges << '\n';

#ifdef DEBUG
            std::vector<double>::iterator dmin = std::min_element(srcEdgeData.begin(), srcEdgeData.end());
            std::vector<double>::iterator dmax = std::max_element(srcEdgeData.begin(), srcEdgeData.end());
            std::cout << "info: min/max src edge data: " << 
                         srcEdgeData[std::distance(srcEdgeData.begin(), dmin)] << "/" <<  
                         srcEdgeData[std::distance(srcEdgeData.begin(), dmax)] << '\n';
#endif

            // apply the weights to the src field to set the dst field
            ier = mnt_regridedges_apply(&rg, &srcEdgeData[0], &dstEdgeData[0]);
            if (ier != 0) {
                std::cerr << "ERROR: could not apply the weights to source edge field \"" << vname << "\" from file \"" << fileAndMeshName << "\"\n";
                return 7;
            }

#ifdef DEBUG
            dmin = std::min_element(dstEdgeData.begin(), dstEdgeData.end());
            dmax = std::max_element(dstEdgeData.begin(), dstEdgeData.end());
            std::cout << "info: min/max dst edge data: " << 
                         dstEdgeData[std::distance(dstEdgeData.begin(), dmin)] << "/" <<  
                         dstEdgeData[std::distance(dstEdgeData.begin(), dmax)] << '\n';
#endif


            // compute loop integrals for each cell
            size_t numDstCells, dstEdgeId;
            int dstEdgeSign;
            mnt_regridedges_getNumDstCells(&rg, &numDstCells);
            int numEdgesPerCell;
            mnt_regridedges_getNumEdgesPerCell(&rg, &numEdgesPerCell);
            std::vector<double> loop_integrals(numDstCells);
            double minAbsLoop = std::numeric_limits<double>::max();
            double maxAbsLoop = - std::numeric_limits<double>::max();
            double avgAbsLoop = 0.0;
            for (size_t dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {
                double loop = 0.0;
                for (int ie = 0; ie < numEdgesPerCell; ++ie) {

                    ier = mnt_grid_getEdgeId(&rg->dstGridObj, dstCellId, ie, &dstEdgeId, &dstEdgeSign);
                    assert(ier == 0);

                    // +1 for ie = 0, 1; -1 for ie = 2, 3
                    int sgn = 1 - 2*(ie/2);
                    loop += sgn * dstEdgeSign * dstEdgeData[dstEdgeId];
                }

                loop_integrals[dstCellId] = loop;
                loop = std::abs(loop);
                minAbsLoop = std::min(loop, minAbsLoop);
                maxAbsLoop = std::max(loop, maxAbsLoop);
                avgAbsLoop += loop;
            }
            avgAbsLoop /= double(numDstCells);
            std::cout << "Min/avg/max cell loop integrals: " << minAbsLoop << "/" << avgAbsLoop << "/" << maxAbsLoop << '\n';

            if (regridFile.size() > 0) {

                // cell by cell data
                std::vector<double> dstCellByCellData(numDstCells * numEdgesPerCell);
                for (size_t dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {
                    for (int ie = 0; ie < 4; ++ie) {
                        ier = mnt_grid_getEdgeId(&rg->dstGridObj, dstCellId, ie, &dstEdgeId, &dstEdgeSign);
                        size_t k = dstCellId*numEdgesPerCell + ie;
                        dstCellByCellData[k] = dstEdgeData[dstEdgeId] * dstEdgeSign;
                    }
                }

                // attach field to grid so we can save the data in file
                //mnt_grid_attach(&rg->dstGridObj, vname.c_str(), numEdgesPerCell, &dstEdgeData[0]);

                std::string loop_integral_varname = std::string("loop_integrals_of_") + vname;
                //mnt_grid_attach(&rg->dstGridObj, loop_integral_varname.c_str(), 1, &loop_integrals[0]);

                std::cout << "info: writing \"" << vname << "\" to VTK file " << regridFile << '\n';
                mnt_grid_dump(&rg->dstGridObj, regridFile.c_str());
            }

            if (dstEdgeDataFile.size() > 0) {
                std::cout << "info: writing \"" << vname << "\" to netCDF file " << dstEdgeDataFile << '\n';
                mnt_regridedges_dumpEdgeField(&rg, dstEdgeDataFile.c_str(), dstEdgeDataFile.size(), 
                                               vname.c_str(), vname.size(), 
                                               numDstEdges, &dstEdgeData[0]);
            }
        }

        // cleanup
        mnt_regridedges_del(&rg);

    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
    }

    return 0;
}
