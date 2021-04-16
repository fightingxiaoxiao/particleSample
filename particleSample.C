/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Chen Xiaoxiao
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    particleSample

Group
    grpPostProcessingUtilities

Description
    //Generate a legacy VTK file of particle tracks for cases that were
    computed using a tracked-parcel-type cloud.//

\*---------------------------------------------------------------------------*/
#include <map>
#include <cmath>
#include <fstream>

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "OFstream.H"
#include "writer.H"

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"

#include "basicKinematicCloud.H"
#define basicKinematicTypeCloud basicKinematicCloud

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote(
        "Generate the transport rate of particles for cases.");

    timeSelector::addOptions();
#include "addRegionOption.H"

#include "setRootCase.H"

#include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#include "createNamedMesh.H"

#include "readDictProp.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //const fileName vtkPath(runTime.rootPath()/runTime.globalCaseName()/"VTK");
    //mkDir(vtkPath);

    // get the particle list length on each processor
    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timeI)
    {
        if (timeI % sampleFrequency != 0)
        {
            continue;
        }
        runTime.setTime(timeDirs[timeI], timeI);
        Info << "Time = " << runTime.timeName() << endl;

        Info << "    Reading particle positions" << endl;

#include "readFields.H"

        Info << "    Read " << returnReduce(kinematicCloud.size(), sumOp<label>())
             << " particles" << endl;

        for (const auto &p : kinematicCloud)
        {
            const label origId = p.origId();
            const label origProc = p.origProc();

            if (origProc >= maxIds.size())
            {
                maxIds.setSize(origProc + 1, -1);
            }

            maxIds[origProc] = max(maxIds[origProc], origId);
        }
    }

    label maxNProcs = returnReduce(maxIds.size(), maxOp<label>());

    Info << "Detected particles originating from " << maxNProcs
         << " processors." << nl << endl;

    maxIds.setSize(maxNProcs, -1);

    Pstream::listCombineGather(maxIds, maxEqOp<label>());
    Pstream::listCombineScatter(maxIds);

    labelList numIds = maxIds + 1;

    Info << nl << "Particle statistics:" << endl;
    forAll(maxIds, proci)
    {
        Info << "    Found " << numIds[proci] << " particles originating"
             << " from processor " << proci << endl;
    }
    Info << nl << endl;

    // Calculate starting ids for particles on each processor
    labelList startIds(numIds.size(), Zero);
    for (label i = 0; i < numIds.size() - 1; i++)
    {
        startIds[i + 1] += startIds[i] + numIds[i];
    }
    label nParticle = startIds.last() + numIds[startIds.size() - 1];

    Info << "\nSampling " << nParticle << " particles for KinematicCloud " << nl << endl;

    std::map<label, vector> allPositionDict;
    std::map<label, vector> _allPositionDict = {};

    std::map<label, scalar> allRhoDict;
    std::map<label, scalar> _allRhoDict = {};

    std::map<label, scalar> allDDict;
    std::map<label, scalar> _allDDict = {};

    std::map<label, label> allnParticleDict;
    std::map<label, label> _allnParticleDict = {};

    std::ofstream totalMassRateOut("./postProcessing/totalMassRate");
    totalMassRateOut << "Time, MassRate" << endl;
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        if (timeI % sampleFrequency != 0)
        {
            continue;
        }
        Info << "Time = " << runTime.timeName() << endl;

        List<pointField> allPositions(Pstream::nProcs());

        List<scalarField> allD(Pstream::nProcs());
        List<scalarField> allRho(Pstream::nProcs());

        List<labelField> allnParticle(Pstream::nProcs());
        List<labelField> allOrigIds(Pstream::nProcs());
        List<labelField> allOrigProcs(Pstream::nProcs());

        // Read particles. Will be size 0 if no particles.
        Info << "    Reading particle properties" << endl;
#include "readFields.H"

        // Collect the track data on all processors that have positions
        allPositions[Pstream::myProcNo()].setSize(
            kinematicCloud.size(),
            point::zero);

        allD[Pstream::myProcNo()].setSize(kinematicCloud.size(), Zero);
        allRho[Pstream::myProcNo()].setSize(kinematicCloud.size(), Zero);

        allnParticle[Pstream::myProcNo()].setSize(kinematicCloud.size(), Zero);
        allOrigIds[Pstream::myProcNo()].setSize(kinematicCloud.size(), Zero);
        allOrigProcs[Pstream::myProcNo()].setSize(kinematicCloud.size(), Zero);

        label i = 0;
        for (const auto &p : kinematicCloud)
        {
            allPositions[Pstream::myProcNo()][i] = p.position();

            allD[Pstream::myProcNo()][i] = p.d();
            allRho[Pstream::myProcNo()][i] = p.rho();

            allnParticle[Pstream::myProcNo()][i] = p.nParticle();
            allOrigIds[Pstream::myProcNo()][i] = p.origId();
            allOrigProcs[Pstream::myProcNo()][i] = p.origProc();

            ++i;
        }

        // Collect the data on the master processor
        Pstream::gatherList(allPositions);
        Pstream::gatherList(allD);
        Pstream::gatherList(allRho);
        Pstream::gatherList(allnParticle);
        Pstream::gatherList(allOrigIds);
        Pstream::gatherList(allOrigProcs);

        if (Pstream::master())
        {

            forAll(allPositions, proci)
            {
                forAll(allPositions[proci], i)
                {
                    label globalId =
                        startIds[allOrigProcs[proci][i]] + allOrigIds[proci][i];

                    allPositionDict[globalId] = allPositions[proci][i];

                    allRhoDict[globalId] = allRho[proci][i];
                    allDDict[globalId] = allD[proci][i];

                    allnParticleDict[globalId] = allnParticle[proci][i];
                }
            }

            scalar particleMassRate = 0.;
            if (_allPositionDict.size() == 0)
            {
                Info << "This is the first data." << endl;
            }
            else
            {
                for (auto &keyTwice : allPositionDict)
                {
                    auto key = keyTwice.first;
                    auto position = allPositionDict[key];
                    auto _position = _allPositionDict[key];

                    if (mag(position - _position) > limitMoveDistanceInOneSample)
                        continue;

                    if (_position[directionIndex] < samplePosition && position[directionIndex] >= samplePosition)
                    {
                        particleMassRate += allnParticleDict[key] * 4 / 3 * M_PI * Foam::pow(allDDict[key] / 2., 3.) * allRhoDict[key];
                    }
                    else if (_position[directionIndex] > samplePosition && position[directionIndex] <= samplePosition)
                    {
                        particleMassRate -= allnParticleDict[key] * 4 / 3 * M_PI * Foam::pow(allDDict[key] / 2., 3.) * allRhoDict[key];
                    }
                }
                Info << "The particle flow rate is " << particleMassRate / runTime.deltaTValue() / sampleFrequency << endl;
                totalMassRateOut << runTime.timeName() << ", " << particleMassRate / runTime.deltaTValue() / sampleFrequency << std::endl;
            }

            Info << "\n\n"
                 << endl;

            _allPositionDict = allPositionDict;
            _allRhoDict = allRhoDict;
            _allDDict = allDDict;
            _allnParticleDict = allnParticleDict;
        }
    }

    totalMassRateOut.close();
    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
