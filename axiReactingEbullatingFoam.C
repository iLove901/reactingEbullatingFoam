/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | 
   \\    /   O peration     | 	OpenFOAM: The Open Source CFD Toolbox
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    reactingEbullatingFoam

Description
	Solver for a system of 2 compressible fluid phases and 1 solid particle phase. 
	Continuous phase is liquid, dispersed phase is gas, and granular phase is solid.

Note:
	190526 -- Move all the phase models to "./multiphaseModels/" (Edited by CHEN Bo)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

/*-------------------ADDING CODE------------test-----------*/
    // Include basic kinematic libaray
    #include "basicKinematicCollidingCloud.H"
    #define basicKinematicTypeCloud basicKinematicCollidingCloud
/*-------------------END ADDING CODE-------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    /*-------------------ADDING CODE-----------------------*/   
        // Specify the cloud name
        // E.g.: myTwoPhaseTestFoam -cloudName <NAME> -help
        argList::addOption
        (
            "cloudName",
            "name",
            "specify alternative cloud name. default is 'kinematicCloud'"
        );
    /*-------------------END ADDING CODE-------------------*/

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

    #include "pUf/createRDeltaTf.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        int nEnergyCorrectors
        (
            pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)
        );

        if (LTS)
        {
            #include "setRDeltaT.H"
            if (faceMomentum)
            {
                #include "setRDeltaTf.H"
            }
        }
        else
        {
            #include "CourantNos.H"
            #include "setDeltaT.H"
        }

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            #include "YEqns.H"

            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
                #include "EEqns.H"
                #include "pUf/pEqn.H"
            }
            else
            {
                #include "pU/UEqns.H"
                #include "EEqns.H"
                #include "pU/pEqn.H"
            }

            fluid.correctKinematics();

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
