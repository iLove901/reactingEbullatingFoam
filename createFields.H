#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<twoPhaseSystem> fluidPtr
(
    twoPhaseSystem::New(mesh)
);
twoPhaseSystem& fluid = fluidPtr();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volScalarField& p = fluid.phase1().thermoRef().p();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());


/*-------------------ADDING CODE-----------------------*/  

    // Look for continuous phase name in "transportProperties"  
    Info << ">>> AxiMeta: Reading name for continuous phase" << endl;
    word continuousPhaseName
    (
        IOdictionary
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ
            )
        ).lookup("continuousPhaseName")
    );

    // Look for velocity field of continuous phase
    Info << ">>> AxiMeta: Reading velocity field of continuous phase" << endl;
    volVectorField Uc
    (
        IOobject
        (
            IOobject::groupName("U", continuousPhaseName),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Reading face flux from continuous phase
    Info << ">>> AxiMeta: Reading/calculating continuous phase face flux field" << endl;
    surfaceScalarField phic
    (
        IOobject
        (
            IOobject::groupName("phi", continuousPhaseName),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Uc) & mesh.Sf()
    );

    // Reading and writing mean rho value 
    // Here rhoc represents the rho for fluid phase, i.e., gas and liquid
    Info << ">>> AxiMeta: Reading and writing density (rho) fro continuous phase" << endl;
    volScalarField rhoc
    (
        IOobject
        (
            "rhom",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
 
    // Calculating dynamic viscosity for continuous phase
    volScalarField muc
    (
        IOobject
        (
            IOobject::groupName("mu", continuousPhaseName),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.phase2().thermo().rho()*fluid.phase2().nu()
    );

/*
Note: 190526-Revised to this part

    // Calculate alpha for fluid phase
    // alphac must be constructed before the cloud
    // so that the drag-models can find it
    Info << ">>> AxiMeta: Creating alpha field for fluid phase" << endl;
    volScalarField alphac
    (
        IOobject
        (
            "alpha.fluid",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    );

    // Get cloud name
    word kinematicCloudName("kinematicCloud");
    args.optionReadIfPresent("cloudName", kinematicCloudName);

    // Constructing kinematic cloud
    Info << ">>> AxiMeta: Constructing kinematic cloud \"" << kinematicCloudName << "\"" << endl;
    basicKinematicTypeCloud kinematicCloud
    (
        kinematicCloudName,
        rhoc,
        Uc,
        muc,
        g
    );

    // Particle fraction upper limit
    scalar alphacMin
    (
        1.0 
        - readScalar
        (
            kinematicCloud.particleProperties().subDict("constantProperties")
        .lookup("alphaMax")
        )
    );

    // Update alphac from the particle locations
    alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
    alphac.correctBoundaryConditions();

    // Correct alphac flux
    surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));
    surfaceScalarField alphaPhic
    (
        IOobject::groupName
        (
            "alphaPhi",
            continuousPhaseName
        ),
        alphacf*phic
    );

*/
/*-------------------END ADDING CODE-------------------*/


