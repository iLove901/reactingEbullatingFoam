        /*-------------------ADDING CODE-----------------------*/
            Info << ">>> AxiMeta: ---------------------CLOUD COLLIDSION SECTION--------------------------------" << endl;
            // Store the particle positions
            kinematicCloud.storeGlobalPositions();

            // fluid.solve();
            // fluid.correct();
            
            // Revise muc
            muc = rhoc*phase2.nu();
            
            // Evolving the cloud
            Info << ">>> AxiMeta: Evolving the kinematic cloud \"" << kinematicCloud.name() << "\"" << endl;
            kinematicCloud.evolve();

            // Update fluid phase volume fraction field
            Info << ">>> AxiMeta: Update fluid phase volume fraction field" << endl;
            alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
            alphac.correctBoundaryConditions();
            alphacf = fvc::interpolate(alphac);
            alphaPhic = alphacf*phic;

            // Correct the fraction for the two fluid phases
            Info << ">>> AxiMeta: Correct volume fractions for gas and liquid phase" << endl;
            // Dimension difference
            // alphac [0 0 0 0 0 0 0]; phase1.alpha() [1 -1 -1 0 0 0 0]
            /* alpha2 = alphac - alpha1;  */


            // Calculate particle and liquid phase interactions
            Info << ">>> AxiMeta: Calculate phase interaction (solid and liquid)" << endl;
            fvVectorMatrix cloudSU(kinematicCloud.SU(Uc));
            volVectorField cloudVolSUSu
            (
                IOobject
                (
                    "cloudVolSUSu",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector
                (
                    "0",
                    cloudSU.dimensions()/dimVolume,
                    Zero
                ),
                zeroGradientFvPatchVectorField::typeName
            );    
            // Calculate SU
            cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
            cloudVolSUSu.correctBoundaryConditions();
            cloudSU.source() = Zero;  
		
            Info << ">>> AxiMeta: -----------------------END CLOUD COLLIDSION-------------------------------" << endl;                  

        /*-------------------END ADDING CODE-------------------*/

