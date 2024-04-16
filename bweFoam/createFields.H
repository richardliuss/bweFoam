Info<< "Reading physicalProperties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Am
(
    "Am",
    dimViscosity,
    physicalProperties.lookup("Am")
);
dimensionedScalar rho
(
    "rho",
    dimDensity,
    physicalProperties.lookup("rho")
);

const dimensionedScalar Beta
(
    "Beta",
    dimensionSet(0, 0, 0, 0, 0),
    physicalProperties.lookup("Beta")
);


Info<< "Creating field nc (Manning coefficient)\n" << endl;
const volScalarField nc
(
    IOobject
    (
        "nc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedScalar dim_s (dimTime, 1);

dimensionedScalar dim_m (dimLength, 1);

Info<< "\nReading gravitationalProperties" << endl;

IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g(gravitationalProperties.lookup("g"));
const Switch rotating(gravitationalProperties.lookup("rotating"));
const dimensionedVector Omega =
    rotating ? gravitationalProperties.lookup("Omega")
             : dimensionedVector("Omega", -dimTime, vector(0,0,0));
const dimensionedScalar magg = mag(g);
const dimensionedVector gHat = g/magg;
const scalar noDimMagG = magg.value();

Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field bath\n" << endl;
const volScalarField bath
(
    IOobject
    (
        "bath",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field hU\n" << endl;
volVectorField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Creating field hTotal for post processing\n" << endl;
volScalarField hTotal
(
    IOobject
    (
        "hTotal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h+bath
);
hTotal.write();

dimensionedScalar hDry
(
    "hDry",
    dimLength,
    physicalProperties.lookup("hDry")
);
dimensionedScalar hDry2
(
    "hDry2",
    dimLength,
    physicalProperties.lookup("hDry2")
);

volScalarField H1
(
    IOobject
    (
        "H1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
(h-hDry)*pos(h-hDry) + hDry        
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
hU/H1
);


Info<< "Reading/calculating velocity face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);


Info<< "Creating Coriolis Force" << endl;
const dimensionedVector F("F", ((2.0*Omega) & gHat)*gHat);

mesh.setFluxRequired(h.name());




Info<< "Creating field wetDry\n" << endl;
volScalarField wetDry
(
    IOobject
    (
        "wetDry",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pos(h-hDry)
);

Info<< "Creating field fircCoeff\n" << endl;
volScalarField fircCoeff
(
    IOobject
    (
        "fircCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (noDimMagG)*mag(U)*pow(nc,2.0)/pow(H1/dim_m,1.0/3.0)/dim_m
);
