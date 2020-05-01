#include <fstream>
#include <iomanip>
#include <limits>
#include "epot_bicgstabsolver.hpp"
//#include "epot_umfpacksolver.hpp"
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "readascii.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"



using namespace std;


const double r0 = 0.75e-3;
const double rplasma = 2.0e-3;

//const double h = 5e-5;
const double h = 1e-4;    // For full geometry
const double Nperh = 500.0;
const uint32_t Npart = Nperh*rplasma/h;

const double q = 1.0;
const double m = 4.0;
const double E0 = 4.0;

const double Tp = 0.0;
const double Tt = 0.5;

const double Te = 5.0;
const double Up = 5.0;
const double Vplasma = 0;
const double Vpuller = -7e3;
const double Veinzel = -1.5e3;
const double Vconv = -7e3;
const double Vgnd = -15e3;
const double Veinzel2 = -20e3;
const double I = 1e-3;
const double J = 1.35*I/(M_PI*r0*r0);

const double sc_alpha = 0.5;
string stamp = "_9";

void simu( int argc, char **argv )
{
    // Input particles
    ReadAscii din( "particles_out.txt", 7 );
    cout << "Reading " << din.rows() << " particles\n";
    
    double xstart = din[2][0];

    double sizereq[3] = { 330.0e-3-xstart,
                          46.0e-3, 
                           0.0e-3 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
                    (int)floor(sizereq[2]/h)+1 );
    Vec3D origo( xstart, 0, 0 );
    Geometry geom( MODE_CYL, meshsize, origo, h );

    MyDXFFile *dxffile = new MyDXFFile( "muokattu9.dxf" );
    dxffile->set_warning_level( 2 );
    MyDXFEntities *e = dxffile->get_entities();
    MyDXFEntitySelection *sel = e->selection_all();
    e->scale( sel, dxffile, 1.0e-3 );

    DXFSolid *s1 = new DXFSolid( dxffile, "plasma" );
    geom.set_solid(  7, s1 );
    DXFSolid *s2 = new DXFSolid( dxffile, "puller" );
    geom.set_solid(  8, s2 );
    DXFSolid *s3 = new DXFSolid( dxffile, "einzel" );
    geom.set_solid(  9, s3 );
    DXFSolid *s4 = new DXFSolid( dxffile, "conv" );
    geom.set_solid( 10, s4 );
    DXFSolid *s5 = new DXFSolid( dxffile, "gnd" );
    geom.set_solid( 11, s5 );
    //DXFSolid *s6 = new DXFSolid( dxffile, "einzel2" );
    //geom.set_solid( 12, s6 );
    
    geom.set_boundary(  1,  Bound(BOUND_DIRICHLET, Vconv) );
    geom.set_boundary(  2,  Bound(BOUND_DIRICHLET,  Vgnd) );
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.0) );

    geom.set_boundary(  7,  Bound(BOUND_DIRICHLET, Vplasma) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET, Vpuller) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET, Veinzel) );
    geom.set_boundary( 10,  Bound(BOUND_DIRICHLET, Vconv) );
    geom.set_boundary( 11,  Bound(BOUND_DIRICHLET, Vgnd) );
    //geom.set_boundary( 12,  Bound(BOUND_DIRICHLET, Veinzel2) );
    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    //EpotUMFPACKSolver solver( geom );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );

    // Define magnetic field
    bool fout[3] = {true, true, false};
    /*
    MeshVectorField bfield( MODE_CYL, fout, 1.0e-3, 1.0, "../bfield_2d.dat" );
    field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO };
    bfield.set_extrapolation( bfldextrpl );
    bfield.translate( Vec3D(-4e-3,0,0) );
    */
    MeshVectorField bfield(geom,fout);

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBaseCyl pdb( geom );
    pdb.set_max_steps( 1000 );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    for( size_t iter = 0; iter < 30; iter++ ) {
	
	solver.solve( epot, scharge_ave );
	int iterc = solver.get_iter();
        if( iterc == 0 ) {
	    ibsimu.message(1) << "Zero iterations, breaking cycle\n";
	    break;
        }
	efield.recalculate();

        pdb.clear(); 
	// Go through all read particles
	// I, t, x, vx, r, vr, w
	for( size_t j = 0; j < din.rows(); j++ ) {
	    double I  = din[0][j];
	    double t  = din[1][j];
	    double x  = din[2][j];
	    double vx = din[3][j];
	    double r  = din[4][j];
	    double vr = din[5][j];
	    double w  = din[6][j];
	    
	    pdb.add_particle( I, 1.0, 4.0, ParticlePCyl(t,x,vx,r,vr,w) );
	}
        pdb.iterate_trajectories( scharge, efield, bfield );

	// Conversion almost completely removes space charge effects
	for( uint32_t j = 0; j < scharge.size(1); j++ ) {
	    for( uint32_t i = 0; i < scharge.size(0); i++ ) {
		double x = scharge.origo(0) + i*scharge.h();
		    if( x > 0.105 )
			scharge(i,j) *= 0.01;
	    }	
	}

	if( iter == 0 ) {
	    scharge_ave = scharge;
	} else {
	    double sc_beta = 1.0-sc_alpha;
            uint32_t nodecount = scharge.nodecount();
            for( uint32_t b = 0; b < nodecount; b++ ) {
                scharge_ave(b) = sc_alpha*scharge(b) + sc_beta*scharge_ave(b);
            }
	}

	// Trajectory diagnostics
        TrajectoryDiagnosticData tdata;
        std::vector<trajectory_diagnostic_e> diagnostics;
        diagnostics.push_back( DIAG_R );
        diagnostics.push_back( DIAG_RP );
        pdb.trajectories_at_plane( tdata, AXIS_X, geom.max(0)-geom.h(), diagnostics );
        Emittance emit( tdata(0).data(), tdata(1).data() );     

        // Output
        ofstream dout( "emittance.txt", ios_base::app );
        dout << emit.alpha() << " "
             << emit.beta() << " "
             << emit.epsilon() << "\n";
        dout.close();

	if( false ) {
	    MeshScalarField tdens( geom );
	    pdb.build_trajectory_density_field( tdens );

	    GTKPlotter plotter( &argc, &argv );
	    plotter.set_geometry( &geom );
	    plotter.set_epot( &epot );
	    plotter.set_bfield( &bfield );
	    plotter.set_scharge( &scharge );
	    plotter.set_trajdens( &tdens );
	    plotter.set_particledatabase( &pdb );
	    plotter.new_geometry_plot_window();
	    plotter.run();
	}
    }

    geom.save( "geom.dat" );
    epot.save( "epot.dat" );
    pdb.save( "pdb.dat" );

    GeomPlotter geomplotter( geom );
    geomplotter.set_size( 1500, 1500 );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    vector<double> eqpotlines;
    eqpotlines.push_back( -4.0 );
    eqpotlines.push_back( -2.0 );
    eqpotlines.push_back( -1.0 );
    eqpotlines.push_back(  0.0 );
    eqpotlines.push_back(  1.0 );
    eqpotlines.push_back(  2.0 );
    eqpotlines.push_back(  4.0 );
    geomplotter.set_eqlines_manual( eqpotlines );
    geomplotter.set_scharge( &scharge_ave );
    geomplotter.set_size( 1500, 300 );
    geomplotter.plot_png( "cont_particle_plot" + stamp + ".png" );

    if( true ) {
	MeshScalarField tdens( geom );
	pdb.build_trajectory_density_field( tdens );
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_bfield( &bfield );
	plotter.set_efield( &efield );
	plotter.set_scharge( &scharge );
	plotter.set_trajdens( &tdens );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
}


int main( int argc, char **argv )
{
    remove( "emittance.txt" );

    try {
	//ibsimu.set_message_output( "ibsimu" + stamp + ".txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
