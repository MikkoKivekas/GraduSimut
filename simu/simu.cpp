
Skip to content
Using University of Jyväskylä Mail with screen readers
MeetNew
Start a meeting
Join a meeting
Chat
1 of 332
MicroAlphatrossin simulaatiot ja gradupalaute
Inbox
	x
Mikko Kivekäs <mikko.m.kivekas@student.jyu.fi>
	
Fri, 27 Mar, 10:22
	
to Taneli
Hei,

oli joskus puhetta niistä MicroAlphatrossin simulaatioista ja etenkin
kuvista, että jakaisit ne gradua varten. Nyt epidemiasulun aikana
ehtisin hyvin tutustumaan niihin

Odottelen vieläkin sinun kommentteja gradun 20%-vaiheeseen. Onko
jotain tiettyjä aiheita ja asioita, joita erityisesti haluat gradussa
nähdä käsiteltävän?


Terveisin
Mikko Kivekäs
Kalvas, Taneli <taneli.v.m.kalvas@jyu.fi>
	
AttachmentsFri, 27 Mar, 11:06
	
to Mikko
Moi!

Sinulla on varmasti jotain uudempaa jo. Laitatko minulle uusimman version gradusta niin luen läpi.

Simulaatiot on tehty IBSimu-softalla (ibsimu.sourceforge.net). Simulaatio koostuu teknisistä syistä kahdesta osasta. Ensimmäistä ajaa simu.cpp ja toista cont.cpp. Kuvat simulaatioista on ohessa. Geometria tulee dxf-tiedostosta, josta tämä versio 9 on se mitä on tällä hetkellä rakennettuna. Magneettikenttä on tiedostossa bfield_2d.dat. Verkossa on jonkun verran opetusmateriaalia noihin simulaatioihin liitttyen. Katso verkkosivut läpi.

Taneli

--
Taneli Kalvas
Ph.D., Senior researcher
Department of Physics, room FL114
P.O. Box 35 (YFL)
40014 University of Jyväskylä, Finland
Mobile: +358-44-314-1602
Fax:    +358-14-617-411
Email:  taneli.kalvas@jyu.fi

________________________________________
From: Mikko Kivekäs <mikko.m.kivekas@student.jyu.fi>
Sent: Friday, March 27, 2020 10:22
To: Kalvas, Taneli
Subject: MicroAlphatrossin simulaatiot ja gradupalaute
...

[Message clipped]  View entire message
6 Attachments
Mikko Kivekäs
	
	AttachmentsFri, 27 Mar, 11:19
Kiitos, perehdyn simulaatioihin ja IBSimuun. Liitteenä nykyinen versio gradusta, siellä on joitain tiedostettuja muutosta vaativia kohtia ja muistiinpanoja seas
Kalvas, Taneli
	
	Fri, 3 Apr, 11:05
Tässä vähän kommentteja. Enemmän alussa ja vähemmän lopussa. Filu dropboxissa. Jostain syystä tuosta tuli valtaisan kokoinen. https://www.dropbox.com/s/euy3axqd
Mikko Kivekäs
	
	Thu, 16 Apr, 14:54
Hei, tässä koitan asennella IBSimua Linux-koneelle, ja IBSimun configurointivaiheessa lyö errorin "required library fontconfig not found". Olen asentanut packag
taneli.v.m.kalvas@jyu.fi
	
	Thu, 16 Apr, 15:57
Moi! Arvaisin ongelman olevan seuraava. Kirjastot on jaettu kahteen osaa. On ns. binääripaketti joka tarvitaan kirjaston käyttämiseen ja sitten on erikseen osa
Mikko Kivekäs
	
	Thu, 16 Apr, 16:14
Kiitos, tällä ratkesi ensimmäinen virheilmoitus. Nyt tuli kuitenkin seuraavaksi ilmoitus, ettei vaadittavaa GSL-versiota ei löydy. Onko tähän ongelmaan valistun
taneli.v.m.kalvas@jyu.fi
	
	Thu, 16 Apr, 16:45
Asensitko gsl:n käsipelillä vai paketinhallinnasta? Jos käsin on arvaukseni se, että profiilistasi puttuu pkg-config -työkalun polusta yksi hakemisto. Myöhemmin
Mikko Kivekäs
	
	Thu, 16 Apr, 17:24
Hei ja kiitos taas avusta, homma pyörähti käyntiin kun asenteli GSLn uudelleen paketinhallinnasta. Kuitenkin make-vaiheessa tulee uusi ongelma vastaan ja antaa
taneli.v.m.kalvas@jyu.fi
	
	Thu, 16 Apr, 17:59
Sinulla on kai Macci? Laita configuren perään optio --disable-sigsegv_stack niin tuon pitäisi mennä pois. Sent: Thursday, April 16, 2020 17:24 Email: taneli.kal
Mikko Kivekäs
	
	Fri, 17 Apr, 11:41
Ei ole Mac, ihan Debian Linux/GNU-distro. Kokeilen miten toimii, jos hylkää virheilmoitukset. Nyt testausvaiheessa en ole kunnolla ymmärtänyt makefilen formaatt
taneli.v.m.kalvas@jyu.fi
	
	Fri, 17 Apr, 16:13
Hylkää virheilmoitukset? Nyt en ymmärrä. Auttoiko --disable-sigsegv_stack. Siis configure ajetaan ensin ja sitten make. IBSimun tapauksessa configure skripti si
Mikko Kivekäs
	
	Fri, 17 Apr, 17:11
Laitoin siis -i make-komentojen perään. Ei auttanut --disable-sigsegv_stack, ei edes tunnistanut tätä komentoa, ja configure on ajettu. Kokeilen noita tutoriaal
taneli.v.m.kalvas@jyu.fi
	
	Fri, 17 Apr, 17:18
Sinun täytyy saada peräkanaa ajettua läpi ./configure --prefix=/home/tvkalvas --disable-sigsegv_stack make make install ilman virheilmoituksia jotta pääset etee
Mikko Kivekäs
	
	Mon, 27 Apr, 16:07
Hei, olen tässä koittanut asennella IBSimua myös toiselle Linux-koneelle ja sain seuraavanlaisen virheilmoituksen libibsimun make-vaiheessa g++ error: unrecogni
Kalvas, Taneli
	
	Wed, 29 Apr, 10:36
Make-vaiheen ongelmaa en ymmärrä. Jälkimmäinen ongelma liittyy siihen, että sinulla ei ole PKG_CONFIG_PATH asetettuna oikein. Kuten aikaisemminkin sanoin on esi
Mikko Kivekäs
	
	Wed, 29 Apr, 10:49
Hei, joo katsotaan vaan huomenna videoyhteydellä kuntoon. 9.30 käy mainiosti.
Kalvas, Taneli
	
	Thu, 30 Apr, 09:32 (13 days ago)
Join Zoom Meeting https://jyufi.zoom.us/j/68849398957 Meeting ID: 688 4939 8957 Password: 743415 One tap mobile +14086528184,,68849398957# US (San Jose) +166990
Mikko Kivekäs
	
	Tue, 5 May, 11:07 (8 days ago)
Hei, olen tässä nyt käynyt IBSimun tutoriaalit läpi ja kävin läpi kesäkoulumateriaalit simulaatioineen. Osa kesäkoulun simulaatioista ei ottanut toimiakseen, en
Kalvas, Taneli
	
	Tue, 5 May, 11:16 (8 days ago)
Jaahas. Kakkosesimerkki sanoo meille siis, että . Layers in dxf-file: 0: '0' 1: 'plasma' 2: 'puller' Koodista taas näemme, että DXF-tiedostosta etsitään layerei
Mikko Kivekäs
	
	Tue, 12 May, 12:55 (20 hours ago)
Hei, kiitos avusta ja nyt kaikki toimii kuten pitäisi. Olen tässä nyt puuhaillut jotain tuon MicroAlphatrossin simulaation kanssa. Millaisia ajatuksia sinulla o
Kalvas, Taneli
	
07:55 (1 hour ago)
	
to Mikko
Moi!

Ajatukset eivät ole kovin kirkkaat siinä mielessä, että en osaa etukäteen sanoa mikä on oleellista ja mielekästä. Tutkimalla selviää. Lähtökohtaisesti kuitenkin sovelluksen kannalta on oleellista, että saamme läpi riittävän määrän ionivirtaa. Ionivirtaan liittyy sellainen seikka, että emme vielä tiedä paljonko ionilähde tuottaa. Jollain ionivirran määrällä meillä on varmasti melko hyvä läpäisy ekstraktiosysteemin läpi. Kun ionikähteestä tulevan virran määrä kasvaa saattaa jossain kohdassa läpäisy heikentyä. Saman tapaista saattaa tapahtua myös suihkun laadulle. Nämä ovat mielenkiintoisia seikkoja. Suihkun laatua kannattanee tarkastella varauksenvaihtokammioon mennessä. Todellisuudessa varauksenvaihtokammiossa tapahtuu myös sirontaa joka tulee heikentämään suihkun laatua. Tämän vuoksi varauksenvaihdon jälkeen simulaatio antaa vain alarajan emittanssille. Simulaatio jonka annoin sinulle on jossakin määrin optimoitu, mutta sitä on myös tietoisesti "huononnettu", sillä tekniset ratkaisut ovat estäneet ionioptiikan kannalta optimaaliset ratkaistut. Nämä ovat siis erityisesti geometriaan liittyviä ratkaisuja.

Tyypillisesti ongelmallisia paikkoja ekstraktiojärjestelmissä ovat paikat, joissa suihku kulkee läheltä elektrodeja. Näissä paikoissa fokusvoima on usein epälineaarinen, joka aiheuttaa suihkun laadun heikentymistä. Kuitenkin tässä tapauksessa suihku on tietoisesti haluttu kulkemaan ahdasta kanavaa pitkin, koska rubidiumin päätymistä ionilähteeseen on haluttu estää.

Tämä fokus on sellainen seikka jota voisi kenties puida vähän teoreettisestikin ja/tai tutkimalla jotain yksittäistä linssiä simulaatiossa. Periaatteessahan ionioptiikka toimii kuten valo-optiikka. Linssille voidaan määrittää fokaalipituus, joka epälineaarisuuden vuoksi on siis säteen funktio. Eli jos hiukkanen kulkee lähellä optista akselia se taittuu fokaalipituuden f mukaan ja jos hiukkanen on kauempana akselista on fokaalipituus jotain muuta. Voisit kenties rakentaa yksinkertaisen simulaation, jossa on vain einzel-linssi, jonka fokuskäyttäytymistä tutkitaan.

Emittanssiarvo jonka IBSimu antaa ulos on puhdas rms-arvo. Geometrisesti tulkittuna ellipsiksi se on pinta-alaltaan puoliakselien tulo eikä siis ellipsin koko pinta-ala. Tästä syystä joskus pii kirjoitetaan yksikköön. Mikäli kontekstista tai merkinnästä on selvää että puhutaan rms-arvosta ei piin merkitsemiselle ole mitään tarvetta minun mielestäni. Toki kaikki eivät ole tätä mieltä. Emittanssi on ongelmallinen koska konventioita on niin monia. Siksi aina emittanssista puhuttaessa on lukijalle tehtävä selväksi mitä juuri sinä kirjoittajana sillä tarkoitat, jotta tieteellisen tekstin täsmällisyys toteutuu.

Taneli

--
Taneli Kalvas
Ph.D., Senior researcher
Department of Physics, room FL114
P.O. Box 35 (YFL)
40014 University of Jyväskylä, Finland
Mobile: +358-44-314-1602
Fax:    +358-14-617-411
Email:  taneli.kalvas@jyu.fi

________________________________________
From: Mikko Kivekäs <mikko.m.kivekas@student.jyu.fi>
Sent: Tuesday, May 12, 2020 12:55
To: Kalvas, Taneli
Subject: Re: MicroAlphatrossin simulaatiot ja gradupalaute

Hei,

kiitos avusta ja nyt kaikki toimii kuten pitäisi.

Olen tässä nyt puuhaillut jotain tuon MicroAlphatrossin simulaation kanssa. Millaisia ajatuksia sinulla oli gradun simulaatioanalyysin suhteen? Mitkä on niitä tärkeitä suihkuparametreja tarkastella, emittanssi ja ionivirta ainakin? Mistä kohtaa emittanssi tulee määrittää, oletan että juuri ennen varauksenvaihtoa vai vasta sen jälkeen? Onko tähän simulaatioon löydetty jo jotkin optimiparametrit?

Lisäkysymys tuosta simulaation antamasta emittanssin yksiköstä: onko pii leivottu numeroarvon sisään vai miten tuo toimii?

Mikko Kivekäs

Email:  taneli.kalvas@jyu.fi<mailto:taneli.kalvas@jyu.fi>

________________________________________
From: Mikko Kivekäs <mikko.m.kivekas@student.jyu.fi<mailto:mikko.m.kivekas@student.jyu.fi>>
Sent: Tuesday, May 5, 2020 11:07
To: Kalvas, Taneli
Subject: Re: MicroAlphatrossin simulaatiot ja gradupalaute

Hei,

olen tässä nyt käynyt IBSimun tutoriaalit läpi ja kävin läpi kesäkoulumateriaalit simulaatioineen. Osa kesäkoulun simulaatioista ei ottanut toimiakseen, eniten ongelmia on dxf-tiedostojen kanssa. Esimerkiksi tällainen error-viesti kesäkoulun toisesta esimerkistä plasma.cpp
   Error in dxf_solid.cpp:73 in DXFSolid(): No entities in layer

Ja einzel3d esimerkissä make antaa virheenä
   Makefile:20 *** puuttuva erotin. Seis.


Mikko Kivekäs



On Thu, 30 Apr 2020 at 09:32, Kalvas, Taneli <taneli.v.m.kalvas@jyu.fi<mailto:taneli.v.m.kalvas@jyu.fi><mailto:taneli.v.m.kalvas@jyu.fi<mailto:taneli.v.m.kalvas@jyu.fi>>> wrote:


--
Taneli Kalvas
Ph.D., Senior researcher
Department of Physics, room FL114
P.O. Box 35 (YFL)
40014 University of Jyväskylä, Finland
Mobile: +358-44-314-1602
Fax:    +358-14-617-411
Email:  taneli.kalvas@jyu.fi<mailto:taneli.kalvas@jyu.fi><mailto:taneli.kalvas@jyu.fi<mailto:taneli.kalvas@jyu.fi>>
Join Zoom Meeting
https://jyufi.zoom.us/j/68849398957

Meeting ID: 688 4939 8957
Password: 743415
One tap mobile
+14086528184,,68849398957# US (San Jose)
+16699006833,,68849398957# US (San Jose)

Dial by your location
        +1 408 652 8184 US (San Jose)
        +1 669 900 6833 US (San Jose)
        +1 929 205 6099 US (New York)
        +1 253 215 8782 US (Tacoma)
        +1 301 715 8592 US (Germantown)
        +1 312 626 6799 US (Chicago)
        +1 346 248 7799 US (Houston)
Meeting ID: 688 4939 8957
Find your local number: https://jyufi.zoom.us/u/cbQCL6zyrK

Join by SIP
68849398957@109.105.112.236<mailto:68849398957@109.105.112.236><mailto:68849398957@109.105.112.236<mailto:68849398957@109.105.112.236>>
68849398957@109.105.112.235<mailto:68849398957@109.105.112.235><mailto:68849398957@109.105.112.235<mailto:68849398957@109.105.112.235>>

Join by H.323
109.105.112.236
109.105.112.235
Meeting ID: 688 4939 8957
Password: 743415



https://jyufi.zoom.us/j/68849398957?pwd=RGw1azVxazFpNHVMclloZ3FTYmJxQT09

--
Taneli Kalvas
Ph.D., Senior researcher
Department of Physics, room FL114
P.O. Box 35 (YFL)
40014 University of Jyväskylä, Finland
Mobile: +358-44-314-1602
Fax:    +358-14-617-411
Email:  taneli.kalvas@jyu.fi<mailto:taneli.kalvas@jyu.fi><mailto:taneli.kalvas@jyu.fi<mailto:taneli.kalvas@jyu.fi>>

________________________________________
From: Mikko Kivekäs <mikko.m.kivekas@student.jyu.fi<mailto:mikko.m.kivekas@student.jyu.fi><mailto:mikko.m.kivekas@student.jyu.fi<mailto:mikko.m.kivekas@student.jyu.fi>>>
Sent: Wednesday, April 29, 2020 10:49
To: Kalvas, Taneli
Subject: Re: MicroAlphatrossin simulaatiot ja gradupalaute

Hei,

joo katsotaan vaan huomenna videoyhteydellä kuntoon. 9.30 käy mainiosti.

Mikko Kivekäs

On Wed, 29 Apr 2020 at 10:36, Kalvas, Taneli <taneli.v.m.kalvas@jyu.fi<mailto:taneli.v.m.kalvas@jyu.fi><mailto:taneli.v.m.kalvas@jyu.fi<mailto:taneli.v.m.kalvas@jyu.fi>><mailto:taneli.v.m.kalvas@jyu.fi<mailto:taneli.v.m.kalvas@jyu.fi><mailto:taneli.v.m.kalvas@jyu.fi<mailto:taneli.v.m.kalvas@jyu.fi>>>> wrote:
Make-vaiheen ongelmaa en ymmärrä. Jälkimmäinen ongelma liittyy siihen, että sinulla ei ole PKG_CONFIG_PATH asetettuna oikein. Kuten aikaisemminkin sanoin on esim ~/.bashrc -tiedostossa oltava määritettynä ylimmäräiset kirjastohakemistot tähän tapaan:

export PKG_CONFIG_PATH="/home/tvkalvas/lib/pkgconfig:$PKG_CONFIG_PATH"

Jos tästä ei tule mitään näin niin hoidetaan IBSimu kuntoon videoyhteyden kanssa. Esim. huominen on kalenterissa tyhjää. Jos aloitetaan 9:30. Käykö?

Taneli

--
Taneli Kalvas
Ph.D., Senior researcher
Department of Physics, room FL114
P.O. Box 35 (YFL)
40014 University of Jyväskylä, Finland
Mobile: +358-44-314-1602
Fax:    +358-14-617-411
Email:  taneli.kalvas@jyu.fi<mailto:taneli.kalvas@jyu.fi><mailto:taneli.kalvas@jyu.fi<mailto:taneli.kalvas@jyu.fi>><mailto:taneli.kalvas@jyu.fi<mailto:taneli.kalvas@jyu.fi><mailto:taneli.kalvas@jyu.fi<mailto:taneli.kalvas@jyu.fi>>>

________________________________________
From: Mikko Kivekäs <mikko.m.kivekas@student.jyu.fi<mailto:mikko.m.kivekas@student.jyu.fi><mailto:mikko.m.kivekas@student.jyu.fi<mailto:mikko.m.kivekas@student.jyu.fi>><mailto:mikko.m.kivekas@student.jyu.fi<mailto:mikko.m.kivekas@student.jyu.fi><mailto:mikko.m.kivekas@student.jyu.fi<mailto:mikko.m.kivekas@student.jyu.fi>>>>
...

[Message clipped]  View entire message
	
	
	

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
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"



using namespace std;


const int nrounds = 15;
const double r0 = 0.75e-3;
const double rplasma = 2.0e-3;

const double h = 4e-5;
const double Nperh = 1000.0;
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
const double Veinzel = -1.375e3;
const double Vconv = -7e3;
const double Vgnd = -15e3;
const double Veinzel2 = -20e3;
const double I = 1e-3;
const double J = 1.35*I/(M_PI*r0*r0);

const double sc_alpha = 0.5;
string stamp = "_10";

/* Veinzel    Ibound2(cont)
 * -1.0e3     0.000147271 A
 * -1.125e3   0.00019536 A
 * -1.25e3    0.000219463 A
 * -1.375e3   0.000204841 A,    0.000203443 A without bfield
 * -1.5e3     0.000172724 A
 * -1.75      0.000123436 A
 * -2.0e3     9.28533e-05 A
 *
 * Total current 0.00103587 A
 */


void simu( int argc, char **argv )
{
    double sizereq[3] = { 71.0e-3,
                          25.0e-3, 
                           0.0e-3 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
                    (int)floor(sizereq[2]/h)+1 );
    Vec3D origo( -1e-3, 0, 0 );
    Geometry geom( MODE_CYL, meshsize, origo, h );

    MyDXFFile *dxffile = new MyDXFFile( "muokattu10.dxf" );
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

    geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  2,  Bound(BOUND_DIRICHLET, Vconv) );
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
    InitialPlasma initp( AXIS_X, 0.2e-3 );
    solver.set_initial_plasma( Up, &initp );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );

    // Define magnetic field
    bool fout[3] = {true, true, false};
    MeshVectorField bfield( MODE_CYL, fout, 1.0e-3, 1.0, "../bfield_2d.dat" );
    field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO };
    bfield.set_extrapolation( bfldextrpl );
    bfield.translate( Vec3D(-4e-3,0,0) );

    //MeshVectorField bfield(geom,fout);
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

    PPlasmaBfieldSuppression psup( epot, 20.0 );
    pdb.set_bfield_suppression( &psup );

    ibsimu.message(1) << "J = " << J << " A/m2\n";

    double rho_tot;
    for( size_t i = 0; i < nrounds; i++ ) {
	
	ibsimu.message(1) << "Iteration round " << i << "\n";

	if( i == 1 ) {
            solver.set_pexp_plasma( rho_tot, Te, Up );
        }
	
	solver.solve( epot, scharge_ave );
	//int iterc = solver.get_iter();
        //if( iterc == 0 ) {
	//ibsimu.message(1) << "Zero iterations, breaking cycle\n";
	//break;
        //}
	efield.recalculate();

        pdb.clear(); 
	ibsimu.message(1) << "J = " << J << " A/m2\n";
	pdb.add_2d_beam_with_energy( Npart, J, q, m, E0, Tp, Tt, 
				     origo[0], 0.0, 
				     origo[0], rplasma );
        pdb.iterate_trajectories( scharge, efield, bfield );
	rho_tot = pdb.get_rhosum();

	if( i == 0 ) {
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

	if( i == nrounds-1 ) {
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

    // Write output file containing all particles
    ofstream fileOut( "particles_out.txt" );
    for( size_t k = 0; k < pdb.size(); k++ ) {

	ParticleCyl &pp = pdb.particle( k );
	
	// Skip ions not at the end
	if( pp(PARTICLE_X) < geom.max(0)-geom.h() )
	    continue;
	
	fileOut << setw(12) << pp.IQ() << " ";
	// t, x, vx, r, vr, w
	for( size_t j = 0; j < 6; j ++ )
	    fileOut << setw(12) << pp(j) << " ";
	fileOut << "\n";
    }
    fileOut.close();
    
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
    geomplotter.plot_png( "particle_plot" + stamp + ".png" );

    if( false ) {
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

simu.cpp
Displaying simu.cpp.
