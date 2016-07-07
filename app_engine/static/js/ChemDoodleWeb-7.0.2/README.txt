ChemDoodle Web Components

=====
ABOUT
=====

The ChemDoodle Web Components library is a pure Javascript chemical graphics and cheminformatics library derived from the ChemDoodle(R) application and produced by iChemLabs. ChemDoodle Web Components allow the wielder to present publication quality 2D and 3D graphics and animations for chemical structures, reactions and spectra. Beyond graphics, this tool provides a framework for user interaction to create dynamic applications through web browsers, desktop platforms and mobile devices such as the iPhone, iPad and Android devices. This library also has complete access to the entire ChemDoodle desktop API through AJAX, allowing for quick access to one of the most robust chemical graphics and informatics packages in existence directly through Javascript.

=======
LICENSE
=======

The ChemDoodle Web Components library is licensed under version 3 of the GNU GENERAL PUBLIC LICENSE.

The complete license is provided in the COPYING.txt file of this download.

Also note that the UI icons provided (in source code as data URIs) are iChemLabs Proprietary and may only be used in projects using ChemDoodle Web Components.

Contact us for alternate licensing options: http://www.ichemlabs.com/contact-us

=========
STRUCTURE
=========

Downloaded zip file contents:

- README.txt                ==> This readme file
- COPYING.txt               ==> GNU GPL Version 3 license

- data                      ==> folder containing chemical data
  - molecules               ==> folder containing structural data
  - spectra                 ==> folder containing spectral data

- install                   ==> folder containing the resources to use ChemDoodle Web Components
  - ChemDoodleWeb.css       ==> Basic CSS file for the library
  - ChemDoodleWeb.js        ==> ChemDoodle Web Components library Javascript (packed)
  - uis                     ==> folder containing the SketcherCanvas and EditorCanvas3D plugin
    - images                        ==> folder containing the corresponding images for the jQuery UI CSS
    - jquery-ui-1.9.2.custom.css    ==> jQuery UI CSS
    - ChemDoodleWeb-uis.js          ==> UIs plugin Javascript (packed)

- samples                   ==> folder containing full sample HTML pages of example components, open them in your browser to view them, view the source to see how they work. View the tutorial on the website for in depth instructions.
  - periodicTable.html                         ==> A PeriodicTableCanvas component.
  - rotator_2D_ballnstick.html                 ==> A RotatorCanvas with a spinning molecule in a 2D-based "Ball and Stick" representation.
  - sketcher_2D_full.html                      ==> A SketcherCanvas component initialized as a Full Sketcher.
  - sketcher_2D_single_molecule.html           ==> A SketcherCanvas component initialized as a Single Molecule Sketcher.
  - spectrum_mass_styled.html                  ==> A PerspectiveCanvas component that has been styled and presents a mass spectrum.
  - transformer_morphine_2D_wireframe.html     ==> A TransformCanvas component. Interact with a 3D model of morphine in a 2D-based "Wireframe" representation.
  - viewer_caffeine_2D_acs1996.html            ==> A ViewerCanvas component with caffeine drawn to ACS Document 1996 style.
  - webgl_transformer_1CRN_3D_Ribbon.html      ==> A TransformCanvas3D component showing an interactive cartoon ribbon model of the PDB ID 1CRN. WebGL is required.
  - webgl_transformer_2OEY_3D_NucleicAcid.html ==> A TransformCanvas3D component showing an interactive model of the Nucleic Acid represented by PDB code 2OEY. WebGL is required.
  - webgl_transformer_DDT_3D_ballnstick.html   ==> A TransformCanvas3D component showing an interactive model of DDT in a "Ball and Stick" representation. WebGL is required.
  - webgl_transformer_MAZ_3D_Crystal.html      ==> A TransformCanvas3D component showing an interactive model of the zeolite MAZ. WebGL is required.

- src                                   ==> folder containing the unpacked Javascript source
  - ChemDoodleWeb-unpacked.js           ==> ChemDoodle Web Components library Javascript (unpacked)
  - ChemDoodleWeb-uis-unpacked.js       ==> SketcherCanvas and EditorCanvas3D plugin library Javascript (unpacked)
  
- tests                       ==> folder containing the qUnit testing suite
  - ChemDoodleWeb-tests.js    ==> ChemDoodle Web Components testing suite, requires qUnit (http://docs.jquery.com/Qunit)

============
INSTRUCTIONS
============

For installation instructions and to obtain the latest version of the ChemDoodle Web Components library, please visit the official ChemDoodle Web Components website: http://web.chemdoodle.com

===========
ATTRIBUTION
===========

If you use ChemDoodle Web Components in your website or product, we request that you provide a link on your site to web.ChemDoodle.com and/or www.ChemDoodle.com. When using the SketcherCanvas or EditorCanvas3D components, you are not allowed to remove any of our links or logos. When requesting support, you must provide the URL to your ChemDoodle Web Components product. Our logos and links must be clearly visible on this page when requesting support. Please give credit where credit is due.

==================
ICHEMLABS SERVICES
==================

ChemDoodle Web Components has built in access to the entire ChemDoodle desktop API through AJAX XMLHttpRequest Level 2 (XHR2) calls to our server. XHR2 requires that our servers recognize your origins (domains). If you are an academic or non-profit organization, we will provide free access to our server resources. Just contact us to register your origins with our server. Otherwise, our rates are very reasonable, please contact us for options.

http://www.ichemlabs.com/contact-us

=========
QUESTIONS
=========

If you have any questions about ChemDoodle Web Components, please contact us at http://www.ichemlabs.com/contact-us .
