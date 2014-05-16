;+
;  NAME: XRDUA
;
;
;  PURPOSE: Analyze powder diffraction data from area detectors
;
;
;  STAND ALONE PROGRAMS:
;       XRDUA(main)     - Analyze and convert 2D XRD powder data manual
;       XRD_BP          - Analyze 2D/1D XRD powder data in batch
;       XRD_CHI         - Analyze 1D XRD powder data manual
;       XRD_XDI         - Edit and process results from XRD_BP
;       ExpSetup        - Help with experimental setup
;       EditOverlap     - Edit area overlap statistics
;       oplotchi        - Overplot chi files
;
;
;  MODIFICATION HISTORY: cfr "history"
;
;
;  BUGS:
;    Error when magnification for dyn view is too high (check routine not sufficient?)
;
;
;  INSTALLATION ISSUES:
;    None so far
;
;  TODO:
;
;    TODO: check Lebail fitting
;    TODO: report memory issues (azimuthal integration)
;    TODO: Constant term when fitting with strip+subtract
;    TODO: Convert spacegroup settings (current version not sufficient?)
;    TODO: detection limit (absolute, relative, information depth)
;    TODO: still getting stuck at peak constraints (carpaint: refined+/-C scatter peak positions)
;    TODO: 2D or 1D: load PDFs and click to identify peak
;    TODO: rigid body
;    TODO: PD --(indexing)--> Pawley --(direct methods)--> Rietveld
;    TODO: check spatial distortion again for using b-splines
;    TODO: Coordination polyhedra + fix some of them during fitting
;    TODO: resize windows
;    TODO: crystal size distribution (histogram of 2D image)
;    TODO: sometimes error with green ROIs
;    TODO: Rietveld -> don't recalculate form factors!
;    TODO: batch az.int. save in 1 tiff
;
;
;  XRDUA CORE:
;
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: idlliboverride
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: redefine IDL library functions
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   function EXTRAC                 - function EXTRAC
;   pro ARROW                       - pro ARROW
;   function Identity               - function Identity
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: updatePlatform
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Update version number
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   function updatePlatform         - Update "version" function
;     
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: xrdua
;
;
;  CALLING SEQUENCE: XRDUA
;
;
;  PURPOSE: Analyze 2D XRD powder data manual
;
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   function PointInsidePolygon     - Is point inside polygon?
;   function PixelsInPolygon        - Find all pixels inside a polygon
;   pro DrawPolygon                 - Draw the maskoff polygon
;   pro AddPolygon_EVENTS           - Maskoff polygon drawing
;   pro ConvertCCDimages            - Convert CCD image format
;   pro RestoreCorrectionev         - Restore corrections
;   function wmean                  - calculate weighted mean + its variance
;   function OverlapnProfiles       - Calculate weighted mean of n profiles
;   function Overlap2Profiles       - Calculate overlap for two profiles
;   pro CreateOverlapTable          - Create the overlap table for all visible areas
;   pro CleanUp_AreaRange           - Cleanup routine for AreaRangeWindow
;   pro AreaRange_event             - Event handler for AreaRangeWindow
;   pro AreaRangeWindow             - Define area by numbers
;   pro CleanUp_PointOptions        - Clean up routine for point options overview window
;   pro PointOptions_event          - Event handler for point options overview window
;   pro PointOptions                - Build the point options overview window
;   function CMethodAllow           - Changing method allowed?
;   function CorrAllow              - Correction and/or restoration allowed?
;   pro CleanUp_Corrections         - Clean up routine for corrections overview window
;   pro CorrOptPaint                - Corrections overview window, set sensitive
;   pro AutoCorrect                 - Automatic correction
;   pro ReadWithAutoStart           - Secure dynamic data in list
;   pro ReadWithAutoStop            - Restore dynamic data in list
;   function ReadWithAutoCorrect    - ReadFile + Autocorrect
;   pro Corrections_event           - Event handler for corrections overview window
;   pro CorrectOptions              - Build the Corrections overview window
;   pro CleanUp_ReadOptions         - ReadOptions cleanup
;   pro ReadOptions_event           - ReadOptions event handler
;   pro ReadOptions                 - Reading options
;   pro BatchProc2DTiff             - Sum/average/superimpose patterns
;   function DoublePresent          - Find double x-y coordinates in an array
;   pro GridNodeCoord               - Convert grid nodes indices to grid node coordinates
;   pro GridNodeCoordInv            - Convert grid node coordinates to indices
;   function GridCoordKeep          - Find output points when grid properties changed
;   function GridCoord              - Find corresponding grid coordinate for found grid-spots
;   function AutoGridCoord          - Find corresponding grid coordinate for found grid-spots (autoconnect mode)
;   pro AutoSortTieq                - Sort input point to connect to regular grid (no point missing)
;   pro MP_PROCESS_EVENTS           - reconnect input point for spatial distortion 'down' and 'moving' event
;   pro MC_PROCESS_EVENTS           - reconnect input point for spatial distortion 'down' event
;   pro MC_DRAWSEL                  - reconnect input point for spatial distortion 'moving' event
;   pro MC_DRAWARROW                - reconnect input point for spatial distortion 'moving' event
;   pro MR_PROCESS_EVENTS           - reconnect input point for spatial distortion 'down' event
;   pro MR_DRAWSEL                  - reconnect input point for spatial distortion 'moving' event
;   pro MR_DRAWARROW                - reconnect input point for spatial distortion 'moving' event
;   function ConnectGrid            - Connect Grid Points for output
;   pro RecalcPixelSize             - Calculate pixel size from grid properties
;   pro ChangeGridParam             - Apply all changes needed when one grid parameter changes
;   pro CleanUp_RealGrid            - Clean Up grid data window
;   pro RealGrid_event              - Edit Grid data
;   pro plotcomp                    - Plot fitted surface and datapoints
;   pro backplot                    - Plot background on image scatterplot
;   pro plotsurf                    - Plot surface
;   pro red_res                     - Reduce fitting results (for bad fits)
;   pro CleanUp_AutoLoad            - Clean Up AutoLoad window
;   pro AutoLoad_event              - Event handler for Auto Loading window
;   pro CleanUp_SelMask             - Clean Up SelMask window
;   pro SelMask_event               - Event handler for selecting fields to load from .msk
;   function eigenns                - Compute eigenvectors and eigenvalues Nonsymmetric Array with n Distinct Real and Complex Eigenvalues
;   function ConicSectionLoci       - convert elliptical parameters to conic loci
;   function ConicTEllipse          - convert from conic to elliptical parameters
;   function fitellipse             - LS Fit ellipse and return center, radii and orientation
;   function RefTilt                - Find tilt parameters and beam center for ellipse using LS-fitting of ellipse
;   pro SortPoints                  - Sort from azimuthal 0 to 2pi
;   pro FitTilt                     - Determine tilt parameters with NLLS
;   pro CleanUp_tilt                - Clean Up for tilt fitting
;   pro tilt_event                  - Event handler for inputting options for tilt fitting
;   pro CleanUp_del_ring            - Clean up for deleting Tilt callibration or Tie rings
;   pro del_ring_event              - Event handler for deleting Tilt callibration or Tie rings
;   pro DTIE_PROCESS_EVENTS         - Change PDF d-spacing value for this ring
;   pro removetiepoints             - Remove tie points from list
;   pro removedoubleTieI            - Remove tiu points with double input
;   pro DELTIE_PROCESS_EVENTS       - Delete Tie point for spatial distortion correction
;   pro TIE_PROCESS_EVENTS          - Select Tie point for spatial distortion correction
;   pro DTILT_PROCESS_EVENTS        - Change PDF d-spacing value for this ring
;   pro DELTILT_PROCESS_EVENTS      - Delete point for tilt correction
;   pro TILT_PROCESS_EVENTS         - Select point for tilt correction
;   function fittdGrid              - Fit 2D data (fit 1 peak at a time)
;   pro CleanUp_GridSearch          - Clean up for finding grid points in selected square
;   pro GridSearch_event            - Event handler for finding grid points in selected square
;   pro surffit_event               - Event handler for fitting selected square mask
;   pro CleanUp_surf_int            - Clean up for selecting square mask to fit
;   pro surf_int_event              - Event handler for selecting square mask to fit
;   pro surf_int                    - Select square mask to fit
;   pro CleanUp_rectangle_sum       - Clean up for rectangle summation
;   pro rectangle_sum_event         - Event handler for rectangle summation
;   pro rectangle_sum               - Rectangle summation
;   function bin_findpeaks          - find peaks in binary image
;   pro easy_grid                   - find spatial distortion grid points
;   pro DeletePreviousAutoSet       - Remove all tie points
;   function AutoSetRing            - Set tie point on Debye ring
;   pro plotoimage                  - Plot 2D oimage after azimuthal integration
;   pro CleanUp_az_int              - Clean up for select mask area to intergrate and plot resulting 2D plot
;   pro az_int_event                - select mask area to intergrate and plot resulting 2D plot
;   pro az_int                      - select mask area to intergrate and plot resulting 2D plot
;   pro CleanUp_auto_set            - Clean up for setting Tilt calibration or Tie points
;   pro AddTi                       - Add Tie points to list
;   pro auto_set_event              - select mask area to set Tilt calibration or Tie points
;   pro auto_set                    - select mask area to set Tilt calibration or Tie points
;   function GetFullArc             - Return the arc that covers the whole image
;   pro DisablePDFrings             - Disable PDF rings that fall of the image before calibrating with them
;   pro CalibFromPDFHandler         - Auto set with PDF data
;   pro CalibFromPDF                - Auto set with PDF data
;   pro OBJECT_PROCESS_EVENTS       - Event handler for 'down' event
;   pro ARC_PROCESS_EVENTS          - Event handler for 'down' event when specify arc mask
;   pro CIRCLE_PROCESS_EVENTS       - Event handler for 'down' event when specify circle mask
;   pro SQUARE_PROCESS_EVENTS       - Event handler for 'down' event when specify circle mask
;   pro ELLIPSEM_PROCESS_EVENTS     - Event handler for 'down' event when specify "mark ellipse"
;   pro ARC_DRAW                    - Event handler for moving events when specify arc mask
;   pro CIRCLE_DRAW                 - Event handler for moving events when specify circle mask
;   pro SQUARE_DRAW                 - Event handler for moving events when specify circle mask
;   pro TTWindowUpdate              - Update window
;   pro ellipsewarning              - Ellipse waring update
;   pro ELLIPSEM_DRAW               - Event handler for moving events when specify "mark ellipse"
;   pro AZIOFF_DRAW_MOVE            - AZIOFF_DRAW helper
;   pro AZIOFF_DRAW_MOVE2           - AZIOFF_DRAW helper
;   pro AZIOFF_DRAW                 - Event handler for shifting the azimuth offset
;   pro usemouse_events             - use mouse to determine center of circle
;   pro circle_events               - add 3 points to determine center circle
;   function calc_middle            - calculate center of the circle determined by three points
;   pro tplot                       - plot things in overview window
;   pro dplot                       - plot things in dynamic window
;   pro RefreshScalingSlider        - Refresh image scaling indicators
;   pro RefreshCircleSlider         - Refresh the circle slider
;   pro RefreshDisplay              - Refresh Display when changes made (Main window)
;   pro MIDDLE_PROCESS_EVENTS       - Point middle with mouse
;   pro mode                        - changing program mode
;   pro CleanUp_Expd                - Clean up window for Experimental geometry window
;   pro Expd_event                  - Event handler for Experimental geometry window
;   pro ZOOM_PROCESS_EVENTS         - Event handler for 'down' event when zooming
;   pro ZOOM_DRAWBOX                - Event handler for moving events when zooming
;   pro ZOOM_CHANGE                 - Resize dynamic window for new zoom
;   pro CleanUpInfo                 - Clean up Info window
;   pro info_event                  - Event handler for "about" window
;   pro XRDUASaveINI                - Save Main Window settings
;   pro CleanUp_XRDUA               - Clean up Main window
;   pro XRDUA_event                 - Main event handler
;   pro XRDUA                       - Main procedure
;
;
;  COMMENTS:
;   Graphical parameters:
;   array:  -> size          = tiffs[1] x tiffs[2]
;              -> zoomsize     = (hrange[1]-hrange[0]+1) x (vrange[1]-vrange[0]+1)
;   left:     -> size            = stath x statv
;              -> scaling: sfh = tiffs[1]/stath (this is fixed!)
;   right:  -> size         = dynh x dynv
;              -> scaling: dynh= (hrange[1]-hrange[0]+1)/sf (sf is editable!)
;    middle of the first pixel in "real coord": [0,0] (other convention would be [0.5,0.5])
;    Tilt angles:
;    a and be are the tilt and rotation angles to convert from circles to ellipses.
;    Geometrical parameters in Fit2D:
;        1. image looks the same:
;            afit2d=-a
;            bfit2d=(90-|b|)*(b/|b|)
;        2. image displayed upside down:
;            afit2d=-a
;            bfit2d=(90-|b|)*(-b/|b|)
;            ycfit2d=ny-1-yc
;    Somethings wrong with fit2d. Try opening and saving tiffs without processing.
;    The result will be flipped vertically.
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: Platform
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group device/platform definitions for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   function Version                - Program version
;   function CompareVersion         - compare version numbers
;   function winidvalid             - Check window
;   function winid                  - Returns next non-used ID under 32 if possible, else above 32
;   pro LoadctX                     - load customized color table
;   pro LoadctRev                   - load reversed color tables
;   function enlarge_dyn            - Enlarge dynamic pixmap/draw_widget
;   function create_dyn             - Create dynamic pixmap/draw_widget
;   function CreateDyn              - New or resize dynamic drawing
;   function test_pixmap_size       - Find largest pixmap on current device
;   pro SetXRDUAEnv                 - Set environment variables
;   pro ControlDevice               - Handle device issues before starting up
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: Files
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group all file read/load functions for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   function ReadCCDInfo            - Read CCD options
;   pro mWRITE_CSV                  - Write CSV file
;   function read_tiff_tagvalue     - Get variable and bytecount for a given tag value type
;   pro reform_tiff_format          - Reform image block
;   function read_tiff_structure    - Read tiff file in a TIFF-structure
;   pro sortformat                  - sort filter for pickfile
;   function ChiXtitle              - Different CHI titles
;   function specblockextract       - Read block in spec file
;   function readspec               - Read spec file
;   function TIFFformat             - Check tiff file 1D or 2D
;   function EDFformat              - Check edf file 1D or 2D
;   function SPEformat              - Check spe file 1D or 2D
;   function ChiFormat              - CHI file format for saving
;   function FileFormat             - Get 1D or 2D fileformats or check whether file has 1D or 2D format
;   function SortListFiles          - Sort list of files
;   function ReadFit2DSpline        - Read Fit2D spline file
;   pro writeims                    - Write ims file
;   function readims                - Read ims file
;   pro write_esg                   - Write esg format
;   function WriteChi               - write .chi file
;   function ReadLOF                - Read list of files to include in batch process
;   pro ReformXDIDatablock          - reform xdi datablock orientation
;   function ReadXDI                - Read .xdi files
;   function Select_Files           - Select multiple files with dialog or commandline
;   function ReadOverlapFile        - Read mask overlap file
;   function SelFile                - Select file to write
;   function PSWINDOW_ASPECT        - get PS aspect ratio
;   function PSWINDOW               - get PS size
;   pro insert_cgm_bin              - add line to cgm-file
;   function SaveImageFilter        - Save image formats
;   pro SaveRGBImage                - Save true image
;   pro SaveImage                   - Save plot as an image
;   function ReadEDF                - Read EDF file
;   function WriteEDF               - Write EDF file
;   function read_d                 - read .d file
;   function read_spe               - read .spe file
;   function write_spe              - write .spe file
;   function read_fio               - read .fio file
;   pro write_fio                   - write .fio file
;   function read_mca               - read .mca file
;   pro writemca                    - write .d file
;   function read_chi               - read .chi file
;   function read_plt               - read .plt file
;   function mread_ascii            - read 2-theta scan in ascii format (2 columns)
;   function read_drn               - read 2-theta scan in drn format (2 columns)
;   function read_x_y               - read .x_y or .asc
;   function read_q                 - read .x_y but instead of (2theta,Int.), (Q-space,Int.)
;   function read_nor               - read .nor (XANES)
;   function read_dat               - read .dat
;   function ReadScan               - data input and save data
;   function CHI_xval               - Convert between d-spacing, two-theta and radial (return all three)
;   function SavePKS                - Save 1D and 2D fitresults
;   function WriteLOF               - Write list of files to include in batch process
;   pro CleanUp_BrowserOptions      - Clean Up BrowserOptions Window
;   pro BrowserOptions_event        - Event handler for setting Browser options
;   pro MakeBOptionsWindow          - Window for setting Browser options
;   function ReadINI                - Read XRDUA Setup file
;   function UpdateINI              - Update XRDUA Setup file
;   function SubscribeToOutID       - New window subscription to OutID
;   function UnSubscribeToOutID     - Remove window subscription to OutID
;   function ReadPDDefault          - Load PDF data from ascii file (first column must be d-spacings)
;   function ReadPDDefault2         - Load PDF data from ascii file
;   function ReadDS                 - Load PDF data from Fit2D format
;   function ReadAID                - Load PDF data from .aid file (NBS*AIDS83)
;   function ReadPDD                - Load PDF data from .pdd file (own format)
;   function ReadBEX                - Load PDF data from .bex file (DHN format)
;   function ReadCif                - Read CIF file
;   function ReadCif2               - Read CIF file (PDF)
;   functiom StorePDD               - Store PDF data
;   function FindNextPDD            - Find next PDF data file
;   pro LoadPDD                     - Load PDF data
;   function SaveCel                - Save cel file
;   function ReadCel                - Read structural data from .cel file
;   function LoadStructure          - Load Structure File
;   function SaveStructure          - Save Structure File
;   function READDM3SIMPLETYPES     - DM3 format helper
;   function READDM3STRING          - DM3 format helper
;   function READDM3STRUCT          - DM3 format helper
;   function READDM3ARRAY           - DM3 format helper
;   function READDM3ELEMENT         - DM3 format helper
;   function READDM3TAGTYPE         - DM3 format helper
;   function READDM3TAGENTRY        - DM3 format helper
;   function READDM3TAGGROUP        - DM3 format helper
;   pro DM3TAGLABEL                 - DM3 format helper
;   pro DM3INFO_EVENT               - DM3 information browser
;   function READDM3                - Read DM3 file
;   function ReadRAW                - Read RAW file
;   function ReadSIFblock           - Read 1 image from multipage SIF
;   function ReadSIF                - Read Andor Technology Multi-Channel Files
;   function function RigakuRaxisHeader - Rigaku R-Axis header format
;   function ReadRigakuRaxis        - Read Rigaku R-Axis format
;   pro ParseSPEHeader              - Parse SPE header information
;   function ReadSPE                - Read Princeton(Photonic Science) data format
;   function ReadCCDtiff            - Read diffraction pattern from tiff
;   function ReadCCDjpg             - Read diffraction pattern from jpg
;   function ReadCCDbmp             - Read diffraction pattern from bmp
;   function ReadMarCCD             - Read MarCCD format
;   function ReadSMART              - Read Bruker SMART format
;   function ReadCCD                - Read CCD data
;   function rebinCCD               - Rebin CCD image
;   function WriteMarCCD            - Write MarCCD format
;   function WriteSMART             - Write Bruker SMART format
;   function WriteCCD               - Write CCD data
;   function CopyFields             - Copy fields of .msk file before updating it
;   pro SaveMask                    - Save .msk file
;   function AdaptDirectory         - Change directory relative to another directory
;   pro LoadMask                    - Load .msk file
;   pro RestoreCorrection           - Restore corrections
;   pro OpenFile                    - Initiate dialog when opening file
;   function DetectFiles            - Detect changes in files in directory
;   pro MainUpdateCore              - Update main program core
;   pro MainUpdate                  - Update main program when selecting file from list
;   pro MainUpdateManual            - Update main program manually
;   pro MainUpdateCoreCalc          - Update main program with average of files core
;   pro az_int_auto                 - perform last integration
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: expsetup
;
;  CALLING SEQUENCE: ExpSetup
;
;  PURPOSE: Help with experimental setup
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    pro RefreshSetupWindow        - Refresh text widgets
;   pro ELLIPSEM2_PROCESS_EVENTS - Event handler for 'down' event when specify "mark ellipse"
;   pro ELLIPSE2_DRAW           - Event handler for moving events when specify "mark ellipse"
;   pro splot                      - Plot things in setup window
;   pro RefreshDisplaySetup     - Refresh Display when changes made (Setup window)
;   pro CleanUp_Exps            - Clean up window for Experimental setup window
;   pro Exps_event              - Event handler for Experimental setup window
;   pro ExpSetup                 - Interface for experimental setup
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: xrd_xdi
;
;  CALLING SEQUENCE: XRD_XDI
;
;  PURPOSE: Edit and process results from XRD_BP
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      function getxdicoord        - get sorted xdi coordinates
;      pro CorrectICR                - ParseXiast helper function
;      function ParseXiast            - Load an ESRF XIA scan and make correction files
;      pro XDIdrawfree                - free draw uvalue
;      pro GetXDIdatadim            - get xdi data dimensions
;    function GetGroupLabel        - Get Group label
;    pro RGBtriangle                - Make RGB Chromaticity Triangle
;    function alpha_blending        - Blend groups to create overlap image
;    pro autofindpixels            - Automatic find pixels to redirect
;    function CheckSubCommand    - Subcommand a number or @-sign
;    function ParseSubCommand    - Parse subcommand for image combination
;    function ParseComboCommand    - Parse command for image combination
;    pro ParseIndices            - Parse command for image addition
;    function ImageOverlap        - Create overlap image of different groups
;    function ImageCombine        - Create sum/superimposed of different groups
;    pro injcorrect                - Cut some pixels
;    function SinAngles            - Calculate angle range for sinogram
;    function RhoFromTomogram    - Get rho for each angle for a tomogram coordinate
;    pro ddistCheck                - Check fitted ddist vs. theoretical ddist
;    pro ddistCheck_event        - Event handler
;    pro DRAWSINUS_EVENTS        - Draw sinus in sinogram
;    function ReconRamlak        - Ramlak reconstruction filter
;    function ReconFilter        - Create tomo reconstruction filter
;    pro PixelJuggle                - Do all pixel juggling
;    function TomoConvergence    - Check convergence for iterative tomographic reconstruction algorithms
;    pro ApplyTomogramConstraints- Put constraints in the tomogram
;    function wobbleshift        - Change in projection center due to elliptical wobble
;    function projcenter            - Calculate wobble for projection center
;    pro FBPrecon                - Reconstruct tomogram with FBP
;    pro MLEMrecon                - Reconstruct tomogram with MLEM
;    pro SIRTrecon                - Algebraic reconstruction
;    function ARTRays            - Rico, intercept and intercept-increment for all angles
;    function ARTweightsArea        - Weights for ART: area coverage
;    function ARTweightsLine        - Weights for ART: line segments
;    function ARTrecon            - Algebraic reconstruction technique
;    pro TOMOrecon                - Tomographic reconstruction
;    function Makewplot3DImage    - Make image (+ tomogram) from 3D XDI
;    function GetGroupindices    - Prompt for group indices
;    pro savecolorbar            - Save a colorbar for the current selected group
;   pro wplot3D_tomosep            - plot only tomo things in XRD_XDI window
;   pro wplot3D                    - plot things in XRD_XDI window
;   function wplot2Dgetxrange    - get xrange for wplot2D
;   function ConvertXDIunit        - Convert mm to um
;   pro wplot2D                    - plot things in XRD_XDI window
;    pro wplot2D1                - wplot2D but all on top of eachother
;   pro RefreshDisplayXDI       - Refresh Display when changes made (XRD_XDI window)
;    pro AnimateTomo                - Animate tomographic reconstruction
;    pro RefreshDisplayXDI_LoadctX- Refresh Display callback for LoadCtX
;    pro DeleteXDIGroups            - Delete groups
;    PlotCorrelateXDIGroups        - Plot result of CorrelateXDIGroups
;    pro CorrelateXDIGroups        - Correlation plot of two groups
;    pro AddXDIGroups            - Add groups
;   pro CleanUp_XDIparam        - Clean up window for param window
;   pro XDIparam_event          - Event handler for param window
;    pro XDIparam                - Create param window
;    pro synchrotron_run            - synchrotron current fit function
;    function synchrotron_run_joints - injection joints in synchrotron run
;    pro synchrotron_run_init    - synchrotron current parameters
;    function synchrotron_run_model - model synchrotron current
;    pro WriteDC                    - Write Doris current in an ims file
;    function ReadDC                - Read Doris current from .spe file
;    pro PlotDCCor                - Plot counters selected from list
;    pro DCNormal_cleanup        - Cleanup handler
;    pro DCNormal_refresh        - Refresh DC display
;    pro DCNormal_DRAW            - Handle draw events
;    pro DCNormal_event            - Handle events
;    pro DCNormal_interactive    - Change DC manually
;    pro DCNormal                - Normalize with DC
;   pro XDI2IMS                    - Convert .xdi file to .ims file
;   pro XDI2TIFF                - Convert .xdi file to .tiff file
;   pro XDI2EDF                    - Convert .xdi file to .edf file
;   pro Type3XDI                   - Save datablock as an .xdi file with type3 groups
;    pro EDF2XDI                    - Convert .edf file to .xdi file
;    pro PyMCADAT2XDI            - Convert PyMCA .dat file to .xdi file
;   pro IMS2XDI                    - Convert .ims file to .xdi file
;    pro TIFF2XDI                - Convert several .tiff files to 1 .xdi file
;    pro Type3CHI                - Save CHI files
;    pro XDI2CHI                    - Save XDI as CHI files
;    pro Type3TIFF                - Save TIFF files
;    pro XDIconvert                - Convert different format and XDI
;    pro DRAW_SINUS                - Draw fitted sinus on sinogram image
;    pro SinZingers                - Detect and remove zingers in a sinogram
;    pro BacklashCorrect            - Correct sinogram for backlash
;    pro SinNormal                - Normalize sinogram rows so there sums equal the mean sum
;    pro CleanUp_FindCOPManual    - Destroy handler
;    pro FindCOPManualPlot        - Plot 0 and 180 projection
;    pro FindCOPManual_Event        - Event handler
;    pro FindCOPManual            - Manually set center of projection
;    pro ModelSinogram            - Calculate backlash shift and center of projection
;    pro WorkingDataCorrect        - Correct working data
;    pro printfphysmask            - Create 'PHYSMASK' block
;    pro UpdateXDITD                - Save re-orientated 2D datablock
;    pro UpdateXDIDD                - Save re-orientated 3D datablock
;    pro UpdateXDI                - Save re-orientated datablock
;    pro xdi_setdim                - Calculate XDI plot dimensions
;    pro xdi_changesize            - Change plotting size
;    pro ChangeXDIcoord_event    - ChangeXDIcoord event handler
;    pro ChangeXDIcoord            - Change of the scan coordinates
;    pro EnlargeDataBlock_Event    - Enlarge datablock event handler
;    pro EnlargeDataBlock        - Enlarge datablock
;    pro ConvertDataBlock        - Re-orientate datablock
;    pro CursorPropTomo            - Cursor selected tomo properties
;    pro CursorProp                - Cursor selected file properties
;    pro SelectGroupMain            - Highlight the current group on diffraction pattern
;   pro DDXDIplot               - Plot .xdi file content in 3D
;    pro INSERT_PROCESS_EVENTS    - event handler for moving event in draw widget
;    pro DEFAULT_PROCESS_EVENTS    - event handler for moving event in draw widget
;    pro DEFAULTTOMO_PROCESS_EVENTS - event handler for moving event in draw widget
;    pro XDILEGEND_EVENTS        - event handler for xdi legend move
;   pro UPDATE_PROCESS_EVENTS   - event handler for 'down' event in draw widget
;    function blockprofile        - get profile from image stack
;    function closeregion        - Close vertices region
;   function AddVert            - Add point to vertices list
;   pro printArray_meanstdev    - Print mean+/-1s and stdev+/-1s
;   function wplot2DgetFWHM        - Get FWHM from 1D plot
;   function AddGaussFit        - Fit gaussian to blockprofile
;   pro updatexdicoord            - Ask for pixel size
;   pro GaussFitROIs            - Fit gaussian to 1D ROIs
;    pro EXTRACT1D_PROCESS_EVENTS - event handler for 'down'+moving event in draw widget
;   pro LOF_PROCESS_EVENTS      - event handler for 'down' event in draw widget
;   pro LOF_DRAWBOX             - event handler for moving event in draw widget
;   pro DRAW_POL                 - draw outline of polygon
;   pro ROI_PROCESS_EVENTS      - Event Handler for draw widget
;    pro RemakePatternBlock        - Re-arrange patterns according to orientation settings
;   pro VoxelPattern            - Get pattern from tomogram pixel by reconstructing for each pattern-pixel
;    pro UPDATE_PROCESS_EVENTS2    - Get averaged pattern from tomogram pixel
;   pro OpenXDI                 - Initiate dialog when opening file
;    pro SaveAllXDI                - wplot2D or wplot3D but all groups 1 window
;    pro RefreshWindowXDI        - Refresh the XRD_XDI window
;    pro modeXDI                    - Switch XDI drawing modes
;    pro XRD_XDI_CleanUp            - Cleanup routine
;    pro CutInsertEvent            - Cut/insert pixel event handler
;    pro ShiftRowEvent            - Shift row/column event handler
;    pro bCorEvent                - Sinogram correction event handler
;    pro CorEvent                - Sinogram correction event handler
;    pro AutoFillEvent            - Autofill pixels event handler
;   pro XRD_XDI_event           - Main event handler
;   pro XRD_XDI                 - Main procedure
;
;
;  COMMENTS:
;    Corrections on list.dataorig:
;        1. DC current
;        2. delta orientation
;        3. insert XDI
;    WorkingDataCorrect: corrections on list.data
;        1. Pixel juggling
;        2. Zinger removal
;        3. Backlash correction
;        4. Normalization
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: xrd_bp
;
;  CALLING SEQUENCE: XRD_BP
;
;  PURPOSE: Process series of patterns
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      pro updatefilecounter        - Update files counter text widget
;    function CalcNrGroups        - Number of groups
;    pro Fillblock                - fill datablock with fitted results
;    pro MakeDatablock            - put processed data in datablock
;    pro Save1Dxdiheader            - dscan xdi header
;    pro Save2Dxdiheader            - mesh xdi header
;   function Makefiles          - Make processing file names
;    function SaveXDI            - save xdi file
;    function CornersAzimuth        - Check whether points are within azimuth
;   function CornersMask        - Find corners of mask and make multiplication mask (optional)
;   function ROIlabels            - Get ROI labels
;    pro SaveModelResult            - Save analyzing model in a readable format
;    function FASTMODE2GetSize    - Subimage x and y size fro grid-mode processing
;   function PrepOut            - Prepare OutWindow for Batch processing
;   pro resetmag_xrdbp            - Reset image magnification
;   function BatchInit          - Initialise batch processing
;    function Process2DGroups    - Process type0 and type2 models
;    pro Process1DROI            - Extract ROI from profile for type1 model
;    function Process1DGroupsS    - Process type1 and type2 models
;   function ProcessFile        - Process current file in Batch Processing
;   pro bbbplot                    - plot things in grid batch processing window
;   pro bbplot                    - plot things in fast batch processing window
;   function bplotgetxrange        - xrange for bplot
;   pro bplot                      - plot things in batch processing window
;   pro RefreshDisplayBatch     - Refresh Display when changes made (Batch processing window)
;   pro MainUpdateDo            - 2D/1D window update
;    pro FILESELECT_PROCESS_EVENTS - Select file in output widget
;    pro ROIFASTMODE_PROCESS_EVENTS - Fastmode window process event
;    pro ROIFASTMODE_DRAWBOX        - Fastmode window process event
;    function EndBatchSequenceSave - Complete Mapping or linescan (only the saving part)
;   function EndBatchSequence   - Complete Mapping or linescan
;    pro FastModelEnd            - Convert list for fast mode
;    pro FastModelBegin            - Reconvert list from fastmode
;   function PauseBatchExe      - Pause Batch Execution for making some changes
;   function ResumeBatchExe     - Resume Batch Execution for making some changes
;   pro UpdateCText             - Sort coordinates of scan
;   pro SSession_XRD_BP         - Make image of XRD_BP enviroment (i.e. variables and display)
;   pro CleanUP_XRD_BP          - Clean Up XRD_BP window
;    pro Refresh_XRD_BP            - Refresh XRD_BP window
;   pro RSession_XRD_BP         - Restore image of XRD_BP enviroment (i.e. variables and display)
;   pro XRD_BP_PROCESSSENSITIVE    - set sensitive; xrd_bp base widget
;   pro SelectGroupMainBP        - update selected group in main window
;   pro XRD_BP_event            - Event handler for batch processing interface
;   pro XRD_BP                  - Interface for batch processing
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: widgets
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Widgets oriented routines
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      function widgetparent        - Find parent widget
;    pro DeleteWBaseChilds        - Delete all children of a widget
;    function wheresibling        - Where IDs with sibling ID
;    function find_by_uname        - Find all widgets with the same uname
;    function widgetsize            - Total widget size
;    pro widgetsetsize            - Resize widget
;    function widgetTLB            - Get Top level base widget
;    function widgettreesize        - Get size of whole widget tree
;    function widgettreesizeinit    - Get size of the top widget in an widget tree
;    pro propagateresize            - Propagate Top-level base resize down a widget tree
;    pro baseresizeevent            - Handle resize events
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: fplot_obj
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: 3D plot applet (object oriented)
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      pro fplot_obj_proplotdefault- Default proplot
;    function NormCoord            - Normalize coordinates
;    function NonOrt2Ort            - Non-orthogonal to orthogonal base transformation
;    function Rotx                - Rotation matrix for Rotx
;    function Roty                - Rotation matrix for Roty
;    function Rotz                - Rotation matrix for Rotz
;    function DefTrans            - Make Default transformation matrix
;    pro CleanUp_flplot_obj        - Clean up routine for flplot_obj
;    pro flplot_obj_refresh        - Refresh display
;    pro flplot_obj_event8        - Event handler for type 8
;    pro flplot_obj_event4        - Event handler for type 4
;    pro flplot_obj_event1        - Event handler for type 1
;    pro flplot_obj_event        - Event handler for flplot_obj
;    pro flplot_obj                - Main procedure
;    pro Example_event            - Event handler for Example
;    pro Example                    - Example of call sequence
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: NLLS
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group Non linear Least-Squares fitting routines for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   function gfunctest_wnoise    - Generate profile with noise
;   function gausstest_wnoise    - Generate Gaussian with noise
;   pro gfunc_stepfunc            - Step function
;      pro gfunc_bellfunc            - Two step functions
;   pro gfunc_gauss                - Gauss function with background
;   function PVconstr           - Construct PV profile
;    pro gfuncBackVar            - Analytical form of partial derivatives and fitfunction:
;                                 (Single 1D Gaussian)
;    pro gfuncSin                - Analytical form of partial derivatives and fitfunction(general sinus)
;   pro gfunc1                  - Analytical form of partial derivatives and fitfunction(gauss)
;                                 (hole spectrum: areas+positions+sigma's+coeff-ortpol)
;   pro gfunc2                  - Analytical form of partial derivatives and fitfunction:
;                                 (Single 1D Gaussian with linear background)
;   pro gfunc3                  - Analytical form of partial derivatives and fitfunction(gauss)
;                                 (hole spectrum: areas+positions+sigma's)
;   pro gfunc4                  - Analytical form of partial derivatives and fitfunction(pseudo-voigt)
;                                 (hole spectrum: areas+positions+sigma's+n)
;   pro gfunc6                  - Analytical form of partial derivatives and fitfunction(pseudo-voigt)
;                                 (hole spectrum: areas+positions+sigma's+n's)
;   pro gfunc6c                 - Analytical form of partial derivatives and fitfunction(pseudo-voigt)
;                                 (hole spectrum: areas+positions+sigma's+n's) with constraints
;    pro gfuncTilt0                - Analytical form of partial derivatives and fitfunction
;                                 optional param: a,b,center,dist
;                                 meaning: RotXc(a)�RotZt(b) (a and b=-angles supplied)
;    pro gfuncTilt1                - Analytical form of partial derivatives and fitfunction
;                                  param: d-spacing
;                                 optional param: a,b,center,dist
;                                 meaning: RotXc(a)�RotZt(b) (a and b=-angles supplied)
;    pro gfuncTilt2                - Analytical form of partial derivatives and fitfunction
;                                  param: wavelength
;                                 optional param: a,b,center,dist
;                                 meaning: RotXc(a)�RotZt(b) (a and b=-angles supplied)
;    pro gfuncTilt3                - Analytical form of partial derivatives and fitfunction
;                                  param: wavelength and d-spacing
;                                 optional param: a,b,center,dist
;                                 meaning: RotXc(a)�RotZt(b) (a and b=-angles supplied)
;    pro gfuncTilt4                - Analytical form of partial derivatives and fitfunction
;                                  param: two-theta
;                                 optional param: a,b,center,dist
;                                 meaning: RotXc(a)�RotZt(b) (a and b=-angles supplied)
;   pro gfunc2dg                   - Analytical form of partial derivatives and fitfunction:
;                                 One 2D Gaussian without constraints on parameters
;   pro gfunc2dPVc              - Analytical form of partial derivatives and fitfunction:
;                                 Sum of 2D PVs with constraints on parameters
;   pro gfunc2dPV               - Analytical form of partial derivatives and fitfunction:
;                                 Sum of 2D PVs without constraints on parameters
;   function Marquardt          - Marquardt-Levenberg non-linear least-squares fitting algorithm
;   function MarquardtC         - Marquardt-Levenberg non-linear least-squares fitting algorithm
;                                 with constraints (chisq-surface modification)
;    pro chiplot                    - Plot chi-square space
;    pro ChiConstruct            - Construct chi-square space
;    function MarquardtMap        - NLLS fit with optional constraints and/or restraints
;    function NLLS_norm            - Vector norm with overflow and underflow handling
;    function NLLS_LMpar            - NLLS helper function
;    pro NLLS_tie                - NLLS helper function
;    function MakeAInfo            - Make mpfit parameter info
;    function NLLSPoissonW        - Weights for Poisson statistics
;    function NLLS                - NLLS fitting (more stable than curvefit)
;    
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: xrd_chi
;
;  CALLING SEQUENCE: XRD_CHI
;
;  PURPOSE: Analyze 1D XRD powder data
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    pro CleanUp_ProfileRange    - Cleanup routine for ProfileRangeWindow
;    pro ProfileRange_event        - Event handler for ProfileRangeWindow
;    pro ProfileRangeWindow        - Define a new x-axis for profile
;    pro MASK1D_PROCESS_EVENT    - event handler for 'down' event in draw widget
;    pro MASK1D_DRAW                - event handler for moving events in draw widget
;    pro MakePeaksout            - Make table for refined parameters of peaks
;    pro MakePeaksin                - Make table for initial parameters of peaks
;    pro RefreshPDD_CHI            - Recalculate PDF hights
;   pro LoadPDD_CHI             - Load PDF data for CHI window
;    pro CleanUp_ParamW            - Clean up routine for Parameter Window
;    pro ParamW_event            - Event handler for Parameter Window
;    pro ParamWindow                - Make parameter input window from XRD_CHI
;   pro peaks                      - detect peaks and mark on plot
;   pro FitScan                 - non-linear least squares fitting of powder pattern
;   pro DRAG_DROP_EVENTS        - Event handler for events in draw widget
;   pro DELM_PROCESS_EVENTS     - Event handler for 'down' event in draw widget
;   pro VERPM_PROCESS_EVENTS    - Event handler for 'down' event in draw widget
;   pro modeCHI                 - change mode -> set flags
;    pro CHANGE_LEGENDCHI_EVENT    - Change legend position
;   pro plotCHI                 - plot XRD_CHI window
;   pro plotres                 - residual plot
;   pro peaksearch              - select peak
;   pro fpeaksearch             - select fitted peak
;   pro Add_Peak_Event           - add a peak
;    pro Add_Move_Event            - motion handler when adding peaks
;   pro Delete_Peak_Event       - delete a peak
;    pro Select_Any_Event        - select any pixel
;   pro Select_Peak_Event       - select a peak
;   pro ROI_PROCESS_EVENTS_CHI  - event handler for 'down' event in draw widget
;   pro ROI_drawbox             - event handler for moving events in draw widget
;   pro ZOOM_PROCESS_EVENTS_CHI - event handler for 'down' event in draw widget
;   pro ZOOM_DRAWBOX_CHI        - event handler for moving events in draw widget
;    pro NewUnits                - Change pattern units
;    pro RefreshUnitsCHI            - Refresh XRD_CHI unit related data
;    pro RefreshDataCHI            - Refresh XRD_CHI data
;    pro RefreshWindowCHI        - Refresh XRD_CHI window (not the draw widget)
;   pro RefreshDisplayCHI       - Refresh XRD_CHI draw widget
;    pro BatchConvertCHI            - Convert 1D profile x-axis
;    pro SaveSingleChi            - Save current loded chi file
;   pro OpenCHI                 - Open .chi file
;   pro XRD_CHI_event           - main event handler
;   pro XRD_CHI                 - Main procedure
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: Coord
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group coordinate transforms for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    function RTDynExt            - From real coord to drawdyn-coord: keep pixel size
;   function RTDyn              - From real coord to drawdyn-coord
;   function DynTR              - From drawdyn-coord to real
;   function RTStatExt            - From real coord to draw-coord: keep pixel size
;   function RTStat              - From real coord to draw-coord
;   function StatTR              - From draw-coord to real coord
;    function ConvertFileNumber    - File number when data order different from Moveh,LR,UB
;   function CurPos             - Cursor position on a magnified image
;
;  COMMENTS:
;    There are four coordinates: Static window, Dynamic window, Real (i.e. image array) and Physical ()
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: CommonP
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group common procedures for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   function XRDUAresourceDIR       - XRDUA resource directory helper function
;   function XRDUAconfigDIR         - XRDUA configurtion directory
;   function XRDUAdlmDIR            - XRDUA DLM build directory
;   function XRDUAworkingDIR        - XRDUA DLM build directory
;   function checkrelease           - release version (=IDL VM), debug version (=IDL)
;   pro printLastError              - print last error
;   pro OpWithDuplicate             - dest[ind] op= addvalue with duplicates in ind
;   function CheckFile              - Check whether a file is still being written
;   pro RestoreStruc                - Copy common fields from one structure to the other and take care of dynamic memory allocation
;   function capitalization         - capitalize an array of strings
;   function MDIALOG_PICKFILE       - dialog_pickfile with catch block
;   function CAPDIALOG_PICKFILE     - dialog_pickfile with catch block and capitalization of the filter
;   function strucbytesize          - size of a structure in bytes
;   function arraybytesize          - array size in bytes
;   pro typetobyte                  - convert array of numbers to bytes
;   pro bytetotype                  - convert bytes to an array of numbers
;   function basetodecimal          - digit string to decimal (mantissa,exp)
;   function decomposefloat         - decompose a floating point number (sign,mantissa,base,exp)
;   function floatpres              - machine precision
;   function floatequal             - check floating point equal
;   function floatequalangle        - check floating point equal for angles
;   function TypeConvert            - Type conversion with type code as input
;   pro SetType                     - Reset variable type
;   pro Heap_Copy                   - copy a data structure with pointers in it
;   function openw_safe             - crash safe openw
;   function openr_safe             - crash safe openr
;   function openu_safe             - crash safe openu
;   pro WRITE_TIFF_SAFE             - crash safe write_tiff
;   function Binwidths              - get binwidths for an array
;   pro ConvertTime                 - Convert time to minutes or hours or days
;   pro ConvertSize                 - Convert size to kB, MB, GB, ...
;   function convertCRange          - X-axis/Y-axis/Z-axis CRange convert
;   pro GetCRange                   - X-axis/Y-axis/Z-axis CRange
;   function MakeRange              - Make Array index range
;   function msymbols               - different symbols
;   function DimSize                - Same as size(...,/dim) but with extended dimension output
;   function DimSizeFull            - dimsize with other "size" return values
;   function rebins                 - rebin for string arrays
;   function subarray               - Extract subarray, expand with zero's where needed
;   function value_locater          - handle reverse order problem for value_locate
;   function chunklindgen           - Chunk index generation
;   function chunkindex             - Chunk indexing
;   function chunktotal             - Partial sums
;   pro Drizzle                     - Drizzling algorithm (double histogram)
;   function SetIntersection        - Set Operations: the intersection
;   pro SetSysVar                   - Restore system variables
;   pro GetSysVar                   - Save system variables
;   function PointSysVar            - Pointer to copy of system variable
;   pro delvar2                     - delvar in code
;   pro delstruct                   - delete members from structure
;   pro addstruct                   - add members to structure
;   pro SRDirectGraph               - Save/Restore direct graphics parameters
;   function CutPath                - Find some thing in a pathstring
;   pro AddRegistry                 - Add something to windows registry
;   pro OpenURL                     - Open web link
;   pro CleanUp_PromptNumber        - Prompt for number clean-up routine
;   pro PromptNumber_Event          - Prompt for number event handler
;   function PromptNumber           - Prompt for number window creation
;   function WeakScl                - Scale image
;   function MakeNumber             - Make numbers in string format
;   function stringr                - Make string and remove all white spaces
;   function ScaleArr               - Scale array
;   function wherer                 - where with rounding
;   pro parseoperator               - parse compare operator
;   function compare                - compare with relational operator as string
;   pro compareself                 - compare in place
;   pro setsens                     - set sensitivity of widgets
;   function UniqValues             - Return only unique values from array
;   function ShrinkArray            - Return array without specified cells
;   pro CleanUp_ChooseFromList      - Choose items clean up routine
;   pro ChooseFromList_Event        - Choose items event handler
;   function ChooseFromList         - Choose between items
;   pro strreplace                  - Replace characters in a string
;   pro ResetDim_cleanup            - Cleanup routine
;   pro ResetDim_event              - ResetDim event handler
;   function ResetDim               - Reform dimensions manually
;   function IndexStruct            - Structure indexing like arrays
;   
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: Colorproc
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Color procedures for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      function gencolor            - Generate RGB color with rainbow-spread
;      pro getnamedcolors            - List of RGB colors and names
;      function namedcolors        - Return color
;      pro plotcolorpickerSV        - Color picker SV image
;      pro plotcolorpickerH        - Color picker H image
;    pro CleanUp_ColorPicker        - Cleanup for Color picker
;    pro ColorPickerCall            - Call external proc on color change
;    pro ColorPicker_Event        - Event handler for Color picker
;    function ColorPicker        - Color picker
;    pro AlphaBlendingWithGradient
;    function DefaultGradient
;    function SetGradientTableValues
;    pro GetGradientTableValues
;    function SelectClipMarker
;    pro GradientClipMarker
;    pro PlotGradientImage
;    pro GradientWidgetCall
;    pro RefreshGradientImage
;    pro GradientWidgetEventDraw
;    pro GradientWidgetEventMotion
;    pro GradientWidgetEventTable
;    pro GradientWidget
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: TDProc
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group common 2D analysis procedures for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      function image_equal        - Two images are equal?
;    function Set2DbkCode        - Make 2D background code
;    function Get2DbkCode        - Extract from 2D background code
;    function MakeStrip2DFilter    - Filter for 2D stripping
;   function Strip2D            - Strip peaks so background is left (2D)
;   function TopHatFilter2D     - 2D Top-hat filter of image
;   function estims             - estimate variance of 1D Gaussian peak (uses gfunc 1)
;   function imgpeaks           - Search for peaks in image and estimate there parameters
;   function fittd              - Fit 2D data
;   function cc                 - cubic convolution
;   function bi                 - bilinear interpolation
;   function mInterpolate       - Find intensities for non-integral coordinates
;   function AzIntWarpBatch        - azimuthal integration of arc or ellipse: batch processing
;   function AzIntWarp          - azimuthal integration of arc or ellipse: with interpolation
;   function inverseAzIntWarp    - inverse azimuthal integration
;    function EqualBin            - Rebin an array so the integrand stays the same
;   pro RemoveSatur              - Remove saturation from image
;    pro RemoveZingers            - Remove zingers from an image
;   function CalcTieO           - Calculate output Tie points for warping to correct Spa.Dist.
;    function bsplineBF            - Calculate 1D b-spline Basis Functions
;    function bsplineint2Dcp        - 2D interpolation with b-spline, control points given
;    pro resample_array            - Fant's resampling algorithm
;    pro WarpSpaDistRunSpline    - Spatial distortion correction with spline
;   function WarpSpaDistPrepSpline - Prepare for WarpSpaDistRun (based on Fit2D's spline file)
;   function WarpSpaDistPrep    - Prepare for WarpSpaDistRun
;   pro WarpSpaDistRun             - Spatial distortion correction
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: ODProc
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group common 1D analysis procedures for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   function Savgol1D_r2        - ie Savgol 1D r=2 (specialy for this case => faster)
;   function Strip1D            - Strip peaks so background is left (1D)
;   function TopHatFilter1D     - 1D Top-hat filter
;    function estpeak            - estimate parameters for one detected peak (G model)
;   function estspe             - detect and estimate areas of 1D peaks (PV model)
;   function ortpol             - calculate background
;   function fitod              - fit 1D data with PV functions
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: OTDProc
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group common 1D/2D analysis procedures for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   pro BackSub                 - Perform background subtraction
;    function minterpol            - IDL's interpol without extrapolation
;    pro SetAnySelect            - Change anyselect value
;    pro MEllipseSet                - Change mellipse value
;    pro UnlinkAreas                - Unlink areas
;    pro LinkAreas                - Link areas in groups
;   pro CleanUp_EditMask        - Clean up window for editing mask areas
;   pro EditMaskHandler         - Event handler for editing mask areas
;   pro EditMaskWindowP         - Make part of window for editing or viewing mask
;    pro ExtractLinkNr            - Extract different link numbers
;   pro EditMaskWindow          - Make window for editing or viewing mask
;   function Savgol             - Return filtered spectrum or coefficients to do so
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: LList
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group common linked list procedures for XRDUA package
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    function PrintLL            - Return the str fields of all entrees of linked list
;    pro DestroyLLEntreeI        - Delete specified entrees in linked list
;   function FindLL             - Find entrees in linked list
;    pro ReplaceStrLL            - Replace part of .str field in linked list
;    pro DeleteLLEntree          - Delete entree in linked list
;   pro DestroyLLEntree         - Delete entree in linked list
;   pro DestroyLL               - Delete linked list
;   function UpdateLL           - Add and delete items from linked list
;    pro GenericListInsertAfter    - Add entree
;    pro GenericListInsertBefore    - Add entree
;    pro GenericListDelete        - Delete entree
;    pro GenericListDeleteI        - Delete entrees
;    pro GenericListDestroy        - Destroy list
;    function GenericListSearch    - Get pointer to structure with specified field
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: editoverlap
;
;  CALLING SEQUENCE: EditOverlap
;
;  PURPOSE: Investigate the mask overlap statistics
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    pro IMG_PROCESS_EVENTS        - Select table entry through image
;    function nFracDigits        - Number of digits after decimal point
;    function FormatnFracDigits    - Format number, using n digits after decimal point
;    function RoundWithSD        - Round a value, using it's standard deviation
;    pro PlotGroupProfile        - Plot group profile
;    pro PlotProfiles            - Plot overlapping profiles
;    pro RefreshEditOverlap        - Refresh EditOverlap draw area
;    function XvarEditField        - Create table for a variable
;    pro CleanUp_EditOverlap        - Clean up main window
;    pro XVarEdit_event            - Main event handler
;    pro editoverlap                - Create Main window
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: editPDF
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Investigate the mask overlap statistics
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    pro ApplyEditPDF            - Apply changes
;   pro CleanUp_EditPDF          - PDF viewing window clean up routine
;    function EditPDFWindowP        - Arrange buttons
;    pro RefreshPDFWindowP        - Rearrange buttons
;    pro ChangeButtonState        - Change state of buttons
;   pro EditPDFHandler             - Event handler for PDF viewing window
;   pro EditPDFWindow            - Window for PDF viewing options
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: outlog
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: replace the IDL output log
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   pro printw                  - Print to text widget instead of IDL output
;                                 (must have uvalue=[0,n] where n is the maximum number of lines)
;    pro makeoutlogcleanup        - Clean up outputlog without parent
;    pro makeoutlogevent            - Handle events from outputlog without parent
;    pro makeoutlog                - Make a log with or without parent
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: fitmodel
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: all modelling routines for powder XRD patterns
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    function ExpData            - Return experimental data
;    pro CleanUp_SetScalingToWeight - SetScalingToWeight    cleanup routine
;    pro SetScalingToWeight_Event - SetScalingToWeight event handler
;    pro SetScalingToWeight        - Change the scaling factors according to the weight fractions
;    function TestTypeObsolete    - Find hkl-peaks again after changing cell parameters
;    function ProfileQuantile    - Evaluate the quantile function of probability density functions
;    pro SavePKSTypeX            - Save fit results for typex (pcsiwin format)
;    pro SendMainRefresh            - Send refresh event to xrd_chi
; ----Editing code----
;    pro ModelType10SetCode        - Set type10 general code
;    pro ModelType11SetCode        - Set type11 general code
;    pro ModelType12SetCode        - Set type12 general code
;    pro ModelTypeSetCode        - Set general code of any type
;    pro ModelTypeSetCodeAll        - Set general code of all groups of any type
;    pro MakeColorsTypeASU        - Make colors for type12 ASUtable from derived code parameters
;    pro MakeColorsType            - Make colors for typex table from derived code parameters
; ----Fill group----
;     pro AddModelTypeX            - Add peaks from another group
;     pro MovePeaksFromGroup        - Move peaks from to another group
;    pro FilterDPeaksFromGroup    - Select duplicated peaks
;    pro FilterUPeaksFromGroupHandler - Event handler
;    pro FilterUPeaksFromGroup    - Select user specified peaks
;    pro FilterVPeaksFromGroup    - Select peaks lower than a threshold
;    pro FilterPeaksFromGroup    - Select peaks with non-physical parameters
;    pro DeletePeaksFromGroup    - Delete peak from group
;    pro FillModelType12            - Add peak to a group from type12 model
;    pro FillModelType11            - Add peak to a group from type11 model
;    pro FillModelType10            - Add peak to a group from type10 model
;    function ChangeSymmetry        - Update space group for type 11 and 12
;    pro EditTypeStructure        - Update structural info for type 11 and 12
;    pro ReCalculateHKL            - Recalculate HKL's for this model in the given 2-theta range
;    pro PeakAutomatic            - Add peaks to a group from type10 model
;    pro PeakFromStruc            - Add peaks to group from type11 and typ12 model
;    pro StrucFromPeak            - Save structure file for type12 models
;    pro PeakFromPDFHandler        - Event handler for PeakFromPDF
;    pro PeakFromPDF                - Add peaks from PDF database to a group from type10 model
;    pro PeakManual_Event        - Event handler for manual peak adding to type10
;    pro PeakManualMove_Event    - Event handler for manual peak adding to type10
; ----Edit groupname----
;    pro CleanUp_PromptModelEntreeName - Cleanup routine for groupname prompt
;    pro PromptModelEntreeName_Event - Event handler for groupname prompt
;    function PromptModelEntreeName - Prompt for a groupname
;    pro ChangeModelEntreeName    - Change a groupname
; ----Model interface----
;    function GetTableID            - Get table identifier
;    pro EditTypeTableASU        - Handle table changes for ASU
;    pro EditTypeTable            - Handle table changes
;    pro EditTypeTableev            - Handle table change-events
;    pro typeTable                - Table of Initial/Refined Peak and Global parameters for typex
;    pro typeTableasu            - Table of Initial/Refined ASU parameters from type12
;    function ExcludedBMP        - Tree bitmap for excluded groups
;    pro RefreshMainTab            - Refresh Main tab
;    pro MainCheckEvent            - Handle event from Main tab
;    pro MainDropEvent            - Handle event from Main tab
;    pro InitGlobal                - Initialize global parameters
;    function Type12nASUParam    - Get number of ASU fitparameters for type12
;    pro toggleASU                - Toggle ASU parameters
;    pro updatespacegroup        - Changed spacegroup
;    pro RefinedToInitialType    - Exchange refined and initial parameters
;    function RefinedToInitial    - Exchange refined and initial parameters
;    pro BranchFileEvent            - Handle events for tree files
;    function TypenParam            - Get number of fitparameters for any type
;    function TotalTypenParam    - Get number of fitparameters for all models included
;    pro KillAtomRadiusMult        - Kill AtomRadiusMult text widget
;    pro BranchFileInfo            - Display info for tree files
;    pro StartRefinement            - Start refinement with model
;    pro BranchFolderEvent        - Handle events for tree folders
;    pro BranchFolderInfo        - Display info for tree folders
;    pro BranchInfo                - Display info for tree files/folders
; ----Adding and deleting model types----
;    pro IncludeModelEntree        - Include/exclude a group
;    pro ExcludeAllButThis        - Exclude all groups except active
;    pro IncludeModelType        - Include/exclude all groups of a type
;   pro DestroyModel             - Free memory for model allocation
;    pro CreateModel                - Allocate memory for model
;    pro SaveFitTypex            - Save type10 in mask
;    pro LoadFitTypex            - Load type10 from mask
;    pro SaveFitModel            - Save model in mask
;    pro LoadFitModel            - Load model from mask
;    pro DeleteModelEntree        - Delete a group
;    pro DeleteModelType            - Delete all groups of a type
;    pro DeleteModelTypeSelect    - Delete all groups of a type with prompt for each group
;    pro InsertModelEntree        - Insert a group in the linked list
;    pro DuplicateModel            - Duplicate model
;    pro ConvertModeltypex_10    - Convert between models
;    pro ConvertModeltype12_11    - Convert between models
;    pro ConvertModeltype        - Convert between models
;    pro CreateModeltype12        - Create model of type12
;    pro CreateModeltype11        - Create model of type11
;    pro CreateModeltype10        - Create model of type10
;    pro MoveModelGroup            - Move model
;    pro AddModeltype            - Add model
; ----Global interface----
;    pro CleanUp_EditModel        - Cleanup routine for model editing window
;    pro EditModel_event            - Event handler for model editing window
;    pro AddModelEntree            - Add existing types to tree
;    pro EditModel                - Make window for fit model editing
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: oplotchi
;
;  CALLING SEQUENCE: oplotchi
;
;  PURPOSE: Plot several chi files on each other
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    function ScaleYCHI            - Extract and scale y for profile
;    pro oplotchi_sum            - Make the sum of the currently displayed patterns
;    pro RefreshXVal                - Recalculate x-values
;    function MakePlotType1Img    - Make image from data
;    pro plotchis3d                - Plot data in 3d
;    pro RefreshPDD_OPLOTCHI        - Recalculate PDF hights
;    pro plotchis                - Plot data
;    pro RefreshDisplayCHIS        - Refresh draw widget
;    pro RefreshDataOPLOTCHI        - Refresh XRD_CHI data
;    pro EditOplotchiParam_cleanup- Cleanup routine for EditOplotchiParam
;    pro EditOplotchiParam_event    - Event handler for EditOplotchiParam
;    pro EditOplotchiParam        - Make window for editing parameters
;    pro ConstructDataList        - Construct datalist
;    pro EditChi_cleanup            - Cleanup routine for EditChi
;    pro EditChiPart                - Make part of edit chi plots window
;    pro EditChi_event            - Event handler for EditChi
;    pro EditChi                    - Edit chi plots window
;    pro ModeOChi                - Draw widget events
;    pro CHANGE_XRANGE_EVENT        - Event handler for changing the x-range
;    pro ChangeXRange            - Change the xrange
;    pro CHANGE_XRANGE_DRAW        - Event handler for changing the x-range
;    pro PLOTTYPE1_EVENT            - Event handler for cursor events
;    pro CHANGE_LEGEND_EVENT        - Event handler for changing the legend position
;    pro oplotchi_cleanup        - CleanUp routine for oplotchi
;    pro oplotchi_event            - Event handler for oplotchi
;    pro oplotchi                - Make window for oplot
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: ellipse
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Ellipse drawing, converting,...
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      pro ArcRecalcCorners        - Calculate invalid corners of an arc

;    function ConicTypeTT        - Get conic section type, using tilt angle and twotheta
;    function azMmtoPixel        - Convert azimuth from mm-detector-frame to pixel-detector-frame
;    function azPixeltoMm        - Convert azimuth from pixel-detector-frame to mm-detector-frame
;    function azTiltToCone        - Transform tilted frame azimuth to cone frame azimuth
;    function azConeToTilt        - Transform cone frame azimuth to tilted frame azimuth
;    function AziMinInc            - Get the azimuth in tilted and rotated frame where the ellipses are closest
;    function MinIncr            - Calc minimum increment for drawing an ellipse
;    function azConeRange        - Azimuthal cone range for drawing an ellipse
;    function AzimuthalShift        - Get azimuthal shift in the cone frame due to 2-angle calibration
;    function azTiltToConeShift    - Transform tilted frame azimuth to cone frame azimuth (with the shift)
;    function azConeShiftToTilt    - Transform cone frame azimuth to tilted frame azimuth
;    function azDetRange            - Azimuthal detector range for drawing an ellipse
;    function RealToCalibAngles  - Convert 3-angle detector rotation to 2-angle detector rotation
;    function CalibToRealAngles    - Convert 2-angle detector rotation to 3-angle detector rotation
;    function Getphipixphi0        - Get phi in the detector frame (in pixels) for phi=0 from a given phid for phi=x
;    pro TTPConvertionConstants    - Constants for convertion tt<>pix
;    function PolarizationFactor    - Calculate the polarization factor for source+optics
;    function AttenuationFactor  - Calculate the attenuation factor for a particular sample
;    function DefaultIntCor2Dinfo - Default 2D intensity correction parameters 
;    function LabelPff            - Pff units
;    function IntCor2D             - Convertion factor to go from J/s/m^2 to J/s/rad^2 + elimination of polarization, attenuation, ...
;    function PixelNegR            - Detector pixel position having a negative R
;    function BraggReform        - Bragg helper function
;    function BraggTtoX            - Convert from two-theta to ...
;    function BraggTtoX_I        - Convert pattern from two-theta to ...
;    function BraggDtoX            - Convert from d-spacing to...
;    function BraggDtoX_I        - Convert pattern from d-spacing to ...
;    function BraggPtoX            - Convert from pixels to ...
;    function BraggPtoX_I        - Convert pattern from pixels to ...
;    function BraggXYtoX            - Convert from (X,Y) to (...,phi)
;    function BraggXYtoX_I        - Convert pattern from (X,Y) to (...,phi)
;    function BraggXtoXY            - Convert from (...,phi) to (X,Y)
;    function BraggXtoXY_I        - Convert pattern from (...,phi) to (X,Y)
;    function SemiAxisCoord        - Calculate semi-axis from 2-theta
;    function EllipsePlotCoord    - BraggXtoXY for plotting (force gaps etc...)
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;
;  NAME: NLLSFitModel
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group Non linear Least-Squares fitting routines for fitmodel
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    function ExtractbitCode32    - Split 32bit code
;    function MakebitCode32        - Make 32bit code
;    function ExtractbitCode        - Get a part of the fit code
;    pro SetbitCode                - Set a part of the fit code
;    
;    function GlobalR            - decay from global parameters
;    pro GlobalRpder                - partial derivatives to decay global parameters
;    pro fitGlobalR                - initialize global decay
;
;    function GlobalAsym            - nmix from global parameters
;    pro GlobalAsympder            - partial derivatives to nmix global parameters
;    pro gfuncGlobalAsym            - Analytical form of partial derivatives and fitfunction for fitGlobalFWHM
;    pro fitGlobalAsym            - initialize global nmix
;    
;    function Globalnmix            - nmix from global parameters
;    pro Globalnmixpder            - partial derivatives to nmix global parameters
;    pro fitGlobalnmix            - initialize global nmix
;    
;    function GlobalFWHMSplitPV    - FWHM from global parameters (for SplitPV)
;    pro GlobalFWHMmodSplitPVpder- partial derivatives to FWHM global parameters (for SplitPV)
;    function GlobalFWHM            - FWHM from global parameters
;    pro GlobalFWHMpder            - partial derivatives to FWHM global parameters
;    function GlobalFWHMntchz    - FWHM and n from global parameters (for TCHZ)
;    pro GlobalFWHMtchzpder        - partial derivatives to FWHM global parameters (for TCHZ)
;    pro gfuncGlobalFWHM            - Analytical form of partial derivatives and fitfunction for fitGlobalFWHM
;    pro fitGlobalFWHM            - initialize global FWHM
;
;    function GlobalPos            - Position from global parameters
;    pro GlobalPospder            - partial derivatives to position global parameters
;
;    function GlobalInt            - Intensity from global parameters
;    pro GlobalIntpder            - partial derivatives to intensity global parameters
;
;    pro gfuncOpvII                - Omega6: pearson VII
;    pro gfuncOtchz                - Omega5: Mod-TCH pseudo-voigt
;    pro gfuncOspvhalf            - gfuncOspv helper function
;    pro gfuncOspv                - Omega7: split pseudo-voigt and mod. split pseudo-voigt
;    pro gfuncOpv                - Omega4: pseudo-voigt
;    pro gfuncOlor                - Omega1,2,3: lorentzian, mod1, mod2
;    pro gfuncOgauss                - Omega0: gaussian
;
;    pro TypexExtractVar            - Extract refined and fixed variables for typex
;    pro gfuncTypex                - Analytical form of partial derivatives and fitfunction of typex
;
;   pro gfuncFitModel              - Analytical form of partial derivatives and fitfunction
;
;    pro MakeLabelsFromTypex        - Labels for all parameters in typex
;    pro MakeAFromTypex            - Make parameter/constr/restr array from model typex
;    pro typexincrease            - Expand Type X for batch processing
;    pro ptrfitincrease            - Expand ptrfit for batch processing
;    pro MakeTypexFromA            - Put parameter/constr/restr array into model typex
;
;    function FitInfoLabels        - Fit statistics labels
;    function FitModelLabels        - Labels for all parameters
;    function InitFitModelShrink    - InitFitModel with only refined parameters as output
;    function InitFitModel        - Prepare fit with model
;    pro FinFitModel                - Put fit results in model
;    pro FitModelReset            - Free refit-info structure
;    pro FitModelDerived            - Put derived fit results in model
;    pro FitModel                - fit with model
;    function CreateFitModelInfo    - Create info structure for fitmodel
;    pro FitModelCHI                - interactive fit with model
;    pro FitModelBatch            - batch fitting with model
;    pro SavefitmodelFitResultBegin_Event - SavefitmodelFitResult helper
;    pro SavefitmodelFitResultBegin - SavefitmodelFitResult helper
;    pro SavefitmodelFitResultEnd - SavefitmodelFitResult helper
;    function MakeYforFitModel    - Convert y from 2-theta to current xtype
;    pro SavefitmodelFitResult    - Save fit as a profile
;    function gfuncPeakparamBegin - gfuncPeakparam helper
;    function gfuncPeakparamEnd    - gfuncPeakparam helper
;    function gfuncPeakparam        - Extract separate peak parameters
;    function ScalingFactorMultiplier - Scaling factor multiplier for fit stability
;    pro FitModelDestabilize        - Reset stability parameters
;    pro FitModelStabilize        - Initialize stability parameters
;    function InitScaling        - Initialize global scaling factor
;    function hkllabels            - Position and label of hkl reflections
;    function fitmodel2D            - Generate 2D profile from Fit
;    pro FitModelplot            - Plot model
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: mmath
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Group strict mathematical concepts, not present in IDL
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      function mtotal                - total but preserving the dimensions
;      function percentile            - Calculate x'th percentile (e.g. 50 = median)
;      function SphericalToCartesian - Convert from spherical to cartesian coordinates
;      function CartesianToSpherical - Convert from cartesian to spherical coordinates
;      function InnerProduct        - Inner product of two vectors (for a given metric tensor)
;      function VectorNorm            - Norm of a vector (for a given metric tensor)
;      function CrossProduct        - Cross product of two vectors
;      function TransformMetricTensor - Transformation of a metric tensor under a change of basis
;      function Combinations        - k combinations out of n elements
;      function Permutations        - n permutations out of n elements
;      function kPermutations        - k permutations out of n elements
;      function mPCA                - Principal Components Analysis
;      function chisqr_test        - Pearson's chi-square test
;      function entropy            - Compute the entropy of a distribution p
;      function joint_entropy        - Compute joint entropy of two random variables
;      function mutual_information    - Compute mutual information of two random variables
;      function pearson_correlate    - Compute the Pearson product-moment correlation coefficient
;      function testnormal            - Test whether observed values come from a normal distribution
;      function hashCRC32            - CRC-32 hash
;      function sign                - get the sign of a number
;      function atantan            - calculate atan(tan(x).a) with quadrant concervation
;    function Equal0ModZ            - x == 0 (mod Z)
;    function seed                - Generate seed for random number generators
;    function randxy                - Random number between x and y
;    function randlong            - Random positive LONG number
;    function digamma            - Evaluate the digamma function
;    function dergamma            - Evaluate the derivative of the gamma function
;    function fround                - Floating point rounding
;    function ffloor                - Floating point floor
;    function fceil                - Floating point ceil
;    function sigdigits            - Get significant digits and truncate/round value to the least significant digit
;    pro rounderror                - Round value occording to error
;    function defaulterror        - Get default precision of a value
;    function fconvol            - Convolution with a kernel given in Fourier space
;    function BlockDiag            - Block diagonalize a matrix
;    function UniqueRows         - Eliminate rows that are a multiple of the others
;    function RealUniqueRows        - Eliminate rows that are equal
;    function UniqueCols            - Eliminate columns that are a multiple of the others
;    function RealUniqueCols        - Eliminate columns that are equal
;    function singular            - Check whether a matrix is singular
;    pro QRfact                    - QR decomposition
;    pro QRsolve                    - Solve QR decomposed system (least squares)
;    pro ReducedRowEchelonFrom     - Reduced row echelon form of a matrix
;    pro RowEchelonFrom            - Row echelon form of a matrix
;    function ElimSolve            - GJSolve helper function (solve with free variables = 1)
;    function GJSolve            - Solve by Gauss–Jordan elimination
;    function MatrixColumnSpace    - Get orthonormal basis for the Column space of a matrix
;    pro mSVDC                    - Custom SVDC
;    function MatrixNullSpace    - Get orthonormal basis for the Null space of a matrix
;    function mSVSOL                - Custom SVSOL
;    function SolveLinSystem        - Solve a linear system of equations
;    pro ScaleRoundColumns        - Scale+round column vectors
;    function JordanDecompHelper - Get generalized eigenvectors
;    function mHQR                - Get eigenvalues of a square matrix
;    function JordanDecomp         - Jordan decomposition
;    function ShurDecomp            - Shur decomposition
;    function COvar                - Covariance matrix
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: mmathrational
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Rational arithmetic
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      function GCD_extended        - Greatest common divisor of two numbers + a pair Bézout coefficients (the ones that are minimal in absolute value)
;      function GCD                - Greatest common divisor of two numbers
;      function GCDmore            - Greatest common divisor of several numbers
;      function LCM                - Least common multiple of two numbers
;      function LCMmore            - Least common multiple of several numbers
;      function prime_sieve        - Find primes less or equal than a number
;      function IntFact            - Integer factorization
;      function isrational            - Check whether matrix is rational
;      function make_ratarray        - Make rational matrix
;      function identityrat        - Make rational identity matrix
;      pro rational                - Convert integer matrix to a rational matrix
;      pro printrational            - print a rational matrix
;      pro ReduceNumDenom            - Reduce numerator and denominator of a rational number
;      pro Float2Rational            - Convert floating point to a rational with maximal denominator "base"
;      pro RationalModZ            - Rational number modulo Z (makes it positive too)
;      pro SignNum                    - Put sign in the numerator of a rational
;      function rationalsign        - Get the sign of rationals
;      function ratarray_equal        - Check whether two rational matrices are the same
;      function ratarray_compare    - Compare two rational matrices
;      pro ReduceNumDenomPrime        - Reduce numerator and denominator of a rational, both prime factorized
;      pro ReduceNumDenomPrime2    - Reduce numerator1,numerator2 and denominator of a rational, all three prime factorized
;      pro ratinv                    - Invert rational numbers
;      pro ratneg                    - Negate rational numbers
;      function ratmultfact        - Multiplication of rational numbers (with integer factorization)
;      function ratmultfast        - Multiplication of rational numbers (without integer factorization)
;      function ratmult            - Multiplication of rational numbers
;      function ratdiv                - Division of rational numbers
;      function ratsumfact            - Summation of rational numbers (with integer factorization)
;      function ratsumfast            - Summation of rational numbers (without integer factorization)
;      function ratsum                - Summation of rational numbers
;      function rattotal            - Summation of rational numbers in an array
;      function ratsub                - Subtraction of rational numbers
;      function ratrebin            - Rebin a rational matrix    
;      function ratmattranspose    - Transpose a rational matrix
;      function ratmatinvert        - Invert rational matrix
;      pro ratmattranspose            - Transpose rational matrix
;      function ratmatdet            - Determinant of a rational matrix
;      function ratmatmult            - Rational amatrix multiplication
;      function equalmodn            - integers are equal modulo n (n in Z)?
;      function SolveLinearCongruence - Solve Linear Congruence Equation
;    function prquot                - Quotient and positive remainder
;    function MinAbsGtZero        - First, non-zero, minimum absolute value
;    function find_non_division    - Smith Normal Form helper
;    function SmithNormalForm    - Smith Normal Form of a matrix
;      pro printpol                - Print rational polynomial
;      function PolDegree            - Degree of rational polynomial
;      pro ReducePol                - Remove the highest zero terms of a rational polynomial
;      function CheckZeroPol        - Check whether the rational polynomial equals to zero
;      function polsum                - Summate rational polynomials
;      function polsub                - Subtract rational polynomials
;      function polmult            - Multiply rational polynomials
;      function poldiv                - Divide rational polynomials
;      function GCDPol                - Greatest common divisor of two rational polynomials
;      function PolFact            - Facorize a rational polynomials
;      pro ESTPermute                - Elementary similarity transformation
;      pro ESTMultiply                - Elementary similarity transformation
;      pro ESTAddMultiple            - Elementary similarity transformation
;      pro ESTAddMultiple2            - Elementary similarity transformation
;      function SearchNonZeroRow    - Find the first non-zero elements in a specific column
;      function SearchNonZeroColumn - Find the first non-zero elements in a specific row
;      function MinimalPolbasevec    - The minimal polynomial (with respect to a rational matrix) of a base vector 
;      function CompanionPol        - The associated polynomial to a companion matrix
;      function MinimalPolvec        - The minimal polynomial (with respect to a rational matrix) of a vector 
;      pro MinimalPolSetvec        - The minimal polynomials (with respect to a rational matrix) of sets of base vectors
;      function canonicalBasisvector - Make a rational canonical base vector
;      pro CheckMinimalPolVec        - Check whether a vector has a given minimal polynomial
;      function MinimalPol            - The minimal polynomial of a rational matrix
;      function vecMinimalPol2        - vecMinimalPol helper
;      function vecMinimalPol        - Find a linear combination of base vectors which has particular minimal polynomial
;      pro FormInterm                - Intermediate Frobenius normal form (correct companion matrices on the diagonal)
;      pro TryRightZeroNewCompanion - Try to set the elements next to a companion submatrix to zero
;      pro RightZeroFormInterm        - Convert the Intermediate Frobenius normal form to the Frobenius normal form
;      function NewCompanion        - Elementary similarity transformations to create a new companion submatrix
;      function IsFrobenius        - Check whether this is a Frobenius normal form
;      function FrobeniusDecomp    - Frobenius normal form of an integer matrix
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: symmetry
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: All functions and operation involving crystal symmetry
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    ----Space- and pointgroup database----
;    function ASUfn                - hkl in ASU?
;    function pgdata                - 85 point group setings (32 real)
;    function SGSFSymbol            - Schoenflies symbol <-> spacegroup IT number
;    function shmsymbol            - Short Hermann-Mauguin symbol <-> spacegroup IT number
;    function sgdata                - 530 spacegroup settings (230 real)
;    ----RT operators from Hall symbol of spacegroup----
;    function MakeSeitz            - Make 4x4 Seitz matrix
;    function MakeSeitzRat        - Make 4x4 rational Seitz matrix
;    function rationalString        - Make string description of a quotient e.g. 0.5 becomes "1/2"
;    function rationalStringToFloat    - Evaluate string description of a quotient e.g. "1/2" becomes 0.5
;    function rationalStringToRational - Evaluate string description of a quotient e.g. "1/2" becomes {num:1,denom:2}
;    function isdigit            - Is character a digit?
;    function StringRTop            - Construct string description from RT operator: ax,by,cz
;    function StringRTopRat        - Construct string description from (rational) RT operator: ax,by,cz
;    function RTopString            - Construct an RT operator from a string description
;    function RTopStringAll        - Construct RT operators from a string descriptions
;    function FindMaxDenomStrops - Maximal denominator for a string of symmetry operations
;    function RTopStringRat        - Construct a rational RT operator from a string description
;    function RTopStringRatAll    - Construct rational RT operators from a string descriptions
;    function Remain                - Corrected floating-point remainder
;    function TrnCode            - Translation part of RTop code
;    function TrnCode64            - 64-bit variant
;    function RotCode            - Rotation part of RTop code
;    function RotCode64            - 64-bit variant
;    function CodeType64            - Return 64bit code data type
;    function SortCodeType64        - Sort code type 64
;    function Symop                - Rotations/translation operation from code
;    function Symop64            - 64-bit variant
;    function SymopAll            - Rotations/translation operation from codes
;    function SymopAll64            - 64-bit variant
;    function Symop_code            - Rotations/translation operation encoded in 32bit integer (long in IDL)
;    function Symop_code64        - 64-bit variant
;    function StringRTopcode64    - Construct string description from RT operator code64: ax,by,cz
;    function RTopcodeString64    - Construct an RT operator code64 from a string description
;    function StringRTopcode        - Construct string description from RT operator code: ax,by,cz
;    function RTopcodeString        - Construct an RT operator code from a string description
;    function SymopAllRat        - Rotations/translation operation from code (rational variant)
;    function isGroupModT        - Check whether a set of isometries is a group (mod T)
;    function isGroup            - Check whether a set of isometries is a group
;    function SysAbs                - Is the hkl abscent or not for the spacegroup
;    function GenEquivalentHKL    - Generate equivalent hkl's
;    function HKLequivalent        - Find equivalent hkl's in two lists
;    function SGFromHall            - List of RT operators from Hall symbol
;    function NTHallFromOps        - Not tabulated Hall symbol from RT operators
;    function SGFromNr            - List of RT operators from IT number
;    function SGFromHM            - List of RT operators from Hermann-Mauguin symbol
;    function SGFromOp            - List of RT operators from symm. op.
;    function SGFromTab            - List of RT operators from tabulated setting index (1..530)
;    function SGFromHash            - List of RT operators from hash
;    ----Generate spacegroup structure----
;    pro ProductOps                - For a series of ops, make all combinations with a second series of ops
;    function ExpandOpsList        - Generate all symmetry operations from Hall matrices
;    function primitive_noninv_ops- Return primitive without inversion symmetry operations
;    function inversion_ops        - Return inversion symmetry operations
;    function primitive_ops        - Return primitive incl inversion symmetry operations
;    function centering_ops        - Return lattice centering symmetry operations
;    function expanded_centering_vecs - Return expanded list of centering vector
;    function sgchangebasis        - Change SG operations to new basis
;    function sgmakeprimitive    - Make spacegroup with primitive cell
;    function generator_ops        - Return minimal list of RT matrices needed to describe all (i.e. spacegroup generators)
;    function laue_ops            - Return Laue symmetry operations
;    pro AddSubgroup                - Add a subgroup if not already present
;    pro MakeClosure                - Close a (incomplete) subgroup of the pointgroup
;    function PGConjugacyClasses - Calculate the conjugacy classes of subgroups for a pointgroup
;    function InequalModZ        - Solve inequality mod Z
;    pro AddOrbit                - Add orbit of an affine subspace
;    function TransformSubspace  - Transform a higher dof subspace to a lower dof subspace
;    function CompareSubspace    - Check whether subspaces are the same or not
;    function ChecRTTable        - DeleteSubspace helper function
;    function DeleteSubspace        - Delete equal or special subspace orbits
;    pro ReduceGorbit            - Remove equivalent G-orbits
;    pro AddGorbit                - Add G-Orbit of an affine subspace
;    function wyckofflabel        - Generate label for wyckoff position
;    function wyckoffcodes        - Wyckoff positions
;    pro SetASUwyckoff            - Set wyckoff code for ASU positions
;    function pgrp_ops            - Return pointgroup symmetry operations
;    function LaueGroupHash        - Get Laue group hash from spacegroup
;    function PointGroupHash        - Get Point group hash from spacegroup
;    function Centrosymmetric    - Space group centrosymmetric?
;    function celparamdefault    - Default crystal system cell parameters
;    function celparamconnect    - Which cell parameters are connected?
;    function celparamfix        - Which cell parameters are fixed?
;    function celparamsetcon        - Set connected cell parameters
;    function MetricTensor        - Calc metric tensor in direct/reciprocal space
;    function UCVolume            - Calc volume direct/reciprocal unit cell
;    function DirectToReciprocalU- Convert direct to reciprocal unit cell
;    function CheckEnantiomorphic- Check whether a spacegroup belongs to an enantiomorphic pair
;    function CheckSymmorphic    - Check whether a spacegroup is symmorphic
;    function BravaisType        - Bravais lattice type
;    function SGInfo                - Return all spacegroup information
;    function SpaceGroupAllops    - Get spacegroup operations and hash
;    function SpaceGroupStructure - Make spacegroup structure
;    function SpaceGroupStructure_event - Make spacegroup structure from widget event
;    function SpacegroupDroplist    - Find droplist selection
;    pro typeTablestruc            - Table of Structure parameters for typex
;    pro selectspacegroup_cleanup- selectspacegroup cleanup handler
;    pro selectspacegroup_event    - selectspacegroup event handler
;    function selectspacegroup    - Select spacegroup (GUI)
;    function UnitCelltoRep         - Get unit cell representation w.r.t. natural basis
;    function ReptoUnitCell         - Get unit cell of its representation w.r.t. natural basis
;    pro TransformCelparam        - Transform cell parameters with basis transform
;    function MatrixGroupsEqual    - Check whether two groups are equal (might be duplicates)
;    function SetofEquivIsometries - Check whether two groups of isometries are conjugate under a particular change of frame
;    function SetofEquivIsometriesWrap - Wrapper for SetofEquivIsometries
;    pro SpaceGroupConvertGetChangeofBasis_Event - Event handler
;    function SpaceGroupConvertGetChangeofBasis - User supply of change of basis (optionally also origin shift)
;    pro AcceptChangeofFrame_Event - Event handler
;    function AcceptChangeofFrame    - Accept/decline the change of frame found by SpaceGroupConvert
;    function FindOriginResolveZeroD - Helper for SpaceGroupConvertFindOrigin
;    function SpaceGroupConvertFindOrigin - Try to find an origin shift the makes two settings equivalent (they already have conjugate point groups and related centering vectors)
;    pro SplitGroupofIsometries    - Split the set of centering vectors from a space group
;    pro print1DLattice            - Print a 1D lattice, based on its code
;    function Gen1DLatticefromcoeff - Coding of a 1D lattice
;    function Simplify1DLattice    - Try to describe a 1D lattice as Z+{...}
;    pro SimplifyCenteringVectors - Remove duplicated and (mod Z) a set of centering vectors
;    function lattice1Dequal        - Verify whether two 1D lattices are identical
;    pro CenteringVectorsUnderChangeofBasis - Convert the list of centering vectors under a change of basis
;    function SpaceGroupConvert    - Convert unconventional space group setting to conventional setting
;    function SpaceGroup            - Generate space group structure
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: cryststruc
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: All functions and operation involving crystal structures
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    function AtomPlotObject::Init    - Initializes an atom object.This function returns a 1
;                                          if initialization is successful, or 0 otherwise.
;    pro AtomPlotObject::Cleanup        - Cleans up all memory associated with the atom
;    pro AtomPlotObject::SetProperty    - Sets the value of properties associated with the atom object
;    pro AtomPlotObject::GetProperty    - Retrieves the value of properties associated with the atom object
;    pro AtomPlotObject::Copy        - Make new object instance with same properties
;    pro AtomPlotObject::Print        - Print atom object member properties
;    pro AtomPlotObject::BuildPoly    - Sets the vertex and connectivity arrays for the
;                                    polygon used to represent the atom
;    pro AtomPlotObject__define    - Defines the object structure for an atom object
;      function getptrFit            - Get ptrFit from a group pointer
;      function indCelParam        - Indices of celparameters
;      function indZeroParam        - Indices of zero-shift
;    pro ZOCompress0                - Compress coordinates between 0 and 1: [0,1[
;    pro ZOCompress1                - Compress coordinates between 0 and 1: [0,1]
;    function dfromhkl            - Calc d-spacing from hkl and celparameters
;    function ttfromhkl            - Calc two-theta(degrees) from hkl and celparameters
;    function ttpdercell            - Partial derivative: tt to cell parameters
;    function GenHKL                - Generate all allowed reflections for the spacegroup with there multiplicity and d-spacing
;    function WyckLabel            - Make wyckoff label for ASUxyz
;    function MakeASUpos            - Make an asymmetric unit position
;    pro AddASUpos                - Add an ASU position to the ASU
;    pro ResetWyckASU            - Reset ASU wyckoff labels
;    pro ResetWyckpASU            - Reset ASU wyckoff labels
;    pro KeepUniXYZ                - Keep only the unique positions
;    pro CompressXYZ                - Compress coordinates between 0 and 1 and keep only the unique ones
;    function CoordPolyhedra        - Get the coordination polyhedra for all orbits
;    pro AppendXYZ                - Append coordinates to those of the unit cell to display structure
;    pro AddUnitCellpos            - Add equivalent positions to unit cell
;    function ASUtoUnitCell        - Build a unit cell from a dynamic ASU
;    pro CopyASUstruc            - Copy ASU info to modeltype12
;    function CopystrucASU        - Copy modeltype12 to ASU
;    function CopystrucASU2        - Copy modeltype12 fitarray to ASU
;    pro printfasu                - Print ASU in file
;    pro readfasu                - Read ASU from file
;    pro PlotASUAddLine            - pro PlotASU helper
;    pro PlotASU                    - function for plotting ASU (or unitcell)
;    pro PlotUC                    - plot ASU (or unitcell)
;    function RefinedToInitialASU- Copy refined ASU values to initial and v.v.
;    function ConvertASUtotable    - Convert ASU info to table value
;    pro EditASUposXYZ            - Edit ASU xyz and change xyz for substituted atoms
;    pro EditASUpos                - Edit an ASU position from the ASU
;    pro AddASUposI                - Add an ASU position to the dynamic ASU
;    pro DeleteASUposI            - Delete ASU position from the dynamic ASU
;    function StructureFactorpder- Calculate unit cell structure factor
;    function Afactorpder        - Calculate Attenuation correction
;    function LGPfactorpder        - Calculate Lorentz-Polarization factor
;    function GetPDRelIntpder    - Get PD relative intensities for model type10 fitarray
;    function GetPawleypeakparam - Get tt and I from Pawley
;    function GetPawleyRelIntpder- Get Pawley relative intensities for model type11 fitarray
;    function GetRietveldpeakparam - Get tt and Fh^2 from Rietveld
;    function GetRietveldRelIntpder    - Get Rietveld relative intensities for model type12 fitarray
;    function UnitCellFundamentals12    - Unit cell constants: volume, mass,...
;    function UnitCellFundamentals11    - Unit cell constants: volume, mass,...
;    function UnitCellFundamentals10    - Unit cell constants: volume, mass,...
;    function UnitCellFundamentals    - Unit cell constants: volume, mass,...
;    pro ScalingFromWeight            - Set the scaling factors according to the given weight percentages
;    function ExperimentalKfactor    - Calculate the experimental intensity multiplier
;    function checkPhysicalPhaseInfo - Use physical phase info?
;    function SetPhysicalPhaseInfo    - Set phase info: weight percentage, density,...
;    function GetPhysicalPhaseInfo    - Get phase info: weight percentage, density,...
;    function GetNumberPhysicalPhaseInfo    - Get number of phase info parameters
;    function CalcStoichiometry    - Calculate stoichiometry
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: fundamentals
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Fundamental constants
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    function speedoflight            - c in m/s
;    function planckconstant            - h in eV.sec
;    function avogadroconstant        - Na in 1/mol
;    function electronradius            - Electron radius in m
;    function electronrestenergy        - Electron rest energy in eV
;    function LambdaEnergy            - Convert wavelength to energy or v.v.
;    function LambdaFreq                - Convert wavelenght to frequency or v.v.
;    pro DefaultIon                    - Default ion and repeats for element
;    function ZElement                - convert from atom number and element and back
;    function ICSDRadii                - ICSD ion radii
;    function AtomColRad                - Colors and Radii as in program CrystalMaker
;    function AtScatDis                - Dispersive part of atomic scattering factor
;    function AtScat                    - Atomic scattering factor and correction for anomalous dispersion
;    function AtMassUse                - Get "relative atomic mass" or "standard atomic weight" (Da) or molar mass (g/mol)
;    function MassAtCoef                - Calculate the mass absorption coefficient for a compound at an energy
;    function CohCrossSection        - Calculate the coherent scattering cross-section for a compound at an energy
;    pro SourceEmission                - X-ray tube emission profile
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;  NAME: Extern
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: 3th party procedures and functions
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      function mrandomn                - Generates random samples from a multivariate normal distribution
;    function hist_nd                - Perform an N-dimensional histogram
;    function Bsort                    - "Stable Sort": keeps order of equal values
;    function frebin                    - rebin/congrid alternative
;    pro ploterror                    - Plot with error bars
;    pro oploterror                    - Oplot with errorbars
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; 
;  NAME: plotroutines
; 
;  CALLING SEQUENCE: None
;
;  PURPOSE: general plotting routines
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;      pro powplot                        - plot with y^pow capabilities
;      pro powoplot                    - oplot with y^pow capabilities
;      pro powplots                    - plots with y^pow capabilities
;      pro powxyouts                    - xyouts with y^pow capabilities
;    function GetBackgroundContrastColor - get color contrasting the background
;   pro imgtvscl                    - tvscl with axis
;   pro TVcolorbar                    - Plot colorbar
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; 
;  NAME: skyscan
; 
;  CALLING SEQUENCE: None
;
;  PURPOSE: routines for skyscan in/output
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;   pro SkyScanLogFile                - Write SkyScan logfile
;   pro SkyScanTIFF_write            - Write tiff file from a TIFF-structure
;   pro SkyScanTIFF_insert            - Insert some bytes somewhere in the tiff
;   pro SkyScanTIFF_patch            - Patch tiff file to SkyScan format
;   function SkyScan_SelectFiles    - Select series of files
;   function SkyScan_BockInterpol    - Interpolate 2D array series
;   pro SkyScan_NRecon                - Prepare raw tomography data to SkyScan's NRecon digestible data
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; 
;  NAME: buildprocedures
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Procedures that are also used for building the project
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;    function ProgramRootDir            - Get root directory
;    pro CDProgramRootDir            - Set dir to root directory
;    
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; 
;  NAME: imagereg
;
;  CALLING SEQUENCE: None
;
;  PURPOSE: Procedures for image registration
;
;  INTERNAL FUNCTIONS AND PROCEDURES:
;  function HXIE_inittransfo        - Prepare transformation structure
;  pro HXie_encrect                    - Get encapsulated
;  function HXIE_inversetransfo        - Calculate inverse transformation
;  function HXIE_combinetransfo        - Combine two transformations
;  function HXie_transformation     - Transform (scale,rotate,shift) grayscale image
;  pro HXie_transform                - Transform (scale,rotate,shift) RGB image
;  function HXie_grayscale            - Make grayscale from RGB
;  function HXie_zerof                - Shift FFT result so it is centered
;  pro HXie_HighPassCenter            - High-pass filter for HXie_LogPolar
;  function HXie_LogPolarBase        - Base for the log-polar transformation of an image
;  pro HXie_LogPolarTrans            - Log-polar transformation of an image
;  pro HXie_ScaleRot2shift            - Convert images so that scal+rot becomes shift
;  function HXie_FourierFreq        - Get Fourier space frequencies in 1 dimension
;  pro HXie_FourierShift            - Find shift between 2 images
;  function HXie_register            - Find scale, rotation and shift between 2 images
;  pro HXie_resizebeforeinsert        - Insert one image in another
;  pro HXie_blend_nan                - Handle NaNs in blending
;  function HXie_blend                - Blend two images
;  function HXie_regandblend        - Register, transform and blend 2 images
;  pro HXie_tvsclblend                - Show blended image
;    
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;-


























; THINGS BELOW MIGHT BE OUTDATED!!!!!!!!!!!!!!!!!!!

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;   Input Formats:
;     *.tif    Tagged Image File Format         - General output of CCD Camera
;     *.mccd   mar ccd                           - Output from a MarCCD Camera
;     *.xrd or *.001 ...                         - SMART Camera output
;     *.pdd    powder diffraction data          - File based on info from PDF database
;     *.plt    PLOTSO file (created by GADDS)   - Azimuthal integrated data
;      *.aid       NBS*AIDS83 format                - PDF data from PCPDFWIN
;      *.bex       DHN format                        - PDF data from DHN
;      *.cel       PowderCell format                - Structure file
;
;   Output Formats:
;     *.pks    peaks file                          - Save 1D or 2D fitresult of a ROI
;     *.jpg    Joint Photographic Experts Group - Save display to this format
;     *.rgn    Region file                      - Regions in XRD-mapping
;
;   Input and Output:
;     *.msk    mask file                           - File to save relevant program status
;     *.chi    FIT2D .chi file                  - Azimuthal integrated data
;     *.xdi    X-ray Diffraction Image          - Batch processing output
;     *.lof    List of Files                       - List of files to include in batch process
;     *.ini    Configuration file                - Some fixed paths and value to start XRDUA
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;listtype:
;    0/1 xrdua
;    10/11: xrd_chi
;    20/21: xrd_xdi
;    30/31: oplotchi
;    40/41: xrd_bp
;    50/51: expsetup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

; Analyzing Model definitions

; Main structure:
fitmodel -> createmodel

; 0: 2D ROI summation
; => cryst. info: none
; => 1 group = a number of 2D ROIs
; => fit parameters: none
; => growth: no
xrd_bp -> batchinit

; 1: 1D ROI summation
; => cryst. info: none
; => 1 group = a number of 1D ROIs
; => fit parameters: none
; => growth: no
xrd_bp -> batchinit

; 2: 2D/1D ROI summation of whole pattern/profile minus the other ROI's (-> rest)
; => cryst. info: none
; => 1 group = the rest image
; => fit parameters: none
; => growth: no
xrd_bp -> batchinit

; 3: General purpose (e.g. XRF)
; => cryst. info: none
; => 1 group = general
; => fit parameters: none
; => growth: no
xrd_xdi -> Type3XDI

; 10: PD (pattern decomposition)
; => cryst. info: none
; => 1 group = a number of peaks
; => fit parameters:
;                peak params: position, I, FWHM, shape, assym.
;                global peak params: zero, isodeform, scaling S, FWHM, shape, assym.
; => growth: change in number of peaks
fitmodel -> CreateModeltype10

; 11: Pawley(/LeBail)
; => cryst. info: Space group (=> hkl and cell param)
; => 1 group = 1 phase
; => fit parameters:
;                peak params: I, FWHM, shape, assym.
;                phase params: unitcell params, scaling S, FWHM, shape, assym.
; => growth: change in number of peaks (when selecting other spacegroup)
fitmodel -> CreateModeltype11

; 12: Rietveld
; => cryst. info: Space group (=> hkl and cell param) and ASU (=> atom positions)
; => 1 group = 1 phase
;                peak params: FWHM, shape, assym.
;                phase params: scale factor, atomic pos., (an)isotropic thermal param, unitcell params, FWHM, shape, assym.
; => growth: change in number of peaks (when selecting other spacegroup)
fitmodel -> CreateModeltype12

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%