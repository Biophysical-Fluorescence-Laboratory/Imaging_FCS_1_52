
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.FileInputStream;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.BorderLayout;
import java.awt.Frame;
import java.awt.Component;
import java.awt.event.MouseListener;
import java.awt.event.KeyListener;
import java.awt.event.ActionListener;
import java.awt.event.ItemListener;
import java.awt.event.AdjustmentListener;
import java.awt.event.MouseEvent;
import java.awt.event.KeyEvent;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.AdjustmentEvent;
import java.util.Map;
import java.util.Collection;
import java.util.Iterator;
import java.util.Date;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.text.SimpleDateFormat;
import javax.swing.JFrame;
import javax.swing.JTextField;
import javax.swing.JComboBox;
import javax.swing.JToggleButton;
import javax.swing.JButton;
import javax.swing.JRadioButton;
import javax.swing.JTable;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTabbedPane;
import javax.swing.JScrollPane;
import javax.swing.SwingWorker;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.Timer;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.HistogramWindow;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.ScrollbarWithLabel;
import ij.plugin.PlugIn;

import imagescience.random.UniformGenerator;
import imagescience.random.GaussianGenerator;
import imagescience.random.PoissonGenerator;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.AbstractCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.NonSquareMatrixException;
import org.apache.commons.math3.linear.NonSymmetricMatrixException;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;

import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import gpufitImFCS.GpufitImFCS;
import gpufitImFCS.GpufitImFCS.*;

/* *
 *  Imaging_FCS
 *  Version 1.52
 *  11.11.2019
 *  
 *  Imaging_FCS is a program to calculate and display spatio-temporal 
 *  correlation functions from a 16 bit stack file.
 *  
 *  Acknowldegements: 
 *  The original program was written by Thorsten Wohland. 
 *  Radek Machan was responsible for many bug fixes in the newer versions and work on the DC-FCCS part. 
 *  Jagadish Sankaran wrote the MSD calculation algorithm and improved the data fitting part of the program. 
 *  Mark Bathe and Syuan-Ming Guo wrote the first GLS and Bayesian Model Selection algorithms in Matlab (http://fcs-bayes.org/) 
 *  and they helped Thorsten with the coding of the algorithm in Java. 
 *  Jagadish Sankaran and Tang Wai Hoh wrote the CUDA codes that enables computations on NVIDIA GPU. 
 *  Many people in the lab checked the program for bugs which improved the program considerably. 
 *  Any feedback or help in improving the program is welcome.
  *  
 *  References:
 *  
 *  Website with derivation of Imaging FCS and FCCS functions and some interactive apps
 *  http://www.dbs.nus.edu.sg/lab/BFL/ (see the Resources tab)
 *  
 *  Software correlator:
 *  Schätzel et al., 1988, J. Mod. Opt., 35: 711 718
 *  Wohland et al., 2001, Biophys. J., 80: 2987–2999
 *  The actual algorithm for the calculation of the correlations was changed according to the MatLab code 
 *  by Guo and Bathe which is available at http://fcs-bayes.org/. The code was checked against the original
 *  of the above papers and they provide identical results.
 *  
 *  Derivation of Imaging FCS models
 *  Sankaran et al., 2009, Biophys. J., 97: 2630–2639
 *  
 *  Derivation of observation volumes for ITIR-FCS and SPIM-FCS
 *  Sankaran et al., 2013, Anal Chem. 85: 3948−3954
 *  
 *  Diffusion law and PSF calibration
 *  Wawrezinieck et al. 2005, Biophys. J., 89: 4029-4042
 *  Bag et al., 2012, ChemPhysChem, 13: 2784–2794
 *  Bag et al., 2013, BBA – Biomembranes, 1838: 802–813
 *  Bag et al., 2016, Methods Appl. Fluoresc. 4: 034003
 *  
 *  SPIM-FCS and SPIM-FCCS
 *  Singh et al., 2013, Opt. Express, 21: 8652
 *  Krieger et al., 2014, Opt. Express, 22: 2358
 *  Krieger et al., 2015, Nat. Protocols, 10(12): 1948-1972 
 *  
 *  Explanation of original Imaging FCS program in Igor
 *  Sankaran et al., 2010, Opt. Express, 18: 25468–25481
 *  
 *  Bayes data fitting (introduction and application)
 *  Raftery 1995, Soc. Meth., 25: 111–164
 *  He et al., 2012, Anal. Chem. 84: 3871–3879
 *  Guo et al., 2012 ,Anal. Chem. 84: 3880–3888
 *  Guo et al., 2014, Biophys. J., 106: 190–200
 *  
 *  Blocking transformation
 *  Flyvbjerg and Petersen, 1989, J. Chem. Phys. 91: 461
 *  
 *  Covariance shrinkage
 *  Schäfer and Strimmer, 2005, Statistical Applications in Genetics and Molecular Biology, 4(1): Article 32
 *  The code for covariance matrix shrinkage is based on their original code http://strimmerlab.org/software.html 
 *  but has been vectorized and simplified by Kevin Murphy. The verison here is derived from the Matlab code
 *  from http://fcs-bayes.org/
 *  
 *  Generalized Least Squares
 *  http://en.wikipedia.org/wiki/Generalized_least_squares or
 *  Seber and Wild, Nonlinear Regression, 1989, Wiley & Sons (and 2005 Wiley Online Library)
 *  
 *  Mean Square Displacement
 *  Shusterman, R. et al., 2004, Phys. Rev. Lett. 92: 048303
 *  
 *  N&B analysis
 *  Digman et al., 2008, Biophys. J., 94: 2310-2332
 *  Unruh et al., 2008, Biophys. J., 95: 5385-5398
 *  
 *  Gpufit: An open-source toolkit for GPU-accelerated curve fitting
 *  Przybylski, A. et al., 2017, Scientific Reports, 7, Article number: 15722
 *  https://github.com/gpufit/Gpufit
 *  https://gpufit.readthedocs.io/en/latest/
 *
 * */
public class Imaging_FCS_1_52 implements PlugIn {

    private static final String VERSION = "v1.52";

    // User definable parameters
    private final String $UFrameTime = "0.001";
    private final String $CorrelQ = "8";
    private final String $UPixelSize = "24";
    private final String $UMagnification = "100";
    private final String $UNA = "1.49";
    private final String $UEmWavelength1 = "515";
    private final String $UEmWavelength2 = "600";
    private final String $ULateralPSF1 = "0.8";
    private final String $ULateralPSF2 = "0.8";
    private final String $UAxialPSF1 = "1000000";
    private final String $UAxialPSF2 = "1000000";
    private final String $workingDir = System.getProperty("user.home");
    private final int decformat = 2;				// Decimal Format for the Fit panel
    private final int decformat2 = 4;				// Second Decimal Format for the Fit panel, mainly used for the parameter G
    private final int fitMaxIterations = 2000;	// maximum number of iterations allowed in the fits; the maximum number is Integer.MAX_VALUE but that can lead to very long evaluations with typically little improvement on the results
    private final int fitMaxEvaluations = 2000;	// maximum number of evaluations allowed in the fits; the maximum number is Integer.MAX_VALUE but that can lead to very long evaluations with typically little improvement on the results
    private final int BCmaxorder = 8;				// highest polynomial order allowed in polynomial bleach correction
    private final int minFrameReq = 100; 			// a warning is given if less than minFrameReq frames are avaialable for any point in the CF
    private final int swMinFrameReq = 20;			// minimum number of frames required for the sliding windows; this is used to calculate a useful correlatorq
    private final int minDLPoints = 4;			// the square of this value determines the minimum number of points needed for diffusion law; e.g for a value of 3, the last point in the dffusion law plot will be calcualted from 3^2 = 9 D values.
    private final int DiffLawMaxPoint = 30; 		// upper limit on number of points in diffusion law plot

    //plotting options
    private boolean plotACFCurves = true;
    private boolean plotSDCurves = true;
    private boolean plotIntensityCurves = true;
    private boolean plotResCurves = true;
    private boolean plotParaHist = true;
    private boolean plotCovmats = false;								// Regularized covariance matrix will not be plotted by default
    private boolean plotBlockingCurve = false;						// blocking curve will not be plotted by default	

    // Image window
    ImagePlus imp;
    ImageCanvas impcan;
    ImageWindow impwin;
    ImageProcessor impip;
    private int width;			// width of loaded stack
    private int height;			// height of loaded stack
    private int frames;			// number of frames in stack
    private int impmin;			// minimum value in the stack
    private double scimp;		// scaling factor to adjust window to an acceptable size for the user
    private String $impTitle;	// title of the window; this will be used in all dependent windows to differentiate them in case multiple instances of the plugin are open

    // Background windows
    ImagePlus bgrimp;			// background window recorded with the same settings and at same position on camera
    ImageProcessor bgrip;
    ImageWindow bgrWin;
//	ImagePlus bgrCorimp;		// the corrected background window containing the StDev of the camera; this does not contain patterns anymore and offset is subtracted (see Hirsch article) 
//	ImageProcessor bgrCorip;
//	ImageWindow bgrCorWin;

    // Parameter map stack window; in this window the fitted parameters are depicted (e.g. diffusion coefficient maps)
    ImagePlus impPara1;
    ImageCanvas impPara1Can;
    ImageWindow impPara1Win;
    private String $impPara1Title;
    boolean impPara1exists = false;

    // Diffusion Law map window
    ImagePlus impDLMap;
    ImageCanvas impDLMapCan;
    ImageWindow impDLMapWin;
    private String $impDLMapTitle;

    // Covariance window; this window is not displayed by default; diplay can be switched on by stting plotCovmats = true
    ImagePlus impCov;
    ImageWindow impCovWin;
    ImageCanvas impCovCan;
    ImageProcessor impCovIp;
    private String $impCovTitle;

    //DCCF Window; note that multiple instances can be open as dCFF can be calculated in different directions
    ImagePlus impDCCF;
    private final int dccfMax = 4;										// How many dCCF windows can the user open? Set to 4 (0...3) since dCCF can be calculated in 4 directions.
    HistogramWindow[] histDCCFWin = new HistogramWindow[dccfMax];		// Windows to display the dCCF histograms
    ImageWindow[] impDCCFWin = new ImageWindow[dccfMax];				// Windows to display the dCCF images
    private String $dccfTitle;											// titles of ImageWindows for dCCF
    private final String[] $histDCCFWinTitle = new String[dccfMax];		// titles of HistogramWindows for dCCF
    private boolean[] dccfCalculated = new boolean[dccfMax];		// an array to remember in which direction has been dCCF calculated
    //private final int[] histDCCFBin = new int[dccfMax];					// This is an array to remember which dCCF window was produced with which binning

    // N&B windows
    ImagePlus impN;
    ImagePlus impB;
    ImageWindow impNWin;
    ImageWindow impBWin;
    private String $impNTitle;
    private String $impBTitle;
    ImagePlus impNum;
    ImagePlus impEps;
    ImageWindow impNumWin;
    ImageWindow impEpsWin;
    private String $impNumTitle;
    private String $impEpsTitle;

    // Plot windows
    PlotWindow acfWindow;		// ACF
    PlotWindow intWindow;		// Intensity Trace
    PlotWindow resWindow;		// Residuals
    PlotWindow sdWindow;		// Standard Devition of the CFs
    PlotWindow msdWindow;		// mean square displacement
    PlotWindow difflawWindow;	// Diffusion Law window
    PlotWindow PSFWindow;		// Window for PSF determination
    PlotWindow paraCorWindow;	// Window for Parameter Scatter Plots
    PlotWindow blockingWindow;	// Window for blocking analysis; not shown by default; to show, set plotBlockingCurve = true
    private String $acfWindowTitle;
    private String $intWindowTitle;
    private String $resWindowTitle;
    private String $sdWindowTitle;
    private String $msdWindowTitle;
    private String $difflawWindowTitle;
    private String $PSFWindowTitle;
    private String $paraCorWindowTitle;
    private String $blockingWindowTitle;

    // Initial dialog box
    private final int splashScreenDimX = 300;
    private final int splashScreenDimY = 200;
    private final int splashScreenDuration = 3000; // milliseconds

    //Window and panel positions and dimensions; this has been adjusted for a 1920x1080 screen
    private final String $panelFont = "SansSerif";							// font and font size of the Panels
    private final int panelFontSize = 12;
    private final int panelPosX = 10;										// control panel, "ImFCS", position and dimensions
    private final int panelPosY = 125;
    private final int panelDimX = 370;
    private final int panelDimY = 370;

    private final int fitPanelPosX = 10;									// control panel, "ImFCS Fitting", position and dimensions
    private final int fitPanelPosY = panelPosY + panelDimY + 10;
    private final int fitPanelDimX = 370;
    private final int fitPanelDimY = 370;

    private final int simPanelPosX = panelPosX + panelDimX + 10;			// control panel, "Simulation", position and dimensions
    private final int simPanelPosY = 125;
    private final int simPanelDimX = 370;
    private final int simPanelDimY = 320;

    private final int difflawPanelPosX = panelPosX + panelDimX + 300;		// control panel, "Diffusion Law", position and dimensions
    private final int difflawPanelPosY = 125;
    private final int difflawPanelDimX = 350;
    private final int difflawPanelDimY = 150;

    private final int NBPosX = panelPosX + panelDimX + 300;		// control panel, "N&B", position and dimensions
    private final int NBPosY = 125;
    private final int NBDimX = 250;
    private final int NBDimY = 150;

    private final int impwinPosX = panelPosX + panelDimX + 10;				// image window positions
    private final int impwinPosY = panelPosY;

    private final int acfWindowPosX = panelPosX + panelDimX + 10;			// ACF window positions and size, with correlation function plots
    private final int acfWindowPosY = panelPosY + 335;
    private final int acfWindowDimX = 200;
    private final int acfWindowDimY = 200;

    private final int intWindowPosX = acfWindowPosX;						// intensity trace window positions and size
    private final int intWindowPosY = acfWindowPosY + acfWindowDimY + 145;
    private final int intWindowDimX = acfWindowDimX;
    private final int intWindowDimY = 50;

    private final int sdWindowPosX = acfWindowPosX + acfWindowDimX + 115;	// fit standard deviation window position and size
    private final int sdWindowPosY = acfWindowPosY;
    private final int sdWindowDimX = acfWindowDimX;
    private final int sdWindowDimY = 50;

    private final int msdWindowPosX = sdWindowPosX + sdWindowDimX + 165;	// msd window position and size
    private final int msdWindowPosY = acfWindowPosY;
    private final int msdWindowDimX = acfWindowDimX;
    private final int msdWindowDimY = acfWindowDimY;

    private final int resWindowPosX = sdWindowPosX;							// fit residuals window position and size
    private final int resWindowPosY = sdWindowPosY + sdWindowDimY + 145;
    private final int resWindowDimX = acfWindowDimX;
    private final int resWindowDimY = 50;

    private final int difflawWindowPosX = acfWindowPosX + 30;				// Diffusion Law window position and size
    private final int difflawWindowPosY = acfWindowPosY + 30;
    private final int difflawWindowDimX = 200;
    private final int difflawWindowDimY = 200;

    private final int PSFWindowPosX = acfWindowPosX + 30;					// PSF window position and size
    private final int PSFWindowPosY = acfWindowPosY + 30;
    private final int PSFWindowDimX = 200;
    private final int PSFWindowDimY = 200;

    private final int paraCorWindowPosX = acfWindowPosX + 30;				// Parameter Scatter Plot Window position and size
    private final int paraCorWindowPosY = acfWindowPosY + 30;
    private final int paraCorWindowDimX = 200;
    private final int paraCorWindowDimY = 200;

    private final int blockingWindowPosX = sdWindowPosX + sdWindowDimX + 110;	// Parameter Blocking Plot Window position and size
    private final int blockingWindowPosY = sdWindowPosY;
    private final int blockingWindowDimX = 200;
    private final int blockingWindowDimY = 100;

    private final int para1PosX = acfWindowPosX + acfWindowDimX + 80;		// parameter map window positions
    private final int para1PosY = panelPosY;

    private final int DCCFPosX = impwinPosX + 50;							// DCCF window positions
    private final int DCCFPosY = impwinPosY + 50;

    private final int NPosX = impwinPosX + 280;								// Number (N&B) window positions
    private final int NPosY = impwinPosY;

    private final int BPosX = NPosX + 280;									// Brightness (N&B) window positions
    private final int BPosY = NPosY;

    private final int NumPosX = impwinPosX + 280;							// Number (N&B) window positions
    private final int NumPosY = impwinPosY + 30;

    private final int EpsPosX = NPosX + 280;								// Brightness (N&B) window positions
    private final int EpsPosY = NPosY + 30;

    private final int filteringPanelPosX = acfWindowPosX + acfWindowDimX + 100;			// control panel, "Thresholds settings", position and dimensions
    private final int filteringPanelPosY = acfWindowPosY;
    private final int filteringPanelDimX = 420;
    private final int filteringPanelDimY = 300;

    private int histDimX = 350;												// histogran windows positions and size
    private int histDimY = 250;												// the maximum number of parameters is defined in histnum
    private int histPosX = para1PosX + 280;
    private int histPosY = para1PosY;

    private final int histDCCFDimX = 350;
    private final int histDCCFDimY = 250;
    private final int histDCCFPosX = DCCFPosX + 280;
    private final int histDCCFPosY = DCCFPosY;

    private final int covPosX = blockingWindowPosX;						// covariance matrix window position
    private final int covPosY = blockingWindowPosY + blockingWindowDimX + 50;

    // zoom factor: The factor is used to determine the size of the windows on screen; this can be adapted 
    private final int zoomFactor = 250;

    // parameters determined form the open image stack and defined by the user
    private int firstframe;		// first frame to be used in the correlation
    private int lastframe;		// last frame to be used in the correlation
    private int binningX = 1;
    private int binningY = 1;
    private int binning = 1;
    private int cfXDistance;	// x,y distance between pixels to be correlated; in DC-FCCS this is the distance between corresponding pixels in red and green channels
    private int cfYDistance;
    private int cfXshift;		// x,y distance used in fitting of FCS correlations between two pixels; (0, 0) is for ACFs
    private int cfYshift;
    private double correlatorp;	// parameters for the correlator structure; correlatorp refers to the number of bins in the first channel group.
    private double correlatorq;	// all higher groups have correlatorp/2 channels; correaltorq refers to the number of higher groups
    private double frametime;	// acquisition time per frame
    private double objmag;		// microscope objective magnification; note the variable 'magnification' is used for the magnification of the stack image in ImageJ
    private double pixelsize;	// pixel size in micrometer on camera chip (variable "a")
    private double pixeldimx;	// pixel size in object space in [m] used in fit
    private double pixeldimy;	// pixel size in object space in [m] used in fit
    private double emlambda;	// emission wavelength
    private double emlambda2;	// emission wavelength 2
    private double NA;			// numerical aperture
    private double sigma;		// lateral PSF factor: sigma * lambda / NA
    private double sigmaZ;		// axial PSF factor: sigma * lambda / NA (for simplicity)
    private double sigma2;		// lateral PSF factor: sigma2 * lambda2 / NA
    private double sigmaZ2;		// axial PSF factor: sigma2 * lambda2 / NA (for simplicity)
    private double psfsize;		// actual lateral PSF in [m] for laser 1, used in fit; this is the 1/e2 radius
    private double psfsize2;	// actual lateral PSF in [m] for laser 2, used in fit; this is the 1/e2 radius
    private double rz;			// offset between the two light sheets in axial direction (of the detection objective);
    private double lsthickness;	// light sheet thickness for SPIM for laser 1 given as 1/e2 radius
    private double lsthickness2; // light sheet thickness for SPIM for laser 2 given as 1/e2 radius
    private double ccfrx;		// distance in x direction in pixels between pixels to be correlated
    private double ccfry;		// distance in y direction in pixels between pixels to be correlated
    private int background;		// background is determined from the smallest count value from the stack; this can be maually changed in the control panel
    private int background2;	// for FCCS background for the second region is determined as the minimum value for that region
    private int initposx;		// x pixel position at which correlation is calculated after pressing "Single"
    private int initposy;		// y pixel position at which correlation is calculated after pressing "Single"
    private double q2;			// fit parameter, ratio of brightness of second component to first component (if it exists)
    private double q3;			// fit parameter, ratio of brightness of third component to first component (if it exists)
    private int slidingWindowLength = 0;	// sliding window size for bleach correction
    private int polyOrder = 0;	// polynomial order for bleach correction
    private int fitstart = 1;	// parameter for settign the range of points of the CF to be fitted
    private int fitend = 1;
    private double fitobsvol;  // normalization factor including the observation volume; for ACFs this is correct; note that N for spatial cross-correaltions has no physical meaning
    private int isgpupresent = 0;
    private int isdlawcalculatedingpu = 0;
    // variables to remember last settings in the ImFCS control panel
    private int noSettings = 31;				// number of individual setting parameters 
    private String[] panelSettings = new String[noSettings];	// array to store the individual settings, the settings are stored in the same order as used to create the results table
    private boolean[] keyParam = {true, true, true, true, true, true, true, true, true, true, false, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, false};
    private int parcormode;
    // whether a change in the respective setting requires resetting of results arrays
    private boolean askOnRewrite = false; 		// whether setParameters should ask before Results arrays
    private boolean expload = false;			// is this an experiment load?
    private boolean checkroi = false;	//determines whether setParameters() checks for the correct setting of the ROI			
    private String bleachCorMem = "none";
    private String filterMem = "none";

    // variables to remember the last settings in Simulation Panel
    private int nosimsettings = 36;								// number of settings in the Simulation panel
    private String[] simSettings = new String[nosimsettings];	// array to store the last used settings

    // parameters for cursor position in window for correlations
    private int cposx;				// actual positions where the current correlation is calculated
    private int cposy;
    private int c2posx;				// actual positions of the second pixel with which the current correlation is calculated
    private int c2posy;
    private int maxcposx;			// max positions allowed in image depending on actual image width/height, binning, and distance between pixels;
    private int maxcposy;
    private int mincposx;			// min positions allowed in image depending on actual image width/height, binning, and distance between pixels;
    private int mincposy;
    private int pixelWidthX;		// determines the number of pixels that can be correlated, depending on whether overlap is allowed for the width and height of the image
    private int pixelHeightY;
    private int pixbinX;			// this is a multiplication factor to determine positions in the image window; it is 1 for overlap and equal to binning for non-overlap
    private int pixbinY;
    private int pixbin;
    private int roi1StartX = 0;		// first pixel of the ROI used for correlation
    private int roi1StartY = 0;
    private int roi1WidthX = 0;		// dimensions of the ROI used for correlations
    private int roi1HeightY = 0;
    private char keyMoveRight = '6';// control keys for cursor position used by this program
    private char keyMoveLeft = '4';	// they keys can be changed if required; they are then used by the key listener for the image window as defined below
    private char keyMoveUp = '8';	// the keys work on the number block of a PC but not on the general keyboard as some of the keys are already used by ImageJ
    private char keyMoveDown = '2';	// potentially this could be removed as the keys are not often used

    // arrays and parameters used for computation of correlations and storage of results
    private double[][] datac;				// temporary array which stores the values along two pixel columns to be correlated through the stack
    private double[][][][] acf; 			// ACF array to store all FCS functions for the pixels in the image; dimensions [ccf: ACF1, ACF2, CCF][width][height][chanum]
    private double[][][][] sdacf; 			// standerd deviation of the ACF; dimensions [ccf: ACF1, ACF2, CCF][width][height][chanum]
    private double[][][][] varacf; 			// standerd deviation of the ACF; dimensions [ccf: ACF1, ACF2, CCF][width][height][chanum]
    private double[][][][] msd; 			// MSD array to store all MSD plots for the pixels in the image; dimensions [ccf: MSD1, MSD2, MSD of CCF][width][height][chanum]
    private double[][] currentCovmats;		// the regularized covariance matrix is at the moment not stored but only updated for the current pixel if "GLS Fit" is selected
    private double aveacf[];				// Array to store average FCS function
    private double varaveacf[];				// Array to store variance of the average FCS function
    private double msdaveacf[];				// Array to store variance of the average FCS function
    private double[] lagtime; 				// lagtime array which stores the correlation times; dimensions [chanum]; defined in setParameters()
    private double[] mtab; 					// number of samples for each correlation channel; dimensions [chanum]; defined in setParameters()
    private double[] samp;					// sampletime (or bin width) of each channel; dimensions [chanum]; defined in setParameters()
    private double[] intTrace1;				// intensity traces
    private double[] intTrace2;
    private double[] intTime;				// time for intensity traces
    private double[][][][] fitacf;			// fit functions; [ccf: ACF1, ACF2, CCF][width][height][chanum]
    private double[][][][] fitres;			// fit parameter results; [ccf: ACF1, ACF2, CCF][width][height][noparam]
    private double[][][][] res;				// fit residuals; [ccf: ACF1, ACF2, CCF][width][height][chanum]
    private double[][][] chi2;				// chi2 values
    private double[][][] dccf;				// array to store dCCF values - in all possible directions of dCCF
    private double[][][] blocked;			// store yes/no whether blocking was succesful
    private double[][][] pixvalid;			// filtering mask, 1 if pixel valid, NaN if not
    private boolean[][][] pixfitted;		// whether pixel has been successfully fitted or not; in the absence of user-defined thresholds, this array determines pixvalid[][][]
    private double[][][] filterThresholds;	// thresholds for filtering parameter values; [ACF1, ACF2, CCF][parameter: noparam + Chi2][min, max] 
    private double[][] CCFq;				// map of cross-correlation amount
    private double fitaveacf[];				// Array to store fit of average FCS function
    private double resaveacf[];				// Array to store residuals of the average FCS function
    private double chi2aveacf;				// chi2 value for fit of average ACF function
    private RealVector transACF;			// transformed data L\y
    private RealMatrix lowerCholDecCovmats;	// lower diagonal matrix from the Cholesky decomposition of the regularized covariance matrix covamts
    private double[] transTheoreticalACF;	// temporary array with transformed theoretical model function L\y(x)
    private double[][] transTheoreticalGradientACF;	// temporary array with transformed theoretical model function L\y(x)

    private int noparam = 11;							// number of fit parameters; at the moment there are two fit models with a maximum of 11 parameters
    private String $param[] = {"N", "D", "vx", "vy", "G", "F2", "D2", "F3", "D3", "Ftrip", "Ttrip", "reduced Chi2", "blocked", "valid pixels", "q map"}; 	// parameter names; has to contain noparam number of parameter names
    private String $channel[] = {"", "g", "r"};			//names of correlation channels in DC-FCCS
    private final int histnum = 10;							// number of parameters for which histograms will be created;
    HistogramWindow histWin;							// define Histogram Windows
    private String $histWindowTitle;

    private double[] initparam;				// store the first fit parameters given by the user
    private double[] paraminitval;			// intitial values for the parameters used in the fit (can be the same as initparam or can be continously set; both options are given but only the first is used momentarily)
    private boolean[] paramfit;				// array with information whether a parameter is fit (true) or hold (false)
    private boolean[] paramfilter;			// array with information whether a filter is applied on the parameter
    private boolean[] userThreshold;		// settings for filtering: whether user has defined any thresholds, whether that was for DC-FCCS and whether the filter has been applied on current file + whether to use the same thresholds for ACFs and CCFs and whether to set q=0 for pixels where ACF valid and CCF invalid
    private int[] lag;						// lag for each correlation channel counted as multiples of the smallest basic time; independent of time; lagtime = lag * frametime
    private final double empty[] = {0.0};			// dummy array used for initilizing plot (should be removed later on)
    private double[][] difflaw = new double[3][1];	// diffusion law values for each pixel binning case for D and its standard deviation
    private double[][][][] difflawarray = new double[1][1][3][1];          // store all points for all subregions of an image for which the diffusion law was calculated 
    private final double[] difflawfit = {0.0, 0.0};	// store the intercept and the slope for the diffusion law
    private final int[] diffLawFitLim = {0, 0};		// store the start and end point of the diff. law fit
    private double[][][] diffLawFitMap = new double[1][1][2];	// store values if the diffusion law is applied to multiple areas in the image
    private int diffLawMapwidth = 1;                            // width of the diffusion law map
    private int diffLawMapheight = 1;                           // height of the diffusion law map
    private int difflawbin = 0;				// store the number of binnings that were calculated for the diffusion law
    private int difflawallbin = 1;				// store the number of binnings that were calculated for the diffusion law
    private int difflawmapbin = 1;				// store the number of binnings that were calculated for the diffusion law
    private double minvalDL;				// minimum value for difflawplot
    private double maxvalDL;				// maximum value for difflawplot
    private int DLbsmem = 0;				// remember input in DLwindow for bin calculations
    private int DLbemem = 0;				// this is used to check for errors
    private double[][][] psfData = new double[1][3][1];	// PSF data
    private int psfmaxbin = 1;				// maximum binning used for PSF calculation
    private int numofpsf = 1;				// this is defined by the user in psfDialogue(); numofPSF = (psfEnd - psfStart) / psfStep
    private float[][] filterArray;			// filter array; intensity filters to select the pixels to be correlated

    private int base; 					// base = number of channels in first group
    private int hbase; 				// hbase = number of channels in all higher groups
    private int chanum; 				// number of total channels of the correlator
    private int lagnum; 				// number of lag groups for the correlator
    private int blockIndS;				// Index whether block Transformation is succesful (1) or maximum blocking is used (0)
    private int nopit = 1;				// points in the shortened intensity traces for plotting
    private double maxsc; 				// scale for ACF plot
    private double minsc;
    private double imaxsc; 				// scale for intensity plot
    private double iminsc;
    private double psfStart;			// values for determining the PSF
    private double psfEnd;
    private double psfStep;
    private int psfBinStart;			// which binnning to be used
    private int psfBinEnd;
    private int filterLL = 0;			// intensity filter values with lower and upper limits LL and UL
    private int filterUL = 65536;
    private String currentmodel;		// remember the current model that is used for fitting; this is not necessarily the cbFitModel setting as in some cases multiple different fits can be used in a single run (e.g. for FCCS or Bayes fitting)

    // background image
    private double[][] bgrVar;			// variance calcualted from background file for each pixel
    private double[][] bgrCoVar;		// covariance calcualted from next neighbour frames of a background file
    private int bgrf;					// bgr number of frames
    private int bgrw;					// bgr width
    private int bgrh;					// bgr height
    private boolean bgrloaded = false;	// flag to indicate whether backgroudn image was loaded by user 
    private String $bgrTitleMem = " ";		// name of the background image file
    private double[][] bgrmean;			// mean values (pixel, row, coumn) for background file
    private double[] bgrrowmean;
    private double[] bgrcolumnmean;
    private double bgrframemean;

    // Number & Brightness
    private double NBslope;					// the slope of intensity vs variance used for correcting N&B
    private boolean NBperformed;			// whether N&B has been calculated
    private boolean NBcorrected;			// whether correction by the slope of intensity vs variance has been used
    private double[][] NBN;					// array to store the Number values in N&B
    private double[][] NBB;					// array to store the Brightness values in N&B
    private double[][] NBNum;				// array to store the corrected number values in N&B
    private double[][] NBEps;				// array to store the corected brightness (epsilon) values in N&B

    // default parameter values determinig Batch processing parameters; user can change them in the Batch Processing dialogue
    private boolean batchCorrelateAll = true;
    private boolean batchFit = false;
    private boolean batchPSF = false;
    private boolean batchDiffLaw = false;
    private boolean batchVerDCCF = false;
    private boolean batchHorDCCF = false;
    private boolean batchDiaUpDCCF = false;
    private boolean batchDiaDownDCCF = false;
    private String $batchSuffix;

    // Parameter for MSD calculation: false is 2D, true is 3D
    private boolean MSDmode = false;

    //Panel Variables for control panel "ImFCS"
    JFrame frame = new JFrame("ImFCS");		// items for ImFCS control panel
    private JTextField tfFirstFrame;		// a detailed description is given in the accompanying documentation
    private JTextField tfLastFrame;
    private JTextField tfFrameTime;
    private JTextField tfBinning;
    private JTextField tfCfXDistance;
    private JTextField tfCfYDistance;
    private JTextField tfBackground;
    private JTextField tfBackground2;
    private JTextField tfCorrelatorQ;
    private JTextField tfPixelSize;
    private JTextField tfMagnification;
    private JTextField tfNA;
    private JTextField tfEmLambda;
    private JTextField tfEmLambda2;
    private JTextField tfSigma;
    private JTextField tfSigmaZ;
    private JTextField tfSigma2;
    private JTextField tfSigmaZ2;
    private JComboBox<String> cbCorrelatorP;
    private JComboBox<String> cbBleachCor;
    private JComboBox<String> cbFilter;
    private JComboBox<String> cbParaCor;
    private JComboBox<String> cbDCCF;
    private boolean setImp = false;			// check whether an image was loaded or an existing one assigned
    private boolean doFit;					// should fits be performed?
    private boolean overlap;				// can pixels overlap if binning is used?
    private boolean doMSD;					// should MSD be performed alongside the ACF?
    private boolean doFiltering;			// should filtering be applied on values in paramter maps?
    private boolean batchSim = false;
    private File simBatchPath;				// directory where simulated files in batch mode can be saved
    private JButton btnBinning = new JButton("Binning");
    private JToggleButton tbFCCSDisplay = new JToggleButton("Off");
    private JToggleButton tbNB = new JToggleButton("N&B Off");
    private JToggleButton tbFit = new JToggleButton("Fit off");
    private JToggleButton tbOverlap = new JToggleButton("Off");
    private JToggleButton tbBGR = new JToggleButton("Bgr NUM");
    private JToggleButton tbMSD = new JToggleButton("MSD Off");
    private JToggleButton tbDL = new JToggleButton("Diff. Law");
    private JToggleButton tbSim = new JToggleButton("Sim off");
    private JToggleButton tbFiltering = new JToggleButton("Threshold");
    private JButton btnSave = new JButton("Save");
    private JButton btnLoad = new JButton("Read");
    private final JButton btnBatch = new JButton("Batch");
    private final JButton btnDCCF = new JButton("dCCF");

    // Parameters for IO
    String $imagePath;
    String $fileName;

    // ImFCS fit panel
    JFrame fitframe = new JFrame("ImFCS Fitting");
    private JComboBox<String> cbFitModel;
    private JTextField tfParama;
    private JTextField tfParamw;
    private JTextField tfParamw2;
    private JTextField tfParamz;
    private JTextField tfParamz2;
    private JTextField tfParamRx;
    private JTextField tfParamRy;
    private JTextField tfParamRz;
    private JTextField tfParamQ2;
    private JTextField tfParamN;
    private JTextField tfParamF2;
    private JTextField tfParamD;
    private JTextField tfParamD2;
    private JTextField tfParamF3;
    private JTextField tfParamQ3;
    private JTextField tfParamD3;
    private JTextField tfParamVx;
    private JTextField tfParamVy;
    private JTextField tfParamG;
    private JTextField tfParamFtrip;
    private JTextField tfParamTtrip;
    private JTextField tfFitStart;
    private JTextField tfFitEnd;
    private JTextField tfFitModel;
    private JTextField tfModProb1;
    private JTextField tfModProb2;
    private JTextField tfModProb3;
    private JToggleButton tbGLS = new JToggleButton("Off");
    private JToggleButton tbBayes = new JToggleButton("Off");
    private JToggleButton tbShow = new JToggleButton("Show", true);
    private final JButton btnTest = new JButton("Test");
    private JToggleButton tbFixPar = new JToggleButton("Free");
    private final JButton btnSetPar = new JButton("Set");
    private final JRadioButton rbtnHoldQ2 = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldN = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldF2 = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldD = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldD2 = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldF3 = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldD3 = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldQ3 = new JRadioButton("Hold");
    private JRadioButton rbtnHoldVx = new JRadioButton("Hold");
    private JRadioButton rbtnHoldVy = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldG = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldFtrip = new JRadioButton("Hold");
    private final JRadioButton rbtnHoldTtrip = new JRadioButton("Hold");

    // Parameter determining whether a simulation, from the program is used instead of a file
    private boolean simFile = false;
    private boolean simDomainFlag = false;
    private boolean simMeshFlag = false;
    private boolean simBlinkFlag = false;

    //Simulation frame
    JFrame simframe = new JFrame("Simulation Panel");
    private JComboBox<String> cbSimMode;
    private JTextField tfSimSeed;
    private JTextField tfSimParticleNum;
    private JTextField tfSimCPS;
    private JTextField tfSimTauBleach;
    private JTextField tfSimPixelNum;
    private JTextField tfSimExtensionFactor;
    private JTextField tfSimTimeStepNum;
    private JTextField tfSimFrameTime;
    private JTextField tfSimStepsPerFrame;
    private JTextField tfSimCurrentStepSize;
    private JTextField tfSimD1;
    private JTextField tfSimDoutDinRatio;
    private JTextField tfSimD2;
    private JTextField tfSimF2;
    private JTextField tfSimD3;
    private JTextField tfSimF3;
    private JTextField tfSimKon;
    private JTextField tfSimKoff;
    private JTextField tfSimCameraOffset;
    private JTextField tfSimCameraNoiseFactor;
    private JTextField tfSimBleachRadius;
    private JTextField tfSimBleachFrame;
    private JTextField tfDomainRadius;
    private JTextField tfDomainDensity;
    private JTextField tfPin;
    private JTextField tfPout;
    private JTextField tfMeshworkSize;
    private JTextField tfHopProbability;
    private JToggleButton tbSimTrip;
    private JButton btnSimulate;
    private JButton btnBatchSim;
    private JButton btnStopSimulation;
    private simulateACFWorker simulateACFInstant;
    private batchSimulateACFWorker batchSimulateACFInstant;

//    private int tmpPartNum;		//variable that identifies whether the coordinates of the simulations of a particle should be plotted in the log window, and if yes which particle that is
    private double batchDStart;
    private double batchDEnd;
    private double batchDStep;
    private double batchD2Start;
    private double batchD2End;
    private double batchD2Step;
    private double batchF2Start;
    private double batchF2End;
    private double batchF2Step;

    //Thresholds settings frame
    JFrame filteringframe;
    private JTextField tfAlFiltN;
    private JTextField tfAhFiltN;
    private JTextField tfAlFiltD;
    private JTextField tfAhFiltD;
    private JTextField tfAlFiltvx;
    private JTextField tfAhFiltvx;
    private JTextField tfAlFiltvy;
    private JTextField tfAhFiltvy;
    private JTextField tfAlFiltG;
    private JTextField tfAhFiltG;
    private JTextField tfAhFiltGrel;
    private JTextField tfAlFiltF2;
    private JTextField tfAhFiltF2;
    private JTextField tfAlFiltD2;
    private JTextField tfAhFiltD2;
    private JTextField tfAlFiltF3;
    private JTextField tfAhFiltF3;
    private JTextField tfAlFiltD3;
    private JTextField tfAhFiltD3;
    private JTextField tfAlFiltFtrip;
    private JTextField tfAhFiltFtrip;
    private JTextField tfAlFiltTtrip;
    private JTextField tfAhFiltTtrip;
    private JTextField tfAlFiltChi2;
    private JTextField tfAhFiltChi2;
    private JTextField tfClFiltN;
    private JTextField tfChFiltN;
    private JTextField tfClFiltD;
    private JTextField tfChFiltD;
    private JTextField tfClFiltvx;
    private JTextField tfChFiltvx;
    private JTextField tfClFiltvy;
    private JTextField tfChFiltvy;
    private JTextField tfClFiltG;
    private JTextField tfChFiltG;
    private JTextField tfChFiltGrel;
    private JTextField tfClFiltF2;
    private JTextField tfChFiltF2;
    private JTextField tfClFiltD2;
    private JTextField tfChFiltD2;
    private JTextField tfClFiltF3;
    private JTextField tfChFiltF3;
    private JTextField tfClFiltD3;
    private JTextField tfChFiltD3;
    private JTextField tfClFiltFtrip;
    private JTextField tfChFiltFtrip;
    private JTextField tfClFiltTtrip;
    private JTextField tfChFiltTtrip;
    private JTextField tfClFiltChi2;
    private JTextField tfChFiltChi2;
    private JRadioButton rbtnFiltN = new JRadioButton();
    private JRadioButton rbtnFiltD = new JRadioButton();
    private JRadioButton rbtnFiltvx = new JRadioButton();
    private JRadioButton rbtnFiltvy = new JRadioButton();
    private JRadioButton rbtnFiltG = new JRadioButton();
    private JRadioButton rbtnFiltGrel = new JRadioButton();
    private JRadioButton rbtnFiltF2 = new JRadioButton();
    private JRadioButton rbtnFiltD2 = new JRadioButton();
    private JRadioButton rbtnFiltF3 = new JRadioButton();
    private JRadioButton rbtnFiltD3 = new JRadioButton();
    private JRadioButton rbtnFiltFtrip = new JRadioButton();
    private JRadioButton rbtnFiltTtrip = new JRadioButton();
    private JRadioButton rbtnFiltChi2 = new JRadioButton();
    private JRadioButton rbtnCuseA = new JRadioButton();
    private JRadioButton rbtnReplaceZero = new JRadioButton();
    private JButton btnFilter = new JButton("Filter");
    private JButton btnReset = new JButton("Reset");

    //N&B frame
    JFrame NBframe = new JFrame("N&B");
    private JComboBox<String> cbNBMode;
    private JTextField tfNBS;
    private JTextField tfNBCalibRatio;
    private JButton btnNB;

    //Diffusion Law frame
    JFrame difflawframe = new JFrame("Diffusion Law Analysis");
    private JTextField tfDLBinStart;
    private JTextField tfDLBinEnd;
    private JTextField tfDLFitStart;
    private JTextField tfDLFitEnd;
    private JTextField tfDLROI;
    private JButton btnDiffLawCalculate;
    private JButton btnDiffLawFit;
    private JToggleButton tbDLROI;

    // ImFCS Results Table
    JFrame resframe = new JFrame("ImFCS Results");
    private final int resTabPosX = panelPosX + 1000;
    private final int resTabPosY = panelPosY;
    private final int resTabDimX = 400;
    private final int resTabDimY = 200;

    private boolean storepar = true;	// will fit parameters be re-used; the parmaeter can be changed in the program in the fit panel
    private int mlnum;					// remember the mouse and keylisteners added in the program
    private int klnum;					// this is used to correctly remove the listeners on exit  

    // main program; starts the FCS panel
    @Override
    public void run(String arg) {
        String tempmsg;
        try {
            if (GpufitImFCS.isCudaAvailable()) {
                isgpupresent = 1;
                tempmsg = "NVIDIA GPU is detected.";
            } else {
                tempmsg = GpufitImFCS.ALERT;
            }
        } catch (Exception e) {
            // NOTE: isgpupresent will not be set to 1. In calculateRoi function, calculations will be done on a CPU instead of GPU.
            tempmsg = "NVIDIA GPU is not detected.";
        }

        /* check screen size and adapt the dimensions of the panels accordingly
        if ( IJ.getScreenSize().getWidth() < 1900 ) {

        }

        if ( IJ.getScreenSize().getWidth() < 1400 ) {

        } */
        // See: https://www.geeksforgeeks.org/java-swing-jdialog-examples/
        // See: https://stackoverflow.com/questions/19274329/automatically-closing-jdialog-without-user-action
        // See: https://stackoverflow.com/questions/26075366/java-how-to-make-a-popup-window-close-automatically
        JFrame f = new JFrame("Loading ImFCS " + VERSION);
        f.setFocusable(true);
        f.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        f.setLayout(new GridLayout(2, 1));
        String modestr = (isgpupresent == 1) ? "GPU mode" : "CPU mode";
        JLabel l1 = new JLabel(modestr); // create a label 
        JLabel l2 = new JLabel("<html><p>" + tempmsg + "</p></html>"); // create a label. Use html tags to wrap long sentences if necessary.
        l1.setVerticalAlignment(JLabel.CENTER);
        l1.setHorizontalAlignment(JLabel.CENTER);
        l2.setVerticalAlignment(JLabel.CENTER);
        l2.setHorizontalAlignment(JLabel.CENTER);
        f.add(l1);
        f.add(l2);
        f.setSize(splashScreenDimX, splashScreenDimY);
        f.setLocation(panelPosX, (int) (panelPosY + panelDimY - splashScreenDimY) / 2);

        Timer timer = new Timer(splashScreenDuration, (ActionEvent e) -> {
            f.setVisible(false);
            f.dispose();

            // set Fonts
            setUIFont(panelFontSize, $panelFont);

            // start plugin and create panels
            SwingUtilities.invokeLater(() -> {
                createImFCSPanel();			// create ImFCS control panel and show it
                createImFCSFit();			// create ImFCS fit panel, but show only when fit is demanded in control panel
                createFiltering();			// create the JFrame for filtering panel and make it visible when "Thresholds" are switched on in the main panel
                createSimPanel();			// create Simulation panel, but show only when simulation is demanded in control panel
                createDiffLawPanel();		// create Diffusion Law panel, but show only when DiffLaw analysis is demanded in control panel
                createNBPanel();			// create NB panel, but show only when NB analysis is demanded in control panel
            });
        }
        );
        timer.setRepeats(false);
        timer.start();
        f.setVisible(true);
    }

    // Font Manager
    public static void setUIFont(int panelFontSize, String $panelFont) {
        UIManager.getLookAndFeelDefaults().put("defaultFont", new java.awt.Font($panelFont, java.awt.Font.PLAIN, panelFontSize));
        UIManager.put("Button.font", new java.awt.Font($panelFont, java.awt.Font.BOLD, panelFontSize));
        UIManager.put("ToggleButton.font", new java.awt.Font($panelFont, java.awt.Font.BOLD, panelFontSize));
        UIManager.put("RadioButton.font", new java.awt.Font($panelFont, java.awt.Font.BOLD, panelFontSize));
        UIManager.put("Label.font", new java.awt.Font($panelFont, java.awt.Font.ITALIC, panelFontSize));
        UIManager.put("ComboBox.font", new java.awt.Font($panelFont, java.awt.Font.PLAIN, panelFontSize));
        UIManager.put("TextField.font", new java.awt.Font($panelFont, java.awt.Font.PLAIN, panelFontSize));
        UIManager.put("ToolTip.font", new java.awt.Font($panelFont, java.awt.Font.PLAIN, panelFontSize));
    }

    /*
	 * CONTROLS FOR IMAGING FCS CONTROL PANEL AND IMAPGEPLUS
	 * 
	 * public void createImFCSPanel(): create a JFrame for the control panel
	 * 
	 * DocumetnListener:
	 * basicValueChanged()
	 * 
	 * MouseListener:
	 * impcanMouseUsed()
	 * 
	 * KeyListener:
	 * impcanKeyUsed(): keys 2, 4, 6, 8 = left right, up, down (definitions can be changed at start of program) 
	 * On PC use number block; on Mac numbers some overlapping usage of keys have nt yet been solved
	 * 
	 * ActionListeners:
	 * btnUseExistingPressed, btnLoadNewPressed, btnROIPressed, btnAllPressed, btnPSFPressed, btnDLPressed, btnParaCorPressed, tbFilteringPressed, 
	 * btnAvePressed, btnDCCFPressed, btnBtfPressed, btnBatchPressed, btnSavePressed, btnLoadPressed, btnExitPressed, btnRTPressed, 
	 * 
	 * ItemListeners:
	 * tbFCCSDisplayPressed, tbNBPressed, tbFitPressed, tbFitPressed, tbMSDPressed, tbOverlapPressed
	 * 
	 * ... and various generic dialogs required in the control panel
	 * fileAlreadyExistsDialogue()
	 * psfDialogue()
	 * filterDialogue()
	 * batchDialogue()
	 * SWSizeDialogue()
	 * PolynomialOrderDialogue()
	 * MSDDialogue()
	 * 
     */
    // create the ImFCS control panel
    public void createImFCSPanel() {

        frame.setFocusable(true);
        frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        frame.setLayout(new GridLayout(18, 4));
        frame.setLocation(new Point(panelPosX, panelPosY));
        frame.setSize(new Dimension(panelDimX, panelDimY));
        frame.setResizable(false);

        tfFirstFrame = new JTextField("1", 8);
        tfFrameTime = new JTextField("0.001", 8);
        tfLastFrame = new JTextField(Integer.toString(frames), 8);
        tfBinning = new JTextField("1 x 1", 8);
        tfBackground = new JTextField(Integer.toString(impmin), 8);
        tfBackground2 = new JTextField(Integer.toString(impmin), 8);
        tfCfXDistance = new JTextField("0", 8);
        tfCfYDistance = new JTextField("0", 8);
        tfPixelSize = new JTextField($UPixelSize, 8);
        tfMagnification = new JTextField($UMagnification, 8);
        tfNA = new JTextField($UNA, 8);
        tfEmLambda = new JTextField($UEmWavelength1, 8);
        tfEmLambda2 = new JTextField($UEmWavelength2, 8);
        tfSigma = new JTextField($ULateralPSF1, 8);
        tfSigmaZ = new JTextField($UAxialPSF1, 8);
        tfSigma2 = new JTextField($ULateralPSF2, 8);
        tfSigmaZ2 = new JTextField($UAxialPSF2, 8);

        cbFitModel = new JComboBox<String>();
        cbFitModel.addItem("FCS");
        cbFitModel.addItem("DC-FCCS");

        cbCorrelatorP = new JComboBox<String>();
        cbCorrelatorP.addItem("16");
        cbCorrelatorP.addItem("32");

        tfCorrelatorQ = new JTextField($CorrelQ, 8);

        cbBleachCor = new JComboBox<String>();
        cbBleachCor.addItem("none");
        cbBleachCor.addItem("Sliding Window");
        cbBleachCor.addItem("Single Exp");
        cbBleachCor.addItem("Double Exp");
        cbBleachCor.addItem("Polynomial");
        cbBleachCor.addItem("Lin Segment");

        cbFilter = new JComboBox<>();
        cbFilter.addItem("none");
        cbFilter.addItem("Intensity");
        cbFilter.addItem("Mean");

        cbParaCor = new JComboBox<>();
        cbParaCor.addItem("N vs D");
        cbParaCor.addItem("N vs F2");
        cbParaCor.addItem("D vs F2");
        cbParaCor.addItem("N*(1-F2) vs D");
        cbParaCor.addItem("N*F2 vs D2");
        cbParaCor.addItem("D vs Sqrt(vx^2+vy^2)");
        cbParaCor.addItem("D2 vs Sqrt(vx^2+vy^2)");

        cbDCCF = new JComboBox<>();
        cbDCCF.addItem("x direction");
        cbDCCF.addItem("y direction");
        cbDCCF.addItem("diagonal /");
        cbDCCF.addItem("diagonal \\");

        tfFrameTime.setToolTipText("Time per frame. NOTE: Changing this value will reinitialize all arrays.");
        tfBinning.setToolTipText("Pixel binning used in the evlauations. NOTE: Changing this value will reinitialize all arrays.");
        tfCfXDistance.setToolTipText("Distance in x-direction for spatial cross-correlation. NOTE: Changing this value will reinitialize all arrays.");
        tfCfYDistance.setToolTipText("Distance in y-direction for spatial cross-correlation. NOTE: Changing this value will reinitialize all arrays.");

        JButton btnUseExisting = new JButton("Use");
        btnUseExisting.setToolTipText("Uses the active existing image in ImageJ.");
        JButton btnLoadNew = new JButton("Load");
        btnLoadNew.setToolTipText("Opens a dialog to open a new image.");
        JButton btnROI = new JButton("ROI");
        btnROI.setToolTipText("Calculates ACFs only in the currently chose ROI.");
        JButton btnAll = new JButton("All");
        btnAll.setToolTipText("Calculates all ACFs.");
        JButton btnPSF = new JButton("PSF");
        btnPSF.setToolTipText("Calculates the calibration for the PSF.");
        JButton btnOptions = new JButton("Options");
        btnOptions.setToolTipText("Select various options regarding the display of results.");
        JButton btnParaCor = new JButton("Scatter");
        btnParaCor.setToolTipText("Calculates a scatter plot for a pair of two parameters from the scroll down menue.");
        JButton btnBtf = new JButton("To Front");
        btnBtf.setToolTipText("Bring all windows of this plugin instance to the front.");
        JButton btnExit = new JButton("Exit");
        JButton btnRT = new JButton("Res. Table");
        btnRT.setToolTipText("Create a results table.");
        JButton btnAve = new JButton("Average");
        btnAve.setToolTipText("Calculate the average ACF from all valid ACFs and fit if fit is switched on; this does not calculate residuals or sd.");

        btnDCCF.setToolTipText("Create a dCCF image to see differences between forward and backward correlation in a direction (see scroll down menue).");
        tbBGR.setToolTipText("Load background file for correction. Has to be the same area recorded under the same conditions as the experimental file");
        tbSim.setToolTipText("Opens/closes Simulation panel.");
        tbDL.setToolTipText("Calculates the Diffusion Law.");
        tbFit.setToolTipText("Switches Fit on/off; opens/closes Fit panel.");
        tbMSD.setToolTipText("Switches Mean Square Displacement calculation and plot on/off.");
        tbFiltering.setToolTipText("Filters the values in parameters maps using user-defined thresholds");
        btnBatch.setToolTipText("Allow to select a list of evaluations to be performed on a range of images.");
        btnSave.setToolTipText("Save the evaluation of the data as binary files. Which data to save can be selected in a dialog.");
        btnLoad.setToolTipText("Load a previously saved experiment. Note that the original image is not automatically loaded along.");

        //row 1
        frame.add(new JLabel("Image"));
        frame.add(btnUseExisting);
        frame.add(btnLoadNew);
        frame.add(btnBatch);

        //row 2
        frame.add(new JLabel("First frame: "));
        frame.add(tfFirstFrame);
        frame.add(new JLabel("Last frame: "));
        frame.add(tfLastFrame);

        //row 3
        frame.add(new JLabel("Frame time: "));
        frame.add(tfFrameTime);
        frame.add(btnBinning);
        frame.add(tfBinning);

        //row 4
        frame.add(new JLabel("CF X distance: "));
        frame.add(tfCfXDistance);
        frame.add(new JLabel("CF Y distance: "));
        frame.add(tfCfYDistance);

        //row 5
        frame.add(new JLabel("Correlator P: "));
        frame.add(cbCorrelatorP);
        frame.add(new JLabel("Correlator Q: "));
        frame.add(tfCorrelatorQ);

        //row 6
        frame.add(new JLabel("Fit Model: "));
        frame.add(cbFitModel);
        frame.add(new JLabel("FCCS display"));
        frame.add(tbFCCSDisplay);

        //row 7
        frame.add(new JLabel("Pixel size [um]:"));
        frame.add(tfPixelSize);
        frame.add(new JLabel("Overlap"));
        frame.add(tbOverlap);

        //row 8
        frame.add(new JLabel("Magnification:"));
        frame.add(tfMagnification);
        frame.add(new JLabel("NA"));
        frame.add(tfNA);

        //row 9
        frame.add(new JLabel("\u03bb\u2081 [nm]:"));
        frame.add(tfEmLambda);
        frame.add(new JLabel("\u03bb\u2082 [nm]:"));
        frame.add(tfEmLambda2);

        //row 10
        frame.add(new JLabel("PSF (xy):"));
        frame.add(tfSigma);
        frame.add(new JLabel("PSF (z):"));
        frame.add(tfSigmaZ);

        //row 11
        frame.add(new JLabel("PSF2 (xy):"));
        frame.add(tfSigma2);
        frame.add(new JLabel("PSF2 (z):"));
        frame.add(tfSigmaZ2);

        //row 12
        frame.add(tbBGR);
        frame.add(tfBackground);
        frame.add(new JLabel("Background 2: "));
        frame.add(tfBackground2);

        //row 13
        frame.add(btnOptions);
        frame.add(tbNB);
        frame.add(tbFiltering);
        frame.add(btnAve);

        //row 14
        frame.add(btnParaCor);
        frame.add(cbParaCor);
        frame.add(new JLabel("Bleach Cor."));
        frame.add(cbBleachCor);

        //row 15
        frame.add(btnDCCF);
        frame.add(cbDCCF);
        frame.add(new JLabel("Filter (All):"));
        frame.add(cbFilter);

        //row 16
        frame.add(btnPSF);
        frame.add(tbDL);
        frame.add(tbFit);
        frame.add(btnAll);

        //row 17
        frame.add(tbSim);
        frame.add(btnRT);
        frame.add(tbMSD);
        frame.add(btnROI);

        //row 18
        frame.add(btnBtf);
        frame.add(btnSave);
        btnSave.setForeground(java.awt.Color.BLUE);
        frame.add(btnLoad);
        btnLoad.setForeground(java.awt.Color.BLUE);
        frame.add(btnExit);
        btnExit.setForeground(java.awt.Color.RED);

        // add listeners
        btnUseExisting.addActionListener(btnUseExistingPressed);
        btnLoadNew.addActionListener(btnLoadNewPressed);
        btnBinning.addActionListener(btnBinningPressed);
        tbBGR.addActionListener(tbBGRPressed);
        btnROI.addActionListener(btnROIPressed);
        btnAll.addActionListener(btnAllPressed);
        btnPSF.addActionListener(btnPSFPressed);
        btnOptions.addActionListener(btnOptionsPressed);
        btnParaCor.addActionListener(btnParaCorPressed);
        btnAve.addActionListener(btnAvePressed);
        btnDCCF.addActionListener(btnDCCFPressed);
        btnBtf.addActionListener(btnBtfPressed);
        btnBatch.addActionListener(btnBatchPressed);
        btnSave.addActionListener(btnSavePressed);
        btnLoad.addActionListener(btnLoadPressed);
        btnExit.addActionListener(btnExitPressed);
        btnRT.addActionListener(btnRTPressed);
        tbDL.addItemListener(tbDLPressed);
        tbFCCSDisplay.addItemListener(tbFCCSDisplayPressed);
        tbNB.addItemListener(tbNBPressed);
        tbFit.addItemListener(tbFitPressed);
        tbOverlap.addItemListener(tbOverlapPressed);
        tbFiltering.addActionListener(tbFilteringPressed);
        tbSim.addItemListener(tbSimPressed);
        tbMSD.addActionListener(tbMSDPressed);
        cbBleachCor.addActionListener(cbBleachCorChanged);
        cbFilter.addActionListener(cbFilterChanged);

        frame.setVisible(true);
    }

    // mouse listener: control mouse click for selection of pixels and subsequent correlation and fit (if fit is on)
    // also checks that parameters were set and that they make sense
    MouseListener impcanMouseUsed = new MouseListener() {
        @Override
        public void mouseClicked(MouseEvent e) {
            int px = e.getX();
            int py = e.getY();

            checkroi = false;

            // read and set parameters form the control panel; this can only be done if image was loaded
            if (setParameters()) {

                // get the mouse coordinates
                initposx = (int) Math.floor(impcan.offScreenX(px) / pixbinX);
                initposy = (int) Math.floor(impcan.offScreenY(py) / pixbinY);

                performCFE(initposx, initposy);
            }
        }

        @Override
        public void mousePressed(MouseEvent e) {
        }	// other mouse events have no action associated yet

        @Override
        public void mouseReleased(MouseEvent e) {
            Roi improi = imp.getRoi();
            if (improi != null) {
                Rectangle rect = improi.getBounds();
                roi1StartX = (int) rect.getX();
                roi1StartY = (int) rect.getY();
                roi1WidthX = (int) rect.getWidth();
                roi1HeightY = (int) rect.getHeight();
            }
            checkroi = true;
        }

        @Override
        public void mouseEntered(MouseEvent e) {
        }

        @Override
        public void mouseExited(MouseEvent e) {
        }
    };

    // key listener for selection of pixels and subsequent correlation and fit (if fit is on)
    // also checks that parameters were set and that they make sense
    // works well on PC but is not working on MAC as key combination relies on number block
    KeyListener impcanKeyUsed = new KeyListener() {
        @Override
        public void keyTyped(KeyEvent e) {
            char keyChar = e.getKeyChar();

            // read and set parameters form the control panel; this can only be done if image was loaded
            if (setParameters()) {

                if (keyChar == keyMoveRight) {
                    cposx += pixbinX;
                    if (cposx >= width) {
                        cposx = width - 1;
                    }
                }

                if (keyChar == keyMoveLeft) {
                    cposx -= pixbinX;
                    if (cposx < 0) {
                        cposx = 0;
                    }
                }

                if (keyChar == keyMoveUp) {
                    cposy -= pixbinY;
                    if (cposy < 0) {
                        cposy = 0;
                    }
                }

                if (keyChar == keyMoveDown) {
                    cposy += pixbinY;
                    if (cposy >= height) {
                        cposy = height - 1;
                    }
                }

                initposx = (int) Math.floor(cposx / pixbinX);
                initposy = (int) Math.floor(cposy / pixbinY);

                performCFE(initposx, initposy);
            }
        }

        @Override
        public void keyPressed(KeyEvent e) {
        }

        @Override
        public void keyReleased(KeyEvent e) {
        }
    };

    // Action Listeners for the buttons on the ImFCS panel
    ActionListener btnExitPressed = (ActionEvent ev) -> {
        exitImFCS();
    };

    ActionListener btnBtfPressed = (ActionEvent ev) -> {
        bringToFront();
    };

    ActionListener btnUseExistingPressed = (ActionEvent ev) -> {
        if (WindowManager.getImageCount() > 0) {
            closeWindows();
            imp = IJ.getImage();
            obtainImage();
        } else {
            IJ.showMessage("No image open.");
        }
    };

    ActionListener btnLoadNewPressed = (ActionEvent ev) -> {
        imp = IJ.openImage();
        if (imp != null) {
            imp.show();
            obtainImage();
            closeWindows();
        }
    };

    ActionListener btnBinningPressed = (ActionEvent ev) -> {
        GenericDialog gd = new GenericDialog("Binning");
        gd.addNumericField("binning x: ", binningX, 0);
        gd.addNumericField("binning y: ", binningY, 0);
        gd.showDialog();
        if (gd.wasOKed()) {
            binningX = (int) gd.getNextNumber();
            binningY = (int) gd.getNextNumber();
            tfBinning.setText(Integer.toString(binningX) + " x " + Integer.toString(binningY));
        }
    };

    ActionListener tbBGRPressed = (ActionEvent ev) -> {
        if (tbBGR.getText().equals("Bgr FILE")) {
            tbBGR.setText("Bgr NUM");
            bgrloaded = false;
        } else {
            tbBGR.setText("Bgr FILE");
            if (loadBGRFile()) {
                bgrloaded = true;
            } else {
                tbBGR.setText("Bgr NUM");
                tbBGR.setSelected(false);
                bgrloaded = false;
            }
        }
    };

    ActionListener btnSavePressed = (ActionEvent ev) -> {
        if (setImp) {
            String xlsxFN = $impTitle;
            int dotind = xlsxFN.lastIndexOf('.');
            if (dotind != -1) {
                xlsxFN = xlsxFN.substring(0, dotind);
            }
            xlsxFN += ".xlsx";
            JFileChooser fc = new JFileChooser($imagePath);
            fc.setSelectedFile(new File(xlsxFN));
            fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
            int returnVal = fc.showSaveDialog(btnSave);
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                writeExperiment(fc.getSelectedFile(), "Failed to write data files", true);
            }
        } else {
            JOptionPane.showMessageDialog(null, "No Image loaded or assigned.");
        }
    };

    ActionListener btnLoadPressed = (ActionEvent ev) -> {
        if (setImp) {
            JFileChooser fc = new JFileChooser($imagePath);
            fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
            int returnVal = fc.showOpenDialog(btnLoad);
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                readExperiment(fc.getSelectedFile(), "Failed to read data files");
            }
        } else {
            JOptionPane.showMessageDialog(null, "No Image loaded or assigned.");
        }
    };

    ActionListener btnBatchPressed = (ActionEvent ev) -> {
        setImp = true;
        if (lastframe < firstframe) {
            firstframe = 1;
            lastframe = 2;
            tfFirstFrame.setText(Integer.toString(firstframe));
            tfLastFrame.setText(Integer.toString(lastframe));
        }
        batchWorker batchInstant = new batchWorker();
        batchInstant.execute();
    };

    ActionListener btnPSFPressed = (ActionEvent ev) -> {
        if (setParameters()) {
            correlatePSFWorker correlatePSFInstant = new correlatePSFWorker();
            correlatePSFInstant.execute();
        }
    };

    ItemListener tbDLPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ev.SELECTED) {
            if (setParameters()) {
                setDLparameters();
                tbDL.setSelected(true);
                difflawframe.setVisible(true);
            }
        } else {
            difflawframe.setVisible(false);
        }
    };

    ActionListener btnOptionsPressed = (ActionEvent ev) -> {
        GenericDialog gd = new GenericDialog("Options");
        gd.addCheckbox("ACF", plotACFCurves);
        gd.addCheckbox("SD", plotSDCurves);
        gd.addCheckbox("Intensity", plotIntensityCurves);
        gd.addCheckbox("Residuals", plotResCurves);
        gd.addCheckbox("Histogram", plotParaHist);
        gd.addCheckbox("Blocking", plotBlockingCurve);
        gd.addCheckbox("Covariance Matrix", plotCovmats);
        gd.hideCancelButton();
        gd.showDialog();
        if (gd.wasOKed()) {
            plotACFCurves = gd.getNextBoolean();
            plotSDCurves = gd.getNextBoolean();
            plotIntensityCurves = gd.getNextBoolean();
            plotResCurves = gd.getNextBoolean();
            plotParaHist = gd.getNextBoolean();
            plotBlockingCurve = gd.getNextBoolean();
            plotCovmats = gd.getNextBoolean();
            if (!plotACFCurves && (plotSDCurves || plotResCurves)) {
                plotSDCurves = false;
                plotResCurves = false;
                IJ.showMessage("Plotting of SD and/or Residuals without the ACF is not supported.");
            }
            if (!plotACFCurves) {
                if (acfWindow != null && acfWindow.isClosed() == false) {	// close ACF window
                    acfWindow.close();
                }
            }
            if (!plotSDCurves) {
                if (sdWindow != null && sdWindow.isClosed() == false) {	// close SD window
                    sdWindow.close();
                }
            }
            if (!plotIntensityCurves) {
                if (intWindow != null && intWindow.isClosed() == false) {	// close intensity trace window
                    intWindow.close();
                }
            }
            if (!plotResCurves) {
                if (resWindow != null && resWindow.isClosed() == false) {	// close fit residuals window
                    resWindow.close();
                }
            }
            if (!plotParaHist) {
                if (histWin != null && histWin.isClosed() == false) {
                    histWin.close();
                }
            }
            if (!plotBlockingCurve) {
                if (blockingWindow != null && blockingWindow.isClosed() == false) {	// close blocking window
                    blockingWindow.close();
                }
            }
            if (!plotCovmats) {
                if (impCovWin != null && impCovWin.isClosed() == false) {	// close covariance window
                    impCovWin.close();
                }
            }
        }
    };

    ItemListener tbNBPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tbNB.setText("N&B On");
            NBframe.setVisible(true);
        } else {
            tbNB.setText("N&B Off");
            NBframe.setVisible(false);
        }
    };

    ActionListener tbFilteringPressed = (ActionEvent ev) -> {
        if (tbFiltering.isSelected()) {
            doFiltering = true;
            tbFiltering.setSelected(true);
            SetThresholds();
            filteringframe.setVisible(true);
        } else {
            doFiltering = false;
            tbFiltering.setSelected(false);
            if (filteringframe != null) {
                filteringframe.setVisible(false);
            }
        }
    };

    ActionListener btnAvePressed = (ActionEvent ev) -> {
        calcAveCF();
    };

    ActionListener btnParaCorPressed = (ActionEvent ev) -> {
        if (setParameters()) {
            plotScatterPlot();
        }
    };

    ActionListener btnDCCFPressed = (ActionEvent ev) -> {
        if (setParameters()) {
            String $dccf = (String) cbDCCF.getSelectedItem();
            int mode;
            int wx = pixelWidthX;
            int hy = pixelHeightY;

            if ("y direction".equals($dccf)) {
                mode = 2;
            } else if ("diagonal /".equals($dccf)) {
                mode = 3;
            } else if ("diagonal \\".equals($dccf)) {
                mode = 4;
            } else {
                mode = 1;
            }
            correlateDCCFWorker correlateDCCFInstant = new correlateDCCFWorker(mode, wx, hy);
            correlateDCCFInstant.execute();
        }
    };

    ActionListener btnRTPressed = (ActionEvent ev) -> {
        createImFCSResultsTable();
    };

    ActionListener btnROIPressed = (ActionEvent ev) -> {
        if (setParameters()) {
            if (impPara1 != null) {		// set the parameter window to the front to avaid a "stack required" error from ImageJ
                WindowManager.setCurrentWindow(impPara1Win);
            }
            Roi improi = imp.getRoi();
            Rectangle rect = improi.getBounds();
            roi1StartX = (int) rect.getX();
            roi1StartY = (int) rect.getY();
            roi1WidthX = (int) rect.getWidth();
            roi1HeightY = (int) rect.getHeight();

            if (improi != null) {
                if (imp.getOverlay() != null) {
                    imp.getOverlay().clear();
                    imp.setOverlay(imp.getOverlay());
                }
                if (cfXDistance != 0 || cfYDistance != 0) {
                    Roi impRoiCCF = (Roi) improi.clone();
                    improi.setLocation(roi1StartX, roi1StartY);
                    improi.setStrokeColor(java.awt.Color.GREEN);
                    imp.setRoi(improi);
                    impRoiCCF.setLocation(roi1StartX + cfXDistance, roi1StartY + cfYDistance);
                    impRoiCCF.setStrokeColor(java.awt.Color.RED);
                    Overlay ccfov = new Overlay(impRoiCCF);
                    imp.setOverlay(ccfov);
                }
                correlateRoiWorker correlateRoiInstant = new correlateRoiWorker(improi);
                correlateRoiInstant.execute();
            } else {
                JOptionPane.showMessageDialog(null, "No ROI chosen.");
            }
        }
    };

    ActionListener btnAllPressed = (ActionEvent ev) -> {
        if (setParameters()) {
            if (imp.getOverlay() != null) {
                imp.getOverlay().clear();
                imp.setOverlay(imp.getOverlay());
            }

            int roi2StartX;
            int roi2WidthX;
            int roi2StartY;
            int roi2HeightY;
            if (cfXDistance > 0) {
                roi1StartX = 0;
                roi2StartX = cfXDistance;
            } else {
                roi1StartX = -cfXDistance;
                roi2StartX = 0;
            }
            if (cfYDistance > 0) {
                roi1StartY = 0;
                roi2StartY = cfYDistance;
            } else {
                roi1StartY = -cfYDistance;
                roi2StartY = 0;
            }
            if (overlap) {
                roi1WidthX = width - Math.abs(cfXDistance);
                roi1HeightY = height - Math.abs(cfYDistance);
            } else {
                roi1WidthX = (int) Math.floor((width - Math.abs(cfXDistance)) / binningX) * binningX;
                roi1HeightY = (int) Math.floor((height - Math.abs(cfYDistance)) / binningY) * binningY;
            }
            roi2WidthX = roi1WidthX;
            roi2HeightY = roi1HeightY;
            Roi impRoi1 = new Roi(roi1StartX, roi1StartY, roi1WidthX, roi1HeightY);
            impRoi1.setStrokeColor(java.awt.Color.GREEN);
            imp.setRoi(impRoi1);
            Roi impRoi2 = new Roi(roi2StartX, roi2StartY, roi2WidthX, roi2HeightY);
            if (cfXDistance != 0 || cfYDistance != 0) {
                impRoi2.setStrokeColor(java.awt.Color.RED);
                Overlay cfov = new Overlay(impRoi2);
                imp.setOverlay(cfov);
            }

            checkroi = true;
            if (setParameters()) {
                correlateRoiWorker correlateRoiInstant = new correlateRoiWorker(impRoi1);
                correlateRoiInstant.execute();
            } else {
                imp.getOverlay().clear();
                imp.deleteRoi();
            }
        }
    };

    ItemListener tbFCCSDisplayPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED && cbFitModel.getSelectedItem() == "DC-FCCS") {
            tbFCCSDisplay.setText("On");
        } else {
            tbFCCSDisplay.setSelected(false);
            tbFCCSDisplay.setText("Off");
        }
    };

    ItemListener tbFitPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            doFit = true;
            tbFit.setText("Fit On");
            fitframe.setVisible(true);
        } else {
            doFit = false;
            tbFit.setText("Fit Off");
            fitframe.setVisible(false);
        }
    };

    ItemListener tbSimPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tbSim.setText("Sim On");
            simframe.setVisible(true);
        } else {
            tbSim.setText("Sim Off");
            simframe.setVisible(false);
        }
    };

    ActionListener tbMSDPressed = (ActionEvent ev) -> {
        if (tbMSD.isSelected()) {
            doMSD = true;
            tbMSD.setText("MSD On");
            tbMSD.setSelected(true);
            MSDDialogue();
        } else {
            doMSD = false;
            tbMSD.setText("MSD Off");
            tbMSD.setSelected(false);
        }
    };

    ItemListener tbOverlapPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            overlap = true;
            tbOverlap.setText("On");
        } else {
            overlap = false;
            tbOverlap.setText("Off");
        }
    };

    // if in saveExperiments the chosen filename already exists ask wjhether to replace or stop
    public boolean fileAlreadyExistsDialog() {
        GenericDialog gd = new GenericDialog("File already exists.");
        gd.addMessage("A file with that name already exists.");
        gd.enableYesNoCancel("Replace", "Cancel");
        gd.hideCancelButton();
        gd.showDialog();
        return gd.wasOKed();
    }

    // generic dialogue to read in start and end values and step size for the PSF calculations
    public boolean psfDialogue() {
        GenericDialog gd = new GenericDialog("PSF");
        psfStart = 0.6;
        psfEnd = 1.0;
        psfStep = 0.1;
        psfBinStart = 1;
        psfBinEnd = 5;
        gd.addNumericField("start values: ", psfStart, 2);
        gd.addNumericField("end value: ", psfEnd, 2);
        gd.addNumericField("step size: ", psfStep, 2);
        gd.addStringField("bin start value: ", Integer.toString(psfBinStart), 2);
        gd.addStringField("bin end value: ", Integer.toString(psfBinEnd), 2);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        psfStart = gd.getNextNumber();
        psfEnd = gd.getNextNumber();
        psfStep = gd.getNextNumber();
        try {
            psfBinStart = Integer.parseInt(gd.getNextString());
            psfBinEnd = Integer.parseInt(gd.getNextString());
        } catch (NumberFormatException nfe) {
            IJ.showMessage("psfBinStart and psfBinEnd have to be Integer.");
            throw new NumberFormatException("Number format error.");
        }

        if (!(psfStart > 0) || !(psfEnd > 0) || !(psfStep > 0) || !(psfBinStart > 0) || !(psfBinEnd > 0) || (psfBinEnd > width) || (psfBinEnd > height)) {
            IJ.showMessage("only positive numbers are accepted");
            return false;
        }
        return true;
    }

    // generic dialogue to read lower and upper limits for the filter functions to determine which pixels should be correlated
    public boolean filterDialogue() {
        double LL;
        double UL;
        GenericDialog gd = new GenericDialog("Filter");
        gd.addMessage("non-negative integers only");
        gd.addNumericField("lower limt: ", filterLL, 2);
        gd.addNumericField("upper limit : ", filterUL, 2);
        gd.showDialog();
        if (gd.wasCanceled()) {
            cbFilter.setSelectedItem(filterMem);
            return false;
        }
        LL = gd.getNextNumber();
        UL = gd.getNextNumber();
        if (!(LL >= 0) || !(UL >= 0)) {
            IJ.showMessage("Illegal filter limits");
            cbFilter.setSelectedItem(filterMem);
            return false;
        } else {
            filterLL = (int) Math.floor(LL);
            filterUL = (int) Math.floor(UL);
        }
        return true;
    }

    // Batch dialogue asking for which evaluations should be done
    public boolean batchDialogue() {
        GenericDialog gd = new GenericDialog("Batch Processing");

        gd.addCheckbox("Correlate All", batchCorrelateAll);
        gd.addCheckbox("Fit", batchFit);
        gd.addCheckbox("PSF calculation", batchPSF);
        gd.addCheckbox("Diffusion Law", batchDiffLaw);
        gd.addCheckbox("verticalDCCF", batchVerDCCF);
        gd.addCheckbox("horizontalDCCF", batchHorDCCF);
        gd.addCheckbox("diaUpDCCF", batchDiaUpDCCF);
        gd.addCheckbox("diaDownDCCF", batchDiaDownDCCF);
        gd.addStringField("File suffix", "");
        gd.addMessage("If empty the date will be used as suffix.");
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }

        // Filter and bleach correction and other settings are automatically taken from FCS panel
        // Fit parameters and fit types are automatically taken from Fit Panel 
        batchCorrelateAll = gd.getNextBoolean();
        batchFit = gd.getNextBoolean();
        batchPSF = gd.getNextBoolean();
        batchDiffLaw = gd.getNextBoolean();
        batchVerDCCF = gd.getNextBoolean();
        batchHorDCCF = gd.getNextBoolean();
        batchDiaUpDCCF = gd.getNextBoolean();
        batchDiaDownDCCF = gd.getNextBoolean();
        $batchSuffix = gd.getNextString();
        return true;
    }

    // dialog asking for the size of the sliding window
    public boolean SWSizeDialogue() {
        GenericDialog gd = new GenericDialog("Window Size");
        gd.addMessage("Only Integers allowed.");
        gd.addNumericField("Size: ", slidingWindowLength, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            cbBleachCor.setSelectedItem(bleachCorMem);
            return false;
        }
        int swtmp = (int) Math.floor(gd.getNextNumber());
        if (swtmp > (lastframe - firstframe) || !(swtmp > 0)) {
            IJ.showMessage("Illegal window size");
            return false;
        } else {
            slidingWindowLength = swtmp;
        }
        return true;
    }

    // Dialog asking for the polynomial order and limiting it to BCmaxorder,  
    // which is defined in the user variables at the start of the plugin
    public boolean PolynomialOrderDialogue() {
        GenericDialog gd = new GenericDialog("Polynomial Order");
        gd.addMessage("Only Integers allowed.");
        gd.addNumericField("Order: ", polyOrder, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            cbBleachCor.setSelectedItem(bleachCorMem);
            return false;
        }
        int potmp = (int) Math.floor(gd.getNextNumber());
        if (!(potmp > 0)) {
            IJ.showMessage("Invalid polynomial order");
            return false;
        }
        if (potmp > BCmaxorder) {
            IJ.showMessage("Order > " + BCmaxorder + "not allowed");
            return false;
        } else {
            polyOrder = potmp;
        }
        return true;
    }

    // Dialog settign the MSD mode (2D or 3D)  
    public boolean MSDDialogue() {
        GenericDialog gd = new GenericDialog("MSD Model");
        gd.addMessage("2D is default (ITIR-FCS)");
        gd.addCheckbox("3D (SPIM-FCS)", MSDmode);
        gd.showDialog();
        if (gd.wasCanceled()) {
            doMSD = false;
            tbMSD.setText("MSD Off");
            tbMSD.setSelected(false);
            return false;
        } else {
            MSDmode = gd.getNextBoolean();
            return true;
        }
    }

    /* 
	 * FCS FIT PANEL
	 * 
	 * public void createImFCSFit(): create the JFrame for fit panel and make it visible if "fit" is switched on in the main panel
	 * 
	 * followed by several listeners:
	 * cbModelChanged()
	 * tbGLSPressed()
	 * tbBayesPressed()
	 * tbShowPressed()
	 * btnTestPressed()
	 * tbFixParPressed()
	 * btnSetParPressed()
	 * cbBleachCorChanged()
	 * cbFilterChanged()
	 * 
     */
    // create the FCS fit panel
    public void createImFCSFit() {

        fitframe.setFocusable(true);
        fitframe.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        fitframe.setLayout(new GridLayout(22, 6));
        fitframe.setLocation(new Point(fitPanelPosX, fitPanelPosY));
        fitframe.setSize(new Dimension(fitPanelDimX, fitPanelDimY));
        fitframe.setResizable(false);

        String $pixelsize = tfPixelSize.getText();
        pixelsize = Double.parseDouble($pixelsize);
        String $objmag = tfMagnification.getText();
        objmag = Double.parseDouble($objmag);
        String strbin = tfBinning.getText();
        int strlen = strbin.length();
        binningX = Integer.parseInt(strbin.substring(0, strbin.indexOf(" x ")));
        binningY = Integer.parseInt(strbin.substring(strbin.indexOf(" x ") + 3), strlen - 1);
        tfParama = new JTextField(Double.toString(pixelsize * 1000 / objmag * binningX), 8);		//XXX

        String $NA = tfNA.getText();
        NA = Double.parseDouble($NA);
        String $sigma = tfSigma.getText();
        sigma = Double.parseDouble($sigma);
        String $emlambda = tfEmLambda.getText();
        emlambda = Double.parseDouble($emlambda);
        tfParamw = new JTextField(IJ.d2s(sigma * emlambda / NA, decformat));

        String $sigma2 = tfSigma2.getText();
        sigma2 = Double.parseDouble($sigma2);
        String $emlambda2 = tfEmLambda2.getText();
        emlambda2 = Double.parseDouble($emlambda2);
        tfParamw2 = new JTextField(IJ.d2s(sigma2 * emlambda2 / NA, decformat));

        String $sigmaZ = tfSigmaZ.getText();
        tfParamz = new JTextField($sigmaZ);
        String $sigmaZ2 = tfSigmaZ2.getText();
        tfParamz2 = new JTextField($sigmaZ2);

        if (cbFitModel.getSelectedItem() == "FCS") {
            ccfrx = Integer.parseInt(tfCfXDistance.getText());
            ccfry = Integer.parseInt(tfCfYDistance.getText());
        } else {
            ccfrx = 0;
            ccfry = 0;
        }
        tfParamRx = new JTextField(IJ.d2s(pixelsize * 1000 / objmag * ccfrx, decformat), 8);
        tfParamRy = new JTextField(IJ.d2s(pixelsize * 1000 / objmag * ccfry, decformat), 8);
        tfParamRz = new JTextField(IJ.d2s(0, decformat), 8); // will be added later as a fit parameter

        tfParamQ2 = new JTextField("1", 8);
        tfParamN = new JTextField("1", 8);
        tfParamF2 = new JTextField("0", 8);
        tfParamD = new JTextField("1", 8);
        tfParamD2 = new JTextField("0", 8);
        tfParamF3 = new JTextField("0", 8);
        tfParamD3 = new JTextField("0", 8);
        tfParamQ3 = new JTextField("1", 8);
        tfParamVx = new JTextField("0", 8);
        tfParamVy = new JTextField("0", 8);
        tfParamG = new JTextField("0", 8);
        tfParamFtrip = new JTextField("0", 8);
        tfParamTtrip = new JTextField("0", 8);
        tfFitStart = new JTextField("1", 8);
        tfFitEnd = new JTextField(Integer.toString(chanum), 8);
        tfFitModel = new JTextField((String) cbFitModel.getSelectedItem());
        tfModProb1 = new JTextField("0", 8);
        tfModProb2 = new JTextField("0", 8);
        tfModProb3 = new JTextField("0", 8);

        // row 1
        fitframe.add(new JLabel("Fit Model: "));
        fitframe.add(tfFitModel);
        tfFitModel.setEditable(false);
        fitframe.add(new JLabel("GLS Fit"));
        fitframe.add(tbGLS);
        fitframe.add(new JLabel("Bayes"));
        fitframe.add(tbBayes);

        // row 2
        fitframe.add(new JLabel("a [nm]: "));
        fitframe.add(tfParama);
        tfParama.setEditable(false);
        fitframe.add(new JLabel("w [nm]: "));
        fitframe.add(tfParamw);
        tfParamw.setEditable(false);
        fitframe.add(new JLabel("z [nm]: "));
        fitframe.add(tfParamz);
        tfParamz.setEditable(false);

        // row 3
        fitframe.add(new JLabel(" "));
        fitframe.add(new JLabel(" "));
        fitframe.add(new JLabel("w2 [nm]: "));
        fitframe.add(tfParamw2);
        tfParamw2.setEditable(false);
        fitframe.add(new JLabel("z2 [nm]: "));
        fitframe.add(tfParamz2);
        tfParamz2.setEditable(false);

        // row 4
        fitframe.add(new JLabel("rx [nm]: "));
        fitframe.add(tfParamRx);
        tfParamRx.setEditable(false);
        fitframe.add(new JLabel("ry [nm]: "));
        fitframe.add(tfParamRy);
        tfParamRy.setEditable(false);
        fitframe.add(new JLabel("rz [nm]: "));
        fitframe.add(tfParamRz);
        tfParamRz.setEditable(false);

        // row 5 (empty)
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));

        // row 6
        fitframe.add(new JLabel("N: "));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("G: "));
        fitframe.add(new JLabel(""));
        fitframe.add(tbShow);
        fitframe.add(btnTest);

        // row 7
        fitframe.add(tfParamN);
        fitframe.add(rbtnHoldN);
        fitframe.add(tfParamG);
        fitframe.add(rbtnHoldG);
        fitframe.add(btnSetPar);
        fitframe.add(tbFixPar);

        // row 8
        fitframe.add(new JLabel("D [um2/s]: "));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("vx [um/s]: "));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("vy [um/s]: "));
        fitframe.add(new JLabel(""));

        // row 9
        fitframe.add(tfParamD);
        fitframe.add(rbtnHoldD);
        fitframe.add(tfParamVx);
        fitframe.add(rbtnHoldVx);
        rbtnHoldVx.setSelected(true);
        fitframe.add(tfParamVy);
        fitframe.add(rbtnHoldVy);
        rbtnHoldVy.setSelected(true);

        // row 10
        fitframe.add(new JLabel("D2 [um2/s]: "));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("Q2: "));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("F2: "));
        fitframe.add(new JLabel(""));

        // row 11
        fitframe.add(tfParamD2);
        fitframe.add(rbtnHoldD2);
        rbtnHoldD2.setSelected(true);
        fitframe.add(tfParamQ2);
        fitframe.add(rbtnHoldQ2);
        rbtnHoldQ2.setSelected(true);
        rbtnHoldQ2.setEnabled(false);
        fitframe.add(tfParamF2);
        fitframe.add(rbtnHoldF2);
        rbtnHoldF2.setSelected(true);

        //row 12
        fitframe.add(new JLabel("D3 [um2/s]: "));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("Q3: "));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("F3: "));
        fitframe.add(new JLabel(""));

        // row 13
        fitframe.add(tfParamD3);
        fitframe.add(rbtnHoldD3);
        rbtnHoldD3.setSelected(true);
        fitframe.add(tfParamQ3);
        fitframe.add(rbtnHoldQ3);
        rbtnHoldQ3.setSelected(true);
        rbtnHoldQ3.setEnabled(false);
        fitframe.add(tfParamF3);
        fitframe.add(rbtnHoldF3);
        rbtnHoldF3.setSelected(true);

        //row 14
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("Ftrip: "));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("Ttrip [us]: "));
        fitframe.add(new JLabel(""));

        // row 15
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(tfParamFtrip);
        fitframe.add(rbtnHoldFtrip);
        rbtnHoldFtrip.setSelected(true);
        fitframe.add(tfParamTtrip);
        fitframe.add(rbtnHoldTtrip);
        rbtnHoldTtrip.setSelected(true);

        // row 16 (empty)
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));

        // row 17
        fitframe.add(new JLabel("Fit start"));
        fitframe.add(tfFitStart);
        fitframe.add(new JLabel("Fit end"));
        fitframe.add(tfFitEnd);
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));

        // row 18 (empty)
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));

        // row 19
        fitframe.add(new JLabel("Bayesian"));
        fitframe.add(new JLabel("Model"));
        fitframe.add(new JLabel("Prob."));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel(""));

        // row 20
        fitframe.add(new JLabel("Model 1"));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("Model 2"));
        fitframe.add(new JLabel(""));
        fitframe.add(new JLabel("Model 3"));
        fitframe.add(new JLabel(""));

        // row 21 (empty)
        fitframe.add(tfModProb1);
        fitframe.add(new JLabel(""));
        fitframe.add(tfModProb2);
        fitframe.add(new JLabel(""));
        fitframe.add(tfModProb3);
        fitframe.add(new JLabel(""));

        // initial paramter settings
        initparam = new double[noparam];
        initparam[0] = 1.0;		// remember starting values
        initparam[1] = 1.0 / Math.pow(10, 12);
        initparam[2] = 0 / Math.pow(10, 6);
        initparam[3] = 0 / Math.pow(10, 6);
        initparam[4] = 0.0;
        initparam[5] = 0;
        initparam[6] = 0;
        initparam[7] = 0;
        initparam[8] = 0;
        initparam[9] = 0;
        initparam[10] = 0;

        cbFitModel.addItemListener(cbModelChanged);
        tbGLS.addItemListener(tbGLSPressed);
        tbBayes.addItemListener(tbBayesPressed);
        tbShow.addItemListener(tbShowPressed);
        btnTest.addActionListener(btnTestPressed);
        tbFixPar.addItemListener(tbFixParPressed);
        btnSetPar.addActionListener(btnSetParPressed);
    }

    ItemListener cbModelChanged = (ItemEvent ev) -> {
        if (cbFitModel.getSelectedItem() == "FCS") {
            tfFitModel.setText("FCS");
            //tfCfXDistance.setText(Integer.toString(0));
            //tfCfYDistance.setText(Integer.toString(0));
            ccfrx = Integer.parseInt(tfCfXDistance.getText());
            tfParamRx.setText(IJ.d2s(pixelsize * 1000 / objmag * ccfrx, decformat));
            ccfry = Integer.parseInt(tfCfYDistance.getText());
            tfParamRy.setText(IJ.d2s(pixelsize * 1000 / objmag * ccfry, decformat));
            tfParamRz.setText(IJ.d2s(0, decformat));
            tbFCCSDisplay.setSelected(false);
            tbFCCSDisplay.setText("Off");
            if (doFiltering) {			// if Thresholds settings panel exists
                doFiltering = false;
                tbFiltering.setSelected(false);
                if (filteringframe != null) {
                    filteringframe.setVisible(false);
                }
            }
        }
        if (cbFitModel.getSelectedItem() == "DC-FCCS") {
            //tfCfXDistance.setText(Integer.toString((int) width/2));
            //tfCfYDistance.setText(Integer.toString(0));
            tfParamRx.setText(IJ.d2s(0, decformat));
            tfParamRy.setText(IJ.d2s(0, decformat));
            tfParamRz.setText(IJ.d2s(0, decformat));
            tfParamVx.setText(IJ.d2s(0, decformat));
            tfParamVy.setText(IJ.d2s(0, decformat));
            rbtnHoldVx.setSelected(true);
            rbtnHoldVy.setSelected(true);
            tfFitModel.setText("DC-FCCS");
            if (doFiltering) {			// if Thresholds settings panel exists, close it
                doFiltering = false;
                tbFiltering.setSelected(false);
                if (filteringframe != null) {
                    filteringframe.setVisible(false);
                }
            }
        }
    };

    ItemListener tbGLSPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tbGLS.setText("On");
        } else {
            tbGLS.setSelected(false);
            tbGLS.setText("Off");
        }
    };

    ItemListener tbBayesPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tbBayes.setText("On");
        } else {
            tbBayes.setSelected(false);
            tbBayes.setText("Off");
        }
    };

    ItemListener tbShowPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tbShow.setSelected(true);
        } else {
            tbShow.setSelected(false);
        }
    };

    ActionListener btnTestPressed = (ActionEvent event) -> {
        prepareFit();
        plotTheoreticalCF();
    };

    ItemListener tbFixParPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            storepar = true;
            tbFixPar.setText("Fix");
        } else {
            storepar = false;
            tbFixPar.setText("Free");
        }
    };

    ActionListener btnSetParPressed = (ActionEvent event) -> {
        initparam[0] = Double.parseDouble(tfParamN.getText());
        initparam[1] = Double.parseDouble(tfParamD.getText()) / Math.pow(10, 12);
        initparam[2] = Double.parseDouble(tfParamVx.getText()) / Math.pow(10, 6);
        initparam[3] = Double.parseDouble(tfParamVy.getText()) / Math.pow(10, 6);
        initparam[4] = Double.parseDouble(tfParamG.getText());
        initparam[5] = Double.parseDouble(tfParamF2.getText());
        initparam[6] = Double.parseDouble(tfParamD2.getText()) / Math.pow(10, 12);
        initparam[7] = Double.parseDouble(tfParamF3.getText());
        initparam[8] = Double.parseDouble(tfParamD3.getText()) / Math.pow(10, 12);
        initparam[9] = Double.parseDouble(tfParamFtrip.getText());
        initparam[10] = Double.parseDouble(tfParamTtrip.getText()) / Math.pow(10, 6);
        tbFixPar.setSelected(true);
        tbFixPar.setText("Fix");
    };

    ActionListener cbBleachCorChanged = (ActionEvent event) -> {
        if (!expload) {
            if (cbBleachCor.getSelectedItem() == "Sliding Window") {
                SWSizeDialogue();
            }
            if (cbBleachCor.getSelectedItem() == "Polynomial") {
                PolynomialOrderDialogue();
            }
            if (cbBleachCor.getSelectedItem() == "Lin Segment") {
                SWSizeDialogue();
            }
        }
        bleachCorMem = (String) cbBleachCor.getSelectedItem();
    };

    ActionListener cbFilterChanged = (ActionEvent event) -> {
        if (!expload) {
            if (cbFilter.getSelectedItem() == "none") {
                filterLL = 0;
                filterUL = 65536;
            }
            if (cbFilter.getSelectedItem() == "Intensity") {
                filterDialogue();
            }
            if (cbFilter.getSelectedItem() == "Mean") {
                filterDialogue();
            }
        }
        filterMem = (String) cbFilter.getSelectedItem();
    };

    /* 
	 * FCS SIMULATION PANEL
	 * 
	 * public void createSimPanel(): create the JFrame for the simulation panel
	 * 
	 * followed by several listeners
	 * btnSimulatePressed()
	 * 
     */
    public void createSimPanel() {
        simframe.setFocusable(true);
        simframe.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        simframe.setLayout(new GridLayout(16, 4));
        simframe.setLocation(new Point(simPanelPosX, simPanelPosY));
        simframe.setSize(new Dimension(simPanelDimX, simPanelDimY));
        simframe.setResizable(false);

        cbSimMode = new JComboBox<>();
        cbSimMode.addItem("2D (free)");
        cbSimMode.addItem("2D (domains)");
        cbSimMode.addItem("2D (mesh)");
        cbSimMode.addItem("2D (dom+mesh");
        cbSimMode.addItem("3D (free)");

        tfSimSeed = new JTextField("1", 8);
        tfSimParticleNum = new JTextField("1000", 8);
        tfSimCPS = new JTextField("10000", 8);
        tfSimTauBleach = new JTextField("0", 8);
        tfSimPixelNum = new JTextField("21", 8);
        tfSimExtensionFactor = new JTextField("1.5", 8);
        tfSimTimeStepNum = new JTextField("50000", 8);
        tfSimFrameTime = new JTextField("0.001", 8);
        tfSimStepsPerFrame = new JTextField("10", 8);
        tfSimCurrentStepSize = new JTextField("20", 8);
        tfSimD1 = new JTextField("1.0", 8);
        tfSimDoutDinRatio = new JTextField("1.0", 8);
        tfSimD2 = new JTextField("0.1", 8);
        tfSimF2 = new JTextField("0.0", 8);
        tfSimD3 = new JTextField("0.01", 8);
        tfSimF3 = new JTextField("0.0", 8);
        tfSimKon = new JTextField("300", 8);
        tfSimKoff = new JTextField("700", 8);
        tfSimCameraOffset = new JTextField("100", 8);
        tfSimCameraNoiseFactor = new JTextField("3", 8);
        tfSimBleachRadius = new JTextField("3.0", 8);
        tfSimBleachFrame = new JTextField("1000000", 8);
        tfDomainRadius = new JTextField("30", 8);
        tfDomainDensity = new JTextField("30", 8);
        tfPin = new JTextField("1.0", 8);
        tfPout = new JTextField("0.6", 8);
        tfMeshworkSize = new JTextField("100", 8);
        tfHopProbability = new JTextField("1", 8);
        btnSimulate = new JButton("Simulate");
        btnBatchSim = new JButton("Batch");
        btnStopSimulation = new JButton("Stop");
        tbSimTrip = new JToggleButton("Triplet off");

        // disbale Stop button; is enabled only when simulations run
        tfSimCurrentStepSize.setEnabled(false);
        tfSimDoutDinRatio.setEnabled(false);
        tfSimKon.setEnabled(false);
        tfSimKoff.setEnabled(false);
        tfDomainRadius.setEnabled(false);
        tfDomainDensity.setEnabled(false);
        tfPin.setEnabled(false);
        tfPout.setEnabled(false);
        tfMeshworkSize.setEnabled(false);
        tfHopProbability.setEnabled(false);
        btnStopSimulation.setEnabled(false);

        // help texts
        JLabel jlSimSeed = new JLabel("Seed");
        jlSimSeed.setToolTipText("Integer: Seed for the random number generator. Using the same seed (>0) leads to reproducible simulations.");
        JLabel jlSimParticleNum = new JLabel("Particle #");
        jlSimParticleNum.setToolTipText("Integer: Seed for the random number generator. Using the same seed (>0) leads to reproducible simulations.");
        JLabel jlSimCPS = new JLabel("CPS");
        jlSimCPS.setToolTipText("Integer: counts per particle per second; brightness of the moleucles");
        JLabel jlSimTauBleach = new JLabel("Bleach time");
        jlSimTauBleach.setToolTipText("Integer: Characteristic bleach time in seconds (based on an exponential). Set to 0 for no bleaching");
        JLabel jlSimPixelNum = new JLabel("Pixel #");
        jlSimPixelNum.setToolTipText("Integer: number of pixels in x and y direction to be simulated. Only square areas are used.");
        JLabel jlSimExtensionFactor = new JLabel("Extension");
        jlSimExtensionFactor.setToolTipText("Double: ratio of simulated to observed region.");
        JLabel jlSimTimeStepNum = new JLabel("Frame #");
        jlSimTimeStepNum.setToolTipText("Integer: Numbers of frames to be simulated.");
        JLabel jlSimFrameTime = new JLabel("Time res");
        jlSimFrameTime.setToolTipText("Double: time per frame in seconds.  Press enter to calcualte new step size.");
        JLabel jlSimStepsPerFrame = new JLabel("Steps per Frame");
        jlSimStepsPerFrame.setToolTipText("Double: number of simulation steps per frame. Press enter to calcualte new step size.");
        JLabel jlSimCurrentStepSize = new JLabel("Step Size [nm]");
        jlSimCurrentStepSize.setToolTipText("Double: Shows the current step size in the simulations based on D1 and time per simualtion step.");
        JLabel jlSimD1 = new JLabel("D1 [um2/s]");
        jlSimD1.setToolTipText("Double: Diffusion coefficient of first species to be simulated.  Press enter to calcualte new step size.");
        JLabel jlSimDoutDinRatio = new JLabel("Dout/Din");
        jlSimDoutDinRatio.setToolTipText("Double: Ratio of diffuion coeeficients of particles outside and inside domains.");
        JLabel jlSimD2 = new JLabel("D2 [um2/s]");
        jlSimD2.setToolTipText("Double: Diffusion coefficient of second species to be simulated (if any).");
        JLabel jlSimF2 = new JLabel("F2");
        jlSimF2.setToolTipText("Double (0 < F2 < 1): Fraction of particles of the total for the second species.");
        JLabel jlSimD3 = new JLabel("D3 [um2/s]");
        jlSimD3.setToolTipText("Double: Diffusion coefficient of third species to be simulated (if any).");
        JLabel jlSimF3 = new JLabel("F3");
        JLabel jlSimKon = new JLabel("kon (triplet)");
        jlSimKon.setToolTipText("Double: on rate for transition to triplet.");
        JLabel jlSimKoff = new JLabel("koff (triplet)");
        jlSimKoff.setToolTipText("Double: off rate for transition from triplet.");
        jlSimF3.setToolTipText("Double (0 < F3 < 1 AND F2 + F3 < 1): Fraction of particles of the total for the second species.");
        JLabel jlSimCameraOffset = new JLabel("Cam Offset");
        jlSimCameraOffset.setToolTipText("Integer: Offset of the CCD camera.");
        JLabel jlSimCameraNoiseFactor = new JLabel("Cam Noise");
        jlSimCameraNoiseFactor.setToolTipText("Integer: noise factor of the camera.");
        JLabel jlSimBleachRadius = new JLabel("FRAP Radius [um]");
        jlSimBleachRadius.setToolTipText("Double: Radius in um within which the particles will be bleached. Only available in 2D.");
        JLabel jlSimBleachFrame = new JLabel("FRAP Frame");
        jlSimBleachFrame.setToolTipText("Integer: Frame at which the assumed bleach pulse happens. Bleachign is assumed to be instantaneous.");
        JLabel jlDomainRadius = new JLabel("Dom Rad [nm]");
        jlDomainRadius.setToolTipText("Radius of domains in nm");
        JLabel jlDomainDensity = new JLabel("Dom Density");
        jlDomainDensity.setToolTipText("Domain Density in numbers/um2");
        JLabel jlPin = new JLabel("Pin");
        jlPin.setToolTipText("Probability to enter a domain");
        JLabel jlPout = new JLabel("Pout");
        jlPout.setToolTipText("Probability to exit a domain");
        JLabel jlMeshworkSize = new JLabel("Mesh Size [nm]");
        jlMeshworkSize.setToolTipText("Mesh size in nm");
        JLabel jlHopProbability = new JLabel("Hop Prob");
        jlHopProbability.setToolTipText("Probablility to hop over a barrier in the meshwork.");
        btnBatchSim.setToolTipText("Run multiple simulations.");
        btnStopSimulation.setToolTipText("Stops running simulation.");

        //row 1
        simframe.add(new JLabel("Mode"));
        simframe.add(cbSimMode);
        simframe.add(new JLabel(""));
        simframe.add(tbSimTrip);

        //row 2
        simframe.add(jlSimSeed);
        simframe.add(tfSimSeed);
        simframe.add(jlSimParticleNum);
        simframe.add(tfSimParticleNum);

        //row 3
        simframe.add(jlSimCPS);
        simframe.add(tfSimCPS);
        simframe.add(jlSimTauBleach);
        simframe.add(tfSimTauBleach);

        //row 4
        simframe.add(jlSimPixelNum);
        simframe.add(tfSimPixelNum);
        simframe.add(jlSimExtensionFactor);
        simframe.add(tfSimExtensionFactor);

        //row 5
        simframe.add(jlSimTimeStepNum);
        simframe.add(tfSimTimeStepNum);
        simframe.add(jlSimFrameTime);
        simframe.add(tfSimFrameTime);

        //row 6
        simframe.add(jlSimStepsPerFrame);
        simframe.add(tfSimStepsPerFrame);
        simframe.add(jlSimCurrentStepSize);
        simframe.add(tfSimCurrentStepSize);

        //row 7
        simframe.add(jlSimD1);
        simframe.add(tfSimD1);
        simframe.add(jlSimDoutDinRatio);
        simframe.add(tfSimDoutDinRatio);

        //row 8
        simframe.add(jlSimD2);
        simframe.add(tfSimD2);
        simframe.add(jlSimF2);
        simframe.add(tfSimF2);

        //row 9
        simframe.add(jlSimD3);
        simframe.add(tfSimD3);
        simframe.add(jlSimF3);
        simframe.add(tfSimF3);

        //row 10
        simframe.add(jlSimKon);
        simframe.add(tfSimKon);
        simframe.add(jlSimKoff);
        simframe.add(tfSimKoff);

        //row 11
        simframe.add(jlSimCameraOffset);
        simframe.add(tfSimCameraOffset);
        simframe.add(jlSimCameraNoiseFactor);
        simframe.add(tfSimCameraNoiseFactor);

        //row 12
        simframe.add(jlSimBleachRadius);
        simframe.add(tfSimBleachRadius);
        simframe.add(jlSimBleachFrame);
        simframe.add(tfSimBleachFrame);

        //rox 13
        simframe.add(jlDomainRadius);
        simframe.add(tfDomainRadius);
        simframe.add(jlDomainDensity);
        simframe.add(tfDomainDensity);

        //rox 14
        simframe.add(jlPin);
        simframe.add(tfPin);
        simframe.add(jlPout);
        simframe.add(tfPout);

        //rox 15
        simframe.add(jlMeshworkSize);
        simframe.add(tfMeshworkSize);
        simframe.add(jlHopProbability);
        simframe.add(tfHopProbability);

        //rox 16
        simframe.add(new JLabel(""));
        simframe.add(btnBatchSim);
        simframe.add(btnStopSimulation);
        simframe.add(btnSimulate);

        cbSimMode.addActionListener(cbSimModeChanged);
        btnSimulate.addActionListener(btnSimulatePressed);
        btnBatchSim.addActionListener(btnBatchSimPressed);
        btnStopSimulation.addActionListener(btnStopSimulationPressed);
        tfSimFrameTime.addActionListener(tfStepSizeChanged);
        tfSimStepsPerFrame.addActionListener(tfStepSizeChanged);
        tfSimD1.addActionListener(tfStepSizeChanged);
        tbSimTrip.addItemListener(tbSimTripPressed);
    }

    ActionListener cbSimModeChanged = (ActionEvent event) -> {
        if (cbSimMode.getSelectedItem().toString().contains("2D")) {
            tfSimBleachRadius.setEnabled(true);
            tfSimBleachFrame.setEnabled(true);
        }
        if (cbSimMode.getSelectedItem().toString().contains("3D")) {
            tfSimBleachRadius.setEnabled(false);
            tfSimBleachFrame.setEnabled(false);
        }
        if (cbSimMode.getSelectedItem().toString().contains("free")) {
            tfSimDoutDinRatio.setEnabled(false);
            tfDomainRadius.setEnabled(false);
            tfDomainDensity.setEnabled(false);
            tfPin.setEnabled(false);
            tfPout.setEnabled(false);
            tfMeshworkSize.setEnabled(false);
            tfHopProbability.setEnabled(false);
            simDomainFlag = false;
            simMeshFlag = false;
        }
        if (cbSimMode.getSelectedItem().toString().contains("dom")) {
            tfSimDoutDinRatio.setEnabled(true);
            tfDomainRadius.setEnabled(true);
            tfDomainDensity.setEnabled(true);
            tfPin.setEnabled(true);
            tfPout.setEnabled(true);
            simDomainFlag = true;
            if (cbSimMode.getSelectedItem().toString().contains("mesh")) {
                tfMeshworkSize.setEnabled(true);
                tfHopProbability.setEnabled(true);
                simMeshFlag = true;
            } else {
                tfMeshworkSize.setEnabled(false);
                tfHopProbability.setEnabled(false);
                simMeshFlag = false;
            }
        } else {
            simDomainFlag = false;
            simMeshFlag = false;
            if (cbSimMode.getSelectedItem().toString().contains("mesh")) {
                tfSimDoutDinRatio.setEnabled(false);
                tfDomainRadius.setEnabled(false);
                tfDomainDensity.setEnabled(false);
                tfPin.setEnabled(false);
                tfPout.setEnabled(false);
                tfMeshworkSize.setEnabled(true);
                tfHopProbability.setEnabled(true);
                simMeshFlag = true;
            }
        }
    };

    ActionListener btnSimulatePressed = (ActionEvent event) -> {
        boolean ask3D;
        ask3D = !cbSimMode.getSelectedItem().toString().contains("2D");
        /*
        GenericDialog gd = new GenericDialog("particle number");
        gd.addNumericField("#: ", -1, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        tmpPartNum = (int) gd.getNextNumber();
         */
        simulateACF(ask3D);
    };

    ActionListener btnBatchSimPressed = (ActionEvent event) -> {
        boolean ask3D;
        ask3D = !cbSimMode.getSelectedItem().toString().contains("2D");
        GenericDialog gd = new GenericDialog("particle number");
        gd.addNumericField("D start ", 1, 1);
        gd.addNumericField("D end ", 10, 1);
        gd.addNumericField("D step ", 1, 1);
        gd.addNumericField("D2 start ", 1, 1);
        gd.addNumericField("D2 end ", 10, 1);
        gd.addNumericField("D2 step ", 1, 1);
        gd.addNumericField("F2 start ", 0, 0);
        gd.addNumericField("F2 end ", 1, 0);
        gd.addNumericField("F2 step ", 0.1, 2);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }

//        tmpPartNum = -1;// for batches we don't put out traces
        batchDStart = gd.getNextNumber();
        batchDEnd = gd.getNextNumber();
        batchDStep = gd.getNextNumber();
        batchD2Start = gd.getNextNumber();
        batchD2End = gd.getNextNumber();
        batchD2Step = gd.getNextNumber();
        batchF2Start = gd.getNextNumber();
        batchF2End = gd.getNextNumber();
        batchF2Step = gd.getNextNumber();

        JFileChooser fc = new JFileChooser($imagePath);
        fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        fc.setMultiSelectionEnabled(false);
        int returnVal = fc.showDialog(null, "Choose directory");
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            simBatchPath = fc.getSelectedFile();
            batchSimulateACFWorker batchSimulateACFInstant = new batchSimulateACFWorker(ask3D);
            batchSimulateACFInstant.execute();
        } else {
            return;
        }
    };

    ActionListener btnStopSimulationPressed = (ActionEvent event) -> {
        if (batchSim) {
            batchSimulateACFInstant.cancel(true);
            simulateACFInstant.cancel(true);
            batchSim = false;
        } else {
            simulateACFInstant.cancel(true);
        }
    };

    ActionListener tfStepSizeChanged = (ActionEvent event) -> {
        double diff = Double.parseDouble(tfSimD1.getText()) / Math.pow(10, 12);
        double timestep = Double.parseDouble(tfSimFrameTime.getText());
        int simsteps = Integer.parseInt(tfSimStepsPerFrame.getText());
        tfSimCurrentStepSize.setText(IJ.d2s(Math.sqrt(4 * diff * timestep / simsteps) * Math.pow(10, 9), decformat));
    };

    ItemListener tbSimTripPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tbSimTrip.setText("Triplet On");
            tfSimKon.setEnabled(true);
            tfSimKoff.setEnabled(true);
            simBlinkFlag = true;
        } else {
            tbSimTrip.setText("Triplet Off");
            tfSimKon.setEnabled(false);
            tfSimKoff.setEnabled(false);
            simBlinkFlag = false;
        }
    };

    /* 
	 * Diffusion Law Panel
	 * public void createDiffLawPanel(): create the JFrame for the Diffusion Law panel
	 * 
	 * followed by two listeners
	 * btnDiffLawCalculatePressed()
	 * btnDiffLawFitPressed()
	 * 
     */
    public void createDiffLawPanel() {
        difflawframe.setFocusable(true);
        difflawframe.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        difflawframe.setLayout(new GridLayout(8, 5));
        difflawframe.setLocation(new Point(difflawPanelPosX, difflawPanelPosY));
        difflawframe.setSize(new Dimension(difflawPanelDimX, difflawPanelDimY));
        difflawframe.setResizable(false);

        tfDLBinStart = new JTextField("1", 8);
        tfDLBinStart.setEditable(true);
        tfDLBinEnd = new JTextField("5", 8);
        tfDLFitStart = new JTextField("1", 8);
        tfDLFitEnd = new JTextField("5", 8);
        tfDLROI = new JTextField("7", 8);
        tfDLROI.setEditable(false);
        btnDiffLawCalculate = new JButton("Cal");
        btnDiffLawFit = new JButton("Fit");
        tbDLROI = new JToggleButton("All");

        // help texts
        btnDiffLawCalculate.setToolTipText("Calculate the diffusion law between the given 'Start-End' range.");
        btnDiffLawFit.setToolTipText("Fit the diffusion law between the given 'Start-End' range..");
        tbDLROI.setToolTipText("Set whether DL is to be claculated over the full image or over smaller ROIs.");

        //row 1
        difflawframe.add(new JLabel("Calc"));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));

        //row 2
        difflawframe.add(new JLabel("Binning:"));
        difflawframe.add(new JLabel("Start"));
        difflawframe.add(tfDLBinStart);
        difflawframe.add(new JLabel("End"));
        difflawframe.add(tfDLBinEnd);

        //row 3
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(btnDiffLawCalculate);

        //row 4
        difflawframe.add(new JLabel("Fitting"));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));

        //row 5
        difflawframe.add(new JLabel("Binning:"));
        difflawframe.add(new JLabel("Start"));
        difflawframe.add(tfDLFitStart);
        difflawframe.add(new JLabel("End"));
        difflawframe.add(tfDLFitEnd);

        //row 6
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(btnDiffLawFit);

        //row 7
        difflawframe.add(new JLabel("ROI"));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));

        //row 8
        difflawframe.add(tbDLROI);
        difflawframe.add(new JLabel(""));
        difflawframe.add(tfDLROI);
        difflawframe.add(new JLabel(""));
        difflawframe.add(new JLabel(""));

        btnDiffLawCalculate.addActionListener(btnDiffLawCalculatePressed);
        btnDiffLawFit.addActionListener(btnDiffLawFitPressed);
        tbDLROI.addItemListener(tbDLROIPressed);
    }

    ActionListener btnDiffLawCalculatePressed = (ActionEvent event) -> {
        int bs;
        int be;
        try {
            bs = Integer.parseInt(tfDLBinStart.getText());
            be = Integer.parseInt(tfDLBinEnd.getText());
        } catch (NumberFormatException nfe) {
            IJ.showMessage("One of the values in the DL window does not have the right format (integer).");
            throw new NumberFormatException("Number format error.");
        }
        DLbsmem = bs;
        DLbemem = be;
        if (setParameters()) {
            if (bs > 0 && be > bs && be < DiffLawMaxPoint) {
                correlateDiffLawWorker correlateDiffLawInstant = new correlateDiffLawWorker();
                correlateDiffLawInstant.execute();
            } else {
                IJ.showMessage("Points are out of range.");
            }
        }
    };

    ActionListener btnDiffLawFitPressed = (ActionEvent event) -> {
        int fs;
        int fe;
        // get and check that values are integer
        try {
            fs = Integer.parseInt(tfDLFitStart.getText());
            fe = Integer.parseInt(tfDLFitEnd.getText());
        } catch (NumberFormatException nfe) {
            IJ.showMessage("One of the values in the DL window does not have the right format (integer).");
            throw new NumberFormatException("Number format error.");
        }

        // if values within range perform DL calculations
        if (fs >= DLbsmem && fe <= DLbemem && fs > 0 && fe > fs) {
            diffLawFitLim[0] = fs;
            diffLawFitLim[1] = fe;
            fitDiffLaw();
        } else {
            IJ.showMessage("Fit points are out of range.");
        }
    };

    ItemListener tbDLROIPressed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tbDLROI.setText("ROI");
            // tfDLROI.setEditable(true); not yet allowed; perhaps determining subregion sizes that is includued later
            tfDLBinEnd.setText("5");
            tfDLFitStart.setText("1");
            tfDLFitEnd.setText("5");
            tfDLBinEnd.setEditable(false);
            tfDLFitStart.setEditable(false);
            tfDLFitEnd.setEditable(false);
        } else {
            tbDLROI.setText("All");
            tfDLROI.setEditable(false);
            tfDLBinEnd.setEditable(true);
            tfDLFitStart.setEditable(true);
            tfDLFitEnd.setEditable(true);
        }
    };

    /*
	 * N&B Panel
	 * public void createNBPanel(): create the JFrame for the NB panel
	 * 
	 * followed by
	 * ActionListener cbNMModechanged()
	 * 
     */
    public void createNBPanel() {
        NBframe.setFocusable(true);
        NBframe.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        NBframe.setLayout(new GridLayout(4, 2));
        NBframe.setLocation(new Point(NBPosX, NBPosY));
        NBframe.setSize(new Dimension(NBDimX, NBDimY));
        NBframe.setResizable(false);

        cbNBMode = new JComboBox<>();
        cbNBMode.addItem("G1");
        cbNBMode.addItem("Calibrated");

        tfNBS = new JTextField(Double.toString(NBslope), 3);
        tfNBCalibRatio = new JTextField("2", 3);

        tfNBCalibRatio.setEnabled(false);
        tfNBS.setEnabled(false);

        btnNB = new JButton("N&B");

        // help texts
        btnNB.setToolTipText("perform N&B analysis");

        NBframe.add(new JLabel("NB mode"));
        NBframe.add(cbNBMode);
        NBframe.add(new JLabel("Calib Ratio"));
        NBframe.add(tfNBCalibRatio);
        NBframe.add(new JLabel("N&B analysis"));
        NBframe.add(btnNB);
        NBframe.add(new JLabel("S value"));
        NBframe.add(tfNBS);

        cbNBMode.addActionListener(cbNBModeChanged);
        btnNB.addActionListener(btnNBPressed);
    }

    ActionListener cbNBModeChanged = (ActionEvent event) -> {
        cbNBMode.setEnabled(false);
        tfNBCalibRatio.setEditable(false);
        btnNB.setEnabled(false);
        tfNBS.setEditable(false);

        if (cbNBMode.getSelectedItem().toString().contains("G1")) {
            tfNBCalibRatio.setEnabled(false);
            tfNBS.setEnabled(false);
        } else {
            tfNBCalibRatio.setEnabled(true);
            tfNBS.setEnabled(true);
            //load bgr file
            if (!bgrloaded) {
                tbBGR.setText("Bgr FILE");
                if (loadBGRFile()) {
                    bgrloaded = true;
                } else {
                    tbBGR.setText("Bgr NUM");
                    tbBGR.setSelected(false);
                    bgrloaded = false;
                    tfNBCalibRatio.setEnabled(false);
                    tfNBS.setEnabled(false);
                    cbNBMode.setSelectedIndex(0);
                }
            }
            //calibrate
            if (setParameters()) {
                camCalibrate();
            }
        }

        cbNBMode.setEnabled(true);
        tfNBCalibRatio.setEditable(true);
        btnNB.setEnabled(true);
        tfNBS.setEditable(true);
    };

    ActionListener btnNBPressed = (ActionEvent event) -> {
        if (setImp == false) {
            IJ.showMessage("No image loaded or assigned.");
            tbNB.setSelected(false); //this makes NB frame invisible and resets the toggle button to off
            return;
        }
        if (setParameters()) {
            performNB(imp, "evaluation");
        }
    };

    /* 
	 * FILTERING PANEL
	 * 
	 * public void createFiltering(): create the JFrame for filtering panel and make it visible when "Thresholds" are switched on in the main panel
	 * 
	 * followed by several functions:
	 * ReadFilteringFrame()
	 * SetThresholds()
	 * listeners:
	 * btnFilterPressed()
	 * btnResetPressed()
	 * and a number of radiobutton listeners
	 * 
     */
    // create the FCS filtering panel
    public void createFiltering() {

        filteringframe = new JFrame("Thresholds settings");

        filteringframe.setFocusable(true);
        filteringframe.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        filteringframe.setLayout(new GridLayout(17, 6));
        filteringframe.setLocation(new Point(filteringPanelPosX, filteringPanelPosY));
        filteringframe.setSize(new Dimension(filteringPanelDimX, filteringPanelDimY));
        filteringframe.setResizable(false);

        // initialize arrays for filtering:
        paramfilter = new boolean[noparam + 2];		// include chi2 value and relative G
        filterThresholds = new double[3][noparam + 2][2];		// include chi2 value and relative G
        userThreshold = new boolean[5];

        userThreshold[0] = false;			//user has not set any threshold

        tfAlFiltN = new JTextField();
        tfAhFiltN = new JTextField();
        tfAlFiltD = new JTextField();
        tfAhFiltD = new JTextField();
        tfAlFiltvx = new JTextField();
        tfAhFiltvx = new JTextField();
        tfAlFiltvy = new JTextField();
        tfAhFiltvy = new JTextField();
        tfAlFiltG = new JTextField();
        tfAhFiltG = new JTextField();
        tfAhFiltGrel = new JTextField();
        tfAlFiltF2 = new JTextField();
        tfAhFiltF2 = new JTextField();
        tfAlFiltD2 = new JTextField();
        tfAhFiltD2 = new JTextField();
        tfAlFiltF3 = new JTextField();
        tfAhFiltF3 = new JTextField();
        tfAlFiltD3 = new JTextField();
        tfAhFiltD3 = new JTextField();
        tfAlFiltFtrip = new JTextField();
        tfAhFiltFtrip = new JTextField();
        tfAlFiltTtrip = new JTextField();
        tfAhFiltTtrip = new JTextField();
        tfAlFiltChi2 = new JTextField();
        tfAhFiltChi2 = new JTextField();
        tfClFiltN = new JTextField();
        tfChFiltN = new JTextField();
        tfClFiltD = new JTextField();
        tfChFiltD = new JTextField();
        tfClFiltvx = new JTextField();
        tfChFiltvx = new JTextField();
        tfClFiltvy = new JTextField();
        tfChFiltvy = new JTextField();
        tfClFiltG = new JTextField();
        tfChFiltG = new JTextField();
        tfChFiltGrel = new JTextField();
        tfClFiltF2 = new JTextField();
        tfChFiltF2 = new JTextField();
        tfClFiltD2 = new JTextField();
        tfChFiltD2 = new JTextField();
        tfClFiltF3 = new JTextField();
        tfChFiltF3 = new JTextField();
        tfClFiltD3 = new JTextField();
        tfChFiltD3 = new JTextField();
        tfClFiltFtrip = new JTextField();
        tfChFiltFtrip = new JTextField();
        tfClFiltTtrip = new JTextField();
        tfChFiltTtrip = new JTextField();
        tfClFiltChi2 = new JTextField();
        tfChFiltChi2 = new JTextField();

        // row 1
        filteringframe.add(new JLabel(" "));
        filteringframe.add(new JLabel(" "));
        filteringframe.add(new JLabel("ACF"));
        filteringframe.add(new JLabel(" "));
        filteringframe.add(new JLabel("CCF"));
        filteringframe.add(new JLabel(" "));

        // row 2
        filteringframe.add(new JLabel(" "));
        filteringframe.add(new JLabel("filter"));
        filteringframe.add(new JLabel("Min."));
        filteringframe.add(new JLabel("Max."));
        filteringframe.add(new JLabel("Min."));
        filteringframe.add(new JLabel("Max."));

        // row 3
        filteringframe.add(new JLabel("N"));
        filteringframe.add(rbtnFiltN);
        filteringframe.add(tfAlFiltN);
        filteringframe.add(tfAhFiltN);
        filteringframe.add(tfClFiltN);
        filteringframe.add(tfChFiltN);

        // row 4
        filteringframe.add(new JLabel("D [um2/s]"));
        filteringframe.add(rbtnFiltD);
        filteringframe.add(tfAlFiltD);
        filteringframe.add(tfAhFiltD);
        filteringframe.add(tfClFiltD);
        filteringframe.add(tfChFiltD);

        // row 5
        filteringframe.add(new JLabel("vx [um/s]"));
        filteringframe.add(rbtnFiltvx);
        filteringframe.add(tfAlFiltvx);
        filteringframe.add(tfAhFiltvx);
        filteringframe.add(tfClFiltvx);
        filteringframe.add(tfChFiltvx);

        // row 6
        filteringframe.add(new JLabel("vy [um/s]"));
        filteringframe.add(rbtnFiltvy);
        filteringframe.add(tfAlFiltvy);
        filteringframe.add(tfAhFiltvy);
        filteringframe.add(tfClFiltvy);
        filteringframe.add(tfChFiltvy);

        // row 7
        filteringframe.add(new JLabel("G"));
        filteringframe.add(rbtnFiltG);
        filteringframe.add(tfAlFiltG);
        filteringframe.add(tfAhFiltG);
        filteringframe.add(tfClFiltG);
        filteringframe.add(tfChFiltG);

        // row 8
        filteringframe.add(new JLabel("rel. G"));
        filteringframe.add(rbtnFiltGrel);
        filteringframe.add(new JLabel("+/- G*N"));
        filteringframe.add(tfAhFiltGrel);
        filteringframe.add(new JLabel(" "));
        filteringframe.add(tfChFiltGrel);

        // row 9
        filteringframe.add(new JLabel("F2"));
        filteringframe.add(rbtnFiltF2);
        filteringframe.add(tfAlFiltF2);
        filteringframe.add(tfAhFiltF2);
        filteringframe.add(tfClFiltF2);
        filteringframe.add(tfChFiltF2);

        // row 10
        filteringframe.add(new JLabel("D2 [um2/s]"));
        filteringframe.add(rbtnFiltD2);
        filteringframe.add(tfAlFiltD2);
        filteringframe.add(tfAhFiltD2);
        filteringframe.add(tfClFiltD2);
        filteringframe.add(tfChFiltD2);

        // row 11
        filteringframe.add(new JLabel("F3"));
        filteringframe.add(rbtnFiltF3);
        filteringframe.add(tfAlFiltF3);
        filteringframe.add(tfAhFiltF3);
        filteringframe.add(tfClFiltF3);
        filteringframe.add(tfChFiltF3);

        // row 12
        filteringframe.add(new JLabel("D3 [um2/s]"));
        filteringframe.add(rbtnFiltD3);
        filteringframe.add(tfAlFiltD3);
        filteringframe.add(tfAhFiltD3);
        filteringframe.add(tfClFiltD3);
        filteringframe.add(tfChFiltD3);

        // row 13
        filteringframe.add(new JLabel("Ftrip"));
        filteringframe.add(rbtnFiltFtrip);
        filteringframe.add(tfAlFiltFtrip);
        filteringframe.add(tfAhFiltFtrip);
        filteringframe.add(tfClFiltFtrip);
        filteringframe.add(tfChFiltFtrip);

        // row 14
        filteringframe.add(new JLabel("Ttrip"));
        filteringframe.add(rbtnFiltTtrip);
        filteringframe.add(tfAlFiltTtrip);
        filteringframe.add(tfAhFiltTtrip);
        filteringframe.add(tfClFiltTtrip);
        filteringframe.add(tfChFiltTtrip);

        // row 15
        filteringframe.add(new JLabel("Chi2"));
        filteringframe.add(rbtnFiltChi2);
        filteringframe.add(tfAlFiltChi2);
        filteringframe.add(tfAhFiltChi2);
        filteringframe.add(tfClFiltChi2);
        filteringframe.add(tfChFiltChi2);

        // row 16
        filteringframe.add(btnFilter);
        filteringframe.add(new JLabel(" "));
        filteringframe.add(new JLabel("CCF "));
        filteringframe.add(new JLabel("same as"));
        filteringframe.add(new JLabel(" ACF"));
        filteringframe.add(rbtnCuseA);

        // row 17
        filteringframe.add(btnReset);
        filteringframe.add(new JLabel(" "));
        filteringframe.add(new JLabel("If CCF"));
        filteringframe.add(new JLabel("invalid"));
        filteringframe.add(new JLabel(" q=0"));
        filteringframe.add(rbtnReplaceZero);

        btnFilter.setToolTipText("Creates filtering mask according to specified thresholds and applies it on the parameter maps");
        btnReset.setToolTipText("Resets the thresholds to their default values and the filtering mask to 1.0 for all fitted pixels");
        rbtnCuseA.setToolTipText("Sets the thresholds for CCF to be the same as for ACFs");
        rbtnReplaceZero.setToolTipText("When calculating cross-correlation amount (q) sets q = 0 for pixels where the ACFs pass the thresholds but the CCF does not");

        btnFilter.addActionListener(btnFilterPressed);
        btnReset.addActionListener(btnResetPressed);
        rbtnCuseA.addItemListener(rbtnCuseAChanged);
        rbtnFiltG.addItemListener(rbtnFiltGChanged);
        rbtnFiltGrel.addItemListener(rbtnFiltGrelChanged);
        rbtnFiltN.addItemListener(rbtnFiltNChanged);
        rbtnFiltD.addItemListener(rbtnFiltDChanged);
        rbtnFiltvx.addItemListener(rbtnFiltvxChanged);
        rbtnFiltvy.addItemListener(rbtnFiltvyChanged);
        rbtnFiltF2.addItemListener(rbtnFiltF2Changed);
        rbtnFiltD2.addItemListener(rbtnFiltD2Changed);
        rbtnFiltF3.addItemListener(rbtnFiltF3Changed);
        rbtnFiltD3.addItemListener(rbtnFiltD3Changed);
        rbtnFiltFtrip.addItemListener(rbtnFiltFtripChanged);
        rbtnFiltTtrip.addItemListener(rbtnFiltTtripChanged);
        rbtnFiltChi2.addItemListener(rbtnFiltChi2Changed);
    }

    public void ReadFilteringFrame() {					// reads the settings from the "thresholds settings" panel and saves them in filterThresholds[][][], paramfilter[] and userThreshold[]

        paramfilter[0] = rbtnFiltN.isSelected();
        paramfilter[1] = rbtnFiltD.isSelected();
        paramfilter[2] = rbtnFiltvx.isSelected();
        paramfilter[3] = rbtnFiltvy.isSelected();
        paramfilter[4] = rbtnFiltG.isSelected();
        paramfilter[5] = rbtnFiltF2.isSelected();
        paramfilter[6] = rbtnFiltD2.isSelected();
        paramfilter[7] = rbtnFiltF3.isSelected();
        paramfilter[8] = rbtnFiltD3.isSelected();
        paramfilter[9] = rbtnFiltFtrip.isSelected();
        paramfilter[10] = rbtnFiltTtrip.isSelected();
        paramfilter[noparam] = rbtnFiltChi2.isSelected();
        paramfilter[noparam + 1] = rbtnFiltGrel.isSelected();

        try {
            filterThresholds[0][0][0] = Double.parseDouble(tfAlFiltN.getText());
            filterThresholds[0][0][1] = Double.parseDouble(tfAhFiltN.getText());
            filterThresholds[0][1][0] = Double.parseDouble(tfAlFiltD.getText()) * Math.pow(10, -12);
            filterThresholds[0][1][1] = Double.parseDouble(tfAhFiltD.getText()) * Math.pow(10, -12);
            filterThresholds[0][2][0] = Double.parseDouble(tfAlFiltvx.getText()) * Math.pow(10, -6);
            filterThresholds[0][2][1] = Double.parseDouble(tfAhFiltvx.getText()) * Math.pow(10, -6);
            filterThresholds[0][3][0] = Double.parseDouble(tfAlFiltvy.getText()) * Math.pow(10, -6);
            filterThresholds[0][3][1] = Double.parseDouble(tfAhFiltvy.getText()) * Math.pow(10, -6);
            filterThresholds[0][4][0] = Double.parseDouble(tfAlFiltG.getText());
            filterThresholds[0][4][1] = Double.parseDouble(tfAhFiltG.getText());
            filterThresholds[0][5][0] = Double.parseDouble(tfAlFiltF2.getText());
            filterThresholds[0][5][1] = Double.parseDouble(tfAhFiltF2.getText());
            filterThresholds[0][6][0] = Double.parseDouble(tfAlFiltD2.getText()) * Math.pow(10, -12);
            filterThresholds[0][6][1] = Double.parseDouble(tfAhFiltD2.getText()) * Math.pow(10, -12);
            filterThresholds[0][7][0] = Double.parseDouble(tfAlFiltF3.getText());
            filterThresholds[0][7][1] = Double.parseDouble(tfAhFiltF3.getText());
            filterThresholds[0][8][0] = Double.parseDouble(tfAlFiltD3.getText()) * Math.pow(10, -12);
            filterThresholds[0][8][1] = Double.parseDouble(tfAhFiltD3.getText()) * Math.pow(10, -12);
            filterThresholds[0][9][0] = Double.parseDouble(tfAlFiltFtrip.getText());
            filterThresholds[0][9][1] = Double.parseDouble(tfAhFiltFtrip.getText());
            filterThresholds[0][10][0] = Double.parseDouble(tfAlFiltTtrip.getText()) * Math.pow(10, -6);
            filterThresholds[0][10][1] = Double.parseDouble(tfAhFiltTtrip.getText()) * Math.pow(10, -6);
            filterThresholds[0][noparam][0] = Double.parseDouble(tfAlFiltChi2.getText());
            filterThresholds[0][noparam][1] = Double.parseDouble(tfAhFiltChi2.getText());
            filterThresholds[0][noparam + 1][1] = Double.parseDouble(tfAhFiltGrel.getText());
            filterThresholds[0][noparam + 1][0] = -1 * filterThresholds[0][noparam + 1][1];

            if (userThreshold[1]) {
                filterThresholds[2][0][0] = Double.parseDouble(tfClFiltN.getText());
                filterThresholds[2][0][1] = Double.parseDouble(tfChFiltN.getText());
                filterThresholds[2][1][0] = Double.parseDouble(tfClFiltD.getText()) * Math.pow(10, -12);
                filterThresholds[2][1][1] = Double.parseDouble(tfChFiltD.getText()) * Math.pow(10, -12);
                filterThresholds[2][2][0] = Double.parseDouble(tfClFiltvx.getText()) * Math.pow(10, -6);
                filterThresholds[2][2][1] = Double.parseDouble(tfChFiltvx.getText()) * Math.pow(10, -6);
                filterThresholds[2][3][0] = Double.parseDouble(tfClFiltvy.getText()) * Math.pow(10, -6);
                filterThresholds[2][3][1] = Double.parseDouble(tfChFiltvy.getText()) * Math.pow(10, -6);
                filterThresholds[2][4][0] = Double.parseDouble(tfClFiltG.getText());
                filterThresholds[2][4][1] = Double.parseDouble(tfChFiltG.getText());
                filterThresholds[2][5][0] = Double.parseDouble(tfClFiltF2.getText());
                filterThresholds[2][5][1] = Double.parseDouble(tfChFiltF2.getText());
                filterThresholds[2][6][0] = Double.parseDouble(tfClFiltD2.getText()) * Math.pow(10, -12);
                filterThresholds[2][6][1] = Double.parseDouble(tfChFiltD2.getText()) * Math.pow(10, -12);
                filterThresholds[2][7][0] = Double.parseDouble(tfClFiltF3.getText());
                filterThresholds[2][7][1] = Double.parseDouble(tfChFiltF3.getText());
                filterThresholds[2][8][0] = Double.parseDouble(tfClFiltD3.getText()) * Math.pow(10, -12);
                filterThresholds[2][8][1] = Double.parseDouble(tfChFiltD3.getText()) * Math.pow(10, -12);
                filterThresholds[2][7][0] = Double.parseDouble(tfClFiltFtrip.getText());
                filterThresholds[2][7][1] = Double.parseDouble(tfChFiltFtrip.getText());
                filterThresholds[2][8][0] = Double.parseDouble(tfClFiltTtrip.getText()) * Math.pow(10, -6);
                filterThresholds[2][8][1] = Double.parseDouble(tfChFiltTtrip.getText()) * Math.pow(10, -6);
                filterThresholds[2][noparam][0] = Double.parseDouble(tfClFiltChi2.getText());
                filterThresholds[2][noparam][1] = Double.parseDouble(tfChFiltChi2.getText());
                filterThresholds[2][noparam + 1][1] = Double.parseDouble(tfChFiltGrel.getText());
                filterThresholds[2][noparam + 1][0] = -1 * filterThresholds[2][noparam + 1][1];

                for (int p = 0; p < noparam + 2; p++) {
                    filterThresholds[1][p][0] = filterThresholds[0][p][0];
                    filterThresholds[1][p][1] = filterThresholds[0][p][1];
                }

                userThreshold[3] = rbtnCuseA.isSelected();
                userThreshold[4] = rbtnReplaceZero.isSelected();
            }
        } catch (NumberFormatException nfe) {
            IJ.showMessage("A value in the thresholds settings has an invalid format.");
            throw new NumberFormatException("Number format error.");
        }
    }

    ;

	public void SetThresholds() {	// sets the Thresholds and other settings to their default values in case no user defined values are vailable or they need to be reset
        double Min;
        double Max;
        double eps = 0.5 * Math.pow(10, -decformat);
        double eps2 = 0.5 * Math.pow(10, -decformat2);
        int nochannels = 3;		// number of correlation channels; 1 for FCS, 3 for DC-FCCS (for cross-correlation and autocorrelation in the 2 channels respectively)

        //if thresholds were set for DC-FCCS and switched to FCS or vice-versa, forget thresholds
        if (((userThreshold[1]) && (cbFitModel.getSelectedItem() != "DC-FCCS")) || ((!userThreshold[1]) && (cbFitModel.getSelectedItem() == "DC-FCCS"))) {
            userThreshold[0] = false;
        }

        //if user has not defined any thresholds, set the default values of all settings
        if (!userThreshold[0]) {
            userThreshold[1] = cbFitModel.getSelectedItem() == "DC-FCCS";
            userThreshold[2] = false;				//the thresholds have not been applied on given file
            userThreshold[3] = false;				//thresholds for ACFs and CCFs are independent
            userThreshold[4] = false;				//do not set q=0

            // if a paremeter map exists, find out min and max of fitted parameters in the current parameter map and use them as starting values of the thresholds
            if (impPara1 != null) {

                for (int m = 0; m < nochannels; m++) {							//loop over individual correlation channels and the cross-correlation
                    for (int p = 0; p < noparam; p++) {
                        Max = -Double.MAX_VALUE;
                        Min = Double.MAX_VALUE;
                        for (int x = 0; x < width; x++) {
                            for (int y = 0; y < height; y++) {
                                if (fitres[m][x][y][p] < Min && fitres[m][x][y][p] != Double.NaN) {
                                    Min = fitres[m][x][y][p];
                                }
                                if (fitres[m][x][y][p] > Max && fitres[m][x][y][p] != Double.NaN) {
                                    Max = fitres[m][x][y][p];
                                }
                                filterThresholds[m][p][0] = Min;
                                filterThresholds[m][p][1] = Max;
                            }
                        }
                    }

                    Max = -Double.MAX_VALUE;
                    Min = Double.MAX_VALUE;
                    for (int x = 0; x < width; x++) {
                        for (int y = 0; y < height; y++) {
                            if (chi2[m][x][y] < Min && chi2[m][x][y] != Double.NaN) {
                                Min = chi2[m][x][y];
                            }
                            if (chi2[m][x][y] > Max && chi2[m][x][y] != Double.NaN) {
                                Max = chi2[m][x][y];
                            }
                            filterThresholds[m][noparam][0] = Min;
                            filterThresholds[m][noparam][1] = Max;
                        }
                    }

                    Max = -Double.MAX_VALUE;
                    Min = Double.MAX_VALUE;
                    for (int x = 0; x < width; x++) {
                        for (int y = 0; y < height; y++) {
                            if (fitres[m][x][y][4] * fitres[m][x][y][0] < Min && fitres[m][x][y][4] != Double.NaN && fitres[m][x][y][0] != Double.NaN) {
                                Min = fitres[m][x][y][4] * fitres[m][x][y][0];
                            }
                            if (fitres[m][x][y][4] * fitres[m][x][y][0] > Max && fitres[m][x][y][4] != Double.NaN && fitres[m][x][y][0] != Double.NaN) {
                                Max = fitres[m][x][y][4] * fitres[m][x][y][0];
                            }
                            if (-1 * Min > Max) {
                                Max = Min;
                            }
                            filterThresholds[m][noparam + 1][0] = -1 * Max;
                            filterThresholds[m][noparam + 1][1] = Max;
                        }
                    }
                }
            }

            //unselect all parameters for filtering and in DC-FCCS mode set the same thresholds for both ACF channels
            for (int p = 0; p < noparam + 2; p++) {
                paramfilter[p] = false;
                if (userThreshold[1]) {
                    if (filterThresholds[0][p][0] > filterThresholds[1][p][0]) {
                        filterThresholds[0][p][0] = filterThresholds[1][p][0];
                    }
                    if (filterThresholds[0][p][1] < filterThresholds[1][p][1]) {
                        filterThresholds[0][p][1] = filterThresholds[1][p][1];
                    }
                }
            }
        }

        // set the values in the panel
        tfAlFiltN.setText(IJ.d2s(filterThresholds[0][0][0] - eps, decformat));
        tfAhFiltN.setText(IJ.d2s(filterThresholds[0][0][1] + eps, decformat));
        tfAlFiltD.setText(IJ.d2s(filterThresholds[0][1][0] * Math.pow(10, 12) - eps, decformat));
        tfAhFiltD.setText(IJ.d2s(filterThresholds[0][1][1] * Math.pow(10, 12) + eps, decformat));
        tfAlFiltvx.setText(IJ.d2s(filterThresholds[0][2][0] * Math.pow(10, 6) - eps, decformat));
        tfAhFiltvx.setText(IJ.d2s(filterThresholds[0][2][1] * Math.pow(10, 6) + eps, decformat));
        tfAlFiltvy.setText(IJ.d2s(filterThresholds[0][3][0] * Math.pow(10, 6) - eps, decformat));
        tfAhFiltvy.setText(IJ.d2s(filterThresholds[0][3][1] * Math.pow(10, 6) + eps, decformat));
        tfAlFiltG.setText(IJ.d2s(filterThresholds[0][4][0] - eps2, decformat2));
        tfAhFiltG.setText(IJ.d2s(filterThresholds[0][4][1] + eps2, decformat2));
        tfAhFiltGrel.setText(IJ.d2s(filterThresholds[0][10][1], decformat));
        tfAlFiltF2.setText(IJ.d2s(filterThresholds[0][5][0] - eps, decformat));
        tfAhFiltF2.setText(IJ.d2s(filterThresholds[0][5][1] + eps, decformat));
        tfAlFiltD2.setText(IJ.d2s(filterThresholds[0][6][0] * Math.pow(10, 12) - eps, decformat));
        tfAhFiltD2.setText(IJ.d2s(filterThresholds[0][6][1] * Math.pow(10, 12) + eps, decformat));
        tfAlFiltF3.setText(IJ.d2s(filterThresholds[0][7][0] - eps, decformat));
        tfAhFiltF3.setText(IJ.d2s(filterThresholds[0][7][1] + eps, decformat));
        tfAlFiltD3.setText(IJ.d2s(filterThresholds[0][8][0] * Math.pow(10, 12) - eps, decformat));
        tfAhFiltD3.setText(IJ.d2s(filterThresholds[0][8][1] * Math.pow(10, 12) + eps, decformat));
        tfAlFiltFtrip.setText(IJ.d2s(filterThresholds[0][9][0] - eps, decformat));
        tfAhFiltFtrip.setText(IJ.d2s(filterThresholds[0][9][1] + eps, decformat));
        tfAlFiltTtrip.setText(IJ.d2s(filterThresholds[0][10][0] * Math.pow(10, 6) - eps, decformat));
        tfAhFiltTtrip.setText(IJ.d2s(filterThresholds[0][10][1] * Math.pow(10, 6) + eps, decformat));
        tfAlFiltChi2.setText(IJ.d2s(filterThresholds[0][11][0] - eps, decformat));
        tfAhFiltChi2.setText(IJ.d2s(filterThresholds[0][11][1] + eps, decformat));
        tfClFiltN.setText(IJ.d2s(filterThresholds[2][0][0] - eps, decformat));
        tfChFiltN.setText(IJ.d2s(filterThresholds[2][0][1] + eps, decformat));
        tfClFiltD.setText(IJ.d2s(filterThresholds[2][1][0] * Math.pow(10, 12) - eps, decformat));
        tfChFiltD.setText(IJ.d2s(filterThresholds[2][1][1] * Math.pow(10, 12) + eps, decformat));
        tfClFiltvx.setText(IJ.d2s(filterThresholds[2][2][0] * Math.pow(10, 6) - eps, decformat));
        tfChFiltvx.setText(IJ.d2s(filterThresholds[2][2][1] * Math.pow(10, 6) + eps, decformat));
        tfClFiltvy.setText(IJ.d2s(filterThresholds[2][3][0] * Math.pow(10, 6) - eps, decformat));
        tfChFiltvy.setText(IJ.d2s(filterThresholds[2][3][1] * Math.pow(10, 6) + eps, decformat));
        tfClFiltG.setText(IJ.d2s(filterThresholds[2][4][0] - eps2, decformat2));
        tfChFiltG.setText(IJ.d2s(filterThresholds[2][4][1] + eps2, decformat2));
        tfChFiltGrel.setText(IJ.d2s(filterThresholds[2][10][1], decformat));
        tfClFiltF2.setText(IJ.d2s(filterThresholds[2][5][0] - eps, decformat));
        tfChFiltF2.setText(IJ.d2s(filterThresholds[2][5][1] + eps, decformat));
        tfClFiltD2.setText(IJ.d2s(filterThresholds[2][6][0] * Math.pow(10, 12) - eps, decformat));
        tfChFiltD2.setText(IJ.d2s(filterThresholds[2][6][1] * Math.pow(10, 12) + eps, decformat));
        tfClFiltF3.setText(IJ.d2s(filterThresholds[2][7][0] - eps, decformat));
        tfChFiltF3.setText(IJ.d2s(filterThresholds[2][7][1] + eps, decformat));
        tfClFiltD3.setText(IJ.d2s(filterThresholds[2][8][0] * Math.pow(10, 12) - eps, decformat));
        tfChFiltD3.setText(IJ.d2s(filterThresholds[2][8][1] * Math.pow(10, 12) + eps, decformat));
        tfClFiltFtrip.setText(IJ.d2s(filterThresholds[2][7][0] - eps, decformat));
        tfChFiltFtrip.setText(IJ.d2s(filterThresholds[2][7][1] + eps, decformat));
        tfClFiltTtrip.setText(IJ.d2s(filterThresholds[2][8][0] * Math.pow(10, 6) - eps, decformat));
        tfChFiltTtrip.setText(IJ.d2s(filterThresholds[2][8][1] * Math.pow(10, 6) + eps, decformat));
        tfClFiltChi2.setText(IJ.d2s(filterThresholds[2][9][0] - eps, decformat));
        tfChFiltChi2.setText(IJ.d2s(filterThresholds[2][9][1] + eps, decformat));

        if (paramfilter[0]) {
            rbtnFiltN.setSelected(true);
            tfAlFiltN.setEnabled(true);
            tfAhFiltN.setEnabled(true);
        } else {
            rbtnFiltN.setSelected(false);
            tfAlFiltN.setEnabled(false);
            tfAhFiltN.setEnabled(false);
        }

        if (paramfilter[1]) {
            rbtnFiltD.setSelected(true);
            tfAlFiltD.setEnabled(true);
            tfAhFiltD.setEnabled(true);
        } else {
            rbtnFiltD.setSelected(false);
            tfAlFiltD.setEnabled(false);
            tfAhFiltD.setEnabled(false);
        }

        if (paramfilter[2]) {
            rbtnFiltvx.setSelected(true);
            tfAlFiltvx.setEnabled(true);
            tfAhFiltvx.setEnabled(true);
        } else {
            rbtnFiltvx.setSelected(false);
            tfAlFiltvx.setEnabled(false);
            tfAhFiltvx.setEnabled(false);
        }

        if (paramfilter[3]) {
            rbtnFiltvy.setSelected(true);
            tfAlFiltvy.setEnabled(true);
            tfAhFiltvy.setEnabled(true);
        } else {
            rbtnFiltvy.setSelected(false);
            tfAlFiltvy.setEnabled(false);
            tfAhFiltvy.setEnabled(false);
        }

        tfAlFiltG.setEnabled(false);
        tfAhFiltG.setEnabled(false);
        if (paramfilter[4]) {
            rbtnFiltG.setSelected(true);
            if (!paramfilter[noparam + 1]) {
                tfAlFiltG.setEnabled(true);
                tfAhFiltG.setEnabled(true);
            } else {
                rbtnFiltG.setEnabled(false);
            }
        } else {
            rbtnFiltG.setSelected(false);
        }

        tfAhFiltGrel.setEnabled(false);
        if (paramfilter[4]) {
            rbtnFiltGrel.setEnabled(true);
            if (paramfilter[noparam + 1]) {
                rbtnFiltGrel.setSelected(true);
                tfAhFiltGrel.setEnabled(true);
            } else {
                rbtnFiltGrel.setSelected(false);
            }
        } else {
            rbtnFiltGrel.setEnabled(false);
            rbtnFiltGrel.setSelected(false);
        }

        if (paramfilter[5]) {
            rbtnFiltF2.setSelected(true);
            tfAlFiltF2.setEnabled(true);
            tfAhFiltF2.setEnabled(true);
        } else {
            rbtnFiltF2.setSelected(false);
            tfAlFiltF2.setEnabled(false);
            tfAhFiltF2.setEnabled(false);
        }

        if (paramfilter[6]) {
            rbtnFiltD2.setSelected(true);
            tfAlFiltD2.setEnabled(true);
            tfAhFiltD2.setEnabled(true);
        } else {
            rbtnFiltD2.setSelected(false);
            tfAlFiltD2.setEnabled(false);
            tfAhFiltD2.setEnabled(false);
        }

        if (paramfilter[7]) {
            rbtnFiltF3.setSelected(true);
            tfAlFiltF3.setEnabled(true);
            tfAhFiltF3.setEnabled(true);
        } else {
            rbtnFiltF3.setSelected(false);
            tfAlFiltF3.setEnabled(false);
            tfAhFiltF3.setEnabled(false);
        }

        if (paramfilter[8]) {
            rbtnFiltD3.setSelected(true);
            tfAlFiltD3.setEnabled(true);
            tfAhFiltD3.setEnabled(true);
        } else {
            rbtnFiltD3.setSelected(false);
            tfAlFiltD3.setEnabled(false);
            tfAhFiltD3.setEnabled(false);
        }

        if (paramfilter[9]) {
            rbtnFiltFtrip.setSelected(true);
            tfAlFiltFtrip.setEnabled(true);
            tfAhFiltFtrip.setEnabled(true);
        } else {
            rbtnFiltFtrip.setSelected(false);
            tfAlFiltFtrip.setEnabled(false);
            tfAhFiltFtrip.setEnabled(false);
        }

        if (paramfilter[10]) {
            rbtnFiltTtrip.setSelected(true);
            tfAlFiltTtrip.setEnabled(true);
            tfAhFiltTtrip.setEnabled(true);
        } else {
            rbtnFiltTtrip.setSelected(false);
            tfAlFiltTtrip.setEnabled(false);
            tfAhFiltTtrip.setEnabled(false);
        }

        if (paramfilter[noparam]) {
            rbtnFiltChi2.setSelected(true);
            tfAlFiltChi2.setEnabled(true);
            tfAhFiltChi2.setEnabled(true);
        } else {
            rbtnFiltChi2.setSelected(false);
            tfAlFiltChi2.setEnabled(false);
            tfAhFiltChi2.setEnabled(false);
        }

        tfClFiltN.setEnabled(false);
        tfChFiltN.setEnabled(false);
        tfClFiltD.setEnabled(false);
        tfChFiltD.setEnabled(false);
        tfClFiltvx.setEnabled(false);
        tfChFiltvx.setEnabled(false);
        tfClFiltvy.setEnabled(false);
        tfChFiltvy.setEnabled(false);
        tfClFiltG.setEnabled(false);
        tfChFiltG.setEnabled(false);
        tfChFiltGrel.setEnabled(false);
        tfClFiltF2.setEnabled(false);
        tfChFiltF2.setEnabled(false);
        tfClFiltD2.setEnabled(false);
        tfChFiltD2.setEnabled(false);
        tfClFiltF3.setEnabled(false);
        tfChFiltF3.setEnabled(false);
        tfClFiltD3.setEnabled(false);
        tfChFiltD3.setEnabled(false);
        tfClFiltFtrip.setEnabled(false);
        tfChFiltFtrip.setEnabled(false);
        tfClFiltTtrip.setEnabled(false);
        tfChFiltTtrip.setEnabled(false);
        tfClFiltChi2.setEnabled(false);
        tfChFiltChi2.setEnabled(false);

        if (userThreshold[1] && !userThreshold[3]) {
            if (paramfilter[0]) {
                tfClFiltN.setEnabled(true);
                tfChFiltN.setEnabled(true);
            }
            if (paramfilter[1]) {
                tfClFiltD.setEnabled(true);
                tfChFiltD.setEnabled(true);
            }
            if (paramfilter[2]) {
                tfClFiltvx.setEnabled(true);
                tfChFiltvx.setEnabled(true);
            }
            if (paramfilter[3]) {
                tfClFiltvy.setEnabled(true);
                tfChFiltvy.setEnabled(true);
            }
            if (!paramfilter[noparam + 1] && paramfilter[4]) {
                tfClFiltG.setEnabled(true);
                tfChFiltG.setEnabled(true);
            }
            if (paramfilter[4] && paramfilter[noparam + 1]) {
                tfChFiltGrel.setEnabled(true);
            }
            if (paramfilter[5]) {
                tfClFiltF2.setEnabled(true);
                tfChFiltF2.setEnabled(true);
            }
            if (paramfilter[6]) {
                tfClFiltD2.setEnabled(true);
                tfChFiltD2.setEnabled(true);
            }
            if (paramfilter[7]) {
                tfClFiltF3.setEnabled(true);
                tfChFiltF3.setEnabled(true);
            }
            if (paramfilter[8]) {
                tfClFiltD3.setEnabled(true);
                tfChFiltD3.setEnabled(true);
            }
            if (paramfilter[9]) {
                tfClFiltFtrip.setEnabled(true);
                tfChFiltFtrip.setEnabled(true);
            }
            if (paramfilter[10]) {
                tfClFiltTtrip.setEnabled(true);
                tfChFiltTtrip.setEnabled(true);
            }
            if (paramfilter[noparam]) {
                tfClFiltChi2.setEnabled(true);
                tfChFiltChi2.setEnabled(true);
            }
        }

        if (userThreshold[1]) {
            rbtnCuseA.setSelected(userThreshold[3]);
            rbtnReplaceZero.setSelected(userThreshold[4]);
            rbtnReplaceZero.setEnabled(true);
            rbtnCuseA.setEnabled(true);
        } else {
            rbtnCuseA.setSelected(false);
            rbtnReplaceZero.setSelected(false);
            rbtnReplaceZero.setEnabled(false);
            rbtnCuseA.setEnabled(false);
        }
    }
    ;

	ActionListener btnFilterPressed = (ActionEvent event) -> {
        if (impPara1 != null) {
            updateParaImp();
        } else {
            JOptionPane.showMessageDialog(null, "No parameter map available. Perform a correlation or load an experiment.");
        }
    };

    ActionListener btnResetPressed = (ActionEvent event) -> {
        userThreshold[0] = false;
        SetThresholds();
        if (impPara1 != null) {
            updateParaImp();
        }
    };

    ItemListener rbtnFiltGChanged = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            rbtnFiltGrel.setEnabled(true);
            tfAlFiltG.setEnabled(true);
            tfAhFiltG.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltG.setEnabled(true);
                tfChFiltG.setEnabled(true);
            }
        } else {
            rbtnFiltGrel.setEnabled(false);
            tfAlFiltG.setEnabled(false);
            tfAhFiltG.setEnabled(false);
            if (userThreshold[1]) {
                tfClFiltG.setEnabled(false);
                tfChFiltG.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltGrelChanged = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAhFiltGrel.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfChFiltGrel.setEnabled(true);
            }
            rbtnFiltG.setEnabled(false);
            tfAlFiltG.setEnabled(false);
            tfAhFiltG.setEnabled(false);
            if (userThreshold[1]) {
                tfClFiltG.setEnabled(false);
                tfChFiltG.setEnabled(false);
            }
        } else {
            tfAhFiltGrel.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfChFiltGrel.setEnabled(false);
            }
            rbtnFiltG.setEnabled(true);
            if (rbtnFiltG.isSelected()) {
                tfAlFiltG.setEnabled(true);
                tfAhFiltG.setEnabled(true);
                if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                    tfClFiltG.setEnabled(true);
                    tfChFiltG.setEnabled(true);
                }
            }
        }
    };

    ItemListener rbtnFiltNChanged = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltN.setEnabled(true);
            tfAhFiltN.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltN.setEnabled(true);
                tfChFiltN.setEnabled(true);
            }
        } else {
            tfAlFiltN.setEnabled(false);
            tfAhFiltN.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltN.setEnabled(false);
                tfChFiltN.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltDChanged = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltD.setEnabled(true);
            tfAhFiltD.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltD.setEnabled(true);
                tfChFiltD.setEnabled(true);
            }
        } else {
            tfAlFiltD.setEnabled(false);
            tfAhFiltD.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltD.setEnabled(false);
                tfChFiltD.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltvxChanged = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltvx.setEnabled(true);
            tfAhFiltvx.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltvx.setEnabled(true);
                tfChFiltvx.setEnabled(true);
            }
        } else {
            tfAlFiltvx.setEnabled(false);
            tfAhFiltvx.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltvx.setEnabled(false);
                tfChFiltvx.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltvyChanged = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltvy.setEnabled(true);
            tfAhFiltvy.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltvy.setEnabled(true);
                tfChFiltvy.setEnabled(true);
            }
        } else {
            tfAlFiltvy.setEnabled(false);
            tfAhFiltvy.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltvy.setEnabled(false);
                tfChFiltvy.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltF2Changed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltF2.setEnabled(true);
            tfAhFiltF2.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltF2.setEnabled(true);
                tfChFiltF2.setEnabled(true);
            }
        } else {
            tfAlFiltF2.setEnabled(false);
            tfAhFiltF2.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltF2.setEnabled(false);
                tfChFiltF2.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltD2Changed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltD2.setEnabled(true);
            tfAhFiltD2.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltD2.setEnabled(true);
                tfChFiltD2.setEnabled(true);
            }
        } else {
            tfAlFiltD2.setEnabled(false);
            tfAhFiltD2.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltD2.setEnabled(false);
                tfChFiltD2.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltF3Changed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltF3.setEnabled(true);
            tfAhFiltF3.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltF3.setEnabled(true);
                tfChFiltF3.setEnabled(true);
            }
        } else {
            tfAlFiltF3.setEnabled(false);
            tfAhFiltF3.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltF3.setEnabled(false);
                tfChFiltF3.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltD3Changed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltD3.setEnabled(true);
            tfAhFiltD3.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltD3.setEnabled(true);
                tfChFiltD3.setEnabled(true);
            }
        } else {
            tfAlFiltD3.setEnabled(false);
            tfAhFiltD3.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltD3.setEnabled(false);
                tfChFiltD3.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltFtripChanged = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltFtrip.setEnabled(true);
            tfAhFiltFtrip.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltFtrip.setEnabled(true);
                tfChFiltFtrip.setEnabled(true);
            }
        } else {
            tfAlFiltFtrip.setEnabled(false);
            tfAhFiltFtrip.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltFtrip.setEnabled(false);
                tfChFiltFtrip.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltTtripChanged = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltTtrip.setEnabled(true);
            tfAhFiltTtrip.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltTtrip.setEnabled(true);
                tfChFiltTtrip.setEnabled(true);
            }
        } else {
            tfAlFiltTtrip.setEnabled(false);
            tfAhFiltTtrip.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltTtrip.setEnabled(false);
                tfChFiltTtrip.setEnabled(false);
            }
        }
    };

    ItemListener rbtnFiltChi2Changed = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfAlFiltChi2.setEnabled(true);
            tfAhFiltChi2.setEnabled(true);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltChi2.setEnabled(true);
                tfChFiltChi2.setEnabled(true);
            }
        } else {
            tfAlFiltChi2.setEnabled(false);
            tfAhFiltChi2.setEnabled(false);
            if (userThreshold[1] && !rbtnCuseA.isSelected()) {
                tfClFiltChi2.setEnabled(false);
                tfChFiltChi2.setEnabled(false);
            }
        }
    };

    ItemListener rbtnCuseAChanged = (ItemEvent ev) -> {
        if (ev.getStateChange() == ItemEvent.SELECTED) {
            tfClFiltN.setEnabled(false);
            tfChFiltN.setEnabled(false);
            tfClFiltD.setEnabled(false);
            tfChFiltD.setEnabled(false);
            tfClFiltvx.setEnabled(false);
            tfChFiltvx.setEnabled(false);
            tfClFiltvy.setEnabled(false);
            tfChFiltvy.setEnabled(false);
            tfClFiltG.setEnabled(false);
            tfChFiltG.setEnabled(false);
            tfClFiltF2.setEnabled(false);
            tfChFiltF2.setEnabled(false);
            tfClFiltD2.setEnabled(false);
            tfChFiltD2.setEnabled(false);
            tfClFiltF3.setEnabled(false);
            tfChFiltF3.setEnabled(false);
            tfClFiltD3.setEnabled(false);
            tfChFiltD3.setEnabled(false);
            tfClFiltFtrip.setEnabled(false);
            tfChFiltFtrip.setEnabled(false);
            tfClFiltTtrip.setEnabled(false);
            tfChFiltTtrip.setEnabled(false);
            tfClFiltChi2.setEnabled(false);
            tfChFiltChi2.setEnabled(false);
            tfChFiltGrel.setEnabled(false);
        } else {
            if (rbtnFiltN.isSelected()) {
                tfClFiltN.setEnabled(true);
            }
            if (rbtnFiltN.isSelected()) {
                tfChFiltN.setEnabled(true);
            }
            if (rbtnFiltD.isSelected()) {
                tfClFiltD.setEnabled(true);
            }
            if (rbtnFiltD.isSelected()) {
                tfChFiltD.setEnabled(true);
            }
            if (rbtnFiltvx.isSelected()) {
                tfClFiltvx.setEnabled(true);
            }
            if (rbtnFiltvx.isSelected()) {
                tfChFiltvx.setEnabled(true);
            }
            if (rbtnFiltvy.isSelected()) {
                tfClFiltvy.setEnabled(true);
            }
            if (rbtnFiltvy.isSelected()) {
                tfChFiltvy.setEnabled(true);
            }
            if (rbtnFiltG.isSelected() && !rbtnFiltGrel.isSelected()) {
                tfClFiltG.setEnabled(true);
            }
            if (rbtnFiltG.isSelected() && !rbtnFiltGrel.isSelected()) {
                tfChFiltG.setEnabled(true);
            }
            if (rbtnFiltF2.isSelected()) {
                tfClFiltF2.setEnabled(true);
            }
            if (rbtnFiltF2.isSelected()) {
                tfChFiltF2.setEnabled(true);
            }
            if (rbtnFiltD2.isSelected()) {
                tfClFiltD2.setEnabled(true);
            }
            if (rbtnFiltD2.isSelected()) {
                tfChFiltD2.setEnabled(true);
            }
            if (rbtnFiltF3.isSelected()) {
                tfClFiltF3.setEnabled(true);
            }
            if (rbtnFiltF3.isSelected()) {
                tfChFiltF3.setEnabled(true);
            }
            if (rbtnFiltD3.isSelected()) {
                tfClFiltD3.setEnabled(true);
            }
            if (rbtnFiltD3.isSelected()) {
                tfChFiltD3.setEnabled(true);
            }
            if (rbtnFiltFtrip.isSelected()) {
                tfClFiltFtrip.setEnabled(true);
            }
            if (rbtnFiltTtrip.isSelected()) {
                tfChFiltTtrip.setEnabled(true);
            }
            if (rbtnFiltChi2.isSelected()) {
                tfClFiltChi2.setEnabled(true);
            }
            if (rbtnFiltChi2.isSelected()) {
                tfChFiltChi2.setEnabled(true);
            }
            if (rbtnFiltG.isSelected() && rbtnFiltGrel.isSelected()) {
                tfChFiltGrel.setEnabled(true);
            }
        }
    };

    /*
	 *  The following are control procedures for the plugin
	 *  createImFCSResultsTable(): create a results JTable with all numerical values inside; this is the same as will be saved as a spreadsheet by the plugin
	 *  writeExperiment(File dir, String $exception): Save essential experiment data (only fit results and window parameters)
	 *  readExperiment(File dir, String $exception): Reads saved experimetn data and reconstitutes parameter and histogram window
	 *  bringToFront(): brings all window of this instance of the plugin to the front
	 *  exitImFCS(): exits the plugin and cleans up windows
	 *  
     */
    // Create a Results Jtable
    public void createImFCSResultsTable() {
        int nov = width * height;
        int dlnov = diffLawMapwidth * diffLawMapheight;
        Object[][] paravalue = new Object[37][2];
        Object[][] thresholds = new Object[noparam + 4][6];
        Object[] thresCol = new Object[6];
        Object[] columnNames = new Object[2];

        // clear frame if it exists to ensure proper updating
        if (resframe != null) {
            resframe.getContentPane().removeAll();
        }

        // parameters to be shown in results table
        columnNames[0] = "Description";
        columnNames[1] = "Value";

        paravalue[0][0] = "First frame";
        paravalue[1][0] = "Last frame";
        paravalue[2][0] = "Frame time";
        paravalue[3][0] = "Binning X";
        paravalue[4][0] = "Binning Y";
        paravalue[5][0] = "CF X distance";
        paravalue[6][0] = "CF Y distance";
        paravalue[7][0] = "Correlator P";
        paravalue[8][0] = "Correlator Q";
        paravalue[9][0] = "Fit model";
        paravalue[10][0] = "FCCS Display";
        paravalue[11][0] = "PixelSize";
        paravalue[12][0] = "Overlap";
        paravalue[13][0] = "Magnification";
        paravalue[14][0] = "NA";
        paravalue[15][0] = "Em Lambda";
        paravalue[16][0] = "Em Lambda2";
        paravalue[17][0] = "Sigma";
        paravalue[18][0] = "SigmaZ";
        paravalue[19][0] = "Sigma2";
        paravalue[20][0] = "SigmaZ2";
        paravalue[21][0] = "Background";
        paravalue[22][0] = "Background2";
        paravalue[23][0] = "Background File";
        paravalue[24][0] = "Bleach Correction";
        paravalue[25][0] = "Sliding Window Length";
        paravalue[26][0] = "Polynomial Order";
        paravalue[27][0] = "Filter";
        paravalue[28][0] = "filterUL";
        paravalue[29][0] = "filterLL";
        paravalue[30][0] = "Fit";
        paravalue[31][0] = "Threshold";
        paravalue[32][0] = "N&B calculated";
        paravalue[33][0] = "N&B corrected by slope";
        paravalue[34][0] = "N&B slope";

        for (int i = 0; i < noSettings; i++) {
            paravalue[i][1] = panelSettings[i];
        }
        paravalue[noSettings][1] = Boolean.toString(userThreshold[2]);
        paravalue[noSettings + 1][1] = Boolean.toString(NBperformed);
        paravalue[noSettings + 3][1] = Boolean.toString(NBcorrected);
        paravalue[noSettings + 4][1] = Double.toString(NBslope);

        JTable table = new JTable(paravalue, columnNames);

        thresCol[0] = "Parameter";
        thresCol[1] = "Filter";
        thresCol[2] = "Min 1";
        thresCol[3] = "Max 1";
        thresCol[4] = "Min CCF";
        thresCol[5] = "Max CCF";

        for (int p = 0; p < noparam + 1; p++) {
            thresholds[p][0] = $param[p];
            thresholds[p][1] = Boolean.toString(paramfilter[p]);
            thresholds[p][2] = filterThresholds[0][p][0];
            thresholds[p][3] = filterThresholds[0][p][1];
            if (userThreshold[1] && !userThreshold[3]) {
                thresholds[p][4] = filterThresholds[2][p][0];
                thresholds[p][5] = filterThresholds[2][p][1];
            }
        }

        thresholds[noparam + 1][0] = "G relative";
        thresholds[noparam + 1][1] = Boolean.toString(paramfilter[noparam + 1]);
        thresholds[noparam + 1][2] = filterThresholds[0][noparam + 1][0];
        thresholds[noparam + 1][3] = filterThresholds[0][noparam + 1][1];
        if (userThreshold[1] && !userThreshold[2]) {
            thresholds[noparam + 1][4] = filterThresholds[2][noparam + 1][0];
            thresholds[noparam + 1][5] = filterThresholds[2][noparam + 1][1];
        }
        if (userThreshold[1]) {
            thresholds[noparam + 2][0] = "CCF thresholds same as ACF";
            thresholds[noparam + 3][0] = "q=0 when CCF invalid and ACF valid";
            thresholds[noparam + 2][1] = Boolean.toString(userThreshold[3]);
            thresholds[noparam + 3][1] = Boolean.toString(userThreshold[4]);
        }

        JTable thresholdtable = new JTable(thresholds, thresCol);

        // Objects that contain array values for creation of tables
        Object[][][] acfvalue = new Object[3][chanum][nov];
        Object[][][] sdvalue = new Object[3][chanum][nov];
        Object[][][] fitacfvalue = new Object[3][chanum][nov];
        Object[][][] resvalue = new Object[3][chanum][nov];
        Object[][][] fitresvalue = new Object[3][noparam + 4][nov + 1];
        Object[][] tauvalue = new Object[chanum][3];
        Object[][] diffLawvalue = new Object[difflawallbin][7];
        Object[][] diffLawMapDatavalue = new Object[difflawmapbin][dlnov + 1];
        Object[][] diffLawMapSDvalue = new Object[difflawmapbin][dlnov + 1];
        Object[][] diffLawMapFitvalue = new Object[2][dlnov + 1];
        Object[][] psfDatavalue = new Object[psfmaxbin][3 * numofpsf];
        Object[] acfcolnames = new Object[nov];
        Object[] fitrescolnames = new Object[nov + 1];
        Object[] taucolnames = new Object[3];
        Object[] diffLawcolnames = new Object[7];
        Object[] diffLawMapcolnames = new Object[dlnov + 1];
        Object[] diffLawMapFitcolnames = new Object[dlnov + 1];
        Object[] psfDatacolnames = new Object[3 * numofpsf];
        Object[][][] msdvalue = new Object[3][chanum][nov];
        Object[][] dccfvalue = new Object[dccfMax][nov + 1];
        Object[] dccfcolnames = new Object[nov + 1];
        Object[][] NBvalue = new Object[4][nov + 1];
        Object[] NBcolnames = new Object[nov + 1];

        // add time
        taucolnames[0] = "S/N";
        taucolnames[1] = "Lagtime";
        taucolnames[2] = "Bin width";
        for (int i = 0; i < chanum; i++) {
            tauvalue[i][0] = i;
            tauvalue[i][1] = lagtime[i];
            tauvalue[i][2] = samp[i];
        }

        // add correlation functions and their data and msd 
        for (int k = 0; k < 3; k++) {
            for (int h = 0; h < height; h++) {
                for (int w = 0; w < width; w++) {
                    acfcolnames[h * width + w] = Integer.toString(h) + "x" + Integer.toString(w);
                    for (int i = 0; i < chanum; i++) {
                        acfvalue[k][i][h * width + w] = acf[k][w][h][i];
                        sdvalue[k][i][h * width + w] = sdacf[k][w][h][i];
                        fitacfvalue[k][i][h * width + w] = fitacf[k][w][h][i];
                        resvalue[k][i][h * width + w] = res[k][w][h][i];
                        msdvalue[k][i][h * width + w] = msd[k][w][h][i];
                    }
                }
            }
        }

        // add fit results
        fitrescolnames[0] = "Parameter";
        for (int k = 0; k < 3; k++) {
            for (int q = 0; q < noparam + 3; q++) {
                fitresvalue[k][q][0] = $param[q];
                if (k == 2) {
                    fitresvalue[k][noparam + 3][0] = $param[noparam + 3];
                }
            }
            for (int h = 0; h < height; h++) {
                for (int w = 0; w < width; w++) {
                    fitrescolnames[h * width + w + 1] = Integer.toString(h) + "x" + Integer.toString(w);
                    for (int i = 0; i < noparam; i++) {
                        fitresvalue[k][i][h * width + w + 1] = fitres[k][w][h][i];
                    }
                    fitresvalue[k][noparam][h * width + w + 1] = chi2[k][w][h];
                    fitresvalue[k][noparam + 1][h * width + w + 1] = blocked[k][w][h];
                    fitresvalue[k][noparam + 2][h * width + w + 1] = pixvalid[k][w][h];
                    if (k == 2) {
                        fitresvalue[k][noparam + 3][h * width + w + 1] = CCFq[w][h];
                    }
                }
            }
        }

        // add diffLaw data
        diffLawcolnames[0] = "Aeff";
        diffLawcolnames[1] = "Time";
        diffLawcolnames[2] = "SD";
        diffLawcolnames[3] = "intercept";
        diffLawcolnames[4] = "slope";
        diffLawcolnames[5] = "fit start";
        diffLawcolnames[6] = "fit end";

        diffLawvalue[0][3] = difflawfit[0];
        diffLawvalue[0][4] = difflawfit[1];
        diffLawvalue[0][5] = diffLawFitLim[0];
        diffLawvalue[0][6] = diffLawFitLim[1];
        for (int i = 0; i < difflawallbin; i++) {
            diffLawvalue[i][0] = difflaw[0][i];
            diffLawvalue[i][1] = difflaw[1][i];
            diffLawvalue[i][2] = difflaw[2][i];
        }

        // add Diff Law Map and SD Data
        diffLawMapcolnames[0] = "Obs Vol";
        for (int k = 0; k < difflawmapbin; k++) {
            diffLawMapDatavalue[k][0] = difflawarray[0][0][0][k];
            diffLawMapSDvalue[k][0] = difflawarray[0][0][0][k];
        }
        for (int h = 0; h < diffLawMapheight; h++) {
            for (int w = 0; w < diffLawMapwidth; w++) {
                diffLawMapcolnames[h * diffLawMapwidth + w + 1] = Integer.toString(h) + "x" + Integer.toString(w);
                for (int k = 0; k < difflawmapbin; k++) {
                    diffLawMapDatavalue[k][h * diffLawMapwidth + w + 1] = difflawarray[w][h][1][k];
                    diffLawMapSDvalue[k][h * diffLawMapwidth + w + 1] = difflawarray[w][h][2][k];
                }
            }
        }

        // add Diff Law Map Fit
        diffLawMapFitcolnames[0] = "Parameters";
        diffLawMapFitvalue[0][0] = "Intercept";
        diffLawMapFitvalue[1][0] = "Slope";
        for (int h = 0; h < diffLawMapheight; h++) {
            for (int w = 0; w < diffLawMapwidth; w++) {
                diffLawMapFitcolnames[h * diffLawMapwidth + w + 1] = Integer.toString(h) + "x" + Integer.toString(w);
                diffLawMapFitvalue[0][h * diffLawMapwidth + w + 1] = diffLawFitMap[w][h][0];
                diffLawMapFitvalue[1][h * diffLawMapwidth + w + 1] = diffLawFitMap[w][h][1];
            }
        }

        // add PSF data
        for (int i = 0; i < numofpsf; i++) {
            psfDatacolnames[0 + i * 3] = "Bin";
            psfDatacolnames[1 + i * 3] = "D (bin " + Integer.toString(i + 1) + ")";
            psfDatacolnames[2 + i * 3] = "SD (bin " + Integer.toString(i + 1) + ")";
        }

        for (int i = 0; i < numofpsf; i++) {
            for (int q = 0; q < 3; q++) {
                for (int k = 0; k < psfmaxbin; k++) {
                    psfDatavalue[k][i * 3] = psfData[i][0][k];
                    psfDatavalue[k][i * 3 + 1] = psfData[i][1][k];
                    psfDatavalue[k][i * 3 + 2] = psfData[i][2][k];
                }
            }
        }

        // add dCCF data
        dccfcolnames[0] = "Direction";
        dccfvalue[0][0] = "horizontal";
        dccfvalue[1][0] = "vertical";
        dccfvalue[2][0] = "diagonal /";
        dccfvalue[3][0] = "diagonal \\";

        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                dccfcolnames[h * width + w + 1] = Integer.toString(h) + "x" + Integer.toString(w);
                for (int i = 0; i < dccfMax; i++) {
                    dccfvalue[i][h * width + w + 1] = dccf[i][w][h];
                }
            }
        }

        // add N&B data
        NBcolnames[0] = " ";
        NBvalue[0][0] = "Number";
        NBvalue[1][0] = "Brightness";
        NBvalue[2][0] = "Num (corrected)";
        NBvalue[3][0] = "Epsilon (corrected)";

        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                NBcolnames[h * width + w + 1] = Integer.toString(h) + "x" + Integer.toString(w);
                NBvalue[0][h * width + w + 1] = NBN[w][h];
                NBvalue[1][h * width + w + 1] = NBB[w][h];
                NBvalue[2][h * width + w + 1] = NBNum[w][h];
                NBvalue[3][h * width + w + 1] = NBEps[w][h];
            }
        }

        // create the tables from arrays; use AUTO_RESIZE_OFF to ensure scroll bar is shown
        JTable tautable = new JTable(tauvalue, taucolnames);
        tautable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        tautable.setCellSelectionEnabled(true);

        JTable acftable0 = new JTable(acfvalue[0], acfcolnames);
        acftable0.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        acftable0.setCellSelectionEnabled(true);
        JTable sdtable0 = new JTable(sdvalue[0], acfcolnames);
        sdtable0.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        sdtable0.setCellSelectionEnabled(true);
        JTable fitacftable0 = new JTable(fitacfvalue[0], acfcolnames);
        fitacftable0.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        fitacftable0.setCellSelectionEnabled(true);
        JTable restable0 = new JTable(resvalue[0], acfcolnames);
        restable0.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        restable0.setCellSelectionEnabled(true);

        JTable acftable1 = new JTable(acfvalue[1], acfcolnames);
        acftable1.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        acftable1.setCellSelectionEnabled(true);
        JTable sdtable1 = new JTable(sdvalue[1], acfcolnames);
        sdtable1.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        sdtable1.setCellSelectionEnabled(true);
        JTable fitacftable1 = new JTable(fitacfvalue[1], acfcolnames);
        fitacftable1.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        fitacftable1.setCellSelectionEnabled(true);
        JTable restable1 = new JTable(resvalue[1], acfcolnames);
        restable1.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        restable1.setCellSelectionEnabled(true);

        JTable acftable2 = new JTable(acfvalue[2], acfcolnames);
        acftable2.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        acftable2.setCellSelectionEnabled(true);
        JTable sdtable2 = new JTable(sdvalue[2], acfcolnames);
        sdtable2.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        sdtable2.setCellSelectionEnabled(true);
        JTable fitacftable2 = new JTable(fitacfvalue[2], acfcolnames);
        fitacftable2.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        fitacftable2.setCellSelectionEnabled(true);
        JTable restable2 = new JTable(resvalue[2], acfcolnames);
        restable2.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        restable2.setCellSelectionEnabled(true);

        JTable fitrestable0 = new JTable(fitresvalue[0], fitrescolnames);
        fitrestable0.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        fitrestable0.setCellSelectionEnabled(true);
        JTable fitrestable1 = new JTable(fitresvalue[1], fitrescolnames);
        fitrestable1.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        fitrestable1.setCellSelectionEnabled(true);
        JTable fitrestable2 = new JTable(fitresvalue[2], fitrescolnames);
        fitrestable2.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        fitrestable2.setCellSelectionEnabled(true);

        JTable diffLawtable = new JTable(diffLawvalue, diffLawcolnames);
        diffLawtable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        diffLawtable.setCellSelectionEnabled(true);

        JTable diffLawMapDatatable = new JTable(diffLawMapDatavalue, diffLawMapcolnames);
        diffLawMapDatatable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        diffLawMapDatatable.setCellSelectionEnabled(true);

        JTable diffLawMapSDtable = new JTable(diffLawMapSDvalue, diffLawMapcolnames);
        diffLawMapSDtable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        diffLawMapSDtable.setCellSelectionEnabled(true);

        JTable diffLawMapFittable = new JTable(diffLawMapFitvalue, diffLawMapFitcolnames);
        diffLawMapFittable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        diffLawMapFittable.setCellSelectionEnabled(true);

        JTable psfDatatable = new JTable(psfDatavalue, psfDatacolnames);
        psfDatatable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        psfDatatable.setCellSelectionEnabled(true);

        JTable msdtable0 = new JTable(msdvalue[0], acfcolnames);
        msdtable0.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        msdtable0.setCellSelectionEnabled(true);
        JTable msdtable1 = new JTable(msdvalue[1], acfcolnames);
        msdtable1.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        msdtable1.setCellSelectionEnabled(true);
        JTable msdtable2 = new JTable(msdvalue[2], acfcolnames);
        msdtable2.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        msdtable2.setCellSelectionEnabled(true);

        JTable dccftable = new JTable(dccfvalue, dccfcolnames);
        dccftable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        dccftable.setCellSelectionEnabled(true);

        JTable NBtable = new JTable(NBvalue, NBcolnames);
        NBtable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        NBtable.setCellSelectionEnabled(true);

        // create and display Result Table
        resframe.setFocusable(true);
        resframe.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        resframe.setLocation(new Point(resTabPosX + 200, resTabPosY + 200));
        resframe.setSize(new Dimension(resTabDimX, resTabDimY));
        resframe.setResizable(true);

        JTabbedPane tabbedPane = new JTabbedPane(JTabbedPane.TOP, JTabbedPane.SCROLL_TAB_LAYOUT);
        tabbedPane.add(new JScrollPane(table), "Panel");
        tabbedPane.add(new JScrollPane(thresholdtable), "Thresholds");
        tabbedPane.add(new JScrollPane(tautable), "Time");

        tabbedPane.add(new JScrollPane(acftable0), "ACF1");
        tabbedPane.add(new JScrollPane(sdtable0), "SD_ACF1");
        tabbedPane.add(new JScrollPane(fitacftable0), "Fit_ACF1");
        tabbedPane.add(new JScrollPane(restable0), "Res_ACF1");
        tabbedPane.add(new JScrollPane(fitrestable0), "Fit_Parameters_ACF1");

        tabbedPane.add(new JScrollPane(acftable1), "ACF2");
        tabbedPane.add(new JScrollPane(sdtable1), "SD_ACF2");
        tabbedPane.add(new JScrollPane(fitacftable1), "Fit_ACF2");
        tabbedPane.add(new JScrollPane(restable1), "Res_ACF2");
        tabbedPane.add(new JScrollPane(fitrestable1), "Fit_Parameters_ACF2");

        tabbedPane.add(new JScrollPane(acftable2), "CCF");
        tabbedPane.add(new JScrollPane(sdtable2), "SD_CCF");
        tabbedPane.add(new JScrollPane(fitacftable2), "Fit_CCF");
        tabbedPane.add(new JScrollPane(restable2), "Res_CCF");
        tabbedPane.add(new JScrollPane(fitrestable2), "Fit_Parameters_CCF");

        tabbedPane.add(new JScrollPane(diffLawtable), "Diffusion Law");
        tabbedPane.add(new JScrollPane(diffLawMapDatatable), "Diffusion Law Map");
        tabbedPane.add(new JScrollPane(diffLawMapSDtable), "Diffusion Law Map SD");
        tabbedPane.add(new JScrollPane(diffLawMapFittable), "Diffusion Law Map Fit");

        tabbedPane.add(new JScrollPane(psfDatatable), "PSF data");

        tabbedPane.add(new JScrollPane(msdtable0), "MSD of ACF1");
        tabbedPane.add(new JScrollPane(msdtable1), "MSD of ACF2");
        tabbedPane.add(new JScrollPane(msdtable2), "MSD of CCF");

        tabbedPane.add(new JScrollPane(dccftable), "dCCF");

        tabbedPane.add(new JScrollPane(NBtable), "N&B");

        resframe.add(tabbedPane, BorderLayout.CENTER);

        resframe.pack();
        resframe.setVisible(true);
    }

    // IO options; the following procedures write and read experimental results
    public void writeExperiment(File sfile, String $exception, boolean showlog) {
        int nov = width * height;
        // There is a limit of 1,048,576 rows and 16,384 (corresponds to 128 x 128) columns in xlsx file as of 2019.
        if (nov > 16384 - 1) { // minus 1 because there is an additional (1st) column in Fit Parameters tab that list respective parameters 
            IJ.log("Unable to save the data. Output data is larger than permissible 16,384 columns in .xlsx file.");
            return;
        }

        int t;
        File file;
        int nofpv = 81;
        String[] paramname = new String[nofpv];
        String[] paramsave = new String[nofpv];

        String $sfile = sfile.toString();
        int dotind = $sfile.lastIndexOf('.');
        if (dotind != -1) {
            $sfile = $sfile.substring(0, dotind);
        }
        file = new File($sfile + ".xlsx");

        if (file.exists()) {
            if (!fileAlreadyExistsDialog()) {
                return;
            }
        }

        XSSFWorkbook wb = new XSSFWorkbook();
        Sheet paraSheet = wb.createSheet("Panel Parameters");
        Sheet thresholdSheet = wb.createSheet("Threshold Settings");
        Sheet tauSheet = wb.createSheet("lagtime");

        Sheet acf0Sheet = wb.createSheet("ACF1");
        Sheet sd0Sheet = wb.createSheet("SD (ACF1)");
        Sheet fitacf0Sheet = wb.createSheet("Fit functions (ACF1)");
        Sheet res0Sheet = wb.createSheet("Residuals (ACF1)");
        Sheet fitres0Sheet = wb.createSheet("Fit Parameters (ACF1)");

        Sheet acf1Sheet = wb.createSheet("ACF2");
        Sheet sd1Sheet = wb.createSheet("SD (ACF2)");
        Sheet fitacf1Sheet = wb.createSheet("Fit functions (ACF2)");
        Sheet res1Sheet = wb.createSheet("Residuals (ACF2)");
        Sheet fitres1Sheet = wb.createSheet("Fit Parameters (ACF2)");

        Sheet acf2Sheet = wb.createSheet("CCFs");
        Sheet sd2Sheet = wb.createSheet("SD (CCF)");
        Sheet fitacf2Sheet = wb.createSheet("Fit functions (CCF)");
        Sheet res2Sheet = wb.createSheet("Residuals (CCF)");
        Sheet fitres2Sheet = wb.createSheet("Fit Parameters (CCF)");

        Sheet diffLawSheet = wb.createSheet("Diffusion Law");
        Sheet diffLawMapDataSheet = wb.createSheet("DiffLaw Map Data");
        Sheet diffLawMapSDSheet = wb.createSheet("DiffLaw Map SD");
        Sheet diffLawMapFitSheet = wb.createSheet("DiffLaw Map Fit");

        Sheet psfDataSheet = wb.createSheet("PSF");

        Sheet msd0Sheet = wb.createSheet("MSD of ACF1");
        Sheet msd1Sheet = wb.createSheet("MSD of ACF2");
        Sheet msd2Sheet = wb.createSheet("MSD of CCF");

        Sheet dCCFDataSheet = wb.createSheet("dCCF");

        Sheet NBDataSheet = wb.createSheet("N&B");

        Sheet ReconParaSheet = wb.createSheet("Reconstruction Control Parameters");

        Row row;
        Row row0;
        Row row1;
        Row row2;
        Row row3;
        Row row4;
        Row row5;

        // write Imaging FCS Panel parameters
        t = 0;
        paramname[t++] = "Image width";
        paramname[t++] = "Image height";
        paramname[t++] = "First frame";
        paramname[t++] = "Last frame";
        paramname[t++] = "Frame time";
        paramname[t++] = "Binning X";
        paramname[t++] = "Binning Y";
        paramname[t++] = "CF X distance";
        paramname[t++] = "CF Y distance";
        paramname[t++] = "Correlator P";
        paramname[t++] = "Correlator Q";
        paramname[t++] = "Fit model";
        paramname[t++] = "FCCS Display";
        paramname[t++] = "PixelSize";
        paramname[t++] = "Overlap";
        paramname[t++] = "Magnification";
        paramname[t++] = "NA";
        paramname[t++] = "Em Lambda";
        paramname[t++] = "Em Lambda2";
        paramname[t++] = "Sigma";
        paramname[t++] = "SigmaZ";
        paramname[t++] = "Sigma2";
        paramname[t++] = "SigmaZ2";
        paramname[t++] = "Background";
        paramname[t++] = "Background2";
        paramname[t++] = "Background File";
        paramname[t++] = "Bleach Correction";
        paramname[t++] = "Sliding Window Length";
        paramname[t++] = "Polynomial Order";
        paramname[t++] = "Filter";
        paramname[t++] = "filterUL";
        paramname[t++] = "filterLL";
        paramname[t++] = "Threshold";
        paramname[t++] = "Fit";
        paramname[t++] = "DL All bins";
        paramname[t++] = "DL Map bins";
        paramname[t++] = "DL Map width";
        paramname[t++] = "DL Map height";
        paramname[t++] = "PSF bin";
        paramname[t++] = "N of PSF";
        paramname[t++] = "N&B calculated";
        paramname[t++] = "N&B corrected by slope";
        paramname[t++] = "N&B slope";

        paramname[t++] = "Simulation";
        paramname[t++] = "Mode";
        paramname[t++] = "Triplet";
        paramname[t++] = "Seed";
        paramname[t++] = "Particle #";
        paramname[t++] = "CPS";
        paramname[t++] = "Bleach Time";
        paramname[t++] = "Pixel #";
        paramname[t++] = "Extension";
        paramname[t++] = "Frame #";
        paramname[t++] = "Time res";
        paramname[t++] = "Steps per frame";
        paramname[t++] = "Step Size [nm] for D1";
        paramname[t++] = "D1";
        paramname[t++] = "Dout/Din";
        paramname[t++] = "D2";
        paramname[t++] = "F2";
        paramname[t++] = "D3";
        paramname[t++] = "F3";
        paramname[t++] = "triplet kon";
        paramname[t++] = "triplet koff";
        paramname[t++] = "Camera offset";
        paramname[t++] = "Camera noise";
        paramname[t++] = "FRAP Radius [um]";
        paramname[t++] = "FRAP Frame";
        paramname[t++] = "Domain Radius [nm]";
        paramname[t++] = "Domain Density [1/um^2]";
        paramname[t++] = "Pin";
        paramname[t++] = "Pout";
        paramname[t++] = "Mesh Size";
        paramname[t++] = "Hop Probability";
        paramname[t++] = "Pixel Size";
        paramname[t++] = "Objestive Magnification";
        paramname[t++] = "Emission Wavelength";
        paramname[t++] = "NA";
        paramname[t++] = "sigma0";
        paramname[t++] = "sigmaZ";

        t = 0;
        paramsave[t++] = Integer.toString(width);
        paramsave[t++] = Integer.toString(height);
        paramsave[t++] = panelSettings[0];
        paramsave[t++] = panelSettings[1];
        paramsave[t++] = panelSettings[2];
        paramsave[t++] = panelSettings[3];
        paramsave[t++] = panelSettings[4];
        paramsave[t++] = panelSettings[5];
        paramsave[t++] = panelSettings[6];
        paramsave[t++] = panelSettings[7];
        paramsave[t++] = panelSettings[8];
        paramsave[t++] = panelSettings[9];
        paramsave[t++] = panelSettings[10];
        paramsave[t++] = panelSettings[11];
        paramsave[t++] = panelSettings[12];
        paramsave[t++] = panelSettings[13];
        paramsave[t++] = panelSettings[14];
        paramsave[t++] = panelSettings[15];
        paramsave[t++] = panelSettings[16];
        paramsave[t++] = panelSettings[17];
        paramsave[t++] = panelSettings[18];
        paramsave[t++] = panelSettings[19];
        paramsave[t++] = panelSettings[20];
        paramsave[t++] = panelSettings[21];
        paramsave[t++] = panelSettings[22];
        paramsave[t++] = panelSettings[23];
        paramsave[t++] = panelSettings[24];
        paramsave[t++] = panelSettings[25];
        paramsave[t++] = panelSettings[26];
        paramsave[t++] = panelSettings[27];
        paramsave[t++] = panelSettings[28];
        paramsave[t++] = panelSettings[29];
        paramsave[t++] = Boolean.toString(userThreshold[2]);
        paramsave[t++] = panelSettings[30];
        paramsave[t++] = Integer.toString(difflawallbin);
        paramsave[t++] = Integer.toString(difflawmapbin);
        paramsave[t++] = Integer.toString(diffLawMapwidth);
        paramsave[t++] = Integer.toString(diffLawMapheight);
        paramsave[t++] = Integer.toString(psfmaxbin);
        paramsave[t++] = Integer.toString(numofpsf);
        paramsave[t++] = Boolean.toString(NBperformed);
        paramsave[t++] = Boolean.toString(NBcorrected);
        paramsave[t++] = Double.toString(NBslope);
        paramsave[t++] = "";
        paramsave[t++] = simSettings[0];
        paramsave[t++] = simSettings[1];
        paramsave[t++] = simSettings[2];
        paramsave[t++] = simSettings[3];
        paramsave[t++] = simSettings[4];
        paramsave[t++] = simSettings[5];
        paramsave[t++] = simSettings[6];
        paramsave[t++] = simSettings[7];
        paramsave[t++] = simSettings[8];
        paramsave[t++] = simSettings[9];
        paramsave[t++] = simSettings[10];
        paramsave[t++] = simSettings[11];
        paramsave[t++] = simSettings[12];
        paramsave[t++] = simSettings[13];
        paramsave[t++] = simSettings[14];
        paramsave[t++] = simSettings[15];
        paramsave[t++] = simSettings[16];
        paramsave[t++] = simSettings[17];
        paramsave[t++] = simSettings[18];
        paramsave[t++] = simSettings[19];
        paramsave[t++] = simSettings[20];
        paramsave[t++] = simSettings[21];
        paramsave[t++] = simSettings[22];
        paramsave[t++] = simSettings[23];
        paramsave[t++] = simSettings[24];
        paramsave[t++] = simSettings[25];
        paramsave[t++] = simSettings[26];
        paramsave[t++] = simSettings[27];
        paramsave[t++] = simSettings[28];
        paramsave[t++] = simSettings[29];
        paramsave[t++] = simSettings[30];
        paramsave[t++] = simSettings[31];
        paramsave[t++] = simSettings[32];
        paramsave[t++] = simSettings[33];
        paramsave[t++] = simSettings[34];
        paramsave[t++] = simSettings[35];

        for (int i = 0; i < nofpv; i++) {
            row = paraSheet.createRow(i);
            row.createCell(0).setCellValue(paramname[i]);
            row.createCell(1).setCellValue(paramsave[i]);
        }

        row = tauSheet.createRow(0);
        row.createCell(0).setCellValue("S/N");
        row.createCell(1).setCellValue("Lagtime");
        row.createCell(2).setCellValue("Bin width");
        for (int i = 0; i < chanum; i++) {
            row = tauSheet.createRow(i + 1);
            row.createCell(0).setCellValue(i);
            row.createCell(1).setCellValue(lagtime[i]);
            row.createCell(2).setCellValue(samp[i]);
        }

        paraSheet.autoSizeColumn(0);
        paraSheet.autoSizeColumn(1);

        // write Threshold settings
        thresholdSheet.createRow(0);
        thresholdSheet.getRow(0).createCell(0).setCellValue("Parameter");
        thresholdSheet.getRow(0).createCell(1).setCellValue("Filter");
        thresholdSheet.getRow(0).createCell(2).setCellValue("Min 1");
        thresholdSheet.getRow(0).createCell(3).setCellValue("Max 1");
        if (userThreshold[1]) {
            thresholdSheet.getRow(0).createCell(4).setCellValue("Min CCF");
            thresholdSheet.getRow(0).createCell(5).setCellValue("Max CCF");
        }

        for (int i = 1; i < noparam + 3; i++) {
            thresholdSheet.createRow(i);
            thresholdSheet.getRow(i).createCell(0).setCellValue($param[i - 1]);
            thresholdSheet.getRow(i).createCell(1).setCellValue(Boolean.toString(paramfilter[i - 1]));
            thresholdSheet.getRow(i).createCell(2).setCellValue(filterThresholds[0][i - 1][0]);
            thresholdSheet.getRow(i).createCell(3).setCellValue(filterThresholds[0][i - 1][1]);
            if (userThreshold[1]) {
                thresholdSheet.getRow(i).createCell(4).setCellValue(filterThresholds[2][i - 1][0]);
                thresholdSheet.getRow(i).createCell(5).setCellValue(filterThresholds[2][i - 1][1]);
            }
        }

        thresholdSheet.getRow(noparam + 2).createCell(0).setCellValue("G relative");
        thresholdSheet.createRow(noparam + 3);
        thresholdSheet.getRow(noparam + 3).createCell(0).setCellValue("User set thresholds");
        thresholdSheet.getRow(noparam + 3).createCell(1).setCellValue(Boolean.toString(userThreshold[0]));
        thresholdSheet.createRow(noparam + 4);
        thresholdSheet.getRow(noparam + 4).createCell(0).setCellValue("Thresholds for DC-FCCS");
        thresholdSheet.getRow(noparam + 4).createCell(1).setCellValue(Boolean.toString(userThreshold[1]));
        thresholdSheet.createRow(noparam + 5);
        thresholdSheet.getRow(noparam + 5).createCell(0).setCellValue("CCF thresholds same as ACF");
        thresholdSheet.getRow(noparam + 5).createCell(1).setCellValue(Boolean.toString(userThreshold[3]));
        thresholdSheet.createRow(noparam + 6);
        thresholdSheet.getRow(noparam + 6).createCell(0).setCellValue("q=0 when CCF invalid and ACF valid");
        thresholdSheet.getRow(noparam + 6).createCell(1).setCellValue(Boolean.toString(userThreshold[4]));

        thresholdSheet.autoSizeColumn(0);

        // write ACF, SD, fitacf, res
        Row rowacf0 = acf0Sheet.createRow(0);
        Row rowacf1 = acf1Sheet.createRow(0);
        Row rowacf2 = acf2Sheet.createRow(0);
        Row rowsd0 = sd0Sheet.createRow(0);
        Row rowsd1 = sd1Sheet.createRow(0);
        Row rowsd2 = sd2Sheet.createRow(0);
        Row rowfitacf0 = fitacf0Sheet.createRow(0);
        Row rowfitacf1 = fitacf1Sheet.createRow(0);
        Row rowfitacf2 = fitacf2Sheet.createRow(0);
        Row rowres0 = res0Sheet.createRow(0);
        Row rowres1 = res1Sheet.createRow(0);
        Row rowres2 = res2Sheet.createRow(0);
        Row rowmsd0 = msd0Sheet.createRow(0);
        Row rowmsd1 = msd1Sheet.createRow(0);
        Row rowmsd2 = msd2Sheet.createRow(0);

        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                String $tmp = "(" + Integer.toString(w) + ", " + Integer.toString(h) + ")";
                rowacf0.createCell(w + h * width).setCellValue($tmp);
                rowacf1.createCell(w + h * width).setCellValue($tmp);
                rowacf2.createCell(w + h * width).setCellValue($tmp);
                rowsd0.createCell(w + h * width).setCellValue($tmp);
                rowsd1.createCell(w + h * width).setCellValue($tmp);
                rowsd2.createCell(w + h * width).setCellValue($tmp);
                rowfitacf0.createCell(w + h * width).setCellValue($tmp);
                rowfitacf1.createCell(w + h * width).setCellValue($tmp);
                rowfitacf2.createCell(w + h * width).setCellValue($tmp);
                rowres0.createCell(w + h * width).setCellValue($tmp);
                rowres1.createCell(w + h * width).setCellValue($tmp);
                rowres2.createCell(w + h * width).setCellValue($tmp);
                rowmsd0.createCell(w + h * width).setCellValue($tmp);
                rowmsd1.createCell(w + h * width).setCellValue($tmp);
                rowmsd2.createCell(w + h * width).setCellValue($tmp);
            }
        }

        for (int i = 0; i < chanum; i++) {
            acf0Sheet.createRow(i + 1);
            acf1Sheet.createRow(i + 1);
            acf2Sheet.createRow(i + 1);
            sd0Sheet.createRow(i + 1);
            sd1Sheet.createRow(i + 1);
            sd2Sheet.createRow(i + 1);
            fitacf0Sheet.createRow(i + 1);
            fitacf1Sheet.createRow(i + 1);
            fitacf2Sheet.createRow(i + 1);
            res0Sheet.createRow(i + 1);
            res1Sheet.createRow(i + 1);
            res2Sheet.createRow(i + 1);
            msd0Sheet.createRow(i + 1);
            msd1Sheet.createRow(i + 1);
            msd2Sheet.createRow(i + 1);
        }

        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                for (int i = 0; i < chanum; i++) {
                    acf0Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(acf[0][w][h][i]);
                    acf1Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(acf[1][w][h][i]);
                    acf2Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(acf[2][w][h][i]);
                    sd0Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(sdacf[0][w][h][i]);
                    sd1Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(sdacf[1][w][h][i]);
                    sd2Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(sdacf[2][w][h][i]);
                    fitacf0Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(fitacf[0][w][h][i]);
                    fitacf1Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(fitacf[1][w][h][i]);
                    fitacf2Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(fitacf[2][w][h][i]);
                    res0Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(res[0][w][h][i]);
                    res1Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(res[1][w][h][i]);
                    res2Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(res[2][w][h][i]);
                    msd0Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(msd[0][w][h][i]);
                    msd1Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(msd[1][w][h][i]);
                    msd2Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(msd[2][w][h][i]);
                }
            }
        }

        // add fit results
        row0 = fitres0Sheet.createRow(0);
        row1 = fitres1Sheet.createRow(0);
        row2 = fitres2Sheet.createRow(0);
        row0.createCell(0).setCellValue("Parameter");
        row1.createCell(0).setCellValue("Parameter");
        row2.createCell(0).setCellValue("Parameter");
        fitres0Sheet.createRow(1).createCell(0).setCellValue("fitted");
        fitres1Sheet.createRow(1).createCell(0).setCellValue("fitted");
        fitres2Sheet.createRow(1).createCell(0).setCellValue("fitted");
        for (int q = 0; q < noparam + 3; q++) {
            fitres0Sheet.createRow(q + 2).createCell(0).setCellValue($param[q]);
            fitres1Sheet.createRow(q + 2).createCell(0).setCellValue($param[q]);
            fitres2Sheet.createRow(q + 2).createCell(0).setCellValue($param[q]);
        }
        fitres2Sheet.createRow(noparam + 5).createCell(0).setCellValue($param[noparam + 3]);

        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                fitres0Sheet.getRow(0).createCell(w + h * width + 1).setCellValue("(" + Integer.toString(w) + ", " + Integer.toString(h) + ")");
                fitres1Sheet.getRow(0).createCell(w + h * width + 1).setCellValue("(" + Integer.toString(w) + ", " + Integer.toString(h) + ")");
                fitres2Sheet.getRow(0).createCell(w + h * width + 1).setCellValue("(" + Integer.toString(w) + ", " + Integer.toString(h) + ")");
                fitres0Sheet.getRow(1).createCell(w + h * width + 1).setCellValue(Boolean.toString(pixfitted[0][w][h]));
                fitres1Sheet.getRow(1).createCell(w + h * width + 1).setCellValue(Boolean.toString(pixfitted[1][w][h]));
                fitres2Sheet.getRow(1).createCell(w + h * width + 1).setCellValue(Boolean.toString(pixfitted[2][w][h]));
                fitres0Sheet.getRow(noparam + 2).createCell(w + h * width + 1).setCellValue(chi2[0][w][h]);
                fitres1Sheet.getRow(noparam + 2).createCell(w + h * width + 1).setCellValue(chi2[1][w][h]);
                fitres2Sheet.getRow(noparam + 2).createCell(w + h * width + 1).setCellValue(chi2[2][w][h]);
                fitres0Sheet.getRow(noparam + 3).createCell(w + h * width + 1).setCellValue(blocked[0][w][h]);
                fitres1Sheet.getRow(noparam + 3).createCell(w + h * width + 1).setCellValue(blocked[1][w][h]);
                fitres2Sheet.getRow(noparam + 3).createCell(w + h * width + 1).setCellValue(blocked[2][w][h]);
                fitres0Sheet.getRow(noparam + 4).createCell(w + h * width + 1).setCellValue(pixvalid[0][w][h]);
                fitres1Sheet.getRow(noparam + 4).createCell(w + h * width + 1).setCellValue(pixvalid[1][w][h]);
                fitres2Sheet.getRow(noparam + 4).createCell(w + h * width + 1).setCellValue(pixvalid[2][w][h]);
                fitres2Sheet.getRow(noparam + 5).createCell(w + h * width + 1).setCellValue(CCFq[w][h]);
            }
        }

        for (int i = 0; i < noparam; i++) {
            for (int h = 0; h < height; h++) {
                for (int w = 0; w < width; w++) {
                    fitres0Sheet.getRow(i + 2).createCell(w + h * width + 1).setCellValue(fitres[0][w][h][i]);
                    fitres1Sheet.getRow(i + 2).createCell(w + h * width + 1).setCellValue(fitres[1][w][h][i]);
                    fitres2Sheet.getRow(i + 2).createCell(w + h * width + 1).setCellValue(fitres[2][w][h][i]);
                }
            }
        }

        fitres0Sheet.autoSizeColumn(0);
        fitres1Sheet.autoSizeColumn(0);
        fitres2Sheet.autoSizeColumn(0);

        // add diffLaw data
        row = diffLawSheet.createRow(0);
        row.createCell(0).setCellValue("Aeff");
        row.createCell(1).setCellValue("Time");
        row.createCell(2).setCellValue("SD");
        row.createCell(3).setCellValue("intercept");
        row.createCell(4).setCellValue("slope");
        row.createCell(5).setCellValue("fit start");
        row.createCell(6).setCellValue("fit end");

        for (int i = 0; i < difflawallbin; i++) {
            row = diffLawSheet.createRow(i + 1);
            row.createCell(0).setCellValue(difflaw[0][i]);
            row.createCell(1).setCellValue(difflaw[1][i]);
            row.createCell(2).setCellValue(difflaw[2][i]);
        }
        row = diffLawSheet.getRow(1);
        row.createCell(3).setCellValue(difflawfit[0]);
        row.createCell(4).setCellValue(difflawfit[1]);
        row.createCell(5).setCellValue(diffLawFitLim[0]);
        row.createCell(6).setCellValue(diffLawFitLim[1]);

        // add diffLaw Map data 
        row = diffLawMapDataSheet.createRow(0);
        row.createCell(0).setCellValue("Obs Vol");

        for (int w = 0; w < diffLawMapwidth; w++) {
            for (int h = 0; h < diffLawMapheight; h++) {
                row.createCell(w + h * diffLawMapwidth + 1).setCellValue("(" + Integer.toString(w) + ", " + Integer.toString(h) + ")");
            }
        }
        for (int k = 0; k < difflawmapbin; k++) {
            row = diffLawMapDataSheet.createRow(k + 1);
            row.createCell(0).setCellValue(difflawarray[0][0][0][k]);
            for (int w = 0; w < diffLawMapwidth; w++) {
                for (int h = 0; h < diffLawMapheight; h++) {
                    row.createCell(w + h * diffLawMapwidth + 1).setCellValue(difflawarray[w][h][1][k]);
                }
            }
        }

        // add diffLaw Map SD
        row = diffLawMapSDSheet.createRow(0);
        row.createCell(0).setCellValue("");
        for (int w = 0; w < diffLawMapwidth; w++) {
            for (int h = 0; h < diffLawMapheight; h++) {
                row.createCell(w + h * diffLawMapwidth).setCellValue("(" + Integer.toString(w) + ", " + Integer.toString(h) + ")");
            }
        }

        for (int k = 0; k < difflawmapbin; k++) {
            row = diffLawMapSDSheet.createRow(k + 1);
            for (int w = 0; w < diffLawMapwidth; w++) {
                for (int h = 0; h < diffLawMapheight; h++) {
                    row.createCell(w + h * diffLawMapwidth).setCellValue(difflawarray[w][h][2][k]);
                }
            }
        }

        // add diffLaw Map Fit data 
        row = diffLawMapFitSheet.createRow(0);
        row.createCell(0).setCellValue("Parameter");
        for (int w = 0; w < diffLawMapwidth; w++) {
            for (int h = 0; h < diffLawMapheight; h++) {
                row.createCell(w + h * diffLawMapwidth + 1).setCellValue("(" + Integer.toString(w) + ", " + Integer.toString(h) + ")");
            }
        }
        row1 = diffLawMapFitSheet.createRow(1);
        row2 = diffLawMapFitSheet.createRow(2);
        row1.createCell(0).setCellValue("intercept");
        row2.createCell(0).setCellValue("slope");
        for (int w = 0; w < diffLawMapwidth; w++) {
            for (int h = 0; h < diffLawMapheight; h++) {
                row1.createCell(w + h * diffLawMapwidth + 1).setCellValue(diffLawFitMap[w][h][0]);
                row2.createCell(w + h * diffLawMapwidth + 1).setCellValue(diffLawFitMap[w][h][1]);
            }
        }

        // add PSF data
        row = psfDataSheet.createRow(0);
        for (int i = 0; i < numofpsf; i++) {
            row.createCell(i * 3).setCellValue("Bin");
            row.createCell(1 + i * 3).setCellValue("D (bin " + Integer.toString(i + 1) + ")");
            row.createCell(2 + i * 3).setCellValue("SD (bin " + Integer.toString(i + 1) + ")");
        }

        for (int k = 0; k < psfmaxbin; k++) {
            row = psfDataSheet.createRow(k + 1);
            for (int i = 0; i < numofpsf; i++) {
                row.createCell(i * 3).setCellValue(psfData[i][0][k]);
                row.createCell(i * 3 + 1).setCellValue(psfData[i][1][k]);
                row.createCell(i * 3 + 2).setCellValue(psfData[i][2][k]);
            }
        }

        row = psfDataSheet.createRow(psfmaxbin + 2);
        row.createCell(0).setCellValue("PSF start");
        row.createCell(1).setCellValue("PSF end");
        row.createCell(2).setCellValue("PSF step");
        row = psfDataSheet.createRow(psfmaxbin + 3);
        row.createCell(0).setCellValue(psfStart);
        row.createCell(1).setCellValue(psfEnd);
        row.createCell(2).setCellValue(psfStep);

        // add dCCF data
        row = dCCFDataSheet.createRow(0);
        row.createCell(0).setCellValue("Direction");
        row.createCell(1).setCellValue("Calculated");
        row = dCCFDataSheet.createRow(1);
        row.createCell(0).setCellValue("horizontal");
        row = dCCFDataSheet.createRow(2);
        row.createCell(0).setCellValue("vertical");
        row = dCCFDataSheet.createRow(3);
        row.createCell(0).setCellValue("diagonal /");
        row = dCCFDataSheet.createRow(4);
        row.createCell(0).setCellValue("diagonal \\");

        for (int i = 0; i < dccfMax; i++) {
            dCCFDataSheet.getRow(i + 1).createCell(1).setCellValue(Boolean.toString(dccfCalculated[i]));
        }

        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                dCCFDataSheet.getRow(0).createCell(w + h * width + 2).setCellValue("(" + Integer.toString(w) + ", " + Integer.toString(h) + ")");
                for (int i = 0; i < dccfMax; i++) {
                    dCCFDataSheet.getRow(i + 1).createCell(w + h * width + 2).setCellValue(dccf[i][w][h]);
                }
            }
        }

        // add N&B data
        row = NBDataSheet.createRow(0);
        row.createCell(0).setCellValue(" ");
        row = NBDataSheet.createRow(1);
        row.createCell(0).setCellValue("Number");
        row = NBDataSheet.createRow(2);
        row.createCell(0).setCellValue("Brightness");
        row = NBDataSheet.createRow(3);
        row.createCell(0).setCellValue("Num (corrected)");
        row = NBDataSheet.createRow(4);
        row.createCell(0).setCellValue("Epsilon (corrected)");

        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                NBDataSheet.getRow(0).createCell(w + h * width + 2).setCellValue("(" + Integer.toString(w) + ", " + Integer.toString(h) + ")");
                NBDataSheet.getRow(1).createCell(w + h * width + 2).setCellValue(NBN[w][h]);
                NBDataSheet.getRow(2).createCell(w + h * width + 2).setCellValue(NBB[w][h]);
                NBDataSheet.getRow(3).createCell(w + h * width + 2).setCellValue(NBNum[w][h]);
                NBDataSheet.getRow(4).createCell(w + h * width + 2).setCellValue(NBEps[w][h]);
            }
        }

        // add Reconstruction Parameters, used to reconstruct the experiments
        t = 0;
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("ACF Window");
        if ((acfWindow != null && acfWindow.isClosed() == false) || batchCorrelateAll) {
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Intensity Trace Window");
        if ((intWindow != null && intWindow.isClosed() == false) || batchCorrelateAll) {	// intensity trace window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Residual Window");
        if ((resWindow != null && resWindow.isClosed() == false) || batchCorrelateAll) {	// fit residuals window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("SD Window");
        if ((sdWindow != null && sdWindow.isClosed() == false) || batchCorrelateAll) {	// SD window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("MSD Window");
        if (msdWindow != null && msdWindow.isClosed() == false) {	// MSD window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Difusion Law Window");
        if ((difflawWindow != null && difflawWindow.isClosed() == false) || batchDiffLaw) {	// diffusion law window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("PSF Window");
        if ((PSFWindow != null && PSFWindow.isClosed() == false) || batchPSF) {	// PSF window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Para Cor Window");
        if (paraCorWindow != null && paraCorWindow.isClosed() == false) {	// paraCor window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Blocking Window");
        if (blockingWindow != null && blockingWindow.isClosed() == false) {	// blocking window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Covariance Window");
        if (impCovWin != null && impCovWin.isClosed() == false) {	// covariance window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("N Window (N&B)");
        if (impNWin != null && impNWin.isClosed() == false) {	// N&B N window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("B Window (N&B)");
        if (impBWin != null && impBWin.isClosed() == false) {	// N&B B window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Num Window (N&B)");
        if (impNumWin != null && impNumWin.isClosed() == false) {	// N&B n window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Epsilon Window (N&B)");
        if (impEpsWin != null && impEpsWin.isClosed() == false) {	// N&B epsilon window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Parameter Window");
        if ((impPara1 != null && impPara1.isVisible()) || batchFit) {							// parameter map window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("Histogram Window");
        if ((histWin != null && histWin.isClosed() == false) || batchFit) {					// histogram windows
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("DL Map Window");
        if (impDLMap != null && impDLMap.isVisible()) {					// DL map window
            row.createCell(1).setCellValue(1);
        } else {
            row.createCell(1).setCellValue(0);
        }

        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("cormode");						// cormode (cormode 1: ACF, cormode 2: ACF or CCF, cormode 3: DC-FCCS, cormode 4:)
        row.createCell(1).setCellValue(parcormode);
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("plotACFCurves");				// parameter whether to plot ACF curves
        row.createCell(1).setCellValue((int) (plotACFCurves ? 1 : 0));
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("plotSDCurves");				// parameter whether to plot ACF curves
        row.createCell(1).setCellValue((int) (plotSDCurves ? 1 : 0));
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("plotIntensityCurves");				// parameter whether to plot ACF curves
        row.createCell(1).setCellValue((int) (plotIntensityCurves ? 1 : 0));
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("plotResCurves");				// parameter whether to plot ACF curves
        row.createCell(1).setCellValue((int) (plotResCurves ? 1 : 0));
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("plotParaHist");				// parameter whether to plot ACF curves
        row.createCell(1).setCellValue((int) (plotParaHist ? 1 : 0));
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("plotBlockingCurve");				// parameter whether to plot ACF curves
        row.createCell(1).setCellValue((int) (plotBlockingCurve ? 1 : 0));
        row = ReconParaSheet.createRow(t++);
        row.createCell(0).setCellValue("plotCovmats");				// parameter whether to plot ACF curves
        row.createCell(1).setCellValue((int) (plotCovmats ? 1 : 0));

        int dccfcount = 0;
        if (batchVerDCCF) {
            dccfcount++;
        }
        if (batchHorDCCF) {
            dccfcount++;
        }
        if (batchDiaUpDCCF) {
            dccfcount++;
        }
        if (batchDiaDownDCCF) {
            dccfcount++;
        }

        for (int x = 0; x < dccfMax; x++) {
            row = ReconParaSheet.createRow(t++);
            row.createCell(0).setCellValue("DCCF Window " + Integer.toString(x));
            if ((impDCCFWin[x] != null && impDCCFWin[x].isClosed() == false) || x < dccfcount) {	// dCCF windows
                row.createCell(1).setCellValue(1);
                row.createCell(2).setCellValue(impDCCFWin[x].getHeight());
                row.createCell(3).setCellValue(impDCCFWin[x].getWidth());
                row.createCell(4).setCellValue($dccfTitle);
            } else {
                row.createCell(1).setCellValue(0);
            }
        }
        for (int x = 0; x < dccfMax; x++) {
            row = ReconParaSheet.createRow(t++);
            row.createCell(0).setCellValue("DCCF Histogram Window " + Integer.toString(x));
            if (histDCCFWin[x] != null && histDCCFWin[x].isClosed() == false) {	// dCCF historgram windows
                row.createCell(1).setCellValue(1);
                row.createCell(2).setCellValue($histDCCFWinTitle[x]);
            } else {
                row.createCell(1).setCellValue(0);
            }
        }

        ReconParaSheet.autoSizeColumn(0);

        // Write the output to a file
        try {
            FileOutputStream fileOut = new FileOutputStream(file);
            wb.write(fileOut);
            fileOut.close();
        } catch (IOException e) {
            throw new RuntimeException($exception, e);
        }
        if (showlog) {
            IJ.log("File(s) saved");
        }
    }

    public void readExperiment(File file, String $exception) {
        XSSFWorkbook wb;
        int PARAM = 0;
        int THRESHOLD = 1;
        int LAGTIME = 2;
        int ACF0 = 3;
        int SD0 = 4;
        int FITACF0 = 5;
        int RES0 = 6;
        int FITRES0 = 7;
        int ACF1 = 8;
        int SD1 = 9;
        int FITACF1 = 10;
        int RES1 = 11;
        int FITRES1 = 12;
        int ACF2 = 13;
        int SD2 = 14;
        int FITACF2 = 15;
        int RES2 = 16;
        int FITRES2 = 17;
        int DIFFLAW = 18;
        int DIFFLAWMAPDATA = 19;
        int DIFFLAWMAPSD = 20;
        int DIFFLAWMAPFIT = 21;
        int PSFDATA = 22;
        int MSD0 = 23;
        int MSD1 = 24;
        int MSD2 = 25;
        int DCCF = 26;
        int NB = 27;
        int RECONPARA = 28;
        int t;
        int rowcount;
        int cellcount;
        int h;
        int w;
        double minval;
        double maxval;
        double tmpread;
        boolean userThresholdM;
        String $direction = "horizontal";
        Iterator<Row> rowIterator;
        Iterator<Cell> cellIterator;
        ArrayList<Integer> reconpara = new ArrayList<>();

        try {
            FileInputStream inp = new FileInputStream(file);
            wb = new XSSFWorkbook(inp);
            //Workbook wb = WorkbookFactory.create(inp);
        } catch (IOException e) {
            throw new RuntimeException($exception, e);
        }

        // the follwoing two stat variables avoid user prompting because parameters or panel items were changed
        // both parameters will be set back after setParmaetrs is invoked after the next read block
        askOnRewrite = false;	// don't prompt the user
        expload = true;			// indicate that this is a load procedure limiting initializations in setParameters() 

        // read panel parameters
        Sheet sheet = wb.getSheetAt(PARAM);
        t = 0;
        // a simple sanity check. The dimensions of input data should match current image file
        int data_width = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        int data_height = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        int data_firstframe = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        int data_lastframe = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        boolean DataOK = (data_width == width) && (data_height == height) && ((data_lastframe - data_firstframe + 1) == frames);
        if (!DataOK) {
            IJ.log("Experimental data does not match the dimensions of current image stack.");
            return;
        }

        t = 0;
        width = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        height = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfFirstFrame.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfLastFrame.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfFrameTime.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        String bx = sheet.getRow(t++).getCell(1).getStringCellValue();
        String by = sheet.getRow(t++).getCell(1).getStringCellValue();
        tfBinning.setText(bx + " x " + by);
        tfCfXDistance.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfCfYDistance.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        cbCorrelatorP.setSelectedItem(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfCorrelatorQ.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        cbFitModel.setSelectedItem(sheet.getRow(t++).getCell(1).getStringCellValue());
        tbFCCSDisplay.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        if (tbFCCSDisplay.getText().equals("On")) {
            tbFCCSDisplay.setSelected(true);
        } else {
            tbFCCSDisplay.setSelected(false);
        }
        tfPixelSize.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tbOverlap.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        if (tbOverlap.getText().equals("On")) {
            tbOverlap.setSelected(true);
        } else {
            tbOverlap.setSelected(false);
        }
        tfMagnification.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfNA.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfEmLambda.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfEmLambda2.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfSigma.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfSigmaZ.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfSigma2.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfSigmaZ2.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfBackground.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        tfBackground2.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        $bgrTitleMem = sheet.getRow(t++).getCell(1).getStringCellValue();
        cbBleachCor.setSelectedItem(sheet.getRow(t++).getCell(1).getStringCellValue());
        slidingWindowLength = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        polyOrder = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        cbFilter.setSelectedItem(sheet.getRow(t++).getCell(1).getStringCellValue());
        filterUL = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        filterLL = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        userThresholdM = Boolean.parseBoolean(sheet.getRow(t++).getCell(1).getStringCellValue());
        tbFit.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        if (tbFit.getText().equals("Fit On")) {
            doFit = true;
            tbFit.setSelected(true);
            fitframe.setVisible(true);
        } else {
            tbFit.setSelected(false);
            fitframe.setVisible(false);
        }
        try {
            difflawallbin = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            cbSimMode.setSelectedItem("");
        }
        difflawmapbin = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        diffLawMapwidth = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        diffLawMapheight = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        psfmaxbin = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        numofpsf = Integer.parseInt(sheet.getRow(t++).getCell(1).getStringCellValue());
        NBperformed = Boolean.parseBoolean(sheet.getRow(t++).getCell(1).getStringCellValue());
        NBcorrected = Boolean.parseBoolean(sheet.getRow(t++).getCell(1).getStringCellValue());
        NBslope = Double.parseDouble(sheet.getRow(t++).getCell(1).getStringCellValue());

        t++;
        try {
            cbSimMode.setSelectedItem(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            cbSimMode.setSelectedItem("");
        }
        t++;
        try {
            tfSimSeed.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimSeed.setText("");
        }
        try {
            tfSimParticleNum.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimParticleNum.setText("");
        }
        try {
            tfSimCPS.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimCPS.setText("");
        }
        try {
            tfSimTauBleach.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimTauBleach.setText("");
        }
        try {
            tfSimPixelNum.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimPixelNum.setText("");
        }
        try {
            tfSimExtensionFactor.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimExtensionFactor.setText("");
        }
        try {
            tfSimTimeStepNum.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimTimeStepNum.setText("");
        }
        try {
            tfSimFrameTime.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimFrameTime.setText("");
        }
        try {
            tfSimStepsPerFrame.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimStepsPerFrame.setText("");
        }
        try {
            tfSimCurrentStepSize.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimCurrentStepSize.setText("");
        }
        try {
            tfSimD1.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimD1.setText("");
        }
        try {
            tfSimDoutDinRatio.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimDoutDinRatio.setText("");
        }
        try {
            tfSimD2.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimD2.setText("");
        }
        try {
            tfSimF2.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimF2.setText("");
        }
        try {
            tfSimD3.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimD3.setText("");
        }
        try {
            tfSimF3.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimF3.setText("");
        }
        try {
            tfSimKon.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimKon.setText("");
        }
        try {
            tfSimKoff.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimKoff.setText("");
        }
        try {
            tfSimCameraOffset.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimCameraOffset.setText("");
        }
        try {
            tfSimCameraNoiseFactor.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimCameraNoiseFactor.setText("");
        }
        try {
            tfSimBleachRadius.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimBleachRadius.setText("");
        }
        try {
            tfSimBleachFrame.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfSimBleachFrame.setText("");
        }
        try {
            tfDomainRadius.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfDomainRadius.setText("");
        }
        try {
            tfDomainDensity.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfDomainDensity.setText("");
        }
        try {
            tfPin.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfPin.setText("");
        }
        try {
            tfPout.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfPout.setText("");
        }
        try {
            tfMeshworkSize.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfMeshworkSize.setText("");
        }
        try {
            tfHopProbability.setText(sheet.getRow(t++).getCell(1).getStringCellValue());
        } catch (NullPointerException e) {
            tfHopProbability.setText("");
        }

        setParameters();
        askOnRewrite = true;
        expload = false;

        // read Threshold settings
        try {
            sheet = wb.getSheetAt(THRESHOLD);

            userThreshold[0] = Boolean.parseBoolean(sheet.getRow(noparam + 3).getCell(1).getStringCellValue());
            userThreshold[1] = Boolean.parseBoolean(sheet.getRow(noparam + 4).getCell(1).getStringCellValue());
            userThreshold[2] = userThresholdM;
            userThreshold[3] = Boolean.parseBoolean(sheet.getRow(noparam + 5).getCell(1).getStringCellValue());
            userThreshold[4] = Boolean.parseBoolean(sheet.getRow(noparam + 6).getCell(1).getStringCellValue());

            for (int i = 1; i < noparam + 3; i++) {
                paramfilter[i - 1] = Boolean.parseBoolean(sheet.getRow(i).getCell(1).getStringCellValue());
                filterThresholds[0][i - 1][0] = sheet.getRow(i).getCell(2).getNumericCellValue();
                filterThresholds[0][i - 1][1] = sheet.getRow(i).getCell(3).getNumericCellValue();
                if (userThreshold[1]) {
                    filterThresholds[2][i - 1][0] = sheet.getRow(i).getCell(4).getNumericCellValue();
                    filterThresholds[2][i - 1][1] = sheet.getRow(i).getCell(5).getNumericCellValue();
                    filterThresholds[1][i - 1][0] = filterThresholds[0][i - 1][0];
                    filterThresholds[1][i - 1][1] = filterThresholds[0][i - 1][1];
                    if (userThreshold[3]) {
                        filterThresholds[2][i - 1][0] = filterThresholds[0][i - 1][0];
                        filterThresholds[2][i - 1][1] = filterThresholds[0][i - 1][1];
                    }
                }
            }

            if (userThreshold[2]) {		// if Threshold has been applied, open Threshold settings panel
                doFiltering = true;
                SetThresholds();
                tbFiltering.setSelected(true);
                filteringframe.setVisible(true);
            }
        } catch (IllegalArgumentException iae) {
        }

        // read ACF0
        try {
            sheet = wb.getSheetAt(ACF0);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        acf[0][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        acf[0][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read SD0
        try {
            sheet = wb.getSheetAt(SD0);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        sdacf[0][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        sdacf[0][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read fitacf0
        try {
            sheet = wb.getSheetAt(FITACF0);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        fitacf[0][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        fitacf[0][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read res0
        try {
            sheet = wb.getSheetAt(RES0);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        res[0][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        res[0][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read FITRES0
        try {
            sheet = wb.getSheetAt(FITRES0);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                if (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                }
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    pixfitted[0][w][h] = Boolean.parseBoolean(cell.getStringCellValue());
                    cellcount++;
                }
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                if (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                }
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        tmpread = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        tmpread = Double.NaN;
                    }
                    if (rowcount < noparam) {
                        fitres[0][w][h][rowcount] = tmpread;
                    }
                    if (rowcount == noparam) {
                        chi2[0][w][h] = tmpread;
                    }
                    if (rowcount == noparam + 1) {
                        blocked[0][w][h] = tmpread;
                    }
                    if (rowcount == noparam + 2) {
                        pixvalid[0][w][h] = tmpread;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read ACF1
        try {
            sheet = wb.getSheetAt(ACF1);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        acf[1][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        acf[1][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read SD1
        try {
            sheet = wb.getSheetAt(SD1);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        sdacf[1][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        sdacf[1][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read fitacf1
        try {
            sheet = wb.getSheetAt(FITACF1);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        fitacf[1][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        fitacf[1][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read res1
        try {
            sheet = wb.getSheetAt(RES1);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        res[1][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        res[1][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read FITRES1
        try {
            sheet = wb.getSheetAt(FITRES1);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                if (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                }
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    pixfitted[1][w][h] = Boolean.parseBoolean(cell.getStringCellValue());
                    cellcount++;
                }
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                if (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                }
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        tmpread = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        tmpread = Double.NaN;
                    }
                    if (rowcount < noparam) {
                        fitres[1][w][h][rowcount] = tmpread;
                    }
                    if (rowcount == noparam) {
                        chi2[1][w][h] = tmpread;
                    }
                    if (rowcount == noparam + 1) {
                        blocked[1][w][h] = tmpread;
                    }
                    if (rowcount == noparam + 2) {
                        pixvalid[1][w][h] = tmpread;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read ACF2
        try {
            sheet = wb.getSheetAt(ACF2);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        acf[2][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        acf[2][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read SD2
        try {
            sheet = wb.getSheetAt(SD2);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        sdacf[2][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        sdacf[2][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read fitacf2
        try {
            sheet = wb.getSheetAt(FITACF2);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        fitacf[2][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        fitacf[2][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read res2
        try {
            sheet = wb.getSheetAt(RES2);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        res[2][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        res[2][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read FITRES2
        try {
            sheet = wb.getSheetAt(FITRES2);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                if (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                }
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    pixfitted[2][w][h] = Boolean.parseBoolean(cell.getStringCellValue());
                    cellcount++;
                }
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                if (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                }
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        tmpread = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        tmpread = Double.NaN;
                    }
                    if (rowcount < noparam) {
                        fitres[2][w][h][rowcount] = tmpread;
                    }
                    if (rowcount == noparam) {
                        chi2[2][w][h] = tmpread;
                    }
                    if (rowcount == noparam + 1) {
                        blocked[2][w][h] = tmpread;
                    }
                    if (rowcount == noparam + 2) {
                        pixvalid[2][w][h] = tmpread;
                    }
                    if (rowcount == noparam + 3) {
                        CCFq[w][h] = tmpread;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (IllegalArgumentException iae) {
        }

        // read DiffLaw
        try {
            difflaw = new double[3][difflawallbin];
            sheet = wb.getSheetAt(DIFFLAW);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                difflaw[0][rowcount] = row.getCell(0).getNumericCellValue();
                difflaw[1][rowcount] = row.getCell(1).getNumericCellValue();
                difflaw[2][rowcount++] = row.getCell(2).getNumericCellValue();
            }
            difflawfit[0] = sheet.getRow(1).getCell(3).getNumericCellValue();
            difflawfit[1] = sheet.getRow(1).getCell(4).getNumericCellValue();
            diffLawFitLim[0] = (int) sheet.getRow(1).getCell(5).getNumericCellValue();
            diffLawFitLim[1] = (int) sheet.getRow(1).getCell(6).getNumericCellValue();
        } catch (Exception e) {
        }

        // read DiffLawMapData
        try {
            difflawarray = new double[diffLawMapwidth][diffLawMapheight][3][difflawmapbin];
            sheet = wb.getSheetAt(DIFFLAWMAPDATA);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                Cell cell = cellIterator.next();
                double tmpobsvol = cell.getNumericCellValue();
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / diffLawMapwidth);
                    w = cellcount - h * diffLawMapwidth;
                    cell = cellIterator.next();
                    difflawarray[w][h][0][rowcount] = tmpobsvol;
                    difflawarray[w][h][1][rowcount] = cell.getNumericCellValue();
                    cellcount++;
                }
                rowcount++;
            }
        } catch (Exception e) {
        }

        // read DiffLawMapSD
        try {
            sheet = wb.getSheetAt(DIFFLAWMAPSD);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / diffLawMapwidth);
                    w = cellcount - h * diffLawMapwidth;
                    Cell cell = cellIterator.next();
                    try {
                        difflawarray[w][h][2][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        difflawarray[w][h][2][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (Exception e) {
        }

        // read DiffLawMapFit
        try {
            diffLawFitMap = new double[diffLawMapwidth][diffLawMapwidth][2];
            sheet = wb.getSheetAt(DIFFLAWMAPFIT);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                if (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                }
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / diffLawMapwidth);
                    w = cellcount - h * diffLawMapwidth;
                    Cell cell = cellIterator.next();
                    if (rowcount == 0) {
                        diffLawFitMap[w][h][0] = cell.getNumericCellValue();
                    }
                    if (rowcount == 1) {
                        diffLawFitMap[w][h][1] = cell.getNumericCellValue();
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (Exception e) {
        }

        // read PSF Data
        try {
            psfData = new double[numofpsf][3][psfmaxbin];
            sheet = wb.getSheetAt(PSFDATA);
            for (int p = 0; p < numofpsf; p++) {
                for (int b = 0; b < psfmaxbin; b++) {
                    psfData[p][0][b] = sheet.getRow(b + 1).getCell(3 * p).getNumericCellValue();
                    psfData[p][1][b] = sheet.getRow(b + 1).getCell(3 * p + 1).getNumericCellValue();
                    psfData[p][2][b] = sheet.getRow(b + 1).getCell(3 * p + 2).getNumericCellValue();
                }
            }
        } catch (Exception e) {
        }

        // read MSD0
        try {
            sheet = wb.getSheetAt(MSD0);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        msd[0][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        msd[0][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (Exception e) {
        }

        // read MSD1
        try {
            sheet = wb.getSheetAt(MSD0);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        msd[1][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        msd[1][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (Exception e) {
        }

        // read MSD2
        try {
            sheet = wb.getSheetAt(MSD0);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        msd[2][w][h][rowcount] = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        msd[2][w][h][rowcount] = Double.NaN;
                    }
                    cellcount++;
                }
                rowcount++;
            }
        } catch (Exception e) {
        }

        //read dCCF
        try {
            sheet = wb.getSheetAt(DCCF);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                if (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                    $direction = cell.getStringCellValue();
                }
                if (cellIterator.hasNext()) {
                    Cell cell1 = cellIterator.next();
                    dccfCalculated[rowcount] = Boolean.parseBoolean(cell1.getStringCellValue());
                    if (dccfCalculated[rowcount]) {
                        cellcount = 0;
                        while (cellIterator.hasNext()) {
                            h = (int) Math.floor(cellcount / width);
                            w = cellcount - h * width;
                            Cell cell = cellIterator.next();
                            try {
                                dccf[rowcount][w][h] = cell.getNumericCellValue();
                            } catch (IllegalStateException ise) {
                                dccf[rowcount][w][h] = Double.NaN;
                            }
                            cellcount++;
                        }
                        createDCCFWindow(width, height, $dccfTitle + " - " + $direction, "Hist - " + $dccfTitle + " - " + $direction, rowcount);
                    }
                }
                rowcount++;
            }
        } catch (Exception e) {
        }

        //read N&B
        try {
            sheet = wb.getSheetAt(NB);
            rowIterator = sheet.iterator();
            if (rowIterator.hasNext()) {
                Row row = rowIterator.next();
            }
            rowcount = 0;
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                cellIterator = row.cellIterator();
                if (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                }
                cellcount = 0;
                while (cellIterator.hasNext()) {
                    h = (int) Math.floor(cellcount / width);
                    w = cellcount - h * width;
                    Cell cell = cellIterator.next();
                    try {
                        tmpread = cell.getNumericCellValue();
                    } catch (IllegalStateException ise) {
                        tmpread = Double.NaN;
                    }
                    if (rowcount == 0) {
                        NBN[w][h] = tmpread;
                    }
                    if (rowcount == 1) {
                        NBB[w][h] = tmpread;
                    }
                    if (rowcount == 2) {
                        NBNum[w][h] = tmpread;
                    }
                    if (rowcount == 3) {
                        NBEps[w][h] = tmpread;
                    }
                    cellcount++;
                }
                rowcount++;
                if (NBperformed) {
                    CreateNBimages(width, height, "evaluation");
                }
            }
        } catch (Exception e) {
        }

        //read reconstruction parameters
        try {
            sheet = wb.getSheetAt(RECONPARA);
            rowIterator = sheet.iterator();
            while (rowIterator.hasNext()) {
                Row row = rowIterator.next();
                Cell cell = row.getCell(1); //Cell cell = row.getCell(1, Row.RETURN_BLANK_AS_NULL);
                reconpara.add((int) cell.getNumericCellValue());
            }

            /* Reconstruction Parameters
		 * 0: ACF Window								17: cormode
		 * 1: Intensity Trace Window					18: plotACFCurves
		 * 2: Residual Window							19: plotSDCurves
		 * 3: SD Window									20: plotIntensityCurves
		 * 4: MSD Window								21: plotResCurves
		 * 5: Difusion Law Window						22: plotParaHist
		 * 6: PSF Window								23: plotBlockingCurve
		 * 7: Para Cor Window							24: plotCovmats
		 * 8: Blocking Window							25: DCCF Window 0
		 * 9: Covariance Window							26: DCCF Window 1
		 * 10: N Window (N&B)							27: DCCF Window 2
		 * 11: B Window (N&B)							28: DCCF Window 3
		 * 12: Num Window (N&B)							29: DCCF Histogram Window 0
		 * 13: Epsilon Window (N&B)						30: DCCF Histogram Window 1
		 * 14: Parameter Window							31: DCCF Histogram Window 2
		 * 15: Histogram Window							32: DCCF Histogram Window 3
		 * 16: DL Map Window
             */
            //RECREATE WINDOWS
            plotACFCurves = (reconpara.get(18) != 0);
            plotSDCurves = (reconpara.get(19) != 0);
            plotIntensityCurves = (reconpara.get(20) != 0);
            plotResCurves = (reconpara.get(21) != 0);
            plotParaHist = (reconpara.get(22) != 0);
            plotBlockingCurve = (reconpara.get(23) != 0);
            plotCovmats = (reconpara.get(24) != 0);

            //create ACF window (Intensity trace, residual window, SD window are not recreated)
            if (reconpara.get(0) == 1) {
                Roi roi1 = new Roi(0, 0, width, height);
                roi1.setStrokeColor(java.awt.Color.BLUE);
                imp.setRoi(roi1);
                if (reconpara.get(14) == 1 && plotACFCurves) {		// if parameter window is open and plot ACF is enabled then plot also the fits
                    plotCF(roi1, reconpara.get(17), true);			// get the cormode from the file 
                } else {
                    plotCF(roi1, reconpara.get(17), false);
                }
            }

            //create MSD window
            if (reconpara.get(4) == 1) {
                Roi roi1 = new Roi(0, 0, width, height);
                roi1.setStrokeColor(java.awt.Color.BLUE);
                imp.setRoi(roi1);
                plotMSD(roi1, reconpara.get(17), true);
            }

            //plot Diffusion Law window
            if (reconpara.get(5) == 1) {
                if (difflawallbin > 1) {
                    minval = 0;
                    maxval = 0;
                    for (int u = 1; u <= difflawallbin; u++) {			// get minimum and maximum values for the plot				
                        if (difflaw[1][u - 1] + difflaw[2][u - 1] > maxval) {
                            maxval = difflaw[1][u - 1] + difflaw[2][u - 1];
                        }
                        if (difflaw[1][u - 1] - difflaw[2][u - 1] < minval) {
                            minval = difflaw[1][u - 1] - difflaw[2][u - 1];
                        }
                    }
                    plotDiffLaw(difflaw, difflawallbin, minvalDL, maxvalDL);
                    fitDiffLaw();
                }
            }

            //plot PSF window
            if (reconpara.get(6) == 1) {
                if (numofpsf > 1 && psfmaxbin > 1) {
                    minval = 0;
                    maxval = 0;
                    for (int x = 0; x < numofpsf; x++) {
                        for (int u = 0; u < psfmaxbin; u++) {
                            if (psfData[x][1][u] + psfData[x][2][u] > maxval) {
                                maxval = psfData[x][1][u] + psfData[x][2][u];
                            }
                            if (psfData[x][1][u] - psfData[x][2][u] < minval) {
                                minval = psfData[x][1][u] - psfData[x][2][u];
                            }
                        }
                    }
                    plotPSF(minval, maxval, numofpsf);
                }
            }

            //create Parameter window
            if (reconpara.get(14) == 1) {
                createParaImp((int) Math.floor((width - cfXDistance) / pixbinX), (int) Math.floor((height - cfYDistance) / pixbinY));
            }

            //plot Diffusion Law Map
            if (reconpara.get(16) == 1) {
                plotDiffLawMaps();
            }

        } catch (Exception e) {
            IJ.showMessage("Reconstruction Parameters not readable");
            throw (e);
        }

        IJ.log("Experiment read");
    }

    // bring all windows of this plugin to the front
    public void bringToFront() {
        Frame tframe = WindowManager.getFrame($impTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($acfWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($intWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($resWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($sdWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($msdWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($difflawWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($PSFWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        for (int i = 0; i < dccfMax; i++) {
            tframe = WindowManager.getFrame($histDCCFWinTitle[i]);
            if (tframe != null) {
                WindowManager.toFront(tframe);
            }
        }
        tframe = WindowManager.getFrame($dccfTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($paraCorWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($impPara1Title);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($histWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($blockingWindowTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($impCovTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($impNTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($impBTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($impNumTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($impEpsTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }
        tframe = WindowManager.getFrame($impDLMapTitle);
        if (tframe != null) {
            WindowManager.toFront(tframe);
        }

        WindowManager.toFront(fitframe);
        WindowManager.toFront(filteringframe);
        WindowManager.toFront(simframe);
        WindowManager.toFront(difflawframe);
        WindowManager.toFront(NBframe);
    }

    // EXIT the program 
    public void exitImFCS() {
        if (impwin != null && impwin.isClosed() == false) {	// remove listeners and overlays from image window
            if (imp.getOverlay() != null) {
                imp.deleteRoi();
                imp.getOverlay().clear();
                imp.setOverlay(imp.getOverlay());
            }
            impcan.removeMouseListener(impcan.getMouseListeners()[mlnum]);
            impcan.removeKeyListener(impcan.getKeyListeners()[klnum]);
        }
        closeWindows();
        frame.setVisible(false);	// close control panel
        frame.dispose();
        fitframe.setVisible(false);	// close fit panel
        fitframe.dispose();
        filteringframe.setVisible(false);	// close thresholding panel
        filteringframe.dispose();
        simframe.setVisible(false);	// close simulation panel
        simframe.dispose();
        difflawframe.setVisible(false);	// close diffusion law panel
        difflawframe.dispose();
        resframe.setVisible(false);	// close result table
        resframe.dispose();
        NBframe.setVisible(false);	// close NB panel
        NBframe.dispose();

    }

    public void closeWindows() {
        if (acfWindow != null && acfWindow.isClosed() == false) {	// close ACF window
            acfWindow.close();
        }
        if (intWindow != null && intWindow.isClosed() == false) {	// close intensity trace window
            intWindow.close();
        }
        if (resWindow != null && resWindow.isClosed() == false) {	// close fit residuals window
            resWindow.close();
        }
        if (sdWindow != null && sdWindow.isClosed() == false) {	// close SD window
            sdWindow.close();
        }
        if (msdWindow != null && msdWindow.isClosed() == false) {	// close MSD window
            msdWindow.close();
        }
        if (difflawWindow != null && difflawWindow.isClosed() == false) {	// close diffusion law window
            difflawWindow.close();
        }
        if (PSFWindow != null && PSFWindow.isClosed() == false) {	// close PSF window
            PSFWindow.close();
        }
        if (paraCorWindow != null && paraCorWindow.isClosed() == false) {	// close paraCor window
            paraCorWindow.close();
        }
        if (blockingWindow != null && blockingWindow.isClosed() == false) {	// close blocking window
            blockingWindow.close();
        }
        if (impCovWin != null && impCovWin.isClosed() == false) {	// close covariance window
            impCovWin.close();
        }
        if (impNWin != null && impNWin.isClosed() == false) {	// close N&B N window
            impNWin.close();
        }
        if (impBWin != null && impBWin.isClosed() == false) {	// close N&B B window
            impBWin.close();
        }
        if (impNumWin != null && impNumWin.isClosed() == false) {	// close N&B n window
            impNumWin.close();
        }
        if (impEpsWin != null && impEpsWin.isClosed() == false) {	// close N&B epsilon window
            impEpsWin.close();
        }
        for (int x = 0; x < dccfMax; x++) {
            if (impDCCFWin[x] != null && impDCCFWin[x].isClosed() == false) {	// close dCCF windows
                impDCCFWin[x].close();
            }
        }
        for (int x = 0; x < dccfMax; x++) {
            if (histDCCFWin[x] != null && histDCCFWin[x].isClosed() == false) {	// close dCCF historgram windows
                histDCCFWin[x].close();
            }
        }
        if (impPara1 != null && impPara1.isVisible()) {		// close parameter map window
            impPara1.close();
        }

        if (histWin != null && histWin.isClosed() == false) {		// close histogram windows
            histWin.close();
        }

        if (impDLMap != null && impDLMap.isVisible()) {				// close parameter map window
            impDLMap.close();
        }
    }

    /*
 *  SwingWorkers for the background threads of longer tasks
 *  
 *  correlateRoiWorker extends SwingWorker<Void, Void>: Correlate ROI or All pixels; perform fit and MSD calculation if chosen
 *  batchWorker extends SwingWorker<Void, Void> : treat a number of files as chosen by the user
 *  correlatePSFWorker extends SwingWorker<Void, Void>: Calcualtes data for PSF plot
 *  correlateDiffLawWorker extends SwingWorker<Void, Void>: Calculates data for Diffusion Law plot
 *  correlateDCCFWorker extends SwingWorker<Void, Void>: Creates one of four possible DCCF (difference of cross-correlation) plots 
 *  
     */
    public class correlateRoiWorker extends SwingWorker<Void, Void> {

        private void failIfInterrupted() throws InterruptedException {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interrupted while calcualting CFs");
            }
        }

        private final Roi currentROI;

        public correlateRoiWorker(final Roi currentROI) {
            this.currentROI = currentROI;
        }

        @Override
        protected Void doInBackground() throws Exception {
            correlateROI(currentROI);
            return null;
        }
    }

    public class batchWorker extends SwingWorker<Void, Void> {

        private void failIfInterrupted() throws InterruptedException {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interrupted while performing batch analysis");
            }
        }

        @Override
        protected Void doInBackground() throws Exception {
            if (batchDialogue()) {
                batchProcessing();
            }
            return null;
        }
    }

    public class correlatePSFWorker extends SwingWorker<Void, Void> {

        private void failIfInterrupted() throws InterruptedException {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interrupted while calcualting CFs");
            }
        }

        @Override
        protected Void doInBackground() throws Exception {
            determinePSF();
            return null;
        }
    }

    public class correlateDiffLawWorker extends SwingWorker<Void, Void> {

        private void failIfInterrupted() throws InterruptedException {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interrupted while calcualting CFs");
            }
        }

        @Override
        protected Void doInBackground() throws Exception {
            determineDiffusionLaw();
            return null;
        }
    }

    public class correlateDCCFWorker extends SwingWorker<Void, Void> {

        private void failIfInterrupted() throws InterruptedException {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interrupted while calcualting CFs");
            }
        }

        private final int mode;
        private final int wx;
        private final int hy;

        public correlateDCCFWorker(final int mode, final int wx, final int hy) {
            this.mode = mode;
            this.wx = wx;
            this.hy = hy;
        }

        @Override
        protected Void doInBackground() throws Exception {
            dccf(mode, wx, hy);
            return null;
        }
    }

    /*
 *  The following procedures check that correct images are open and paramters are meaningful and partly define reuqired parameters for the correlation
 *  
 *  public void obtainImage(): use existing image or load new image and set essential variables and arrays for this image; sets parameter setImp to true infroming the program that an appropriate image is available
 *  public void check(ImagePlus image): checks whether the open file is a stack with 16 bit gray values
 *  public void checkImp(): checks whether an image was already loaded or assigned for the plugin
 *  public boolean loadBGRFile(): load a bckground file that was recorded under the same conditions and in the same ROI
 *  public boolean setParameters(): checks and sets the necessary parameters for the plugin from the ControlPanel; this will be called each time an action on the panel is performed
 *  public void initializeArrays(): initializes arrays if essential parameters have been changed or image loaded; is called from setParameters() and obtainImage()
 *  public void initializeFitres(int a, int b, int c, int d): initialize Fitres and CCFq with NaN
 *  public void initializepixvalid(int a, int b, int c): initialize filtering mask with Double.NaN
 *  public void initializedCCF(int a, int b, int c): initializes array for dCCF with NaN and closes dCCF windows and histograms if they exist
 *  public void setDLparameters(): set the parametesr for the diffusion law
 *  public void performCFE(int px, int py): performs correlation and fit on a particular pixel as chosen by coordinates in the panel, mouse click, or key events
 *  public void performNB(ImagePlus tmpip): performs N&B analysis
 *  public void CreateNBimages(int x, int y): creates Number and Brightness maps from the results of N&B analysis
 *  public void camCalibrate(): adds a backgroud image and calculates the parameter S from a B vs epsilon' plot.
 *  public int minDetermination(ImagePlus image): determines minimum value in an image
 *  public void determinePSF(): calculates the diffusion law plots for vrarious PSF sizes to determine the PSF
 *  public void plotPSF(double minval, double maxval, int numofpsf)
 *  public void determineDiffusionLaw(): calculate the FCS diffusion law plot for an image
 *  public void plotDiffLaw(double[][] dldat, in dlbin, double minval, double maxval): plot the diffusion law
 *  public void correlateROI(ROI improi): correlate multiple pixels which were chose by an ROI or correlate all pixels (ROI=whole image)
 *  public void batchProcessing(): starts the batch processing of several files according to the evlauations chose in the batchDialog
 *  public void plotScatterPlot(): plot a scatter plot for a selected parameter pair (not all pairs possible at the moment)
 *  public void dccf(int mode, int wx, int hy): calcualtes on of four possible dCCF images
 *  public void createDCCFWindow(int dw, int dh, String $imagetitle, String $histtitle): creates the windows for dccf; is also used to recreate an experiment after readExperiment
 *  
     */
    public void obtainImage() {

        if (impwin != null && impwin.isClosed() == false) {	// remove listeners and overlays from image window if a previous one existed
            if (imp.getOverlay() != null) {
                imp.deleteRoi();
                imp.getOverlay().clear();
                imp.setOverlay(imp.getOverlay());
            }
            impcan.removeMouseListener(impcan.getMouseListeners()[mlnum]);
            impcan.removeKeyListener(impcan.getKeyListeners()[klnum]);
        }

        if (!simFile) {	// remember image path if it's not a simulation
            try {
                $imagePath = imp.getOriginalFileInfo().directory;
                $fileName = imp.getOriginalFileInfo().fileName;
            } catch (NullPointerException ne) {
                IJ.log("No file path resolved");
            } // avoids exception if a simulation file is re-used
        } else {
            simFile = false; // reset simFile for the next loaded picture
        }

        // check file format
        check(imp);

        // get image parameters
        impwin = imp.getWindow();		// get image window
        //impip = imp.getStack().getProcessor(1);		// processor
        impcan = imp.getCanvas();		// and canvas
        impcan.setFocusable(true);		// set focusable
        $impTitle = imp.getTitle();		// get title of image

        // add listeners
        impcan.addMouseListener(impcanMouseUsed);
        mlnum = impcan.getMouseListeners().length - 1;	// remember the listener
        impcan.addKeyListener(impcanKeyUsed);
        klnum = impcan.getKeyListeners().length - 1;	// remember the listener

        // get width, height, number of frames in stack, and magnification of the image window
        width = imp.getWidth();
        height = imp.getHeight();
        frames = imp.getStackSize();

        // check whether a background file has been loaded and whether this is to be used 
        if (bgrloaded && tbBGR.getText().equals("Bgr FILE")) {
            if (width != bgrw || height != bgrh || frames != bgrf) {
                if (impwin != null && impwin.isClosed() == false) {	// remove listeners and overlays from image window if a previous one existed
                    if (imp.getOverlay() != null) {
                        imp.deleteRoi();
                        imp.getOverlay().clear();
                        imp.setOverlay(imp.getOverlay());
                    }
                    impcan.removeMouseListener(impcan.getMouseListeners()[mlnum]);
                    impcan.removeKeyListener(impcan.getKeyListeners()[klnum]);
                }
                impwin.close();
                IJ.showMessage("Image is not the same size as Background. Image not loaded");
                return;
            }
            impmin = 0;
        } else {
            impmin = minDetermination(imp); //calculate the minimum of the image stack; this will be used as the default background value
        }

        // set position of image window
        impwin.setLocation(impwinPosX, impwinPosY);

        // enlarge image to better see pixels
        if (width >= height) {
            scimp = zoomFactor / width;	//adjustable: zoomFactor is by default 250 (see parameter definitions), a value chosen as it produces a good size on the screen
        } else {
            scimp = zoomFactor / height;
        }
        if (scimp < 1.0) {
            scimp = 1.0;
        }
        scimp *= 100;			// transfrom this into %tage to run ImageJ command
        IJ.run(imp, "Original Scale", "");
        IJ.run(imp, "Set... ", "zoom=" + scimp + " x=" + (int) Math.floor(width / 2) + " y=" + (int) Math.floor(height / 2));
        IJ.run("In [+]", ""); 	// This needs to be used since ImageJ 1.48v to set the window to the right size; 
        // this might be a bug and is an ad hoc solution for the moment; before only the "Set" command was necessary

        // set file parameters in control window
        tfLastFrame.setText(Integer.toString(frames));
        tfBackground.setText(Integer.toString(impmin));
        tfBackground2.setText(Integer.toString(impmin));
        slidingWindowLength = (int) Math.floor(frames / 20);

        // indicate that image was loaded
        setImp = true;

        //no need to warn the user that results will be reset
        askOnRewrite = false;

        // define array for this picture which will be used to store the results; as these arrays should stay with the window they are defined here; 
        // other changeable arrays are defined in setParameters(); however, if parameters are changed the arrays are re-initialized.
        //  chanum = (int) (16 + (Integer.parseInt($CorrelQ) - 1) * 8 + 1); 	// value of chanum based on default Q and P, used for arrays initialization
        // NOTE: calculation of chanum has been revised. It is now based on values/settings on GUI. Instead of default values as the commented line above.
        // Otherwise, there may be a downstream error, ie. array out of bound error, example in plotCF function, cormode 2, fitacf array
        String $bcmode = (String) cbBleachCor.getSelectedItem();
        if ("Sliding Window".equals($bcmode)) {
            lagnum = (int) Math.floor((Math.log(slidingWindowLength / (swMinFrameReq + correlatorp)) + 1) / Math.log(2));
            if (lagnum < correlatorq) {	// allow smaller correlatorq values as minimum but not larger
                correlatorq = lagnum;
                tfCorrelatorQ.setText(Integer.toString(lagnum));
            } else {
                lagnum = (int) correlatorq;
            }
            chanum = (int) (correlatorp + (lagnum - 1) * correlatorp / 2 + 1);
        } else {
            chanum = (int) (correlatorp + (correlatorq - 1) * correlatorp / 2 + 1);
            lagnum = (int) correlatorq;
        }
        tfFitEnd.setText(Integer.toString(chanum - 1));
        fitstart = Integer.parseInt(tfFitStart.getText());
        fitend = Integer.parseInt(tfFitEnd.getText());
        initializeArrays();

        //initialize arrays for fitting
        paramfit = new boolean[noparam];
        paraminitval = new double[noparam];

        //initialize arrays and variables for N&B
        NBN = new double[width][height];
        NBB = new double[width][height];
        NBNum = new double[width][height];
        NBEps = new double[width][height];
        NBslope = 1.0;
        NBperformed = false;
        NBcorrected = false;

        // re-initialize arrays for PSF data, diffusion law and diff. law fits; close respective windows if they exist
        psfData = new double[1][3][1];
        psfmaxbin = 1;
        numofpsf = 1;
        if (PSFWindow != null && !PSFWindow.isClosed()) {
            PSFWindow.close();
        }
        difflaw = new double[3][1];
        difflawbin = 1;
        if (difflawWindow != null && !difflawWindow.isClosed()) {
            difflawWindow.close();
        }
        ///YYY(difflawmap init)
        if (impDLMap != null) {
            impDLMap.close();		// close diffusion law map window if it exists
        }
        difflawarray = new double[1][1][3][1];
        diffLawFitMap = new double[1][1][2];
        diffLawMapwidth = 1;
        diffLawMapheight = 1;
        difflawmapbin = 1;
        difflawallbin = 1;
        for (int k = 0; k < 2; k++) {
            difflawfit[k] = 0.0;
            diffLawFitLim[k] = 0;
        }

        // change ImFCS panel title
        frame.setTitle("ImFCS - " + $impTitle);
        fitframe.setTitle("ImFCS Fitting - " + $impTitle);

        // set window titles
        $acfWindowTitle = "CF - " + $impTitle;
        $intWindowTitle = "Int - " + $impTitle;
        $resWindowTitle = "Res - " + $impTitle;
        $sdWindowTitle = "StdDev - " + $impTitle;
        $msdWindowTitle = "MSD - " + $impTitle;
        $difflawWindowTitle = "DiffLaw - " + $impTitle;
        $PSFWindowTitle = "PSF - " + $impTitle;
        $paraCorWindowTitle = "Scatter - " + $impTitle;
        $impPara1Title = "Maps - " + $impTitle;
        $impDLMapTitle = "DL Map - " + $impTitle;
        $dccfTitle = "dCCF - " + $impTitle;
        $histWindowTitle = $impTitle;
        $blockingWindowTitle = "Blocking Transform - " + $impTitle;
        $impCovTitle = "Covariance - " + $impTitle;
        $impNTitle = "Number (N&B) - " + $impTitle;
        $impBTitle = "Brightness (N&B) - " + $impTitle;
        $impNumTitle = "Num - " + $impTitle;
        $impEpsTitle = "Epsilon - " + $impTitle;
    }

    // Check image type
    public void check(ImagePlus image) {
        // slice numbers start with 1 for historical reasons
        int type = image.getType();
        if (type != ImagePlus.GRAY16) {
            IJ.showMessage("Wrong Image Format", "Only GRAY16 Tiff stacks supported");
            throw new RuntimeException("Only GRAY16 Tiff stacks supported");
        }
    }

    // check whether an image was either loaded or an existing image assigned to the plugin for treatment
    public void checkImp() {
        if (setImp == false) {
            IJ.showMessage("Image not loaded or assigned.\nPlease use \"Load\" or \"Use\" button.");
            throw new RuntimeException("No Image loaded/assigned");
        }
    }

    // get background characteristics that were recorded under the same conditions and in the same ROI
    public boolean loadBGRFile() {
        // open a background image
        String $bgrTitle = $bgrTitleMem;
        int returnVal;
        double[][] mean;
        double[][] mean2;

        // get the background file
        JFileChooser fc = new JFileChooser($imagePath);
        fc.setSelectedFile(new File($bgrTitle));
        fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
        returnVal = fc.showOpenDialog(btnLoad);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            bgrimp = IJ.openImage(fc.getSelectedFile().getAbsolutePath());
            if (bgrimp == null) {
                JOptionPane.showMessageDialog(null, "Selected file does not exist or it is not an image");
                return false;
            }
            check(bgrimp);
            bgrWin = bgrimp.getWindow();
            $bgrTitle = fc.getName(fc.getSelectedFile());
            // set the size variables for the bgr window
            bgrf = bgrimp.getStackSize();
            bgrw = bgrimp.getWidth();
            bgrh = bgrimp.getHeight();
            if (impwin != null && impwin.isClosed() == false) {
                if (width != bgrw || height != bgrh || frames != bgrf) {
                    IJ.showMessage("Background image is not the same size as Image. Background image not loaded.");
                    if (bgrWin != null && bgrWin.isClosed() == false) {	// close N&B B window
                        bgrWin.close();
                    }
                    return false;
                }
            }
            //bgrimp.show();
        } else {
            JOptionPane.showMessageDialog(null, "No background image loaded.");
            return false;
        }

        //arrays for the calculation of row and column means
        bgrrowmean = new double[bgrh];
        bgrcolumnmean = new double[bgrw];

        //get background statistics for the fits
        ImageStatistics bgrStat = bgrimp.getStatistics();
        bgrframemean = bgrStat.mean;

        // calculate the row means
        for (int z = 1; z <= bgrf; z++) {
            bgrip = bgrimp.getStack().getProcessor(z);
            for (int y = 0; y < bgrh; y++) {
                double sum = 0;
                for (int x = 0; x < bgrw; x++) {
                    sum += bgrip.getPixel(x, y);
                }
                bgrrowmean[y] += sum / (bgrw * bgrf);
            }
        }

        // calculate column means
        for (int z = 1; z <= bgrf; z++) {
            bgrip = bgrimp.getStack().getProcessor(z);
            for (int x = 0; x < bgrw; x++) {
                double sum = 0;
                for (int y = 0; y < bgrh; y++) {
                    sum += bgrip.getPixel(x, y);
                }
                bgrcolumnmean[x] += sum / (bgrh * bgrf);
            }
        }

        // calculate covariance
        ImageProcessor bgrCorip;
        ImageProcessor bgrCorip2;
        mean = new double[bgrw][bgrh];
        mean2 = new double[bgrw][bgrh];
        bgrCoVar = new double[bgrw][bgrh];
        for (int z = 1; z < bgrf; z++) {
            bgrCorip = bgrimp.getStack().getProcessor(z);
            bgrCorip2 = bgrimp.getStack().getProcessor(z + 1);
            for (int x = 0; x < bgrw; x++) {
                for (int y = 0; y < bgrh; y++) {
                    mean[x][y] += bgrCorip.getPixel(x, y);
                    mean2[x][y] += bgrCorip2.getPixel(x, y);
                    bgrCoVar[x][y] += bgrCorip.getPixel(x, y) * bgrCorip2.getPixel(x, y);
                }
            }
        }

        for (int x = 0; x < bgrw; x++) {	// caluclate mean and variance: E(x) = Sum(x)/(n-1), E(x*x') = Sum(x*x')/(n-1) and coVar = E(x*x') - E(x)E(x') 
            for (int y = 0; y < bgrh; y++) {
                mean[x][y] /= (bgrf - 1);
                mean2[x][y] /= (bgrf - 1);
                bgrCoVar[x][y] = bgrCoVar[x][y] / (bgrf - 1) - mean[x][y] * mean2[x][y];
            }
        }

        // calculate variance
        bgrmean = new double[bgrw][bgrh];
        bgrVar = new double[bgrw][bgrh];
        for (int z = 1; z <= bgrf; z++) {
            bgrCorip = bgrimp.getStack().getProcessor(z);
            for (int x = 0; x < bgrw; x++) {
                for (int y = 0; y < bgrh; y++) {
                    bgrmean[x][y] += bgrCorip.getPixel(x, y);
                    bgrVar[x][y] += Math.pow(bgrCorip.getPixel(x, y), 2.0);
                }
            }
        }

        for (int x = 0; x < bgrw; x++) {	// caluclate mean and variance: E(x) = Sum(x)/n, E(x^2) = Sum(x^2)/n and Var = E(x^2) - E(x)^2 
            for (int y = 0; y < bgrh; y++) {
                bgrmean[x][y] /= bgrf;
                bgrVar[x][y] = bgrVar[x][y] / bgrf - Math.pow(bgrmean[x][y], 2.0);
            }
        }

        // set BGR values
        bgrStat = bgrimp.getStatistics();
        tfBackground.setText("0");
        $bgrTitleMem = $bgrTitle;
        background = 0;
        return true;
    }

    // checks and sets the necessary parameters for the plugin from the ControlPanel; 
    // this will be called each time an action on the panel is performed; as it resets
    // note that setParameters initializes all arrays and thus deletes all previous results 
    public boolean setParameters() {
        int index = 0;
        boolean resetResults = false; 								// whether Result arrays need to be reset
        boolean proceed = true;									// whether program should proceed resetting the Results 
        boolean onlySigmaOrBinChanged = true;							// whether sigma0 is the only parameter changed in the panel - PSF calibration is not reset in that case
        boolean onlyBinChanged = true;							// whether binning is the only parameter changed in the panel - diff law is not reset in that case
        String[] newPanelSettings = new String[noSettings];		// an array to temporarily hold the settings from the Panel

        checkImp();		// check if image was loaded or assigned

        // read settings from the panel and store them temporarily in newPanelSettings
        int tmpct = 0;
        newPanelSettings[tmpct++] = tfFirstFrame.getText();		// index 0
        newPanelSettings[tmpct++] = tfLastFrame.getText();
        newPanelSettings[tmpct++] = tfFrameTime.getText();
        String strbin = tfBinning.getText();
        int strlen = strbin.length();
        newPanelSettings[tmpct++] = strbin.substring(0, strbin.indexOf(" x "));
        newPanelSettings[tmpct++] = strbin.substring(strbin.indexOf(" x ") + 3);
        newPanelSettings[tmpct++] = tfCfXDistance.getText();
        newPanelSettings[tmpct++] = tfCfYDistance.getText();
        newPanelSettings[tmpct++] = (String) cbCorrelatorP.getSelectedItem();
        newPanelSettings[tmpct++] = tfCorrelatorQ.getText();
        newPanelSettings[tmpct++] = (String) cbFitModel.getSelectedItem();	// index 10
        newPanelSettings[tmpct++] = tbFCCSDisplay.getText();
        newPanelSettings[tmpct++] = tfPixelSize.getText();
        newPanelSettings[tmpct++] = tbOverlap.getText();
        newPanelSettings[tmpct++] = tfMagnification.getText();
        newPanelSettings[tmpct++] = tfNA.getText();
        newPanelSettings[tmpct++] = tfEmLambda.getText();
        newPanelSettings[tmpct++] = tfEmLambda2.getText();
        newPanelSettings[tmpct++] = tfSigma.getText();
        newPanelSettings[tmpct++] = tfSigmaZ.getText();
        newPanelSettings[tmpct++] = tfSigma2.getText();		// index = 20
        newPanelSettings[tmpct++] = tfSigmaZ2.getText();
        newPanelSettings[tmpct++] = tfBackground.getText();
        newPanelSettings[tmpct++] = tfBackground2.getText();
        newPanelSettings[tmpct++] = $bgrTitleMem;
        newPanelSettings[tmpct++] = (String) cbBleachCor.getSelectedItem();
        newPanelSettings[tmpct++] = Integer.toString(slidingWindowLength);
        newPanelSettings[tmpct++] = Integer.toString(polyOrder);
        newPanelSettings[tmpct++] = (String) cbFilter.getSelectedItem();
        newPanelSettings[tmpct++] = Integer.toString(filterUL);
        newPanelSettings[tmpct++] = Integer.toString(filterLL);		// index = 30
        newPanelSettings[tmpct++] = tbFit.getText();

        // check whether any settings have changed
        for (int i = 0; i < noSettings; i++) {
            if (!(newPanelSettings[i].equals(panelSettings[i])) && keyParam[i]) {
                resetResults = true;
                if (i != 17 && i != 19 && i != 3 && i != 4) {
                    onlySigmaOrBinChanged = false;
                }
                if (i != 3) {
                    onlyBinChanged = false;
                }
                if (i != 4) {
                    onlyBinChanged = false;
                }
            }
        }

        if (resetResults && askOnRewrite) {
            GenericDialog gd = new GenericDialog("Delete the Results and start new?");
            gd.addMessage("Some of the parameter settings in the main panel have been changed");
            gd.addMessage("Continuing will results in deleting some Results");
            gd.addMessage("Do you wish to proceed?");
            gd.showDialog();
            if (gd.wasOKed()) {
                proceed = true;
            }
            if (gd.wasCanceled()) {
                proceed = false;
            }
        }

        if (proceed) {
            // assign the values to the variables used in the calculations

            try {
                firstframe = Integer.parseInt(newPanelSettings[0]);
                lastframe = Integer.parseInt(newPanelSettings[1]);
                frametime = Double.parseDouble(newPanelSettings[2]);
                binningX = Integer.parseInt(newPanelSettings[3]);
                binningY = Integer.parseInt(newPanelSettings[4]);
                cfXDistance = Integer.parseInt(newPanelSettings[5]);
                cfYDistance = Integer.parseInt(newPanelSettings[6]);
                correlatorp = Integer.parseInt(newPanelSettings[7]);
                correlatorq = Integer.parseInt(newPanelSettings[8]);
                pixelsize = Double.parseDouble(newPanelSettings[11]);
                objmag = Double.parseDouble(newPanelSettings[13]);
                NA = Double.parseDouble(newPanelSettings[14]);
                emlambda = Double.parseDouble(newPanelSettings[15]);
                emlambda2 = Double.parseDouble(newPanelSettings[16]);
                sigma = Double.parseDouble(newPanelSettings[17]);
                sigmaZ = Double.parseDouble(newPanelSettings[18]);
                sigma2 = Double.parseDouble(newPanelSettings[19]);
                sigmaZ2 = Double.parseDouble(newPanelSettings[20]);
                background = Integer.parseInt(newPanelSettings[21]);
                background2 = Integer.parseInt(newPanelSettings[22]);
            } catch (NumberFormatException nfe) {
                IJ.showMessage("A parameter in the panel has an invalid format");
                throw new NumberFormatException("Number format error.");
            }

            if (cbFitModel.getSelectedItem() == "FCS") {
                cfXshift = cfXDistance;
                cfYshift = cfYDistance;
            } else {
                cfXshift = 0;
                cfYshift = 0;
            }

            // set maximum, and minimum cursor positions possible in the image, depending on binning and whether overlap is allowed
            // parameters need to be checked according to these settings as allowed parameter ranges differ in overlap and non-overlap mode
            if (overlap) {
                pixelWidthX = width - binningX;		// these values are the correct maximum if counting from 0
                pixelHeightY = height - binningY;
                pixbinX = 1;
                pixbinY = 1;
            } else {
                pixelWidthX = (int) Math.floor(width / binningX) - 1;
                pixelHeightY = (int) Math.floor(height / binningY) - 1;
                pixbinX = binningX;
                pixbinY = binningY;
            }

            // check parameter settings
            // check that numbers are not out of bounds and make sense
            if (firstframe < 1 || firstframe > frames || lastframe < 1 || lastframe > frames || firstframe >= lastframe) {
                JOptionPane.showMessageDialog(null, "Frames set incorrectly");
                return false;
            }

            if (binningX < 1 || binningY < 1) {
                JOptionPane.showMessageDialog(null, "Parameter \"binning\" is smaller than 1.");
                return false;
            }

            if (binningX > width || binningY > height) {
                JOptionPane.showMessageDialog(null, "Parameter \"binning\" is larger than image size.");
                return false;
            }

            // check that the distance between pixels to be correlated is smaller than the size of the image
            if (Math.abs(cfXDistance) > width - binningX || Math.abs(cfYDistance) > height - binningY) {
                JOptionPane.showMessageDialog(null, "Correlation distance is larger than image size.");
                return false;
            }

            // check that pixel values are within image
            if (checkroi) {
                if ((cfXDistance < 0 && Math.ceil((double) roi1StartX / pixbinX) * pixbinX < cfXDistance * (-1)) || (cfYDistance < 0 && Math.ceil((double) roi1StartY / pixbinY) * pixbinY < cfYDistance * (-1))) {
                    JOptionPane.showMessageDialog(null, "Correlation points are not within image.");
                    return false;
                }

                if ((cfXDistance >= 0 && Math.floor(((double) roi1StartX + roi1WidthX - binningX) / pixbinX) * pixbinX + cfXDistance > width - binningX) || (cfYDistance >= 0 && Math.floor(((double) roi1StartY + roi1HeightY - binningY) / pixbinY) * pixbinY + cfYDistance > height - binningY)) {
                    JOptionPane.showMessageDialog(null, "Correlation points are not within image.");
                    return false;
                }

                if (cbFitModel.getSelectedItem() == "DC-FCCS") { //this applies only for dual-color cross-correlations
                    // check that the correlation areas don't overlap
                    if (Math.abs(cfXDistance) < roi1WidthX && Math.abs(cfYDistance) < roi1HeightY) {
                        JOptionPane.showMessageDialog(null, "Cross-correlation areas overlap.");
                        return false;
                    }
                }
            }

            if (cfXDistance % pixbinX != 0 || cfXDistance % pixbinX != 0) {
                JOptionPane.showMessageDialog(null, "Warning: CF distance is not a multiple of pixel bin.");
            }

            // check that background1 value is sensible
            if (background < 0) {
                JOptionPane.showMessageDialog(null, "Invalid background 1 value");
                return false;
            } else if (background > impmin) {
                JOptionPane.showMessageDialog(null, "Warning: Background 1 value is larger than smallest pixel value");
            }

            // check that background2 value is sensible
            if (background2 < 0) {
                JOptionPane.showMessageDialog(null, "Invalid background 2 value");
                return false;
            } else if (background2 > impmin) {
                JOptionPane.showMessageDialog(null, "Warning: Background 2 value is larger than smallest pixel value");
            }

            // check whether correlator Q is sensible
            if (correlatorq < 1) {
                JOptionPane.showMessageDialog(null, "Correlator Q is smaller than 1");
                return false;
            }

            // check whether there are enough frames for the correlation and give a warning if there are less than minFrameReq frames for the last correlation channel
            if ((lastframe - firstframe + 1 - Math.pow(2, correlatorq - 1) * correlatorp / 2 - correlatorp / 4 * Math.pow(2, correlatorq)) / Math.pow(2, correlatorq - 1) < 1) {
                JOptionPane.showMessageDialog(null, "Not enough frames");
                return false;
            } else if ((lastframe - firstframe + 1 - Math.pow(2, correlatorq - 1) * correlatorp / 2 - correlatorp / 4 * Math.pow(2, correlatorq)) / Math.pow(2, correlatorq - 1) < minFrameReq) {
                JOptionPane.showMessageDialog(null, "Warning: Less than " + minFrameReq + " data point for last correlation channel.");
            }

            // set common arrays and parameters according to user settings in the panel
            if ((lastframe - firstframe + 1) < 1000) {	// use 1000 points for the intensity, except when less than 1000 frames are present
                nopit = (lastframe - firstframe + 1);
            } else {
                nopit = 1000;
            }

            // if sliding window correction is needed then the correlator structure has to be adapted to the smaller number of lagtimes
            int num; 	// total number of frames to be correlated; sliding window length or lastframe-firstframe+1
            String $bcmode = (String) cbBleachCor.getSelectedItem();
            if ("Sliding Window".equals($bcmode)) {
                lagnum = (int) Math.floor((Math.log(slidingWindowLength / (swMinFrameReq + correlatorp)) + 1) / Math.log(2));
                if (lagnum < correlatorq) {	// allow smaller correlatorq values as minimum but not larger
                    correlatorq = lagnum;
                    tfCorrelatorQ.setText(Integer.toString(lagnum));
                } else {
                    lagnum = (int) correlatorq;
                }
                chanum = (int) (correlatorp + (lagnum - 1) * correlatorp / 2 + 1);
                num = slidingWindowLength;
            } else {
                chanum = (int) (correlatorp + (correlatorq - 1) * correlatorp / 2 + 1);
                lagnum = (int) correlatorq;
                num = (lastframe - firstframe + 1);
            }

            // initialize arrays required for calculations; they change with each new paramter setting and are thus re-initialized
            intTrace1 = new double[nopit];
            intTrace2 = new double[nopit];
            intTime = new double[nopit];

            // the arrays for CFs and results of fitting are reinitialized only if key parameters have changed to prevent mixing results for different settings in a single parameter map
            if (resetResults) {
                initializeArrays();
                tfFitEnd.setText(Integer.toString(chanum - 1));			// reset the fitting range for data
                fitstart = Integer.parseInt(tfFitStart.getText());
                fitend = Integer.parseInt(tfFitEnd.getText());
                if (!onlyBinChanged && !expload) {					// if other parameter than binning has been changed, reset diff. law data and close diff. law window if it exists; but don't do anything if htis is an experimetn load
                    difflaw = new double[3][1];
                    difflawbin = 1;
                    if (difflawWindow != null && !difflawWindow.isClosed()) {
                        difflawWindow.close();
                    }
                    //YYY(close difflawmap and initializearrays)
                    if (impDLMap != null) {
                        impDLMap.close();		// close diffusion law map window if it exists
                    }
                    difflawarray = new double[1][1][3][1];
                    diffLawFitMap = new double[1][1][2];
                    diffLawMapwidth = 1;
                    diffLawMapheight = 1;
                    difflawmapbin = 1;
                    difflawallbin = 1;
                    for (int k = 0; k < 2; k++) {
                        difflawfit[k] = 0.0;
                        diffLawFitLim[k] = 0;
                    }
                    if (!onlySigmaOrBinChanged) {				// if other parameter than sigma or bin have been changed, reset PSF data and close PSF window if it exists
                        psfData = new double[1][3][1];
                        psfmaxbin = 1;
                        numofpsf = 1;
                        if (PSFWindow != null && !PSFWindow.isClosed()) {
                            PSFWindow.close();
                        }
                    }
                }
            }

            // initialize arrays required for calculations; they change with each new paramter setting				
            base = (int) correlatorp; 		// base = number of channels in first group
            hbase = (int) correlatorp / 2; 	// hbase = number of channels in all higher groups
            mtab = new double[chanum]; 		// number of samples for each correlation channel
            lag = new int[chanum];			// lag for each correlation channel; indepenednet of time
            samp = new double[chanum];		// sampletime (or bin width) of each channel
            lagtime = new double[chanum];	// lagtime = lag*frametime; this is the actual lagtime in seconds for each channel

            for (int x = 0; x <= hbase; x++) {	// calculate lag and lagtimes for the 0 lagtime channel and the first 8 channels
                lag[x] = x;
                lagtime[x] = lag[x] * frametime;
            }

            for (int x = 1; x <= lagnum; x++) {	// calculate lag and lagtimes for all higher channels
                for (int y = 1; y <= hbase; y++) {
                    lag[x * hbase + y] = (int) (Math.pow(2, x - 1) * y + (base / 4) * Math.pow(2, x));
                    lagtime[x * hbase + y] = lag[x * hbase + y] * frametime;
                }
            }

            for (int x = 0; x <= base; x++) {	// calculate sampletimes (bin width) for the 0 lagtime channel and the first 8 channels
                samp[x] = 1;
            }

            for (int x = 2; x <= lagnum; x++) {	// calculate sampletimes (bin width) for all higher channels
                for (int y = 1; y <= hbase; y++) {
                    samp[x * hbase + y] = Math.pow(2, x - 1);
                }
            }

            // calculate the number of samples for each channel including 0 lagtime; this differs for sliding window correction
            // the variable num takes care of this
            for (int x = 0; x <= (chanum - 1); x++) {
                mtab[x] = (int) Math.floor((num - lag[x]) / samp[x]);
            }

            // set initial, maximum, and minimum cursor positions possible in the image
            if (cfXDistance >= 0) {
                maxcposx = pixelWidthX - (int) Math.ceil(((double) cfXDistance - (width - (pixelWidthX * pixbinX + binningX))) / pixbinX);
                mincposx = 0;
            } else {
                maxcposx = pixelWidthX;
                mincposx = -(int) Math.floor((double) cfXDistance / pixbinX);
            }

            if (cfYDistance >= 0) {
                maxcposy = pixelHeightY - (int) Math.ceil(((double) cfYDistance - (height - (pixelHeightY * pixbinY + binningY))) / pixbinY);
                mincposy = 0;
            } else {
                maxcposy = pixelHeightY;
                mincposy = -(int) Math.floor((double) cfYDistance / pixbinY);
            }

            // set values in the fit window
            tfParama.setText(IJ.d2s(pixelsize * 1000 / objmag * binningX));	// XXX
            tfParamw.setText(IJ.d2s(sigma * emlambda / NA, decformat));
            tfParamw2.setText(IJ.d2s(sigma2 * emlambda2 / NA, decformat));
            tfParamz.setText(IJ.d2s(sigmaZ * emlambda / NA, decformat));
            tfParamz2.setText(IJ.d2s(sigmaZ2 * emlambda2 / NA, decformat));
            pixeldimx = (pixelsize * 1000 / objmag * binningX) / Math.pow(10, 9);
            pixeldimy = (pixelsize * 1000 / objmag * binningY) / Math.pow(10, 9);
            psfsize = (sigma * emlambda / NA) / Math.pow(10, 9);
            psfsize2 = (sigma2 * emlambda2 / NA) / Math.pow(10, 9);
            lsthickness = (sigmaZ * emlambda / NA) / Math.pow(10, 9);
            lsthickness2 = (sigmaZ2 * emlambda2 / NA) / Math.pow(10, 9);

            tfParamRx.setText(IJ.d2s(pixelsize * 1000 / objmag * cfXshift, decformat)); // set the value of rx
            tfParamRy.setText(IJ.d2s(pixelsize * 1000 / objmag * cfYshift, decformat)); // set the value of ry

            // Rz is set to 0 at the moment, no fit possible. Can be added later.
            // read initial values from Fit window only when they are to be refreshed
            if (tbFixPar.getText().equals("Free")) {
                initparam[0] = Double.parseDouble(tfParamN.getText());
                initparam[1] = Double.parseDouble(tfParamD.getText()) / Math.pow(10, 12);
                initparam[2] = Double.parseDouble(tfParamVx.getText()) / Math.pow(10, 6);
                initparam[3] = Double.parseDouble(tfParamVy.getText()) / Math.pow(10, 6);
                initparam[4] = Double.parseDouble(tfParamG.getText());
                initparam[5] = Double.parseDouble(tfParamF2.getText());
                initparam[6] = Double.parseDouble(tfParamD2.getText()) / Math.pow(10, 12);
                initparam[7] = Double.parseDouble(tfParamF3.getText());
                initparam[8] = Double.parseDouble(tfParamD3.getText()) / Math.pow(10, 12);
                initparam[9] = Double.parseDouble(tfParamFtrip.getText());
                initparam[10] = Double.parseDouble(tfParamTtrip.getText()) / Math.pow(10, 6);
            }

            // save the new settings in panelSettings
            System.arraycopy(newPanelSettings, 0, panelSettings, 0, noSettings);
            askOnRewrite = true;

            return true;
        } else {
            return false;
        }
    }

    // initializes arrays if certain parameters are changed and is called from setParameters() 
    // the plugin thus does not allow to have reuslts from different settings in one parameter map
    public void initializeArrays() {
        acf = new double[3][width][height][chanum];
        varacf = new double[3][width][height][chanum];
        sdacf = new double[3][width][height][chanum];
        fitacf = new double[3][width][height][chanum];
        res = new double[3][width][height][chanum];
        msd = new double[3][width][height][chanum];
        blocked = new double[3][width][height];
        aveacf = new double[chanum];
        fitaveacf = new double[chanum];
        varaveacf = new double[chanum];
        msdaveacf = new double[chanum];
        resaveacf = new double[chanum];
        initializedCCF(dccfMax, width, height);		//initialize array for dCCF
        initializeFitres(3, width, height, noparam);	// initialize fitres and pixfiltered
        initializepixvalid(3, width, height);			// initialize filtering mask
        userThreshold[2] = false;						// tresholds haven't been applied on the data
        if (impPara1 != null && impPara1.isVisible()) {		//close parameter and histograms windows if they exists
            impPara1.close();
        }
        if (histWin != null && !histWin.isClosed()) {
            histWin.close();
        }
    }

    public void initializeFitres(int a, int b, int c, int d) {    //initializes arrays for fit parameters, chi2 and cross-correlation amount with NaN
        fitres = new double[a][b][c][d];
        chi2 = new double[a][b][c];
        CCFq = new double[b][c];
        pixfitted = new boolean[a][b][c];
        for (int q = 0; q < a; q++) {
            for (int r = 0; r < b; r++) {
                for (int s = 0; s < c; s++) {
                    chi2[q][r][s] = Double.NaN;
                    pixfitted[q][r][s] = false;
                    for (int t = 0; t < d; t++) {
                        fitres[q][r][s][t] = Double.NaN;
                    }
                }
            }
        }
        for (int r = 0; r < b; r++) {
            for (int s = 0; s < c; s++) {
                CCFq[r][s] = Double.NaN;
            }
        }
    }

    public void initializepixvalid(int a, int b, int c) {
        pixvalid = new double[a][b][c];
        for (int q = 0; q < a; q++) {
            for (int r = 0; r < b; r++) {
                for (int s = 0; s < c; s++) {
                    pixvalid[q][r][s] = Double.NaN;
                }
            }
        }
    }

    public void initializedCCF(int a, int b, int c) {
        dccf = new double[a][b][c];
        for (int q = 0; q < a; q++) {
            if (impDCCFWin[q] != null && !impDCCFWin[q].isClosed()) {
                impDCCFWin[q].close();
            }
            if (histDCCFWin[q] != null && !histDCCFWin[q].isClosed()) {
                histDCCFWin[q].close();
            }
            dccfCalculated[q] = false;
            for (int r = 0; r < b; r++) {
                for (int s = 0; s < c; s++) {
                    dccf[q][r][s] = Double.NaN;
                }
            }
        }
    }

    public void setDLparameters() {
        if (overlap) {						// determine the number of points in the diffusion law plot
            if (width <= height) {
                tfDLBinEnd.setText(Integer.toString(width + 1 - minDLPoints));
                tfDLFitEnd.setText(Integer.toString(width + 1 - minDLPoints));
            } else {
                tfDLBinEnd.setText(Integer.toString(height + 1 - minDLPoints));
                tfDLFitEnd.setText(Integer.toString(height + 1 - minDLPoints));
            }
        } else {
            if (width <= height) {
                tfDLBinEnd.setText(Integer.toString((int) Math.floor(width / minDLPoints)));
                tfDLFitEnd.setText(Integer.toString((int) Math.floor(width / minDLPoints)));
            } else {
                tfDLBinEnd.setText(Integer.toString((int) Math.floor(height / minDLPoints)));
                tfDLFitEnd.setText(Integer.toString((int) Math.floor(height / minDLPoints)));
            }
        }
    }

    //perform correlations and fits on a particular pixel as chosen by coordinates in the panel, mouse click, or key events
    public void performCFE(int px, int py) {
        // check that pixel values given by the user (as read in setParameters()) are within image
        // px, and py are the pixel positions in the imp
        double q1;
        double q2;

        if (cbFitModel.getSelectedItem() == "DC-FCCS" && (Math.abs(cfXDistance) < binningX && Math.abs(cfYDistance) < binningY)) {		//in DC-FCCS mode check that cross-correlation areas don't overlap
            JOptionPane.showMessageDialog(null, "Cross-correlation areas overlap.");
        } else {

            if (px <= maxcposx && px >= mincposx && py <= maxcposy && py >= mincposy) {
                cposx = px * pixbinX;
                cposy = py * pixbinY;
                c2posx = px * pixbinX + cfXDistance;
                c2posy = py * pixbinY + cfYDistance;

                //set the ROI in the image
                if (imp.getOverlay() != null) {
                    imp.getOverlay().clear();
                    imp.setOverlay(imp.getOverlay());
                }
                Roi roi1 = new Roi(cposx, cposy, binningX, binningY);
                roi1.setStrokeColor(java.awt.Color.BLUE);
                imp.setRoi(roi1);

                //set a second ROI in cross-correlation mode
                if (cfXDistance + cfYDistance != 0 || cbFitModel.getSelectedItem() == "DC-FCCS") {
                    Roi roi2 = new Roi(c2posx, c2posy, binningX, binningY);
                    roi2.setStrokeColor(java.awt.Color.RED);
                    Overlay impov = new Overlay(roi2);
                    imp.setOverlay(impov);
                }

                if (cbFitModel.getSelectedItem() == "FCS") {
                    // calculate intensities and correlations for the point indicated in the control panel
                    // note that calcIntensityTrace() has to come before correlate() as this trace is also used for bleach correction
                    calcIntensityTrace(imp, cposx, cposy, c2posx, c2posy, firstframe, lastframe);
                    correlate(imp, cposx, cposy, c2posx, c2posy, 0, firstframe, lastframe);

                    //if fit is selected, perform the fit
                    if (doFit) {
                        prepareFit();
                        fit(px, py, 0, "FCS");
                        if (impPara1 != null) {		//update parameter window if it exists
                            updateParaImp();
                        } else {
                            createParaImp(maxcposx - mincposx + 1, maxcposy - mincposy + 1);
                        }
                    }

                    //plot the correlation function and the intensity
                    plotCF(roi1, 1, false);
                    plotIntensityTrace(roi1, 1);
                    if (doMSD) { 					// plot MSD if selected
                        plotMSD(roi1, 1, false);
                    }
                }

                if (cbFitModel.getSelectedItem() == "DC-FCCS" && !tbFCCSDisplay.isSelected()) {

                    // calculate intensities and correlations for the point indicated in the control panel
                    // note that calcIntensityTrace() has to come before correlate() as this trace is also used for blaech correction
                    calcIntensityTrace(imp, cposx, cposy, c2posx, c2posy, firstframe, lastframe);
                    correlate(imp, cposx, cposy, c2posx, c2posy, 2, firstframe, lastframe);

                    //if fit is selected, perform the fit
                    if (doFit) {
                        prepareFit();
                        fit(px, py, 2, "DC-FCCS");
                        if (impPara1 != null) {		//update parameter window if it exists
                            updateParaImp();
                        } else {
                            createParaImp(maxcposx - mincposx + 1, maxcposy - mincposy + 1);
                        }
                    }

                    //plot the correlation function and the intensity
                    plotCF(roi1, 1, false);
                    plotIntensityTrace(roi1, 1);
                    if (doMSD) { 					// plot MSD if selected
                        plotMSD(roi1, 1, false);
                    }
                }

                if (cbFitModel.getSelectedItem() == "DC-FCCS" && tbFCCSDisplay.isSelected()) {
                    // calculate intensities and correlations for the point indicated in the control panel
                    // note that calcIntensityTrace() has to come before correlate() as this trace is also used for exponential bleach correction
                    calcIntensityTrace(imp, cposx, cposy, cposx, cposy, firstframe, lastframe);
                    correlate(imp, cposx, cposy, cposx, cposy, 0, firstframe, lastframe);
                    calcIntensityTrace(imp, c2posx, c2posy, c2posx, c2posy, firstframe, lastframe);
                    correlate(imp, c2posx, c2posy, c2posx, c2posy, 1, firstframe, lastframe);
                    calcIntensityTrace(imp, cposx, cposy, c2posx, c2posy, firstframe, lastframe);
                    correlate(imp, cposx, cposy, c2posx, c2posy, 2, firstframe, lastframe);

                    //if fit is selected, perform the fit
                    if (doFit) {
                        prepareFit();
                        fit(px, py, 0, "FCS");
                        prepareFit();
                        fit(px, py, 1, "FCS");
                        prepareFit();
                        fit(px, py, 2, "DC-FCCS");

                        //calculate q value
                        q1 = fitres[0][px][py][0] / fitres[2][px][py][0];
                        q2 = fitres[1][px][py][0] / fitres[2][px][py][0];
                        if (q1 > q2) {
                            CCFq[px][py] = q1;
                        } else {
                            CCFq[px][py] = q2;
                        }
                        if (CCFq[px][py] > 1) {
                            CCFq[px][py] = Double.NaN;	//cross-correlation > 1 is non-physical, most likely due to very noisy CFs, better discarded
                        }
                        if (impPara1 != null) {		//update parameter window if it exists
                            updateParaImp();
                        } else {
                            createParaImp(maxcposx - mincposx + 1, maxcposy - mincposy + 1);
                        }
                    }

                    //plot the correlation function and the intensity
                    plotCF(roi1, 3, false);
                    plotIntensityTrace(roi1, 1);
                    if (doMSD) { 					// plot MSD if selected
                        plotMSD(roi1, 3, false);
                    }
                }

            } else {
                JOptionPane.showMessageDialog(null, "Pixel coordinates are out of image.");
            }
        }
    }

    //perform number and brightnes (N&B) analysis
    public void performNB(ImagePlus tmpimp, String $mode) {
        if (!bgrloaded) {
            IJ.showMessage("Background not loaded");
            return;
        }

        // close existing windows
        if (impNWin != null && impNWin.isClosed() == false) {	// close N&B N window
            impNWin.close();
        }
        if (impBWin != null && impBWin.isClosed() == false) {	// close N&B B window
            impBWin.close();
        }

        // close existing windows
        if (impNumWin != null && impNumWin.isClosed() == false) {	// close N&B Num window
            impNumWin.close();
        }
        if (impEpsWin != null && impEpsWin.isClosed() == false) {	// close N&B Eps window
            impEpsWin.close();
        }

        int x = tmpimp.getWidth();
        int y = tmpimp.getHeight();
        int z = lastframe - firstframe + 1; //tmpimp.getStackSize();
        double[] data;
        double[][] mean;
        double[][] mean2;
        double[][] variance;
        double[][] covariance;
        double N;
        double B;

        //determine filter if applicable
        filterArray = new float[width][height];
        if (cbFilter.getSelectedItem() == "none") {
            for (int x1 = 0; x1 < x; x1++) {
                for (int y1 = 0; y1 < y; y1++) {
                    filterArray[x1][y1] = 1;
                }
            }
        }

        if (cbFilter.getSelectedItem() == "Intensity") {
            for (int x1 = 0; x1 < x; x1++) {
                for (int y1 = 0; y1 < y; y1++) {
                    filterArray[x1][y1] = imp.getStack().getProcessor(firstframe).get(x1, y1);
                    if (!(filterArray[x1][y1] > filterLL && filterArray[x1][y1] < filterUL)) {
                        filterArray[x1][y1] = Float.NaN;
                    }
                }
            }
        }

        if (cbFilter.getSelectedItem() == "Mean") {
            for (int x1 = 0; x1 < x; x1++) {
                for (int y1 = 0; y1 < y; y1++) {
                    for (int z1 = firstframe; z1 <= lastframe; z1++) {
                        filterArray[x1][y1] += imp.getStack().getProcessor(z1).get(x1, y1);
                    }
                    filterArray[x1][y1] /= z;
                    if (!(filterArray[x1][y1] > filterLL && filterArray[x1][y1] < filterUL)) {
                        filterArray[x1][y1] = Float.NaN;
                    }
                }
            }
        }

        boolean IsGPUCalculationOK = true;
        if (isgpupresent == 1) {
            roi1StartX = 0;
            roi1StartY = 0;

            // Object to store some of input values for GPU calculations
            GpufitImFCS.ACFParameters GPUparams = new ACFParameters();
            IsGPUCalculationOK = GPU_Initialize_GPUparams(GPUparams, true);

            float[] pixels = new float[GPUparams.w_temp * GPUparams.h_temp * GPUparams.framediff];
            if (IsGPUCalculationOK) {
                IsGPUCalculationOK = GPU_get_pixels(GPUparams, pixels, true);
            }

            IJ.showProgress(10, 100);

            // ********************************************************************************************************************
            // Bleach correction
            // ********************************************************************************************************************   
            double[] bleachcorr_params = new double[GPUparams.w_temp * GPUparams.h_temp * GPUparams.bleachcorr_order];
            if (GPUparams.bleachcorr_gpu && IsGPUCalculationOK) {
                IsGPUCalculationOK = GPU_Calculate_BleachCorrection(GPUparams, pixels, bleachcorr_params);
                IJ.showProgress(30, 100);
            }
            // pixels1 is the output array in which the auto and cross-calculations values are stored.
            double[] pixels1 = new double[GPUparams.width * GPUparams.height * GPUparams.chanum];
            double[] blockvararray = new double[GPUparams.width * GPUparams.height * GPUparams.chanum];
            if (IsGPUCalculationOK) {

                // N&B Calculations. Note cfXdistance and cfYdistance should be zero.
                double[] NBmeanGPU = new double[GPUparams.width * GPUparams.height];
                double[] NBcovarianceGPU = new double[GPUparams.width * GPUparams.height];

                try {
                    IsGPUCalculationOK = GPU_Calculate_ACF(pixels, pixels1, blockvararray, NBmeanGPU, NBcovarianceGPU, bleachcorr_params, GPUparams);

                    if (!IsGPUCalculationOK) {
                        throw new Exception("calculation ACF not OK.");
                    }

                    for (int i = 0; i < GPUparams.width; i++) {
                        for (int j = 0; j < GPUparams.height; j++) {
                            if (filterArray[i][j] == filterArray[i][j]) {
                                NBB[i][j] = (NBcovarianceGPU[j * GPUparams.width + i] - bgrCoVar[i][j]) / NBmeanGPU[j * GPUparams.width + i]; // B = (var -var0) / (mean - offset)
                                NBN[i][j] = NBmeanGPU[j * GPUparams.width + i] * NBmeanGPU[j * GPUparams.width + i] / (NBcovarianceGPU[j * GPUparams.width + i] - bgrCoVar[i][j]); // N= (mean - offset)^2/(var -var0)
                            } else {
                                NBB[i][j] = Float.NaN;
                                NBN[i][j] = Float.NaN;
                            }
                        }
                    }
                    IJ.showProgress(80, 100);

                } catch (Exception e) {
                    e.printStackTrace(System.out);
                }

                IJ.showProgress(100, 100);

                try {
                    GpufitImFCS.resetGPU();
                } catch (Exception e) {
                    System.out.println("Unable to reset GPU");
                }
            }
        }

        if (!IsGPUCalculationOK) {
            System.out.println("An error was encountered while performing calculations on GPU. Calculations will be done on CPU instead.");
        }

        // CPU calculations
        if (isgpupresent != 1 || !IsGPUCalculationOK) {

            //mean and variance of the image
            mean = new double[x][y];
            variance = new double[x][y];
            mean2 = new double[x][y];
            covariance = new double[x][y];
            data = new double[z + 1];

            if (cbNBMode.getSelectedItem().equals("G1")) {
                //calculate mean and covariance of next neighbour frames of the image stack
                for (int i = 0; i < x; i++) {		// calculate Sum(x) and Sum(x^2)
                    for (int j = 0; j < y; j++) {
                        calcIntensityTrace(tmpimp, i, j, i, j, firstframe, lastframe);
                        data = getIntensity(tmpimp, i, j, 1, firstframe, lastframe);
                        for (int k = 1; k < z; k++) {
                            mean[i][j] += data[k];
                            mean2[i][j] += data[k + 1];
                            covariance[i][j] += data[k] * data[k + 1];
                        }
                    }
                }

                for (int i = 0; i < x; i++) {	// caluclate mean and variance: E(x) = Sum(x)/(n-1), E(x*x') = Sum(x*x')/(n-1) and coVar = E(x*x') - E(x)E(x') 
                    for (int j = 0; j < y; j++) {
                        mean[i][j] /= (z - 1);
                        mean2[i][j] /= (z - 1);
                        covariance[i][j] = covariance[i][j] / (z - 1) - mean[i][j] * mean2[i][j];
                    }
                }

                // calculate the N and B images 
                for (int i = 0; i < x; i++) {
                    for (int j = 0; j < y; j++) {
                        if (filterArray[i][j] == filterArray[i][j]) {
                            NBB[i][j] = (covariance[i][j] - bgrCoVar[i][j]) / (mean[i][j]);			// epsilon' = G gamma epsilon = B
                            NBN[i][j] = (mean[i][j]) / NBB[i][j];									// n' = (mean - offset) / epsilon' = (mean - offset) / B
                        } else {
                            NBB[i][j] = Float.NaN;
                            NBN[i][j] = Float.NaN;
                        }
                    }																			// note that offset has already been corrected for in getIntensity()
                }
            } else {
                //calculate mean and variance of the image
                z = tmpimp.getStackSize();
                for (int i = 0; i < x; i++) {		// calculate Sum(x) and Sum(x^2)
                    for (int j = 0; j < y; j++) {
                        calcIntensityTrace(tmpimp, i, j, i, j, 1, z);
                        data = getIntensity(tmpimp, i, j, 1, 1, z);
                        for (int k = 1; k <= z; k++) {
                            mean[i][j] += data[k];
                            variance[i][j] += Math.pow(data[k], 2.0);
                        }
                    }
                }

                for (int i = 0; i < x; i++) {	// calculate mean and variance: E(x) = Sum(x)/n, E(x^2) = Sum(x^2)/n and Var = E(x^2) - E(x)^2 
                    for (int j = 0; j < y; j++) {
                        mean[i][j] /= z;
                        variance[i][j] = variance[i][j] / z - Math.pow(mean[i][j], 2.0);
                    }
                }

                // calculate the N and B values 
                for (int i = 0; i < x; i++) {
                    for (int j = 0; j < y; j++) {
                        if (filterArray[i][j] == filterArray[i][j]) {
                            NBB[i][j] = (variance[i][j] - bgrVar[i][j]) / (mean[i][j]);					// B = (var -var0) / (mean - offset)
                            NBN[i][j] = Math.pow(mean[i][j], 2.0) / (variance[i][j] - bgrVar[i][j]);	// N= (mean - offset)^2/(var -var0)
                        } else {
                            NBB[i][j] = Float.NaN;
                            NBN[i][j] = Float.NaN;
                        }
                    }																				// note that offset has already been corrected for in getIntensity()
                }

                //calculate corrected images using S
                if ($mode.equals("evaluation")) {
                    double S = Double.parseDouble(tfNBS.getText());
                    for (int i = 0; i < x; i++) {
                        for (int j = 0; j < y; j++) {
                            if (filterArray[i][j] == filterArray[i][j]) {
                                NBNum[i][j] = (mean[i][j]) / (NBB[i][j] - S);					// n = (mean - offset) / (B - S)
                                NBEps[i][j] = NBB[i][j] / S - 1;								// epsilon = B/S - 1	
                            } else {
                                NBNum[i][j] = Float.NaN;
                                NBEps[i][j] = Float.NaN;
                            }
                        }																	// note that offset has already been corrected for in getIntensity()
                    }
                }
            }
        }
        NBperformed = true;

        CreateNBimages(x, y, $mode);
    }

    public void CreateNBimages(int x, int y, String $mode) {
        impN = IJ.createImage($impNTitle, "GRAY32", x, y, 1);
        impB = IJ.createImage($impBTitle, "GRAY32", x, y, 1);
        ImageProcessor ipN = impN.getStack().getProcessor(1);
        ImageProcessor ipB = impB.getStack().getProcessor(1);

        for (int i = 0; i < x; i++) {
            for (int j = 0; j < y; j++) {
                ipN.putPixelValue(i, j, NBN[i][j]);
                ipB.putPixelValue(i, j, NBB[i][j]);
            }
        }

        if ($mode.equals("evaluation")) {
            impN.show();
            IJ.run(impN, "Cyan Hot", "");						// apply LUT
            IJ.run(impN, "Enhance Contrast", "saturated=0.35");	// enhance contrast
            impNWin = impN.getWindow();
            impNWin.setLocation(NPosX, NPosY);
            IJ.run(impN, "Set... ", "zoom=" + scimp + " x=" + 0 + " y=" + 0); //then zoom to fit within application
            IJ.run("In [+]", ""); 	// This needs to be used since ImageJ 1.48v to set the window to the right size;

            impB.show();
            IJ.run(impB, "Yellow Hot", "");						// apply LUT
            IJ.run(impB, "Enhance Contrast", "saturated=0.35");	// enhance contrast
            impBWin = impB.getWindow();
            impBWin.setLocation(BPosX, BPosY);
            IJ.run(impB, "Set... ", "zoom=" + scimp + " x=" + 0 + " y=" + 0); //then zoom to fit within application
            IJ.run("In [+]", ""); 	// This needs to be used since ImageJ 1.48v to set the window to the right size;

            if (cbNBMode.getSelectedItem().equals("Calibrated")) {
                impNum = IJ.createImage($impNumTitle, "GRAY32", x, y, 1);
                impEps = IJ.createImage($impEpsTitle, "GRAY32", x, y, 1);
                ImageProcessor ipNum = impNum.getStack().getProcessor(1);
                ImageProcessor ipEps = impEps.getStack().getProcessor(1);

                for (int i = 0; i < x; i++) {
                    for (int j = 0; j < y; j++) {
                        ipNum.putPixelValue(i, j, NBNum[i][j]);
                        ipEps.putPixelValue(i, j, NBEps[i][j]);
                    }
                }

                impNum.show();
                IJ.run(impNum, "Cyan Hot", "");						// apply LUT
                IJ.run(impNum, "Enhance Contrast", "saturated=0.35");	// enhance contrast
                impNumWin = impNum.getWindow();
                impNumWin.setLocation(NumPosX, NumPosY);
                IJ.run(impNum, "Set... ", "zoom=" + scimp + " x=" + 0 + " y=" + 0); //then zoom to fit within application
                IJ.run("In [+]", ""); 	// This needs to be used since ImageJ 1.48v to set the window to the right size;

                impEps.show();
                IJ.run(impEps, "Yellow Hot", "");						// apply LUT
                IJ.run(impEps, "Enhance Contrast", "saturated=0.35");	// enhance contrast
                impEpsWin = impEps.getWindow();
                impEpsWin.setLocation(EpsPosX, EpsPosY);
                IJ.run(impEps, "Set... ", "zoom=" + scimp + " x=" + 0 + " y=" + 0); //then zoom to fit within application
                IJ.run("In [+]", ""); 	// This needs to be used since ImageJ 1.48v to set the window to the right size;
            }
        }
    }

    public void camCalibrate() {
        // open a background image
        String $cd = $imagePath;
        int returnVal;
        ImagePlus tmpimp;
        ImageProcessor tmpip;
        int tmpwidth;
        int tmpheight;
        int tmpframes;
        double[] Bmean;
        double[] BStdDev;
        double[] Nmean;
        double[] NStdDev;

        // check whether the image is already background corrected
        if (!bgrloaded) {
            IJ.showMessage("Background file was not loaded.");
            return;
        }

        // open a set of calibration measurements
        String $suffix;
        String $sfile;

        Bmean = new double[2];
        BStdDev = new double[2];
        Nmean = new double[2];
        NStdDev = new double[2];

        for (int ct = 0; ct < 2; ct++) {
            FileNameExtensionFilter filter = new FileNameExtensionFilter("Tiff files", "tif");
            JFileChooser batchfc = new JFileChooser($cd);
            batchfc.setMultiSelectionEnabled(false);
            batchfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
            batchfc.addChoosableFileFilter(filter);
            returnVal = batchfc.showOpenDialog(null);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = batchfc.getSelectedFile();
                tmpimp = IJ.openImage(file.getAbsolutePath());
                if (tmpimp == null) {
                    JOptionPane.showMessageDialog(null, "Selected file does not exist or it is not an image");
                    return;
                }
                check(tmpimp);
                // tmpimp.show();
                tmpwidth = imp.getWidth();
                tmpheight = imp.getHeight();
                tmpframes = imp.getStackSize();

                if (bgrf != tmpframes || bgrw != tmpwidth || bgrh != tmpheight) {
                    IJ.showMessage("Experiment and background file are not of the same dimensions");
                    tmpimp.changes = false;
                    tmpimp.close();
                    return;
                }
                // correct with the background file and perform N&B
                performNB(tmpimp, "calibration");
                Bmean[ct] = impB.getStatistics().mean;

                tmpimp.changes = false;
                tmpimp.close();
                impB.changes = false;
                impB.close();
                impN.changes = false;
                impN.close();
            }
        }

        double calibratio = Double.parseDouble(tfNBCalibRatio.getText());
        double S = Math.round(100 * (calibratio * Bmean[0] - Bmean[1]) / (calibratio - 1)) / 100;
        tfNBS.setText(Double.toString(S));
        NBslope = S;			// store the slope for saving in output file
        NBcorrected = true;
    }

    // determine minimum value in stack
    public int minDetermination(ImagePlus image) {
        int min;
        min = image.getStack().getProcessor(1).get(1, 1);
        for (int z = 1; z <= frames; z++) {
            for (int x = 1; x < width; x++) {
                for (int y = 1; y < height; y++) {
                    if (image.getStack().getProcessor(z).get(x, y) < min) {
                        min = image.getStack().getProcessor(z).get(x, y);
                    }
                }
            }
        }

        return min;
    }

    public void determinePSF() {
        double maxval = 0;
        double minval = 1;

        // check whether the image is at least 20x20 pixels
        if (width < 20 || height < 20) {
            JOptionPane.showMessageDialog(null, "Image is too small to provide good PSF statistics. At least 20x20 pixels are required.");
            return;
        }

        // get input for estimated PSF range for user
        if (psfDialogue()) {
            psfmaxbin = psfBinEnd - psfBinStart + 1; // calculate for binning as defined by user
            numofpsf = (int) Math.floor((psfEnd - psfStart) / psfStep + 1.1); // how many calculations are necessary?

            IJ.showStatus("Calculating PSF plots");

            if (numofpsf >= 10) {
                JOptionPane.showMessageDialog(null, "PSF range and stepsize lead to too many calculations. Please chose a smaller range.");
                return;
            }

            double[][][][] results;
            psfData = new double[numofpsf][3][psfmaxbin];

            // perform the calculations for all PSFs
            for (int x = 0; x < numofpsf; x++) {
                psfsize = ((psfStart + x * psfStep) * emlambda / NA) / Math.pow(10, 9); // calculate PSF size for each case
                results = dlfit(psfmaxbin, x);
                for (int u = 0; u < psfmaxbin; u++) {
                    psfData[x][0][u] = u + psfBinStart;	// x-value: area size
                    psfData[x][1][u] = results[0][0][1][u];	// y-value: diffusion coefficient
                    psfData[x][2][u] = results[0][0][2][u] / Math.sqrt(width / (u + 1) * height / (u + 1));	// error bars: SEM of diffusion coefficient
                    if (psfData[x][1][u] + psfData[x][2][u] > maxval) {
                        maxval = psfData[x][1][u] + psfData[x][2][u];
                    }
                    if (psfData[x][1][u] - psfData[x][2][u] < minval) {
                        minval = psfData[x][1][u] - psfData[x][2][u];
                    }
                }
            }
            minval /= 2; // make margin for plot label
            plotPSF(minval, maxval, numofpsf);
            IJ.showStatus("Done");
        }
    }

    public void plotPSF(double minval, double maxval, int numofpsf) {
        // Plot PSF calibration and provide fit values in title
        // minval and maxval are determined inthe calling program
        // numofpsf determines how many theoretical values were assumed for the PSF calcualtion 
        double[] labelPosX = {0.1, 0.25, 0.4, 0.55, 0.7, 0.1, 0.25, 0.4, 0.55, 0.7};
        double[] labelPosY = {0.8, 0.8, 0.8, 0.8, 0.8, 0.95, 0.95, 0.95, 0.95, 0.95};
        java.awt.Color[] colors = {java.awt.Color.BLUE, java.awt.Color.CYAN, java.awt.Color.GREEN, java.awt.Color.ORANGE, java.awt.Color.PINK, java.awt.Color.MAGENTA, java.awt.Color.RED, java.awt.Color.LIGHT_GRAY, java.awt.Color.GRAY, java.awt.Color.BLACK};
        Plot plot = new Plot($PSFWindowTitle, "binning", "D [um2/s]");
        plot.setFrameSize(PSFWindowDimX, PSFWindowDimY);
        plot.setLimits(0, psfData[0][0][psfmaxbin - 1] * 1.1, minval * 0.9, maxval * 1.1);
        plot.setJustification(Plot.CENTER);
        for (int x = 0; x < numofpsf; x++) {
            plot.setColor(colors[x]);
            plot.addPoints(psfData[x][0], psfData[x][1], psfData[x][2], Plot.LINE);
            plot.addLabel(labelPosX[x], labelPosY[x], IJ.d2s(psfStart + x * psfStep, decformat));
        }
        plot.draw();

        if (PSFWindow == null || PSFWindow.isClosed() == true) {
            PSFWindow = plot.show();
            PSFWindow.setLocation(PSFWindowPosX, PSFWindowPosY);
        } else {
            PSFWindow.drawPlot(plot);
            PSFWindow.setTitle($PSFWindowTitle);
        }
    }

    // calculate the FCS diffusion law plot for an image
    public void determineDiffusionLaw() {
        maxvalDL = 0;
        minvalDL = 0;

        IJ.showStatus("Calculating FCS Diffusion Law plot");
        difflawbin = Integer.parseInt(tfDLBinEnd.getText());

        if (difflawbin > DiffLawMaxPoint) {
            difflawbin = DiffLawMaxPoint; //upper limit on number of points in diffsion law plot
        }
        double[][][][] results;
        difflaw = new double[3][difflawbin];

        results = dlfit(difflawbin, -1);

        if (tbDLROI.getText() == "All") {
            difflawallbin = difflawbin;							// remember the bins used for DL for the whole image
            for (int u = 1; u <= difflawbin; u++) {			// convert the fitted D into the diffusion law plot								
                difflaw[0][u - 1] = results[0][0][0][u - 1];
                difflaw[1][u - 1] = results[0][0][0][u - 1] / results[0][0][1][u - 1];
                difflaw[2][u - 1] = results[0][0][0][u - 1] / Math.pow(results[0][0][1][u - 1], 2) * Math.sqrt(results[0][0][2][u - 1]);
                if (difflaw[1][u - 1] + difflaw[2][u - 1] > maxvalDL) {
                    maxvalDL = difflaw[1][u - 1] + difflaw[2][u - 1];
                }
                if (difflaw[1][u - 1] - difflaw[2][u - 1] < minvalDL) {
                    minvalDL = difflaw[1][u - 1] - difflaw[2][u - 1];
                }
            }
            plotDiffLaw(difflaw, difflawbin, minvalDL, maxvalDL);
        } else {
            difflawmapbin = difflawbin;									// remember the bins used for DL for the whole image
            diffLawMapwidth = (int) Math.floor(width / Integer.parseInt(tfDLROI.getText()));
            diffLawMapheight = (int) Math.floor(height / Integer.parseInt(tfDLROI.getText()));
            difflawarray = new double[diffLawMapwidth][diffLawMapheight][3][difflawbin];

            for (int i = 0; i < results.length; i++) {
                for (int j = 0; j < results[0].length; j++) {
                    for (int u = 1; u <= difflawbin; u++) {			// convert the fitted D into the diffusion law plot				
                        difflawarray[i][j][0][u - 1] = results[i][j][0][u - 1];
                        difflawarray[i][j][1][u - 1] = results[i][j][0][u - 1] / results[i][j][1][u - 1];
                        difflawarray[i][j][2][u - 1] = results[i][j][0][u - 1] / Math.pow(results[i][j][1][u - 1], 2) * Math.sqrt(results[i][j][2][u - 1]);
                        if (difflawarray[i][j][1][u - 1] + difflawarray[i][j][2][u - 1] > maxvalDL) {
                            maxvalDL = difflawarray[i][j][1][u - 1] + difflawarray[i][j][2][u - 1];
                        }
                        if (difflawarray[i][j][1][u - 1] - difflawarray[i][j][2][u - 1] < minvalDL) {
                            minvalDL = difflawarray[i][j][1][u - 1] - difflawarray[i][j][2][u - 1];
                        }
                    }

                    diffLawFitLim[0] = 1;
                    diffLawFitLim[1] = 5;
                    difflaw = difflawarray[i][j];
                    fitDiffLaw();
                    diffLawFitMap[i][j][0] = difflawfit[0];
                    diffLawFitMap[i][j][1] = difflawfit[1];
                }
            }
            plotDiffLawMaps();
        }
        IJ.showStatus("Done");
    }

    // Fit the Diffusion Law
    public void fitDiffLaw() {
        double[] result;
        int lower = diffLawFitLim[0] - 1;
        int upper = diffLawFitLim[1] - 1;
        int range = upper - lower + 1;

        if (lower >= upper || upper < 0 || lower < 0) {
            IJ.showMessage("Fit start/end not correctly set.");
            diffLawFitLim[0] = 0;
            diffLawFitLim[1] = 0;
            return;
        }

        double[][] difflawsegment = new double[3][range];

        for (int x = 0; x < range; x++) {
            difflawsegment[0][x] = difflaw[0][x + lower];
            difflawsegment[1][x] = difflaw[1][x + lower];
            difflawsegment[2][x] = difflaw[2][x + lower];
        }

        LineFit lfit = new LineFit();
        result = lfit.doFit(difflawsegment[0], difflawsegment[1], difflawsegment[2], range);
        difflawfit[0] = result[0];
        difflawfit[1] = result[1];

        double[][] fitfunc = new double[2][range];
        for (int x = 0; x < range; x++) {
            fitfunc[0][x] = difflaw[0][x + lower];
            fitfunc[1][x] = difflawfit[0] + difflawfit[1] * difflaw[0][x + lower];
        }

        plotDiffLaw(difflaw, difflawbin, minvalDL, maxvalDL);

        Plot plot = difflawWindow.getPlot();
        plot.addLabel(0.3, 0, Math.floor(difflawfit[0] * 100 + 0.5) / 100 + " + " + Math.floor(difflawfit[1] * 100 + 0.5) / 100 + " * Aeff");
        plot.setColor(java.awt.Color.RED);
        plot.addPoints(fitfunc[0], fitfunc[1], Plot.LINE);
        plot.draw();

        difflawWindow = plot.show();
        difflawWindow.setLocation(difflawWindowPosX, difflawWindowPosY);
    }

    // Plot the Diffusion Law
    public void plotDiffLaw(double[][] dldata, int dlbin, double minval, double maxval) {
        //Plot diffusion law and provide fit values in title
        Plot plot = new Plot($difflawWindowTitle, "Aeff [um^2]", "Aeff/D [s]");
        plot.setFrameSize(difflawWindowDimX, difflawWindowDimY);
        plot.setLimits(0, dldata[0][dlbin - 1] * 1.1, minval * 0.9, maxval * 1.1);
        plot.setColor(java.awt.Color.BLUE);
        plot.setJustification(Plot.CENTER);
        plot.addPoints(dldata[0], dldata[1], dldata[2], Plot.CIRCLE);
        plot.draw();

        if (difflawWindow == null || difflawWindow.isClosed() == true) {
            difflawWindow = plot.show();
            plot.setLimitsToFit(true);
            difflawWindow.setLocation(difflawWindowPosX, difflawWindowPosY);
        } else {
            difflawWindow.drawPlot(plot);
            plot.setLimitsToFit(true);
            difflawWindow.setTitle($difflawWindowTitle);
        }
    }

    // Plot Diffusion Law Maps
    public void plotDiffLawMaps() {
        if (impDLMap != null) {		// close diffusion law map window if it exists
            impDLMap.close();
        }

        impDLMap = IJ.createImage($impDLMapTitle, "GRAY32", diffLawMapwidth, diffLawMapheight, 2);	// create a stack for the fit parameters plus chi2, blocking status, filtering mask plus q map
        for (int u = 0; u < 2; u++) {
            ImageProcessor ipDLMap = impDLMap.getStack().getProcessor(u + 1);
            for (int x = 0; x < diffLawMapwidth; x++) {
                for (int y = 0; y < diffLawMapheight; y++) {
                    ipDLMap.putPixelValue(x, y, diffLawFitMap[x][y][u]);
                }
            }
        }

        impDLMap.show();
        impDLMapWin = impDLMap.getWindow();
        impDLMapWin.setLocation(para1PosX, para1PosY);
        impDLMap.setSlice(1);
        IJ.run("Set Label...", "label=" + "tau0");
        impDLMap.setSlice(2);
        IJ.run("Set Label...", "label=" + "slope");

        IJ.run(impDLMap, "Red Hot", "");	// apply "Fire" LUT
        IJ.run(impDLMap, "Original Scale", ""); 	//first set image to original scale
        IJ.run(impDLMap, "Set... ", "zoom=" + scimp + " x=" + (int) Math.floor(diffLawMapwidth / 2) + " y=" + (int) Math.floor(diffLawMapheight / 2)); //then zoom to fit within application
        IJ.run("In [+]", ""); 	// This needs to be used since ImageJ 1.48v to set the window to the right size; 
        // this might be a bug and is an ad hoc solution for the moment; before only the "Set" command was necessary

        impDLMap.setSlice(1);				// set back to slice 1 for viewing
        IJ.run(impDLMap, "Enhance Contrast", "saturated=0.35");	//autoscaling the contrast for slice 1 

        Component[] impDLMapcomp = impDLMapWin.getComponents();	// check which component is the scrollbar and add an AdjustmentListener
        ScrollbarWithLabel impDLMapscrollbar;
        for (int i = 0; i < impDLMapcomp.length; i++) {
            if (impDLMapcomp[i] instanceof ScrollbarWithLabel) {
                impDLMapscrollbar = (ScrollbarWithLabel) impDLMapWin.getComponent(i);
                impDLMapscrollbar.addAdjustmentListener(impDLMapAdjusted);
            }
        }

        // add a mouse listener to retrieve DL plot for each pixel
        impDLMapCan = impDLMap.getCanvas();		// get canvas
        impDLMapCan.setFocusable(true);			// set focusable
        impDLMapCan.addMouseListener(DLMapMouseClicked);
    }

    // Check whether scrollbar was changed
    AdjustmentListener impDLMapAdjusted = (AdjustmentEvent e) -> {
        int slice = impDLMap.getSlice();
        impDLMap.setSlice(slice);
        IJ.run(impDLMap, "Enhance Contrast", "saturated=0.35");
    };

    // plot diffusion law for selected pixel
    MouseListener DLMapMouseClicked = new MouseListener() {
        @Override
        public void mouseClicked(MouseEvent e) {
            int px = e.getX(); // get event coordinates
            int py = e.getY();
            int cpx = (int) Math.floor(impDLMapCan.offScreenX(px)); //convert to pixel numbers
            int cpy = (int) Math.floor(impDLMapCan.offScreenY(py));

            pixeldimx = (pixelsize * 1000 / objmag * Integer.parseInt(tfDLBinEnd.getText())) / Math.pow(10, 9); // set to maximum binning
            pixeldimy = (pixelsize * 1000 / objmag * Integer.parseInt(tfDLBinEnd.getText())) / Math.pow(10, 9); // set to maximum binning
            double tmpx = obsvolFCS_ST2D1p(2) * Math.pow(10, 12);
            pixeldimx = (pixelsize * 1000 / objmag * binningX) / Math.pow(10, 9); // reset to actual value
            pixeldimy = (pixelsize * 1000 / objmag * binningY) / Math.pow(10, 9); // reset to actual value

            double[][] fitfunc = new double[2][2];

            fitfunc[0][0] = 0;
            fitfunc[1][0] = diffLawFitMap[cpx][cpy][0];
            fitfunc[0][1] = tmpx;
            fitfunc[1][1] = diffLawFitMap[cpx][cpy][0] + diffLawFitMap[cpx][cpy][1] * tmpx;

            for (double[][][] difflawarray1 : difflawarray) {
                for (int j = 0; j < difflawarray[0].length; j++) {
                    for (int u = 1; u <= difflawarray[0][0][0].length; u++) {
                        // convert the fitted D into the diffusion law plot
                        if (difflawarray1[j][1][u - 1] + difflawarray1[j][2][u - 1] > maxvalDL) {
                            maxvalDL = difflawarray1[j][1][u - 1] + difflawarray1[j][2][u - 1];
                        }
                        if (difflawarray1[j][1][u - 1] - difflawarray1[j][2][u - 1] < minvalDL) {
                            minvalDL = difflawarray1[j][1][u - 1] - difflawarray1[j][2][u - 1];
                        }
                    }
                }
            }

            plotDiffLaw(difflawarray[cpx][cpy], difflawmapbin, minvalDL, maxvalDL);

            Plot plot = difflawWindow.getPlot();
            plot.addLabel(0.3, 0, Math.floor(diffLawFitMap[cpx][cpy][0] * 100 + 0.5) / 100 + " + " + Math.floor(diffLawFitMap[cpx][cpy][1] * 100 + 0.5) / 100 + " * Aeff");
            plot.setColor(java.awt.Color.RED);
            plot.addPoints(fitfunc[0], fitfunc[1], Plot.LINE);
            plot.draw();
        }

        @Override
        public void mousePressed(MouseEvent e) {
        }	// other mouse events have no action associated yet

        @Override
        public void mouseReleased(MouseEvent e) {
        }

        @Override
        public void mouseEntered(MouseEvent e) {
        }

        @Override
        public void mouseExited(MouseEvent e) {
        }
    };

    // Calculate correlations for a ROI
    public void correlateROI(Roi improi) {

        IJ.showStatus("Correlating ROI");

        boolean IsGPUCalculationOK = true;

        if (isgpupresent == 1) {
            // The user can only perform 2 particle fits in the GPU version.
            prepareFit();
            IsGPUCalculationOK = GPU_Calculate_ACF_All(improi);
            if (IsGPUCalculationOK) {
                IJ.showStatus("Done");
            }
        }

        if (!IsGPUCalculationOK) {
            System.out.println("An error was encountered while performing calculations on GPU. Calculations will be done on CPU instead.");
        }

        // CPU calculations
        if (isgpupresent != 1 || !IsGPUCalculationOK) {
            Rectangle imprect = improi.getBounds();
            int startXmap = (int) Math.ceil(imprect.getX() / pixbinX);
            int startYmap = (int) Math.ceil(imprect.getY() / pixbinY);
            int startX = startXmap * pixbinX;
            int startY = startYmap * pixbinY;
            int endXmap = (int) Math.floor((imprect.getX() + imprect.getWidth() - binningX) / pixbinX);
            int endYmap = (int) Math.floor((imprect.getY() + imprect.getHeight() - binningY) / pixbinY);
            int endX = endXmap * pixbinX;
            int endY = endYmap * pixbinY;
            int pixcount = 0;
            double q1;
            double q2;

            //switch the FCCS display off; this would result in too many functions to be displayed
            tbFCCSDisplay.setSelected(false);
            tbFCCSDisplay.setText("Off");
            filterArray = new float[width][height]; // calculate the mean image of the stack

            if (cbFilter.getSelectedItem() == "Mean" || cbFilter.getSelectedItem() == "Intensity") {
                initializeFitres(3, width, height, noparam); // reset the fitresult array and thus the parameter window
                if (cbFilter.getSelectedItem() == "Mean") {
                    for (int x1 = startX; x1 <= endX; x1 = x1 + pixbinX) {
                        for (int x2 = startY; x2 <= endY; x2 = x2 + pixbinY) {
                            for (int x3 = firstframe; x3 <= lastframe; x3++) {
                                for (int x4 = 0; x4 < binningX; x4++) {
                                    for (int x5 = 0; x5 < binningY; x5++) {
                                        if (improi.contains(x1, x2) && improi.contains(x1, x2 + binningY - 1) && improi.contains(x1 + binningX - 1, x2) && improi.contains(x1 + binningX - 1, x2 + binningY - 1)) {
                                            filterArray[x1][x2] += imp.getStack().getProcessor(x3).get(x1 + x4, x2 + x5);
                                            pixcount++;
                                        } else {
                                            filterArray[x1][x2] = Float.NaN;
                                        }
                                    }
                                }
                            }
                            filterArray[x1][x2] /= (lastframe - firstframe + 1);
                        }
                    }
                } else { // if "Intensity" was selcted, then get the first frame and bin if necessary
                    for (int x1 = startX; x1 <= endX; x1 = x1 + pixbinX) {
                        for (int x2 = startY; x2 <= endY; x2 = x2 + pixbinY) {
                            for (int x3 = 0; x3 < binningX; x3++) {
                                for (int x4 = 0; x4 < binningY; x4++) {
                                    if (improi.contains(x1, x2) && improi.contains(x1, x2 + binningY - 1) && improi.contains(x1 + binningX - 1, x2) && improi.contains(x1 + binningX - 1, x2 + binningY - 1)) {
                                        filterArray[x1][x2] += imp.getStack().getProcessor(firstframe).get(x1 + x3, x2 + x4);
                                        pixcount++;
                                    } else {
                                        filterArray[x1][x2] = Float.NaN;
                                    }
                                }
                            }
                        }
                    }
                }
            } else { // if "none" was selected
                for (int x1 = startX; x1 <= endX; x1 = x1 + pixbinX) {
                    for (int x2 = startY; x2 <= endY; x2 = x2 + pixbinY) {
                        for (int x3 = 0; x3 < binningX; x3++) {
                            for (int x4 = 0; x4 < binningY; x4++) {
                                if (improi.contains(x1, x2) && improi.contains(x1, x2 + binningY - 1) && improi.contains(x1 + binningX - 1, x2) && improi.contains(x1 + binningX - 1, x2 + binningY - 1)) {
                                    filterArray[x1][x2] += imp.getStack().getProcessor(firstframe).get(x1 + x3, x2 + x4);
                                    pixcount++;
                                } else {
                                    filterArray[x1][x2] = Float.NaN;
                                }
                            }
                        }
                    }
                }
            }

            //do the FCS or DC-FCCS evaluation
            if (cbFitModel.getSelectedItem() == "FCS") {
                if (doFit) {
                    prepareFit();
                    for (int x = startXmap; x <= endXmap; x++) {
                        for (int y = startYmap; y <= endYmap; y++) {
                            if (filterArray[x * pixbinX][y * pixbinY] >= filterLL * binningX * binningY && filterArray[x * pixbinX][y * pixbinY] <= filterUL * binningX * binningY) {
                                calcIntensityTrace(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, firstframe, lastframe);
                                correlate(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, 0, firstframe, lastframe);
                                if (!fit(x, y, 0, "FCS")) {
                                    return;
                                }
                            }
                        }
                        IJ.showProgress(x - startXmap, startXmap - endXmap);
                    }
                } else {
                    for (int x = startXmap; x <= endXmap; x++) {
                        for (int y = startYmap; y <= endYmap; y++) {
                            if (!Float.isNaN(filterArray[x * pixbinX][y * pixbinY])) {
                                calcIntensityTrace(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, firstframe, lastframe);
                                correlate(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, 0, firstframe, lastframe);
                            }
                            IJ.showProgress(x - startXmap, startXmap - endXmap);
                        }
                    }
                }
            }

            if (cbFitModel.getSelectedItem() == "DC-FCCS") {
                if (doFit) {
                    prepareFit();
                    for (int x = startXmap; x <= endXmap; x++) {
                        for (int y = startYmap; y <= endYmap; y++) {
                            if (filterArray[x * pixbinX][y * pixbinY] >= filterLL * binningX * binningY && filterArray[x * pixbinX][y * pixbinY] <= filterUL * binningX * binningY) {
                                calcIntensityTrace(imp, x * pixbinX, y * pixbinY, x * pixbinX, y * pixbinY, firstframe, lastframe);
                                correlate(imp, x * pixbinX, y * pixbinY, x * pixbinX, y * pixbinY, 0, firstframe, lastframe);
                                if (!fit(x, y, 0, "FCS")) {
                                    return;
                                }
                                calcIntensityTrace(imp, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, firstframe, lastframe);
                                correlate(imp, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, 1, firstframe, lastframe);
                                if (!fit(x, y, 1, "FCS")) {
                                    return;
                                }
                                calcIntensityTrace(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, firstframe, lastframe);
                                correlate(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, 2, firstframe, lastframe);
                                if (!fit(x, y, 2, "DC-FCCS")) {
                                    return;
                                }
                                //calculate q value
                                q1 = fitres[0][x][y][0] / fitres[2][x][y][0];
                                q2 = fitres[1][x][y][0] / fitres[2][x][y][0];
                                if (q1 > q2) {
                                    CCFq[x][y] = q1;
                                } else {
                                    CCFq[x][y] = q2;
                                }
                                if (CCFq[x][y] > 1) {
                                    CCFq[x][y] = Double.NaN;	//cross-correlation > 1 is non-physical, most likely due to very noisy CFs, better discarded
                                }
                            }
                        }
                        IJ.showProgress(x - startXmap, startXmap - endXmap);
                    }
                } else {
                    for (int x = startXmap; x <= endXmap; x++) {
                        for (int y = startYmap; y <= endYmap; y++) {
                            if (filterArray[x * pixbinX][y * pixbinY] >= 0) {
                                calcIntensityTrace(imp, x * pixbinX, y * pixbinY, x * pixbinX, y * pixbinY, firstframe, lastframe);
                                correlate(imp, x * pixbinX, y * pixbinY, x * pixbinX, y * pixbinY, 0, firstframe, lastframe);
                                calcIntensityTrace(imp, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, firstframe, lastframe);
                                correlate(imp, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, 1, firstframe, lastframe);
                                calcIntensityTrace(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, firstframe, lastframe);
                                correlate(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance, 2, firstframe, lastframe);
                            }
                        }
                        IJ.showProgress(x - startXmap, startXmap - endXmap);
                    }
                }
            }

            if (pixcount > 0) {
                IJ.showStatus("Plotting data.");
                plotCF(improi, 2, false);
                calcAverageIntensityTrace(filterArray, startX, startY, endX, endY, firstframe, lastframe);
                plotIntensityTrace(improi, 2);
                if (doMSD) { 					// plot MSD if selected
                    plotMSD(improi, 2, false);
                }
                if (doFit) {	// create parameter map window
                    createParaImp(maxcposx - mincposx + 1, maxcposy - mincposy + 1);
                }
            } else {
                JOptionPane.showMessageDialog(null, "ROI does not cover a whole single pixel in the binned image.");
            }

            IJ.showStatus("Done");
        }
    }

    // batch processing
    public void batchProcessing() {
        int numOfFiles;
        String $suffix;
        String $sfile;
        askOnRewrite = false;
        JFileChooser batchfc = new JFileChooser($workingDir);
        batchfc.setMultiSelectionEnabled(true);
        batchfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
        int returnVal = batchfc.showOpenDialog(null);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File[] files = batchfc.getSelectedFiles();
            numOfFiles = files.length;
            for (File currentFile : files) {
                imp = IJ.openImage(currentFile.getAbsolutePath());
                imp.show();
                obtainImage();
                setParameters();

                // perform selected functions on the current image
                if (batchCorrelateAll) {
                    if (batchFit) {
                        doFit = true;
                    }
                    if (imp.getOverlay() != null) {
                        imp.getOverlay().clear();
                        imp.setOverlay(imp.getOverlay());
                    }

                    int roi2StartX;
                    int roi2WidthX;
                    int roi2StartY;
                    int roi2HeightY;
                    if (cfXDistance > 0) {
                        roi1StartX = 0;
                        roi2StartX = cfXDistance;
                    } else {
                        roi1StartX = -cfXDistance;
                        roi2StartX = 0;
                    }
                    if (cfYDistance > 0) {
                        roi1StartY = 0;
                        roi2StartY = cfYDistance;
                    } else {
                        roi1StartY = -cfYDistance;
                        roi2StartY = 0;
                    }
                    if (overlap) {
                        roi1WidthX = width - Math.abs(cfXDistance);
                        roi1HeightY = height - Math.abs(cfYDistance);
                    } else {
                        roi1WidthX = (int) Math.floor((width - Math.abs(cfXDistance)) / binningX) * binningX;
                        roi1HeightY = (int) Math.floor((height - Math.abs(cfYDistance)) / binningY) * binningY;
                    }
                    roi2WidthX = roi1WidthX;
                    roi2HeightY = roi1HeightY;
                    Roi impRoi1 = new Roi(roi1StartX, roi1StartY, roi1WidthX, roi1HeightY);
                    impRoi1.setStrokeColor(java.awt.Color.GREEN);
                    imp.setRoi(impRoi1);
                    Roi impRoi2 = new Roi(roi2StartX, roi2StartY, roi2WidthX, roi2HeightY);
                    if (cfXDistance != 0 || cfYDistance != 0) {
                        impRoi2.setStrokeColor(java.awt.Color.RED);
                        Overlay cfov = new Overlay(impRoi2);
                        imp.setOverlay(cfov);
                    }
                    correlateRoiWorker correlateRoiInstant = new correlateRoiWorker(impRoi1);
                    correlateRoiInstant.execute();
                    while (correlateRoiInstant.isDone() == false) {
                    }
                }

                if (batchPSF) {
                    correlatePSFWorker correlatePSFInstant = new correlatePSFWorker();
                    correlatePSFInstant.execute();
                    while (correlatePSFInstant.isDone() == false) {
                    }
                }

                if (batchDiffLaw) {
                    correlateDiffLawWorker correlateDiffLawInstant = new correlateDiffLawWorker();
                    correlateDiffLawInstant.execute();
                    while (correlateDiffLawInstant.isDone() == false) {
                    }
                    diffLawFitLim[0] = 1;
                    diffLawFitLim[1] = 5;
                    fitDiffLaw();
                }

                if (batchHorDCCF) {
                    String $dccf = (String) cbDCCF.getSelectedItem();
                    int mode;
                    int wx = pixelWidthX;
                    int hy = pixelHeightY;
                    mode = 2;
                    correlateDCCFWorker correlateDCCFInstant = new correlateDCCFWorker(mode, wx, hy);
                    correlateDCCFInstant.execute();
                    while (correlateDCCFInstant.isDone() == false) {
                    }
                }

                if (batchVerDCCF) {
                    String $dccf = (String) cbDCCF.getSelectedItem();
                    int mode;
                    int wx = pixelWidthX;
                    int hy = pixelHeightY;
                    mode = 1;
                    correlateDCCFWorker correlateDCCFInstant = new correlateDCCFWorker(mode, wx, hy);
                    correlateDCCFInstant.execute();
                    while (correlateDCCFInstant.isDone() == false) {
                    }
                }

                if (batchDiaUpDCCF) {
                    String $dccf = (String) cbDCCF.getSelectedItem();
                    int mode;
                    int wx = pixelWidthX;
                    int hy = pixelHeightY;
                    mode = 3;
                    correlateDCCFWorker correlateDCCFInstant = new correlateDCCFWorker(mode, wx, hy);
                    correlateDCCFInstant.execute();
                    while (correlateDCCFInstant.isDone() == false) {
                    }
                }

                if (batchDiaDownDCCF) {
                    String $dccf = (String) cbDCCF.getSelectedItem();
                    int mode;
                    int wx = pixelWidthX;
                    int hy = pixelHeightY;
                    mode = 4;
                    correlateDCCFWorker correlateDCCFInstant = new correlateDCCFWorker(mode, wx, hy);
                    correlateDCCFInstant.execute();
                    while (correlateDCCFInstant.isDone() == false) {
                    }
                }

                // clean up and close all windows
                if (acfWindow != null && acfWindow.isClosed() == false) {	// close ACF window
                    acfWindow.close();
                }
                if (intWindow != null && intWindow.isClosed() == false) {	// close intensity trace window
                    intWindow.close();
                }
                if (resWindow != null && resWindow.isClosed() == false) {	// close fit residuals window
                    resWindow.close();
                }
                if (sdWindow != null && sdWindow.isClosed() == false) {	// close SD window
                    sdWindow.close();
                }
                if (msdWindow != null && msdWindow.isClosed() == false) {	// close SD window
                    msdWindow.close();
                }
                if (difflawWindow != null && difflawWindow.isClosed() == false) {	// close diffusion law window
                    difflawWindow.close();
                }
                if (PSFWindow != null && PSFWindow.isClosed() == false) {	// close PSF window
                    PSFWindow.close();
                }
                if (paraCorWindow != null && paraCorWindow.isClosed() == false) {	// close paraCor window
                    paraCorWindow.close();
                }
                if (blockingWindow != null && blockingWindow.isClosed() == false) {	// close blocking window
                    blockingWindow.close();
                }
                if (impCovWin != null && impCovWin.isClosed() == false) {	// close covariance window
                    impCovWin.close();
                }
                if (impNWin != null && impNWin.isClosed() == false) {	// close covariance window
                    impNWin.close();
                }
                if (impBWin != null && impBWin.isClosed() == false) {	// close covariance window
                    impBWin.close();
                }
                for (int x = 0; x < dccfMax; x++) {
                    if (impDCCFWin[x] != null && impDCCFWin[x].isClosed() == false) {	// close DCCF windows
                        impDCCFWin[x].close();
                    }
                }
                for (int x = 0; x < dccfMax; x++) {
                    if (histDCCFWin[x] != null && histDCCFWin[x].isClosed() == false) {	// close dCCF historgram windows
                        histDCCFWin[x].close();
                    }
                }
                if (doFit) {	// close parameter map window
                    impPara1.close();
                }
                // close histogram windows
                if (histWin != null && histWin.isClosed() == false) {
                    histWin.close();
                }

                imp.changes = false;
                imp.close();

                // save the results
                if ($batchSuffix.equals("")) {
                    $suffix = new SimpleDateFormat("yyyy_MM_dd-HH_mm_ss").format(new Date());
                } else {
                    $suffix = $batchSuffix;
                }

                $sfile = currentFile.toString();
                int dotind = $sfile.lastIndexOf('.');
                if (dotind != -1) {
                    $sfile = $sfile.substring(0, dotind);
                }

                File newFile = new File($sfile + $suffix + ".xlsx");

                if (newFile.exists()) {
                    JOptionPane.showMessageDialog(null, "A file with name " + newFile + " already exists.");
                } else {
                    writeExperiment(newFile, "Failed to write data files", false);
                }
                IJ.log(currentFile.toString() + " processed. Data saved in " + newFile.toString());
            }
            IJ.showMessage("Batch processing done.");
        } else {
            setImp = false;
        }
    }

    // prepare the scatter plot for two parameters as selected from the FCS panel menue
    public void plotScatterPlot() {
        String $paracor = (String) cbParaCor.getSelectedItem();
        String $xlabel = "";
        String $ylabel = "";
        int x, y;
        double minx = 0;
        double maxx = 0;
        double miny = 0;
        double maxy = 0;
        double[][] scplot = new double[2][(maxcposx - mincposx + 1) * (maxcposy - mincposy + 1)];

        if (null != $paracor) {
            switch ($paracor) {
                case "N vs F2":
                    for (x = mincposx; x <= maxcposx; x++) {
                        for (y = mincposy; y <= maxcposy; y++) {
                            scplot[0][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][0];
                            scplot[1][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][5];
                            $xlabel = "F2";
                            $ylabel = "N";
                        }
                    }
                    break;
                case "D vs F2":
                    for (x = mincposx; x <= maxcposx; x++) {
                        for (y = mincposy; y <= maxcposy; y++) {
                            scplot[0][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][1] * Math.pow(10, 12);
                            scplot[1][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][5];
                            $xlabel = "F2";
                            $ylabel = "D";
                        }
                    }
                    break;
                case "N*(1-F2) vs D":
                    for (x = mincposx; x <= maxcposx; x++) {
                        for (y = mincposy; y <= maxcposy; y++) {
                            scplot[0][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][0] * (1 - fitres[0][x][y][5]);
                            scplot[1][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][1] * Math.pow(10, 12);
                            $xlabel = "D";
                            $ylabel = "N*(1-F2)";
                        }
                    }
                    break;
                case "N*F2 vs D2":
                    for (x = mincposx; x <= maxcposx; x++) {
                        for (y = mincposy; y <= maxcposy; y++) {
                            scplot[0][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][0] * fitres[0][x][y][5];
                            scplot[1][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][6] * Math.pow(10, 12);
                            $xlabel = "D2";
                            $ylabel = "N*F2";
                        }
                    }
                    break;
                case "D vs Sqrt(vx^2+vy^2)":
                    for (x = mincposx; x <= maxcposx; x++) {
                        for (y = mincposy; y <= maxcposy; y++) {
                            scplot[0][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][1] * Math.pow(10, 12);
                            scplot[1][x * (maxcposy - mincposy + 1) + y] = Math.sqrt(Math.pow(fitres[0][x][y][2], 2) + Math.pow(fitres[0][x][y][3], 2)) * Math.pow(10, 6);
                            $xlabel = "Sqrt(vx^2+vy^2)";
                            $ylabel = "D";
                        }
                    }
                    break;
                case "D2 vs Sqrt(vx^2+vy^2)":
                    for (x = mincposx; x <= maxcposx; x++) {
                        for (y = mincposy; y <= maxcposy; y++) {
                            scplot[0][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][6] * Math.pow(10, 12);
                            scplot[1][x * (maxcposy - mincposy + 1) + y] = Math.sqrt(Math.pow(fitres[0][x][y][2], 2) + Math.pow(fitres[0][x][y][3], 2)) * Math.pow(10, 6);
                            $xlabel = "Sqrt(vx^2+vy^2)";
                            $ylabel = "D2";
                        }
                    }
                    break;
                default:
                    // "N vs D"
                    for (x = mincposx; x <= maxcposx; x++) {
                        for (y = mincposy; y <= maxcposy; y++) {
                            scplot[0][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][0];
                            scplot[1][x * (maxcposy - mincposy + 1) + y] = fitres[0][x][y][1] * Math.pow(10, 12);
                            $xlabel = "D";
                            $ylabel = "N";
                        }
                    }
                    break;
            }
        }

        // find min and max values for x and y axes
        for (x = mincposx; x <= maxcposx; x++) {
            for (y = mincposy; y <= maxcposy; y++) {
                if (minx > scplot[1][x * (maxcposy - mincposy + 1) + y]) {
                    minx = scplot[1][x * (maxcposy - mincposy + 1) + y];
                }
                if (maxx < scplot[1][x * (maxcposy - mincposy + 1) + y]) {
                    maxx = scplot[1][x * (maxcposy - mincposy + 1) + y];
                }
                if (miny > scplot[0][x * (maxcposy - mincposy + 1) + y]) {
                    miny = scplot[0][x * (maxcposy - mincposy + 1) + y];
                }
                if (maxy < scplot[0][x * (maxcposy - mincposy + 1) + y]) {
                    maxy = scplot[0][x * (maxcposy - mincposy + 1) + y];
                }
            }
        }

        // plot
        Plot plot = new Plot($paraCorWindowTitle, $xlabel, $ylabel, empty, empty);
        plot.setFrameSize(paraCorWindowDimX, paraCorWindowDimY);
        plot.setLimits(minx, maxx, miny, maxy);
        plot.setColor(java.awt.Color.BLUE);
        plot.setJustification(Plot.CENTER);
        plot.draw();

        // create plot label 
        plot.addLabel(0.5, 0, $paracor);

        // either create a new plot window or plot within the existing window
        if (paraCorWindow == null || paraCorWindow.isClosed() == true) {
            paraCorWindow = plot.show();
            paraCorWindow.setLocation(paraCorWindowPosX, paraCorWindowPosY);
        } else {
            paraCorWindow.drawPlot(plot);
            paraCorWindow.setTitle($paraCorWindowTitle);
        }

        plot.addPoints(scplot[1], scplot[0], Plot.CIRCLE);
        paraCorWindow.drawPlot(plot);
    }

    // calculate dCCF image
    public void dccf(int mode, int wx, int hy) {
        // mode determines whether dccf is calculated horizontally, vertically, diagonally up or diagonally down
        // wx, hy stand for width and height of the area to be used.
        int sx; // start and end values for x and y
        int ex;
        int sy;
        int ey;
        int dx;	// determines the direction for the dccf, as determined by mode
        int dy;
        double[][][][] acfMem = new double[2][width][height][chanum];	// an array to remember the values in acf[0][x][y][c], fited acf, residuals and standard deviation and return them there after the calculation is finished
        String $imagetitle = $dccfTitle + " - " + (String) cbDCCF.getSelectedItem();
        String $histtitle = "Hist - " + $dccfTitle + " - " + (String) cbDCCF.getSelectedItem();
        IJ.showStatus("Correlating all pixels");

        // crating a copy of acf[0][x][y][c] and related arrays
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int c = 0; c < chanum; c++) {
                    acfMem[0][x][y][c] = acf[0][x][y][c];
                    acfMem[1][x][y][c] = sdacf[0][x][y][c];
                }
            }
        }

        //set beginning and end points
        if (mode == 1) {
            sx = 0;
            ex = wx - 1;
            sy = 0;
            ey = hy;
            dx = 1;
            dy = 0;
        } else if (mode == 2) {
            sx = 0;
            ex = wx;
            sy = 0;
            ey = hy - 1;
            dx = 0;
            dy = 1;
        } else if (mode == 3) {
            sx = 0;
            ex = wx - 1;
            sy = 1;
            ey = hy;
            dx = 1;
            dy = -1;
        } else {
            sx = 0;
            ex = wx - 1;
            sy = 0;
            ey = hy - 1;
            dx = 1;
            dy = 1;
        }

        // width and height of the dCCF image 
        int dw = ex - sx + 1;
        int dh = ey - sy + 1;

        for (int x = sx; x <= ex; x++) {
            for (int y = sy; y <= ey; y++) {
                dccf[mode - 1][x - sx][y - sy] = 0.0;		//initialize given pixel with 0
                calcIntensityTrace(imp, x * pixbinX, y * pixbinY, (x + dx) * pixbinX, (y + dy) * pixbinY, firstframe, lastframe);
                correlate(imp, x * pixbinX, y * pixbinY, (x + dx) * pixbinX, (y + dy) * pixbinY, 0, firstframe, lastframe);
                calcIntensityTrace(imp, (x + dx) * pixbinX, (y + dy) * pixbinY, x * pixbinX, y * pixbinY, firstframe, lastframe);
                correlate(imp, (x + dx) * pixbinX, (y + dy) * pixbinY, x * pixbinX, y * pixbinY, 0, firstframe, lastframe);
                for (int i = 1; i < chanum; i++) {
                    dccf[mode - 1][x - sx][y - sy] += acf[0][x][y][i] - acf[0][x + dx][y + dy][i];
                }
            }
            IJ.showProgress(x - sx, ex - sx + 1);
        }

        createDCCFWindow(dw, dh, $imagetitle, $histtitle, mode - 1);
        dccfCalculated[mode - 1] = true;

        // reverting changes to acf[0][x][y][c] and related arrays done by dlfit
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int c = 0; c < chanum; c++) {
                    acf[0][x][y][c] = acfMem[0][x][y][c];
                    sdacf[0][x][y][c] = acfMem[1][x][y][c];
                }
            }
        }

        IJ.showStatus("Done");
    }

    public void createDCCFWindow(int dw, int dh, String $imagetitle, String $histtitle, int mode) {
        // dw, dh are width and height for the window to be plotted
        // the strings are the title for the image and histograms
        if (impDCCFWin[mode] != null && impDCCFWin[mode].isClosed() == false) {	// close DCCF window for the given direction if it exists
            impDCCFWin[mode].close();
        }
        if (histDCCFWin[mode] != null && histDCCFWin[mode].isClosed() == false) {	// close dCCF historgram window for the given direction if it exists
            histDCCFWin[mode].close();
        }

        impDCCF = IJ.createImage($imagetitle, "GRAY32", dw, dh, 1);	// create ImagePlus
        ImageProcessor ipDCCF = impDCCF.getStack().getProcessor(1);
        for (int x = 0; x < dw; x++) {
            for (int y = 0; y < dh; y++) {
                ipDCCF.putPixelValue(x, y, dccf[mode][x][y]);
            }
        }

        impDCCF.show();
        impDCCFWin[mode] = impDCCF.getWindow();
        impDCCFWin[mode].setLocation(DCCFPosX, DCCFPosY);
        IJ.run(impDCCF, "Original Scale", ""); 	//first set image to original scale
        IJ.run(impDCCF, "Set... ", "zoom=" + scimp + " x=" + (int) Math.floor(dw / 2) + " y=" + (int) Math.floor(dh / 2)); //then zoom to fit within application
        IJ.run("In [+]", ""); 	// This needs to be used since ImageJ 1.48v to set the window to the right size; 
        // this might be a bug and is an ad hoc solution for the moment; before only the "Set" command was necessary
        IJ.run(impDCCF, "Enhance Contrast", "saturated=0.35");	//autoscale the contrast
        double histMin = impDCCF.getStatistics().min;
        double histMax = impDCCF.getStatistics().max;
        int histYMax = impDCCF.getStatistics().histYMax;
        int pixelCount = impDCCF.getStatistics().pixelCount;
        double stdDev = impDCCF.getStatistics().stdDev;
        int nBins = (int) Math.floor(Math.cbrt(pixelCount) * (histMax - histMin) / (4 * stdDev)) + 1;
        $histDCCFWinTitle[mode] = $histtitle;
        histDCCFWin[mode] = new HistogramWindow($histtitle, impDCCF, nBins, histMin, histMax, histYMax);
        //histDCCFBin[dccfCount] = binning;
        histDCCFWin[mode].setLocationAndSize(histDCCFPosX, histDCCFPosY, histDCCFDimX, histDCCFDimY);
    }

    /* 
 * The main programs start here. This includes correlations, intensity calculations, bleach correction and data fitting.
 * 
 * public void correlate(ImagePlus image, int px1, int py1, int px2, int py2, int kcf, int initialframe, int finalframe): correlate any two pixels with coordinates (px1, py1) and (px2, py2)
 * public Map correlator(double[][] intcor, int initialframe, int finalframe): first does blocking and then starts the correlation
 * public int blockTransform(double[][] intcor, int num, int blocklag): perform blocking on the data to see whether a good estimate of the SD can be obtained
 * public Map calculateCF(double[][] intcor, int num, int ind, int blocklag): performs the actual calculation of the correlation, SD, and covariance matrix
 * public double[] getIntensity(ImagePlus image, int px, int py, int mode, int initialframe, int finalframe): calculates the intensity for the correlations and performs a bleach correction if selected
 * public void calcIntensityTrace(ImagePlus image, int ipx1, int ipy1, int ipx2, int ipy2, int initialframe, int finalframe): calculates the intensity trace for the intensity plot; 
 *    this has a reduced number of pixels compared to the full intensity trace
 * public void calcAverageIntensityTrace(Roi introi, int initialframe, int finalframe): calcualtes the average intensity for the whole stack
 * 
     */
    //correlate one pixel with itself or two pixels with each other
    public void correlate(ImagePlus image, int px1, int py1, int px2, int py2, int kcf, int initialframe, int finalframe) {
        // image: the imp to be used
        // px1, py1, px2, py2: pixel cooredinates for pixel 1 and pixel 2 which are to be correalted
        // if px1 = px2 AND py1 = py2: then a autocorrelation is calculated
        // kcf (0, 1, or 2) determines whether ACF1, ACF2, or CCF is calculated
        // initialframe and finalframe provide the range of frames to be used for the correlation
        int num = (finalframe - initialframe + 1); 	// total number of frames to be correlated
        int numofsw; 							// number of sliding windows
        int swinitialframe; 					// if sliding window (bleach) correction is selected these are the initial and final frames of the sub-windows
        int swfinalframe;
        int pxm1;								// pixel coordinates on the binned grid used to store the output and map it to the parameter map
        int pym1;

        if (kcf == 1 && px1 == px2 && py1 == py2) {	// if red channel in the DC-FCCS mode, map the output to the corresponding green-channel pixels on the binned grid
            pxm1 = (int) (px1 - cfXDistance) / pixbinX;
            pym1 = (int) (py1 - cfYDistance) / pixbinY;
        } else {										// otherwise map to the pixel on the pixel on the binned grid
            pxm1 = (int) px1 / pixbinX;
            pym1 = (int) py1 / pixbinY;
        }

        String $bcmode = (String) cbBleachCor.getSelectedItem();
        if ("Sliding Window".equals($bcmode)) {						// if sliding window correction was selected then ...
            numofsw = (int) Math.floor(num / slidingWindowLength);
            datac = new double[2][slidingWindowLength + 1];			// get the intensity data for the correlation
            for (int x = 0; x < chanum; x++) {				// initialize acf and sdacf
                acf[kcf][pxm1][pym1][x] = 0;
                varacf[kcf][pxm1][pym1][x] = 0;
                sdacf[kcf][pxm1][pym1][x] = 0;
            }
            for (int i = 0; i < numofsw; i++) {						// loop over the number of sliding windows to calculate the correlations of the sub-intervals of the intensity trace
                swinitialframe = i * slidingWindowLength + initialframe;
                swfinalframe = (i + 1) * slidingWindowLength + initialframe - 1;
                datac[0] = getIntensity(image, px1, py1, 1, swinitialframe, swfinalframe);		// getIntensity for first pixel; performs a bleach correction if indicated in the panel
                if (px1 != px2 || py1 != py2) {												// if the two pixels are not equal (i.e. for a cross-correaltion)
                    datac[1] = getIntensity(image, px2, py2, 2, swinitialframe, swfinalframe);	// getIntensity for second pixel
                } else {
                    datac[1] = datac[0];			// otherwise perform an autocorrelation
                }

                Map result;
                result = correlator(datac, swinitialframe, swfinalframe);		// correlate the data; no normalization yet

                for (int x = 0; x <= (chanum - 1); x++) {	// calculate the normalized acf and its standard deviation
                    acf[kcf][pxm1][pym1][x] += ((double[]) result.get("corav"))[x] / numofsw;
                    varacf[kcf][pxm1][pym1][x] += ((double[]) result.get("blockvar"))[x] / numofsw;
                    sdacf[kcf][pxm1][pym1][x] += ((double[]) result.get("blocksd"))[x] / numofsw;
                }
            }
        } else {	// if sliding window is not selected, correlate the full intensity trace
            datac = new double[2][num + 1];							// get the intensity data for the correlation 
            datac[0] = getIntensity(image, px1, py1, 1, initialframe, finalframe);		// getIntensity for first pixel; performs a bleach correction if indicated in the panel
            if (px1 != px2 || py1 != py2) {						// if the two pixels are not equal (i.e. for a cross-correlation)
                datac[1] = getIntensity(image, px2, py2, 2, initialframe, finalframe);	// getIntensity for second pixel
            } else {
                datac[1] = datac[0];			// otherwise perform an autocorrelation
            }

            Map result;
            result = correlator(datac, initialframe, finalframe);		// correlate the data
            acf[kcf][pxm1][pym1] = (double[]) result.get("corav");			// acf
            varacf[kcf][pxm1][pym1] = (double[]) result.get("blockvar");		// variance of the ACF; blocked
            sdacf[kcf][pxm1][pym1] = (double[]) result.get("blocksd");		// standard deviation of the ACF; blocked
            currentCovmats = (double[][]) result.get("covmats");
        }

        blocked[kcf][pxm1][pym1] = blockIndS;			// store whether blocking worked successfully for the pixel

        // calculate MSD if switched on
        if (doMSD) {
            if (!MSDmode) { // 2D if MSDmode is false, otherwise 3D
                msd[kcf][pxm1][pym1] = correlationToMSD(acf[kcf][pxm1][pym1], pixeldimx * Math.pow(10, 6), psfsize * Math.pow(10, 6));
            } else {
                msd[kcf][pxm1][pym1] = correlationToMSD3D(acf[kcf][pxm1][pym1], pixeldimx * Math.pow(10, 6), psfsize * Math.pow(10, 6), lsthickness * Math.pow(10, 6));
            }
        }
    }

    // correlator calculates correlation functions
    public Map correlator(double[][] intcor, int initialframe, int finalframe) {
        // intcor contains the array of intensity values to be correlated for pixels 1 and 2
        // initialframe and finalframe provide the range of frames to be used for the correlation
        int num = (finalframe - initialframe + 1); 			// total number of frames to be correlated
        int blockIndex;									// index at which optimal blocking is reached; if it fails maximal blocking is used

        blockIndex = blockTransform(intcor, num, 1);			// perform blocking on the first channel to determine when intensity bins are independent

        Map result;
        result = calculateCF(intcor, num, blockIndex, 1);			// perform optimal blocking and return the CF, SD and covariance matrix

        return result;
    }

    public int blockTransform(double[][] intcor, int num, int blocklag) {
        // intcor is the array of intensity values for the two traces which are correlated
        // num is the number of frames which are correlated
        // blocklag defines for which lag the blocking will be done; typically we use the smalles, i.e. 1
        int blocknum = (int) Math.floor(Math.log(mtab[blocklag]) / Math.log(2)) - 2;	// number of blocking operations that can be performed given blocklag
        int numbin = num;		// number of data points when they are binned
        int del;				// delay or correlation time expressed in lags
        int currentIncrement;
        int crwin = 2;			// 3 points that fit the error bar overlap criterion 
        double sumprod = 0.0;		// sum of all intensity products; divide by num to get the average <i(n)i(n+del)>
        double sumprod2 = 0.0;	// sum of all intensity products squared; divide by num to get the average <(i(n)i(n+del))^2>
        double directm = 0.0;		// direct monitor required for ACF normalization
        double delayedm = 0.0;	// delayed monitor required for ACF normalization
        double[][] intblock;
        double[] prod = new double[num];
        double[][] varblock;
        double[] upper;
        double[] lower;
        double[][] blockpoints;
        double[] blocksd;
        int[] crt;
        int[] cr12;			// do neighbouring points have overlapping error bars; together with crwin=2 this tests for three points that have overlapping erro bars
        int[] cr3;
        int[] diffpos;
        int last0 = 0;
        int ind = 0;
        double[] prodnum;
        double minblock;
        double maxblock;

        varblock = new double[3][blocknum];
        prodnum = new double[blocknum];
        intblock = new double[2][num];
        blocksd = new double[chanum];
        upper = new double[blocknum];
        lower = new double[blocknum];
        crt = new int[blocknum - 1];
        cr12 = new int[blocknum - 2];
        cr3 = new int[blocknum - 2];
        diffpos = new int[blocknum - 1];
        blockpoints = new double[3][3];

        for (int x = 0; x < numbin; x++) {
            intblock[0][x] = intcor[0][x];
            intblock[1][x] = intcor[1][x];
        }

        currentIncrement = blocklag;		// at the moment we always do blocking for smallest lag which is 1 but in general it can be used freely

        for (int x = 1; x < chanum; x++) {									// run over all channels
            if (currentIncrement != samp[x]) {								// check whether the channel width has changed
                currentIncrement = (int) samp[x];								// set the currentIncrement accordingly
                numbin = (int) Math.floor(numbin / 2);							// and correct the number of actual data points accordingly
                for (int y = 0; y < numbin; y++) {								// if yes, bin the data according to the width of the current channel
                    intblock[0][y] = (intblock[0][2 * y] + intblock[0][2 * y + 1]);
                    intblock[1][y] = (intblock[1][2 * y] + intblock[1][2 * y + 1]);
                }

            }

            if (x == blocklag) {										// if the channel number is equal to the blocklag ...
                del = lag[x] / currentIncrement;							// calculate the delay, i.e. the correlation time
                for (int y = 0; y < numbin - del; y++) {				// calculate the ...
                    directm += intblock[0][y];							// direct and ...
                    delayedm += intblock[1][y + del];					// delayed monitor
                }
                prodnum[0] = numbin - del; 								// number of correlation products
                directm /= prodnum[0];									// calculate average of direct and delayed monitor, 
                delayedm /= prodnum[0];									// i.e. the average intesity <n(0)> and <n(tau)>

                for (int y = 0; y < prodnum[0]; y++) {					// calculate the correlation
                    prod[y] = intblock[0][y] * intblock[1][y + del] - delayedm * intblock[0][y] - directm * intblock[1][y + del] + delayedm * directm;
                    sumprod += prod[y];									// calculate the sum of prod, i.e. the raw correlation value ...
                    sumprod2 += Math.pow(prod[y], 2);					// ... and the sum of the squares
                }

                varblock[0][0] = currentIncrement * frametime;			// the time of the block curve
                varblock[1][0] = (sumprod2 / prodnum[0] - Math.pow(sumprod / prodnum[0], 2)) / (prodnum[0] * Math.pow(directm * delayedm, 2));	// value of the block curve

                for (int y = 1; y < blocknum; y++) {					// perform blocking operations
                    prodnum[y] = (int) Math.floor(prodnum[y - 1] / 2);	// the number of samples for the blocking curve decreases by a factor 2 with every step
                    sumprod = 0;
                    sumprod2 = 0;
                    for (int z = 0; z < prodnum[y]; z++) {			// bin the correlation data and calculate the blocked values for the SD
                        prod[z] = (prod[2 * z] + prod[2 * z + 1]) / 2;
                        sumprod += prod[z];
                        sumprod2 += Math.pow(prod[z], 2);
                    }
                    varblock[0][y] = (currentIncrement * Math.pow(2, y)) * frametime;	// the time of the block curve
                    varblock[1][y] = (sumprod2 / prodnum[y] - Math.pow(sumprod / prodnum[y], 2)) / (prodnum[y] * Math.pow(directm * delayedm, 2));	// value of the block curve
                }
            }
        }

        for (int x = 0; x < blocknum; x++) {
            varblock[1][x] = Math.sqrt(varblock[1][x]);							// calculate the standard deviation
            varblock[2][x] = varblock[1][x] / Math.sqrt(2 * (prodnum[x] - 1));	// calculate the error 
            upper[x] = varblock[1][x] + varblock[2][x];							// upper and lower quartile
            lower[x] = varblock[1][x] - varblock[2][x];
        }

        // determine index where blocking criteria are fulfilled
        for (int x = 0; x < blocknum - 1; x++) {							// do neighboring points have overlapping error bars?
            if (upper[x] > lower[x + 1] && upper[x + 1] > lower[x]) {
                crt[x] = 1;
            }
        }

        for (int x = 0; x < blocknum - 2; x++) {							// do three adjacent points have overlapping error bars?
            if (crt[x] * crt[x + 1] == 1) {
                cr12[x] = 1;
            }
        }

        for (int x = 0; x < blocknum - 1; x++) {							// do neighboring points have a positive difference (increasing SD)?
            if (varblock[1][x + 1] - varblock[1][x] > 0) {
                diffpos[x] = 1;
            }
        }

        for (int x = 0; x < blocknum - 2; x++) {							// do three neighboring points monotonically increase?
            if (diffpos[x] * diffpos[x + 1] == 1) {
                cr3[x] = 1;
            }
        }

        for (int x = 0; x < blocknum - 2; x++) {							// find the last triple of points with monotonically increasing differences and non-overlapping error bars
            if ((cr3[x] == 1 && cr12[x] == 0)) {
                last0 = x;
            }
        }

        for (int x = 0; x <= last0; x++) {								// indices of two pairs that pass criterion 1 an 2
            cr12[x] = 0;
        }

        cr12[blocknum - 3] = 0;												// criterion 3, the last two points can't be part of the blocking triple
        cr12[blocknum - 4] = 0;

        for (int x = blocknum - 5; x > 0; x--) {							// index of triplet with overlapping error bars and after which no other triplet has a significant monotonic increase
            if (cr12[x] == 1) {											// or 4 increasing points
                ind = x + 1;												// take the middle of the three points as the blocking limit
            }
        }

        if (ind == 0) {													// if optimal blocking is not possible, use maximal blocking
            blockIndS = 0;
            if (blocknum - 3 > 0) {
                ind = blocknum - 3; 											// maximal blocking is performed for the 3rd last point in the blocking curve if that exists
            } else {
                ind = blocknum - 1;
            }
        } else {
            blockIndS = 1;
        }

        ind = (int) Math.max(ind, correlatorq - 1);				// block at least until maximum sample time

        // Plot the blocking curve if selected
        if (plotBlockingCurve == true) {
            minblock = varblock[1][0];
            maxblock = varblock[1][0];
            for (int x = 0; x < blocknum; x++) {
                if (varblock[1][x] > maxblock) {
                    maxblock = varblock[1][x];
                }
                if (varblock[1][x] < minblock) {
                    minblock = varblock[1][x];
                }
            }
            minblock *= 0.9;
            maxblock *= 1.1;

            Plot plot = new Plot("blocking", "x", "SD", varblock[0], varblock[1]);
            plot.setFrameSize(blockingWindowDimX, blockingWindowDimY);
            plot.setLogScaleX();
            plot.setLimits(varblock[0][0] / 2, 2 * varblock[0][blocknum - 1], minblock, maxblock);
            plot.setColor(java.awt.Color.BLUE);
            plot.setJustification(Plot.CENTER);
            plot.addPoints(varblock[0], varblock[1], varblock[2], Plot.CIRCLE);
            plot.draw();
            if (ind != 0) {								// Plot the points where blocking is succesful in red
                blockpoints[0][0] = varblock[0][ind - 1];
                blockpoints[1][0] = varblock[1][ind - 1];
                blockpoints[2][0] = varblock[2][ind - 1];
                blockpoints[0][1] = varblock[0][ind];
                blockpoints[1][1] = varblock[1][ind];
                blockpoints[2][1] = varblock[2][ind];
                blockpoints[0][2] = varblock[0][ind + 1];
                blockpoints[1][2] = varblock[1][ind + 1];
                blockpoints[2][2] = varblock[2][ind + 1];
                plot.setColor(java.awt.Color.RED);
                plot.addPoints(blockpoints[0], blockpoints[1], blockpoints[2], Plot.CIRCLE);
            }
            plot.draw();
            // either create a new plot window or plot within the existing window
            if (blockingWindow == null || blockingWindow.isClosed() == true) {
                blockingWindow = plot.show();
                blockingWindow.setLocation(blockingWindowPosX, blockingWindowPosY);
            } else {
                blockingWindow.drawPlot(plot);
            }
        }

        return ind;	// return the blockIndex
    }

    // calculate the standard deviation by blocking
    public Map calculateCF(double[][] intcor, int num, int ind, int blocklag) {
        // intcor is the array of intensity values for the two traces which are correlated
        // num is the number of frames which are correlated
        // ind is the blockindex, at which the SD has converged, previously found in blockSD()
        // blocklag defines for which lag the blocking will be done; typically we use the smalles, i.e. 1
        int numbin = num;		// number of data points when they are binned
        int del;				// delay or correlation time expressed in lags
        int currentIncrement;
        int ctbin = 0;			// count how often the data was binned
        int binct;
        int pnum = (int) Math.floor(mtab[chanum - 1] / Math.pow(2, Math.max(ind - Math.log(samp[chanum - 1]) / Math.log(2), 0)));	// minimum number of prodcuts given the used correlator structure and the blockIndex ind
        int[] prodnum = new int[chanum];
        double sumprod;		// sum of all intensity products; divide by num to get the average <i(n)i(n+del)>
        double sumprod2;	// sum of all intensity products squared; divide by num to get the average <(i(n)i(n+del))^2>
        double[] directm = new double[chanum];		// direct monitor required for ACF normalization
        double[] delayedm = new double[chanum];	// delayed monitor required for ACF normalization
        double[] blockvar;
        double[] blocksd;
        double[][] intblock;
        double[][] prod = new double[chanum][num];
        double[] corav = new double[chanum];
        double[] mcov = new double[chanum];
        double[] diagcovmat = new double[chanum];
        double[][] covmat = new double[chanum][chanum];
        double[][] covmats = new double[chanum - 1][chanum - 1];	//the final results does not contain information about the zero lagtime channel
        double[][] cormat = new double[chanum][chanum];
        double[][] denomshrink = new double[chanum][chanum];
        double numerator = 0;
        double denominator = 0;
        double median = 0;
        double lamvar;
        double lamcov;

        intblock = new double[2][num];
        blockvar = new double[chanum];
        blocksd = new double[chanum];

        for (int x = 0; x < numbin; x++) {	//re-initialize the intensities
            intblock[0][x] = intcor[0][x + 1];
            intblock[1][x] = intcor[1][x + 1];
        }

        currentIncrement = 1;		// at the moment we always do blocking for smallest lag
        blocksd[0] = 0;				// we do not calcualte the SD for the 0 lagtime as it is not used for fitting (shot noise)

        for (int x = 0; x < chanum; x++) {						// run over all channels except the 0 lag time

            if (currentIncrement != samp[x]) {					// check whether the channel width has changed
                numbin = (int) Math.floor(numbin / 2);				// if yes, correct the number of actual data points
                currentIncrement = (int) samp[x];					// set the currentIncrement accordingly
                ctbin++;											// count how often the data was binned
                for (int y = 0; y < numbin; y++) {					// and bin the data according to the width of the current channel
                    intblock[0][y] = (intblock[0][2 * y] + intblock[0][2 * y + 1]);
                    intblock[1][y] = (intblock[1][2 * y] + intblock[1][2 * y + 1]);
                }
            }

            del = lag[x] / currentIncrement;							// calculate the delay, i.e. the correlation time ...
            prodnum[x] = numbin - del;								// and the number of products for that delay; //(int) (mtab[chanum-1]*(samp[chanum-1]/samp[x]));//IJ.log(Double.toString(prodnum[x])); 
            for (int y = 0; y < prodnum[x]; y++) {				// calculate the ...
                directm[x] += intblock[0][y];						// direct and ...
                delayedm[x] += intblock[1][y + del];				// delayed monitor
            }
            directm[x] /= prodnum[x];		// calculate average of direct and delayed monitor, i.e. the average intensity <n(0)> and <n(tau)>
            delayedm[x] /= prodnum[x];

            sumprod = 0;
            sumprod2 = 0;

            for (int y = 0; y < prodnum[x]; y++) {					// calculate the correlation
                prod[x][y] = intblock[0][y] * intblock[1][y + del] - delayedm[x] * intblock[0][y] - directm[x] * intblock[1][y + del] + delayedm[x] * directm[x];
                sumprod += prod[x][y];								// calculate the sum of prod, i.e. the raw correlation value ...
                sumprod2 += Math.pow(prod[x][y], 2);				// ... and the sum of the squares
            }

            corav[x] = sumprod / (prodnum[x] * directm[x] * delayedm[x]);	// calculate the ACF, i.e. the mean for the later calculations of the variance-covariance matrix

            binct = ind - ctbin; // determine whether data needs to be further binned or is already exceeding the blocking number
            sumprod = 0;
            sumprod2 = 0;

            for (int y = 1; y <= binct; y++) {					// bin the data until block time is reached
                prodnum[x] = (int) Math.floor(prodnum[x] / 2);		// for each binning the number of data points is halfed
                for (int z = 0; z < prodnum[x]; z++) {			// do the binning and divide by 2 so that average value does not change
                    prod[x][z] = (prod[x][2 * z] + prod[x][2 * z + 1]) / 2;
                }
            }

            prodnum[x] = pnum;										// use only the minimal number of products to achieve a symmetric variance matrix
            for (int z = 0; z < prodnum[x]; z++) {
                sumprod += prod[x][z];								// calculate the sum of prod, i.e. the raw correlation value ...
                sumprod2 += Math.pow(prod[x][z], 2);				// ... and the sum of the squares
            }

            blockvar[x] = (sumprod2 / prodnum[x] - Math.pow(sumprod / prodnum[x], 2)) / ((prodnum[x] - 1) * Math.pow(directm[x] * delayedm[x], 2));	// variance after blocking; extra division by prodnum to obtain SEM
            blocksd[x] = Math.sqrt(blockvar[x]);																									// standard deviation after blocking
        }

        // if GLS is selected then calulate the regularized covariance matrix
        if (tbGLS.isSelected()) {
            // Calculate the mean of the products used for the variance-covariance matrix
            for (int x = 1; x < chanum; x++) {
                for (int z = 0; z < pnum; z++) {
                    mcov[x] += prod[x][z] / (directm[x] * delayedm[x]);
                }
                mcov[x] /= pnum;		// normalize by the number of products
            }

            // Calculate the variance-covariance matrix
            for (int x = 1; x < chanum; x++) {
                for (int y = 1; y <= x; y++) {	// calculate only the upper triangular part as the matrix is symmetric
                    for (int z = 0; z < pnum; z++) {
                        covmat[x][y] += (prod[x][z] / (directm[x] * delayedm[x]) - mcov[x]) * (prod[y][z] / (directm[y] * delayedm[y]) - mcov[y]);
                    }
                    covmat[x][y] /= (pnum - 1);		// normalize by the number of products
                    covmat[y][x] = covmat[x][y];	// lower triangular part is equal to upper triangular part
                }
            }

            // Regularize variance-covariance matrix
            // first determine the shrinkage weight for the variance 
            for (int x = 1; x < chanum; x++) {			// get the variance (diagonal of covariance matrix) ...
                diagcovmat[x] = covmat[x][x];
            }

            Arrays.sort(diagcovmat);						// ... and determine the median
            double pos1 = Math.floor((diagcovmat.length - 1.0) / 2.0);
            double pos2 = Math.ceil((diagcovmat.length - 1.0) / 2.0);
            if (pos1 == pos2) {
                median = diagcovmat[(int) pos1];
            } else {
                median = (diagcovmat[(int) pos1] + diagcovmat[(int) pos2]) / 2.0;
            }

            double tmpnum;									// determine the variance of the variance
            for (int x = 1; x < chanum; x++) {
                tmpnum = 0;
                for (int z = 0; z < pnum; z++) {
                    tmpnum += Math.pow((Math.pow(prod[x][z] / (directm[x] * delayedm[x]) - mcov[x], 2) - covmat[x][x]), 2);
                }
                tmpnum *= (pnum) / Math.pow(pnum - 1, 3);
                numerator += tmpnum;
                denominator += Math.pow(covmat[x][x] - median, 2);
            }
            lamvar = Math.min(1, numerator / denominator);		// shrinkage weight for the variance	
            lamvar = Math.max(lamvar, 0);

            // determine the shrinkage weight for the covariance
            for (int x = 1; x < chanum; x++) {						// calculate the sample correlation matrix
                for (int y = 1; y < chanum; y++) {
                    cormat[x][y] = covmat[x][y] / Math.sqrt(covmat[x][x] * covmat[y][y]);
                }
            }

            numerator = 0;
            denominator = 0;
            double cmx;										// tmp variables to simplify ... 
            double cmy;										// ... in the loop
            for (int x = 1; x < chanum; x++) {			// determine the variance of the covariance
                tmpnum = 0;
                for (int y = 1; y < x; y++) {				//sum only over the upper triangle as the matrix is symmetric
                    for (int z = 0; z < pnum; z++) {
                        cmx = (prod[x][z] / (directm[x] * delayedm[x]) - mcov[x]) / Math.sqrt(covmat[x][x]);
                        cmy = (prod[y][z] / (directm[y] * delayedm[y]) - mcov[y]) / Math.sqrt(covmat[y][y]);
                        tmpnum += Math.pow(cmx * cmy - cormat[x][y], 2);
                    }
                    tmpnum *= (pnum) / Math.pow(pnum - 1, 3);
                    numerator += tmpnum;
                    denominator += Math.pow(cormat[x][y], 2);		// sum of squares of off-diagonal elements of correlation matrix
                }
            }
            lamcov = Math.min(1, numerator / denominator);		// shrinkage weight for the covariance
            lamcov = Math.max(lamcov, 0);

            // calculate the off-diagonal elements of the regularized variance-covariance matrix
            for (int x = 1; x < chanum; x++) {		// do not include zero lagtime channel as we don't use it for fitting				
                for (int y = 1; y < x; y++) {
                    cmx = lamvar * median + (1 - lamvar) * covmat[x][x];
                    cmy = lamvar * median + (1 - lamvar) * covmat[y][y];
                    covmats[x - 1][y - 1] = (1 - lamcov) * cormat[x][y] * Math.sqrt(cmx * cmy) / pnum;
                    covmats[y - 1][x - 1] = covmats[x - 1][y - 1];
                }
            }
            for (int x = 1; x < chanum; x++) {	// diagonal elements of the regularized variance-covariance matrix
                covmats[x - 1][x - 1] = (lamvar * median + (1 - lamvar) * covmat[x][x]) / pnum;
            }

            // Plot covariance matrix if selected
            if (plotCovmats == true) {
                if (impCov != null) {		// close covariance window if it exists
                    impCov.close();
                }
                impCov = IJ.createImage($impCovTitle, "GRAY32", chanum - 1, chanum - 1, 1);
                impCovIp = impCov.getProcessor();
                for (int x = 1; x < chanum; x++) {		// calculate the covariance
                    for (int y = 1; y < chanum; y++) {
                        impCovIp.putPixelValue(x - 1, y - 1, covmats[x - 1][y - 1]);			// regularized var-cov matrix
                        //impCovIp.putPixelValue(x, y, Math.sqrt(covmat[x][y]));	// non-regularized var-cov matrix
                        //impCovIp.putPixelValue(x, y, Math.sqrt(cormat[x][y]));	// sample correlation matrix
                    }
                }

                impCov.show();
                IJ.run(impCov, "Spectrum", "");							// apply "Spectrum" LUT
                IJ.run(impCov, "Enhance Contrast", "saturated=0.35");	// enhance contrast
                impCovWin = impCov.getWindow();
                impCovWin.setLocation(covPosX, covPosY);
                IJ.run(impCov, "Set... ", "zoom=" + 200 + " x=" + 0 + " y=" + 0); //then zoom to fit within application
                IJ.run("In [+]", ""); 	// This needs to be used since ImageJ 1.48v to set the window to the right size; 
            }
        } // end of if GLS selected statement

        Map<String, Object> map = new HashMap<>();
        if (tbGLS.isSelected()) {	//hand over either the correlation function corav or the actual function used to calcualte the covariance matrix; they differ only slightly
            map.put("corav", mcov);
        } else {
            map.put("corav", corav);
        }
        map.put("blockvar", blockvar);
        map.put("blocksd", blocksd);
        map.put("covmats", covmats);
        return map;

    }

    // get intensity data for correlate() and correct for bleaching if required; note that you should call calcIntensityTrace() before to obtain intTrace1 and 2
    public double[] getIntensity(ImagePlus image, int px, int py, int mode, int initialframe, int finalframe) {
        // image: imp form which intensity will be taken
        // px, py: coordinates of pixel within image
        // mode: determines whether intensity for pixel 1 or pixel 2 is read, and in case of DC-FCCS, 
        // whether background1 or background 2 is to be subtracted from the intensity trace
        // initialframe and finalframe provide the range of frames to be used
        int num = (finalframe - initialframe + 1);
        double[] intdat = new double[num + 1];
        double[] res = new double[5];
        int bckg;

        if (mode == 2 && cbFitModel.getSelectedItem() == "DC-FCCS") {
            bckg = background2;
        } else {
            bckg = background;
        }

        for (int x = 1; x <= num; x++) {	//read data from all relevant pixels, depending on the selected frames and binning
            for (int i = 0; i < binningX; i++) {
                for (int k = 0; k < binningY; k++) {
                    if (bgrloaded) {
                        bckg = (int) Math.round(bgrmean[px + i][py + k]);
                    }
                    intdat[x] += image.getStack().getProcessor(initialframe + x - 1).get(px + i, py + k) - bckg;
                }
            }

        }

        String $bcmode = (String) cbBleachCor.getSelectedItem();	// perform single or double exponential bleach corrections if selected 

        if ("Single Exp".equals($bcmode)) {
            SingleExpFit efit = new SingleExpFit();
            if (mode == 1) {
                res = efit.doFit(intTrace1);			// note that the bleach correction is performed on the averaged intensity traces to make it faster
            } else {									// while you may have 20,000 intensity points, intTrace1 and 2 contain only 1,000 points
                res = efit.doFit(intTrace2);			// see definition in setParameters()
            }
            if (res[0] * res[1] != 0) {				// correct the full intensity trace if the fit was succesful
                for (int x = 1; x <= num; x++) {
                    intdat[x] = intdat[x] / Math.sqrt((res[0] * Math.exp(-frametime * (x + 0.5) / res[1]) + res[2]) / (res[0] + res[2])) + (res[0] + res[2]) * (1 - Math.sqrt((res[0] * Math.exp(-frametime * (x + 0.5) / res[1]) + res[2]) / (res[0] + res[2])));
                }
                if (mode == 1) {
                    for (int x = 0; x < nopit; x++) {
                        intTrace1[x] = intTrace1[x] / Math.sqrt((res[0] * Math.exp(-intTime[x] / res[1]) + res[2]) / (res[0] + res[2])) + (res[0] + res[2]) * (1 - Math.sqrt((res[0] * Math.exp(-intTime[x] / res[1]) + res[2]) / (res[0] + res[2])));
                    }
                }
                if (mode == 2) {
                    for (int x = 0; x < nopit; x++) {
                        intTrace2[x] = intTrace2[x] / Math.sqrt((res[0] * Math.exp(-intTime[x] / res[1]) + res[2]) / (res[0] + res[2])) + (res[0] + res[2]) * (1 - Math.sqrt((res[0] * Math.exp(-intTime[x] / res[1]) + res[2]) / (res[0] + res[2])));
                    }
                }
            } else {
                IJ.log("Exponential Fit not successful for (" + px + ", " + py + ")");
            }
        }

        if ("Double Exp".equals($bcmode)) {		// same as single exponential fit only for double exponential
            DoubleExpFit defit = new DoubleExpFit();
            if (mode == 1) {
                res = defit.doFit(intTrace1);
            } else {
                res = defit.doFit(intTrace2);
            }
            if (res[0] * res[1] * res[2] * res[3] != 0) {
                for (int x = 1; x <= num; x++) {
                    intdat[x] = intdat[x] / Math.sqrt((res[0] * Math.exp(-frametime * (x + 0.5) / res[1]) + res[2] * Math.exp(-frametime * (x + 0.5) / res[3]) + res[4]) / (res[0] + res[2] + res[4])) + (res[0] + res[2] + res[4]) * (1 - Math.sqrt((res[0] * Math.exp(-frametime * (x + 0.5) / res[1]) + res[2] * Math.exp(-frametime * (x + 0.5) / res[3]) + res[4]) / (res[0] + res[2] + res[4])));
                }
                if (mode == 1) {
                    for (int x = 0; x < nopit; x++) {
                        intTrace1[x] = intTrace1[x] / Math.sqrt((res[0] * Math.exp(-intTime[x] / res[1]) + res[2] * Math.exp(-intTime[x] / res[3]) + res[4]) / (res[0] + res[2] + res[4])) + (res[0] + res[2] + res[4]) * (1 - Math.sqrt((res[0] * Math.exp(-intTime[x] / res[1]) + res[2] * Math.exp(-intTime[x] / res[3]) + res[4]) / (res[0] + res[2] + res[4])));
                    }
                }
                if (mode == 2) {
                    for (int x = 0; x < nopit; x++) {
                        intTrace2[x] = intTrace2[x] / Math.sqrt((res[0] * Math.exp(-intTime[x] / res[1]) + res[2] * Math.exp(-intTime[x] / res[3]) + res[4]) / (res[0] + res[2] + res[4])) + (res[0] + res[2] + res[4]) * (1 - Math.sqrt((res[0] * Math.exp(-intTime[x] / res[1]) + res[2] * Math.exp(-intTime[x] / res[3]) + res[4]) / (res[0] + res[2] + res[4])));
                    }
                }
            } else {
                IJ.log("Double Exponential Fit not successful for (" + px + ", " + py + ")");
            }
        }

        if ("Polynomial".equals($bcmode)) {		//fitting with a polynomial of selected order
            PolynomFit polfit = new PolynomFit();
            double corfunc;
            int maxord = polyOrder;
            if (mode == 1) {
                res = polfit.doFit(intTrace1);			// note that the bleach correction is performed on the averaged intensity traces to make it faster
            } else {									// while you may have 20,000 intensity points, intTrace1 and 2 contain only 1,000 points
                res = polfit.doFit(intTrace2);			// see definition in setParameters()
            }

            for (int x = 1; x <= num; x++) {
                corfunc = 0;
                for (int i = 0; i <= maxord; i++) {
                    corfunc += res[i] * Math.pow(frametime * (x + 0.5), i);
                }
                intdat[x] = intdat[x] / Math.sqrt(corfunc / res[0]) + res[0] * (1 - Math.sqrt(corfunc / res[0]));
            }
            if (mode == 1) {
                for (int x = 0; x < nopit; x++) {
                    corfunc = 0;
                    for (int i = 0; i <= maxord; i++) {
                        corfunc += res[i] * Math.pow(intTime[x], i);
                    }
                    intTrace1[x] = intTrace1[x] / Math.sqrt(corfunc / res[0]) + res[0] * (1 - Math.sqrt(corfunc / res[0]));
                }
            }
            if (mode == 2) {
                for (int x = 0; x < nopit; x++) {
                    corfunc = 0;
                    for (int i = 0; i <= maxord; i++) {
                        corfunc += res[i] * Math.pow(intTime[x], i);
                    }
                    intTrace2[x] = intTrace2[x] / Math.sqrt(corfunc / res[0]) + res[0] * (1 - Math.sqrt(corfunc / res[0]));
                }
            }
        }

        if ("Lin Segment".equals($bcmode)) {		// approximating bleaching by a partially linear function
            int ave = (int) Math.floor((finalframe - initialframe + 1) / nopit); // number of points averaged in intTrace
            int bcNopit = (int) num / slidingWindowLength; // number of linear segments
            int bcAve = (int) Math.floor((finalframe - initialframe + 1) / bcNopit);
            double[] bcTrace = new double[bcNopit];
            int nf;
            double bcInt0;
            double sum;

            for (int x = 0; x < bcNopit; x++) {		// calculating the average intensity in each segment
                sum = 0;	// initialize arrays with 0
                for (int y = 1 + x * bcAve; y < (x + 1) * bcAve; y++) {
                    sum += intdat[y];
                }
                bcTrace[x] = sum / bcAve;
            }

            bcInt0 = bcTrace[0] + (bcTrace[0] - bcTrace[1]) / 2;	// Initial intensity obtained by extrapolating the line between average intensities in 1st and second segment to 0

            for (int x = 1; x < (int) Math.floor(bcAve / 2); x++) {
                intdat[x] = intdat[x] / Math.sqrt((bcTrace[0] + (bcTrace[0] - bcTrace[1]) / bcAve * (bcAve / 2 - x)) / bcInt0) + bcInt0 * (1 - Math.sqrt((bcTrace[0] + (bcTrace[0] - bcTrace[1]) / bcAve * (bcAve / 2 - x)) / bcInt0));
            }

            for (int x = 1; x < bcNopit; x++) {
                for (int y = 0; y < bcAve; y++) {
                    nf = (x - 1) * bcAve + (int) Math.floor(bcAve / 2) + y;
                    intdat[nf] = intdat[nf] / Math.sqrt((bcTrace[x - 1] + (bcTrace[x] - bcTrace[x - 1]) * y / bcAve) / bcInt0) + bcInt0 * (1 - Math.sqrt((bcTrace[x - 1] + (bcTrace[x] - bcTrace[x - 1]) * y / bcAve) / bcInt0));
                }
            }

            for (int x = (bcNopit - 1) * bcAve + (int) Math.floor(bcAve / 2); x <= num; x++) {
                nf = x - (bcNopit - 1) * bcAve + (int) Math.floor(bcAve / 2) + bcAve;
                intdat[x] = intdat[x] / Math.sqrt((bcTrace[bcNopit - 2] + (bcTrace[bcNopit - 1] - bcTrace[bcNopit - 2]) * nf / bcAve) / bcInt0) + bcInt0 * (1 - Math.sqrt((bcTrace[bcNopit - 2] + (bcTrace[bcNopit - 1] - bcTrace[bcNopit - 2]) * nf / bcAve) / bcInt0));
            }

            if (mode == 1) {
                for (int x = 0; x < nopit; x++) {
                    intTrace1[x] = intdat[(int) (x + 0.5) * ave + 1];
                }
            }
            if (mode == 2) {
                for (int x = 0; x < nopit; x++) {
                    intTrace2[x] = intdat[(int) (x + 0.5) * ave + 1];
                }
            }
        }

        return (intdat);
    }

    // calculate reduced intensity traces for plotting average and use not more than 'nopit' (defined in setParameters()) points for a trace
    public void calcIntensityTrace(ImagePlus image, int ipx1, int ipy1, int ipx2, int ipy2, int initialframe, int finalframe) {
        // image: imp form which intensity will be taken
        // px1, py1, px2, py2: coordinates of pixels to be correlated
        // initialframe and finalframe provide the range of frames to be used
        int ave = (int) Math.floor((finalframe - initialframe + 1) / nopit); // calculate number of data points which are averaged
        int sum1;
        int sum2;
        int bckg1 = background;
        int bckg2 = background;	//needs to be adapted for background2 once available

        if (cbFitModel.getSelectedItem() == "DC-FCCS") {
            bckg2 = background2;
        }

        for (int x = 0; x < nopit; x++) {
            sum1 = 0;	// initialize arrays with 0
            sum2 = 0;
            for (int i = 0; i < binningX; i++) {
                for (int k = 0; k < binningY; k++) {
                    for (int y = initialframe + x * ave; y <= initialframe + (x + 1) * ave - 1; y++) {
                        if (bgrloaded) { // if a background image is loaded, then subtract the mean of the background image for each pixel
                            bckg1 = (int) Math.round(bgrmean[ipx1 + i][ipy1 + k]);
                            bckg2 = (int) Math.round(bgrmean[ipx1 + i][ipy1 + k]);
                        }
                        sum1 += image.getStack().getProcessor(y).get(ipx1 + i, ipy1 + k) - bckg1;
                        sum2 += image.getStack().getProcessor(y).get(ipx2 + i, ipy2 + k) - bckg2;
                    }
                }
            }
            intTime[x] = frametime * (x + 0.5) * ave;
            intTrace1[x] = sum1 / ave;	// calculate average intensity for the 'ave' points
            intTrace2[x] = sum2 / ave;
        }
    }

    // calcualte the average intensity for the image
    public void calcAverageIntensityTrace(float[][] filterArray, int startX, int startY, int endX, int endY, int initialframe, int finalframe) {
        // introi: roi over which average is to be determined
        // initialframe and finalframe provide the range of frames to be used

        if (!plotIntensityCurves) {
            return;	// perform only if intensities are to be plotted
        }
        int ave;
        int pixcount = 0;
        ave = (int) Math.floor((finalframe - initialframe + 1) / nopit);
        intTrace1 = new double[nopit]; // reset before updating.

        for (int i = 0; i < nopit; i++) {
            for (int j = firstframe + i * ave; j <= firstframe + (i + 1) * ave - 1; j++) {
                pixcount = 0;
                for (int x1 = startX; x1 <= endX; x1 = x1 + pixbinX) {
                    for (int x2 = startY; x2 <= endY; x2 = x2 + pixbinY) {
                        for (int x4 = 0; x4 < binningX; x4++) {
                            for (int x5 = 0; x5 < binningY; x5++) {
                                if (!Float.isNaN(filterArray[x1][x2])) {
                                    intTrace1[i] += imp.getStack().getProcessor(j).get(x1 + x4, x2 + x5) - background;
                                    pixcount++;
                                    intTime[i] = frametime * (i + 0.5) * ave;
                                }
                            }
                        }
                    }
                }
            }
            intTrace1[i] /= (ave * pixcount);	// calculate average intensity for the 'ave' points
        }
    }

    // calculate and fit average ACF
    public void calcAveCF() {
        int count = 0;
        double correlatedwidth = 0; //parameters to account for space correlations and binning
        double correlatedheight = 0;
        Roi tmproi = new Roi(0, 0, 0, 0);

        if (!setImp) {
            IJ.showMessage("No image open.");
            return;
        }

        try {
            for (int k = 0; k < chanum; k++) {
                aveacf[k] = 0.0;
                varaveacf[k] = 0.0;
                correlatedwidth = width - Math.abs(cfXshift);
                correlatedheight = height - Math.abs(cfYshift);
                if (!overlap) {
                    correlatedwidth /= binningX;
                    correlatedheight /= binningY;
                }

                if (doFit) {
                    for (int i = 0; i < (correlatedwidth); i++) {
                        for (int j = 0; j < (correlatedheight); j++) {
                            if (acf[0][i][j][k] != 0 && !Double.isNaN(acf[0][i][j][k])) { // use this if only fitted pixels to be averaged:  pixvalid[0][i][j] == 1.0
                                aveacf[k] += acf[0][i][j][k];
                                varaveacf[k] += Math.pow(acf[0][i][j][k], 2.0);
                                count++;
                            }
                        }
                    }
                } else {
                    for (int i = 0; i < (correlatedwidth); i++) {
                        for (int j = 0; j < (correlatedheight); j++) {
                            if (acf[0][i][j][k] != 0 && !Double.isNaN(acf[0][i][j][k])) {
                                aveacf[k] += acf[0][i][j][k];
                                varaveacf[k] += Math.pow(acf[0][i][j][k], 2.0);
                                count++;
                            }
                        }
                    }
                }

                if (count > 0) {
                    aveacf[k] /= count;
                    varaveacf[k] /= count;

                }
                varaveacf[k] = varaveacf[k] - Math.pow(aveacf[k], 2.0);
                count = 0;
            }

            if (doMSD) {
                if (!MSDmode) { // 2D if MSDmode is false, otherwise 3D
                    msdaveacf = correlationToMSD(aveacf, pixeldimx * Math.pow(10, 6), psfsize * Math.pow(10, 6));
                } else {
                    msdaveacf = correlationToMSD3D(aveacf, pixeldimx * Math.pow(10, 6), psfsize * Math.pow(10, 6), lsthickness * Math.pow(10, 6));
                }
                plotMSD(tmproi, 4, false);
            }

            if (doFit) {
                prepareFit();
                if (initparam[1] < 0 || initparam[6] < 0 || initparam[8] < 0 || initparam[5] < 0 || initparam[7] < 0) {
                    JOptionPane.showMessageDialog(null, "Fit parameters are out of range. Either D < 0 or F < 0.");
                    return;
                }
                averageFit afit = new averageFit();
                afit.doFit(cbFitModel.getSelectedItem().toString());
            }
            plotCF(tmproi, 4, false);
        } catch (Exception e) {/* just ignore the command */        }
    }

    /* 
 * Fitting procedures:
 * 
 * public void prepareFit(): Reads current status of the parameters given, which are used as initial values in the fits and determines which ones are free and which ones are held constant.
 * public void fit(int cpx, int cpy, int kcf, String model): Performs the actual fit for correlation data
 * public void bayesFit(int cpx, int cpy, int kcf, String model)
 * public class standardFit extends AbstractCurveFitter: perform least suqares fit using the blocked SD
 * public class GLSFit extends AbstractCurveFitter: perform generalized least squares fit using the regularized covarinace matrix
 * public class LineFit extends AbstractCurveFitter: Line fit used for the diffusion law plot
 * public class SingleExpFit extends AbstractCurveFitter: single exponential fit for bleach correction
 * public class DoubleExpFit extends AbstractCurveFitter: double exponential fit for bleach correction
 * public class PolynomFit extends AbstractCurveFitter: perform a polynomial fit for bleach correction
 * public double[][] dlfit(int numbin, int numpsf): a procedure that performs the full diffusion law calculations for the PSF determination
 * public double[] correlationToMSD(double[] corrfunc, double psize, double psfwidth): Calculate the MSD form the ACF or CCF
 * public double cubicrootFunction(double a, double b, double c): solver required in correlationToMSD
 * public double[] quarticFunction(double a, double b, double c, double d): solver required in correlationToMSD
 * 
     */
    // data fitting
    public void prepareFit() {

        String $fitmode = (String) cbFitModel.getSelectedItem();

        // always start with the same initial values; the initial values are first defined in createImFCSFit(); 
        // subsequently they are read in setParameters() from the fit window. One thus should set the prameters to 
        // realistic values, otherwise the fit might not converge and hang the plugin for a long time
        paraminitval[0] = initparam[0];
        paramfit[0] = !rbtnHoldN.isSelected();

        paraminitval[1] = initparam[1];
        paramfit[1] = !rbtnHoldD.isSelected();

        paraminitval[2] = initparam[2];
        paramfit[2] = !rbtnHoldVx.isSelected();

        paraminitval[3] = initparam[3];
        paramfit[3] = !rbtnHoldVy.isSelected();

        paraminitval[4] = initparam[4];
        paramfit[4] = !rbtnHoldG.isSelected();

        paraminitval[5] = initparam[5];
        paramfit[5] = !rbtnHoldF2.isSelected();

        paraminitval[6] = initparam[6];
        paramfit[6] = !rbtnHoldD2.isSelected();

        paraminitval[7] = initparam[7];
        paramfit[7] = !rbtnHoldF3.isSelected();

        paraminitval[8] = initparam[8];
        paramfit[8] = !rbtnHoldD3.isSelected();

        paraminitval[9] = initparam[9];
        paramfit[9] = !rbtnHoldFtrip.isSelected();

        paraminitval[10] = initparam[10];
        paramfit[10] = !rbtnHoldTtrip.isSelected();

        String $q2 = tfParamQ2.getText();
        q2 = Double.parseDouble($q2);

        String $q3 = tfParamQ3.getText();
        q3 = Double.parseDouble($q3);

        pixeldimx = (pixelsize * 1000 / objmag * binningX) / Math.pow(10, 9);
        pixeldimy = (pixelsize * 1000 / objmag * binningY) / Math.pow(10, 9);
        if (Double.parseDouble(tfSigmaZ.getText()) <= 0.0 || Double.parseDouble(tfSigmaZ.getText()) > 100) {
            fitobsvol = obsvolFCS_ST2D1p(2);
        } else {
            fitobsvol = obsvolFCS_ST2D1p(3);
        }
    }

    public boolean fit(int cpx, int cpy, int kcf, String model) {
        pixfitted[kcf][cpx][cpy] = false;
        if (tbGLS.isSelected()) {						// select whether GLS or OLS fits are used and whether Bayes is used or not
            if (tbBayes.isSelected()) {					// Bayes fitting will call fitting routines multiple times for model comparisons
                String $testq2 = tfParamQ2.getText();
                double testq2 = Double.parseDouble($testq2);
                String $testq3 = tfParamQ3.getText();
                double testq3 = Double.parseDouble($testq3);
                if (testq2 * testq3 <= 0) {
                    JOptionPane.showMessageDialog(null, "Q2 and/or Q3 are not set correctly; both have to be larger than 0.");
                    return false;
                }
                if (initparam[1] < 0 || initparam[6] < 0 || initparam[8] < 0 || initparam[5] < 0 || initparam[7] < 0) {
                    JOptionPane.showMessageDialog(null, "Fit parameters are out of range. Either D < 0 or F < 0.");
                    return false;
                }
                bayesFit(cpx, cpy, kcf, model);
            } else {
                if (initparam[1] < 0 || initparam[6] < 0 || initparam[8] < 0 || initparam[5] < 0 || initparam[7] < 0) {
                    JOptionPane.showMessageDialog(null, "Fit parameters are out of range. Either D < 0 or F < 0.");
                    return false;
                }
                GLSFit glsfit = new GLSFit();
                glsfit.doFit(cpx, cpy, kcf, model);
            }
        } else {
            if (tbBayes.isSelected()) {
                String $testq2 = tfParamQ2.getText();
                double testq2 = Double.parseDouble($testq2);
                String $testq3 = tfParamQ3.getText();
                double testq3 = Double.parseDouble($testq3);
                if (testq2 * testq3 <= 0) {
                    JOptionPane.showMessageDialog(null, "Q2 and/or Q3 are not set correctly; both have to be larger than 0.");
                    return false;
                }
                if (initparam[1] < 0 || initparam[6] < 0 || initparam[8] < 0 || initparam[5] < 0 || initparam[7] < 0) {
                    JOptionPane.showMessageDialog(null, "Fit parameters are out of range. Either D < 0 or F < 0.");
                    return false;
                }
                bayesFit(cpx, cpy, kcf, model);
            } else {
                if (initparam[1] < 0 || initparam[6] < 0 || initparam[8] < 0 || initparam[5] < 0 || initparam[7] < 0) {
                    JOptionPane.showMessageDialog(null, "Fit parameters are out of range. Either D < 0 or F < 0.");
                    return false;
                }
                standardFit sfit = new standardFit();
                sfit.doFit(cpx, cpy, kcf, model);
            }
        }
        return true;
    }

    // fit procedure for Bayes Hypothesis testing
    public void bayesFit(int cpx, int cpy, int kcf, String model) { // at the moment we will allow onle one and two-component fits
        // cpx, cpy: pixel coordinates
        // kcf (0, 1, 2): ACF1, ACF2, CCF
        // model: which fit model to be used; the models are defined below
        ParametricUnivariateFunction fitfunction;
        ParametricUnivariateFunction plotfunction;
        double[] paramem = new double[noparam];
        boolean[] paramfitmem = new boolean[noparam];
        double[][] covpar;
        double[] residuals;
        double[] sigma;
        double[] modprob = new double[2];
        double normprob;
        double residsum;
        double logresid;
        double prodsigma;
        double det;
        double pi = 3.14159265359;
        int boxsize = 200;
        int numel;

        // take the existing fitparameters and do a one or two component fit; if the D2 etc. are not specified, specify them and take D2 as 10 time D1 etc.
        // remember the fit parameter settings so they can be reset after the Bayes fitting
        paramem[0] = Double.parseDouble(tfParamN.getText());
        paramem[1] = Double.parseDouble(tfParamD.getText()) / Math.pow(10, 12);
        paramem[2] = Double.parseDouble(tfParamVx.getText()) / Math.pow(10, 6);
        paramem[3] = Double.parseDouble(tfParamVy.getText()) / Math.pow(10, 6);
        paramem[4] = Double.parseDouble(tfParamG.getText());
        paramem[5] = Double.parseDouble(tfParamF2.getText());
        paramem[6] = Double.parseDouble(tfParamD2.getText()) / Math.pow(10, 12);

        System.arraycopy(paramfit, 0, paramfitmem, 0, noparam);

        // set initial values for one-component fit
        paramfit[0] = true;
        paramfit[1] = true;
        paramfit[2] = false;
        paramfit[3] = false;
        paramfit[4] = true;
        paramfit[5] = false;
        paramfit[6] = false;
        paramfit[7] = false;
        paramfit[8] = false;
        initparam[0] = paramem[0];
        rbtnHoldN.setSelected(!paramfit[0]);
        initparam[1] = paramem[1];
        rbtnHoldD.setSelected(!paramfit[1]);
        initparam[2] = 0;
        rbtnHoldVx.setSelected(true);
        initparam[3] = 0;
        rbtnHoldVy.setSelected(true);
        initparam[4] = paramem[4];
        rbtnHoldG.setSelected(!paramfit[4]);
        initparam[5] = 0;
        rbtnHoldF2.setSelected(true);
        initparam[6] = 0;
        rbtnHoldD2.setSelected(true);
        initparam[7] = 0;
        rbtnHoldF3.setSelected(true);
        initparam[8] = 0;
        rbtnHoldD3.setSelected(true);
        initparam[9] = 0;
        rbtnHoldFtrip.setSelected(true);
        initparam[10] = 0;
        rbtnHoldTtrip.setSelected(true);

        residsum = 0;
        if (!tbGLS.isSelected()) {
            standardFit sfit = new standardFit();
            Map map;
            map = sfit.doFit(cpx, cpy, kcf, model);

            covpar = (double[][]) map.get("covariance");
            residuals = (double[]) map.get("residuals");
            sigma = (double[]) map.get("sigma");

            for (int i = 1; i < residuals.length; i++) {
                residsum += Math.pow(residuals[i], 2);
            }

            logresid = -0.5 * residsum;
        } else {
            GLSFit glsfit = new GLSFit();
            Map map;
            map = glsfit.doFit(cpx, cpy, kcf, model);

            covpar = (double[][]) map.get("covariance");
            residuals = (double[]) map.get("residuals");
            sigma = (double[]) map.get("sigma");

            RealMatrix mat = MatrixUtils.createRealMatrix(currentCovmats);
            RealMatrix matT = mat.transpose();
            DecompositionSolver solver = new LUDecomposition(matT).getSolver();
            RealVector resvec = new ArrayRealVector(residuals);
            RealVector solution = solver.solve(resvec);
            logresid = -0.5 * solution.dotProduct(resvec);
        }

        numel = 3;
        det = determinant(covpar, numel);

        prodsigma = 1;
        for (int i = 0; i < sigma.length; i++) {
            prodsigma *= sigma[i];
        }

        modprob[0] = Math.exp(0.5 * numel * Math.log(2 * pi) + 0.5 * Math.log(det) + logresid - Math.log(prodsigma * Math.pow(2 * boxsize, numel)));

        // set initial values for two-component fit
        paramfit[0] = true;
        paramfit[1] = true;
        paramfit[2] = false;
        paramfit[3] = false;
        paramfit[4] = true;
        paramfit[5] = true;
        paramfit[6] = true;
        paramfit[7] = false;
        paramfit[8] = false;
        initparam[0] = paramem[0];
        rbtnHoldN.setSelected(!paramfit[0]);
        initparam[1] = paramem[1];
        rbtnHoldD.setSelected(!paramfit[1]);
        initparam[2] = 0;
        rbtnHoldVx.setSelected(true);
        initparam[3] = 0;
        rbtnHoldVy.setSelected(true);
        initparam[4] = paramem[4];
        rbtnHoldG.setSelected(!paramfit[4]);
        initparam[5] = 0.5;
        rbtnHoldF2.setSelected(!paramfit[5]);
        initparam[6] = 0.1 * paramem[1];
        rbtnHoldD2.setSelected(!paramfit[6]);
        initparam[7] = 0;
        rbtnHoldF3.setSelected(true);
        initparam[8] = 0;
        rbtnHoldD3.setSelected(true);
        initparam[9] = 0;
        rbtnHoldFtrip.setSelected(true);
        initparam[10] = 0;
        rbtnHoldTtrip.setSelected(true);

        residsum = 0;
        if (!tbGLS.isSelected()) {
            standardFit sfit = new standardFit();
            Map map;
            map = sfit.doFit(cpx, cpy, kcf, model);

            covpar = (double[][]) map.get("covariance");
            residuals = (double[]) map.get("residuals");
            sigma = (double[]) map.get("sigma");

            for (int i = 1; i < residuals.length; i++) {
                residsum += Math.pow(residuals[i], 2);
            }

            logresid = -0.5 * residsum;
        } else {
            GLSFit glsfit = new GLSFit();
            Map map;
            map = glsfit.doFit(cpx, cpy, kcf, model);

            covpar = (double[][]) map.get("covariance");
            residuals = (double[]) map.get("residuals");
            sigma = (double[]) map.get("sigma");

            RealMatrix mat = MatrixUtils.createRealMatrix(currentCovmats);
            RealMatrix matT = mat.transpose();
            DecompositionSolver solver = new LUDecomposition(matT).getSolver();
            RealVector resvec = new ArrayRealVector(residuals);
            RealVector solution = solver.solve(resvec);
            logresid = -0.5 * solution.dotProduct(resvec);
        }

        numel = 5;
        det = determinant(covpar, numel);

        prodsigma = 1;
        for (int i = 0; i < sigma.length; i++) {
            prodsigma *= sigma[i];
        }

        modprob[1] = Math.exp(0.5 * numel * Math.log(2 * pi) + 0.5 * Math.log(det) + logresid - Math.log(prodsigma * Math.pow(2 * boxsize, numel)));

        normprob = 0;
        for (int i = 0; i < modprob.length; i++) // calculate the normalization for the modle probabilities
        {
            normprob += modprob[i];
        }

        tfModProb1.setText(IJ.d2s(modprob[0] / normprob));
        tfModProb2.setText(IJ.d2s(modprob[1] / normprob));
    }

    // standard Nonlinear Least Squares fit
    public class standardFit extends AbstractCurveFitter {

        private int numfreefitpar;

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess;

            int i = 0;				// add data 
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                weights[i] = point.getWeight();
                i++;
            }

            int numparfit = 0;		// count how many paramters will be fit
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    numparfit++;
                }
            }

            numfreefitpar = numparfit;

            initialGuess = new double[numfreefitpar];	// setup the initial guesses for the parameters
            int num = 0;
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    initialGuess[num] = paraminitval[i];
                    num++;
                }
            }

            ParametricUnivariateFunction fitfunction;
            if (currentmodel.equals("FCS")) {	// select the fit model to be used; extra fit models can be added here
                fitfunction = new FCS_3p();
            } else {
                fitfunction = new FCCS_2p();
            }

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(fitfunction, points);

            return new LeastSquaresBuilder().
                    maxEvaluations(fitMaxIterations).
                    maxIterations(fitMaxEvaluations).
                    start(initialGuess).
                    target(target).
                    weight(new DiagonalMatrix(weights)).
                    model(model.getModelFunction(), model.getModelFunctionJacobian()).
                    build();
        }

        public Map doFit(int cpx, int cpy, int kcf, String model) {
            //standardFit fitter = new standardFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();

            currentmodel = model;	//set the presently used model

            ParametricUnivariateFunction fitfunction;
            if (currentmodel.equals("FCS")) {	// select the fit model to be used; extra fit models can be added here
                fitfunction = new FCS_3p();
            } else {
                fitfunction = new FCCS_2p();
            }

            fitstart = Integer.parseInt(tfFitStart.getText());
            fitend = Integer.parseInt(tfFitEnd.getText());

            if (fitstart < 1) {
                IJ.showMessage("Illegal Value for Fit start");
                tfFitStart.setText("1");
                fitstart = 1;
            }
            if (fitend >= chanum) {
                IJ.showMessage("Illegal Value for Fit end");
                tfFitEnd.setText(Integer.toString(chanum - 1));
                fitend = chanum - 1;
            }

            // Add points here
            for (int i = fitstart; i <= fitend; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1 / varacf[kcf][cpx][cpy][i], lagtime[i], acf[kcf][cpx][cpy][i]);
                points.add(point);
            }

            Map<String, Object> map = new HashMap<>();
            try {
                final Optimum topt = getOptimizer().optimize(getProblem(points));
                double result[] = topt.getPoint().toArray();
                double tmpres[] = topt.getResiduals().toArray();
                double[] tres = new double[chanum - 1];

                // store results and create fit and residual function for the plot
                for (int i = 1; i < chanum; i++) {
                    if (i >= fitstart && i <= fitend) {
                        fitacf[kcf][cpx][cpy][i] = fitfunction.value(lagtime[i], result);  // calculate the fit function

                        if (i == 0) {
                            res[kcf][cpx][cpy][i] = acf[kcf][cpx][cpy][i] - fitacf[kcf][cpx][cpy][i];	// calculate the residuals
                        } else {
                            res[kcf][cpx][cpy][i] = tmpres[i - fitstart];							// use the weighted residuals
                            tres[i - 1] = acf[kcf][cpx][cpy][i] - fitacf[kcf][cpx][cpy][i];			// unweighted residuals for model probability calculations
                        }
                    } else {
                        fitacf[kcf][cpx][cpy][i] = 0;
                        res[kcf][cpx][cpy][i] = 0;
                        if (i > 0) {
                            tres[i - 1] = 0;
                        }
                    }
                }

                chi2[kcf][cpx][cpy] = 0; // initialize chi2 value for this pixel	
                for (int i = fitstart; i <= fitend; i++) {
                    chi2[kcf][cpx][cpy] += Math.pow(res[kcf][cpx][cpy][i], 2) / ((fitend - fitstart) - numfreefitpar - 1);	// calculate chi2 value; do not include the 0 lagtime channel which contains shot noise
                }

                int num = 0;										// store the fit results 
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        fitres[kcf][cpx][cpy][i] = result[num];		// for free parameters use the fit results
                        num++;
                    } else {
                        fitres[kcf][cpx][cpy][i] = paraminitval[i];	// for fixed parameters use the initial values
                    }
                }
                pixfitted[kcf][cpx][cpy] = true;

                // update parameters in fit window
                if (tbShow.isSelected()) {
                    tfParamN.setText(IJ.d2s(fitres[kcf][cpx][cpy][0], decformat));
                    tfParamD.setText(IJ.d2s(fitres[kcf][cpx][cpy][1] * Math.pow(10, 12), decformat));
                    tfParamVx.setText(IJ.d2s(fitres[kcf][cpx][cpy][2] * Math.pow(10, 6), decformat));
                    tfParamVy.setText(IJ.d2s(fitres[kcf][cpx][cpy][3] * Math.pow(10, 6), decformat));
                    tfParamG.setText(IJ.d2s(fitres[kcf][cpx][cpy][4], decformat2));
                    tfParamF2.setText(IJ.d2s(fitres[kcf][cpx][cpy][5], decformat));
                    tfParamD2.setText(IJ.d2s(fitres[kcf][cpx][cpy][6] * Math.pow(10, 12), decformat));
                    tfParamF3.setText(IJ.d2s(fitres[kcf][cpx][cpy][7], decformat));
                    tfParamD3.setText(IJ.d2s(fitres[kcf][cpx][cpy][8] * Math.pow(10, 12), decformat));
                    tfParamFtrip.setText(IJ.d2s(fitres[kcf][cpx][cpy][9], decformat));
                    tfParamTtrip.setText(IJ.d2s(fitres[kcf][cpx][cpy][10] * Math.pow(10, 6), decformat));
                }

                map.put("results", result);
                map.put("covariance", topt.getCovariances(1).getData());
                map.put("residuals", tres);
                map.put("sigma", topt.getSigma(1).toArray());

            } catch (Exception e) {
                IJ.log(e.getClass().getName() + " at pixel " + Integer.toString(cpy) + " - " + Integer.toString(cpx));
                ArrayList<Double> result = new ArrayList<>();
                double[] tres = new double[chanum - 1];
                Arrays.fill(tres, 1);
                double[] tsig = new double[chanum - 1];
                Arrays.fill(tsig, 1);
                int num = 0;
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        fitres[kcf][cpx][cpy][i] = Double.NaN;
                        result.add(initparam[i]);		// return the initial values
                        num++;
                    }
                }
                for (int i = 1; i < chanum; i++) {
                    fitacf[kcf][cpx][cpy][i] = 0.0;
                    res[kcf][cpx][cpy][i] = 0.0;	// calculate the residuals
                }
                double[][] tcov = new double[num][num];
                for (double[] tcov1 : tcov) {
                    Arrays.fill(tcov1, 1);
                }
                map.put("results", result.toArray());
                map.put("covariance", tcov);
                map.put("residuals", tres);
                map.put("sigma", tsig);
                pixfitted[kcf][cpx][cpy] = false;
            }

            return map;
        }
    }

    // Nonlinear Least Squares fit for the average ACF
    public class averageFit extends AbstractCurveFitter {

        private int numfreefitpar;

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess;

            int i = 0;				// add data 
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                weights[i] = point.getWeight();
                i++;
            }

            int numparfit = 0;		// count how many paramters will be fit
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    numparfit++;
                }
            }

            numfreefitpar = numparfit;

            initialGuess = new double[numfreefitpar];	// setup the initial guesses for the parameters
            int num = 0;
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    initialGuess[num] = paraminitval[i];
                    num++;
                }
            }

            ParametricUnivariateFunction fitfunction;
            if (currentmodel.equals("FCS")) {	// select the fit model to be used; extra fit models can be added here
                fitfunction = new FCS_3p();
            } else {
                fitfunction = new FCCS_2p();
            }

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(fitfunction, points);

            return new LeastSquaresBuilder().
                    maxEvaluations(fitMaxIterations).
                    maxIterations(fitMaxEvaluations).
                    start(initialGuess).
                    target(target).
                    weight(new DiagonalMatrix(weights)).
                    model(model.getModelFunction(), model.getModelFunctionJacobian()).
                    build();
        }

        public Map doFit(String model) {
            //standardFit fitter = new standardFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();

            currentmodel = model;	//set the presently used model

            ParametricUnivariateFunction fitfunction;
            if (currentmodel.equals("FCS")) {	// select the fit model to be used; extra fit models can be added here
                fitfunction = new FCS_3p();
            } else {
                fitfunction = new FCCS_2p();
            }

            // Add points here
            for (int i = 1; i < chanum; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1 / varaveacf[i], lagtime[i], aveacf[i]);
                points.add(point);
            }

            Map<String, Object> map = new HashMap<>();
            try {
                final Optimum topt = getOptimizer().optimize(getProblem(points));
                double result[] = topt.getPoint().toArray();
                double tmpres[] = topt.getResiduals().toArray();

                // store results and create fit and residual function for the plot
                for (int i = 0; i < chanum; i++) {
                    fitaveacf[i] = fitfunction.value(lagtime[i], result);			// calculate the fit function
                    if (i == 0) {
                        resaveacf[i] = aveacf[i] - fitaveacf[i];	// calculate the residuals
                    } else {
                        resaveacf[i] = tmpres[i - 1];									// use the weighted residuals
                    }
                }

                chi2aveacf = 0; // initialize chi2 value for this pixel	
                for (int i = 1; i < chanum; i++) {
                    chi2aveacf += Math.pow(resaveacf[i], 2) / (chanum - numfreefitpar - 1);	// calculate chi2 value; do not include the 0 lagtime channel which contains shot noise
                }

                int num = 0;										// store the fit results
                double[] avefitres = new double[noparam];
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        avefitres[i] = result[num];		// for free parameters use the fit results
                        num++;
                    } else {
                        avefitres[i] = paraminitval[i];	// for fixed parameters use the initial values
                    }
                }

                // update parameters in fit window
                if (tbShow.isSelected()) {
                    tfParamN.setText(IJ.d2s(avefitres[0], decformat));
                    tfParamD.setText(IJ.d2s(avefitres[1] * Math.pow(10, 12), decformat));
                    tfParamVx.setText(IJ.d2s(avefitres[2] * Math.pow(10, 6), decformat));
                    tfParamVy.setText(IJ.d2s(avefitres[3] * Math.pow(10, 6), decformat));
                    tfParamG.setText(IJ.d2s(avefitres[4], decformat2));
                    tfParamF2.setText(IJ.d2s(avefitres[5], decformat));
                    tfParamD2.setText(IJ.d2s(avefitres[6] * Math.pow(10, 12), decformat));
                    tfParamF3.setText(IJ.d2s(avefitres[7], decformat));
                    tfParamD3.setText(IJ.d2s(avefitres[8] * Math.pow(10, 12), decformat));
                    tfParamFtrip.setText(IJ.d2s(avefitres[9], decformat));
                    tfParamTtrip.setText(IJ.d2s(avefitres[10] * Math.pow(10, 6), decformat));
                }

                map.put("results", result);
                map.put("covariance", topt.getCovariances(1).getData());

            } catch (Exception e) {
                IJ.log(e.getClass().getName() + " for average ACF ");
                ArrayList<Double> result = new ArrayList<>();
                int num = 0;
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        result.add(initparam[i]);		// return the initial values
                        num++;
                    }
                }
                for (int i = 0; i < chanum; i++) {
                    fitaveacf[i] = 0.0;
                    resaveacf[i] = 0.0;	// calculate the residuals

                }

                map.put("results", result.toArray());
            }

            return map;
        }
    }

    // Generalized Least Squares fit; requires covariance matrix as calculated in calculateCF()
    public class GLSFit extends AbstractCurveFitter {

        public int numfreefitpar;

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess;

            int i = 0;				// add data 
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                weights[i] = point.getWeight();
                i += 1;
            }

            int numparfit = 0;		// count how many paramters will be fit
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    numparfit++;
                }
            }

            numfreefitpar = numparfit;

            initialGuess = new double[numparfit];	// setup the initial guesses for the parameters
            int num = 0;
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    initialGuess[num] = paraminitval[i];
                    num++;
                }
            }

            ParametricUnivariateFunction glsfitfunction;
            glsfitfunction = new GLS_fitFunction();

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(glsfitfunction, points);

            return new LeastSquaresBuilder().
                    maxEvaluations(fitMaxIterations).
                    maxIterations(fitMaxEvaluations).
                    start(initialGuess).
                    target(target).
                    weight(new DiagonalMatrix(weights)).
                    model(model.getModelFunction(), model.getModelFunctionJacobian()).
                    build();
        }

        public Map doFit(int cpx, int cpy, int kcf, String model) {
            //GLSFit fitter = new GLSFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();

            dataTransform(acf[kcf][cpx][cpy], currentCovmats);	// transform data
            transTheoreticalACF = new double[chanum - 1];
            transTheoreticalGradientACF = new double[noparam][chanum - 1];

            // Add points here
            for (int i = 1; i < chanum; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1, lagtime[i], transACF.getEntry(i - 1));
                points.add(point);
            }
            currentmodel = model;
            ParametricUnivariateFunction plotfunction;
            if (currentmodel.equals("FCS")) {	// select the fit model to be used; extra fit models can be added here
                plotfunction = new FCS_3p();
            } else {
                plotfunction = new FCCS_2p();
            }

            Map<String, Object> map = new HashMap<>();
            try {
                final Optimum topt = getOptimizer().optimize(getProblem(points));
                double result[] = topt.getPoint().toArray();
                double tmpres[] = topt.getResiduals().toArray();
                double[] tres = new double[chanum - 1];

                // store results and create fit and residual function for the plot
                for (int i = 0; i < chanum; i++) {
                    fitacf[kcf][cpx][cpy][i] = plotfunction.value(lagtime[i], result);			// calculate the fit function
                    if (i == 0) {
                        res[kcf][cpx][cpy][i] = acf[kcf][cpx][cpy][i] - fitacf[kcf][cpx][cpy][i];	// calculate the residuals
                    } else {
                        res[kcf][cpx][cpy][i] = tmpres[i - 1];
                        tres[i - 1] = acf[kcf][cpx][cpy][i] - fitacf[kcf][cpx][cpy][i];			// unweighted residuals for model probability calculations
                    }
                }

                chi2[kcf][cpx][cpy] = 0; // initialize chi2 value for this pixel
                for (int i = 1; i < chanum; i++) {
                    chi2[kcf][cpx][cpy] += Math.pow(res[kcf][cpx][cpy][i], 2) / (chanum - numfreefitpar - 1);	// calculate chi2 value; do not include the 0 lagtime channel which contains shot nosie
                }

                int num = 0;							// store the fit results 
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        fitres[kcf][cpx][cpy][i] = result[num];		// for free parameters use the fit results
                        num++;
                    } else {
                        fitres[kcf][cpx][cpy][i] = paraminitval[i];	// for fixed parameters use the initial values
                    }
                }
                pixfitted[kcf][cpx][cpy] = true;

                // update parameters in fit window
                if (tbShow.isSelected()) {
                    tfParamN.setText(IJ.d2s(fitres[kcf][cpx][cpy][0], decformat));
                    tfParamD.setText(IJ.d2s(fitres[kcf][cpx][cpy][1] * Math.pow(10, 12), decformat));
                    tfParamVx.setText(IJ.d2s(fitres[kcf][cpx][cpy][2] * Math.pow(10, 6), decformat));
                    tfParamVy.setText(IJ.d2s(fitres[kcf][cpx][cpy][3] * Math.pow(10, 6), decformat));
                    tfParamG.setText(IJ.d2s(fitres[kcf][cpx][cpy][4], decformat2));
                    tfParamF2.setText(IJ.d2s(fitres[kcf][cpx][cpy][5], decformat));
                    tfParamD2.setText(IJ.d2s(fitres[kcf][cpx][cpy][6] * Math.pow(10, 12), decformat));
                    tfParamF3.setText(IJ.d2s(fitres[kcf][cpx][cpy][7], decformat));
                    tfParamD3.setText(IJ.d2s(fitres[kcf][cpx][cpy][8] * Math.pow(10, 12), decformat));
                    tfParamFtrip.setText(IJ.d2s(fitres[kcf][cpx][cpy][9], decformat));
                    tfParamTtrip.setText(IJ.d2s(fitres[kcf][cpx][cpy][10] * Math.pow(10, 6), decformat));
                }

                map.put("results", result);
                map.put("covariance", topt.getCovariances(1).getData());
                map.put("residuals", tres);
                map.put("sigma", topt.getSigma(1).toArray());

            } catch (Exception e) {
                IJ.log(e.getClass().getName() + " at pixel " + Integer.toString(cpy) + " - " + Integer.toString(cpx));
                ArrayList<Double> result = new ArrayList<>();
                double[] tres = new double[chanum - 1];
                Arrays.fill(tres, 1);
                double[] tsig = new double[chanum - 1];
                Arrays.fill(tsig, 1);
                int num = 0;
                /*for (int i = 0; i < noparam; i++) { 
					if ( paramfit[i] == true ) {
						result.add(initparam[i]);		// return the initial values
						num++;
					}
				}
                 */

                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        fitres[kcf][cpx][cpy][i] = Double.NaN;
                        result.add(initparam[i]);		// return the initial values
                        num++;
                    }
                }
                for (int i = 0; i < chanum; i++) {
                    fitacf[kcf][cpx][cpy][i] = 0.0;
                    res[kcf][cpx][cpy][i] = 0.0;	// calculate the residuals

                }
                double[][] tcov = new double[num][num];
                for (double[] tcov1 : tcov) {
                    Arrays.fill(tcov1, 1);
                }
                map.put("results", result.toArray());
                map.put("covariance", tcov);
                map.put("residuals", tres);
                map.put("sigma", tsig);
                pixfitted[kcf][cpx][cpy] = false;
            }

            return map;
        }
    }

    // Linear fit
    public class LineFit extends AbstractCurveFitter {	// simple line fit mainly used for diffusion law fit

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] xtarget = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess = new double[2];

            int i = 0;
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                xtarget[i] = point.getX();
                weights[i] = point.getWeight();
                i += 1;
            }

            // initial guesses
            initialGuess[0] = target[0];						// use first point as intercept estimate
            initialGuess[1] = (target[1] - target[0]) / (xtarget[1] - xtarget[0]);	// use slope calculated from first two points as slope estimate		

            ParametricUnivariateFunction function;
            function = new Line();

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(function, points);

            return new LeastSquaresBuilder().
                    maxEvaluations(Integer.MAX_VALUE).
                    maxIterations(Integer.MAX_VALUE).
                    start(initialGuess).
                    target(target).
                    weight(new DiagonalMatrix(weights)).
                    model(model.getModelFunction(), model.getModelFunctionJacobian()).
                    build();
        }

        public double[] doFit(double[] xtrace, double[] trace, double[] sdtrace, int num) {
            LineFit fitter = new LineFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();

            // Add points here
            for (int i = 0; i < num; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1 / sdtrace[i] / sdtrace[i], xtrace[i], trace[i]);
                points.add(point);
            }

            double result[] = fitter.fit(points);
            return (result);
        }
    }

    // Single exponential fit for bleach correction
    public class SingleExpFit extends AbstractCurveFitter {

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess = new double[3];

            int i = 0;
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                weights[i] = point.getWeight();
                i += 1;
            }

            // initial guesses
            initialGuess[0] = target[1];							// use first point as intercept estimate
            initialGuess[1] = intTime[(int) Math.floor(len / 2)];	// use half the intensity trace time as first estimate for the exponential decay	
            initialGuess[2] = 0;

            ParametricUnivariateFunction function;
            function = new SingleExp();

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(function, points);

            return new LeastSquaresBuilder().
                    maxEvaluations(Integer.MAX_VALUE).
                    maxIterations(Integer.MAX_VALUE).
                    start(initialGuess).
                    target(target).
                    weight(new DiagonalMatrix(weights)).
                    model(model.getModelFunction(), model.getModelFunctionJacobian()).
                    build();
        }

        public double[] doFit(double[] itrace) {
            SingleExpFit fitter = new SingleExpFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();
            int num = itrace.length;

            // Add points here
            for (int i = 0; i < num; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1, intTime[i], itrace[i]);
                points.add(point);
            }

            double result[] = fitter.fit(points);
            return (result);
        }
    }

    // Double exponential fit for bleach correction
    public class DoubleExpFit extends AbstractCurveFitter {

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess = new double[5];

            int i = 0;
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                weights[i] = point.getWeight();
                i += 1;
            }

            // initial guesses
            initialGuess[0] = target[1] / 2;							// amplitude for first and second decay are set equal; estimated from half of the first point
            initialGuess[1] = intTime[(int) Math.floor(len / 10)];	// use a tenth of the intensity trace time as first estimate for the first exponential decay
            initialGuess[2] = target[1] / 2;							// amplitude estimated from half of the first point
            initialGuess[3] = intTime[(int) Math.floor(len / 2)];						// use half the intensity trace time as first estimate for the second exponential decay
            initialGuess[4] = 0;

            ParametricUnivariateFunction function;
            function = new DoubleExp();

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(function, points);

            return new LeastSquaresBuilder().
                    maxEvaluations(Integer.MAX_VALUE).
                    maxIterations(Integer.MAX_VALUE).
                    start(initialGuess).
                    target(target).
                    weight(new DiagonalMatrix(weights)).
                    model(model.getModelFunction(), model.getModelFunctionJacobian()).
                    build();
        }

        public double[] doFit(double[] itrace) {
            DoubleExpFit fitter = new DoubleExpFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();
            int num = itrace.length;

            // Add points here
            for (int i = 0; i < num; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1, intTime[i], itrace[i]);
                points.add(point);
            }

            double result[] = fitter.fit(points);
            return (result);
        }
    }

    // Polynomial fit for bleach correction
    public class PolynomFit extends AbstractCurveFitter {

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess = new double[polyOrder + 1];

            int i = 0;
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                weights[i] = point.getWeight();
                i += 1;
            }

            // initial guesses
            initialGuess[0] = target[len - 1];							// use the last point as offset estimate
            for (int j = 1; j <= polyOrder; j++) {						// use a straight line as the first estimate
                initialGuess[j] = 0;
            }

            ParametricUnivariateFunction function;
            function = new Polynomial();

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(function, points);

            return new LeastSquaresBuilder().
                    maxEvaluations(Integer.MAX_VALUE).
                    maxIterations(Integer.MAX_VALUE).
                    start(initialGuess).
                    target(target).
                    weight(new DiagonalMatrix(weights)).
                    model(model.getModelFunction(), model.getModelFunctionJacobian()).
                    build();
        }

        public double[] doFit(double[] itrace) {
            PolynomFit fitter = new PolynomFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();
            int num = itrace.length;

            // Add points here
            for (int i = 0; i < num; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1, intTime[i], itrace[i]);
                points.add(point);
            }

            double result[] = fitter.fit(points);
            return (result);
        }
    }

    // Diffusion law fit
    public double[][][][] dlfit(int numbin, int numpsf) {					// perform calculations of average diffusion coefficients for different binning
        double[][][][] results;												// this is used by the PSF and the Diffusion Law calculations by using a range
        double[][][] fitresMem = new double[width][height][noparam + 3];	// an array to remember the values in fitres[0][x][y][q], Chi2, blocking success and filtering mask and return them there after the calculation is finished
        boolean[][] pixfitMem = new boolean[width][height];				// an array for temporary storage of pixfitted[0][x][y]
        double[][][][] acfMem = new double[4][width][height][chanum];		// an array to remember the values in acf[0][x][y][c], fited acf, residuals and standard deviation and return them there after the calculation is finished
        int bxmem = binningX;
        int bymem = binningY;
        int bpxmem = pixbinX;
        int bpymem = pixbinY;
        int maxs = 0;														// of binnings and caluclating the average D each time from fits to all pixels
        int max = 0;
        int count;
        int itw;
        int ith;
        int ROIdim = Integer.parseInt(tfDLROI.getText());

        doMSD = false;									// do not calculate msd when calculating diffusion law

        // creating a copy of fitres[0][x][y][q] and related arrays
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                System.arraycopy(fitres[0][x][y], 0, fitresMem[x][y], 0, noparam);
                fitresMem[x][y][noparam] = chi2[0][x][y];
                fitresMem[x][y][noparam + 1] = blocked[0][x][y];
                fitresMem[x][y][noparam + 2] = pixvalid[0][x][y];
                pixfitMem[x][y] = pixfitted[0][x][y];
            }
        }

        // creating a copy of acf[0][x][y][c] and related arrays
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int c = 0; c < chanum; c++) {
                    acfMem[0][x][y][c] = acf[0][x][y][c];
                    acfMem[1][x][y][c] = sdacf[0][x][y][c];
                    acfMem[2][x][y][c] = fitacf[0][x][y][c];
                    acfMem[3][x][y][c] = res[0][x][y][c];
                }
            }
        }

        if ("All".equals(tbDLROI.getText()) || numpsf >= 0) { 			// if All is selected or if PSF is to be calculated
            itw = 1;
            ith = 1;
            results = new double[1][1][3][numbin];
            if (overlap) {													// counting in overlap mode is different from non-overlap mode as different number of pixels will be calculated
                for (int x = 1; x <= numbin; x++) {
                    maxs += width - (x - 1);
                }
            } else {
                for (int x = 1; x <= numbin; x++) {
                    maxs += (int) Math.floor(width / x);
                }
            }
            if (numpsf >= 0) { 												// counting in status bar is different for diffusion law and PSF calculation
                max = maxs * (int) Math.floor((psfEnd - psfStart) / psfStep + 1.1);
                count = numpsf * maxs;
            } else {
                max = maxs;
                count = 0;
            }
        } else {																// if ROI is selected
            diffLawMapwidth = (int) Math.floor(width / ROIdim);
            diffLawMapheight = (int) Math.floor(height / ROIdim);
            itw = diffLawMapwidth;
            ith = diffLawMapheight;
            results = new double[diffLawMapwidth][diffLawMapheight][3][numbin];
            diffLawFitMap = new double[diffLawMapwidth][diffLawMapheight][2];		// define the array for diff law map
            for (int n = 1; n < numbin; n++) {
                maxs += n;
            }
            max = itw * ith * (ROIdim * numbin - maxs);
            count = 0;
        }

        int u_counter = 1; // for updating progress bar when in GPU mode
        for (int u = numbin; u >= 1; u--) {									// run for different bin sizes					
            binning = u;														// the value of binning (global variable) influences also calcIntensity and correlate; so it has to be correctly set here
            binningX = binning;
            binningY = binning;
            pixeldimx = (pixelsize * 1000 / objmag * binning) / Math.pow(10, 9);
            pixeldimy = (pixelsize * 1000 / objmag * binning) / Math.pow(10, 9);

            if (overlap || "ROI".equals(tbDLROI.getText())) {				// set pixel number according to whether overlap is allowed; in the case of diffusion law calcualtions
                pixelWidthX = width - binning;								// for sub-areas overalp is automatically allowed
                pixelHeightY = height - binning;
                pixbin = 1;
            } else {
                pixelWidthX = (int) Math.floor(width / binning) - 1;
                pixelHeightY = (int) Math.floor(height / binning) - 1;
                pixbin = binning;
            }
            pixbinX = pixbin;
            pixbinY = pixbin;

            maxcposx = pixelWidthX;						// define the maxposition
            mincposx = 0;
            maxcposy = pixelHeightY;
            mincposy = 0;

            initparam[0] = Double.parseDouble(tfParamN.getText());
            initparam[1] = Double.parseDouble(tfParamD.getText()) / Math.pow(10, 12);
            initparam[2] = Double.parseDouble(tfParamVx.getText()) / Math.pow(10, 6);
            initparam[3] = Double.parseDouble(tfParamVy.getText()) / Math.pow(10, 6);
            initparam[4] = Double.parseDouble(tfParamG.getText());
            initparam[5] = Double.parseDouble(tfParamF2.getText());
            initparam[6] = Double.parseDouble(tfParamD2.getText()) / Math.pow(10, 12);
            initparam[7] = Double.parseDouble(tfParamF3.getText());
            initparam[8] = Double.parseDouble(tfParamD3.getText()) / Math.pow(10, 12);
            initparam[9] = Double.parseDouble(tfParamFtrip.getText());
            initparam[10] = Double.parseDouble(tfParamTtrip.getText()) / Math.pow(10, 6);

            /*rbtnHoldN.setSelected(false);				// set hold values for parameters
			rbtnHoldD.setSelected(false);
			rbtnHoldVx.setSelected(true);
			rbtnHoldVy.setSelected(true);
			rbtnHoldG.setSelected(false);
			rbtnHoldF2.setSelected(true);
			rbtnHoldD2.setSelected(true);
			rbtnHoldF3.setSelected(true);
			rbtnHoldD3.setSelected(true);
			rbtnHoldFtrip.setSelected(true);
			rbtnHoldTtrip.setSelected(true);*/
            prepareFit();

            int[][] count_dlaw = new int[itw][ith];

            for (int r = 0; r < itw; r++) {
                for (int s = 0; s < ith; s++) {
                    results[r][s][0][u - 1] = obsvolFCS_ST2D1p(2) * Math.pow(10, 12);			// observation area in um^2

                    if ("ROI".equals(tbDLROI.getText())) {
                        maxcposx = (r + 1) * ROIdim - binning;
                        mincposx = ROIdim * r;
                        maxcposy = (s + 1) * ROIdim - binning;
                        mincposy = ROIdim * s;
                    }

                    boolean IsGPUCalculationOK = false; // Runs the calculation in CPU if error is encountered in GPU.
                    if (isgpupresent == 1) {
                        boolean temp_doFit = doFit;
                        doFit = true;
                        isdlawcalculatedingpu = 1;

                        roi1StartX = (cfXDistance > 0) ? 0 : -cfXDistance;
                        roi1StartY = (cfYDistance > 0) ? 0 : -cfYDistance;

                        if (overlap) {
                            roi1WidthX = width - Math.abs(cfXDistance);
                            roi1HeightY = height - Math.abs(cfYDistance);
                        } else {
                            roi1WidthX = (int) Math.floor((width - Math.abs(cfXDistance)) / binningX) * binningX;
                            roi1HeightY = (int) Math.floor((height - Math.abs(cfYDistance)) / binningY) * binningY;
                        }

                        Roi impRoi1 = new Roi(roi1StartX, roi1StartY, roi1WidthX, roi1HeightY);

                        IsGPUCalculationOK = GPU_Calculate_ACF_All(impRoi1);
                        doFit = temp_doFit;
                    }

                    for (int x = mincposx; x <= maxcposx; x++) {
                        for (int y = mincposy; y <= maxcposy; y++) {
                            if (isgpupresent != 1 || !IsGPUCalculationOK) {
                                calcIntensityTrace(imp, x * pixbin, y * pixbin, x * pixbin, y * pixbin, firstframe, lastframe);
                                correlate(imp, x * pixbin, y * pixbin, x * pixbin, y * pixbin, 0, firstframe, lastframe);	// correlate the pixel
                                fit(x, y, 0, "FCS");										// fit the correlation function with a simple 1 component fit
                            }
                            if (pixfitted[0][x][y]) {						// take only data that has been fitted
                                results[r][s][1][u - 1] += fitres[0][x][y][1];					// sum all results
                                results[r][s][2][u - 1] += Math.pow(fitres[0][x][y][1], 2);
                                count_dlaw[r][s]++;
                            }
                            // and their square to calcualte averages and standard deviations
                        }
                        if (isgpupresent != 1 || !IsGPUCalculationOK) {
                            IJ.showProgress(count++, max);
                        }
                    }

                    results[r][s][1][u - 1] *= Math.pow(10, 12) / count_dlaw[r][s];	// average D in um2/s
                    results[r][s][2][u - 1] *= Math.pow(10, 24) / count_dlaw[r][s];	// normalize the average of the square
                    results[r][s][2][u - 1] = results[r][s][2][u - 1] - Math.pow(results[r][s][1][u - 1], 2);
                }
            }

            if (isgpupresent == 1) {
                int pcnt = (u_counter == numbin) ? 100 : (int) Math.floor(u_counter * 100 / numbin);
                IJ.showProgress(pcnt, 100);
                u_counter += 1;
            }
        }

        // reverting changes to fitres[0][x][y][q] and related arrays done by dlfit
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                System.arraycopy(fitresMem[x][y], 0, fitres[0][x][y], 0, noparam);
                chi2[0][x][y] = fitresMem[x][y][noparam];
                blocked[0][x][y] = fitresMem[x][y][noparam + 1];
                pixvalid[0][x][y] = fitresMem[x][y][noparam + 2];
                pixfitted[0][x][y] = pixfitMem[x][y];
            }
        }

        // reverting changes to acf[0][x][y][c] and related arrays done by dlfit
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int c = 0; c < chanum; c++) {
                    acf[0][x][y][c] = acfMem[0][x][y][c];
                    sdacf[0][x][y][c] = acfMem[1][x][y][c];
                    fitacf[0][x][y][c] = acfMem[2][x][y][c];
                    res[0][x][y][c] = acfMem[3][x][y][c];
                }
            }
        }

        binningX = bxmem;
        binningY = bymem;
        pixbinX = bpxmem;
        pixbinY = bpymem;

        doMSD = tbMSD.isSelected();		// set MSD selection again to selected value 
        isdlawcalculatedingpu = 0;
        return results;
    }

    // calculate mean square displacement from ITIR-FCS correlation function
    public double[] correlationToMSD(double[] corrfunc, double psize, double psfwidth) {
        // corrfunc: correlation fucntion to be inverted
        // pszie: pixel size in um
        // psfwidth: pixelwidth in um
        double pi = 3.14159265359;
        double a = -17.0 * 1260.0 / 29.0 / 180.0;
        double b = 1260.0 / 29.0 / 3.0;
        double c = -1260.0 / 29.0;
        double d;
        int cutoffInd = chanum - 1;
        double[] msdarray = new double[chanum]; // array is defined as full length even if it is not completely used

        for (int i = chanum - 1; i > 1; i--) {
            if (corrfunc[i] / corrfunc[1] < 0.1) {
                cutoffInd = i;
            }
        }

        for (int i = 1; i < cutoffInd; i++) {
            d = pi * corrfunc[i] / corrfunc[1] / (obsvolFCS_ST2D1p(2) * Math.pow(10, 12)) * Math.pow(psize, 2) * 1260.0 / 29.0;
            double[] g = quarticFunction(a, b, c, d);
            msdarray[i] = (psize * psize / (g[1]) - psfwidth * psfwidth);
        }

        for (int i = 0; i < msdarray.length; i++) {
            if (Double.isNaN(msdarray[i])) {
                msdarray[i] = 0;
            }
        }
        return msdarray;
    }

    // calculate mean square displacement from SPIM correlation function
    public double[] correlationToMSD3D(double[] corrfunc, double psize, double psfwidth, double psfwidthz) {
        // corrfunc: correlation function to be inverted
        // psize: pixel size in um
        // psfwidth: psf width in um
        // psfwidthz: z psf width in um	

        int cutoffInd = chanum - 1;
        double[] msdarray = new double[chanum]; // array is defined as full length even if it is not completely used

        msdarray[0] = 0.0;
        msdarray[1] = 0.0;

        for (int i = chanum - 1; i > 1; i--) {
            if (corrfunc[i] / corrfunc[1] < 0.1) {
                cutoffInd = i;
            }
        }

        final double temp1_msd_3d = corrfunc[1];
        final double temp_psfwidth = psfwidth;
        final double temp_psfwidthz = psfwidthz;
        final double temp_psize = psize;

        for (int i = 2; i < cutoffInd; i++) {
            final double temp_msd_3d = corrfunc[i];
            UnivariateFunction f = (double x) -> {
                double pi = 3.14159265359;
                double sqrpi = Math.sqrt(pi);

                double p = 1 - Math.pow(temp_psfwidth, 2.0) / Math.pow(temp_psfwidthz, 2.0);
                double q = Math.pow(temp_psize, 2.0) / Math.pow(temp_psfwidthz, 2.0);
                double p1 = -15.0 * Math.pow(p, 2.0) + 20.0 * p * q + 2.0 * Math.pow(q, 2.0);
                double p2 = 945.0 * Math.pow(p, 3.0) - 630.0 * Math.pow(p, 2.0) * q;

                double coef_0 = -1.0;
                double coef_2 = 1 / 63.0 / q * (p2 - 420.0 * p * Math.pow(q, 2.0) - 64.0 * Math.pow(q, 3.0)) / p1;
                double coef_3 = -1 / pi / Math.sqrt(q);
                double coef_4 = 1 / 7560.0 / Math.pow(q, 2.0) / p1 * (14175.0 * Math.pow(p, 4.0) + 37800.0 * Math.pow(p, 3.0) * q - 30240.0 * Math.pow(p, 2.0) * Math.pow(q, 2.0) - 3840.0 * p * Math.pow(q, 3.0) - 1132.0 * Math.pow(q, 4.0));
                double coef_5 = 1 / 126.0 / pi / Math.pow(q, 1.5) * (p2 + 126.0 * p * Math.pow(q, 2.0) - 44.0 * Math.pow(q, 3.0)) / p1;

                double corrp = temp_msd_3d * sqrpi * temp_psize * temp_psize * temp_psfwidthz / (obsvolFCS_ST2D1p(3) * Math.sqrt(pi) * temp_psfwidthz * Math.pow(10, 12)) / temp1_msd_3d;
                return coef_5 * Math.pow(x, 5.0) - coef_4 * corrp * Math.pow(x, 4.0) + coef_3 * Math.pow(x, 3.0) - coef_2 * corrp * Math.pow(x, 2.0) - coef_0 * corrp;
            };

            try {
                UnivariateSolver solver = new BrentSolver();
                double result = solver.solve(100, f, 0, psize / psfwidth);
                msdarray[i] = 1.5 * (Math.pow(psize, 2.0) / Math.pow(result, 2.0) - Math.pow(psfwidth, 2.0));
            } catch (RuntimeException ex) {
                IJ.log(ex.getMessage());
                throw ex;
            }

        }

        for (int i = 0; i < msdarray.length; i++) {
            if (Double.isNaN(msdarray[i])) {
                msdarray[i] = 0;
            }
        }

        return msdarray;
    }

    //This calculates the solution to a third order polynomial which is obtained from the fourth order polynomial. The cubic which appears is called resolvent cubic equation
    public double cubicrootFunction(double a, double b, double c) {
        double pi = 3.14159265359;
        double root = 0.0;
        double qq, pp;
        double Qa, y1 = 0.0, y2, y3 = 0.0;
        double AA = 0.0, BB = 0.0, sQ = 0.0;
        double cosA, alpha;

        pp = b - (a * a) / 3.0;
        qq = c + 2.0 * (a / 3.0) * (a / 3.0) * (a / 3.0) - a * b / 3.0;
        Qa = (pp / 3.0) * (pp / 3.0) * (pp / 3.0) + (qq / 2.0) * (qq / 2.0);
        a = a / 3.0;
        if (Qa >= 0.0) {
            sQ = Math.sqrt(Qa);
            if (sQ > qq / 2.0) {
                AA = Math.pow(sQ - qq / 2.0, 0.33333333333333333);
            } else {
                AA = -Math.pow(qq / 2.0 - sQ, 0.33333333333333333);
            }
            if (-sQ - qq / 2.0 > 0.0) {
                BB = Math.pow(-sQ - qq / 2.0, 0.33333333333333333);
            } else {
                BB = -Math.pow(sQ + qq / 2.0, 0.33333333333333333);
            }
            y1 = AA + BB - a;
        } else {
            pp = pp / 3;
            cosA = -qq / (2 * Math.sqrt(-pp * pp * pp));
            alpha = Math.acos(cosA) / 3.0;
            if (alpha == 0) {
                alpha = 2.0 * pi / 3.0;
            }
            y1 = 2 * Math.sqrt(-pp) * Math.cos(alpha) - a;
        }
        return y1;
    }

    //Fourth order polynomial solver
    public double[] quarticFunction(double a, double b, double c, double d) {
        double AA = -b;
        double BB = a * c - 4.0 * d;
        double CC = -a * a * d + 4.0 * b * d - c * c;
        double y1 = cubicrootFunction(AA, BB, CC);
        double RR = Math.sqrt(a * a / 4.0 - b + y1);
        double DD, EE, x1, x2, x3, x4;

        if (RR == 0) {
            DD = Math.sqrt(3.0 * a * a / 4.0 - 2.0 * b + 2.0 * Math.sqrt(y1 * y1 - 4.0 * d));
            EE = Math.sqrt(3.0 * a * a / 4.0 - 2.0 * b - 2.0 * Math.sqrt(y1 * y1 - 4.0 * d));
        } else {
            DD = Math.sqrt(3.0 * a * a / 4.0 - RR * RR - 2.0 * b + (4.0 * a * b - 8.0 * c - a * a * a) / (4.0 * RR));
            EE = Math.sqrt(3.0 * a * a / 4.0 - RR * RR - 2.0 * b - (4.0 * a * b - 8.0 * c - a * a * a) / (4.0 * RR));
        }

        a = a / 4.0;
        RR = RR / 2.0;
        DD = DD / 2.0;
        x1 = -a + RR + DD;
        x2 = -a + RR - DD;
        EE = EE / 2.0;
        x3 = -a - RR + EE;
        x4 = -a - RR - EE;
        double[] returningarray = {x1, x2, x3, x4};
        return returningarray;
    }

    /*  
 * Following are the output procedures for the program:
 * 
 * public void plotCF(Roi plotroi, int cormode, boolean map): plots the correlation functions (no fit necessary)
 * public void plotTheoreticalCF(): calcualtes and plots a theoretical CF from the parameters in the fit window
 * public void plotIntensityTrace(Roi plotroi, int cormode): plots the intensity traces (no fit necessary)
 * public void plotResiduals(int rpx1, int rpy1): plots the fit residuals
 * public void plotSD(int sdpx1, int sdpy1, int cormode): plots the SD
 * public void plotMSD(Roi MSDroi, int cormode): plots the MSD
 * public void createParaImp(int wx, int hy): create the window for the parameter maps
 * public void impPara1Adjusted(): Listener for the Parameter window; if the map is changed the histograms are 
 * adapted to the presently shown maps
 * MouseListener para1MouseClicked: MouseListener which allows the user to click into the parameter window and 
 * obtain its values in the Fit Panel
 * public void updateParaImp(): update parameter maps if new fits or filtering were performed 
 * public boolean filterPix (int m, int x, int y): decides whether a pixel is valid - whether the fot parameters fall within the thresholds
 * 
     */
    public void plotCF(Roi plotroi, int cormode, boolean map) {
        // plotroi: roi for which CFs are to be plotted
        // cormode = 1: plot a single function; correlation calculated between (cpx1, cpy1) to (cpx2, cpy2) 
        // cormode = 2: plot all CFs in a ROI; note that this can be ACFs or CCFs depending on the setting
        // of cfXDistance and cfYDistance
        // cormode = 3, plot ACF and CCFs for FCCS
        // map: is roi selected in parameter window (reproduction of values only) or is it from a new calculation?
        // cormode = 4, plot average ACF
        if (!plotACFCurves) {
            return;
        }

        parcormode = cormode;

        int ct = 0;
        int cpx1;
        int cpy1;
        int cpxf;
        int cpyf;
        int cpx2;
        int cpy2;

        if (cbFitModel.getSelectedItem() == "DC-FCCS" && cormode == 1) { 	// index ct ensures that in cormode 1 when DC-FCCS is selected the CCF is plotted
            ct = 2;
        }

        Rectangle plotrect = plotroi.getBounds();
        if (map == true) {										// if the ROI is selected in the parameter map
            cpx1 = (int) plotrect.getX();
            cpy1 = (int) plotrect.getY();
            if (cfXDistance < 0) {
                cpx1 = cpx1 - (int) Math.floor(cfXDistance / pixbinX);
            }
            if (cfYDistance < 0) {
                cpy1 = cpy1 - (int) Math.floor(cfYDistance / pixbinY);
            }
            cpxf = (int) (cpx1 + plotrect.getWidth());
            cpyf = (int) (cpy1 + plotrect.getHeight());
            cpx2 = (cpx1 * pixbinX + cfXDistance);
            cpy2 = (cpy1 * pixbinY + cfYDistance);
        } else {
            cpx1 = (int) Math.ceil(plotrect.getX() / pixbinX);
            cpy1 = (int) Math.ceil(plotrect.getY() / pixbinY);
            cpxf = (int) Math.floor((plotrect.getX() + plotrect.getWidth() - binningX) / pixbinX);
            cpyf = (int) Math.floor((plotrect.getY() + plotrect.getHeight() - binningY) / pixbinY);
            cpx2 = (cpx1 * pixbinX + cfXDistance);
            cpy2 = (cpy1 * pixbinY + cfYDistance);
        }

        maxsc = acf[ct][cpx1][cpy1][1];		// minimum and maximum setting for plot window
        minsc = 1000;						// set minsc to a high value to make sure it is not 0

        if (cormode == 1) {					// if single ACF or CCF
            for (int z = 1; z <= (chanum - 1); z++) {	// find minimum and maximum values in correlation functions that will be plotted
                if (maxsc < acf[ct][cpx1][cpy1][z]) {
                    maxsc = acf[ct][cpx1][cpy1][z];
                }
                if (minsc > acf[ct][cpx1][cpy1][z]) {
                    minsc = acf[ct][cpx1][cpy1][z];
                }
            }

            // maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            minsc -= minsc * 0.1;
            maxsc += maxsc * 0.1;

            // plot
            Plot plot = new Plot($acfWindowTitle, "tau [s]", "G(tau)", lagtime, acf[ct][cpx1][cpy1]);
            plot.setFrameSize(acfWindowDimX, acfWindowDimY);
            plot.setLogScaleX();
            plot.setLimits(lagtime[1], 2 * lagtime[chanum - 1], minsc, maxsc);
            plot.setColor(java.awt.Color.BLUE);
            plot.setJustification(Plot.CENTER);

            // create plot label for ACF of CCF
            if (cfXDistance + cfYDistance != 0) {
                plot.addLabel(0.5, 0, " CCF of (" + cpx1 * pixbinX + ", " + cpy1 * pixbinY + ")" + " and (" + cpx2 + ", " + cpy2 + ")" + " at (" + binningX + "x" + binningY + "binning");
            } else if (cfXDistance + cfYDistance != 0 && cbFitModel.getSelectedItem() == "DC-FCCS") {
                plot.addLabel(0.5, 0, " CCF of (" + cpx1 * pixbinX + ", " + cpy1 * pixbinY + ")" + " and (" + cpx2 + ", " + cpy2 + ")" + " at " + binningX + "x" + binningY + "binning");
            } else {
                plot.addLabel(0.5, 0, " ACF of (" + cpx1 * pixbinX + ", " + cpy1 * pixbinY + ") at " + binningX + "x" + binningY + "binning");
            }
            plot.draw();
            int num1 = 0;
            if (isgpupresent == 1) {
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i]) {
                        num1++;
                    }
                }
            }
            if (isgpupresent == 1) {
                ParametricUnivariateFunction theofunction = new FCS_3p();
                double[] parameters = new double[num1];
                num1 = 0;
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i]) {
                        parameters[num1] = fitres[0][cpx1][cpy1][i];
                        num1++;
                    }
                    for (int j = 1; j < chanum; j++) {
                        fitacf[0][cpx1][cpy1][j] = theofunction.value(lagtime[j], parameters);
                    }
                }
            }

            // if fit has been performed add the fit to the plot
            if (doFit) {
                plot.setColor(java.awt.Color.RED);
                plot.addPoints(Arrays.copyOfRange(lagtime, fitstart, fitend + 1), Arrays.copyOfRange(fitacf[ct][cpx1][cpy1], fitstart, fitend + 1), Plot.LINE);
                //plot.setColor(java.awt.Color.BLUE);
                plot.draw();
                if (plotResCurves) {
                    plotResiduals(cpx1, cpy1); // and add a residual window
                }
            }

            // either create a new plot window or plot within the existing window
            if (acfWindow == null || acfWindow.isClosed() == true) {
                acfWindow = plot.show();
                acfWindow.setLocation(acfWindowPosX, acfWindowPosY);
            } else {
                acfWindow.drawPlot(plot);
                acfWindow.setTitle($acfWindowTitle);
            }

            // plot the standard deviation of the CF
            if (plotSDCurves) {
                plotSD(cpx1, cpy1, ct);
            }

        }

        if (cormode == 2) {					// if multiple ACF or CCF
            for (int x = cpx1; x <= cpxf; x++) {		// find minimum and maximum values in correlation functions that will be plotted
                for (int y = cpy1; y <= cpyf; y++) {
                    for (int z = 1; z <= (chanum - 1); z++) {
                        if ((!map && plotroi.contains(x * pixbinX, y * pixbinY) && plotroi.contains(x * pixbinX, y * pixbinY + binningY - 1) && plotroi.contains(x * pixbinX + binningX - 1, y * pixbinY) && plotroi.contains(x * pixbinX + binningX - 1, y * pixbinY + binningY - 1)) || (map == true && plotroi.contains(x, y))) {
                            if (maxsc < acf[0][x][y][z]) {
                                maxsc = acf[0][x][y][z];
                            }
                            if (acf[0][x][y][z] != 0 && minsc > acf[0][x][y][z]) { // make sure minsc is not set to 0 because of a missing CF
                                minsc = acf[0][x][y][z];
                            }
                        }
                    }
                }
            }

            // maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            minsc -= minsc * 0.1;
            maxsc += maxsc * 0.1;

            //create empty plot
            Plot plot = new Plot($acfWindowTitle, "tau [s]", "G(tau)", empty, empty);
            plot.setFrameSize(acfWindowDimX, acfWindowDimY);
            plot.setLogScaleX();
            plot.setLimits(lagtime[1], 2 * lagtime[chanum - 1], minsc, maxsc);
            plot.setColor(java.awt.Color.BLUE);
            plot.setJustification(Plot.CENTER);

            // If this is a DC-FCCS plot where all functions are supposed to be plotted then plot the cross-correlations only
            int cftype;
            if (cbFitModel.getSelectedItem() == "DC-FCCS") {
                cftype = 2;		// cross-correlations are stored in acf[2]; the autocorrelations in acf[0] and acf[1] are not shown in this view
            } else {
                cftype = 0;		// autocorrelations are stored in acf[0];
            }

            // create plot label for ACF of CCF
            if (cfXDistance + cfYDistance != 0 && cftype != 2) {
                plot.addLabel(0.5, 0, " CCFs of pixels in the ROIs at " + binningX + "x" + binningY + "binning");
            } else if (cfXDistance + cfYDistance != 0 && cftype == 2) {
                plot.addLabel(0.5, 0, " CCFs of pixels in the ROIs at " + binningX + "x" + binningY + "binning");
            } else {
                plot.addLabel(0.5, 0, " ACFs from the ROI at " + binningX + "x" + binningY + "binning");
            }

            // plot all CFs and if fit has been performed, add fits
            if (doFit) {
                for (int y = cpy1; y <= cpyf; y++) {
                    for (int x = cpx1; x <= cpxf; x++) {
                        if ((!map && plotroi.contains(x * pixbinX, y * pixbinY) && plotroi.contains(x * pixbinX, y * pixbinY + binningY - 1) && plotroi.contains(x * pixbinX + binningX - 1, y * pixbinY) && plotroi.contains(x * pixbinX + binningX - 1, y * pixbinY + binningY - 1)) || (map == true && plotroi.contains(x, y))) {
                            plot.setColor(java.awt.Color.BLUE);
                            plot.addPoints(lagtime, acf[cftype][x][y], Plot.LINE);
                            plot.setColor(java.awt.Color.RED);
                            int num1 = 0;
                            if (isgpupresent == 1) {
                                for (int i = 0; i < noparam; i++) {
                                    if (paramfit[i]) {
                                        num1++;
                                    }
                                }
                            }

                            if (isgpupresent == 1) {
                                double[] parameters = new double[num1];
                                ParametricUnivariateFunction theofunction;
                                theofunction = new FCS_3p();
                                num1 = 0;
                                for (int i = 0; i < noparam; i++) {
                                    if (paramfit[i]) {
                                        parameters[num1] = fitres[0][x][y][i];
                                        num1++;
                                    }
                                }

                                for (int i = 1; i < chanum; i++) {
                                    fitacf[0][x][y][i] = theofunction.value(lagtime[i], parameters);
                                }
                            }
                            plot.addPoints(Arrays.copyOfRange(lagtime, fitstart, fitend + 1), Arrays.copyOfRange(fitacf[cftype][x][y], fitstart, fitend + 1), Plot.LINE);
                        }
                    }
                }
            } else {
                for (int y = cpy1; y <= cpyf; y++) {
                    for (int x = cpx1; x <= cpxf; x++) {
                        if ((!map && plotroi.contains(x * pixbinX, y * pixbinY) && plotroi.contains(x * pixbinX, y * pixbinY + binningY - 1) && plotroi.contains(x * pixbinX + binningX - 1, y * pixbinY) && plotroi.contains(x * pixbinX + binningX - 1, y * pixbinY + binningY - 1)) || (map == true && plotroi.contains(x, y))) {
                            plot.setColor(java.awt.Color.BLUE);
                            plot.addPoints(lagtime, acf[cftype][x][y], Plot.LINE);
                        }
                    }
                }
            }

            // either create a new plot window or plot within the existing window
            if (acfWindow == null || acfWindow.isClosed() == true) {
                acfWindow = plot.show();
                acfWindow.setLocation(acfWindowPosX, acfWindowPosY);
            } else {
                acfWindow.drawPlot(plot);
                acfWindow.setTitle($acfWindowTitle);
            }
        }

        if (cormode == 3) {					// if FCCS plot with ACF1, ACF2 and CCF
            for (int z = 1; z <= (chanum - 1); z++) {	// find minimum and maximum values in correlation functions that will be plotted
                if (maxsc < acf[0][cpx1][cpy1][z]) {
                    maxsc = acf[0][cpx1][cpy1][z];
                }
                if (maxsc < acf[1][cpx1][cpy1][z]) {
                    maxsc = acf[1][cpx1][cpy1][z];
                }
                if (maxsc < acf[2][cpx1][cpy1][z]) {
                    maxsc = acf[2][cpx1][cpy1][z];
                }
                if (minsc > acf[0][cpx1][cpy1][z]) {
                    minsc = acf[0][cpx1][cpy1][z];
                }
                if (minsc > acf[1][cpx1][cpy1][z]) {
                    minsc = acf[1][cpx1][cpy1][z];
                }
                if (minsc > acf[2][cpx1][cpy1][z]) {
                    minsc = acf[2][cpx1][cpy1][z];
                }
            }

            // maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            minsc -= minsc * 0.1;
            maxsc += maxsc * 0.1;

            // plot
            Plot plot = new Plot($acfWindowTitle, "tau [s]", "G(tau)");
            plot.setFrameSize(acfWindowDimX, acfWindowDimY);
            plot.setLogScaleX();
            plot.setLimits(lagtime[1], 2 * lagtime[chanum - 1], minsc, maxsc);
            plot.setColor(java.awt.Color.GREEN);
            plot.addPoints(lagtime, acf[0][cpx1][cpy1], Plot.LINE);
            plot.setJustification(Plot.CENTER);
            plot.setColor(java.awt.Color.RED);
            plot.addPoints(lagtime, acf[1][cpx1][cpy1], Plot.LINE);
            plot.setColor(java.awt.Color.BLUE);
            plot.addPoints(lagtime, acf[2][cpx1][cpy1], Plot.LINE);
            plot.addLabel(0.5, 0, " ACF1, ACF2, CCF of (" + cpx1 * pixbinX + ", " + cpy1 * pixbinY + ")" + " and (" + cpx2 * pixbinX + ", " + cpy2 * pixbinY + ")" + " at (" + binningX + "x" + binningY + "binning");
            plot.draw();

            // if fit has been performed add the fit to the plot
            if (doFit) {
                plot.setColor(java.awt.Color.BLACK);
                plot.addPoints(Arrays.copyOfRange(lagtime, fitstart, fitend + 1), Arrays.copyOfRange(fitacf[0][cpx1][cpy1], fitstart, fitend + 1), Plot.LINE);
                plot.addPoints(Arrays.copyOfRange(lagtime, fitstart, fitend + 1), Arrays.copyOfRange(fitacf[1][cpx1][cpy1], fitstart, fitend + 1), Plot.LINE);
                plot.addPoints(Arrays.copyOfRange(lagtime, fitstart, fitend + 1), Arrays.copyOfRange(fitacf[2][cpx1][cpy1], fitstart, fitend + 1), Plot.LINE);
                plot.draw();
                plot.setColor(java.awt.Color.GREEN);
                if (plotResCurves) {
                    plotResiduals(cpx1, cpy1); // and add a residual window
                }
            }

            // either create a new plot window or plot within the existing window
            if (acfWindow == null || acfWindow.isClosed() == true) {
                acfWindow = plot.show();
                acfWindow.setLocation(acfWindowPosX, acfWindowPosY);
            } else {
                acfWindow.drawPlot(plot);
                acfWindow.setTitle($acfWindowTitle);
            }

            // plot the standard deviation of the CF
            if (plotSDCurves) {
                plotSD(cpx1, cpy1, 2);
            }
        }

        if (cormode == 4) {					// if average ACF
            maxsc = aveacf[1];		// minimum and maximum setting for plot window
            minsc = 1000;
            for (int z = 1; z <= (chanum - 1); z++) {	// find minimum and maximum values in correlation functions that will be plotted
                if (maxsc < aveacf[z]) {
                    maxsc = aveacf[z];
                }
                if (minsc > aveacf[z]) {
                    minsc = aveacf[z];
                }
            }

            // maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            minsc -= minsc * 0.1;
            maxsc += maxsc * 0.1;

            // plot
            Plot plot = new Plot($acfWindowTitle, "tau [s]", "G(tau)", lagtime, aveacf);
            plot.setFrameSize(acfWindowDimX, acfWindowDimY);
            plot.setLogScaleX();
            plot.setLimits(lagtime[1], 2 * lagtime[chanum - 1], minsc, maxsc);
            plot.setColor(java.awt.Color.BLUE);
            plot.setJustification(Plot.CENTER);

            // create plot label for ACF of CCF
            plot.addLabel(0.5, 0, " Average ACF ");

            plot.draw();

            // if fit has been performed add the fit to the plot
            if (doFit) {
                plot.setColor(java.awt.Color.RED);
                plot.addPoints(Arrays.copyOfRange(lagtime, fitstart, fitend + 1), Arrays.copyOfRange(fitaveacf, fitstart, fitend + 1), Plot.LINE);
                plot.setColor(java.awt.Color.BLUE);
                plot.draw();
                /*if (plotResCurves) {
					plotResiduals(cpx1, cpy1); // and add a residual window
				}*/
            }

            // either create a new plot window or plot within the existing window
            if (acfWindow == null || acfWindow.isClosed() == true) {
                acfWindow = plot.show();
                acfWindow.setLocation(acfWindowPosX, acfWindowPosY);
            } else {
                acfWindow.drawPlot(plot);
                acfWindow.setTitle($acfWindowTitle);
            }

            // plot the standard deviation of the CF
            /*if (plotSDCurves) {
				plotSD(cpx1, cpy1, ct);
			}*/
        }

    }

    // plot a theretical ACF usng the parameters in the Fit Panel
    public void plotTheoreticalCF() {
        if (acfWindow == null || acfWindow.isClosed() == true) {
            JOptionPane.showMessageDialog(null, "No correlation window open");
            return;
        }

        ParametricUnivariateFunction theofunction;
        ArrayList<Double> paralist = new ArrayList<>();
        int cpx1 = initposx;
        int cpy1 = initposy;
        int cpx2 = cpx1 * pixbinX + cfXDistance;
        int cpy2 = cpy1 * pixbinY + cfYDistance;
        if (tfFitModel.getText().equals("FCS")) {	// select the fit model to be used; extra fit models can be added here
            theofunction = new FCS_3p();
            if (paramfit[0] == true) {
                paralist.add(Double.parseDouble(tfParamN.getText()));
            }
            if (paramfit[1] == true) {
                paralist.add(Double.parseDouble(tfParamD.getText()) / Math.pow(10, 12));
            }
            if (paramfit[2] == true) {
                paralist.add(Double.parseDouble(tfParamVx.getText()) / Math.pow(10, 6));
            }
            if (paramfit[3] == true) {
                paralist.add(Double.parseDouble(tfParamVy.getText()) / Math.pow(10, 6));
            }
            if (paramfit[4] == true) {
                paralist.add(Double.parseDouble(tfParamG.getText()));
            }
            if (paramfit[5] == true) {
                paralist.add(Double.parseDouble(tfParamF2.getText()));
            }
            if (paramfit[6] == true) {
                paralist.add(Double.parseDouble(tfParamD2.getText()) / Math.pow(10, 12));
            }
            if (paramfit[7] == true) {
                paralist.add(Double.parseDouble(tfParamF3.getText()));
            }
            if (paramfit[8] == true) {
                paralist.add(Double.parseDouble(tfParamD3.getText()) / Math.pow(10, 12));
            }
            if (paramfit[9] == true) {
                paralist.add(Double.parseDouble(tfParamFtrip.getText()));
            }
            if (paramfit[10] == true) {
                paralist.add(Double.parseDouble(tfParamTtrip.getText()) / Math.pow(10, 6));
            }
        } else {
            theofunction = new FCCS_2p();
            if (paramfit[0] == true) {
                paralist.add(Double.parseDouble(tfParamN.getText()));
            }
            if (paramfit[1] == true) {
                paralist.add(Double.parseDouble(tfParamD.getText()) / Math.pow(10, 12));
            }
            paramfit[2] = false;
            paraminitval[2] = 0;
            paramfit[3] = false;
            paraminitval[3] = 0;
            if (paramfit[4] == true) {
                paralist.add(Double.parseDouble(tfParamG.getText()));
            }
            if (paramfit[5] == true) {
                paralist.add(Double.parseDouble(tfParamF2.getText()));
            }
            if (paramfit[6] == true) {
                paralist.add(Double.parseDouble(tfParamD2.getText()) / Math.pow(10, 12));
            }
            paramfit[7] = false;
            paraminitval[7] = 0;
            paramfit[8] = false;
            paraminitval[8] = 0;
            if (paramfit[9] == true) {
                paralist.add(Double.parseDouble(tfParamFtrip.getText()));
            }
            if (paramfit[10] == true) {
                paralist.add(Double.parseDouble(tfParamTtrip.getText()) / Math.pow(10, 6));
            }
        }

        double[] plotdat = new double[chanum];
        double[] parameters = new double[paralist.size()];
        for (int i = 0; i < paralist.size(); i++) {
            parameters[i] = paralist.toArray(new Double[paralist.size()])[i];
        }
        for (int i = 0; i < chanum; i++) {
            plotdat[i] = theofunction.value(lagtime[i], parameters);
        }

        int ct = 0;

        if (cbFitModel.getSelectedItem() == "DC-FCCS") { 	// index ct ensures that in cormode 1 when DC-FCCS is selected the CCF is plotted
            ct = 2;
        }

        maxsc = acf[ct][cpx1][cpy1][1];		// minimum and maximum setting for plot window
        minsc = 1000;				// set minsc to a high value to make sure it is not 0

        for (int z = 1; z <= (chanum - 1); z++) {	// find minimum and maximum values in correlation functions that will be plotted
            if (maxsc < acf[ct][cpx1][cpy1][z]) {
                maxsc = acf[ct][cpx1][cpy1][z];
            }
            if (maxsc < plotdat[z]) {
                maxsc = plotdat[z];
            }
            if (minsc > acf[ct][cpx1][cpy1][z]) {
                minsc = acf[ct][cpx1][cpy1][z];
            }
            if (minsc > plotdat[z]) {
                minsc = plotdat[z];
            }
        }

        // maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
        minsc -= minsc * 0.1;
        maxsc += maxsc * 0.1;

        // plot
        Plot plot = new Plot($acfWindowTitle, "tau [s]", "G(tau)", lagtime, acf[ct][cpx1][cpy1]);
        plot.setFrameSize(acfWindowDimX, acfWindowDimY);
        plot.setLogScaleX();
        plot.setLimits(lagtime[1], 2 * lagtime[chanum - 1], minsc, maxsc);
        plot.setColor(java.awt.Color.BLUE);
        plot.setJustification(Plot.CENTER);
        plot.draw();

        // create plot label for ACF of CCF
        if (cfXDistance + cfYDistance != 0) {
            plot.addLabel(0.5, 0, " CCF of (" + cpx1 * pixbinX + ", " + cpy1 * pixbinY + ")" + " and (" + cpx2 + ", " + cpy2 + ")" + " at (" + binningX + "x" + binningY + "binning");
        } else if (cfXDistance + cfYDistance != 0 && cbFitModel.getSelectedItem() == "DC-FCCS") {
            plot.addLabel(0.5, 0, " CCF of (" + cpx1 * pixbinX + ", " + cpy1 * pixbinY + ")" + " and (" + cpx2 + ", " + cpy2 + ")" + " at " + binningX + "x" + binningY + "binning");
        } else {
            plot.addLabel(0.5, 0, " ACF of (" + cpx1 * pixbinX + ", " + cpy1 * pixbinY + ") at " + binningX + "x" + binningY + "binning");
        }

        plot.setColor(java.awt.Color.RED);
        plot.addPoints(lagtime, plotdat, Plot.LINE);
        plot.setColor(java.awt.Color.BLUE);
        plot.draw();

        acfWindow.drawPlot(plot);
        acfWindow.setTitle($acfWindowTitle);
    }

    // plot the intensity traces
    public void plotIntensityTrace(Roi plotroi, int cormode) {
        // plotroi: roi for which Intensities are to be plotted
        // cormode = 1: plot all intensity traces for a ROI for ACF 
        // cormode = 2: plot all intensity traces in a ROI; ACFs or CCFs depending on the setting of cfXDistance and cfYDistance
        // cormode = 3, plot all intensity traces for a ROI for ACF and CCFs for FCCS

        if (!plotIntensityCurves) {
            return;
        }

        Rectangle plotrect = plotroi.getBounds();
        int ipx1 = (int) Math.ceil(plotrect.getX() / pixbinX) * pixbinX;
        int ipy1 = (int) Math.ceil(plotrect.getY() / pixbinY) * pixbinY;
        int ipx2 = (ipx1 + cfXDistance);
        int ipy2 = (ipy1 + cfYDistance);

        iminsc = intTrace1[1];		// minimum and maximum setting for plot window
        imaxsc = intTrace1[1];

        if (cormode == 1) {
            for (int x = 0; x < nopit; x++) {	// find minimum and maximum values in intensity trace 1 that will be plotted
                if (imaxsc < intTrace1[x]) {
                    imaxsc = intTrace1[x];
                }
                if (iminsc > intTrace1[x]) {
                    iminsc = intTrace1[x];
                }
            }

            if ((ipx1 - ipx2) != 0 || (ipy1 - ipy2) != 0) { // in the case of a CCF find also the max and min value in the intensity trace 2
                for (int x = 0; x < nopit; x++) {	// find minimum and maximum values in intensity trace 1 that will be plotted
                    if (imaxsc < intTrace2[x]) {
                        imaxsc = intTrace2[x];
                    }
                    if (iminsc > intTrace2[x]) {
                        iminsc = intTrace2[x];
                    }
                }
            }

            iminsc -= iminsc * 0.1;		// maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            imaxsc += imaxsc * 0.1;

            //plot intensity traces
            Plot iplot = new Plot($intWindowTitle, "time [s]", "Intensity", intTime, intTrace1);
            iplot.setFrameSize(intWindowDimX, intWindowDimY);
            iplot.setLimits(intTime[1], intTime[nopit - 1], iminsc, imaxsc);

            if ((ipx1 - ipx2) != 0 || (ipy1 - ipy2) != 0) {
                iplot.setColor(java.awt.Color.GREEN);
            } else {
                iplot.setColor(java.awt.Color.BLUE);
            }

            iplot.setJustification(Plot.CENTER);

            if ((ipx1 - ipx2) != 0 || (ipy1 - ipy2) != 0) {
                iplot.addLabel(0.5, 0, $intWindowTitle + " (" + ipx1 + ", " + ipy1 + ") and (" + ipx2 + ", " + ipy2 + ")");
            } else {
                iplot.addLabel(0.5, 0, $intWindowTitle + " (" + ipx1 + ", " + ipy1 + ")");
            }

            iplot.draw();

            if ((ipx1 - ipx2) != 0 || (ipy1 - ipy2) != 0) {
                iplot.setColor(java.awt.Color.RED);
                iplot.addPoints(intTime, intTrace2, Plot.LINE);
                iplot.setColor(java.awt.Color.GREEN);
            }

            if (intWindow == null || intWindow.isClosed() == true) {	// create new plot if window doesn't exist, or reuse existing window
                intWindow = iplot.show();
                intWindow.setLocation(intWindowPosX, intWindowPosY);
            } else {
                intWindow.drawPlot(iplot);
                intWindow.setTitle($intWindowTitle);
            }
        }

        if (cormode == 2) {
            for (int x = 0; x < nopit; x++) {	// find minimum and maximum values in intensity trace 1 that will be plotted
                if (imaxsc < intTrace1[x]) {
                    imaxsc = intTrace1[x];
                }
                if (iminsc > intTrace1[x]) {
                    iminsc = intTrace1[x];
                }
            }

            iminsc -= iminsc * 0.1;		// maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            imaxsc += imaxsc * 0.1;

            //plot intensity traces
            Plot iplot = new Plot($intWindowTitle, "time [s]", "Intensity", intTime, intTrace1);
            iplot.setFrameSize(intWindowDimX, intWindowDimY);
            iplot.setLimits(intTime[1], intTime[nopit - 1], iminsc, imaxsc);
            iplot.addLabel(0, 0, "Average intensity trace of whole ROI");
            iplot.setJustification(Plot.LEFT);
            iplot.setColor(java.awt.Color.BLUE);
            iplot.draw();

            if (intWindow == null || intWindow.isClosed() == true) {	// create new plot if window doesn't exist, or reuse existing window
                intWindow = iplot.show();
                intWindow.setLocation(intWindowPosX, intWindowPosY);
            } else {
                intWindow.drawPlot(iplot);
                intWindow.setTitle($intWindowTitle);
            }
        }

    }

    // plot the residuals
    public void plotResiduals(int rpx1, int rpy1) {
        // rpxq, rpy1: pixel coordinates

        if (!plotResCurves) {
            return;
        }

        int ct = 0;
        if (cbFitModel.getSelectedItem() == "DC-FCCS") {
            ct = 2; 	// index ct ensures that when DC-FCCS is selected residual of CCF fit is plotted
        }

        if (isgpupresent == 1) {
            for (int i = 1; i < chanum; i++) {
                res[ct][rpx1][rpy1][i] = acf[ct][rpx1][rpy1][i] - fitacf[ct][rpx1][rpy1][i];
            }
        }
        double rminsc = res[ct][rpx1][rpy1][1];		// minimum and maximum setting for plot window
        double rmaxsc = res[ct][rpx1][rpy1][1];

        if (!tbFCCSDisplay.isSelected()) {
            for (int x = 1; x < chanum; x++) {	// find minimum and maximum values in residuals that will be plotted
                if (rmaxsc < res[ct][rpx1][rpy1][x]) {
                    rmaxsc = res[ct][rpx1][rpy1][x];
                }
                if (rminsc > res[ct][rpx1][rpy1][x]) {
                    rminsc = res[ct][rpx1][rpy1][x];
                }
            }
        }

        if (cbFitModel.getSelectedItem() == "DC-FCCS" && tbFCCSDisplay.isSelected()) {
            for (int z = 1; z <= (chanum - 1); z++) {	// find minimum and maximum values in residuals that will be plotted
                if (rmaxsc < res[0][rpx1][rpy1][z]) {
                    rmaxsc = res[0][rpx1][rpy1][z];
                }
                if (rmaxsc < res[1][rpx1][rpy1][z]) {
                    rmaxsc = res[1][rpx1][rpy1][z];
                }
                if (rmaxsc < res[2][rpx1][rpy1][z]) {
                    rmaxsc = res[2][rpx1][rpy1][z];
                }
                if (rminsc > res[0][rpx1][rpy1][z]) {
                    rminsc = res[0][rpx1][rpy1][z];
                }
                if (rminsc > res[1][rpx1][rpy1][z]) {
                    rminsc = res[1][rpx1][rpy1][z];
                }
                if (rminsc > res[2][rpx1][rpy1][z]) {
                    rminsc = res[2][rpx1][rpy1][z];
                }
            }
        }

        rminsc -= rminsc * 0.1;		// maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
        rmaxsc += rmaxsc * 0.1;

        //plot residuals
        Plot rplot = new Plot($resWindowTitle, "time [s]", "Res", lagtime, res[ct][rpx1][rpy1]);
        rplot.setFrameSize(resWindowDimX, resWindowDimY);
        rplot.setLogScaleX();
        rplot.setLimits(lagtime[1], lagtime[chanum - 1], rminsc, rmaxsc);
        rplot.setColor(java.awt.Color.BLUE);
        rplot.setJustification(Plot.CENTER);
        rplot.addLabel(0.5, 0, " Residuals (" + rpx1 * pixbinX + ", " + rpy1 * pixbinY + ")");
        rplot.draw();

        if (cbFitModel.getSelectedItem() == "DC-FCCS" && tbFCCSDisplay.isSelected()) {
            rplot.setColor(java.awt.Color.GREEN);
            rplot.addPoints(lagtime, res[0][rpx1][rpy1], Plot.LINE);
            rplot.setColor(java.awt.Color.RED);
            rplot.addPoints(lagtime, res[1][rpx1][rpy1], Plot.LINE);
            rplot.setColor(java.awt.Color.BLUE);
            rplot.addPoints(lagtime, res[2][rpx1][rpy1], Plot.LINE);
        }

        if (resWindow == null || resWindow.isClosed() == true) {	// create new plot if window doesn't exist, or reuse existing window
            resWindow = rplot.show();
            resWindow.setLocation(resWindowPosX, resWindowPosY);
        } else {
            resWindow.drawPlot(rplot);
            resWindow.setTitle($resWindowTitle);
        }
    }

    // plot the standard deviation
    public void plotSD(int sdpx1, int sdpy1, int cormode) {
        // sdpx1, sdpy1: pixel coordinates
        // cormode 0: ACF, cormode 2: CCF
        double sdminsc = sdacf[cormode][sdpx1][sdpy1][1];		// minimum and maximum setting for plot window
        double sdmaxsc = sdacf[cormode][sdpx1][sdpy1][1];

        if (cbFitModel.getSelectedItem() == "DC-FCCS" && tbFCCSDisplay.isSelected()) {
            for (int z = 1; z <= (chanum - 1); z++) {	// find minimum and maximum values in residuals that will be plotted
                if (sdmaxsc < sdacf[0][sdpx1][sdpy1][z]) {
                    sdmaxsc = sdacf[0][sdpx1][sdpy1][z];
                }
                if (sdmaxsc < sdacf[1][sdpx1][sdpy1][z]) {
                    sdmaxsc = sdacf[1][sdpx1][sdpy1][z];
                }
                if (sdmaxsc < sdacf[2][sdpx1][sdpy1][z]) {
                    sdmaxsc = sdacf[2][sdpx1][sdpy1][z];
                }
                if (sdminsc > sdacf[0][sdpx1][sdpy1][z]) {
                    sdminsc = sdacf[0][sdpx1][sdpy1][z];
                }
                if (sdminsc > sdacf[1][sdpx1][sdpy1][z]) {
                    sdminsc = sdacf[1][sdpx1][sdpy1][z];
                }
                if (sdminsc > sdacf[2][sdpx1][sdpy1][z]) {
                    sdminsc = sdacf[2][sdpx1][sdpy1][z];
                }
            }
        } else {
            for (int x = 1; x < chanum; x++) {	// find minimum and maximum values in sdacf that will be plotted
                if (sdmaxsc < sdacf[cormode][sdpx1][sdpy1][x]) {
                    sdmaxsc = sdacf[cormode][sdpx1][sdpy1][x];
                }
                if (sdminsc > sdacf[cormode][sdpx1][sdpy1][x]) {
                    sdminsc = sdacf[cormode][sdpx1][sdpy1][x];
                }
            }
        }

        sdminsc -= sdminsc * 0.1;		// maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
        sdmaxsc += sdmaxsc * 0.1;

        //plot standard deviation
        Plot sdplot = new Plot($sdWindowTitle, "time [s]", "SD", lagtime, sdacf[cormode][sdpx1][sdpy1]);
        sdplot.setFrameSize(sdWindowDimX, sdWindowDimY);
        sdplot.setLogScaleX();
        sdplot.setLimits(lagtime[1], lagtime[chanum - 1], sdminsc, sdmaxsc);
        sdplot.setColor(java.awt.Color.BLUE);
        sdplot.setJustification(Plot.CENTER);
        sdplot.addLabel(0.5, 0, " StdDev (" + sdpx1 * pixbinX + ", " + sdpy1 * pixbinY + ")");
        sdplot.draw();

        if (cbFitModel.getSelectedItem() == "DC-FCCS" && tbFCCSDisplay.isSelected()) {
            sdplot.setColor(java.awt.Color.GREEN);
            sdplot.addPoints(lagtime, sdacf[0][sdpx1][sdpy1], Plot.LINE);
            sdplot.setColor(java.awt.Color.RED);
            sdplot.addPoints(lagtime, sdacf[1][sdpx1][sdpy1], Plot.LINE);
        }

        if (sdWindow == null || sdWindow.isClosed() == true) {	// create new plot if window doesn't exist, or reuse existing window
            sdWindow = sdplot.show();
            sdWindow.setLocation(sdWindowPosX, sdWindowPosY);
        } else {
            sdWindow.drawPlot(sdplot);
            sdWindow.setTitle($sdWindowTitle);
        }
    }

    // MSD plots; MSD calcualtion requires prior ACF calculation
    public void plotMSD(Roi MSDroi, int cormode, boolean map) {
        // MSDroi: ROI for which MSD is to be calculated
        // cormode 1: ACF, cormode 2: ACF or CCF, cormode 3: DC-FCCS
        // map: is roi selected in parameter window (reproduction of values only) or is it from a new calculation?

        int cutoff = chanum - 1;
        int ct = 0;

        if (cbFitModel.getSelectedItem() == "DC-FCCS") {
            ct = 2; //in cormode 1, when DC-FCCS is selected, MSD of the cross-correlation is ploted
        }
        Rectangle MSDrect = MSDroi.getBounds();

        int msdpx1 = (int) Math.ceil(MSDrect.getX() / pixbinX);
        int msdpy1 = (int) Math.ceil(MSDrect.getY() / pixbinY);
        int msdpxf = (int) Math.floor((MSDrect.getX() + MSDrect.getWidth() - binningX) / pixbinX);
        int msdpyf = (int) Math.floor((MSDrect.getY() + MSDrect.getHeight() - binningY) / pixbinY);
        if (map == true) {										// if the ROI is selected in the parameter map
            msdpx1 = (int) MSDrect.getX();
            msdpy1 = (int) MSDrect.getY();
            if (cfXDistance < 0) {
                msdpx1 = msdpx1 - (int) Math.floor(cfXDistance / pixbinX);
            }
            if (cfYDistance < 0) {
                msdpy1 = msdpy1 - (int) Math.floor(cfYDistance / pixbinY);
            }
            msdpxf = (int) (msdpx1 + MSDrect.getWidth());
            msdpyf = (int) (msdpy1 + MSDrect.getHeight());
        }
        double msdminsc = msd[ct][msdpx1][msdpy1][1];		// minimum and maximum setting for plot window
        double msdmaxsc = msd[ct][msdpx1][msdpy1][1];

        if (cormode == 1) {

            int i = chanum - 1;
            while ((msd[ct][msdpx1][msdpy1][i] == 0) && (i > 1)) {
                i--;
            }
            cutoff = i + 1;

            double[] msdtime = new double[cutoff];
            msdtime = Arrays.copyOfRange(lagtime, 0, cutoff);
            double[] msdvalue = new double[cutoff];
            msdvalue = Arrays.copyOfRange(msd[ct][msdpx1][msdpy1], 0, cutoff);

            for (int x = 1; x < cutoff; x++) {
                if (msdmaxsc < msd[ct][msdpx1][msdpy1][x]) {
                    msdmaxsc = msd[ct][msdpx1][msdpy1][x];
                }
                if (msdminsc > msd[ct][msdpx1][msdpy1][x]) {
                    msdminsc = msd[ct][msdpx1][msdpy1][x];
                }
            }

            msdminsc -= msdminsc * 0.1;		// maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            msdmaxsc += msdmaxsc * 0.1;

            //plot MSD
            Plot msdplot = new Plot($msdWindowTitle, "time [s]", "MSD (um^2)", msdtime, msdvalue);
            msdplot.setFrameSize(msdWindowDimX, msdWindowDimY);
            msdplot.setLimits(msdtime[0], msdtime[cutoff - 1], msdminsc, msdmaxsc);
            msdplot.setColor(java.awt.Color.BLUE);
            msdplot.setJustification(Plot.CENTER);
            msdplot.addLabel(0.5, 0, " MSD (" + msdpx1 * pixbinX + ", " + msdpy1 * pixbinY + ")");
            msdplot.draw();

            if (msdWindow == null || msdWindow.isClosed() == true) {	// create new plot if window doesn't exist, or reuse existing window
                msdWindow = msdplot.show();
                msdWindow.setLocation(msdWindowPosX, msdWindowPosY);
            } else {
                msdWindow.drawPlot(msdplot);
                msdWindow.setTitle($msdWindowTitle);
            }
        }

        if (cormode == 2) {

            int i = chanum;						// find cutoff channel for the whole plot 
            boolean zerofound = true;
            while (zerofound && i > 2) {
                i--;
                zerofound = false;
                for (int x = msdpx1; x <= msdpxf; x++) {
                    for (int y = msdpy1; y <= msdpyf; y++) {
                        if (((!map && MSDroi.contains(x * pixbinX, y * pixbinY) && MSDroi.contains(x * pixbinX, y * pixbinY + binningX - 1) && MSDroi.contains(x * pixbinX + binningX - 1, y * pixbinY) && MSDroi.contains(x * pixbinX + binningX - 1, y * pixbinY + binningY - 1)) || (map && MSDroi.contains(x, y))) && (msd[ct][x][y][i] == 0.0)) {
                            zerofound = true;
                        }
                    }
                }
            }
            cutoff = i + 1;

            double[] msdtime = new double[cutoff];
            msdtime = Arrays.copyOfRange(lagtime, 0, cutoff);
            double[] msdvalue = new double[cutoff];

            for (int x = msdpx1; x <= msdpxf; x++) {		// find minimum and maximum values in MSD plots that will be plotted
                for (int y = msdpy1; y <= msdpyf; y++) {
                    for (int z = 1; z < cutoff; z++) {
                        if ((!map && MSDroi.contains(x * pixbinX, y * pixbinY) && MSDroi.contains(x * pixbinX, y * pixbinY + binningY - 1) && MSDroi.contains(x * pixbinX + binningX - 1, y * pixbinY) && MSDroi.contains(x * pixbinX + binningX - 1, y * pixbinY + binningY - 1)) || (map == true && MSDroi.contains(x, y))) {
                            if (msdmaxsc < msd[ct][x][y][z]) {
                                msdmaxsc = msd[ct][x][y][z];
                            }
                            if (msdminsc > msd[ct][x][y][z]) {
                                msdminsc = msd[ct][x][y][z];
                            }
                        }
                    }
                }
            }

            msdminsc -= msdminsc * 0.1;		// maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            msdmaxsc += msdmaxsc * 0.1;

            //plot MSD
            Plot msdplot = new Plot($msdWindowTitle, "time [s]", "MSD (um^2)", empty, empty);
            msdplot.setFrameSize(msdWindowDimX, msdWindowDimY);
            msdplot.setLimits(msdtime[0], msdtime[cutoff - 1], msdminsc, msdmaxsc);
            msdplot.setColor(java.awt.Color.BLUE);
            msdplot.setJustification(Plot.CENTER);
            msdplot.addLabel(0.5, 0, "MSDs of pixels in the ROIs at " + binningX + "x" + binningY + "binning");
            msdplot.draw();

            if (msdWindow == null || msdWindow.isClosed() == true) {	// create new plot if window doesn't exist, or reuse existing window
                msdWindow = msdplot.show();
                msdWindow.setLocation(msdWindowPosX, msdWindowPosY);
            } else {
                msdWindow.drawPlot(msdplot);
                msdWindow.setTitle($msdWindowTitle);
            }

            for (int y = msdpy1; y <= msdpyf; y++) {
                for (int x = msdpx1; x <= msdpxf; x++) {
                    if ((!map && MSDroi.contains(x * pixbinX, y * pixbinY) && MSDroi.contains(x * pixbinX, y * pixbinY + binningY - 1) && MSDroi.contains(x * pixbinX + binningX - 1, y * pixbinY) && MSDroi.contains(x * pixbinX + binningX - 1, y * pixbinY + binningY - 1)) || (map == true && MSDroi.contains(x, y))) {
                        msdvalue = Arrays.copyOfRange(msd[ct][x][y], 0, cutoff);
                        msdplot.setColor(java.awt.Color.BLUE);
                        msdplot.addPoints(msdtime, msdvalue, Plot.LINE);
                    }
                }
            }
            msdWindow.drawPlot(msdplot);
        }

        if (cormode == 3) {

            int i = chanum - 1;
            while ((msd[0][msdpx1][msdpy1][i] == 0 || msd[1][msdpx1][msdpy1][i] == 0 || msd[2][msdpx1][msdpy1][i] == 0) && (i > 2)) {
                i--;
            }
            cutoff = i + 1;

            double[] msdtime = new double[cutoff];
            msdtime = Arrays.copyOfRange(lagtime, 0, cutoff);
            double[][] msdvalue = new double[3][cutoff];
            msdvalue[0] = Arrays.copyOfRange(msd[0][msdpx1][msdpy1], 0, cutoff);
            msdvalue[1] = Arrays.copyOfRange(msd[1][msdpx1][msdpy1], 0, cutoff);
            msdvalue[2] = Arrays.copyOfRange(msd[2][msdpx1][msdpy1], 0, cutoff);

            for (int z = 1; z < cutoff; z++) { // find minimum and maximum values in msd that will be plotted
                if (msdmaxsc < msd[0][msdpx1][msdpy1][z]) {
                    msdmaxsc = msd[0][msdpx1][msdpy1][z];
                }
                if (msdmaxsc < msd[1][msdpx1][msdpy1][z]) {
                    msdmaxsc = msd[1][msdpx1][msdpy1][z];
                }
                if (msdmaxsc < msd[2][msdpx1][msdpy1][z]) {
                    msdmaxsc = msd[2][msdpx1][msdpy1][z];
                }
                if (msdminsc > msd[0][msdpx1][msdpy1][z]) {
                    msdminsc = msd[0][msdpx1][msdpy1][z];
                }
                if (msdminsc > msd[1][msdpx1][msdpy1][z]) {
                    msdminsc = msd[1][msdpx1][msdpy1][z];
                }
                if (msdminsc > msd[2][msdpx1][msdpy1][z]) {
                    msdminsc = msd[2][msdpx1][msdpy1][z];
                }
            }

            msdminsc -= msdminsc * 0.1;		// maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            msdmaxsc += msdmaxsc * 0.1;

            //plot MSD
            Plot msdplot = new Plot($msdWindowTitle, "time [s]", "MSD (um^2)", empty, empty);
            msdplot.setFrameSize(msdWindowDimX, msdWindowDimY);
            msdplot.setLimits(msdtime[0], msdtime[cutoff - 1], msdminsc, msdmaxsc);
            msdplot.setColor(java.awt.Color.BLUE);
            msdplot.setJustification(Plot.CENTER);
            msdplot.addLabel(0.5, 0, " MSD (" + msdpx1 * pixbinX + ", " + msdpy1 * pixbinY + ")");
            msdplot.draw();
            msdplot.setColor(java.awt.Color.GREEN);
            msdplot.addPoints(msdtime, msdvalue[0], Plot.LINE);
            msdplot.setColor(java.awt.Color.RED);
            msdplot.addPoints(msdtime, msdvalue[1], Plot.LINE);
            msdplot.setColor(java.awt.Color.BLUE);
            msdplot.addPoints(msdtime, msdvalue[2], Plot.LINE);

            if (msdWindow == null || msdWindow.isClosed() == true) {	// create new plot if window doesn't exist, or reuse existing window
                msdWindow = msdplot.show();
                msdWindow.setLocation(msdWindowPosX, msdWindowPosY);
            } else {
                msdWindow.drawPlot(msdplot);
                msdWindow.setTitle($msdWindowTitle);
            }
        }

        if (cormode == 4) {

            int i = chanum - 1;
            while ((msdaveacf[i] == 0) && (i > 1)) {
                i--;
            }
            cutoff = i + 1;

            double[] msdtime = new double[cutoff];
            msdtime = Arrays.copyOfRange(lagtime, 0, cutoff);
            double[] msdvalue = new double[cutoff];
            msdvalue = Arrays.copyOfRange(msdaveacf, 0, cutoff);

            for (int x = 1; x < cutoff; x++) {
                if (msdmaxsc < msdaveacf[x]) {
                    msdmaxsc = msdaveacf[x];
                }
                if (msdminsc > msdaveacf[x]) {
                    msdminsc = msdaveacf[x];
                }
            }

            msdminsc -= msdminsc * 0.1;		// maximum scales are to be 10% larger than maximum value and 10% smaller than minimum value
            msdmaxsc += msdmaxsc * 0.1;

            //plot MSD
            Plot msdplot = new Plot($msdWindowTitle, "time [s]", "MSD (um^2)", msdtime, msdvalue);
            msdplot.setFrameSize(msdWindowDimX, msdWindowDimY);
            msdplot.setLimits(msdtime[0], msdtime[cutoff - 1], msdminsc, msdmaxsc);
            msdplot.setColor(java.awt.Color.BLUE);
            msdplot.setJustification(Plot.CENTER);
            msdplot.addLabel(0.5, 0, " MSD of average ACF");
            msdplot.draw();

            if (msdWindow == null || msdWindow.isClosed() == true) {	// create new plot if window doesn't exist, or reuse existing window
                msdWindow = msdplot.show();
                msdWindow.setLocation(msdWindowPosX, msdWindowPosY);
            } else {
                msdWindow.drawPlot(msdplot);
                msdWindow.setTitle($msdWindowTitle);
            }
        }
    }

    // create parameter maps
    public void createParaImp(int wx, int hy) {
        // wx, hy: image width and height
        int nochannels = 1;		// number of correlation channels; 1 for FCS, 3 for DC-FCCS (for cross-correlation and autocorrelation in the 2 channels respectively)
        int xm = mincposx;		//introducing pixel shifts in order to align all channels
        int ym = mincposy;
        int cshift = 0;			// index shift for renumbering the correlation channels so that cross-correlation parameter maps are in the upper slices of the stack
        int cm;					// index of the correlation channel in the loop
        int cq = 0;				// additional one frame for q map in DC-FCCS mode

        impPara1exists = true;

        if (cbFitModel.getSelectedItem() == "DC-FCCS") {
            nochannels = 3;
            cshift = 2;
            cq = 1;
        }

        if (impPara1 != null) {		// close parameter window if it exists
            impPara1.close();
        }

        if (histWin != null && histWin.isClosed() == false) {		// close histogram window if it exists
            histWin.close();
        }

        if (doFiltering) {
            userThreshold[0] = true;
            userThreshold[2] = true;
            ReadFilteringFrame(); 									// reads the current Thresholding settings
        } else {
            userThreshold[2] = false;
        }

        impPara1 = IJ.createImage($impPara1Title, "GRAY32", wx, hy, nochannels * (noparam + 3) + cq);	// create a stack for the fit parameters plus chi2, blocking status, filtering mask plus q map

        for (int m = 0; m < nochannels; m++) {							//loop over individual correlation channels and the cross-correlation
            cm = (m + cshift) % 3;
            for (int x = 0; x < wx; x++) {
                for (int y = 0; y < hy; y++) {
                    if (doFiltering) {											//if thresholding on, apply the thresholds	
                        if (filterPix(cm, x + xm, y + ym)) {
                            pixvalid[cm][x + xm][y + ym] = 1.0;
                        } else {
                            pixvalid[cm][x + xm][y + ym] = Double.NaN;
                        }
                    } else {
                        if (pixfitted[cm][x + xm][y + ym]) {
                            pixvalid[cm][x + xm][y + ym] = 1.0;
                        } else {
                            pixvalid[cm][x + xm][y + ym] = Double.NaN;
                        }
                    }
                }
            }
        }

        if (cq == 1) {
            ImageProcessor ipPara1 = impPara1.getStack().getProcessor(nochannels * (noparam + 3) + cq);	// if q map calculated set the stack to the last frame
            for (int x = 0; x < wx; x++) {					// fill the frame with q map
                for (int y = 0; y < hy; y++) {
                    ipPara1.putPixelValue(x, y, CCFq[x + xm][y + ym] * pixvalid[0][x + xm][y + ym] * pixvalid[1][x + xm][y + ym] * pixvalid[2][x + xm][y + ym]);
                    if (userThreshold[4]) {
                        if ((pixvalid[2][x + xm][y + ym] != 1.0) && (pixvalid[0][x + xm][y + ym] * pixvalid[1][x + xm][y + ym] == 1)) {
                            ipPara1.putPixelValue(x, y, 0.0);
                        }
                    }
                }
            }
        }

        for (int m = 0; m < nochannels; m++) {							//loop over individual correlation channels and the cross-correlation
            cm = (m + cshift) % 3;

            for (int p = 0; p < noparam; p++) {					// put pixel values for fit parameters in the maps
                ImageProcessor ipPara1 = impPara1.getStack().getProcessor((m * noparam) + p + 1);
                for (int x = 0; x < wx; x++) {
                    for (int y = 0; y < hy; y++) {
                        ipPara1.putPixelValue(x, y, fitres[cm][x + xm][y + ym][p] * pixvalid[cm][x + xm][y + ym]);
                    }
                }
            }

            ImageProcessor ipPara1 = impPara1.getStack().getProcessor(nochannels * noparam + m + 1);	// set the stack to the frame for chi2
            for (int x = 0; x < wx; x++) {					// fill the frame with chi2 values
                for (int y = 0; y < hy; y++) {
                    ipPara1.putPixelValue(x, y, chi2[cm][x + xm][y + ym] * pixvalid[cm][x + xm][y + ym]);
                }
            }

            ipPara1 = impPara1.getStack().getProcessor(nochannels * (noparam + 1) + m + 1);	// set the stack to the frame for blocked values
            for (int x = 0; x < wx; x++) {					// fill the frame with blocked values
                for (int y = 0; y < hy; y++) {
                    ipPara1.putPixelValue(x, y, blocked[cm][x + xm][y + ym]);
                }
            }

            ipPara1 = impPara1.getStack().getProcessor(nochannels * (noparam + 2) + m + 1);	// set the stack to the frame for filtering mask
            for (int x = 0; x < wx; x++) {					// fill the frame with filtering mask
                for (int y = 0; y < hy; y++) {
                    ipPara1.putPixelValue(x, y, pixvalid[cm][x + xm][y + ym]);
                }
            }

            impPara1.show();
            impPara1Win = impPara1.getWindow();
            impPara1Win.setLocation(para1PosX, para1PosY);

            for (int i = noparam; i >= 1; i--) {		// set label for each parameter map
                impPara1.setSlice(i + m * noparam);
                IJ.run("Set Label...", "label=" + $param[i - 1] + $channel[m]);
            }

            impPara1.setSlice(nochannels * noparam + m + 1);		// set label for Chi2
            IJ.run("Set Label...", "label=" + $param[noparam] + $channel[m]);

            impPara1.setSlice(nochannels * (noparam + 1) + m + 1);		// set label for blocking success
            IJ.run("Set Label...", "label=" + $param[noparam + 1] + $channel[m]);

            impPara1.setSlice(nochannels * (noparam + 3) + cq);		// set label for q map
            IJ.run("Set Label...", "label=" + $param[noparam + 3]);

            impPara1.setSlice(nochannels * (noparam + 2) + m + 1);		// set label for filtering mask
            IJ.run("Set Label...", "label=" + $param[noparam + 2] + $channel[m]);
        }

        IJ.run(impPara1, "Red Hot", "");	// apply "Fire" LUT
        IJ.run(impPara1, "Original Scale", ""); 	//first set image to original scale
        IJ.run(impPara1, "Set... ", "zoom=" + scimp + " x=" + (int) Math.floor(wx / 2) + " y=" + (int) Math.floor(hy / 2)); //then zoom to fit within application
        IJ.run("In [+]", ""); 	// This needs to be used since ImageJ 1.48v to set the window to the right size; 
        // this might be a bug and is an ad hoc solution for the moment; before only the "Set" command was necessary

        impPara1.setSlice(1);				// set back to slice 1 for viewing
        IJ.run(impPara1, "Enhance Contrast", "saturated=0.35");	//autoscaling the contrast for slice 1 
        Component[] impPara1comp = impPara1Win.getComponents();	// check which component is the scrollbar and add an AdjustmentListener
        ScrollbarWithLabel impPara1scrollbar;
        for (int i = 0; i < impPara1comp.length; i++) {
            if (impPara1comp[i] instanceof ScrollbarWithLabel) {
                impPara1scrollbar = (ScrollbarWithLabel) impPara1Win.getComponent(i);
                impPara1scrollbar.addAdjustmentListener(impPara1Adjusted);
            }
        }

        // create histogram window
        impPara1.setSlice(1);
        double histMin = impPara1.getStatistics().min;
        double histMax = impPara1.getStatistics().max;
        int histYMax = impPara1.getStatistics().histYMax;
        int pixelCount = impPara1.getStatistics().pixelCount;
        double stdDev = impPara1.getStatistics().stdDev;
        int q1 = 0;			// determine first quartile
        int countQ = 0;
        while (countQ < Math.ceil(pixelCount / 4.0)) {
            countQ += impPara1.getStatistics().getHistogram()[q1++];
        }
        int q3 = 0;			// determine third quartile
        countQ = 0;
        while (countQ < Math.ceil(3.0 * pixelCount / 4.0)) {
            countQ += impPara1.getStatistics().getHistogram()[q3++];
        }
        double iqr = (q3 - q1) * impPara1.getStatistics().binSize;		// calculate interquartile distance
        int nBins;
        if (iqr > 0) {
            nBins = (int) Math.ceil(Math.cbrt(pixelCount) * (histMax - histMin) / (2.0 * iqr)); // Freedman-Diaconis rule for number of bins in histogram
        } else {
            nBins = 10;
        }

        $histWindowTitle = $param[0] + " - " + $impTitle;
        histWin = new HistogramWindow($histWindowTitle, impPara1, nBins, histMin, histMax, histYMax);
        histWin.setLocationAndSize(histPosX, histPosY, histDimX, histDimY);

        impPara1.setSlice(1);
        impPara1Can = impPara1.getCanvas();		// get canvas
        impPara1Can.setFocusable(true);			// set focusable

        // add listeners
        impPara1Can.addMouseListener(para1MouseClicked);
    }

    // Check whether scrollbar was changed
    AdjustmentListener impPara1Adjusted = (AdjustmentEvent e) -> {
        int slice = impPara1.getSlice();
        int nochannels = 1;
        int chan;
        int par;
        if (cbFitModel.getSelectedItem() == "DC-FCCS") {
            nochannels = 3;
        }
        if (slice <= noparam * nochannels) {
            chan = (int) Math.floor((slice - 1) / noparam);
            par = slice - 1 - chan * noparam;
        } else if (slice <= (noparam + 1) * nochannels) {
            chan = slice - noparam * nochannels - 1;	// Chi2
            par = noparam;
        } else if (slice <= (noparam + 2) * nochannels) {
            chan = slice - (noparam + 1) * nochannels - 1;	// blocking success
            par = noparam + 1;
        } else if (slice <= (noparam + 3) * nochannels) {
            chan = slice - (noparam + 2) * nochannels - 1;	// filtering mask
            par = noparam + 2;
        } else {
            chan = 0;	// q map
            par = noparam + 3;
        }
        impPara1.setSlice(slice);
        IJ.run(impPara1, "Enhance Contrast", "saturated=0.35");
        double histMin = impPara1.getStatistics().min;
        double histMax = impPara1.getStatistics().max;
        int histYMax = impPara1.getStatistics().histYMax;
        int pixelCount = impPara1.getStatistics().pixelCount;
        double stdDev = impPara1.getStatistics().stdDev;
        int q1 = 0;			// determine first quartile
        int countQ = 0;
        while (countQ < Math.ceil(pixelCount / 4.0)) {
            countQ += impPara1.getStatistics().getHistogram()[q1++];
        }
        int q3 = 0; // determine third quartile
        countQ = 0;
        while (countQ < Math.ceil(3.0 * pixelCount / 4.0)) {
            countQ += impPara1.getStatistics().getHistogram()[q3++];
        }
        double iqr = (q3 - q1) * impPara1.getStatistics().binSize; // calculate interquartile distance
        int nBins;
        if (iqr > 0) {
            nBins = (int) Math.ceil(Math.cbrt(pixelCount) * (histMax - histMin) / (2.0 * iqr)); // Freedman-Diaconis rule for number of bins in histogram
        } else {
            nBins = 10;
        }
        if (plotParaHist) {
            $histWindowTitle = $param[par] + $channel[chan] + " - " + $impTitle;
            if (histWin.isClosed() == true) {
                histWin = new HistogramWindow($histWindowTitle, impPara1, nBins, histMin, histMax, histYMax);
                histWin.setLocationAndSize(histPosX, histPosY, histDimX, histDimY);
            } else {
                histWin.showHistogram(impPara1, nBins, histMin, histMax);
                histWin.setTitle($histWindowTitle);
            }
        }
    };

    // retrieve fit values from parameter maps window and display in Fit panel if user clicks a pixel in the parameter map window
    MouseListener para1MouseClicked = new MouseListener() {
        @Override
        public void mouseClicked(MouseEvent e) {

            int px = e.getX(); // get event coordinates
            int py = e.getY();
            int cpx = (int) Math.floor(impPara1Can.offScreenX(px)); //convert to pixel numbers
            int cpy = (int) Math.floor(impPara1Can.offScreenY(py));
            int cpxpar = cpx;
            int cpypar = cpy;
            if (cfXDistance < 0) {
                cpxpar = cpxpar - (int) Math.floor(cfXDistance / pixbinX);
            }
            if (cfYDistance < 0) {
                cpypar = cpypar - (int) Math.floor(cfYDistance / pixbinY);
            }

            if (e.getClickCount() == 1) {
                // update the fit window with the parameters form the selected pixel
                tfParamN.setText(IJ.d2s(fitres[0][cpxpar][cpypar][0], decformat));
                tfParamD.setText(IJ.d2s(fitres[0][cpxpar][cpypar][1] * Math.pow(10, 12), decformat));
                tfParamVx.setText(IJ.d2s(fitres[0][cpxpar][cpypar][2] * Math.pow(10, 6), decformat));
                tfParamVy.setText(IJ.d2s(fitres[0][cpxpar][cpypar][3] * Math.pow(10, 6), decformat));
                tfParamG.setText(IJ.d2s(fitres[0][cpxpar][cpypar][4], decformat2));
                tfParamF2.setText(IJ.d2s(fitres[0][cpxpar][cpypar][5], decformat));
                tfParamD2.setText(IJ.d2s(fitres[0][cpxpar][cpypar][6] * Math.pow(10, 12), decformat));
                tfParamF3.setText(IJ.d2s(fitres[0][cpxpar][cpypar][7], decformat));
                tfParamD3.setText(IJ.d2s(fitres[0][cpxpar][cpypar][8] * Math.pow(10, 12), decformat));
                tfParamFtrip.setText(IJ.d2s(fitres[0][cpxpar][cpypar][9], decformat));
                tfParamTtrip.setText(IJ.d2s(fitres[0][cpxpar][cpypar][10] * Math.pow(10, 6), decformat));
            } else if (e.getClickCount() == 2) {
                tfParamN.setText("NaN");
                tfParamD.setText("NaN");
                tfParamVx.setText("NaN");
                tfParamVy.setText("NaN");
                tfParamG.setText("NaN");
                tfParamF2.setText("NaN");
                tfParamD2.setText("NaN");
                tfParamF3.setText("NaN");
                tfParamD3.setText("NaN");
                tfParamFtrip.setText("NaN");
                tfParamTtrip.setText("NaN");
                CCFq[cpxpar][cpypar] = Double.NaN;
                for (int x = 0; x < 3; x++) {
                    for (int q = 0; q < noparam; q++) {
                        fitres[x][cpxpar][cpypar][q] = Double.NaN;
                    }
                    chi2[x][cpxpar][cpypar] = Double.NaN;
                    blocked[x][cpxpar][cpypar] = Double.NaN;
                    pixfitted[x][cpxpar][cpypar] = false;
                    for (int y = 0; y < chanum; y++) {
                        acf[x][cpxpar][cpypar][y] = 0;
                        sdacf[x][cpxpar][cpypar][y] = 0;
                        varacf[x][cpxpar][cpypar][y] = 0;
                        msd[x][cpxpar][cpypar][y] = 0;
                        fitacf[x][cpxpar][cpypar][y] = 0;
                        res[x][cpxpar][cpypar][y] = 0;
                    }
                }
                updateParaImp();
            }

            if (impPara1.getOverlay() != null) {	//set a 1-pixel ROI
                impPara1.getOverlay().clear();
                impPara1.setOverlay(imp.getOverlay());
            }
            Roi pararoi1 = new Roi(cpx, cpy, 1, 1);
            pararoi1.setStrokeColor(java.awt.Color.BLUE);
            impPara1.setRoi(pararoi1);

            // plot the correlation function of the pixel
            if ((cbFitModel.getSelectedItem() == "FCS") || (cbFitModel.getSelectedItem() == "DC-FCCS" && !tbFCCSDisplay.isSelected())) {
                plotCF(pararoi1, 1, true);
                calcIntensityTrace(imp, cpxpar, cpypar, cpxpar, cpypar, firstframe, lastframe);
                plotIntensityTrace(pararoi1, 1);
                if (doMSD) {
                    plotMSD(pararoi1, 1, true);
                }
            } else {
                if (cbFitModel.getSelectedItem() == "DC-FCCS" && tbFCCSDisplay.isSelected()) {
                    plotCF(pararoi1, 3, true);
                    calcIntensityTrace(imp, cpxpar, cpypar, cpxpar, cpypar, firstframe, lastframe);
                    plotIntensityTrace(pararoi1, 3);
                    if (doMSD) {
                        plotMSD(pararoi1, 3, true);
                    }
                }
            }
            impPara1.deleteRoi();	// delete the ROI again so that image statistics uses full image
        }

        @Override
        public void mousePressed(MouseEvent e) {
        }	// other mouse events have no action associated yet

        @Override
        public void mouseReleased(MouseEvent e) {
            Roi pararoi = impPara1.getRoi();
            if (pararoi != null && pararoi.getFeretsDiameter() > 1) {
                pararoi.setStrokeColor(java.awt.Color.BLUE);
                plotCF(pararoi, 2, true);
                if (doMSD) {
                    plotMSD(pararoi, 2, true);
                }
            }
        }

        @Override
        public void mouseEntered(MouseEvent e) {
        }

        @Override
        public void mouseExited(MouseEvent e) {
        }
    };

    // update single parameter value in parameter window
    public void updateParaImp() {
        int wx = maxcposx - mincposx + 1;		// current image dimensions; this takes acount of binning
        int hy = maxcposy - mincposy + 1;
        int slice = impPara1.getSlice();		// current slice position
        int nochannels = 1;						// number of correlation channels; 1 for FCS, 3 for DC-FCCS (for cross-correlation and autocorrelation in the 2 channels respectively)
        int xm = mincposx;						// introducing pixel shifts in order to align all channels
        int ym = mincposy;
        int cshift = 0;							// index shift for renumbering the correlation channels so that cross-correlation parameter maps are in the upper slices of the stack
        int cm;									// index of the correlation channel in the loop
        int chan;
        int par;
        int cq = 0;								// equals to 1 is q map is calculated to ensure a slice is added for the map

        if (cbFitModel.getSelectedItem() == "DC-FCCS") {
            nochannels = 3;
            cshift = 2;
            cq = 1;
        }

        if (slice <= noparam * nochannels) {
            chan = (int) Math.floor((slice - 1) / noparam);
            par = slice - 1 - chan * noparam;
        } else if (slice <= (noparam + 1) * nochannels) {
            chan = slice - noparam * nochannels - 1;	// Chi2
            par = noparam;
        } else if (slice <= (noparam + 2) * nochannels) {
            chan = slice - (noparam + 1) * nochannels - 1;	// blocking success
            par = noparam + 1;
        } else if (slice <= (noparam + 3) * nochannels) {
            chan = slice - (noparam + 2) * nochannels - 1;	// filtering mask
            par = noparam + 2;
        } else {
            chan = 0;	// q map
            par = noparam + 3;
        }

        if (!impPara1.isVisible()) {
            createParaImp(maxcposx - mincposx + 1, maxcposy - mincposy + 1);
            WindowManager.setCurrentWindow(impPara1Win);
        }

        if (doFiltering) {
            userThreshold[0] = true;
            userThreshold[2] = true;
            ReadFilteringFrame();
        } else {
            userThreshold[2] = false;
        }

        for (int m = 0; m < nochannels; m++) {							//loop over individual correlation channels and the cross-correlation
            cm = (m + cshift) % 3;

            for (int x = 0; x < wx; x++) {
                for (int y = 0; y < hy; y++) {
                    if (doFiltering) {											//if thresholding on, apply the thresholds	
                        if (filterPix(cm, x + xm, y + ym)) {
                            pixvalid[cm][x + xm][y + ym] = 1.0;
                        } else {
                            pixvalid[cm][x + xm][y + ym] = Double.NaN;
                        }
                    } else {
                        if (pixfitted[cm][x + xm][y + ym]) {
                            pixvalid[cm][x + xm][y + ym] = 1.0;
                        } else {
                            pixvalid[cm][x + xm][y + ym] = Double.NaN;
                        }
                    }
                }
            }

            for (int p = 0; p < noparam; p++) {
                ImageProcessor ipPara1 = impPara1.getStack().getProcessor((m * noparam) + p + 1);
                for (int x = 0; x < wx; x++) {
                    for (int y = 0; y < hy; y++) {
                        ipPara1.putPixelValue(x, y, fitres[cm][x + xm][y + ym][p] * pixvalid[cm][x + xm][y + ym]);
                    }
                }
            }

            ImageProcessor ipPara1 = impPara1.getStack().getProcessor(nochannels * noparam + m + 1);	// set the stack to the frame for chi2
            for (int x = 0; x < wx; x++) {					// fill the frame with chi2 values
                for (int y = 0; y < hy; y++) {
                    ipPara1.putPixelValue(x, y, chi2[cm][x + xm][y + ym] * pixvalid[cm][x + xm][y + ym]);
                }
            }

            ipPara1 = impPara1.getStack().getProcessor(nochannels * (noparam + 1) + m + 1);	// set the stack to the frame for blocked values
            for (int x = 0; x < wx; x++) {					// fill the frame with blocked values
                for (int y = 0; y < hy; y++) {
                    ipPara1.putPixelValue(x, y, blocked[cm][x + xm][y + ym]);
                }
            }

            ipPara1 = impPara1.getStack().getProcessor(nochannels * (noparam + 2) + m + 1);	// set the stack to the frame for filtering mask
            for (int x = 0; x < wx; x++) {					// fill the frame with filtering mask
                for (int y = 0; y < hy; y++) {
                    ipPara1.putPixelValue(x, y, pixvalid[cm][x + xm][y + ym]);
                }
            }
        }

        if (cq == 1) {
            ImageProcessor ipPara1 = impPara1.getStack().getProcessor(nochannels * (noparam + 3) + cq);	// if q map calculated set the stack to the last frame
            for (int x = 0; x < wx; x++) {					// fill the frame with q map
                for (int y = 0; y < hy; y++) {
                    ipPara1.putPixelValue(x, y, CCFq[x + xm][y + ym] * pixvalid[0][x + xm][y + ym] * pixvalid[1][x + xm][y + ym] * pixvalid[2][x + xm][y + ym]);
                    if (userThreshold[4]) {
                        if ((pixvalid[2][x + xm][y + ym] != 1.0) && (pixvalid[0][x + xm][y + ym] * pixvalid[1][x + xm][y + ym] == 1)) {
                            ipPara1.putPixelValue(x, y, 0.0);
                        }
                    }
                }
            }
        }

        if (((par < noparam) && paramfit[par]) || par == noparam || par == noparam + 3) {	// create histograms for parameters which were not held fixed including chi2 and q

            impPara1.setSlice(slice);
            double histMin = impPara1.getStatistics().min;
            double histMax = impPara1.getStatistics().max;
            int histYMax = impPara1.getStatistics().histYMax;
            int pixelCount = impPara1.getStatistics().pixelCount;
            double stdDev = impPara1.getStatistics().stdDev;
            int q1 = 0;			// determine first quartile
            int countQ = 0;
            while (countQ < Math.ceil(pixelCount / 4.0)) {
                countQ += impPara1.getStatistics().getHistogram()[q1++];
            }
            int q3 = 0;			// determine third quartile
            countQ = 0;
            while (countQ < Math.ceil(3.0 * pixelCount / 4.0)) {
                countQ += impPara1.getStatistics().getHistogram()[q3++];
            }
            double iqr = (q3 - q1) * impPara1.getStatistics().binSize;		// calculate interquartile distance
            int nBins;
            if (iqr > 0) {
                nBins = (int) Math.ceil(Math.cbrt(pixelCount) * (histMax - histMin) / (2.0 * iqr)); // Freedman-Diaconis rule for number of bins in histogram
            } else {
                nBins = 10;
            }
            if (plotParaHist) {
                $histWindowTitle = $param[par] + $channel[chan] + " - " + $impTitle;
                if (histWin.isClosed() == true) {
                    histWin = new HistogramWindow($histWindowTitle, impPara1, nBins, histMin, histMax, histYMax);
                    histWin.setLocationAndSize(histPosX, histPosY, histDimX, histDimY);
                } else {
                    histWin.showHistogram(impPara1, nBins, histMin, histMax);
                    histWin.setTitle($histWindowTitle);
                }
            }
        }
        impPara1.setSlice(slice);
        IJ.run(impPara1, "Enhance Contrast", "saturated=0.35");
    }

    public boolean filterPix(int m, int x, int y) {		//boolean function which decides whether a pixel passed the user-defined thresholds (true)
        boolean valid = true;
        int cm = m;
        double[][] thresholds = new double[noparam + 1][2];		//local copy of the thresholds for given channel

        if (!pixfitted[m][x][y]) {
            valid = false;
            return valid;
        }

        if (userThreshold[3]) {
            cm = 0; 					//use the thresholds for 1st channel for all channels if the same thresholds for ACF and CCF are used
        }
        for (int p = 0; p < noparam + 1; p++) {
            System.arraycopy(filterThresholds[cm][p], 0, thresholds[p], 0, 2);
        }

        if (paramfilter[noparam + 1]) {											//find the absolute thresholds for Ginf if filtered relative to N
            thresholds[4][0] = filterThresholds[cm][noparam + 1][0] / fitres[cm][x][y][0];
            thresholds[4][1] = filterThresholds[cm][noparam + 1][1] / fitres[cm][x][y][0];
        }

        for (int p = 0; p < noparam; p++) {					// loop over all parameters and if they have been fitted and selected for filtering, apply thresholds
            if (paramfilter[p]) {
                if (fitres[m][x][y][p] < thresholds[p][0] || fitres[m][x][y][p] > thresholds[p][1]) {
                    valid = false;
                }
            }
        }

        if (paramfilter[noparam]) {
            if (chi2[m][x][y] < thresholds[noparam][0] || chi2[m][x][y] > thresholds[noparam][1]) {
                valid = false;
            }
        }

        return valid;
    }

    /* 
 * 	The following classes contain the fitting functions required for Imaging_FCS
 * 	
 * 	class Line implements ParametricUnivariateFunction: linear function
 * 	class SingleExp implements ParametricUnivariateFunction: single exponential function for bleach correction
 * 	class DoubleExp implements ParametricUnivariateFunction: double exponential function for bleach correction
 * 	class Polynomial implements ParametricUnivariateFunction: polynomial fucntion for the bleach correction
 * 	
     */
    // linear fit for diffusion law plot
    class Line implements ParametricUnivariateFunction {

        @Override
        public double[] gradient(double x, double[] params) {
            double t0 = params[0];
            double b = params[1];

            double[] grad = new double[]{
                1,
                x
            };

            return grad;
        }

        @Override
        public double value(double x, double[] params) {
            double t0 = params[0];
            double b = params[1];
            return t0 + b * x;
        }
    }

    // single exponential for bleach correction
    class SingleExp implements ParametricUnivariateFunction {

        @Override
        public double[] gradient(double x, double[] params) {
            double A = params[0];
            double t1 = params[1];
            double C = params[2];

            double[] grad = new double[]{
                Math.exp(-x / t1),
                A * x * Math.exp(-x / t1) / (Math.pow(t1, 2)),
                1
            };

            return grad;
        }

        @Override
        public double value(double x, double[] params) {
            double A = params[0];
            double t1 = params[1];
            double C = params[2];
            return A * Math.exp(-x / t1) + C;
        }
    }

    // double exponential for bleach correction
    class DoubleExp implements ParametricUnivariateFunction {

        @Override
        public double[] gradient(double x, double[] params) {
            double A = params[0];
            double t1 = params[1];
            double B = params[2];
            double t2 = params[3];
            double C = params[4];

            double[] grad = new double[]{
                Math.exp(-x / t1),
                A * x * Math.exp(-x / t1) / (Math.pow(t1, 2)),
                Math.exp(-x / t2),
                B * x * Math.exp(-x / t2) / (Math.pow(t2, 2)),
                1
            };

            return grad;
        }

        @Override
        public double value(double x, double[] params) {
            double A = params[0];
            double t1 = params[1];
            double B = params[2];
            double t2 = params[3];
            double C = params[4];
            return A * Math.exp(-x / t1) + B * Math.exp(-x / t2) + C;
        }
    }

    // polynomial for bleach correction
    class Polynomial implements ParametricUnivariateFunction {

        @Override
        public double[] gradient(double x, double[] params) {
            int maxord = polyOrder;

            double[] grad = new double[maxord + 1];
            for (int i = 0; i <= maxord; i++) {
                grad[i] = Math.pow(x, i);
            }

            return grad;
        }

        @Override
        public double value(double x, double[] params) {
            int maxord = polyOrder;
            double[] A = new double[maxord + 1];
            for (int i = 0; i <= maxord; i++) {
                A[i] = params[i];
            }

            double val = 0;
            for (int i = 0; i <= maxord; i++) {
                val += A[i] * Math.pow(x, i);
            }

            return val;
        }
    }

    /*
 * Fit functions
 * 
 * class FCS_3p implements ParametricUnivariateFunction: FCS fit including diffusion, flow and spatial cross-correlation for ITIR-FCS and SPIM-FCS; accepts up to 3 components/particles with areas of rectangular shape
 * class FCCS_2p implements ParametricUnivariateFunction: DC-FCCS fit; assumes only diffusion; accepts up to 2 components/particles
 * class GLS_fitFunction implements ParametricUnivariateFunction: generalized least squares fit; transforms the model functions using the regularized covarince matrix
 * public double obsvolFCS_ST2D1p(int dim): Calculate the observation area/volume  
 * public void dataTransform(double[] corav, double[][] covmats): transforms the ACF using the regularized covarince matrix for the GLS fit
 * 
     */
    // FCS model: 3D fit assuming 3 components and flow in x and y direction; this is the general fit formula; parameters can be set to 0 to obtain simpler models
    // the models and their derivation are provided on our website in CDF files (http://staff.science.nus.edu.sg/~chmwt/)
    class FCS_3p implements ParametricUnivariateFunction {
        // general parameters

        double pi = 3.14159265359;
        double sqrpi = Math.sqrt(pi);
        double ax = pixeldimx;
        double ay = pixeldimy;
        double s = psfsize;
        double sz = lsthickness;
        double psfz = 2 * emlambda / Math.pow(10, 9.0) * 1.33 / Math.pow(NA, 2.0); // size of PSF in axial direction
        //double szeff = Math.sqrt( 1 / ( Math.pow(sz, -2.0) + Math.pow(psfz, -2.0) ) ); // convolution of two Gaussians depending on illumination profile and detection PSF
        double szeff = sz;
        double rx = ax * cfXshift / binningX;
        double ry = ay * cfYshift / binningY;

        @Override
        public double[] gradient(double x, double[] params) {
            double[] pareq = new double[noparam];
            int num = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i]) {
                    pareq[i] = params[num];
                    num++;
                } else {
                    pareq[i] = paraminitval[i];
                }
            }

            // note that x is used here as the time variable instead of t; that can be come confusing as x and y are used in the names for the paramaters to indicate spatial directions
            // pareq[0] = N
            // pareq[1] = D
            // pareq[2] = vx
            // pareq[3] = vy
            // pareq[4] = G
            // pareq[5] = F2
            // pareq[6] = D2
            // pareq[7] = F3
            // pareq[8] = D3
            // pareq[9] = Ftrip
            // pareq[10] = Ttrip
            //COMPONENT1
            // help variables, which are dependent on time, to write the full function
            double p0t = Math.sqrt(4 * pareq[1] * x + Math.pow(s, 2));
            double p1xt = ax + rx - pareq[2] * x;
            double p2xt = ax - rx + pareq[2] * x;
            double p3xt = rx - pareq[2] * x;
            double p4xt = 2 * Math.pow(ax, 2) + 3 * Math.pow(rx, 2) - 6 * x * rx * pareq[2] + 3 * Math.pow(x * pareq[2], 2);
            double p5xt = Math.pow(p3xt, 2) + Math.pow(p1xt, 2);
            double p6xt = Math.pow(p3xt, 2) + Math.pow(p2xt, 2);
            double p7xt = 2 * (Math.pow(ax, 2) + Math.pow(rx, 2) - 2 * x * rx * pareq[2] + Math.pow(x * pareq[2], 2));
            double p1yt = ay + ry - pareq[3] * x;
            double p2yt = ay - ry + pareq[3] * x;
            double p3yt = ry - pareq[3] * x;
            double p4yt = 2 * Math.pow(ay, 2) + 3 * Math.pow(ry, 2) - 6 * x * ry * pareq[3] + 3 * Math.pow(x * pareq[3], 2);
            double p5yt = Math.pow(p3yt, 2) + Math.pow(p1yt, 2);
            double p6yt = Math.pow(p3yt, 2) + Math.pow(p2yt, 2);
            double p7yt = 2 * (Math.pow(ay, 2) + Math.pow(ry, 2) - 2 * x * ry * pareq[3] + Math.pow(x * pareq[3], 2));
            double pexpxt = Math.exp(-Math.pow(p1xt / p0t, 2)) + Math.exp(-Math.pow(p2xt / p0t, 2)) - 2 * Math.exp(-Math.pow(p3xt / p0t, 2));
            double perfxt = p1xt * Erf.erf(p1xt / p0t) + p2xt * Erf.erf(p2xt / p0t) - 2 * p3xt * Erf.erf(p3xt / p0t);
            double dDpexpxt = 2 * Math.exp(-p4xt / Math.pow(p0t, 2)) * (Math.exp(p5xt / Math.pow(p0t, 2)) + Math.exp(p6xt / Math.pow(p0t, 2)) - 2 * Math.exp(p7xt / Math.pow(p0t, 2)));
            double dvxperfxt = (Erf.erf(p2xt / p0t) + 2 * Erf.erf(p3xt / p0t) - Erf.erf(p1xt / p0t)) * x;
            double pexpyt = Math.exp(-Math.pow(p1yt / p0t, 2)) + Math.exp(-Math.pow(p2yt / p0t, 2)) - 2 * Math.exp(-Math.pow(p3yt / p0t, 2));
            double dDpexpyt = 2 * Math.exp(-p4yt / Math.pow(p0t, 2)) * (Math.exp(p5yt / Math.pow(p0t, 2)) + Math.exp(p6yt / Math.pow(p0t, 2)) - 2 * Math.exp(p7yt / Math.pow(p0t, 2)));
            double dvyperfyt = (Erf.erf(p2yt / p0t) + 2 * Erf.erf(p3yt / p0t) - Erf.erf(p1yt / p0t)) * x;
            double perfyt = p1yt * Erf.erf(p1yt / p0t) + p2yt * Erf.erf(p2yt / p0t) - 2 * p3yt * Erf.erf(p3yt / p0t);

            //CF for the lateral dimension (x, y) and its derivative for D
            double plat = (p0t / sqrpi * pexpxt + perfxt) * (p0t / sqrpi * pexpyt + perfyt) / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double dDplat = (1 / (sqrpi * p0t)) * (dDpexpyt * x * (p0t / sqrpi * pexpxt + perfxt) + dDpexpxt * x * (p0t / sqrpi * pexpyt + perfyt)) / (4 * Math.pow(ax * ay, 2) / fitobsvol);

            //CF for the axial dimension (z) and its derivative for D
            double pspim = 1 / Math.sqrt(1 + (4 * pareq[1] * x) / Math.pow(szeff, 2));
            double dDpspim = -4 * x / (2 * Math.pow(szeff, 2) * Math.pow(Math.sqrt(1 + (4 * pareq[1] * x) / Math.pow(szeff, 2)), 3));

            double acf1 = plat * pspim;

            //COMPONENT2
            // help variables, which are dependent on time, to write the full function
            double p0t2 = Math.sqrt(4 * pareq[6] * x + Math.pow(s, 2));
            double p1xt2 = ax + rx - pareq[2] * x;
            double p2xt2 = ax - rx + pareq[2] * x;
            double p3xt2 = rx - pareq[2] * x;
            double p4xt2 = 2 * Math.pow(ax, 2) + 3 * Math.pow(rx, 2) - 6 * x * rx * pareq[2] + 3 * Math.pow(x * pareq[2], 2);
            double p5xt2 = Math.pow(p3xt2, 2) + Math.pow(p1xt2, 2);
            double p6xt2 = Math.pow(p3xt2, 2) + Math.pow(p2xt2, 2);
            double p7xt2 = 2 * (Math.pow(ax, 2) + Math.pow(rx, 2) - 2 * x * rx * pareq[2] + Math.pow(x * pareq[2], 2));
            double p1yt2 = ay + ry - pareq[3] * x;
            double p2yt2 = ay - ry + pareq[3] * x;
            double p3yt2 = ry - pareq[3] * x;
            double p4yt2 = 2 * Math.pow(ay, 2) + 3 * Math.pow(ry, 2) - 6 * x * ry * pareq[3] + 3 * Math.pow(x * pareq[3], 2);
            double p5yt2 = Math.pow(p3yt2, 2) + Math.pow(p1yt2, 2);
            double p6yt2 = Math.pow(p3yt2, 2) + Math.pow(p2yt2, 2);
            double p7yt2 = 2 * (Math.pow(ay, 2) + Math.pow(ry, 2) - 2 * x * ry * pareq[3] + Math.pow(x * pareq[3], 2));
            double pexpxt2 = Math.exp(-Math.pow(p1xt2 / p0t2, 2)) + Math.exp(-Math.pow(p2xt2 / p0t2, 2)) - 2 * Math.exp(-Math.pow(p3xt2 / p0t2, 2));
            double perfxt2 = p1xt2 * Erf.erf(p1xt2 / p0t2) + p2xt2 * Erf.erf(p2xt2 / p0t2) - 2 * p3xt2 * Erf.erf(p3xt2 / p0t2);
            double dDpexpxt2 = 2 * Math.exp(-p4xt2 / Math.pow(p0t2, 2)) * (Math.exp(p5xt2 / Math.pow(p0t2, 2)) + Math.exp(p6xt2 / Math.pow(p0t2, 2)) - 2 * Math.exp(p7xt2 / Math.pow(p0t2, 2)));
            double dvxperfxt2 = (Erf.erf(p2xt2 / p0t2) + 2 * Erf.erf(p3xt2 / p0t2) - Erf.erf(p1xt2 / p0t2)) * x;
            double pexpyt2 = Math.exp(-Math.pow(p1yt2 / p0t2, 2)) + Math.exp(-Math.pow(p2yt2 / p0t2, 2)) - 2 * Math.exp(-Math.pow(p3yt2 / p0t2, 2));
            double dDpexpyt2 = 2 * Math.exp(-p4yt2 / Math.pow(p0t2, 2)) * (Math.exp(p5yt2 / Math.pow(p0t2, 2)) + Math.exp(p6yt2 / Math.pow(p0t2, 2)) - 2 * Math.exp(p7yt2 / Math.pow(p0t2, 2)));
            double dvyperfyt2 = (Erf.erf(p2yt2 / p0t2) + 2 * Erf.erf(p3yt2 / p0t2) - Erf.erf(p1yt2 / p0t2)) * x;
            double perfyt2 = p1yt2 * Erf.erf(p1yt2 / p0t2) + p2yt2 * Erf.erf(p2yt2 / p0t2) - 2 * p3yt2 * Erf.erf(p3yt2 / p0t2);

            //CF for the lateral dimension (x, y) and its derivative for D
            double plat2 = (p0t2 / sqrpi * pexpxt2 + perfxt2) * (p0t2 / sqrpi * pexpyt2 + perfyt2) / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double dDplat2 = (1 / (sqrpi * p0t2)) * (dDpexpyt2 * x * (p0t2 / sqrpi * pexpxt2 + perfxt2) + dDpexpxt2 * x * (p0t2 / sqrpi * pexpyt2 + perfyt2)) / (4 * Math.pow(ax * ay, 2) / fitobsvol);

            //CF for the axial dimension (z) and its derivative for D
            double pspim2 = 1 / Math.sqrt(1 + (4 * pareq[6] * x) / Math.pow(szeff, 2));
            double dDpspim2 = -4 * x / (2 * Math.pow(szeff, 2) * Math.pow(Math.sqrt(1 + (4 * pareq[6] * x) / Math.pow(szeff, 2)), 3));

            double acf2 = plat2 * pspim2;

            //COMPONENT3
            // help variables, which are dependent on time, to write the full function
            double p0t3 = Math.sqrt(4 * pareq[8] * x + Math.pow(s, 2));
            double p1xt3 = ax + rx - pareq[2] * x;
            double p2xt3 = ax - rx + pareq[2] * x;
            double p3xt3 = rx - pareq[2] * x;
            double p4xt3 = 2 * Math.pow(ax, 2) + 3 * Math.pow(rx, 2) - 6 * x * rx * pareq[2] + 3 * Math.pow(x * pareq[2], 2);
            double p5xt3 = Math.pow(p3xt2, 2) + Math.pow(p1xt2, 2);
            double p6xt3 = Math.pow(p3xt2, 2) + Math.pow(p2xt2, 2);
            double p7xt3 = 2 * (Math.pow(ax, 2) + Math.pow(rx, 2) - 2 * x * rx * pareq[2] + Math.pow(x * pareq[2], 2));
            double p1yt3 = ay + ry - pareq[3] * x;
            double p2yt3 = ay - ry + pareq[3] * x;
            double p3yt3 = ry - pareq[3] * x;
            double p4yt3 = 2 * Math.pow(ay, 2) + 3 * Math.pow(ry, 2) - 6 * x * ry * pareq[3] + 3 * Math.pow(x * pareq[3], 2);
            double p5yt3 = Math.pow(p3yt3, 2) + Math.pow(p1yt3, 2);
            double p6yt3 = Math.pow(p3yt3, 2) + Math.pow(p2yt3, 2);
            double p7yt3 = 2 * (Math.pow(ay, 2) + Math.pow(ry, 2) - 2 * x * ry * pareq[3] + Math.pow(x * pareq[3], 2));
            double pexpxt3 = Math.exp(-Math.pow(p1xt3 / p0t3, 2)) + Math.exp(-Math.pow(p2xt3 / p0t3, 2)) - 2 * Math.exp(-Math.pow(p3xt3 / p0t3, 2));
            double perfxt3 = p1xt3 * Erf.erf(p1xt3 / p0t3) + p2xt3 * Erf.erf(p2xt3 / p0t3) - 2 * p3xt3 * Erf.erf(p3xt3 / p0t3);
            double dDpexpxt3 = 2 * Math.exp(-p4xt3 / Math.pow(p0t3, 2)) * (Math.exp(p5xt3 / Math.pow(p0t3, 2)) + Math.exp(p6xt3 / Math.pow(p0t3, 2)) - 2 * Math.exp(p7xt3 / Math.pow(p0t3, 2)));
            double dvxperfxt3 = (Erf.erf(p2xt3 / p0t3) + 2 * Erf.erf(p3xt3 / p0t3) - Erf.erf(p1xt3 / p0t3)) * x;
            double pexpyt3 = Math.exp(-Math.pow(p1yt3 / p0t3, 2)) + Math.exp(-Math.pow(p2yt3 / p0t3, 2)) - 2 * Math.exp(-Math.pow(p3yt3 / p0t3, 2));
            double dDpexpyt3 = 2 * Math.exp(-p4yt3 / Math.pow(p0t3, 2)) * (Math.exp(p5yt3 / Math.pow(p0t3, 2)) + Math.exp(p6yt3 / Math.pow(p0t3, 2)) - 2 * Math.exp(p7yt3 / Math.pow(p0t3, 2)));
            double dvyperfyt3 = (Erf.erf(p2yt3 / p0t3) + 2 * Erf.erf(p3yt3 / p0t3) - Erf.erf(p1yt3 / p0t3)) * x;
            double perfyt3 = p1yt3 * Erf.erf(p1yt3 / p0t3) + p2yt3 * Erf.erf(p2yt3 / p0t3) - 2 * p3yt3 * Erf.erf(p3yt3 / p0t3);

            // TRIPLET
            double triplet = 1 + pareq[9] / (1 - pareq[9]) * Math.exp(-x / pareq[10]);
            double dtripletFtrip = Math.exp(-x / pareq[10]) * (1 / (1 - pareq[9]) + pareq[9] / Math.pow(1 - pareq[9], 2));
            double dtripletTtrip = Math.exp(-x / pareq[10]) * (pareq[9] * x) / ((1 - pareq[9]) * Math.pow(pareq[10], 2));

            //CF for the lateral dimension (x, y) and its derivative for D
            double plat3 = (p0t3 / sqrpi * pexpxt3 + perfxt3) * (p0t3 / sqrpi * pexpyt3 + perfyt3) / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double dDplat3 = (1 / (sqrpi * p0t3)) * (dDpexpyt3 * x * (p0t3 / sqrpi * pexpxt3 + perfxt3) + dDpexpxt3 * x * (p0t3 / sqrpi * pexpyt3 + perfyt3)) / (4 * Math.pow(ax * ay, 2) / fitobsvol);

            //CF for the axial dimension (z) and its derivative for D
            double pspim3 = 1 / Math.sqrt(1 + (4 * pareq[8] * x) / Math.pow(szeff, 2));
            double dDpspim3 = -4 * x / (2 * Math.pow(szeff, 2) * Math.pow(Math.sqrt(1 + (4 * pareq[8] * x) / Math.pow(szeff, 2)), 3));

            double acf3 = plat3 * pspim3;

            double pf1 = (1 - pareq[5] - pareq[7]) / (1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7]);
            double pf2 = (Math.pow(q2, 2) * pareq[5]) / (1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7]);
            double pf3 = (Math.pow(q3, 2) * pareq[7]) / (1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7]);
            double dfnom = Math.pow(1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7], 3);
            double df21 = 1 - pareq[5] - pareq[7] + q2 * pareq[5] - q3 * pareq[7] + 2 * q2 * pareq[7] - 2 * q2;
            double df22 = Math.pow(q2, 2) * (1 + pareq[5] - pareq[7] - q2 * pareq[5] + q3 * pareq[7]);
            double df23 = 2 * pareq[7] * Math.pow(q3, 2) * (1 - q2);
            double df31 = 1 - pareq[5] - pareq[7] - q2 * pareq[5] + 2 * q3 * pareq[5] - 2 * q3 + q3 * pareq[7];
            double df32 = 2 * pareq[5] * Math.pow(q2, 2) * (1 - q3);
            double df33 = Math.pow(q3, 2) * (1 - pareq[5] + pareq[7] + q2 * pareq[5] - q3 * pareq[7]);

            double pacf = (1 / pareq[0]) * ((1 - pareq[5] - pareq[7]) * acf1 + Math.pow(q2, 2) * pareq[5] * acf2 + Math.pow(q3, 2) * pareq[7] * acf3) / Math.pow(1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7], 2) * triplet + pareq[4];

            double[] grad = new double[]{
                (-1 / Math.pow(pareq[0], 2)) * (pf1 * acf1 + pf2 * acf2 + pf3 * acf3) * triplet,
                (1 / pareq[0]) * pf1 * (plat * dDpspim + pspim * dDplat),
                (1 / pareq[0]) * (pf1 * ((p0t / sqrpi * pexpyt + perfyt) * dvxperfxt) * pspim / (4 * Math.pow(ax * ay, 2) / fitobsvol) + pf2 * ((p0t2 / sqrpi * pexpyt2 + perfyt2) * dvxperfxt2) * pspim2 / (4 * Math.pow(ax * ay, 2) / fitobsvol) + pf3 * ((p0t3 / sqrpi * pexpyt3 + perfyt3) * dvxperfxt3) * pspim3 / (4 * Math.pow(ax * ay, 2) / fitobsvol)) * triplet,
                (1 / pareq[0]) * (pf1 * ((p0t / sqrpi * pexpxt + perfxt) * dvyperfyt) * pspim / (4 * Math.pow(ax * ay, 2) / fitobsvol) + pf2 * ((p0t2 / sqrpi * pexpxt2 + perfxt2) * dvyperfyt2) * pspim2 / (4 * Math.pow(ax * ay, 2) / fitobsvol) + pf3 * ((p0t3 / sqrpi * pexpxt3 + perfxt3) * dvyperfyt3) * pspim3 / (4 * Math.pow(ax * ay, 2) / fitobsvol)) * triplet,
                1,
                (1 / pareq[0]) * (1 / dfnom) * (df21 * acf1 + df22 * acf2 + df23 * acf3) * triplet,
                (1 / pareq[0]) * pf2 * (plat2 * dDpspim2 + pspim2 * dDplat2) * triplet,
                (1 / pareq[0]) * (1 / dfnom) * (df31 * acf1 + df32 * acf2 + df33 * acf3) * triplet,
                (1 / pareq[0]) * pf3 * (plat3 * dDpspim3 + pspim3 * dDplat3) * triplet,
                dtripletFtrip * pacf,
                dtripletTtrip * pacf
            };

            double[] gradret = new double[num]; // return the gradients of the fit model in respect to the fit parameters
            num = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    gradret[num] = grad[i];
                    num++;
                }
            }

            return gradret;
        }

        @Override
        public double value(double x, double[] params) {
            double[] pareq = new double[noparam];
            int num = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i]) {
                    pareq[i] = params[num];
                    num++;
                } else {
                    pareq[i] = paraminitval[i];
                }
            }

            // note that x is used here as the time variable instead of t; that can be come confusing as x and y are used in the names for the paramaters to indicate spatial directions
            // pareq[0] = N
            // pareq[1] = D
            // pareq[2] = vx
            // pareq[3] = vy
            // pareq[4] = G
            // pareq[5] = F2
            // pareq[6] = D2
            // pareq[7] = F3
            // pareq[8] = D3
            // pareq[9] = Ftrip
            // pareq[10] = Dtrip
            //q2 and q3, the brightness of the second and third components are fixed parameters and have been globaly defined; see prepareFit()
            // COMPONENT 1
            // help variables, which are dependent on time, to write the full function
            double p0t = Math.sqrt(4 * pareq[1] * x + Math.pow(s, 2));
            double p1xt = ax + rx - pareq[2] * x;
            double p2xt = ax - rx + pareq[2] * x;
            double p3xt = rx - pareq[2] * x;
            double p1yt = ay + ry - pareq[3] * x;
            double p2yt = ay - ry + pareq[3] * x;
            double p3yt = ry - pareq[3] * x;
            double pexpxt = Math.exp(-Math.pow(p1xt / p0t, 2)) + Math.exp(-Math.pow(p2xt / p0t, 2)) - 2 * Math.exp(-Math.pow(p3xt / p0t, 2));
            double perfxt = p1xt * Erf.erf(p1xt / p0t) + p2xt * Erf.erf(p2xt / p0t) - 2 * p3xt * Erf.erf(p3xt / p0t);
            double pexpyt = Math.exp(-Math.pow(p1yt / p0t, 2)) + Math.exp(-Math.pow(p2yt / p0t, 2)) - 2 * Math.exp(-Math.pow(p3yt / p0t, 2));
            double perfyt = p1yt * Erf.erf(p1yt / p0t) + p2yt * Erf.erf(p2yt / p0t) - 2 * p3yt * Erf.erf(p3yt / p0t);

            double pplane1 = (p0t / sqrpi * pexpxt + perfxt) * (p0t / sqrpi * pexpyt + perfyt) / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double pspim1 = 1 / Math.sqrt(1 + (4 * pareq[1] * x) / Math.pow(szeff, 2));
            double acf1 = pplane1 * pspim1;

            // COMPONENT 2
            // help variables, which are dependent on time, to write the full function
            double p0t2 = Math.sqrt(4 * pareq[6] * x + Math.pow(s, 2));
            double p1xt2 = ax + rx - pareq[2] * x;
            double p2xt2 = ax - rx + pareq[2] * x;
            double p3xt2 = rx - pareq[2] * x;
            double p1yt2 = ay + ry - pareq[3] * x;
            double p2yt2 = ay - ry + pareq[3] * x;
            double p3yt2 = ry - pareq[3] * x;
            double pexpxt2 = Math.exp(-Math.pow(p1xt2 / p0t2, 2)) + Math.exp(-Math.pow(p2xt2 / p0t2, 2)) - 2 * Math.exp(-Math.pow(p3xt2 / p0t2, 2));
            double perfxt2 = p1xt * Erf.erf(p1xt2 / p0t2) + p2xt2 * Erf.erf(p2xt2 / p0t2) - 2 * p3xt2 * Erf.erf(p3xt2 / p0t2);
            double pexpyt2 = Math.exp(-Math.pow(p1yt2 / p0t2, 2)) + Math.exp(-Math.pow(p2yt2 / p0t2, 2)) - 2 * Math.exp(-Math.pow(p3yt2 / p0t2, 2));
            double perfyt2 = p1yt2 * Erf.erf(p1yt2 / p0t2) + p2yt2 * Erf.erf(p2yt2 / p0t2) - 2 * p3yt2 * Erf.erf(p3yt2 / p0t2);

            double pplane2 = (p0t2 / sqrpi * pexpxt2 + perfxt2) * (p0t2 / sqrpi * pexpyt2 + perfyt2) / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double pspim2 = 1 / Math.sqrt(1 + (4 * pareq[6] * x) / Math.pow(szeff, 2));
            double acf2 = pplane2 * pspim2;

            // COMPONENT 3
            // help variables, which are dependent on time, to write the full function
            double p0t3 = Math.sqrt(4 * pareq[8] * x + Math.pow(s, 2));
            double p1xt3 = ax + rx - pareq[2] * x;
            double p2xt3 = ax - rx + pareq[2] * x;
            double p3xt3 = rx - pareq[2] * x;
            double p1yt3 = ay + ry - pareq[3] * x;
            double p2yt3 = ay - ry + pareq[3] * x;
            double p3yt3 = ry - pareq[3] * x;
            double pexpxt3 = Math.exp(-Math.pow(p1xt3 / p0t3, 2)) + Math.exp(-Math.pow(p2xt3 / p0t3, 2)) - 2 * Math.exp(-Math.pow(p3xt3 / p0t3, 2));
            double perfxt3 = p1xt * Erf.erf(p1xt3 / p0t3) + p2xt3 * Erf.erf(p2xt3 / p0t3) - 2 * p3xt3 * Erf.erf(p3xt3 / p0t3);
            double pexpyt3 = Math.exp(-Math.pow(p1yt3 / p0t3, 2)) + Math.exp(-Math.pow(p2yt3 / p0t3, 2)) - 2 * Math.exp(-Math.pow(p3yt3 / p0t3, 2));
            double perfyt3 = p1yt3 * Erf.erf(p1yt3 / p0t3) + p2yt3 * Erf.erf(p2yt3 / p0t3) - 2 * p3yt3 * Erf.erf(p3yt3 / p0t3);

            double pplane3 = (p0t3 / sqrpi * pexpxt3 + perfxt3) * (p0t3 / sqrpi * pexpyt3 + perfyt3) / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double pspim3 = 1 / Math.sqrt(1 + (4 * pareq[8] * x) / Math.pow(szeff, 2));
            double acf3 = pplane3 * pspim3;

            // TRIPLET
            double triplet = 1 + pareq[9] / (1 - pareq[9]) * Math.exp(-x / pareq[10]);

            return (1 / pareq[0]) * ((1 - pareq[5] - pareq[7]) * acf1 + Math.pow(q2, 2) * pareq[5] * acf2 + Math.pow(q3, 2) * pareq[7] * acf3) / Math.pow(1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7], 2) * triplet + pareq[4];
        }
    }

    // DC-FCCS model, applicabale for ITIR-FCCS and SPIM-FCCS
    // the models and their derivation are provided on our website in CDF files (http://staff.science.nus.edu.sg/~chmwt/)
    class FCCS_2p implements ParametricUnivariateFunction {

        // general parameters
        double pi = 3.14159265359;
        double sqrpi = Math.sqrt(pi);
        double a = pixeldimx;
        double s1 = psfsize;
        double s2 = psfsize2;
        double sz1 = lsthickness;
        double sz2 = lsthickness2;
        double psfz1 = 2 * emlambda / Math.pow(10, 9.0) * 1.33 / Math.pow(NA, 2.0); // size of PSF in axial direction
        double psfz2 = 2 * emlambda2 / Math.pow(10, 9.0) * 1.33 / Math.pow(NA, 2.0); // size of PSF in axial direction
        double szeff1 = Math.sqrt(1 / (Math.pow(sz1, -2.0) + Math.pow(psfz1, -2.0))); // convolution of two Gaussians depending on illumination profile and detection PSF
        double szeff2 = Math.sqrt(1 / (Math.pow(sz2, -2.0) + Math.pow(psfz2, -2.0))); // convolution of two Gaussians depending on illumination profile and detection PSF
        double rz = 0; // shift in light sheet position; can be introduced if necessary

        @Override
        public double[] gradient(double x, double[] params) {
            double[] pareq = new double[noparam];
            int num = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    pareq[i] = params[num];
                    num++;
                } else {
                    pareq[i] = paraminitval[i];
                }
            }

            // note that x is used here as the time variable instead of t; that can become confusing as x and y are used in the names for the paramaters to indicate spatial directions
            // pareq[0] = N
            // pareq[1] = D
            // pareq[2] = 0 (no flow in x direction)
            // pareq[3] = 0 (no flow in y direction)
            // pareq[4] = G
            // pareq[5] = F2
            // pareq[6] = D2
            // pareq[7] = 0
            // pareq[8] = 0
            // pareq[9] = 0
            // pareq[10] = 0
            //COMPONENT1
            // help variables, which are dependent on time, to write the full function
            double p1t = Math.sqrt(4 * pareq[1] * x + Math.pow(s1, 2) / 2 + Math.pow(s2, 2) / 2) / a;
            double p10 = Math.sqrt(Math.pow(s1, 2) / 2 + Math.pow(s2, 2) / 2) / a;
            double p2t = Math.sqrt(1 + 8 * pareq[1] * x / (Math.pow(szeff1, 2) + Math.pow(szeff2, 2)));
            double p3t = Math.exp(-1 / Math.pow(p1t, 2)) - 1;
            double p30 = Math.exp(-1 / Math.pow(p10, 2)) - 1;
            double p4t = Erf.erf(1 / p1t);
            double p40 = Erf.erf(1 / p10);
            double acfnum = Math.pow(p1t / sqrpi * p3t + p4t, 2);
            double acfden = Math.pow(p10 / sqrpi * p30 + p40, 2) * p2t;
            double acf1 = acfnum / acfden;

            // COMPONENT 2
            // help variables, which are dependent on time, to write the full function
            double p1t2 = Math.sqrt(4 * pareq[6] * x + Math.pow(s1, 2) / 2 + Math.pow(s2, 2) / 2) / a;
            double p102 = Math.sqrt(Math.pow(s1, 2) / 2 + Math.pow(s2, 2) / 2) / a;
            double p2t2 = Math.sqrt(1 + 8 * pareq[6] * x / (Math.pow(szeff1, 2) + Math.pow(szeff2, 2)));
            double p3t2 = Math.exp(-1 / Math.pow(p1t2, 2)) - 1;
            double p302 = Math.exp(-1 / Math.pow(p102, 2)) - 1;
            double p4t2 = Erf.erf(1 / p1t2);
            double p402 = Erf.erf(1 / p102);
            double acfnum2 = Math.pow(p1t2 / sqrpi * p3t2 + p4t2, 2);
            double acfden2 = Math.pow(p102 / sqrpi * p302 + p402, 2) * p2t2;
            double acf2 = acfnum2 / acfden2;

            double acfc = ((1 - pareq[5]) * acf1 + Math.pow(q2, 2) * pareq[5] * acf2) / Math.pow(1 - pareq[5] + q2 * pareq[5], 2);
            double dDacfa1 = 4 * x / pareq[0] * (p4t + p3t * p1t / sqrpi) / Math.pow(p40 + p30 * p10 / sqrpi, 2) / p2t;
            double dDacfb1 = p3t / Math.pow(a, 2) / sqrpi / p1t - (p4t + p3t * p1t / sqrpi) / (Math.pow(szeff1, 2) + Math.pow(szeff2, 2)) / Math.pow(p2t, 2);
            double dDacf1 = dDacfa1 * dDacfb1;
            double dDacfa2 = 4 * x / pareq[0] * (p4t2 + p3t2 * p1t2 / sqrpi) / Math.pow(p402 + p302 * p102 / sqrpi, 2) / p2t2;
            double dDacfb2 = p3t2 / Math.pow(a, 2) / sqrpi / p1t2 - (p4t2 + p3t2 * p1t2 / sqrpi) / (Math.pow(szeff1, 2) + Math.pow(szeff2, 2)) / Math.pow(p2t2, 2);
            double dDacf2 = dDacfa2 * dDacfb2;

            double dF2acf = (1 / pareq[0]) * ((Math.pow(q2, 2) * acf2 - acf1) / (1 - pareq[5] + q2 * pareq[5]) - 2 * (q2 - 1) * ((1 - pareq[5]) * acf1 + Math.pow(q2, 2) * pareq[5] * acf2) / Math.pow(1 - pareq[5] + q2 * pareq[5], 3));

            double[] grad = new double[]{
                (-1 / Math.pow(pareq[0], 2)) * acfc,
                dDacf1,
                0,
                0,
                1,
                dF2acf,
                dDacf2,
                0,
                0
            };

            double[] gradret = new double[num]; // return the gradients of the fit model in respect to the fit parameters
            num = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    gradret[num] = grad[i];
                    num++;
                }
            }

            return gradret;
        }

        @Override
        public double value(double x, double[] params) {
            double[] pareq = new double[noparam];
            int num = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    pareq[i] = params[num];
                    num++;
                } else {
                    pareq[i] = paraminitval[i];
                }
            }

            // note that x is used here as the time variable instead of t; that can be come confusing as x and y are used in the names for the paramaters to indicate spatial directions
            // pareq[0] = N
            // pareq[1] = D
            // pareq[2] = 0 (no flow in x direction)
            // pareq[3] = 0 (no flow in y direction)
            // pareq[4] = G
            // pareq[5] = F2
            // pareq[6] = D2
            // pareq[7] = 0
            // pareq[8] = 0
            // pareq[9] = 0
            // pareq[10] = 0
            //COMPONENT1
            // help variables, which are dependent on time, to write the full function
            double p1t = Math.sqrt(4 * pareq[1] * x + Math.pow(s1, 2) / 2 + Math.pow(s2, 2) / 2) / a;
            double p10 = Math.sqrt(Math.pow(s1, 2) / 2 + Math.pow(s2, 2) / 2) / a;
            double p2t = Math.sqrt(1 + 8 * pareq[1] * x / (Math.pow(szeff1, 2) + Math.pow(szeff2, 2)));
            double p3t = Math.exp(-1 / Math.pow(p1t, 2)) - 1;
            double p30 = Math.exp(-1 / Math.pow(p10, 2)) - 1;
            double p4t = Erf.erf(1 / p1t);
            double p40 = Erf.erf(1 / p10);
            double acfnum = Math.pow(p1t / sqrpi * p3t + p4t, 2);
            double acfden = Math.pow(p10 / sqrpi * p30 + p40, 2) * p2t;
            double acf1 = acfnum / acfden;

            // COMPONENT 2
            // help variables, which are dependent on time, to write the full function
            double p1t2 = Math.sqrt(4 * pareq[6] * x + Math.pow(s1, 2) / 2 + Math.pow(s2, 2) / 2) / a;
            double p102 = Math.sqrt(Math.pow(s1, 2) / 2 + Math.pow(s2, 2) / 2) / a;
            double p2t2 = Math.sqrt(1 + 8 * pareq[6] * x / (Math.pow(szeff1, 2) + Math.pow(szeff2, 2)));
            double p3t2 = Math.exp(-1 / Math.pow(p1t2, 2)) - 1;
            double p302 = Math.exp(-1 / Math.pow(p102, 2)) - 1;
            double p4t2 = Erf.erf(1 / p1t2);
            double p402 = Erf.erf(1 / p102);
            double acfnum2 = Math.pow(p1t2 / sqrpi * p3t2 + p4t2, 2);
            double acfden2 = Math.pow(p102 / sqrpi * p302 + p402, 2) * p2t2;
            double acf2 = acfnum2 / acfden2;

            return (1 / pareq[0]) * ((1 - pareq[5]) * acf1 + Math.pow(q2, 2) * pareq[5] * acf2) / Math.pow(1 - pareq[5] + q2 * pareq[5], 2) + pareq[4];
        }
    }

    // generalized least square fit; uses FCS_3p, FCS_SX_3p, and FCCS_2p to tranform the data via the covariance matrix
    class GLS_fitFunction implements ParametricUnivariateFunction {

        @Override
        public double[] gradient(double x, double[] params) {
            int num = params.length;							// number of values in the gradient of FCS_3p, which depends on the number of parameters to be fit			
            double[][] gradtau = new double[num][chanum - 1];
            double[] tmpgrad = new double[num]; // temporary variable;
            double[] finalgrad = new double[num];
            int sol = 0;						// get the index of the solution for the particular tau

            for (int g = 1; g < chanum; g++) {	// determine which element is required by checking the lagtime
                if (lagtime[g] == x) {
                    sol = g - 1;
                }
            }

            if (x == lagtime[1]) {
                RealVector[] solution = new RealVector[num];				// vector containing solutions
                ParametricUnivariateFunction function;
                if ("FCS".equals((String) cbFitModel.getSelectedItem())) {	// select the fit model to be used; extra fit models can be added here
                    function = new FCS_3p();
                } else {
                    function = new FCCS_2p();
                }

                for (int y = 1; y < chanum; y++) {							// get all gradients for all lag times
                    tmpgrad = function.gradient(lagtime[y], params);
                    for (int z = 0; z < num; z++) {
                        gradtau[z][y - 1] = tmpgrad[z];
                    }
                }

                for (int z = 0; z < num; z++) {							// solve for a new correlation vector with independent elements
                    DecompositionSolver solver = new LUDecomposition(lowerCholDecCovmats).getSolver();
                    RealVector constants = new ArrayRealVector(gradtau[z]);
                    solution[z] = solver.solve(constants);
                    transTheoreticalGradientACF[z] = solution[z].toArray();			// remember the transformed theoretical function
                }

                for (int z = 0; z < num; z++) {
                    finalgrad[z] = transTheoreticalGradientACF[z][sol];
                }
            } else {
                for (int g = 2; g < chanum; g++) {		// determine which element is required by checking the lagtime
                    if (lagtime[g] == x) {
                        sol = g - 1;
                    }
                }
                for (int z = 0; z < num; z++) {
                    finalgrad[z] = transTheoreticalGradientACF[z][sol];
                }
            }

            return finalgrad;
        }

        @Override
        public double value(double x, double[] params) {
            double[] valtau = new double[chanum - 1];	// array for theoretical ACF before transormation
            double retval;								// return value
            int sol = 0;								// index of the solution for a particular tau

            if (x == lagtime[1]) {
                ParametricUnivariateFunction function;
                if ("FCS".equals((String) cbFitModel.getSelectedItem())) {							// select the fit model to be used; extra fit models can be added here
                    function = new FCS_3p();
                } else {
                    function = new FCCS_2p();
                }

                // calculate the correlation function for this particular set of parameters; do not take the zero lagtime into account
                for (int y = 0; y < chanum - 1; y++) {
                    valtau[y] = function.value(lagtime[y + 1], params);
                }

                // use the regularized covariance matrix to transform the data
                DecompositionSolver solver = new LUDecomposition(lowerCholDecCovmats).getSolver(); // solve for a new correlation vector with independent elements
                RealVector constants = new ArrayRealVector(valtau);
                RealVector solution = solver.solve(constants);

                for (int y = 0; y < chanum - 1; y++) {												// remember the transformed theoretical function
                    transTheoreticalACF[y] = solution.getEntry(y);
                }

                retval = transTheoreticalACF[sol];

            } else {
                for (int g = 2; g < chanum; g++) {
                    if (lagtime[g] == x) {
                        sol = g - 1;
                    }
                }
                retval = transTheoreticalACF[sol];
            }

            return retval;
        }
    }

    // calculation of the observation area; this is used in the Diffusion Law Plot as the y-axis
    // the calculation of the observation area/volume is provided on our website in CDF files (http://www.dbs.nus.edu.sg/lab/BFL/index.html)
    public double obsvolFCS_ST2D1p(int dim) {
        // general parameters
        double pi = 3.14159265359;
        double sqrpi = Math.sqrt(pi);
        double ax = pixeldimx;
        double ay = pixeldimy;
        double s = psfsize;
        double sz = lsthickness;
        double psfz = 2 * emlambda / Math.pow(10, 9.0) * 1.33 / Math.pow(NA, 2.0); // size of PSF in axial direction
        double szeff = Math.sqrt(1 / (Math.pow(sz, -2.0) + Math.pow(psfz, -2.0))); // convolution of two Gaussians depending on illumination profile and detection PSF
        double rx = ax * cfXshift / binningX;
        double ry = ay * cfYshift / binningY;

        // help variables, for t = 0, to write the full fit function
        double p00 = s;
        double p1x0 = ax;
        double p2x0 = ax;
        double p1y0 = ay;
        double p2y0 = ay;
        double pexpx0 = 2 * Math.exp(-Math.pow(p1x0 / p00, 2)) - 2;
        double perfx0 = 2 * p1x0 * Erf.erf(p1x0 / p00);
        double pexpy0 = 2 * Math.exp(-Math.pow(p1y0 / p00, 2)) - 2;
        double perfy0 = 2 * p1y0 * Erf.erf(p1y0 / p00);

        //return (p00/sqrpi * pexpx0 + perfx0) * (p00/sqrpi * pexpy0 + perfy0) * Math.pow(sz, 2);
        if (dim == 2) {
            return 4 * Math.pow(ax * ay, 2) / ((p00 / sqrpi * pexpx0 + perfx0) * (p00 / sqrpi * pexpy0 + perfy0));
        } else {
            //return sqrpi * szeff * 4 * Math.pow(ax*ay, 2)/( (p00/sqrpi * pexpx0 + perfx0) * (p00/sqrpi * pexpy0 + perfy0) );
            return 4 * Math.pow(ax * ay, 2) / ((p00 / sqrpi * pexpx0 + perfx0) * (p00 / sqrpi * pexpy0 + perfy0));
        }

    }

    // Transform data using the regularized covariance matrix for a generalized least squares fit
    public void dataTransform(double[] corav, double[][] covmats) { // this needs to be run before the Bayes fit or any generalized least suqares 
        // corav: correlation function
        // covmats: covariance matrix
        double[] cortmp = new double[chanum - 1];
        for (int x = 0; x < chanum - 1; x++) {					// remove zero lagtime channel from the correlation as it is not fit along
            cortmp[x] = corav[x + 1];
        }

        RealMatrix mat = MatrixUtils.createRealMatrix(covmats);	// lower triangular matrix of the CholeskyDecomposition of covmats
        RealMatrix matL;

        try {													// perfrom the Cholesky decomposition and determine the lower triangular matrix
            CholeskyDecomposition CD = new CholeskyDecomposition(mat);
            matL = CD.getL();
        } catch (NonSquareMatrixException | NonSymmetricMatrixException | NonPositiveDefiniteMatrixException ex) {
            IJ.log(ex.getMessage());
            throw ex;
        }

        DecompositionSolver solver = new LUDecomposition(matL).getSolver(); // solve for a new correlation vector with independent elements
        RealVector constants = new ArrayRealVector(cortmp);
        RealVector solution = solver.solve(constants);
        lowerCholDecCovmats = matL;
        transACF = solution;
    }

    /*
 * FCS simulations
 * 
 *  public void simulateACF(): Decide whether 2D or 3D simulation
 *  public void simulateACF2D(): program for 2D simulations
 *  public void simulateACF3D(): program for 3D simulations
 *  
     */
    public void simulateACF(boolean ask3D) {
        //start simulation; simulateACFInstant keeps a reference so the simualtions can be cancelled
        simulateACFInstant = new simulateACFWorker(ask3D);
        simulateACFInstant.execute();
    }

    public void batchSimulateACF(boolean ask3D) {
        // run a set of simulations
        int Drange;
        int D2range;
        int F2range;

        if (batchDStep == 0 || (batchDEnd - batchDStart) == 0) {
            Drange = 0;
        } else {
            Drange = (int) Math.floor((batchDEnd - batchDStart) / batchDStep);
        }

        if (batchD2Step == 0 || (batchD2End - batchD2Start) == 0) {
            D2range = 0;
        } else {
            D2range = (int) Math.floor((batchD2End - batchD2Start) / batchD2Step);
        }

        if (batchF2Step == 0 || (batchF2End - batchF2Start) == 0) {
            F2range = 0;
        } else {
            F2range = (int) Math.floor((batchF2End - batchF2Start) / batchF2Step);
        }
        for (int i = 0; i <= Drange; i++) {
            for (int j = 0; j <= D2range; j++) {
                for (int k = 0; k <= F2range; k++) {
                    tfSimD1.setText(Double.toString(batchDStart + i * batchDStep));
                    tfSimD2.setText(Double.toString(batchD2Start + j * batchD2Step));
                    tfSimF2.setText(Double.toString(batchF2Start + k * batchF2Step));
                    simulateACFInstant = new simulateACFWorker(ask3D);
                    simulateACFInstant.execute();
                    while (simulateACFInstant.isDone() == false) {
                    }
                }
            }
        }
    }

    public class batchSimulateACFWorker extends SwingWorker<Void, Void> {

        private void failIfInterrupted() throws InterruptedException {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interrupted while simulating ACFs");
            }
        }

        private final boolean ask3D;

        public batchSimulateACFWorker(final boolean ask3D) {
            this.ask3D = ask3D;
        }

        @Override
        protected Void doInBackground() throws Exception {
            if (ask3D) {
                if (Double.parseDouble(tfSigmaZ.getText()) <= 0.0 || Double.parseDouble(tfSigmaZ.getText()) > 100) {
                    IJ.showMessage("LightSheetThickness is either <= 0 or is too large.");
                    return null;
                }
                if (Double.parseDouble(tfNA.getText()) >= 1.33) {
                    IJ.showMessage("For 3D simulations NA has to be smaller than 1.33.");
                    return null;
                }
            }
            batchSim = true;
            btnStopSimulation.setEnabled(true);
            batchSimulateACF(ask3D);
            btnStopSimulation.setEnabled(false);
            batchSim = false;
            return null;
        }
    }

    public class simulateACFWorker extends SwingWorker<Void, Void> {

        private void failIfInterrupted() throws InterruptedException {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interrupted while simulating ACFs");
            }
        }

        private final boolean ask3D;

        public simulateACFWorker(final boolean ask3D) {
            this.ask3D = ask3D;
        }

        @Override
        protected Void doInBackground() throws Exception {
            if (ask3D) {
                if (Double.parseDouble(tfSigmaZ.getText()) <= 0.0 || Double.parseDouble(tfSigmaZ.getText()) > 100) {
                    IJ.showMessage("LightSheetThickness is either <= 0 or is too large.");
                    return null;
                }
                if (Double.parseDouble(tfNA.getText()) >= 1.33) {
                    IJ.showMessage("For 3D simulations NA has to be smaller than 1.33.");
                    return null;
                }
                btnStopSimulation.setEnabled(true);
                simulateACF3D();
            } else {
                btnStopSimulation.setEnabled(true);
                simulateACF2D();
            }
            btnStopSimulation.setEnabled(false);
            return null;
        }
    }

    public void simulateACF2D() {
        int simSeed = 0;
        int simNoParticles = 1000; 						// number of simulated particles
        int simCPS = 10000;								// average count rate per particle per second
        double simTauBleach = 100000;					// bleach time in seconds
        int simPixelnum = 21; 							// width of image in pixels
        double simExtFactor = 1.5; 						// factor by which the simulated area is bigger than the observed area
        int simNoTStep = 50000; 						// number of frames to be simulated
        double simFrameTime = 0.001;					// time resolution of the camera in second
        int simStepsPerFrame = 10;						// simulation steps per frame
        double simD1 = 1.0 / Math.pow(10, 12);			// particle 1 diffusion coefficient
        double simDoutDinRatio = 1.0;					// ratio of diffusion coefficients outside over inside of domains
        double simD2 = 0.1 / Math.pow(10, 12);			// particle 2 diffusion coefficient
        double simD3 = 0.01 / Math.pow(10, 12);			// particle 3 diffusion coefficient
        double simF2 = 0.0;								// fraction of particle 2
        double simF3 = 0.0;								// fraction of particle 3
        double simKon = 1.0;							// on-rate for triplet
        double simKoff = 0.0;							// off-rate for triplet
        int simCameraOffset = 100;						// offset of CCD camera
        double simCameraNoiseFactor = 3.0;				// noise of CCD camera
        double simBleachRadius = 3.0;					// bleach radius
        int simBleachFrame = 10000000;					// frame at which bleach happens
        double simDomainRadius = 100.0;					// Radius of domains
        double simDomainDensity = 0.0;					// Density of domains in number/um2
        double simPin = 1.0;							// Probability to enter domain
        double simPout = 1.0;							// Probability to exit domain
        double simMeshworkSize = 100.0;					// Size of meshes
        double simHopProbability = 1.0;					// hop probability over meshwork barriers
        double simPixelSizeRS = 24 / Math.pow(10, 6);		// pixel size in real space
        double simMag = 63.0;							// objective maginification
        double simWavelength = 614.0 / Math.pow(10, 9);	// observation wavelegnth
        double simNA = 1.0; 							// NA of the objective
        double simSigma0 = 0.8;							// actual resolution
        String[] newSimSettings = new String[nosimsettings];	// an array for reading out the settings in Simulation panel

        newSimSettings[0] = (String) cbSimMode.getSelectedItem();
        newSimSettings[1] = tbSimTrip.getText();
        newSimSettings[2] = tfSimSeed.getText();
        newSimSettings[3] = tfSimParticleNum.getText();
        newSimSettings[4] = tfSimCPS.getText();
        newSimSettings[5] = tfSimTauBleach.getText();
        newSimSettings[6] = tfSimPixelNum.getText();
        newSimSettings[7] = tfSimExtensionFactor.getText();
        newSimSettings[8] = tfSimTimeStepNum.getText();
        newSimSettings[9] = tfSimFrameTime.getText();
        newSimSettings[10] = tfSimStepsPerFrame.getText();
        newSimSettings[11] = tfSimCurrentStepSize.getText();
        newSimSettings[12] = tfSimD1.getText();
        newSimSettings[13] = tfSimDoutDinRatio.getText();
        newSimSettings[14] = tfSimD2.getText();
        newSimSettings[15] = tfSimF2.getText();
        newSimSettings[16] = tfSimD3.getText();
        newSimSettings[17] = tfSimF3.getText();
        newSimSettings[18] = tfSimKon.getText();
        newSimSettings[19] = tfSimKoff.getText();
        newSimSettings[20] = tfSimCameraOffset.getText();
        newSimSettings[21] = tfSimCameraNoiseFactor.getText();
        newSimSettings[22] = tfSimBleachRadius.getText();
        newSimSettings[23] = tfSimBleachFrame.getText();
        newSimSettings[24] = tfDomainRadius.getText();
        newSimSettings[25] = tfDomainDensity.getText();
        newSimSettings[26] = tfPin.getText();
        newSimSettings[27] = tfPout.getText();
        newSimSettings[28] = tfMeshworkSize.getText();
        newSimSettings[29] = tfHopProbability.getText();
        newSimSettings[30] = tfPixelSize.getText();
        newSimSettings[31] = tfMagnification.getText();
        newSimSettings[32] = tfEmLambda.getText();
        newSimSettings[33] = tfNA.getText();
        newSimSettings[34] = tfSigma.getText();

        try {
            simSeed = Integer.parseInt(newSimSettings[2]);
            simNoParticles = Integer.parseInt(newSimSettings[3]); 					// number of simulated particles
            simCPS = Integer.parseInt(newSimSettings[4]); 							// average count rate per particle per second
            simTauBleach = Double.parseDouble(newSimSettings[5]);					// bleach time in seconds
            simPixelnum = Integer.parseInt(newSimSettings[6]); 						// width of image in pixels
            simExtFactor = Double.parseDouble(newSimSettings[7]); 					// factor by which the simulated area is bigger than the observed area
            simNoTStep = Integer.parseInt(newSimSettings[8]); 						// number of frames to be simulated
            simFrameTime = Double.parseDouble(newSimSettings[9]);					// time resolution of the camera in second
            simStepsPerFrame = Integer.parseInt(newSimSettings[10]);				// steps of simulations doen for each frame
            simD1 = Double.parseDouble(newSimSettings[12]) / Math.pow(10, 12);		// particle 1 diffusion coefficient
            simDoutDinRatio = Double.parseDouble(newSimSettings[13]);				// ratio of diffusion coefficients outside over inside of domains
            simD2 = Double.parseDouble(newSimSettings[14]) / Math.pow(10, 12);		// particle 2 diffusion coefficient
            simD3 = Double.parseDouble(newSimSettings[16]) / Math.pow(10, 12);		// particle 3 diffusion coefficient
            simF2 = Double.parseDouble(newSimSettings[15]);							// fraction of particle 2
            simF3 = Double.parseDouble(newSimSettings[17]);							// fraction of particle 3
            simKon = Double.parseDouble(newSimSettings[18]);						// on-rate for triplet
            simKoff = Double.parseDouble(newSimSettings[19]);						// off-rate for triplet
            simCameraOffset = Integer.parseInt(newSimSettings[20]);					// offset of CCD camera
            simCameraNoiseFactor = Integer.parseInt(newSimSettings[21]);			// noise of CCD camera
            simBleachRadius = Double.parseDouble(newSimSettings[22]) / Math.pow(10, 6);		// bleach radius
            simBleachFrame = Integer.parseInt(newSimSettings[23]);							// frame at which bleach happens
            simDomainRadius = Double.parseDouble(newSimSettings[24]) / Math.pow(10, 9);		// Radius of domains
            simDomainDensity = Double.parseDouble(newSimSettings[25]);						// Density of domains in number/um2
            simPin = Double.parseDouble(newSimSettings[26]);								// Probability to enter domain
            simPout = Double.parseDouble(newSimSettings[27]);								// Probability to exit domain
            simMeshworkSize = Double.parseDouble(newSimSettings[28]) / Math.pow(10, 9);		// Size of meshes
            simHopProbability = Double.parseDouble(newSimSettings[29]);						// hop probability over meshwork barriers
            simPixelSizeRS = Double.parseDouble(newSimSettings[30]) / Math.pow(10, 6);		// pixel size in real space
            simMag = Double.parseDouble(newSimSettings[31]);								// objective maginification
            simWavelength = Double.parseDouble(newSimSettings[32]) / Math.pow(10, 9);			// observation wavelegnth
            simNA = Double.parseDouble(newSimSettings[33]); 								// NA of the objective
            simSigma0 = Double.parseDouble(newSimSettings[34]);								// actual resolution
        } catch (NumberFormatException nfe) {
            IJ.showMessage("One of the values in the simulation window does not have the right format (integer or double).");
            throw new NumberFormatException("Number format error.");
        }

        double simTStep = simFrameTime / simStepsPerFrame;
        double simDarkF = simKoff / (simKoff + simKon);		//fraction of molecules in the dark state
        if (simDoutDinRatio <= 0) {
            IJ.showMessage("Dout/Din <= 0 is not allowed");
            return;
        }

        int numOfSeeds = 18;
        int[] simSeedArray = new int[numOfSeeds];							// array of simulation Seeds so that the random number generators are different
        double simPixelSize = simPixelSizeRS / simMag;						// pixel size in object space
        double simPSFSize = 0.5 * simSigma0 * simWavelength / simNA;			// PSF size
        double simGridSize = simPixelnum * simPixelSize;					// gridsize; i.e. the size of the pixel area of the detector
        double simMidPos = simGridSize / 2.0;  									// middel position
        double simSizeLL = -simExtFactor * simGridSize; 					// lower limit of the simulation area
        double simSizeUL = simExtFactor * simGridSize; 						// upper limit of the simulation area
        double simDetectorSize = simGridSize / 2.0; 							// the detector extends from -simDetectorSize to simDetectorSize, i.e. it is the same as simGridSize
        double simPhotonsPerStep = Math.floor(simCPS * simTStep + 0.5); 	// number of photons per particle and time step
        double bleachFactor = 2.0;											// the 2.0 ensures that no bleaching happens if simTauBleach is 0
        if (simTauBleach != 0) {
            bleachFactor = Math.exp(-simTStep / simTauBleach);
        }
        double[] blinkFactor = new double[2];
        blinkFactor[0] = Math.exp(-simTStep * simKon);
        blinkFactor[1] = Math.exp(-simTStep * simKoff);

        double gridLength = (simSizeUL - simSizeLL);									// length of full simualtion grid
        double gridMidPos = gridLength / 2;												// half length of full simulation grid
        int numberofdomains = (int) Math.ceil(Math.pow(gridLength * Math.pow(10, 6), 2) * simDomainDensity);
        double[][] domains = new double[numberofdomains][3];
        int subgridnum = (int) Math.ceil(gridLength / (simDomainRadius * 10)) + 1;		// number N of elements in a NxN array 
        double subgridsize = gridLength / (subgridnum - 1);					 			// gridsize in the NxN array
        int maxdomainperarea = (int) Math.ceil(Math.pow(subgridsize / (simDomainRadius * 0.5), 2.0));
        int maxct = 0;
        int[][][] domainsorted = new int[1][1][1];										// subgridsize defines a grid with sizes larger than the largest domain radius
        int[][][] dsortmp = new int[subgridnum][subgridnum][maxdomainperarea];
        int[][] dctr = new int[subgridnum][subgridnum];									// temporary counter
        int num1 = (int) Math.round(simNoParticles * (1 - simF2 - simF3));				// divide particle into their types according to their fractions (1- F2 - F3), F2, F3
        int num2 = (int) Math.round(simNoParticles * simF2);
        int num3 = (int) Math.round(simNoParticles * simF3);
        double[][] particles = new double[simNoParticles][5]; 							// array for particle positions (0:x, 1:y), whether particle is bleached (2) or in dark state (4) and if particle is in domain then (3) contains domain number
        ImagePlus impSim = IJ.createImage("2D Simulation", "GRAY16", simPixelnum, simPixelnum, simNoTStep);

        if (simSeed == 0) {
            Arrays.fill(simSeedArray, 0);
        } else {
            for (int x = 0; x < numOfSeeds; x++) {
                simSeedArray[x] = simSeed + (int) Math.pow(x, 2.0);
            }
        }

        int cs = 0;
        UniformGenerator rugxpos = new UniformGenerator(simSizeLL, simSizeUL, simSeedArray[cs++]);
        UniformGenerator rugypos = new UniformGenerator(simSizeLL, simSizeUL, simSeedArray[cs++]);
        UniformGenerator ruig = new UniformGenerator(simSeedArray[cs++]);
        UniformGenerator rugbf = new UniformGenerator(simSeedArray[cs++]);
        UniformGenerator rugpin = new UniformGenerator(simSeedArray[cs++]);
        UniformGenerator rugpout = new UniformGenerator(simSeedArray[cs++]);
        UniformGenerator rugphop = new UniformGenerator(simSeedArray[cs++]);
        UniformGenerator rugblink = new UniformGenerator(simSeedArray[cs++]);
        GaussianGenerator rgg1 = new GaussianGenerator(0, Math.sqrt(2 * simD1 * simTStep), simSeedArray[cs++]);
        GaussianGenerator rgg2 = new GaussianGenerator(0, Math.sqrt(2 * simD2 * simTStep), simSeedArray[cs++]);
        GaussianGenerator rgg3 = new GaussianGenerator(0, Math.sqrt(2 * simD3 * simTStep), simSeedArray[cs++]);
        GaussianGenerator rggpsf = new GaussianGenerator(0, simPSFSize, simSeedArray[cs++]);
        PoissonGenerator rpgphoton = new PoissonGenerator(simTStep * simCPS, simSeedArray[cs++]);
        //PoissonGenerator rpgnoise = new PoissonGenerator(simCameraNoiseFactor, simSeedArray[cs++]);
        GaussianGenerator rggnoise = new GaussianGenerator(0, Math.sqrt(simCameraNoiseFactor), simSeedArray[cs++]);
        GaussianGenerator rggdom1 = new GaussianGenerator(0, Math.sqrt(2 * simD1 / simDoutDinRatio * simTStep), simSeedArray[cs++]);
        GaussianGenerator rggdom2 = new GaussianGenerator(0, Math.sqrt(2 * simD2 / simDoutDinRatio * simTStep), simSeedArray[cs++]);
        GaussianGenerator rggdom3 = new GaussianGenerator(0, Math.sqrt(2 * simD3 / simDoutDinRatio * simTStep), simSeedArray[cs++]);
        GaussianGenerator rggdrad = new GaussianGenerator(simDomainRadius, simDomainRadius / 10, simSeedArray[cs++]);

        if (simDomainFlag) {	// create domains
            int counter = 1;    // the 0 position is not used
            int totcount = 0;
            IJ.showStatus("creating non-overlapping domains");

            while (counter < numberofdomains && totcount < 10 * numberofdomains) {
                if (Thread.currentThread().isInterrupted()) {
                    IJ.showStatus("Simulation Interrupted");
                    IJ.showProgress(1);
                    return;
                }
                domains[counter][0] = rugxpos.next();
                domains[counter][1] = rugypos.next();
                domains[counter][2] = rggdrad.next();
                for (int x = 0; x < counter; x++) {	// check that domains do not overlap; if there is overlap create a new domain
                    if (Math.pow(domains[counter][0] - domains[x][0], 2) + Math.pow(domains[counter][1] - domains[x][1], 2) < Math.pow(domains[counter][2] + domains[x][2], 2)) {
                        x = counter--;
                    }
                }
                counter++;
                totcount++;
                IJ.showProgress(counter, numberofdomains);
            }

            if (totcount >= 10 * numberofdomains) {
                IJ.showMessage("Domains too dense, cannot place them without overlap.");
                IJ.showStatus("Simulation Error");
                IJ.showProgress(1);
                return;
            }

            maxct = 0;
            for (int x = 1; x < numberofdomains; x++) {
                int xt = (int) Math.floor((domains[x][0] + gridMidPos) / subgridsize);
                int yt = (int) Math.floor((domains[x][1] + gridMidPos) / subgridsize);
                dsortmp[xt][yt][dctr[xt][yt]++] = x;		// dsortmp stores the number of each domain whose centre is in a particular grid area
                if (dctr[xt][yt] > maxct) {
                    maxct = dctr[xt][yt];	// maximum number of domains detected in any grid
                }
            }

            maxct *= 9;	// the domains of 9 neighbouring pixels are combined into one grid, so the maximum number increases accordingly

            domainsorted = new int[subgridnum][subgridnum][maxct];	// domains will be sorted into a grid for faster testing whether particles are in domains

            for (int x = 0; x < subgridnum; x++) {		// domainsorted contains for each grid position all domains in that and all directly surrounding grid positions
                for (int y = 0; y < subgridnum; y++) {	// as the grid is larger than a domain radius, a particle in that grid can be only in any of these domains if at all
                    int dct = 0;
                    int dtmp = 0;
                    while (dsortmp[x][y][dtmp] > 0) {
                        domainsorted[x][y][dct++] = dsortmp[x][y][dtmp++];
                    }
                    dtmp = 0;
                    if ((x + 1) < subgridnum) {
                        while (dsortmp[x + 1][y][dtmp] > 0) {
                            domainsorted[x][y][dct++] = dsortmp[x + 1][y][dtmp++];
                        }
                    }
                    dtmp = 0;
                    if ((y + 1) < subgridnum) {
                        while (dsortmp[x][y + 1][dtmp] > 0) {
                            domainsorted[x][y][dct++] = dsortmp[x][y + 1][dtmp++];
                        }
                    }
                    dtmp = 0;
                    if ((x + 1) < subgridnum && (y + 1) < subgridnum) {
                        while (dsortmp[x + 1][y + 1][dtmp] > 0) {
                            domainsorted[x][y][dct++] = dsortmp[x + 1][y + 1][dtmp++];
                        }
                    }
                    dtmp = 0;
                    if ((x - 1) >= 0) {
                        while (dsortmp[x - 1][y][dtmp] > 0) {
                            domainsorted[x][y][dct++] = dsortmp[x - 1][y][dtmp++];
                        }
                    }
                    dtmp = 0;
                    if ((y - 1) >= 0) {
                        while (dsortmp[x][y - 1][dtmp] > 0) {
                            domainsorted[x][y][dct++] = dsortmp[x][y - 1][dtmp++];
                        }
                    }
                    dtmp = 0;
                    if ((x - 1) >= 0 && (y - 1) >= 0) {
                        while (dsortmp[x - 1][y - 1][dtmp] > 0) {
                            domainsorted[x][y][dct++] = dsortmp[x - 1][y - 1][dtmp++];
                        }
                    }
                    dtmp = 0;
                    if ((x + 1) < subgridnum && (y - 1) >= 0) {
                        while (dsortmp[x + 1][y - 1][dtmp] > 0) {
                            domainsorted[x][y][dct++] = dsortmp[x + 1][y - 1][dtmp++];
                        }
                    }
                    dtmp = 0;
                    if ((x - 1) >= 0 && (y + 1) < subgridnum) {
                        while (dsortmp[x - 1][y + 1][dtmp] > 0) {
                            domainsorted[x][y][dct++] = dsortmp[x - 1][y + 1][dtmp++];
                        }
                    }
                }
            }
        }

        // determine intial positions for all particles
        for (int m = 0; m < simNoParticles; m++) {
            particles[m][0] = rugxpos.next();
            particles[m][1] = rugypos.next();
            particles[m][2] = 1.0;
            if (simBlinkFlag) {
                if ((int) ((m + 1) * simDarkF) > (int) (m * simDarkF)) {
                    particles[m][4] = 0.0;
                } else {
                    particles[m][4] = 1.0;
                }
            } else {
                particles[m][4] = 1.0;
            }
        }

        if (simDomainFlag) {	// check for each particle whether it is in a domain
            for (int m = 0; m < simNoParticles; m++) {
                particles[m][3] = simCheckInDomain(particles[m][0], particles[m][1], domains, domainsorted, maxct, gridMidPos, subgridsize);
            }
        }

        IJ.showStatus("Simulating ...");
        for (int n = 0; n < simNoTStep; n++) {	// run over all time steps/frames
            // check for interruption and stop excecution
            if (Thread.currentThread().isInterrupted()) {
                IJ.showStatus("Simulation Interrupted");
                IJ.showProgress(1);
                return;
            }

            ImageProcessor ipSim = impSim.getStack().getProcessor(n + 1); // get the image processor

            for (int dx = 0; dx < simPixelnum; dx++) { // add the camera offset and a noise term to each pixel
                for (int dy = 0; dy < simPixelnum; dy++) {
                    ipSim.putPixelValue(dx, dy, simCameraOffset + rggnoise.next());
                }
            }

            if (n == simBleachFrame) {	// if bleach frame is reached, bleach all particles within the bleach region
                for (int m = 0; m < simNoParticles; m++) {
                    if (Math.sqrt(Math.pow(particles[m][0], 2.0) + Math.pow(particles[m][1], 2.0)) < simBleachRadius) {
                        particles[m][2] = 0.0;
                    }
                }
            }

            for (int s = 0; s < simStepsPerFrame; s++) {		// there will be simStepsPerFrame for each frame simulated
                for (int m = 0; m < simNoParticles; m++) {	// change positions of all particles

                    double dx = 0;	// step sizes
                    double dy = 0;

                    if (!simDomainFlag && !simMeshFlag) {

                        if (m < num1) {	// if there are no domains and no mesh then diffuse freely
                            dx = rgg1.next();
                            dy = rgg1.next();
                        } else if (m >= num1 && m < num1 + num2) {
                            dx = rgg2.next();
                            dy = rgg2.next();
                        } else if (m >= num1 + num2) {
                            dx = rgg3.next();
                            dy = rgg3.next();
                        }

                    } else if (!simDomainFlag && simMeshFlag) {	// simulate diffusion on a simple meshwork grid

                        if (m < num1) {	// if there are no domains and no mesh then diffuse freely
                            dx = rgg1.next();
                            dy = rgg1.next();
                        } else if (m >= num1 && m < num1 + num2) {
                            dx = rgg2.next();
                            dy = rgg2.next();
                        } else if (m >= num1 + num2) {
                            dx = rgg3.next();
                            dy = rgg3.next();
                        }

                        boolean hoptrue = false;

                        if (simHopProbability > rugphop.next() || simHopProbability == 1) {
                            hoptrue = true;
                        }

                        if (!hoptrue) {	// if hop is not true, step inside the mesh only
                            while ((Math.floor(particles[m][0] / simMeshworkSize) != Math.floor((particles[m][0] + dx) / simMeshworkSize))) {
                                if (m < num1) {
                                    dx = rgg1.next();
                                } else if (m >= num1 && m < num1 + num2) {
                                    dx = rgg2.next();
                                } else if (m >= num1 + num2) {
                                    dx = rgg3.next();
                                }
                            }
                            while ((Math.floor(particles[m][1] / simMeshworkSize) != Math.floor((particles[m][1] + dy) / simMeshworkSize))) {
                                if (m < num1) {
                                    dy = rgg1.next();
                                } else if (m >= num1 && m < num1 + num2) {
                                    dy = rgg2.next();
                                } else if (m >= num1 + num2) {
                                    dy = rgg3.next();
                                }
                            }
                        }

                    } else if (simDomainFlag && !simMeshFlag) {	// if there are domains, determine for each particle whether it is in a domain and have it diffuse accordingly

                        int domnum = (int) particles[m][3];

                        if (domnum != 0) {			// if the particle is in a domain, the next step size is determined by domain diffusion
                            if (m < num1) {		// if there is no border crossing this would be the 1. case of diffusion within a domain (in->in)
                                dx = rggdom1.next();
                                dy = rggdom1.next();
                            } else if (m >= num1 && m < num1 + num2) {
                                dx = rggdom2.next();
                                dy = rggdom2.next();
                            } else if (m >= num1 + num2) {
                                dx = rggdom3.next();
                                dy = rggdom3.next();
                            }
                        } else {					// if the particle is not in a domain, the next step size is determined by diffusion of the surrounding matrix
                            if (m < num1) {		// if there is no border crossing this would be the 2. case of diffusion outside domains (out->out)
                                dx = rgg1.next();
                                dy = rgg1.next();
                            } else if (m >= num1 && m < num1 + num2) {
                                dx = rgg2.next();
                                dy = rgg2.next();
                            } else if (m >= num1 + num2) {
                                dx = rgg3.next();
                                dy = rgg3.next();
                            }
                        }

                        boolean crossinout = false;		// are crossings allowed?
                        boolean crossoutin = false;
                        if (simPout > rugpout.next() || simPout == 1.0) {
                            crossinout = true;
                        }
                        if (simPin > rugpin.next() || simPin == 1.0) {
                            crossoutin = true;
                        }

                        if (domnum != 0 && crossinout) {	// if inside domain and in-out allowed
                            int domcheck = simCheckInDomain(particles[m][0] + dx, particles[m][1] + dy, domains, domainsorted, maxct, gridMidPos, subgridsize); // is new position in domain
                            if (domcheck == 0) {	// act only if new position is actually outside domain
                                double domx = domains[domnum][0];	// domain coordinates
                                double domy = domains[domnum][1];
                                double domr = domains[domnum][2];
                                double px = particles[m][0] - domx;
                                double py = particles[m][1] - domy;
                                double sol = (-(px * dx + py * dy) + Math.sqrt(-Math.pow(py * dx - px * dy, 2.0) + Math.pow(dx * domr, 2.0) + Math.pow(dy * domr, 2.0))) / (Math.pow(dx, 2.0) + Math.pow(dy, 2.0));
                                dx = (sol + (1 - sol) * Math.sqrt(simDoutDinRatio)) * dx; // move the particle to the border with Din and outside with Dout
                                dy = (sol + (1 - sol) * Math.sqrt(simDoutDinRatio)) * dy;
                            }
                        }

                        if (domnum != 0 && !crossinout) {	// if inside domain and in-out not allowed
                            double domx = domains[domnum][0];	// domain coordinates
                            double domy = domains[domnum][1];
                            double domr = domains[domnum][2];
                            while (Math.pow(particles[m][0] + dx - domx, 2.0) + Math.pow(particles[m][1] + dy - domy, 2.0) > Math.pow(domr, 2.0)) {	// find (dx, dy) to stay in domain
                                if (m < num1) {
                                    dx = rggdom1.next();
                                    dy = rggdom1.next();
                                } else if (m >= num1 && m < num1 + num2) {
                                    dx = rggdom2.next();
                                    dy = rggdom2.next();
                                } else if (m >= num1 + num2) {
                                    dx = rggdom3.next();
                                    dy = rggdom3.next();
                                }
                            }
                        }

                        if (domnum == 0 && crossoutin) {	// if outside domain and out-in allowed
                            int domcheck = simCheckInDomain(particles[m][0] + dx, particles[m][1] + dy, domains, domainsorted, maxct, gridMidPos, subgridsize); // is new position in domain
                            if (domcheck != 0) { // act only if step brings particle into domain
                                double domx = domains[domcheck][0];	// domain coordinates
                                double domy = domains[domcheck][1];
                                double domr = domains[domcheck][2];
                                double px = particles[m][0] - domx;
                                double py = particles[m][1] - domy;
                                double sol = (-(px * dx + py * dy) - Math.sqrt(-Math.pow(py * dx - px * dy, 2.0) + Math.pow(dx * domr, 2.0) + Math.pow(dy * domr, 2.0))) / (Math.pow(dx, 2.0) + Math.pow(dy, 2.0));
                                dx = (sol + (1 - sol) / Math.sqrt(simDoutDinRatio)) * dx; // move the particle to the border with Dout and inside with Din
                                dy = (sol + (1 - sol) / Math.sqrt(simDoutDinRatio)) * dy;
                            }
                        }

                        if (domnum == 0 && !crossoutin) {	// if outside domain and out-in not allowed
                            while (simCheckInDomain(particles[m][0] + dx, particles[m][1] + dy, domains, domainsorted, maxct, gridMidPos, subgridsize) != 0) { // find (dx, dy) to stay outside domain
                                if (m < num1) {
                                    dx = rgg1.next();
                                    dy = rgg1.next();
                                } else if (m >= num1 && m < num1 + num2) {
                                    dx = rgg2.next();
                                    dy = rgg2.next();
                                } else if (m >= num1 + num2) {
                                    dx = rgg3.next();
                                    dy = rgg3.next();
                                }
                            }
                        }

                    } else if (simDomainFlag && simMeshFlag) {
                        IJ.showMessage("Mesh and Domain diffusion has not been implemented yet");
                        IJ.showStatus("Done");
                        return;
                    }

                    particles[m][0] += dx;	// finalize step
                    particles[m][1] += dy;

                    if (rugbf.next() > bleachFactor) {
                        particles[m][2] = 0.0;
                    }

                    int index = (int) particles[m][4];				//blinking of particles
                    if (simBlinkFlag) {
                        if (rugblink.next() > blinkFactor[index]) {
                            particles[m][4] = Math.abs(1.0 - particles[m][4]);
                        }
                    }

                    if (particles[m][0] > simSizeUL || particles[m][1] > simSizeUL || particles[m][0] < simSizeLL || particles[m][1] < simSizeLL) {
                        //Reset particle on border if particle left the simulation region
                        int tmp1 = (int) Math.floor(ruig.next() + 0.5);
                        int tmp2 = (int) (1 - 2 * Math.floor(ruig.next() + 0.5));
                        particles[m][0] = tmp1 * rugxpos.next() + (1 - tmp1) * tmp2 * simSizeUL;
                        particles[m][1] = (1 - tmp1) * rugypos.next() + tmp1 * tmp2 * simSizeUL;
                        particles[m][2] = 1.0;
                    }

                    if (simDomainFlag) {    // check the domain location of the particle
                        particles[m][3] = simCheckInDomain(particles[m][0], particles[m][1], domains, domainsorted, maxct, gridMidPos, subgridsize);
                    }
                    /*
                    if (m == tmpPartNum) {
                        IJ.log("{" + Double.toString(particles[m][0]) + ", " + Double.toString(particles[m][1]) + "},");	// provide a list of coordinates for a particular particle in the log window		
                    }
                     */
                    int nop = (int) Math.round(rpgphoton.next() * particles[m][2] * particles[m][4]); // create photons if the particle is fluorescent
                    for (int p = 0; p < nop; p++) { // run over emitted photons
                        double cordx = particles[m][0] + rggpsf.next();
                        double cordy = particles[m][1] + rggpsf.next();
                        if (cordx < simDetectorSize && cordy < simDetectorSize && cordx > -simDetectorSize && cordy > -simDetectorSize) {
                            int tpx = (int) Math.floor((cordx + simMidPos) / simPixelSize);
                            int tpy = (int) Math.floor((cordy + simMidPos) / simPixelSize);
                            int tmp = (int) (ipSim.getPixelValue(tpx, tpy) + 1);
                            ipSim.putPixelValue(tpx, tpy, tmp);
                        }
                    } // end photon loop (p)
                } // end particle loop (m)
            } // end step per frame loop (s)
            IJ.showProgress(n, simNoTStep);
        } // end frame loop (n)

        if (!batchSim) {
            // show the simulation file
            System.arraycopy(newSimSettings, 0, simSettings, 0, nosimsettings); // save the settings used for the simulation in simSettings
            imp = (ImagePlus) impSim.clone();
            imp.show();
            IJ.run(imp, "Enhance Contrast", "saturated=0.35");	//autoscaling the contrast 
            simFile = true;
            obtainImage();
            closeWindows();
        } else {
            // save the simulation file
            String $fs = simBatchPath.getAbsolutePath().toString() + "/" + "sim" + Double.toString(simD1 * Math.pow(10, 12)) + "-" + Double.toString(simD2 * Math.pow(10, 12)) + "-" + Double.toString(simF2) + ".tif";
            IJ.saveAsTiff(impSim, $fs);
        }
    }

    public int simCheckInDomain(double px, double py, double[][] domains, int[][][] domainsorted, int maxct, double gridMidPos, double subgridsize) {
        // px and py are particle positions
        // domains has the position and size of the different existing domains
        // domainsorted is essentially a grid in which the domains are sorted according to their place in the simulated area; this makes sure not all domains but only domians in the vicinity of the particle are searched
        // maxct is the maximum number of domains in any subgrid area (max length of domainsorted[][]) 
        // gridMidPos and subgridsize define the subgrid on which the domains were sorted
        int maxxt = domainsorted.length;
        int maxyt = domainsorted[0].length;
        int xt = (int) Math.floor((px + gridMidPos) / subgridsize);
        int yt = (int) Math.floor((py + gridMidPos) / subgridsize);
        int result = 0;
        if (xt > maxxt || xt < 0 || yt > maxyt || yt < 0) {
            return result;
        }
        boolean indomain = false;
        boolean onborder = false;
        int ct = 0;
        while (domainsorted[xt][yt][ct] > 0 && ct < maxct && !indomain) {	// check whether particle is in domain, and if yes, remember the domain number
            if (Math.pow(px - domains[domainsorted[xt][yt][ct]][0], 2.0) + Math.pow(py - domains[domainsorted[xt][yt][ct]][1], 2.0) - Math.pow(domains[domainsorted[xt][yt][ct]][2], 2.0) <= 0.0) {
                result = domainsorted[xt][yt][ct];	//remember number of domain in which particle resides
                indomain = true;
                if (Math.pow(px - domains[domainsorted[xt][yt][ct]][0], 2.0) + Math.pow(py - domains[domainsorted[xt][yt][ct]][1], 2.0) - Math.pow(domains[domainsorted[xt][yt][ct]][2], 2.0) == 0.0) {
                    onborder = true;
                }
            }
            ct++;
        }
        return result;
    }

    public void simulateACF3D() {
        int simSeed = 0;
        int simNoParticles = 1000; 						// number of simulated particles
        int simCPS = 10000;								// average count rate per particle per second
        double simTauBleach = 100000;					// bleach time in seconds
        int simPixelnum = 21; 							// width of image in pixels
        double simExtFactor = 1.5; 						// factor by which the simulated area is bigger than the observed area
        int simNoTStep = 50000; 						// number of frames to be simulated
        double simFrameTime = 0.001;					// time resolution of the camera in second
        int simStepsPerFrame = 10;						// simulation steps per frame
        double simD1 = 1.0 / Math.pow(10, 12);			// particle 1 diffusion coefficient
        double simD2 = 0.1 / Math.pow(10, 12);			// particle 2 diffusion coefficient
        double simD3 = 0.01 / Math.pow(10, 12);			// particle 3 diffusion coefficient
        double simF2 = 0.0;								// fraction of particle 2
        double simF3 = 0.0;								// fraction of particle 3
        double simKon = 1.0;							// on-rate for triplet
        double simKoff = 0.0;							// off-rate for triplet
        int simCameraOffset = 100;						// offset of CCD camera
        double simCameraNoiseFactor = 3.0;				// noise of CCD camera
        double simBleachRadius = 3.0;					// bleach radius
        int simBleachFrame = 10000000;					// frame at which bleach happens
        double simPixelSizeRS = 24 / Math.pow(10, 6);		// pixel size in real space
        double simMag = 63.0;							// objective maginification
        double simWavelength = 614.0 / Math.pow(10, 9);	// observation wavelegnth
        double simNA = 1.0; 							// NA of the objective
        double simSigma0 = 0.8;							// actual resolution
        double simSigmaZ = 1.0;							// axial resolution
        String[] newSimSettings = new String[nosimsettings];	// an array for reading out the settings in Simulation panel

        newSimSettings[0] = (String) cbSimMode.getSelectedItem();
        newSimSettings[1] = tbSimTrip.getText();
        newSimSettings[2] = tfSimSeed.getText();
        newSimSettings[3] = tfSimParticleNum.getText();
        newSimSettings[4] = tfSimCPS.getText();
        newSimSettings[5] = tfSimTauBleach.getText();
        newSimSettings[6] = tfSimPixelNum.getText();
        newSimSettings[7] = tfSimExtensionFactor.getText();
        newSimSettings[8] = tfSimTimeStepNum.getText();
        newSimSettings[9] = tfSimFrameTime.getText();
        newSimSettings[10] = tfSimStepsPerFrame.getText();
        newSimSettings[11] = tfSimCurrentStepSize.getText();
        newSimSettings[12] = tfSimD1.getText();
        newSimSettings[13] = tfSimDoutDinRatio.getText();
        newSimSettings[14] = tfSimD2.getText();
        newSimSettings[15] = tfSimF2.getText();
        newSimSettings[16] = tfSimD3.getText();
        newSimSettings[17] = tfSimF3.getText();
        newSimSettings[18] = tfSimKon.getText();
        newSimSettings[19] = tfSimKoff.getText();
        newSimSettings[20] = tfSimCameraOffset.getText();
        newSimSettings[21] = tfSimCameraNoiseFactor.getText();
        newSimSettings[22] = tfSimBleachRadius.getText();
        newSimSettings[23] = tfSimBleachFrame.getText();
        newSimSettings[24] = tfDomainRadius.getText();
        newSimSettings[25] = tfDomainDensity.getText();
        newSimSettings[26] = tfPin.getText();
        newSimSettings[27] = tfPout.getText();
        newSimSettings[28] = tfMeshworkSize.getText();
        newSimSettings[29] = tfHopProbability.getText();
        newSimSettings[30] = tfPixelSize.getText();
        newSimSettings[31] = tfMagnification.getText();
        newSimSettings[32] = tfEmLambda.getText();
        newSimSettings[33] = tfNA.getText();
        newSimSettings[34] = tfSigma.getText();
        newSimSettings[35] = tfSigmaZ.getText();

        try {
            simSeed = Integer.parseInt(newSimSettings[2]);
            simNoParticles = Integer.parseInt(newSimSettings[3]); 					// number of simulated particles
            simCPS = Integer.parseInt(newSimSettings[4]); 							// average count rate per particle per second
            simTauBleach = Double.parseDouble(newSimSettings[5]);					// bleach time in seconds
            simPixelnum = Integer.parseInt(newSimSettings[6]); 						// width of image in pixels
            simExtFactor = Double.parseDouble(newSimSettings[7]); 					// factor by which the simulated area is bigger than the observed area
            simNoTStep = Integer.parseInt(newSimSettings[8]); 						// number of frames to be simulated
            simFrameTime = Double.parseDouble(newSimSettings[9]);					// time resolution of the camera in second
            simStepsPerFrame = Integer.parseInt(newSimSettings[10]);				// steps of simulations doen for each frame
            simD1 = Double.parseDouble(newSimSettings[12]) / Math.pow(10, 12);		// particle 1 diffusion coefficient		
            simD2 = Double.parseDouble(newSimSettings[14]) / Math.pow(10, 12);		// particle 2 diffusion coefficient
            simD3 = Double.parseDouble(newSimSettings[16]) / Math.pow(10, 12);		// particle 3 diffusion coefficient
            simF2 = Double.parseDouble(newSimSettings[15]);							// fraction of particle 2
            simF3 = Double.parseDouble(newSimSettings[17]);							// fraction of particle 
            simKon = Double.parseDouble(newSimSettings[18]);						// on-rate for triplet
            simKoff = Double.parseDouble(newSimSettings[19]);						// off-rate for triplet
            simCameraOffset = Integer.parseInt(newSimSettings[20]);					// offset of CCD camera
            simCameraNoiseFactor = Integer.parseInt(newSimSettings[21]);			// noise of CCD camera
            simBleachRadius = Double.parseDouble(newSimSettings[22]) / Math.pow(10, 6);		// bleach radius
            simBleachFrame = Integer.parseInt(newSimSettings[23]);							// frame at which bleach happens
            simPixelSizeRS = Double.parseDouble(newSimSettings[30]) / Math.pow(10, 6);		// pixel size in real space
            simMag = Double.parseDouble(newSimSettings[31]);								// objective maginification
            simWavelength = Double.parseDouble(newSimSettings[32]) / Math.pow(10, 9);			// observation wavelegnth
            simNA = Double.parseDouble(newSimSettings[33]); 								// NA of the objective
            simSigma0 = Double.parseDouble(newSimSettings[34]);								// actual resolution
            simSigmaZ = Double.parseDouble(newSimSettings[35]);
        } catch (NumberFormatException nfe) {
            IJ.showMessage("One of the values in the simulation window does not have the right format (integer or double).");
            throw new NumberFormatException("Number format error.");
        }

        double simTStep = simFrameTime / simStepsPerFrame;
        double simDarkF = simKoff / (simKoff + simKon);		//fraction of molecules in the dark state

        int numOfSeeds = 50;
        int[] simSeedArray = new int[numOfSeeds];							// array of simulation Seeds so that the random number generators are different
        double simPixelSize = simPixelSizeRS / simMag;					// pixel size in object space
        double simPSFSize = 0.5 * simSigma0 * simWavelength / simNA;		// PSF size
        double simGridSize = simPixelnum * simPixelSize;				// gridsize
        double simMidPos = simGridSize / 2.0;  								// middel position
        double simSizeLL = -simExtFactor * simGridSize; 				// lower limit of the simulation area
        double simSizeUL = simExtFactor * simGridSize; 					// upper limit of the simulation area
        double simDetectorSize = simGridSize / 2.0; 						// size of the observation areas
        double simPhotonsPerStep = Math.floor(simCPS * simTStep + 0.5); // number of photons per particle and time step
        double bleachFactor = 2.0;									// the 2.0 ensures that no bleaching happens
        if (simTauBleach != 0) {
            bleachFactor = Math.exp(-simTStep / simTauBleach);
        }
        double[] blinkFactor = new double[2];
        blinkFactor[0] = Math.exp(-simTStep * simKon);
        blinkFactor[1] = Math.exp(-simTStep * simKoff);

        double lightSheetThickness = simSigmaZ * simWavelength / simNA / 2.0; // division by 2 to yield the 1/sqrt(e) radius

        double simThicknessLL = -10.0 * lightSheetThickness;			// lower and upper limit of the zdirection of the simulation volume
        double simThicknessUL = 10.0 * lightSheetThickness;			// z dimension is 20 times the light sheet thickness

        // divide particle in to their types according to their fractions (1- F2 - F3), F2, F3
        int num1 = (int) Math.ceil(simNoParticles * (1 - simF2 - simF3));
        int num2 = (int) Math.ceil(simNoParticles * simF2);
        int num3 = simNoParticles - num1 - num2;
        double[][] particles = new double[simNoParticles][5];	// array for particle positions and whether particle is bleached or in the dark state
        double zcor;
        //	double zfac = Math.tan( Math.asin(simNA/1.333) );			// factor describing the spread of the PSF cross-section on the camera if the particle is not in the focal plane
        double zfac = simNA / Math.sqrt(1.776889 - Math.pow(simNA, 2));
        ImagePlus impSim = IJ.createImage("3D Simulation", "GRAY16", simPixelnum, simPixelnum, simNoTStep);

        if (simSeed == 0) {
            Arrays.fill(simSeedArray, 0);
        } else {
            for (int x = 0; x < numOfSeeds; x++) {
                simSeedArray[x] = simSeed + (int) Math.pow(x, 2.0);
            }
        }

        int cs = 0;
        UniformGenerator rugxpos = new UniformGenerator(simSizeLL, simSizeUL, simSeedArray[cs++]);
        UniformGenerator rugypos = new UniformGenerator(simSizeLL, simSizeUL, simSeedArray[cs++]);
        UniformGenerator rugzpos = new UniformGenerator(simThicknessLL, simThicknessUL, simSeedArray[cs++]);
        UniformGenerator ruig = new UniformGenerator(simSeedArray[cs++]);
        UniformGenerator rugblink = new UniformGenerator(simSeedArray[cs++]);
        GaussianGenerator rgg1 = new GaussianGenerator(0, Math.sqrt(2 * simD1 * simTStep), simSeedArray[cs++]);
        GaussianGenerator rgg2 = new GaussianGenerator(0, Math.sqrt(2 * simD2 * simTStep), simSeedArray[cs++]);
        GaussianGenerator rgg3 = new GaussianGenerator(0, Math.sqrt(2 * simD3 * simTStep), simSeedArray[cs++]);
        GaussianGenerator rggpsf = new GaussianGenerator(0, simPSFSize, simSeedArray[cs++]);
        PoissonGenerator rpgphoton = new PoissonGenerator(simTStep * simCPS, simSeedArray[cs++]);
        GaussianGenerator rggnoise = new GaussianGenerator(0, Math.sqrt(simCameraNoiseFactor), simSeedArray[cs++]);
        UniformGenerator dug1 = new UniformGenerator(0, 3, simSeedArray[cs++]);
        UniformGenerator dug2 = new UniformGenerator(0, 2, simSeedArray[cs++]);
        UniformGenerator rugbf = new UniformGenerator(0, 1, simSeedArray[cs++]);
        UniformGenerator BMU1 = new UniformGenerator(0, 1, simSeedArray[cs++]);
        UniformGenerator BMU2 = new UniformGenerator(0, 1, simSeedArray[cs++]);
        UniformGenerator BMU3 = new UniformGenerator(0, 1, simSeedArray[cs++]);
        UniformGenerator BMU4 = new UniformGenerator(0, 1, simSeedArray[cs++]);
        // determine intial positions for all particles
        for (int m = 0; m < simNoParticles; m++) {
            particles[m][0] = rugxpos.next();
            particles[m][1] = rugypos.next();
            //	particles[m][0] = 0;
            //	particles[m][1] = 0;
            particles[m][2] = rugzpos.next();
            //	particles[m][2] = 2*lightSheetThickness;
            particles[m][3] = 1;
            if (simBlinkFlag) {
                if ((int) ((m + 1) * simDarkF) > (int) (m * simDarkF)) {
                    particles[m][4] = 0.0;
                } else {
                    particles[m][4] = 1.0;
                }
            } else {
                particles[m][4] = 1.0;
            }
        }

        for (int n = 0; n < simNoTStep; n++) {	// run over all time steps/frames

            // check for interruption and stop excecution
            if (Thread.currentThread().isInterrupted()) {
                IJ.showStatus("Simulation Interrupted");
                IJ.showProgress(1);
                return;
            }

            ImageProcessor ipSim = impSim.getStack().getProcessor(n + 1);

            for (int dx = 0; dx < simPixelnum; dx++) { // add the camera offset and a noise term to each pixel
                for (int dy = 0; dy < simPixelnum; dy++) {
                    ipSim.putPixelValue(dx, dy, simCameraOffset + rggnoise.next());
                }
            }
            for (int s = 0; s < simStepsPerFrame; s++) {
                for (int m = 0; m < simNoParticles; m++) {	// change positions of all particles for each time step
                    int numOfSeeds1 = m + 1;
                    int[] simSeedArray1 = new int[numOfSeeds1];
                    int ks = 0;
                    if (m < num1) {
                        particles[m][0] += rgg1.next();
                        particles[m][1] += rgg1.next();
                        //	particles[m][0] = 0;
                        //	particles[m][1] = 0;
                        particles[m][2] += rgg1.next();
                        //	particles[m][2] = 2*lightSheetThickness;
                    } else if (m >= num1 && m < num1 + num2) {
                        particles[m][0] += rgg2.next();
                        particles[m][1] += rgg2.next();
                        particles[m][2] += rgg2.next();
                    } else if (m >= num1 + num2) {
                        particles[m][0] += rgg3.next();
                        particles[m][1] += rgg3.next();
                        particles[m][2] += rgg3.next();
                    }

                    if (particles[m][3] != 0.0) {				// bleaching of particles
                        if (rugbf.next() > bleachFactor) {
                            particles[m][3] = 0.0;
                        }
                    }

                    int index = (int) particles[m][4];			//blinking of particles
                    if (simBlinkFlag) {
                        if (rugblink.next() > blinkFactor[index]) {
                            particles[m][4] = Math.abs(1.0 - particles[m][4]);
                        }
                    }

                    if (particles[m][0] > simSizeUL || particles[m][1] > simSizeUL || particles[m][2] > simThicknessUL || particles[m][0] < simSizeLL || particles[m][1] < simSizeLL || particles[m][2] < simThicknessLL) {
                        //Reset particle on border
                        int tmp1 = (int) Math.ceil(dug1.next());
                        int tmp2 = (int) Math.ceil(dug2.next());
                        if (tmp1 == 1) {
                            particles[m][0] = rugxpos.next();
                            particles[m][1] = rugypos.next();
                            if (tmp2 == 1) {
                                particles[m][2] = simThicknessLL;
                            } else {
                                particles[m][2] = simThicknessUL;
                            }
                        } else if (tmp1 == 2) {
                            particles[m][0] = rugxpos.next();
                            particles[m][2] = rugzpos.next();
                            if (tmp2 == 1) {
                                particles[m][1] = simSizeLL;
                            } else {
                                particles[m][1] = simSizeUL;
                            }
                        } else {
                            particles[m][1] = rugypos.next();
                            particles[m][2] = rugzpos.next();
                            if (tmp2 == 1) {
                                particles[m][0] = simSizeLL;
                            } else {
                                particles[m][0] = simSizeUL;
                            }
                        }
                        particles[m][3] = 1.0; // * new particle is fluorescent
                    }
                    zcor = (simPSFSize + (Math.abs(particles[m][2]) * (zfac / 2))); // factor describing the increase of the PSF at the focal plane for a particle not situated in the focal plane
                    //	zcor=simPSFSize;
                    //	double photoncount=simTStep * simCPS* Math.exp(-0.5 * Math.pow(particles[m][2] / lightSheetThickness, 2.0));
                    //	PoissonGenerator rpgphoton = new PoissonGenerator(photoncount, simSeedArray1[ks++]);
                    int nop = (int) Math.round((Math.abs(rpgphoton.next() * Math.exp(-0.5 * Math.pow(particles[m][2] / lightSheetThickness, 2.0)) * particles[m][3] * particles[m][4])));
                    //	int nop = (int)((Math.round((Math.abs(rpgphoton.next())))) * particles[m][3] * particles[m][4]);
                    //	int nop= (int)nop1;
                    if (nop < 0) {						// no negative photons
                        nop = 0;
                    }
                    /*
                    if (m == tmpPartNum) {
                        IJ.log(Double.toString(particles[m][2]));	// provide a list of coordinates for a particular particle in the log window		
                    }
                     */
                    for (int p = 0; p < nop; p++) { // run over emitted photons
                        //	double cordx = particles[m][0] + (rggpsf.next() );
                        //	double cordy = particles[m][1] + (rggpsf.next() );
                        double cordx = particles[m][0] + zcor * (Math.sqrt(-2 * Math.log(BMU1.next()))) * Math.cos(2 * 3.1415 * BMU2.next());
                        double cordy = particles[m][1] + zcor * (Math.sqrt(-2 * Math.log(BMU3.next()))) * Math.cos(2 * 3.1415 * BMU4.next());
                        if (cordx < simDetectorSize && cordy < simDetectorSize && cordx > -simDetectorSize && cordy > -simDetectorSize) {
                            int tpx = (int) Math.floor((cordx + simMidPos) / simPixelSize);
                            int tpy = (int) Math.floor((cordy + simMidPos) / simPixelSize);
                            int tmp = (int) (ipSim.getPixelValue(tpx, tpy) + 1);
                            ipSim.putPixelValue(tpx, tpy, tmp);
                        }
                    } // end photon loop (p)
                } // end particle loop (m)
            } // end steps per frame loop (s)
            IJ.showProgress(n, simNoTStep);
        } // end frame loop (n)
        System.arraycopy(newSimSettings, 0, simSettings, 0, nosimsettings); // save the settings used for the simulation in simSettings

        imp = (ImagePlus) impSim.clone();
        imp.show();
        IJ.run(imp, "Enhance Contrast", "saturated=0.35");	//autoscaling the contrast 
        simFile = true;
        obtainImage();
    }

    /* Miscellaneous functions
 *  
 *  public double determinant(double A[][], int N): calculate determinant
 *  public double[][] generateSubArray (double A[][], int N, int j1): generate a subarray 
 *  
     */
    // calculate determinant of a subarray; required for the Bayes model probabilities
    // this was adapted from http://stackoverflow.com/questions/16602350/calculating-matrix-determinant
    public double determinant(double A[][], int N) {
        double res;
        double[][] m;

        // 1x1 matrix
        if (N == 1) {
            res = A[0][0];
        } // 2x2 matrix
        else if (N == 2) {
            res = A[0][0] * A[1][1] - A[1][0] * A[0][1];
        } // NxN matrix
        else {
            res = 0;
            for (int j1 = 0; j1 < N; j1++) {
                m = generateSubArray(A, N, j1);
                res += Math.pow(-1.0, 1.0 + j1 + 1.0) * A[0][j1] * determinant(m, N - 1);
            }
        }
        return res;
    }

    public double[][] generateSubArray(double A[][], int N, int j1) {
        double[][] m = new double[N - 1][];
        for (int k = 0; k < (N - 1); k++) {
            m[k] = new double[N - 1];
        }
        for (int i = 1; i < N; i++) {
            int j2 = 0;
            for (int j = 0; j < N; j++) {
                if (j == j1) {
                    continue;
                }
                m[i - 1][j2] = A[i][j];
                j2++;
            }
        }
        return m;
    }

    private void getVoxels(int this_width, int this_height, int this_framediff, float[] pixels, boolean IsNBCalculation) {
        int startX = IsNBCalculation ? 0 : roi1StartX;
        int startY = IsNBCalculation ? 0 : roi1StartY;
        int startFrame = firstframe - 1;

        try {
            imp.getStack().getVoxels(startX, startY, startFrame, this_width, this_height, this_framediff, pixels);
        } catch (Exception e) {
            for (int t = startFrame; t < this_framediff; t++) {
                for (int j = startY; j < this_height; j++) {
                    for (int i = startX; i < this_width; i++) {
                        pixels[t * this_height * this_width + j * this_width + i] = (float) imp.getStack().getVoxel(i, j, t);
                    }
                }
            }
        }
    }

    /*
    arrays pixels1, blockvararray, NBmeanGPU, and NBcovarianceGPU are passed by reference.
     */
    private boolean GPU_Calculate_ACF(float pixels[], double pixels1[], double blockvararray[],
            double NBmeanGPU[], double NBcovarianceGPU[],
            double bleachcorr_params[], ACFParameters GPUparams) {

        boolean IsGPUCalculationOK = true;

        double[] blocked1D = new double[GPUparams.width * GPUparams.height * GPUparams.chanum];

        // run GPU code
        try {
            GpufitImFCS.calcACF(pixels, pixels1, blockvararray, NBmeanGPU, NBcovarianceGPU, blocked1D, bleachcorr_params, samp, lag, GPUparams);
        } catch (Exception e) {
            IsGPUCalculationOK = false;
            e.printStackTrace(System.out);
        }

        //reset blocked
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                blocked[0][i][j] = 0;
            }
        }

        // copy result to blocked
        for (int y = 0; y < GPUparams.height; y++) {
            for (int x = 0; x < GPUparams.width; x++) {
                blocked[0][x + roi1StartX][y + roi1StartY] = (int) blocked1D[GPUparams.width * GPUparams.height + y * GPUparams.width + x];
            }
        }
        return IsGPUCalculationOK;
    }

    private boolean GPU_Initialize_GPUparams(GpufitImFCS.ACFParameters GPUparams, boolean isNBcalculation) {
        //pixbinX and pixbinY in the cpu plugin are strides in width and height respectively. They  are determined automatically depending on the overlap in the plugin.
        int w_out = isNBcalculation ? width : (int) Math.floor((Math.min(roi1WidthX + cfXDistance, width) - cfXDistance) / pixbinX);
        int h_out = isNBcalculation ? height : (int) Math.floor((Math.min(roi1HeightY + cfYDistance, height) - cfYDistance) / pixbinY);

        if (w_out <= 0 || h_out <= 0) {
            return false;
        } else {
            // win_star and hin_star determines the necessary area/pixels of all required data for GPU calculations. See function getVoxels too.
            int win_star = w_out * pixbinX + cfXDistance;
            int hin_star = h_out * pixbinY + cfYDistance;

            // w_temp and h_temp are 'intermediate' dimensions after accounting for binning.
            // example, if binning is not 1x1, w_temp x h_temp is an area smaller than win_star x hin_star.
            // it is essentially binning with stride 1 (even if overlap is false), 
            // so that we have all the necessary data for all scenarios of our GPU calculations. It simplifies/generalizes the indexing task in our CUDA kernels. 
            // else if binning is 1x1, w_temp = win_star and h_temp = hin_star.
            int w_temp = win_star - binningX + 1;
            int h_temp = hin_star - binningY + 1;

            // stores the current/actual image width and height
            // perform single or double exponential bleach corrections if selected
            String $bcmode = (String) cbBleachCor.getSelectedItem();

            if ("Polynomial".equals($bcmode)) {
                GPUparams.bleachcorr_gpu = true;
                GPUparams.bleachcorr_order = polyOrder + 1;
            } else {
                GPUparams.bleachcorr_gpu = false;
                GPUparams.bleachcorr_order = 0;
            }

            GPUparams.width = w_out;  //output width
            GPUparams.height = h_out;  //output height
            GPUparams.win_star = isNBcalculation ? width : win_star;
            GPUparams.hin_star = isNBcalculation ? height : hin_star;
            GPUparams.w_temp = isNBcalculation ? width : w_temp; // win_star - pixbinX + 1
            GPUparams.h_temp = isNBcalculation ? height : h_temp; // hin_start - pixbinY + 1
            GPUparams.pixbinX = isNBcalculation ? 1 : pixbinX;
            GPUparams.pixbinY = isNBcalculation ? 1 : pixbinY;
            GPUparams.binningX = isNBcalculation ? 1 : binningX; // binning in X axis
            GPUparams.binningY = isNBcalculation ? 1 : binningY; // binning in Y axis
            GPUparams.firstframe = firstframe;
            GPUparams.lastframe = lastframe;
            GPUparams.framediff = lastframe - firstframe + 1;
            GPUparams.cfXDistance = isNBcalculation ? 0 : cfXDistance;
            GPUparams.cfYDistance = isNBcalculation ? 0 : cfYDistance;
            GPUparams.correlatorp = correlatorp;
            GPUparams.correlatorq = correlatorq;
            GPUparams.frametime = frametime;
            GPUparams.background = background;
            GPUparams.mtab1 = mtab[1]; // mtab[1], used to calculate blocknumgpu.
            GPUparams.mtabchanumminus1 = mtab[chanum - 1]; // mtab[chanum-1], used to calculate pnumgpu[counter_indexarray]
            GPUparams.sampchanumminus1 = samp[chanum - 1]; // samp[chanum-1], used to calculate pnumgpu[counter_indexarray]
            GPUparams.chanum = chanum;
            GPUparams.isNBcalculation = isNBcalculation; // true;
            GPUparams.nopit = nopit;
            GPUparams.ave = (int) Math.floor((lastframe - firstframe + 1) / nopit);

            return true;
        }

    }

    private boolean GPU_get_pixels(GpufitImFCS.ACFParameters GPUparams, float[] pixels, boolean isNBcalculation) {
        // pixels is the input intensity array on which the auto and cross-correlation will be calculated.

        boolean IsGPUCalculationOK = true;

        if (GPUparams.binningY == 1 && GPUparams.binningX == 1) {
            try {
                getVoxels(GPUparams.w_temp, GPUparams.h_temp, GPUparams.framediff, pixels, isNBcalculation);
            } catch (Exception e) {
                IsGPUCalculationOK = false;
                e.printStackTrace(System.out);
            }

        } else {
            float[] pixels_approach2 = new float[GPUparams.win_star * GPUparams.hin_star * GPUparams.framediff];
            getVoxels(GPUparams.win_star, GPUparams.hin_star, GPUparams.framediff, pixels_approach2, isNBcalculation);

            // We found that the JNI function SetFloatArrayRegion fails when the output array is too huge. Tentatively, we set a limit of 96*96*50000.
            boolean WithinSizeLimit = (GPUparams.w_temp * GPUparams.h_temp * GPUparams.framediff) < 96 * 96 * 50000;

            if (WithinSizeLimit && GpufitImFCS.isBinningMemorySufficient(GPUparams)) {
                // Binning on GPU                
                try {
                    GpufitImFCS.calcBinning(pixels_approach2, pixels, GPUparams);
                } catch (Exception e) {
                    IsGPUCalculationOK = false;
                    e.printStackTrace(System.out);
                }
            } else {
                // Binning on CPU
                try {
                    float sum;
                    for (int z = 0; z < GPUparams.framediff; z++) {
                        for (int y = 0; y < GPUparams.h_temp; y++) {
                            for (int x = 0; x < GPUparams.w_temp; x++) {
                                sum = 0;
                                for (int k = 0; k < GPUparams.binningY; k++) {
                                    for (int i = 0; i < GPUparams.binningX; i++) {
                                        sum += (float) pixels_approach2[z * GPUparams.win_star * GPUparams.hin_star + (y + k) * GPUparams.win_star + x + i];
                                    }
                                }
                                pixels[z * GPUparams.w_temp * GPUparams.h_temp + y * GPUparams.w_temp + x] = sum;
                            }
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace(System.out);
                }
            }
        }

        if (bgrloaded) {
            float[] bgrmeangpu = new float[GPUparams.w_temp * GPUparams.h_temp];

            for (int j = 0; j < GPUparams.h_temp; j++) {
                for (int i = 0; i < GPUparams.w_temp; i++) {
                    bgrmeangpu[j * GPUparams.w_temp + i] = (float) (int) Math.round(bgrmean[i][j]);
                    //bgrmeangpu[j * GPUparams.w_temp + i] = (float) bgrmean[i][j];
                    //bckg = (int) Math.round(bgrmean[px + i][py + k]);
                }

            }

            for (int k = 0; k < GPUparams.framediff; k++) {
                for (int j = 0; j < GPUparams.h_temp; j++) {
                    for (int i = 0; i < GPUparams.w_temp; i++) {
                        pixels[k * GPUparams.h_temp * GPUparams.w_temp + j * GPUparams.w_temp + i] -= bgrmeangpu[j * GPUparams.w_temp + i];
                    }
                }
            }
        } else {
            int background_correction = background * GPUparams.binningX * GPUparams.binningY;
            for (int i = 0; i < GPUparams.w_temp * GPUparams.h_temp * GPUparams.framediff; i++) {
                pixels[i] -= background_correction;
            }
        }

        return IsGPUCalculationOK;
    }

    private boolean GPU_Calculate_BleachCorrection(GpufitImFCS.ACFParameters GPUparams, float[] pixels, double[] bleachcorr_params) {
        boolean IsGPUCalculationOK = true;

        try {
            float[] datableach_correction = new float[GPUparams.w_temp * GPUparams.h_temp * nopit];
            GpufitImFCS.calcDataBleachCorrection(pixels, datableach_correction, GPUparams);
            int numberFitsbleach_correction = GPUparams.w_temp * GPUparams.h_temp, numberPoints_bleach_correction = nopit;
            float tolerance_bleachcorrection = 0.0000000000000001f; // similar to GPU_ACF_Fit
            int maxNumberIterations_bleachcorrection = fitMaxIterations;
            Model model_bleachcorrection = Model.LINEAR_1D;
            Estimator estimator_bleachcorrection = Estimator.LSE;

            Boolean[] parameters_to_fit_bleachcorrection = new Boolean[model_bleachcorrection.numberParameters];
            parameters_to_fit_bleachcorrection[0] = true;
            for (int i = 1; i < model_bleachcorrection.numberParameters; i++) {
                parameters_to_fit_bleachcorrection[i] = i < GPUparams.bleachcorr_order;
            }

            // NOTE: initialization of 0-th term is different from CPU code, where in class PolynomFit (extends AbstractCurveFitter), the last point in target array is used as offset estimate.
            float[] initialParameters_bleachcorrection = new float[numberFitsbleach_correction * model_bleachcorrection.numberParameters];
            for (int i = 0; i < numberFitsbleach_correction; i++) {
                int offset = i * model_bleachcorrection.numberParameters;
                for (int j = 0; j < model_bleachcorrection.numberParameters; j++) {
                    initialParameters_bleachcorrection[offset + j] = (j == 0) ? datableach_correction[(i + 1) * nopit - 1] : (float) 0; // last value of every nopit points. This works the best.
                }
            }

            float[] intTime_bleachcorrection = new float[nopit];
            for (int z1 = 0; z1 < nopit; z1++) {
                intTime_bleachcorrection[z1] = (float) (frametime * (z1 + 0.5) * GPUparams.ave);
            }

            float[] weights_bleachcorr = new float[numberFitsbleach_correction * numberPoints_bleach_correction];
            for (int i = 0; i < numberFitsbleach_correction * numberPoints_bleach_correction; i++) {
                weights_bleachcorr[i] = (float) 1.0;
            }

            FitModel fitModel_bleachcorrection = new FitModel(numberFitsbleach_correction, numberPoints_bleach_correction, true, model_bleachcorrection, tolerance_bleachcorrection, maxNumberIterations_bleachcorrection, GPUparams.bleachcorr_order, parameters_to_fit_bleachcorrection, estimator_bleachcorrection, nopit * Float.SIZE / 8);
            fitModel_bleachcorrection.weights.clear();
            fitModel_bleachcorrection.weights.put(weights_bleachcorr);
            fitModel_bleachcorrection.userInfo.clear();
            fitModel_bleachcorrection.userInfo.put(intTime_bleachcorrection);
            fitModel_bleachcorrection.data.clear();
            fitModel_bleachcorrection.data.put(datableach_correction);
            fitModel_bleachcorrection.initialParameters.clear();
            fitModel_bleachcorrection.initialParameters.put(initialParameters_bleachcorrection);
            FitResult fitResult_bleachcorrection = GpufitImFCS.fit(fitModel_bleachcorrection);
            // boolean[] converged_bleachcorrection = new boolean[numberFitsbleach_correction];

            int numberConverged_bleachcorrection = 0, numberMaxIterationExceeded_bleachcorrection = 0, numberSingularHessian_bleachcorrection = 0, numberNegativeCurvatureMLE_bleachcorrection = 0;
            for (int i = 0; i < numberFitsbleach_correction; i++) {
                FitState fitState_bleachcorrection = FitState.fromID(fitResult_bleachcorrection.states.get(i));
                // converged_bleachcorrection[i] = fitState_bleachcorrection.equals(FitState.CONVERGED);
                switch (fitState_bleachcorrection) {
                    case CONVERGED:
                        numberConverged_bleachcorrection++;
                        break;
                    case MAX_ITERATIONS:
                        numberMaxIterationExceeded_bleachcorrection++;
                        break;
                    case SINGULAR_HESSIAN:
                        numberSingularHessian_bleachcorrection++;
                        break;
                    case NEG_CURVATURE_MLE:
                        numberNegativeCurvatureMLE_bleachcorrection++;
                }
            }

            int counter4 = 0;
            for (int y1 = 0; y1 < GPUparams.h_temp; y1++) {
                for (int x1 = 0; x1 < GPUparams.w_temp; x1++) {
                    for (int ii = 0; ii < GPUparams.bleachcorr_order; ii++) {
                        bleachcorr_params[counter4] = fitResult_bleachcorrection.parameters.get((y1 * GPUparams.w_temp + x1) * model_bleachcorrection.numberParameters + ii);
                        counter4 += 1;
                    }
                }
            }

            fitModel_bleachcorrection.reset();
            fitResult_bleachcorrection.reset();
        } catch (Exception e) {
            IsGPUCalculationOK = false;
            e.printStackTrace(System.out);
        }

        return IsGPUCalculationOK;
    }

    private boolean GPU_ACF_Fit(GpufitImFCS.ACFParameters GPUparams, double[] pixels1, double[] blockvararray) {
        boolean IsGPUCalculationOK = true;

        try {
            prepareFit();

            int numberFits = GPUparams.width * GPUparams.height;
            int numberPoints = chanum - 1;
            float tolerance = 0.0000000000000001f;
            int maxNumberIterations = fitMaxIterations;
            Model model = Model.ACF_1D;
            Estimator estimator = Estimator.LSE;
            float[] trueParameters;

            trueParameters = new float[]{
                (float) Double.parseDouble(tfParamN.getText()),
                (float) Double.parseDouble(tfParamD.getText()) / (float) Math.pow(10, 12),
                (float) Double.parseDouble(tfParamVx.getText()) / (float) Math.pow(10, 6),
                (float) Double.parseDouble(tfParamVy.getText()) / (float) Math.pow(10, 6),
                (float) Double.parseDouble(tfParamG.getText()),
                (float) Double.parseDouble(tfParamF2.getText()),
                (float) Double.parseDouble(tfParamD2.getText()) / (float) Math.pow(10, 12),
                (float) Double.parseDouble(tfParamF3.getText()),
                (float) Double.parseDouble(tfParamD3.getText()) / (float) Math.pow(10, 12),
                (float) Double.parseDouble(tfParamFtrip.getText()),
                (float) Double.parseDouble(tfParamTtrip.getText()) / (float) Math.pow(10, 6),
                (float) pixeldimx,
                (float) pixeldimy,
                (float) psfsize,
                (float) lsthickness,
                (float) pixeldimx * cfXshift / binningX,
                (float) pixeldimy * cfYshift / binningY,
                (float) fitobsvol,
                (float) Double.parseDouble(tfParamQ2.getText()),
                (float) Double.parseDouble(tfParamQ3.getText())};

            float[] initialParameters = new float[numberFits * model.numberParameters];

            for (int i = 0; i < numberFits; i++) {
                int offset = i * model.numberParameters;
                // System.arraycopy(trueParameters, 0, initialParameters, offset, model.numberParameters);
                for (int j = 0; j < model.numberParameters; j++) {
                    initialParameters[offset + j] = trueParameters[j];
                }
            }

            float[] user_info = new float[numberPoints];
            for (int i = 0; i < GPUparams.chanum - 1; i++) {
                user_info[i] = (float) lagtime[i + 1];
            }

            float[] data = new float[numberFits * numberPoints];
            int counter = 0;

            // Out of the 20 parameters passed to the GPUfit, only the first 11 are fitting parameters.
            Boolean[] parameters_to_fit = new Boolean[model.numberParameters];
            parameters_to_fit[0] = !rbtnHoldN.isSelected();
            parameters_to_fit[1] = !rbtnHoldD.isSelected();
            parameters_to_fit[2] = !rbtnHoldVx.isSelected();
            parameters_to_fit[3] = !rbtnHoldVy.isSelected();
            parameters_to_fit[4] = !rbtnHoldG.isSelected();
            parameters_to_fit[5] = !rbtnHoldF2.isSelected();
            parameters_to_fit[6] = !rbtnHoldD2.isSelected();
            parameters_to_fit[7] = !rbtnHoldF3.isSelected();
            parameters_to_fit[8] = !rbtnHoldD3.isSelected();
            parameters_to_fit[9] = !rbtnHoldFtrip.isSelected();
            parameters_to_fit[10] = !rbtnHoldTtrip.isSelected();
            for (int i = 11; i < model.numberParameters; i++) {
                parameters_to_fit[i] = false;
            }

            for (int y = 0; y < GPUparams.height; y++) {
                for (int x = 0; x < GPUparams.width; x++) {
                    for (int z = 0; z < GPUparams.chanum - 1; z++) {
                        data[counter] = (float) pixels1[(z + 1) * GPUparams.width * GPUparams.height + y * GPUparams.width + x];
                        counter += 1;
                    }
                }
            }

            float[] weights = new float[numberFits * numberPoints];
            int counter1 = 0;
            for (int y = 0; y < GPUparams.height; y++) {
                for (int x = 0; x < GPUparams.width; x++) {
                    for (int z = 0; z < GPUparams.chanum - 1; z++) {
                        weights[counter1] = (float) 1.0 / (float) blockvararray[(z + 1) * GPUparams.width * GPUparams.height + y * GPUparams.width + x];
                        varacf[0][x][y][z] = blockvararray[(z + 1) * GPUparams.width * GPUparams.height + y * GPUparams.width + x];
                        sdacf[0][x][y][z] = Math.sqrt((double) varacf[0][x][y][z]);
                        counter1 = counter1 + 1;
                    }
                }
            }

            float[] userinfo2 = new float[numberPoints];
            // System.arraycopy(user_info, 0, userinfo2, 0, numberPoints);
            for (int i = 0; i < numberPoints; i++) {
                userinfo2[i] = user_info[i];
            }

            FitModel fitModel = new FitModel(numberFits, numberPoints, true, model, tolerance, maxNumberIterations, GPUparams.bleachcorr_order, parameters_to_fit, estimator, (GPUparams.chanum - 1) * Float.SIZE / 8);
            fitModel.weights.clear();
            fitModel.weights.put(weights);
            fitModel.userInfo.clear();
            fitModel.userInfo.put(userinfo2);
            fitModel.data.clear();
            fitModel.data.put(data);
            fitModel.initialParameters.clear();
            fitModel.initialParameters.put(initialParameters);

            FitResult fitResult = GpufitImFCS.fit(fitModel);

            boolean[] converged = new boolean[numberFits];
            int numberConverged = 0, numberMaxIterationExceeded = 0, numberSingularHessian = 0, numberNegativeCurvatureMLE = 0;
            for (int i = 0; i < numberFits; i++) {
                FitState fitState = FitState.fromID(fitResult.states.get(i));
                converged[i] = fitState.equals(FitState.CONVERGED);
                switch (fitState) {
                    case CONVERGED:
                        numberConverged++;
                        break;
                    case MAX_ITERATIONS:
                        numberMaxIterationExceeded++;
                        break;
                    case SINGULAR_HESSIAN:
                        numberSingularHessian++;
                        break;
                    case NEG_CURVATURE_MLE:
                        numberNegativeCurvatureMLE++;
                }
            }

            float[] convergedParameterMean = new float[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            // float[] convergedParameterStd = new float[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

            for (int i = 0; i < numberFits; i++) {
                for (int j = 0; j < model.numberParameters; j++) {
                    if (converged[i]) {
                        convergedParameterMean[j] += fitResult.parameters.get(i * model.numberParameters + j);
                    }
                }
            }

            int parfitcounters = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    parfitcounters++;
                }
            }

            double[][][] gpuresultarray = new double[width][height][model.numberParameters];

            // reset
            for (int i = 0; i < width; i++) {
                for (int j = 0; j < height; j++) {
                    chi2[0][i][j] = Double.NaN;
                    pixfitted[0][i][j] = false;
                    pixvalid[0][i][j] = Double.NaN;

                    for (int k = 0; k < model.numberParameters; k++) {
                        gpuresultarray[i][j][k] = Double.NaN;
                    }

                    for (int k = 0; k < noparam; k++) {
                        fitres[0][i][j][k] = Double.NaN;
                    }
                }
            }

            int numfreefitpar = parfitcounters;
            for (int i1 = 0; i1 < GPUparams.width; i1++) {
                for (int i2 = 0; i2 < GPUparams.height; i2++) {
                    chi2[0][i1 + roi1StartX][i2 + roi1StartY] = fitResult.chiSquares.get(i2 * GPUparams.width + i1) / ((fitend - fitstart) - numfreefitpar - 1);

                    if (converged[i2 * GPUparams.width + i1]) {
                        pixfitted[0][i1 + roi1StartX][i2 + roi1StartY] = true;
                        pixvalid[0][i1 + roi1StartX][i2 + roi1StartY] = 1.0;

                        for (int j = 0; j < model.numberParameters; j++) {
                            gpuresultarray[i1 + roi1StartX][i2 + roi1StartY][j] = (double) fitResult.parameters.get((i2 * GPUparams.width + i1) * model.numberParameters + j);
                        }
                    }
                }
            }

            int LBoundX = roi1StartX;
            int UBoundX = roi1StartX + GPUparams.width;
            int LBoundY = roi1StartY;
            int UBoundY = roi1StartY + GPUparams.height;

            for (int i = 0; i < width; i++) {
                for (int j = 0; j < height; j++) {
                    for (int k = 0; k < noparam; k++) {
                        if (i >= LBoundX && i < UBoundX && j >= LBoundY && j < UBoundY) {
                            if (filterArray[i * pixbinX][j * pixbinY] >= filterLL * binningX * binningY && filterArray[i * pixbinX][j * pixbinY] <= filterUL * binningX * binningY) {
                                fitres[0][i][j][k] = gpuresultarray[i][j][k];
                            }
                        }
                    }
                }
            }

            fitModel.reset();
            fitResult.reset();
        } catch (NumberFormatException e) {
            IsGPUCalculationOK = false;
            e.printStackTrace(System.out);
        }

        return IsGPUCalculationOK;
    }

    public boolean GPU_Calculate_ACF_All(Roi improi) {

        boolean IsGPUCalculationOK;

        if (!"none".equals(bleachCorMem) && !"Polynomial".equals(bleachCorMem)) {
            System.out.println("Calculation on GPU is currently supported only for Polynomial bleach correction.");
            return false;
        }

        if (cfXDistance < 0 || cfYDistance < 0) {
            System.out.println("Calculations on GPU current supports only non negative CF X and CF Y distance currently.");
            return false;
        }

        // ********************************************************************************************************************
        // Parameters
        // ********************************************************************************************************************
        int framediff = lastframe - firstframe + 1;

        Rectangle imprect = improi.getBounds();
        int startXmap = (int) Math.ceil(imprect.getX() / pixbinX);
        int startYmap = (int) Math.ceil(imprect.getY() / pixbinY);
        int endXmap = (int) Math.floor((imprect.getX() + imprect.getWidth() - binningX) / pixbinX);
        int endYmap = (int) Math.floor((imprect.getY() + imprect.getHeight() - binningY) / pixbinY);

        int startX = startXmap * pixbinX;
        int startY = startYmap * pixbinY;
        int endX = endXmap * pixbinX;
        int endY = endYmap * pixbinY;

        tbFCCSDisplay.setSelected(false);
        tbFCCSDisplay.setText("Off");
        filterArray = new float[width][height]; // calculate the mean image of the stack

        if (cbFilter.getSelectedItem() == "Mean" || cbFilter.getSelectedItem() == "Intensity") {
            initializeFitres(3, width, height, noparam); // reset the fitresult array and thus the parameter window
            if (cbFilter.getSelectedItem() == "Mean") {
                for (int x1 = startX; x1 <= endX; x1 = x1 + pixbinX) {
                    for (int x2 = startY; x2 <= endY; x2 = x2 + pixbinY) {
                        for (int x3 = firstframe; x3 <= lastframe; x3++) {
                            for (int x4 = 0; x4 < binningX; x4++) {
                                for (int x5 = 0; x5 < binningY; x5++) {
                                    if (improi.contains(x1, x2) && improi.contains(x1, x2 + binningY - 1) && improi.contains(x1 + binningX - 1, x2) && improi.contains(x1 + binningX - 1, x2 + binningY - 1)) {
                                        filterArray[x1][x2] += imp.getStack().getProcessor(x3).get(x1 + x4, x2 + x5);
                                    } else {
                                        filterArray[x1][x2] = Float.NaN;
                                    }
                                }
                            }
                        }
                        filterArray[x1][x2] /= framediff;
                    }
                }
            } else { // if "Intensity" was selcted, then get the first frame and bin if necessary
                for (int x1 = startX; x1 <= endX; x1 = x1 + pixbinX) {
                    for (int x2 = startY; x2 <= endY; x2 = x2 + pixbinY) {
                        for (int x3 = 0; x3 < binningX; x3++) {
                            for (int x4 = 0; x4 < binningY; x4++) {
                                if (improi.contains(x1, x2) && improi.contains(x1, x2 + binningY - 1) && improi.contains(x1 + binningX - 1, x2) && improi.contains(x1 + binningX - 1, x2 + binningY - 1)) {
                                    filterArray[x1][x2] += imp.getStack().getProcessor(firstframe).get(x1 + x3, x2 + x4);
                                } else {
                                    filterArray[x1][x2] = Float.NaN;
                                }
                            }
                        }
                    }
                }
            }
        } else { // if "none" was selected
            for (int x1 = startX; x1 <= endX; x1 = x1 + pixbinX) {
                for (int x2 = startY; x2 <= endY; x2 = x2 + pixbinY) {
                    for (int x3 = 0; x3 < binningX; x3++) {
                        for (int x4 = 0; x4 < binningY; x4++) {
                            if (improi.contains(x1, x2) && improi.contains(x1, x2 + binningY - 1) && improi.contains(x1 + binningX - 1, x2) && improi.contains(x1 + binningX - 1, x2 + binningY - 1)) {
                                filterArray[x1][x2] += imp.getStack().getProcessor(firstframe).get(x1 + x3, x2 + x4);
                            } else {
                                filterArray[x1][x2] = Float.NaN;
                            }
                        }
                    }
                }
            }
        }

        if (isdlawcalculatedingpu == 0) {
            IJ.showProgress(5, 100); // 5%
        }

        // Object to store some of input values for GPU calculations
        GpufitImFCS.ACFParameters GPUparams = new ACFParameters();
        IsGPUCalculationOK = GPU_Initialize_GPUparams(GPUparams, false);

        assert GPUparams.width == (endXmap - startXmap + 1) : "Invalid GPUparams.width";
        assert GPUparams.height == (endYmap - startYmap + 1) : "Invalid GPUparams height";

        if (!GpufitImFCS.isACFmemorySufficient(GPUparams)) {
            IsGPUCalculationOK = false;
            System.out.println("Insufficient GPU memory. Cannot perform ACF calculations on GPU.");
        }

        // ********************************************************************************************************************
        // Get pixels.
        // ******************************************************************************************************************** 
        float[] pixels = new float[GPUparams.w_temp * GPUparams.h_temp * GPUparams.framediff];
        if (IsGPUCalculationOK) {
            IsGPUCalculationOK = GPU_get_pixels(GPUparams, pixels, false);
        }

        if (isdlawcalculatedingpu == 0) {
            IJ.showProgress(10, 100);
        }

        // ********************************************************************************************************************
        // Bleach correction
        // ********************************************************************************************************************   
        double[] bleachcorr_params = new double[GPUparams.w_temp * GPUparams.h_temp * GPUparams.bleachcorr_order];
        if (GPUparams.bleachcorr_gpu && IsGPUCalculationOK) {
            IsGPUCalculationOK = GPU_Calculate_BleachCorrection(GPUparams, pixels, bleachcorr_params);

            if (isdlawcalculatedingpu == 0) {
                IJ.showProgress(40, 100);
            }
        }

        // ********************************************************************************************************************
        // Calculate ACF
        // ********************************************************************************************************************  
        // pixels1 is the output array in which the auto and cross-calculations values are stored.
        double[] pixels1 = new double[GPUparams.width * GPUparams.height * GPUparams.chanum];
        double[] blockvararray = new double[GPUparams.width * GPUparams.height * GPUparams.chanum];

        if (IsGPUCalculationOK) {

            // N&B Calculations. Note cfXdistance and cfYdistance should be zero.
            double[] NBmeanGPU = new double[GPUparams.width * GPUparams.height];
            double[] NBcovarianceGPU = new double[GPUparams.width * GPUparams.height];

            try {
                IsGPUCalculationOK = GPU_Calculate_ACF(pixels, pixels1, blockvararray, NBmeanGPU, NBcovarianceGPU, bleachcorr_params, GPUparams);

                if (!IsGPUCalculationOK) {
                    throw new Exception("Error in ACF calculation.");
                }

                // clear acf array
                for (int i = 0; i < width; i++) {
                    for (int j = 0; j < height; j++) {
                        for (int k = 0; k < chanum; k++) {
                            acf[0][i][j][k] = 0.0;
                        }
                    }
                }

                for (int z = 0; z < GPUparams.chanum; z++) {
                    for (int y = 0; y < GPUparams.height; y++) {
                        for (int x = 0; x < GPUparams.width; x++) {
//                            acf[0][x][y][z] = pixels1[z * GPUparams.width * GPUparams.height + y * GPUparams.width + x];
                            acf[0][x + roi1StartX][y + roi1StartY][z] = pixels1[z * GPUparams.width * GPUparams.height + y * GPUparams.width + x];
                        }
                    }
                }

                if (doMSD) {
                    for (int x = 0; x < width; x++) {
                        for (int y = 0; y < height; y++) {
                            if (!MSDmode) { // 2D if MSDmode is false, otherwise 3D                
                                msd[0][x][y] = correlationToMSD(acf[0][x][y], pixeldimx * Math.pow(10, 6), psfsize * Math.pow(10, 6));
                            } else {
                                msd[0][x][y] = correlationToMSD3D(acf[0][x][y], pixeldimx * Math.pow(10, 6), psfsize * Math.pow(10, 6), lsthickness * Math.pow(10, 6));
                            }
                        }
                    }
                }

                if (isdlawcalculatedingpu == 0) {
                    IJ.showProgress(70, 100);
                }

            } catch (Exception e) {
                e.printStackTrace(System.out);
            }

        }

        // ********************************************************************************************************************
        // ACF_1D fit
        // ********************************************************************************************************************       
        if (IsGPUCalculationOK) {
            if (doFit) {
                IsGPUCalculationOK = GPU_ACF_Fit(GPUparams, pixels1, blockvararray);
            }

            if (isdlawcalculatedingpu == 0) {
                IJ.showProgress(100, 100);
            }

            if (IsGPUCalculationOK) {

                if (isdlawcalculatedingpu == 0) {
                    IJ.showStatus("Plotting data.");
                    plotCF(improi, 2, false);
                    calcAverageIntensityTrace(filterArray, startX, startY, endX, endY, firstframe, lastframe);
                    plotIntensityTrace(improi, 2);
                    // plot MSD if selected
                    if (doMSD) {
                        plotMSD(improi, 2, false);
                    }
                    // create parameter map window
                    if (doFit) {
                        createParaImp(maxcposx - mincposx + 1, maxcposy - mincposy + 1);
//                        createParaImp(GPUparams.width, GPUparams.height);
                    }
                }
            }
        }

        try {
            GpufitImFCS.resetGPU();
        } catch (Exception e) {
            System.out.println("Unable to reset GPU");
        }

        return IsGPUCalculationOK;
    }
}
