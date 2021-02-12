
Imaging FCS 1.52 is a basic ImageJ plugin to calculate and view spatio-temporal correlation functions from 16 bit gray tiff stack files. It was written under FIJI (ImageJ 1.51f; Java 1.8.0_102) and requires Imagescience for statistics (simulator) and Apache Poi for file reading and writing.

The details about the sofware are also provided here. 
https://www.dbs.nus.edu.sg/lab/BFL/imfcs_image_j_plugin.html

ImagingFCS 1.52 tries to provide a comprehensive software tool to calculate and evaluate spatiotemporal correlation functions. It includes the calculation of all auto- or cross-correlation functions for arbitrary pixel binning and regions of interest within an image, provides fit functions for total internal reflection fluorescence (TIRF) and single plane illumination microscopy (SPIM) based FCS measurements, can calculate the FCS diffusion laws and contains a basic simulator to create simulated data for different diffusive modes.

ImagingFCS runs under ImageJ, FIJI and Micromanager, and it runs on PC as well as on Mac OS. We will use in the following text always FIJI but it should be understood that the same is true for ImageJ and Micromanager.

Currently we provide two ways to install the plugin.

# Option 1: By using the Update sites of ImageJ/Fiji
Click on Help -> Update. Later click on Manage update sites. Check the boxes next to ImagingFCS and ImageScience. Click on Close. Now click on the Apply Changes button. The plugin will be downloaded. As instructed, please restart ImageJ. The plugin can now be used.

# Option 2: By downloading the files provided here in GITHUB
The following files are needed:
1. Imaging_FCS_1_52.jar : Put this file in the plugin folder of FIJI (“Fiji.app\plugins”). 
2. gpufitImFCS.jar : Put this file in the jars folder of FIJI (“Fiji.app\jars”). Start FIJI. You will find Imaging_FCS_1_52 under the plugin tab.
3. Imagescience: Either install imagescience.jar in the jar folder within FIJI or link the update side to imagescience. This supports the probability distributions used in the simulator.
  http://www.imagescience.org/meijering/software/imagescience/

4. Apache POI :  You need to install Apache poi-3.17 (the latest stable release at the time of writing). The Apache Poi provides the necessary code for the writing and reading of .xlsx spreadsheet files, which are used to store and read experimental data. You can copy the whole poi-3.17 folder into the jars folder of Fiji (\Fiji.app\jars). It has also been found that sometimes there are errors in reading the jar files inside the poi folder. In case, if there are errors while trying to run the plugin and if the error is associated with poi files, one suggestion is to place all the jar files inside the poi folder directly under the jars folder. In total, there must be 13 jar files as per poi 3.17. Six of them are found just inside the folder. Five of them are in the lib folder and another two of them are in the ooxml-lib folder.
  http://poi.apache.org/download.html


If you want to compile the program yourself:
You will then need Imaging_FCS_1_52.java : This is the Java code from which the .jar file was produced. You can open it in FIJI and compile it yourself.
For manual compilation of GpufitImFCS.java and the CUDA code, agpufitjni.cu , please refer to section 1.3.3 in the Imaging FCS_1_52 Manual
1. PC: If you want to compile the program on a PC in ImageJ, you need to install a JDK (e.g. I used Jave SE - jdk1.8.0_102). For details see http://forum.imagej.net/t/no-javac-jar-found/2340 and http://stackoverflow.com/questions/18455732/play-framework-cant-find-javac .
2. Mac: If you would to compile the program yourself you need to install a Java IDE, i.e. an Integrative Development Environment. Netbeans, for instance, is free and worked fine for us. But there are other free IDEs, e.g. Eclipse, JSource, IntelliJ IDEA etc. OR drag and drop Imaging_FCS_1_52.java onto FIJI control panel and make sure gpufitImFCS.jar is in "jars" folder of FIJI. Note that manual compilation of GpufitImFCS.java requires Mac OS with supported NVIDIA graphic cards. We recommend manual compilation of GpufitImFCS.java through Windows or Linux.

# An important note:
Switch off SCIFIO for file opening in ImageJ2 as Imaging_FCS_1_52 does not work with SCIFIO yet. To switch SCIFIO off, go to Edit\Options\ImageJ2 in the Fiji control bar. A dialog will appear in which you can untick the option “Use SCIFIO when opening files”.

# Example File
Bilayer.tif: This is an example tiff stack to test the program. In general a tiff stack should contain at least 20,000 frames which were recorded with a time resolution of 1 ms or less. This is sufficient to resolve the dynamics in lipid bilayers. However, we recommend to take at least 50,000 frames for better statistics. For faster processes, shorter frame times and more frames are required (see Sankaran et al. Analytical Chemistry 2013).

# Manual
This manual contains the basic instructions how to use the program, the definition of all items in the control and fit panels, the file formats of the saved data, and the theoretical functions used for fitting.

# References:
# Correlator Scheme
1. Meseth, U., T. Wohland, R. Rigler, and H. Vogel. 1999. "Resolution of fluorescence correlation measurements." Biophys J. 76: 1619–1631.
2. Shi, X., and T. "Wohland. 2010. Fluorescence Correlation Spectroscopy." In: Diaspro A, editor. Nanoscopy and Multidimensional Optical Fluorescence Microscopy. Boca Raton: CRC Press. pp. 6.1–6.34.
3. Sankaran, J., X. Shi, L.Y. Ho, E.H.K. Stelzer, and T. Wohland. 2010. "ImFCS: a software for imaging FCS data analysis and visualization." Opt Express. 18: 25468–25481.
# Imaging FCS publications
1. Kannan, B., J. Har, P. Liu, I. Maruyama, J. Ding, et al. 2006. "Electron multiplying charge-coupled device camera based fluorescence correlation spectroscopy." Anal Chem. 78: 3444–3451.
2. Kannan, B., L. Guo, T. Sudhaharan, S. Ahmed, I. Maruyama, et al. 2007. "Spatially resolved total internal reflection fluorescence correlation microscopy using an electron multiplying charge-coupled device camera." Anal Chem. 79: 4463–4470.
3 Guo, L., J.Y. Har, J. Sankaran, Y. Hong, B. Kannan, et al. 2008. "Molecular Diffusion Measurement in Lipid Bilayers over Wide Concentration Ranges: A Comparative Study." ChemPhysChem. 9: 721–728.
4. Sankaran, J., M. Manna, L. Guo, R. Kraut, and T. Wohland. 2009. "Diffusion, Transport, and Cell Membrane Organization Investigated by Imaging Fluorescence Cross-Correlation Spectroscopy." Biophys J. 97: 2630–2639.
5. Wohland, T., X. Shi, J. Sankaran, and E. Stelzer. 2010. "Single plane illumination fluorescence correlation spectroscopy (SPIM-FCS) probes inhomogeneous three-dimensional environments." Opt Express. 18: 10627–10641.
6. Kraut, R., N. Bag, and T. Wohland. 2012. "Methods in Cell Biology." Elsevier.
7. Bag, N., J. Sankaran, A. Paul, R.S. Kraut, and T. Wohland. 2012. "Calibration and Limits of Camera-Based Fluorescence Correlation Spectroscopy: A Supported Lipid Bilayer Study." ChemPhysChem. 13: 2784–2794.
8. Sankaran, J., N. Bag, R.S. Kraut, and T. Wohland. 2013. "Accuracy and Precision in Camera-Based Fluorescence Correlation Spectroscopy Measurements." Anal Chem. : 130404142328008.
9. Bag, N., A. Ali, V.S. Chauhan, T. Wohland, and A. Mishra. 2013. "Membrane destabilization by monomeric hIAPP observed by imaging fluorescence correlation spectroscopy." Chem. Commun. 49: 9155.
10. Singh, A.P., J.W. Krieger, J. Buchholz, E. Charbon, J. Langowski, et al. 2013. "The performance of 2D array detectors for light sheet based fluorescence correlation spectroscopy." Opt Express. 21: 8652.
11. Bag N., D.H.X. Yap, T. Wohland, "Temperature dependence of diffusion in model and live cell membranes characterized by imaging fluorescence correlation spectroscopy." BBA Biomembranes 1838 (2014) 802–813.
12. Guo, S.M., N. Bag, A. Mishra, T. Wohland, and M. Bathe, "Bayesian Total Internal Reflection Fluorescence Correlation Spectroscopy Reveals hIAPP-Induced Plasma Membrane Domain Organization in Live Cells." Biophy. J. 2014, (106) 190–200.
13. Krieger, J.W., A.P. Singh, C.S. Garbe, T. Wohland and J. Langowski, “Dual-Color Fluorescence Cross-Correlation Spectroscopy on a Single Plane Illumination Microscope (SPIM-FCCS)”, Optics Express (2014), 22(3), 2358–2375.
14. Bag, N; Huang, S; Wohland, T. “Plasma Membrane Organization of Epidermal Growth Factor Receptor in Resting and Ligand-Bound States.” Biophys J. 2015; 109(9):1925-36.
15. Machan R, Foo YH, Wohland T. "On the Equivalence of FCS and FRAP: Simultaneous Lipid Membrane Measurements." Biophys. J. 2016, (111) 152–161.
16. Bag N, Ng XW, Sankaran J, Wohland T. "Spatiotemporal mapping of diffusion dynamics and organization in plasma membranes." Methods Appl. Fluoresc. 4 (2016) 034003.
17. Ng XW, Teh C., Korzh V, and Wohland T. "The Secreted Signaling Protein Wnt3 Is Associated with Membrane Domains In Vivo: A SPIM-FCS Study." Biophy.J. 2016, (111) 418–429.

# Imaging FCS reviews
1. N. Bag and T. Wohland, “Imaging Fluorescence Fluctuation Spectroscopy: New Tools for Quantitative Bioimaging”. Ann. Rev. Phys. Chem. (2014) 65:225–248.
2. A.P. Singh and T. Wohland, “Applications of imaging fluorescence correlation spectroscopy”, Curr. Op. Chem. Biol. (2014) 20:29–35.
Gpufit
3. A. Przybylski, B.X.R. Thiel, J. Keller-Findeisen, B. Stock, and M. Bates, “Gpufit: An open-source toolkit for GPU-accelerated curve fitting.”, Sci. Rep. (2017) 1–9.

# Disclaimer
The software and data on this site are provided for personal or academic use only and may not be used in any commercial venture or distributions. All files have been virus scanned, however, for your own protection; you should scan these files again. You assume the entire risk related to your use of this software and data. By using the software and data on this site your expressly assume all risks of data loss or damage alleged to have been caused by the software and data. The Biophysical Fluorescence Laboratory at NUS is providing this data "as is," and disclaims any and all warranties, whether express or implied, including (without limitation) any implied warranties of merchantability or fitness for a particular purpose. In no event will the Biophysical Fluorescence Laboratory at NUS and/or NUS be liable to you or to any third party for any direct, indirect, incidental, consequential, special or exemplary damages or lost profit resulting from any use or misuse of this software and data.


