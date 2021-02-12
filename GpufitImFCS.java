package gpufitImFCS;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Field;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.nio.file.Files;

/**
 * NOTE: Code amended from original gpufit code to put all files into a single
 * Java file for ease of distribution or publishing the code online.
 *
 * Java binding for Gpufit, a Levenberg Marquardt curve fitting library written
 * in CUDA See https://github.com/gpufit/Gpufit,
 * http://gpufit.readthedocs.io/en/latest/bindings.html#java
 *
 * Mirror of the C interface of Gpufit. Loads the Gpufit library and calls it
 * via JNI. On the Java side the arguments to Gpufit are bundled in a FitModel,
 * the results of the fit are bundled in a FitResult. The data arrays in
 * FitModel and FitResult are pre-allocated by us and need to be filled by the
 * user of this library accordingly.
 *
 * See the documentation of FitModel and FitResult as well as the examples.
 *
 * The ACF_1D function is a generalized function which can be used for fitting
 * auto and cross-correlation functions. The LINEAR_1D function is for
 * polynomial bleach correction fitting.
 */
public class GpufitImFCS {

    // NOTE: 
    // VERSION of the used Gpufit library.
    // VERSION must be updated when .dll/.so files are changed so that they are placed in a new sub-folder named after this VERSION num in Fiji.App > jars.
    public static final String VERSION = "v1_0_0";
    public static boolean IsOperatingSystemOK = false;
    public static boolean IsLoadingLibraryOK = false;
    private static final String CUBLASDLLFN = "cublas64_92.dll";
    public static String ALERT;

    static {
        boolean IsWindows = false;
        boolean proceed = true;
        boolean hasError = false;
        String thisAlert = "NVIDIA GPU is not detected.";
        String libName = null;
        boolean writeToFijiAppJars = true; // For testing. Switch to false if we want the files only in temporary folder.

        try {
            // determine OS
            String osname = System.getProperty("os.name");
            String osnamelc = osname.toLowerCase();
            boolean WriteToTempDir = true;

            // reference1: https://stackoverflow.com/questions/4691095/java-loading-dlls-by-a-relative-path-and-hide-them-inside-a-jar
            // reference2: https://stackoverflow.com/questions/2937406/how-to-bundle-a-native-library-and-a-jni-library-inside-a-jar 
            // NOTE: the "/" before the file name. It indicates the current directory.  
            String libDirPath = null;

            // reference: https://www.mkyong.com/java/how-to-detect-os-in-java-systemgetpropertyosname/                        
            if (osnamelc.contains("win")) {
                libName = "agpufit.dll";
                IsWindows = true;
            } else if ((osnamelc.contains("nix") || osnamelc.contains("nux")) || osnamelc.contains("aix")) {
                libName = "libagpufit.so";
            } else {
                proceed = false;
                thisAlert = "GPU mode is currently supported on Windows and Linux, with NVIDIA GPU.";
            }

            if (proceed) {
                IsOperatingSystemOK = true;

                File curr_dir = new File(System.getProperty("java.class.path"));
                File Fiji_jars_dir = curr_dir.getAbsoluteFile().getParentFile();

                if (writeToFijiAppJars && Fiji_jars_dir.canRead()) {

                    // check if gpufitImFCS-cublas folder is present
                    File tempdir = new File(Fiji_jars_dir.toString() + "/gpufitImFCS-cublas");
                    if (!(tempdir.exists() && tempdir.isDirectory())) {
                        tempdir.mkdir();
                    }

                    File gpufitImFCS_cublas_dir = new File(tempdir.toString() + "/" + VERSION);
//                    File gpufitImFCS_cublas_dir = new File(Fiji_jars_dir.toString() + "/gpufitImFCS-cublas/" + VERSION);
                    libDirPath = gpufitImFCS_cublas_dir.toString();

                    boolean Write_Cublas = IsWindows;
                    boolean Write_library = true;

                    if (gpufitImFCS_cublas_dir.exists() && gpufitImFCS_cublas_dir.isDirectory()) {
                        // Directory exists. Check if the libraries are in the folder.
                        if (IsWindows) {
                            Write_Cublas = !(new File(gpufitImFCS_cublas_dir.toString() + "/" + CUBLASDLLFN).exists());
                        }
                        Write_library = !(new File(gpufitImFCS_cublas_dir.toString() + "/" + libName).exists());
                    } else {
                        gpufitImFCS_cublas_dir.mkdir();
                    }

                    if (Fiji_jars_dir.canWrite()) {
                        WriteToTempDir = false;
                        if (Write_Cublas) {
                            writeLibraryFile(gpufitImFCS_cublas_dir, CUBLASDLLFN, false);
                        }

                        if (Write_library) {
                            writeLibraryFile(gpufitImFCS_cublas_dir, libName, false);
                        }
                    } else {
                        WriteToTempDir = Write_Cublas || Write_library;
                    }

                }

                if (WriteToTempDir) {
                    // write files to temporary folder.
                    File tmpDir = Files.createTempDirectory("gpufitImFCS-lib").toFile();
                    tmpDir.deleteOnExit();

                    libDirPath = tmpDir.toString();

                    if (IsWindows) {
                        writeLibraryFile(tmpDir, CUBLASDLLFN, true);
                    }
                    writeLibraryFile(tmpDir, libName, true);
                }

                if (IsWindows) {
                    try {
                        System.load(libDirPath + "/" + CUBLASDLLFN);
                    } catch (Exception e) {
                        thisAlert = "Unable to load " + CUBLASDLLFN + ". Please install CUDA Toolkit 9.2 or newer if it is not present on the computer.";
                        System.out.println(thisAlert);
                    }
                }
                // Proceed to load gpufit.dll. It may still work if user has Cuda Toolkit installed.
                System.load(libDirPath + "/" + libName);
                IsLoadingLibraryOK = true;
            }

        } catch (FileNotFoundException ex) {
            hasError = true;
            System.out.println("FileNotFoundException Error.");
//            ex.printStackTrace(System.out);
//            Logger.getLogger(GpufitImFCS.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            hasError = true;
            System.out.println("IOException Error");
//            ex.printStackTrace(System.out);
//            Logger.getLogger(GpufitImFCS.class.getName()).log(Level.SEVERE, null, ex);
        } catch (UnsatisfiedLinkError ex) {
            hasError = true;
            System.out.println("Unsatisfied Link Error");
//            ex.printStackTrace(System.out);
//            Logger.getLogger(GpufitImFCS.class.getName()).log(Level.SEVERE, null, ex);
        }

        // if OS is OK but encountered error loading libraries.
        if (proceed && hasError) {
            thisAlert = "Unable to load " + libName + ".";
            System.out.println(thisAlert);
        }
        ALERT = thisAlert;
    }

    private static void writeLibraryFile(File Directory, String LibName, Boolean DeleteOnExit) {
        try {
            //NOTE: include package name, which becomes the folder name in .jar file.'
            InputStream in = ClassLoader.class.getResourceAsStream("/gpufitImFCS/" + LibName);
            if (in == null) {
                throw new FileNotFoundException("Library " + LibName + " is not available");
            }

            File temp = new File(Directory, LibName);
            if (DeleteOnExit) {
                temp.deleteOnExit();
            }

            byte[] buffer = new byte[1024];
            try (FileOutputStream fos = new FileOutputStream(temp)) {
                int read;
                //NOTE: other methods didn't work. the read = -1 seems important. Examples of other methods:
                // 1) length = is.read(buffer); os.write(buffer, 0, length);
                // 2) while ((readBytes = stream.read(buffer)) > 0) { resStreamOut.write(buffer, 0, readBytes);}
                // 3) try (InputStream in = url.openStream()) { Files.copy(in, nativeLibTmpFile.toPath());}
                // 4) InputStream source = ClassLoader.class.getResourceAsStream("/libagpufit.so"); File outfile = new File(tempdir, "libagpufit.so"); FileOutputStream fos = new java.io.FileOutputStream(outfile); while (source.available() > 0) {  // write contents of 'is' to 'fos' fos.write(source.read());}
                while ((read = in.read(buffer)) != -1) {
                    fos.write(buffer, 0, read);
                }
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    /**
     * -------------------------------------------------------------------------------------------------
     * from GpufitUtils.java START /**
     * -------------------------------------------------------------------------------------------------
     */
    /**
     * Creates a direct ByteBuffer with the native Byte order, which is exactly
     * what we need to send and receive data arrays via JNI.
     *
     * @param length Desired number of elements in ByteBuffer
     * @return Direct ByteBuffer
     */
    public static ByteBuffer allocateDirectByteBuffer(int length) {
        assertNotNegative(length, "Parameter length must be non-negative");
        ByteBuffer buffer = ByteBuffer.allocateDirect(length);
        buffer.order(ByteOrder.nativeOrder());
        return buffer;
    }

    /**
     * Creates a FloatBuffer of a certain length backed up by a direct
     * ByteBuffer.
     *
     * @param length Desired number of elements in FloatBuffer
     * @return FloatBuffer backed by direct ByteBuffer
     */
    public static FloatBuffer allocateDirectFloatBuffer(int length) {
        return allocateDirectByteBuffer(Float.BYTES * length).asFloatBuffer();
    }

    /**
     * Creates an IntBuffer of a certain length backed up by a direct
     * ByteBuffer.
     *
     * @param length Desired number of elements in IntBuffer
     * @return IntBuffer backed by direct ByteBuffer
     */
    public static IntBuffer allocateDirectIntBuffer(int length) {
        return allocateDirectByteBuffer(Integer.BYTES * length).asIntBuffer();
    }

    /**
     * Ensures that an integer value is not negative. Throws a runtime exception
     * otherwise.
     *
     * @param value Integer value
     * @param message Error message
     */
    public static void assertNotNegative(int value, String message) {
        if (!(value >= 0)) {
            throw new RuntimeException(message);
        }
    }

    /**
     * Ensures that a boolean is true. Throws a runtime exception otherwise.
     *
     * @param value Boolean value
     * @param message Error message
     */
    public static void assertTrue(boolean value, String message) {
        if (!value) {
            throw new RuntimeException(message);
        }
    }

    /**
     * Ensures that an object reference passed as a parameter to the calling
     * method is not null.
     *
     * @param <T> Type of the reference to be checked
     * @param reference Object reference
     * @return Non-null reference that was checked
     */
    public static <T> T verifyNotNull(T reference) throws NullPointerException {
        if (null == reference) {
            throw new NullPointerException();
        }
        return reference;
    }

    /**
     * Adds a path to the java.library.path programmatically.
     *
     * See also:
     * https://stackoverflow.com/questions/11783632/how-do-i-load-and-use-native-library-in-java
     *
     * @param path Path String to add to the java.library.path
     */
    public static void addPathToJavaLibraryPath(String path) {

        // add path to system property
        String libraryPath = System.getProperty("java.library.path");
        libraryPath += File.pathSeparator + path;
        System.setProperty("java.library.path", libraryPath);

        // clear field in class loader
        try {
            Field fieldSysPath = ClassLoader.class.getDeclaredField("sys_paths");
            fieldSysPath.setAccessible(true);
            fieldSysPath.set(null, null);
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new RuntimeException("Java Library Path addition failed.");
        }
    }

    /**
     * -------------------------------------------------------------------------------------------------
     * from GpufitUtils.java START /**
     * -------------------------------------------------------------------------------------------------
     */
    /**
     * -------------------------------------------------------------------------------------------------
     * from Gpufit.java START /**
     * -------------------------------------------------------------------------------------------------
     */
    /**
     * Native method. More of less calls gpufit() in the gpufit C interface
     * directly. Used only internally.
     */
    private static native int fit(int numberFits, int numberPoints, FloatBuffer data, FloatBuffer weights, int model_id, FloatBuffer initialParameters, float tolerance, int maxNumberIterations, int num_valid_coefs, IntBuffer parametersToFit, int estimatorID, int userInfoSize, FloatBuffer userInfo, FloatBuffer outputParameters, IntBuffer outputStates, FloatBuffer outputChiSquares, IntBuffer outputNumberIterations);

    /**
     * Use this method to perform a parallel fit of many single fits of the same
     * Function model and the same fit data size in parallel. The input is given
     * as a {@link FitModel}, the output as a {@link FitResult}. If null is
     * passed for fitResult, a new output structure is created, otherwise the
     * results overwrite the existent content (size is checked for compatibility
     * before).
     *
     * A call to Gpufit is performed. If an error is encountered, a runtime
     * exception is thrown.
     *
     * @param fitModel Fit data including the model
     * @param fitResult Fit result (could be old one which is reused) or null
     * @return Fit result
     */
    public static FitResult fit(FitModel fitModel, FitResult fitResult) {

        // Should we reuse fitResult?
        if (null == fitResult) {
            fitResult = new FitResult(fitModel.numberFits, fitModel.model.numberParameters);
        } else {
            // check sizes
            fitResult.isCompatible(fitModel.numberFits, fitModel.model.numberParameters);
        }

        // call into native code to perform a fit
        long t0 = System.currentTimeMillis();
        int status = GpufitImFCS.fit(fitModel.numberFits, fitModel.numberPoints, fitModel.data, fitModel.weights, fitModel.model.id, fitModel.initialParameters, fitModel.tolerance, fitModel.maxNumberIterations, fitModel.num_valid_coefs, fitModel.parametersToFit, fitModel.estimator.id, fitModel.userInfo.capacity(), fitModel.userInfo, fitResult.parameters, fitResult.states, fitResult.chiSquares, fitResult.numberIterations);
        long t1 = System.currentTimeMillis();
        fitResult.fitDuration = (float) (t1 - t0) / 1000;

        // check status
        if (status != Status.OK.id) {
            String message = getLastError();
            throw new RuntimeException(String.format("status = %s, message = %s", status, message));
        }

        // return results
        return fitResult;
    }

    /**
     * Convenience method. Calls fit with only a FitModel. Returns the result.
     *
     * @param fitModel Fit data including the model
     * @return Fit result
     */
    public static FitResult fit(FitModel fitModel) {
        return fit(fitModel, null);
    }

    /**
     * Native method. Returns a string representing the last error message from
     * Gpufit.
     *
     * @return The last error message from Gpufit.
     */
    public static native String getLastError();

    /**
     * Native method. Indicates if CUDA capability is available.
     *
     * @return True if available, false otherwise.
     */
    private static native boolean isCudaAvailableInt();

    public static boolean isCudaAvailable() {
        if (IsOperatingSystemOK && IsLoadingLibraryOK) {
            return isCudaAvailableInt();
        } else {
            return false;
        }
    }

    /**
     * Native method. Gets the CUDA runtime and driver version as ints. Used
     * only internally.
     */
    private static native int[] getCudaVersionAsArray();

    /**
     * Gets the CUDA runtime and driver version as strings.
     *
     * Throws an exception if CUDA is not available.
     *
     * @return The CUDA version.
     */
    public static CudaVersion getCudaVersion() {
        //String[] version = getOpenCLVersionAsArray();
        int[] version = getCudaVersionAsArray();
        if (null == version) {
            String message = getLastError();
            throw new RuntimeException(message);
        }

        //For OpenCL
        //return new CudaVersion(version[0], version[1]);
        //For CUDA
        return new CudaVersion(versionAsString(version[0]), versionAsString(version[1]));
    }

    /**
     * Special conversion for our versions from Integer to String. The
     * convention is that the integer is major version * 1000 + minor version *
     * 10.
     *
     * @param version An integer version.
     * @return A version string.
     */
    private static String versionAsString(int version) {
        return String.format("%d.%d", version / 1000, (version % 1000) / 10);
    }

    /**
     * The status of a call to the Gpufit routines (see Gpufit documentation for
     * details). Used only internally.
     */
    private enum Status {

        OK(0), ERROR(-1);

        final int id;

        Status(int id) {
            this.id = id;
        }
    }

    public static class ACFParameters {

        public int width; //output width
        public int height; //output height
        public int win_star; // w output x pixbinX + cfXdistance
        public int hin_star; // h output x pixbinY + cfYdistance
        public int w_temp; // win_star - pixbinX + 1
        public int h_temp; // hin_start - pixbinY + 1
        public int pixbinX;
        public int pixbinY;
        public int binningX; // binning in X axis
        public int binningY; // binning in Y axis
        public int firstframe;
        public int lastframe;
        public int framediff;
        public int cfXDistance;
        public int cfYDistance;
        public double correlatorp;
        public double correlatorq;
        public double frametime;
        public int background;
        public double mtab1; // mtab[1], used to calculate blocknumgpu.
        public double mtabchanumminus1; // mtab[chanum-1], used to calculate pnumgpu[counter_indexarray]
        public double sampchanumminus1; // samp[chanum-1], used to calculate pnumgpu[counter_indexarray]
        public int chanum;
        public boolean isNBcalculation;
        public boolean bleachcorr_gpu;
        public int bleachcorr_order;
        public int nopit;
        public int ave;

        public ACFParameters() {
        }
    }

    /**
     * Native method. Reset gpu by running cudaDeviceReset function
     */
    public static native void resetGPU();

    /**
     * Native method. Get data bleach correction given the 1D pixels values.
     *
     * @param pixels 1D float array of pixels values
     * @param outdata output data, by reference.
     * @param ACFInputParams Values nopit, ave, cfXDistance, cfYDistance, width
     * and height are required
     */
    public static native void calcDataBleachCorrection(float[] pixels, float[] outdata, ACFParameters ACFInputParams);

    /**
     * Native method. Find if GPU memory is sufficient to perform binning.
     *
     * @param ACFInputParams Object that encapsulates / stores the necessary
     * input values.
     * @return true is memory is sufficient, false otherwise.
     */
    public static native boolean isBinningMemorySufficient(ACFParameters ACFInputParams);

    /**
     * Native method. Do binning given a 1D array.
     *
     * @param indata 1D float array. Dimension win_star = w_out x pixbinX +
     * cfXDistance, hin_star = h_out x pixbinY + cfYDistance
     * @param outdata output data, by reference.
     * @param ACFInputParams Values pixbinX, pixbinY are required.
     */
    public static native void calcBinning(float[] indata, float[] outdata, ACFParameters ACFInputParams);

    /**
     * Native method. Indicates if memory is sufficient to run calcacf2 and
     * calcacf3 calcacf3 will be executed first, then some pointers will be
     * freed before running calcacf2.
     *
     * @param ACFInputParams Object that encapsulates/stores the necessary input
     * values.
     * @return true if memory is sufficient, false otherwise.
     */
    public static native boolean isACFmemorySufficient(ACFParameters ACFInputParams);

    public static native void calcACF(float[] pixels, double[] pixels1, double[] blockvararray,
            double[] NBmeanGPU, double[] NBcovarianceGPU,
            double[] blocked1D, double[] bleachcorr_params,
            double[] samp, int[] lag, ACFParameters ACFInputParams);

    public static native void resetDevice();

    /**
     * -------------------------------------------------------------------------------------------------
     * from Gpufit.java END /**
     * -------------------------------------------------------------------------------------------------
     */
    /**
     * -------------------------------------------------------------------------------------------------
     * from CudaVersion.java START /**
     * -------------------------------------------------------------------------------------------------
     */
    public static class CudaVersion {

        public final String runtime, driver;

        public CudaVersion(String runtime, String driver) {
            this.runtime = runtime;
            this.driver = driver;
        }
    }

    /**
     * -------------------------------------------------------------------------------------------------
     * from CudaVersion.java END /**
     * -------------------------------------------------------------------------------------------------
     */
    /**
     * -------------------------------------------------------------------------------------------------
     * from Enumerator.java START /**
     * -------------------------------------------------------------------------------------------------
     */
    public enum Estimator {

        /**
         * Least-squares estimator
         */
        LSE(0),
        /**
         * Poisson maximum likelihood estimator
         */
        MLE(1);

        public final int id;

        Estimator(int id) {
            this.id = id;
        }

    }

    /**
     * -------------------------------------------------------------------------------------------------
     * from Enumerator.java END /**
     * -------------------------------------------------------------------------------------------------
     */
    /**
     * -------------------------------------------------------------------------------------------------
     * from FitModel.java START /**
     * -------------------------------------------------------------------------------------------------
     *
     * @param numberFits
     * @param numberPoints
     * @param withWeights
     * @param model
     * @param tolerance
     * @param maxNumberIterations
     * @param num_valid_coefs
     * @param parametersToFit
     * @param estimator
     * @param userInfoSize
     * @return *
     */
    // To make public FitModel function.
    public static FitModel FitModel(int numberFits, int numberPoints, boolean withWeights, Model model, Float tolerance, Integer maxNumberIterations, Integer num_valid_coefs, Boolean[] parametersToFit, Estimator estimator, int userInfoSize) {
        return new FitModel(numberFits, numberPoints, withWeights, model, tolerance, maxNumberIterations, num_valid_coefs, parametersToFit, estimator, userInfoSize);
    }

    public static class FitModel {

        /**
         * Number of fits, i.e. number of independent data sets
         */
        public final int numberFits;
        /**
         * Number of data points per fit
         */
        public final int numberPoints;
        /**
         * Buffer holding the data
         */
        public final FloatBuffer data;
        /**
         * Buffer holding the data weights (or null)
         */
        public final FloatBuffer weights;
        /**
         * Fit function model enum
         */
        public final Model model;
        /**
         * Initial value of parameters for each data set
         */
        public final FloatBuffer initialParameters;
        /**
         * Minimal fit tolerance.
         */
        public final float tolerance;
        /**
         * Maximal number of iterations per data set.
         */
        public final int maxNumberIterations;
        /**
         * Max number of coefficients to fit. Used mainly for linear_1d
         * polynomial fit.
         */
        public final int num_valid_coefs;
        /**
         * Indication which parameters should be fitted (value 1) and which
         * should be kept constant (value 0).
         */
        public final IntBuffer parametersToFit;
        /**
         * Fit estimator enum
         */
        public final Estimator estimator;
        /**
         * Additional user info (optional).
         */
        public final FloatBuffer userInfo;

        /**
         * Provide a number of input arguments for the fit.
         *
         * Indicate if weights or userInfo is needed by passing withWeights true
         * and userInfoSize is larger than zero.
         *
         * Sets default values for tolerance, maxNumberIterations,
         * parametersToFit and estimator if these parameters are passed as null.
         *
         * Fill data, weights (if needed), initialParameters and userInfo (if
         * needed) afterwards with values.
         *
         * @param numberFits Number of fits
         * @param numberPoints Number of data points per fit
         * @param withWeights If True, weights will be pre-allocated, otherwise
         * not
         * @param model Fit model enum
         * @param tolerance Fit tolerance (if null, a default value is chosen)
         * @param maxNumberIterations Maximal number of iterations per fit (if
         * null, a default value is chosen)
         * @param num_valid_coefs Max number of coefficients to fit. Used mainly
         * for linear_1d polynomial fit.
         * @param parametersToFit For each parameter indicates if it should be
         * fit (true) or kept constant (false). If null all parameters are fit
         * by default.
         * @param estimator Fit estimator enum (if null, a default value is
         * chosen)
         * @param userInfoSize If positive, userInfo is pre-allocated with
         * userInfoSize as capacity, otherwise not
         */
        public FitModel(int numberFits, int numberPoints, boolean withWeights, Model model, Float tolerance, Integer maxNumberIterations, Integer num_valid_coefs, Boolean[] parametersToFit, Estimator estimator, int userInfoSize) {

            this.numberFits = numberFits;
            this.numberPoints = numberPoints;
            this.data = allocateDirectFloatBuffer(numberFits * numberPoints);
            this.weights = withWeights ? allocateDirectFloatBuffer(numberFits * numberPoints) : null;
            this.model = verifyNotNull(model);
            this.initialParameters = allocateDirectFloatBuffer(numberFits * model.numberParameters);
            this.tolerance = tolerance == null ? 1e-4f : tolerance;
            this.maxNumberIterations = maxNumberIterations == null ? 25 : maxNumberIterations;
            this.num_valid_coefs = num_valid_coefs;
            this.parametersToFit = allocateDirectIntBuffer(model.numberParameters);
            if (null == parametersToFit) {
                // fill with ones
                for (int i = 0; i < model.numberParameters; i++) {
                    this.parametersToFit.put(1);
                }
            } else {
                // fill with given values
                for (int i = 0; i < model.numberParameters; i++) {
                    this.parametersToFit.put(parametersToFit[i] ? 1 : 0);
                }
            }
            this.estimator = estimator == null ? Estimator.LSE : estimator;
            this.userInfo = allocateDirectFloatBuffer(Math.max(0, userInfoSize));
        }

        public void reset() {
            data.clear();
            weights.clear();
            initialParameters.clear();
            parametersToFit.clear();
            userInfo.clear();
        }
    }

    /**
     * -------------------------------------------------------------------------------------------------
     * from FitModel.java END /**
     * -------------------------------------------------------------------------------------------------
     */
    /**
     * -------------------------------------------------------------------------------------------------
     * from FitResult.java START /**
     * -------------------------------------------------------------------------------------------------
     */
    public static class FitResult {

        /**
         * Values of resulting fit parameters for each data set.
         */
        public final FloatBuffer parameters;
        /**
         * Ids of the fit results (@see FitResult)
         */
        public final IntBuffer states;
        /**
         * Final Chi² for each fit
         */
        public final FloatBuffer chiSquares;
        /**
         * Used number of iterations for each fit
         */
        public final IntBuffer numberIterations;
        /**
         * Duration of fit in seconds.
         */
        public float fitDuration;

        /**
         * Given the number of fits and the number of parameters of the fit
         * model, pre-allocates memory for the fit results.
         *
         * @param numberFits Number of fits in the call to Gpufit.fit
         * @param numberParameters Number of parameters of the model
         */
        public FitResult(int numberFits, int numberParameters) {
            parameters = allocateDirectFloatBuffer(numberFits * numberParameters);
            states = allocateDirectIntBuffer(numberFits);
            chiSquares = allocateDirectFloatBuffer(numberFits);
            numberIterations = allocateDirectIntBuffer(numberFits);
        }

        /**
         * Checks where the sizes of the existing member variables is consistent
         * with a given number of fits and number of parameters of the fit
         * model. Required when re-using an instance between different runs of
         * the fit.
         *
         * @param numberFits Number of fits in the call to Gpufit.fit
         * @param numberParameters Number of parameters of the model
         */
        public void isCompatible(int numberFits, int numberParameters) {
            assertTrue(parameters.capacity() == numberFits * numberParameters, "Expected different size of parameters");
            assertTrue(states.capacity() == numberFits, "Expected different size of states");
            assertTrue(chiSquares.capacity() == numberFits, "Expected different size of chiSquares");
            assertTrue(numberIterations.capacity() == numberFits, "Expected different size of numberIterations");
        }

        public void reset() {
            parameters.clear();
            states.clear();
            chiSquares.clear();
            numberIterations.clear();
        }

    }

    /**
     * -------------------------------------------------------------------------------------------------
     * from FitResult.java END /**
     * -------------------------------------------------------------------------------------------------
     */
    /**
     * -------------------------------------------------------------------------------------------------
     * from FitState.java START /**
     * -------------------------------------------------------------------------------------------------
     */
    public enum FitState {

        /**
         * Converged within the given tolerance
         */
        CONVERGED(0),
        /**
         * Maximum number of iterations exceeded
         */
        MAX_ITERATIONS(1),
        /**
         * Singular Hessian
         */
        SINGULAR_HESSIAN(2),
        /**
         * Negative curvature in the MLE
         */
        NEG_CURVATURE_MLE(3),
        /**
         * GPU not ready
         */
        GPU_NOT_READY(4);

        /**
         * Id is the same as the output of the Gpufit fit.
         */
        private final int id;

        FitState(int id) {
            this.id = id;
        }

        /**
         * Retrieves the enum corresponding to a certain id.
         *
         * Throws a runtime exception if the id is unknown.
         *
         * @param id An id of a FitState (for example the output of the fit
         * routines).
         * @return An FitState enum member as defined above.
         */
        public static FitState fromID(int id) {
            for (FitState fitState : values()) {
                if (fitState.id == id) {
                    return fitState;
                }
            }
            throw new RuntimeException("Unknown id");
        }
    }

    /**
     * -------------------------------------------------------------------------------------------------
     * from FitState.java END /**
     * -------------------------------------------------------------------------------------------------
     */
    /**
     * -------------------------------------------------------------------------------------------------
     * from Model.java START /**
     * -------------------------------------------------------------------------------------------------
     * *
     */
    public enum Model {

        // ACF_1D function is a generalized function which can be used for fitting auto and cross-correlation functions.
        // LINEAR_1D function is for polynomial bleach correction fitting.
        GAUSS_2D(1, 5), ACF_1D(2, 20), LINEAR_1D(3, 11);

        public final int id, numberParameters;

        Model(int id, int numberParameters) {
            this.id = id;
            this.numberParameters = numberParameters;
        }
    }
    /**
     * -------------------------------------------------------------------------------------------------
     * from Model.java END /**
     * -------------------------------------------------------------------------------------------------
     */
}
