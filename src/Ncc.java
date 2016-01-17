import auxiliary.Interpolator;
import ij.ImageStack;

/** */
public class Ncc extends Thread {

    private int begN, endN;

    // INPUT
    public static float[][][][] xy2;                        // N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)
    public static float[][][][] vxy2;                        // N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)

    public static float tanSampling = Float.NaN;            // tangential sampling step (inherited from Delineator)
    public static float[][] inimg_xy;                    // input image (necessary for feature extraction)

    public static float cross_sigma_ratio = Float.NaN;            // necessary for the fitter template defeinition
    public static int patchRadLen = Integer.MIN_VALUE;    // radial
    public static int patchTanLen = Integer.MIN_VALUE;    // cross
    public static int patchTanLenHalf = Integer.MIN_VALUE;// cross

    // OUTPUT
    public static float[][] scores;                    // N(foreground locs.) x 4(max. threads) patch fitting ncc score

    public Ncc(int n0, int n1) {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(
            float[][][][] _xy2,
            float[][][][] _vxy2,
            int _patchRadLen,
            int _patchTanLen,
            float _tanSampling,
            float[][] _inimg_xy,
            float _sigma_ratio
    ) {
        xy2 = _xy2;
        vxy2 = _vxy2;

        tanSampling = _tanSampling;
        inimg_xy = _inimg_xy;

        cross_sigma_ratio = _sigma_ratio;
        patchRadLen = _patchRadLen;   // how many samples will be taken radially
        patchTanLen = _patchTanLen;   // how many along the cross section
        patchTanLenHalf = (_patchTanLen - 1) / 2; // patchTanLen = 2 * patchTanLenHalf + 1

        // allocate output
        scores = new float[xy2.length][4];            // N(foreground locs.)   x 4 (branches)

    }

    public static ImageStack getTemplates()
    {
        Fitter fitter  			= new Fitter(patchRadLen, patchTanLen, cross_sigma_ratio, false);
        return fitter.getTemplates();
    }

    public void run() {

        //*** auxiliary
        float[] ptch = new float[patchRadLen*patchTanLen];
        Fitter fitter  			= new Fitter(patchRadLen, patchTanLen, cross_sigma_ratio, false); // dimRadial, dimTangen = patch dimensions, verbose = false
//		new ImagePlus("thread"+this.getId(), fitter.getTemplates()).show();
        //*** auxiliary

        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {


            /**************************************************************************
             scores[locationIdx]
             **************************************************************************/
            if (xy2[locationIdx]!=null) {

                for (int b = 0; b<xy2[locationIdx].length; b++) {

                    if (xy2[locationIdx][b]!=null) {

                        // fill the values of the patch for branch b
                        for (int radialIdx=0; radialIdx<patchRadLen; radialIdx++) {  // loop radially, will correspond to row
                            for (int shift=-patchTanLenHalf; shift<=patchTanLenHalf; shift++) {

                                // sample the values using refined locs and refined vecs
                                float at_x = xy2[locationIdx][b][0][radialIdx] + shift * tanSampling * vxy2[locationIdx][b][0][radialIdx];
                                float at_y = xy2[locationIdx][b][1][radialIdx] + shift * tanSampling * vxy2[locationIdx][b][1][radialIdx];

                                // radialIdx and tangenIdx mark the xy position within the patch
                                int tangenIdx = shift + patchTanLenHalf;
                                ptch[radialIdx*patchTanLen+tangenIdx] = Interpolator.interpolateAt(at_x, at_y, inimg_xy);
                            }
                        }

                        // patch[] is formed, 1d array stores 2d patch row-by-row

                        // normalize it before fitting
                        //Stat.normalize(ptch); // stays in ptch variable

                        // normalize row by row (due to the patch)  patchRadLen x 2*patchTanLenHalf+1
                        for (int radIdx = 0; radIdx < patchRadLen; radIdx++) {

                            // max per tangential row
                            float row_max = Float.NEGATIVE_INFINITY;
                            for (int tanIdx = 0; tanIdx < patchTanLen; tanIdx++) {
                                if (ptch[radIdx*patchTanLen+tanIdx]>row_max) {
                                    row_max = ptch[radIdx*patchTanLen+tanIdx];
                                }
                            }

                            // normalize
                            for (int tanIdx = 0; tanIdx < patchTanLen; tanIdx++) {
                                ptch[radIdx*patchTanLen+tanIdx] /= row_max;
                            }

                        }

                        // calculate the ncc fit

                        float[]	fit_idx_score = fitter.fit(ptch); // returns [0] - profile index, [1] - fit score
                        scores[locationIdx][b] = fit_idx_score[1];

                    }
                    else {
                        scores[locationIdx][b] = Float.NaN;
                    }
                }
            }
            else {
                for (int i = 0; i < scores[locationIdx].length; i++) scores[locationIdx][i] = Float.NaN;

            }

        }

    }

}
